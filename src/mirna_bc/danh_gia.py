"""Khung đánh giá mô hình sống còn cho ND4 (quyết định D1, 11/2026). Chỉ chạy ở ND4, không dùng ở ND1.

- Đánh giá nội bộ chính: bootstrap hiệu chỉnh lạc quan theo Harrell (200 lần, khớp lại TOÀN BỘ pipeline mỗi lần).
- Phân tích phụ: repeated stratified 5-fold CV (20 lần), chỉ số tính trên dự đoán ngoài fold gộp theo từng lần lặp.
- Chỉ số: Harrell C, Uno C cắt tại tau = 3 năm, AUC phụ thuộc thời gian tại 2 và 3 năm.
- Mô hình chính: CoxNet (elastic-net, l1_ratio = 0,5), chọn alpha theo quy tắc 1-SE trên CV lặp ở vòng trong;
  chuẩn hóa z chỉ nằm bên trong pipeline.
- Đối chứng null: cặp ngẫu nhiên (giữ miRNA, thay gen bằng gen không phải đích cùng thập phân vị biểu hiện).
- QC miRBind2: đích dự đoán có lệch về phía tương quan âm so với gen không phải đích (KS một phía).

Thời gian tính bằng NGÀY (giống TCGA-CDR). y là mảng có cấu trúc của scikit-survival (event, time).
"""
import warnings
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
from statsmodels.stats.multitest import multipletests
from sklearn.base import BaseEstimator, clone
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored, concordance_index_ipcw, cumulative_dynamic_auc
from sksurv.util import Surv

NAM = 365.25
TAU = 3 * NAM
MOC_AUC = (2 * NAM, 3 * NAM)


def tao_y(su_kien, thoi_gian):
    """Mảng y cho scikit-survival từ cột sự kiện (0/1) và thời gian (ngày)."""
    return Surv.from_arrays(np.asarray(su_kien).astype(bool), np.asarray(thoi_gian, float))


# ---------------------------------------------------------------------------
# Chỉ số
# ---------------------------------------------------------------------------
def chi_so(y, diem, tau=TAU, moc_auc=MOC_AUC):
    """Harrell C, Uno C (cắt tại tau), AUC(t) tại các mốc. diem: điểm nguy cơ (lớn = nguy cơ cao).
    Phân bố kiểm duyệt cho IPCW ước lượng từ chính tập được đánh giá. Chỉ số không tính được -> NaN."""
    e, t = y[y.dtype.names[0]], y[y.dtype.names[1]]
    kq = {"harrell_c": concordance_index_censored(e, t, diem)[0]}
    try:
        kq["uno_c"] = concordance_index_ipcw(y, y, diem, tau=tau)[0]
    except Exception:
        kq["uno_c"] = np.nan
    for m in moc_auc:
        try:
            kq[f"auc_{m / NAM:g}n"] = cumulative_dynamic_auc(y, y, diem, [m])[0][0]
        except Exception:
            kq[f"auc_{m / NAM:g}n"] = np.nan
    return kq


# ---------------------------------------------------------------------------
# CoxNet chọn alpha theo quy tắc 1-SE
# ---------------------------------------------------------------------------
class CoxNet1SE(BaseEstimator):
    """CoxNet elastic-net, alpha chọn bằng CV lặp (phân tầng theo sự kiện) với quy tắc 1-SE:
    alpha lớn nhất (mô hình đơn giản nhất) có Harrell C trung bình ≥ C tốt nhất − 1 SE."""

    def __init__(self, l1_ratio=0.5, alpha_min_ratio=0.01, n_alphas=50, k=5, n_lap=3, seed=0):
        self.l1_ratio, self.alpha_min_ratio, self.n_alphas = l1_ratio, alpha_min_ratio, n_alphas
        self.k, self.n_lap, self.seed = k, n_lap, seed

    def _mo_hinh(self, **kw):
        return CoxnetSurvivalAnalysis(l1_ratio=self.l1_ratio, alpha_min_ratio=self.alpha_min_ratio, **kw)

    def fit(self, X, y):
        X = np.asarray(X, float)
        e = y[y.dtype.names[0]]
        duong = self._mo_hinh(n_alphas=self.n_alphas).fit(X, y).alphas_
        diem_c = []                                            # (lần chia) × alpha
        for lap in range(self.n_lap):
            cv = StratifiedKFold(self.k, shuffle=True, random_state=self.seed + lap)
            for tr, te in cv.split(X, e):
                hang = np.full(len(duong), np.nan)
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        m = self._mo_hinh(alphas=duong).fit(X[tr], y[tr])
                    et, tt = y[te][y.dtype.names[0]], y[te][y.dtype.names[1]]
                    if et.sum() > 0:
                        for i, a in enumerate(m.alphas_):
                            hang[i] = concordance_index_censored(et, tt, m.predict(X[te], alpha=a))[0]
                except (ArithmeticError, ValueError):
                    pass
                diem_c.append(hang)
        D = np.array(diem_c)
        tb = np.nanmean(D, axis=0)
        se = np.nanstd(D, axis=0, ddof=1) / np.sqrt(np.sum(~np.isnan(D), axis=0))
        i_tot = int(np.nanargmax(tb))
        hop_le = np.where(tb >= tb[i_tot] - se[i_tot])[0]
        self.alpha_ = float(duong[hop_le.min()])               # đường alpha giảm dần: chỉ số nhỏ = alpha lớn
        self.cv_c_ = pd.DataFrame({"alpha": duong, "c_tb": tb, "se": se})
        self.mo_hinh_ = self._mo_hinh(alphas=[self.alpha_], fit_baseline_model=True).fit(X, y)
        self.coef_ = self.mo_hinh_.coef_[:, 0]
        return self

    def predict(self, X):
        return self.mo_hinh_.predict(np.asarray(X, float))


def tao_coxnet_1se(**kw):
    """Pipeline chuẩn của ND4: chuẩn hóa z (khớp trên dữ liệu huấn luyện) + CoxNet 1-SE."""
    return make_pipeline(StandardScaler(), CoxNet1SE(**kw))


# ---------------------------------------------------------------------------
# Bootstrap hiệu chỉnh lạc quan (Harrell)
# ---------------------------------------------------------------------------
def bootstrap_lac_quan(pipeline, X, y, B=200, seed=2026, tau=TAU, moc_auc=MOC_AUC):
    """Mỗi lần: lấy mẫu có hoàn lại, khớp TOÀN BỘ pipeline (gồm chọn siêu tham số ở vòng trong),
    tính chỉ số trên chính mẫu bootstrap và trên dữ liệu gốc.
    optimism = trung bình (chỉ số_boot − chỉ số_gốc); hiệu chỉnh = biểu kiến − optimism.
    Trả về (bảng tóm tắt theo chỉ số, bảng từng lần bootstrap)."""
    X = np.asarray(X, float)
    rng = np.random.default_rng(seed)
    bieu_kien = chi_so(y, clone(pipeline).fit(X, y).predict(X), tau, moc_auc)
    lan = []
    for b in range(B):
        idx = rng.integers(0, len(y), len(y))
        if y[idx][y.dtype.names[0]].sum() < 2:
            continue
        try:
            m = clone(pipeline).fit(X[idx], y[idx])
        except Exception:
            continue
        cb = chi_so(y[idx], m.predict(X[idx]), tau, moc_auc)
        cg = chi_so(y, m.predict(X), tau, moc_auc)
        lan.append({"lan": b, **{f"{k}_boot": v for k, v in cb.items()}, **{f"{k}_goc": v for k, v in cg.items()}})
    lan = pd.DataFrame(lan)
    tom = []
    for k, v in bieu_kien.items():
        lq = np.nanmean(lan[f"{k}_boot"] - lan[f"{k}_goc"]) if len(lan) else np.nan
        tom.append({"chi_so": k, "bieu_kien": v, "lac_quan": lq, "hieu_chinh": v - lq,
                    "so_lan_hop_le": int((lan[f"{k}_boot"] - lan[f"{k}_goc"]).notna().sum()) if len(lan) else 0})
    return pd.DataFrame(tom), lan


# ---------------------------------------------------------------------------
# Repeated stratified CV
# ---------------------------------------------------------------------------
def cv_lap(pipeline, X, y, n_lap=20, k=5, seed=2026, tau=TAU, moc_auc=MOC_AUC):
    """Repeated stratified k-fold CV (phân tầng theo sự kiện). Mỗi lần lặp: gộp dự đoán ngoài fold của mọi fold
    rồi tính chỉ số một lần. Lưu ý: điểm nguy cơ của các fold khác nhau có thang hơi khác nhau, nên gộp có sai
    lệch nhẹ (thường làm C thấp đi). Trả về (tóm tắt: trung vị và khoảng 2,5–97,5%, bảng từng lần lặp)."""
    X = np.asarray(X, float)
    e = y[y.dtype.names[0]]
    lan = []
    for l in range(n_lap):
        du_doan = np.full(len(y), np.nan)
        for tr, te in StratifiedKFold(k, shuffle=True, random_state=seed + l).split(X, e):
            du_doan[te] = clone(pipeline).fit(X[tr], y[tr]).predict(X[te])
        lan.append({"lan": l, **chi_so(y, du_doan, tau, moc_auc)})
    lan = pd.DataFrame(lan)
    q = lan.drop(columns="lan").quantile([0.025, 0.5, 0.975]).T
    q.columns = ["q2_5", "trung_vi", "q97_5"]
    return q.rename_axis("chi_so").reset_index(), lan


# ---------------------------------------------------------------------------
# Đối chứng null
# ---------------------------------------------------------------------------
def _thap_phan_vi(trung_binh_bieu_hien, n_bin):
    return pd.qcut(trung_binh_bieu_hien.rank(method="first"), n_bin, labels=False)


def cap_ngau_nhien(cap, diem_mirbind2, trung_binh_bieu_hien, rng, n_bin=10):
    """Đối chứng null (ii). Giữ nguyên miRNA; thay gen bằng một gen KHÔNG phải đích
    (điểm miRBind2 dưới trung vị của miRNA đó), cùng thập phân vị biểu hiện trung bình, chưa dùng cho miRNA đó.

    cap: DataFrame cột mirna, gen. diem_mirbind2: DataFrame miRNA × gen (điểm). trung_binh_bieu_hien: Series theo gen.
    rng: np.random.Generator. Trả về DataFrame cùng dạng cap (cột gen đã thay)."""
    bin_gen = _thap_phan_vi(trung_binh_bieu_hien.loc[diem_mirbind2.columns], n_bin)
    moi = []
    for mir, g in cap[["mirna", "gen"]].itertuples(index=False):
        hang = diem_mirbind2.loc[mir]
        da_dung = set(cap.gen[cap.mirna == mir]) | {x for m, x in moi if m == mir}
        ung_vien = hang.index[(hang < hang.median()) & (bin_gen == bin_gen[g]) & ~hang.index.isin(da_dung)]
        if len(ung_vien) == 0:
            raise ValueError(f"Không còn gen không phải đích cùng thập phân vị cho {mir}|{g}")
        moi.append((mir, rng.choice(np.asarray(ung_vien))))
    return pd.DataFrame(moi, columns=["mirna", "gen"])


def p_thuc_nghiem(c_that, c_null):
    """p = (1 + #{C_null ≥ C_thật}) / (1 + số lần null)."""
    c_null = np.asarray(c_null, float)
    c_null = c_null[~np.isnan(c_null)]
    return (1 + np.sum(c_null >= c_that)) / (1 + len(c_null))


def chay_null_cap_ngau_nhien(ham_c, cap, diem_mirbind2, trung_binh_bieu_hien, c_that, n=200, seed=2026, n_bin=10):
    """Chạy lại toàn bộ pipeline `n` lần với cặp ngẫu nhiên.
    ham_c(cap) -> C hiệu chỉnh (vd. gọi bootstrap_lac_quan hoặc cv_lap trên đặc trưng tính từ cap).
    Trả về (p thực nghiệm, mảng C null)."""
    rng = np.random.default_rng(seed)
    c_null = [ham_c(cap_ngau_nhien(cap, diem_mirbind2, trung_binh_bieu_hien, rng, n_bin)) for _ in range(n)]
    return p_thuc_nghiem(c_that, c_null), np.asarray(c_null)


def dac_trung_rieng_le(mir_log2, mr_log2, cap):
    """Đối chứng (i): cùng miRNA và gen của các cặp, nhưng đưa vào thành đặc trưng riêng (mẫu × phân tử)."""
    mir = mir_log2.loc[pd.unique(cap.mirna)]
    mr = mr_log2.loc[pd.unique(cap.gen)]
    return pd.concat([mir.T.add_prefix("mir:"), mr.T.add_prefix("gen:")], axis=1)


def bo_doi_chung(dac_trung_cap, rieng_le, lam_sang_mien_dich):
    """Bốn bộ đặc trưng đối chứng bắt buộc (mẫu × đặc trưng, cùng thứ tự mẫu):
    cap (chính), (i) rieng_le, (iii) lam_sang_mien_dich (tuổi, giai đoạn, hạch, điểm miễn dịch, purity),
    (iv) lam_sang_mien_dich + cap (giá trị cộng thêm). Đối chứng (ii) dùng chay_null_cap_ngau_nhien."""
    assert list(dac_trung_cap.index) == list(rieng_le.index) == list(lam_sang_mien_dich.index), "Lệch thứ tự mẫu"
    return {"cap": dac_trung_cap, "i_rieng_le": rieng_le, "iii_lam_sang_mien_dich": lam_sang_mien_dich,
            "iv_lam_sang_mien_dich_cong_cap": pd.concat([lam_sang_mien_dich, dac_trung_cap], axis=1)}


# ---------------------------------------------------------------------------
# QC miRBind2
# ---------------------------------------------------------------------------
def kiem_tra_ks_mirbind2(rho_tung_phan, diem_mirbind2, trung_binh_bieu_hien, top=0.05, n_bin=10, k_ghep=5, seed=2026):
    """Với mỗi miRNA: so phân bố rho từng phần của đích top `top` (theo điểm miRBind2 của miRNA đó) với gen không
    phải đích (điểm < trung vị) được ghép cùng thập phân vị biểu hiện (k_ghep gen/đích, có hoàn lại).
    ks_2samp(x_dich, x_khong_dich, alternative="greater"): giả thuyết đối là CDF của đích nằm trên, tức đích lệch
    về phía âm. Trả về (bảng theo miRNA, tóm tắt: tỷ lệ miRNA có FDR < 0,05 và trung vị chênh lệch rho)."""
    rng = np.random.default_rng(seed)
    gen = rho_tung_phan.columns.intersection(diem_mirbind2.columns)
    bin_gen = _thap_phan_vi(trung_binh_bieu_hien.loc[gen], n_bin)
    dong = []
    for mir in rho_tung_phan.index.intersection(diem_mirbind2.index):
        d, r = diem_mirbind2.loc[mir, gen], rho_tung_phan.loc[mir, gen]
        dich = gen[d >= d.quantile(1 - top)]
        khong = gen[d < d.median()]
        ghep = []
        for b, so in bin_gen[dich].value_counts().items():
            nhom = khong[bin_gen[khong] == b]
            if len(nhom):
                ghep += list(rng.choice(np.asarray(nhom), so * k_ghep))
        if len(dich) < 5 or len(ghep) < 5:
            continue
        x_d, x_k = r[dich].values, r[ghep].values
        kq = ks_2samp(x_d, x_k, alternative="greater")
        dong.append({"mirna": mir, "so_dich": len(dich), "so_khong_dich_ghep": len(ghep),
                     "trung_vi_rho_dich": np.median(x_d), "trung_vi_rho_khong_dich": np.median(x_k),
                     "chenh_trung_vi": np.median(x_d) - np.median(x_k), "ks": kq.statistic, "p": kq.pvalue})
    bang = pd.DataFrame(dong)
    if len(bang):
        bang["fdr"] = multipletests(bang.p, method="fdr_bh")[1]
    tom = {"so_mirna": len(bang),
           "ty_le_fdr_05": float((bang.fdr < 0.05).mean()) if len(bang) else np.nan,
           "trung_vi_chenh_rho": float(bang.chenh_trung_vi.median()) if len(bang) else np.nan}
    return bang, tom
