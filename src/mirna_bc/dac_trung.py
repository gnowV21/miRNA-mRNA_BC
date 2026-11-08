"""Đặc trưng cặp miRNA–mRNA và tương quan từng phần (ND2). Không dùng dữ liệu sống còn.

Ba dạng đặc trưng (quyết định D5, 11/2026):
- f_rank (chính): hiệu phân vị trong mẫu, bất biến với mọi biến đổi đơn điệu theo từng nền tảng.
- f_raw (độ nhạy): log2(miRNA+1) − log2(mRNA+1), đầu vào đã ở thang log2.
- f_ref (độ nhạy, dạng qPCR): log2 trừ trung bình log2 của phân tử tham chiếu trong cùng mẫu (giống ΔCt).
Bộ đặc trưng phụ: điểm hoạt động miRNA (D10). Đối chứng thăm dò: ΔPCC so với mô lành (D11, Lee et al., 2023).
"""
import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as phan_bo_t
from statsmodels.stats.multitest import multipletests


def f_raw(mir_log2, mr_log2, cap):
    """log2(miRNA+1) − log2(mRNA+1), đầu vào đã ở thang log2. cap: DataFrame cột mirna, gen.
    Trả về cặp ("mirna|gen") × mẫu."""
    return pd.DataFrame({f"{a}|{b}": mir_log2.loc[a] - mr_log2.loc[b]
                         for a, b in zip(cap.mirna, cap.gen)}).T


def chuan_hoa_tham_chieu(L, tham_chieu):
    """Trừ trung bình log2 của các phân tử tham chiếu trong cùng mẫu (giống ΔCt)."""
    return L - L.loc[list(tham_chieu)].mean(axis=0)


def f_ref(mir_log2, mr_log2, cap, ref_mir, ref_gen):
    return f_raw(chuan_hoa_tham_chieu(mir_log2, ref_mir), chuan_hoa_tham_chieu(mr_log2, ref_gen), cap)


def f_rank(mir, mr, cap):
    """Hiệu phân vị trong mẫu. Bất biến với mọi biến đổi đơn điệu theo từng nền tảng.
    BẮT BUỘC: mir và mr có cùng tập miRNA/gen ở mọi cohort (giao nền tảng)."""
    return f_raw(mir.rank(axis=0, pct=True), mr.rank(axis=0, pct=True), cap)


def diem_hoat_dong_mirna(mir, mr, cap, toi_thieu=5):
    """D10 (08/11/2026). Điểm hoạt động miRNA = trung bình f_rank trên các đích đã qua bộ lọc ND2 của miRNA đó.

    mir, mr: phân tử × mẫu (cùng tập phân tử ở mọi cohort, như f_rank). cap: DataFrame cột mirna, gen
    (cặp đã qua bộ lọc ND2). Chỉ tính cho miRNA có ≥ `toi_thieu` đích. Trọng số cố định (trung bình),
    không học từ kết cục. Trả về miRNA × mẫu; bất biến với biến đổi đơn điệu như f_rank.
    """
    F = f_rank(mir, mr, cap)                                  # cặp × mẫu, index "mirna|gen"
    nhom = F.index.str.split("|").str[0]
    dem = nhom.value_counts()
    giu = nhom.isin(dem.index[dem >= toi_thieu])
    return F[giu].groupby(nhom[giu]).mean()


def delta_pcc(mir_mau, mr_mau, mir_ref, mr_ref, cap, tuyet_doi=False):
    """D11 (08/11/2026) – ΔPCC theo Lee et al. (2023), đối chứng thăm dò.

    *_ref: phân tử × n mẫu mô lành; *_mau: phân tử × mẫu u. cap: DataFrame cột mirna, gen.
    Trả về cặp ("mirna|gen") × mẫu u: PCC(n mẫu lành + 1 mẫu u) − PCC(n mẫu lành), tính bằng công thức
    cập nhật tổng (không lặp qua từng mẫu).

    Quy ước của D11:
    - Đầu vào ở **thang log2** (cùng đầu vào với f_raw: log2(RPM + 1), log2(TPM + 1)), khác Lee 2023 vốn dùng
      CPM thô; lý do: Pearson trên CPM thô bị vài giá trị cực đoan chi phối.
    - Biến thể chính là **ΔPCC có dấu**; |ΔPCC| như Lee 2023 là biến thể phụ (`tuyet_doi=True`).
    - Dùng **cùng tập cặp** đã qua bộ lọc ND2 với f_rank, để so sánh công bằng.
    Xấp xỉ bậc một (hàm ảnh hưởng của hệ số tương quan): ΔPCC ≈ [z_x·z_y − r(z_x² + z_y²)/2] / (n + 1),
    với z chuẩn hóa theo trung bình và độ lệch chuẩn của mẫu lành, r = PCC trên mẫu lành.
    """
    X = mir_ref.loc[cap.mirna].to_numpy(float)                # cặp × n
    Y = mr_ref.loc[cap.gen].to_numpy(float)
    x = mir_mau.loc[cap.mirna].to_numpy(float)                # cặp × mẫu u
    y = mr_mau.loc[cap.gen].to_numpy(float)
    # Trừ trung bình mẫu lành của từng cặp cho cả hai phía (PCC không đổi) để giảm sai số làm tròn
    mx, my = X.mean(1, keepdims=True), Y.mean(1, keepdims=True)
    X, Y, x, y = X - mx, Y - my, x - mx, y - my
    n = X.shape[1]
    Sx, Sy = X.sum(1), Y.sum(1)
    Sxx, Syy, Sxy = (X * X).sum(1), (Y * Y).sum(1), (X * Y).sum(1)

    def pcc(m, sx, sy, sxx, syy, sxy):
        return (m * sxy - sx * sy) / np.sqrt((m * sxx - sx**2) * (m * syy - sy**2))

    r0 = pcc(n, Sx, Sy, Sxx, Syy, Sxy)[:, None]
    r1 = pcc(n + 1, Sx[:, None] + x, Sy[:, None] + y, Sxx[:, None] + x**2,
             Syy[:, None] + y**2, Sxy[:, None] + x * y)
    d = r1 - r0
    return pd.DataFrame(np.abs(d) if tuyet_doi else d,
                        index=[f"{a}|{b}" for a, b in zip(cap.mirna, cap.gen)], columns=mir_mau.columns)


def spearman_tung_phan(X, Y, Z):
    """X: mẫu × p, Y: mẫu × q, Z: mẫu × k hiệp biến (không có NaN; k = 0 cho tương quan thô).
    Trả về rho (p × q), p-value, FDR (BH trên toàn ma trận)."""
    Z = np.asarray(Z, float).reshape(len(X), -1)
    n, k = Z.shape
    R = lambda A: np.apply_along_axis(rankdata, 0, np.asarray(A, float))
    Zr = np.column_stack([np.ones(n), R(Z)]) if k else np.ones((n, 1))
    def phan_du(A):
        A = R(A); E = A - Zr @ np.linalg.lstsq(Zr, A, rcond=None)[0]
        return (E - E.mean(0)) / E.std(0, ddof=1)
    rho = phan_du(X).T @ phan_du(Y) / (n - 1)
    df = n - 2 - k
    tstat = rho * np.sqrt(df / np.clip(1 - rho**2, 1e-12, None))
    p = 2 * phan_bo_t.sf(np.abs(tstat), df)
    fdr = multipletests(p.ravel(), method="fdr_bh")[1].reshape(p.shape)
    return rho, p, fdr


def chon_tham_chieu(L, so=4, Z=None, nguong_tuong_quan=0.2, top_bieu_hien=0.3):
    """Chọn phân tử tham chiếu không dùng kết cục.

    L: phân tử × mẫu, thang log2. Z: mẫu × điểm thành phần (tùy chọn).
    Tiêu chí: biểu hiện cao (top 30% theo trung bình), |rho Spearman| với mọi điểm thành phần < ngưỡng,
    rồi lấy `so` phân tử ổn định nhất (SD nhỏ nhất của L − trung bình mẫu).
    Trả về DataFrame (đã sắp) gồm trung_binh, sd_tuong_doi và rho_max_thanh_phan (nếu có Z).
    """
    tb = L.mean(axis=1)
    ung_vien = L.loc[tb >= tb.quantile(1 - top_bieu_hien)]
    sd = (ung_vien - L.mean(axis=0)).std(axis=1)
    bang = pd.DataFrame({"trung_binh": tb.loc[ung_vien.index], "sd_tuong_doi": sd})
    if Z is not None:
        Zr = Z.loc[L.columns].rank()
        bang["rho_max_thanh_phan"] = ung_vien.T.rank().apply(lambda c: Zr.corrwith(c).abs().max())
        bang = bang[bang.rho_max_thanh_phan < nguong_tuong_quan]
    return bang.sort_values("sd_tuong_doi").head(so)
