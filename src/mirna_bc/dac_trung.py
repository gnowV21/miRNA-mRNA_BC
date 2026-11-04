"""Đặc trưng cặp miRNA–mRNA và tương quan từng phần (ND2). Không dùng dữ liệu sống còn.

Ba dạng đặc trưng (quyết định D5, 11/2026):
- f_rank (chính): hiệu phân vị trong mẫu, bất biến với mọi biến đổi đơn điệu theo từng nền tảng.
- f_raw (độ nhạy): log2(miRNA+1) − log2(mRNA+1), đầu vào đã ở thang log2.
- f_ref (độ nhạy, dạng qPCR): log2 trừ trung bình log2 của phân tử tham chiếu trong cùng mẫu (giống ΔCt).
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
