import numpy as np, pandas as pd
from mirna_bc.dac_trung import (f_raw, f_rank, f_ref, spearman_tung_phan, chon_tham_chieu,
                                diem_hoat_dong_mirna, delta_pcc)

rng = np.random.default_rng(1)
MAU = [f"s{i}" for i in range(40)]
MIR = pd.DataFrame(rng.normal(8, 2, (30, 40)), index=[f"miR{i}" for i in range(30)], columns=MAU)
MR = pd.DataFrame(rng.normal(5, 2, (80, 40)), index=[f"G{i}" for i in range(80)], columns=MAU)
CAP = pd.DataFrame({"mirna": ["miR0", "miR1", "miR2"], "gen": ["G0", "G5", "G9"]})


def test_f_rank_bat_bien_don_dieu():
    a = f_rank(MIR, MR, CAP)
    pd.testing.assert_frame_equal(a, f_rank(MIR, np.exp2(MR), CAP))
    pd.testing.assert_frame_equal(a, f_rank(MIR * 3 + 1, MR ** 3, CAP))


def test_f_raw_dich_hang_so():
    a = f_raw(MIR, MR, CAP)
    MR2 = MR.copy(); MR2.loc["G0"] += 2.5
    b = f_raw(MIR, MR2, CAP)
    assert np.allclose(a.loc["miR0|G0"] - b.loc["miR0|G0"], 2.5)
    pd.testing.assert_frame_equal(a.drop("miR0|G0"), b.drop("miR0|G0"))


def test_f_ref_bat_bien_lech_theo_mau():
    lech = pd.Series(rng.normal(0, 3, 40), index=MAU)        # lệch cộng theo mẫu (vd. lượng RNA đầu vào)
    a = f_ref(MIR, MR, CAP, ["miR10", "miR11"], ["G20", "G21"])
    b = f_ref(MIR + lech, MR - lech, CAP, ["miR10", "miR11"], ["G20", "G21"])
    pd.testing.assert_frame_equal(a, b)


def test_spearman_tung_phan_loai_yeu_to_an():
    n = 300
    an = rng.normal(size=n)                                   # yếu tố ẩn (vd. thâm nhiễm miễn dịch)
    X = pd.DataFrame({"mir": an + 0.3 * rng.normal(size=n)})
    Y = pd.DataFrame({"gen": -an + 0.3 * rng.normal(size=n)})
    rho_tho, _, _ = spearman_tung_phan(X, Y, np.empty((n, 0)))
    rho_tp, _, _ = spearman_tung_phan(X, Y, an[:, None])
    assert rho_tho[0, 0] < -0.8
    assert abs(rho_tp[0, 0]) < 0.15


def test_chon_tham_chieu():
    L = MIR.copy(); L.loc["miR29"] = L.mean(axis=0) + 5     # bám sát trung bình mẫu, biểu hiện cao
    bang = chon_tham_chieu(L, so=3)
    assert bang.index[0] == "miR29" and len(bang) == 3


def _tuong_quan(n_ref, n_mau, so_cap, r=0.5, seed=2):
    """Dữ liệu giả: so_cap cặp, mỗi cặp có miRNA và gen tương quan r (thang log2)."""
    g = np.random.default_rng(seed)
    C = np.array([[1, r], [r, 1]])
    def sinh(n, ten):
        z = g.multivariate_normal([0, 0], C, size=(so_cap, n))     # cặp × n × 2
        return (pd.DataFrame(8 + 2 * z[..., 0], index=[f"miR{k}" for k in range(so_cap)], columns=[f"{ten}{i}" for i in range(n)]),
                pd.DataFrame(5 + 1.5 * z[..., 1], index=[f"G{k}" for k in range(so_cap)], columns=[f"{ten}{i}" for i in range(n)]))
    mir_ref, mr_ref = sinh(n_ref, "lanh")
    mir_mau, mr_mau = sinh(n_mau, "u")
    cap = pd.DataFrame({"mirna": [f"miR{k}" for k in range(so_cap)], "gen": [f"G{k}" for k in range(so_cap)]})
    return mir_mau, mr_mau, mir_ref, mr_ref, cap


def test_delta_pcc_khop_tinh_truc_tiep():
    mir_mau, mr_mau, mir_ref, mr_ref, cap = _tuong_quan(60, 7, 4)
    D = delta_pcc(mir_mau, mr_mau, mir_ref, mr_ref, cap)
    assert D.shape == (4, 7) and list(D.columns) == list(mir_mau.columns)
    for a, b in zip(cap.mirna, cap.gen):
        x, y = mir_ref.loc[a].to_numpy(), mr_ref.loc[b].to_numpy()
        r0 = np.corrcoef(x, y)[0, 1]
        for s in mir_mau.columns:
            r1 = np.corrcoef(np.append(x, mir_mau.loc[a, s]), np.append(y, mr_mau.loc[b, s]))[0, 1]
            assert abs(D.loc[f"{a}|{b}", s] - (r1 - r0)) < 1e-10
    pd.testing.assert_frame_equal(delta_pcc(mir_mau, mr_mau, mir_ref, mr_ref, cap, tuyet_doi=True), D.abs())


def test_delta_pcc_xap_xi_ham_anh_huong():
    n = 500
    mir_mau, mr_mau, mir_ref, mr_ref, cap = _tuong_quan(n, 200, 3, r=0.4, seed=3)
    D = delta_pcc(mir_mau, mr_mau, mir_ref, mr_ref, cap).to_numpy()
    X, Y = mir_ref.loc[cap.mirna].to_numpy(), mr_ref.loc[cap.gen].to_numpy()
    r = np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(X, Y)])[:, None]
    zx = (mir_mau.loc[cap.mirna].to_numpy() - X.mean(1, keepdims=True)) / X.std(1, keepdims=True)
    zy = (mr_mau.loc[cap.gen].to_numpy() - Y.mean(1, keepdims=True)) / Y.std(1, keepdims=True)
    xap_xi = (zx * zy - r * (zx**2 + zy**2) / 2) / (n + 1)
    assert np.linalg.norm(D - xap_xi) / np.linalg.norm(D) < 0.05


def test_diem_hoat_dong_mirna():
    so_dich = {"miR0": 6, "miR1": 5, "miR2": 4}                # miR2 chỉ có 4 đích: bị loại
    cap = pd.DataFrame([(m, f"G{10 * k + j}") for k, (m, d) in enumerate(so_dich.items()) for j in range(d)],
                       columns=["mirna", "gen"])
    A = diem_hoat_dong_mirna(MIR, MR, cap)
    assert sorted(A.index) == ["miR0", "miR1"] and list(A.columns) == MAU
    F = f_rank(MIR, MR, cap)
    assert np.allclose(A.loc["miR1"], F.loc[[f"miR1|G{10 + j}" for j in range(5)]].mean())
    pd.testing.assert_frame_equal(A, diem_hoat_dong_mirna(MIR, np.exp2(MR), cap))
