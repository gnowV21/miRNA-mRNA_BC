import numpy as np, pandas as pd
from mirna_bc.dac_trung import f_raw, f_rank, f_ref, spearman_tung_phan, chon_tham_chieu

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
