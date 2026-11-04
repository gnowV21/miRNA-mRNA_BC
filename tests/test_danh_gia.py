"""Test khung đánh giá trên dữ liệu giả (không dùng dữ liệu bệnh nhân)."""
import numpy as np, pandas as pd
from mirna_bc.danh_gia import (tao_y, chi_so, tao_coxnet_1se, bootstrap_lac_quan, cv_lap, cap_ngau_nhien,
                               p_thuc_nghiem, kiem_tra_ks_mirbind2, dac_trung_rieng_le, bo_doi_chung, NAM)


def _song_con_gia(n=150, p=6, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    nguy_co = 0.9 * X[:, 0]                                   # chỉ đặc trưng 0 mang tín hiệu
    t_su_kien = rng.exponential(3 * NAM * np.exp(-nguy_co))
    t_kiem_duyet = rng.uniform(0.5 * NAM, 6 * NAM, n)
    return X, tao_y(t_su_kien <= t_kiem_duyet, np.minimum(t_su_kien, t_kiem_duyet))


def test_chi_so_xep_hang_dung():
    X, y = _song_con_gia()
    tot, ngau_nhien = chi_so(y, X[:, 0]), chi_so(y, np.random.default_rng(1).normal(size=len(y)))
    assert tot["harrell_c"] > 0.65 and abs(ngau_nhien["harrell_c"] - 0.5) < 0.1
    assert set(tot) == {"harrell_c", "uno_c", "auc_2n", "auc_3n"}


def test_coxnet_1se_chon_dac_trung_co_tin_hieu():
    X, y = _song_con_gia()
    m = tao_coxnet_1se(k=3, n_lap=1).fit(X, y)
    he_so = m[-1].coef_
    assert abs(he_so[0]) == np.max(np.abs(he_so)) and he_so[0] > 0


def test_bootstrap_lac_quan():
    X, y = _song_con_gia()
    tom, lan = bootstrap_lac_quan(tao_coxnet_1se(k=3, n_lap=1), X, y, B=8, seed=3)
    tom = tom.set_index("chi_so")
    assert len(lan) >= 6
    assert tom.loc["harrell_c", "hieu_chinh"] <= tom.loc["harrell_c", "bieu_kien"] + 0.02
    assert tom.loc["harrell_c", "hieu_chinh"] > 0.6


def test_cv_lap():
    X, y = _song_con_gia()
    tom, lan = cv_lap(tao_coxnet_1se(k=3, n_lap=1), X, y, n_lap=2, k=3)
    assert len(lan) == 2
    r = tom.set_index("chi_so").loc["harrell_c"]
    assert r.q2_5 <= r.trung_vi <= r.q97_5 and r.trung_vi > 0.6


def _mirbind_gia(seed=0):
    rng = np.random.default_rng(seed)
    mir, gen = [f"miR{i}" for i in range(4)], [f"G{j}" for j in range(400)]
    diem = pd.DataFrame(rng.uniform(size=(4, 400)), index=mir, columns=gen)
    tb = pd.Series(rng.normal(5, 2, 400), index=gen)
    return diem, tb


def test_cap_ngau_nhien_giu_mirna_doi_gen_khong_dich_cung_bin():
    diem, tb = _mirbind_gia()
    cap = pd.DataFrame({"mirna": ["miR0", "miR1"], "gen": [diem.loc["miR0"].idxmax(), diem.loc["miR1"].idxmax()]})
    moi = cap_ngau_nhien(cap, diem, tb, np.random.default_rng(0))
    assert list(moi.mirna) == list(cap.mirna)
    bin_ = pd.qcut(tb.rank(method="first"), 10, labels=False)
    for (m, g0), g1 in zip(cap.itertuples(index=False), moi.gen):
        assert diem.loc[m, g1] < diem.loc[m].median() and bin_[g1] == bin_[g0] and g1 != g0
    pd.testing.assert_frame_equal(moi, cap_ngau_nhien(cap, diem, tb, np.random.default_rng(0)))


def test_p_thuc_nghiem():
    assert p_thuc_nghiem(0.7, [0.5] * 199 + [0.8]) == 2 / 201
    assert p_thuc_nghiem(0.5, [0.6] * 200) == 1.0


def test_ks_mirbind2_phat_hien_dich_lech_am():
    diem, tb = _mirbind_gia()
    rng = np.random.default_rng(5)
    rho = pd.DataFrame(rng.normal(0, 0.1, diem.shape), index=diem.index, columns=diem.columns)
    co_tin_hieu = rho.copy()
    for m in ["miR0", "miR1", "miR2"]:
        dich = diem.columns[diem.loc[m] >= diem.loc[m].quantile(0.9)]
        co_tin_hieu.loc[m, dich] -= 0.15
    bang, tom = kiem_tra_ks_mirbind2(co_tin_hieu, diem, tb, top=0.1)
    assert bang.set_index("mirna").loc[["miR0", "miR1", "miR2"], "fdr"].lt(0.05).all()
    assert tom["trung_vi_chenh_rho"] < 0
    _, tom0 = kiem_tra_ks_mirbind2(rho, diem, tb, top=0.1)
    assert tom0["ty_le_fdr_05"] <= 0.25


def test_dac_trung_rieng_le_va_bo_doi_chung():
    mau = [f"s{i}" for i in range(5)]
    mir = pd.DataFrame(np.arange(10).reshape(2, 5), index=["miR0", "miR1"], columns=mau)
    mr = pd.DataFrame(np.arange(15).reshape(3, 5), index=["G0", "G1", "G2"], columns=mau)
    cap = pd.DataFrame({"mirna": ["miR0", "miR0"], "gen": ["G0", "G2"]})
    rl = dac_trung_rieng_le(mir, mr, cap)
    assert list(rl.columns) == ["mir:miR0", "gen:G0", "gen:G2"] and list(rl.index) == mau
    ls = pd.DataFrame({"tuoi": range(5)}, index=mau)
    bo = bo_doi_chung(rl.iloc[:, :1], rl, ls)
    assert bo["iv_lam_sang_mien_dich_cong_cap"].shape == (5, 2)
