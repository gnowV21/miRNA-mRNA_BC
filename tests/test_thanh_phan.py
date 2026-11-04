import numpy as np, pandas as pd
from mirna_bc.thanh_phan import diem_thanh_phan, BO_GEN


def _du_lieu(seed=0):
    rng = np.random.default_rng(seed)
    gen = sorted({g for ds in BO_GEN.values() for g in ds}) + [f"G{i}" for i in range(300)]
    return pd.DataFrame(rng.gamma(2, 50, (len(gen), 12)), index=gen, columns=[f"m{j}" for j in range(12)])


def test_bat_bien_nhan_hang_so_va_log():
    X = _du_lieu()
    d0 = diem_thanh_phan(X)
    pd.testing.assert_frame_equal(d0, diem_thanh_phan(X * 7.3))
    pd.testing.assert_frame_equal(d0, diem_thanh_phan(np.log2(X + 1)))
    pd.testing.assert_frame_equal(d0, diem_thanh_phan(np.exp2(np.log2(X + 1))))


def test_loai_gen_va_so_gen():
    X = _du_lieu()
    d = diem_thanh_phan(X, loai_gen=["CD68"])
    assert d["dai_thuc_bao_so_gen"].iloc[0] == len(BO_GEN["dai_thuc_bao"]) - 1
    assert d[["mien_dich", "mo_dem", "bieu_mo"]].stack().between(0, 1).all()
