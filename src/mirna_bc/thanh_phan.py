"""Điểm thành phần mô tính trong từng mẫu (single-sample), dùng chung cho TCGA (TPM) và METABRIC (microarray log2).

Điểm = trung bình phân vị (trong cùng mẫu) của các gen chỉ dấu, nên không phụ thuộc thang đo:
mọi biến đổi đơn điệu tăng áp lên cả mẫu (nhân hằng số, log, exp) cho cùng một điểm.
Lưu ý: phân vị phụ thuộc tập gen nền. Muốn so sánh giữa hai cohort thì đưa vào cùng một tập gen.
"""
import pandas as pd

BO_GEN = {
    "mien_dich": ["PTPRC", "CD2", "CD3D", "CD3E", "CD247", "CD48", "CD52", "CORO1A", "LCP2",
                  "IKZF1", "CCL5", "GZMA", "GZMK", "CD8A", "MS4A1", "CD79A", "CXCL9", "CXCL10"],
    "mo_dem":    ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A3", "DCN", "LUM", "FAP",
                  "PDGFRB", "SPARC", "THY1", "FBN1", "VCAN", "MMP2"],
    "bieu_mo":   ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "CLDN4", "CLDN7", "ELF3", "SPINT2"],
    "LAR":       ["AR", "FOXA1", "PIP", "AGR2", "TFF1", "SPDEF", "KRT18"],
    "dai_thuc_bao": ["CD68", "CD163", "CSF1R", "C1QA", "C1QB", "MS4A4A"],
}


def diem_thanh_phan(bieu_hien: pd.DataFrame, bo_gen=BO_GEN, loai_gen=()) -> pd.DataFrame:
    """bieu_hien: gen (Hugo symbol) × mẫu, thang bất kỳ đơn điệu. Trả về mẫu × điểm (0–1).

    loai_gen: gen không dùng khi tính điểm (vd. gen bị định lượng sai do chồng lấn trong STAR).
    Cột "<tên>_so_gen" ghi số gen thực sự dùng.
    """
    pv = bieu_hien.rank(axis=0, pct=True)                  # phân vị trong từng mẫu
    loai = set(loai_gen)
    kq = {}
    for ten, ds in bo_gen.items():
        co = [g for g in ds if g in pv.index and g not in loai]
        kq[ten] = pv.loc[co].mean(axis=0)
        kq[f"{ten}_so_gen"] = len(co)
    return pd.DataFrame(kq)
