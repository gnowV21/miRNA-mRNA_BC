"""Gắn cờ gen có exon chồng lấn với gen khác (bất kể mạch), thay cho pyranges (không build được trên Windows).

Lý do: STAR --quantMode GeneCounts (GDC "STAR - Counts") bỏ các read rơi vào exon của nhiều gen
(dòng N_ambiguous). Với thư viện unstranded như TCGA-BRCA, gen chồng lấn ở mạch nào cũng bị ảnh hưởng.
Ví dụ: CD68 và EIF4A1 (chr17) có số đếm thấp bất thường.
"""
import gzip
import numpy as np
import pandas as pd


def doc_exon_gtf(path):
    """Đọc exon từ GTF (.gtf hoặc .gtf.gz). Trả về DataFrame: chrom, start (0-based), end, gene_id, gene_name, gene_type."""
    mo = gzip.open if str(path).endswith(".gz") else open
    hang = []
    with mo(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t", 8)
            if p[2] != "exon":
                continue
            a = p[8]
            lay = lambda k: a.split(f'{k} "', 1)[1].split('"', 1)[0]
            hang.append((p[0], int(p[3]) - 1, int(p[4]), lay("gene_id"), lay("gene_name"), lay("gene_type")))
    return pd.DataFrame(hang, columns=["chrom", "start", "end", "gene_id", "gene_name", "gene_type"])


def _gop(df, khoa):
    """Gộp các đoạn chồng nhau trong cùng nhóm `khoa`. Đoạn nửa mở [start, end)."""
    df = df.sort_values([khoa, "start"]).reset_index(drop=True)
    cmax = df.groupby(khoa).end.cummax()
    truoc = cmax.groupby(df[khoa]).shift()
    khoi = (truoc.isna() | (df.start >= truoc)).cumsum()
    return (df.assign(_k=khoi.values).groupby("_k")
              .agg(**{khoa: (khoa, "first"), "chrom": ("chrom", "first"), "start": ("start", "min"), "end": ("end", "max")})
              .reset_index(drop=True))


def gen_chong_lan(exon):
    """exon: kết quả doc_exon_gtf. Trả về mỗi gen một dòng:
    gene_id, so_gen_chong (số gen khác có exon chồng), bp_chong (bp exon bị chồng, đã gộp, không đếm trùng),
    bp_exon (tổng bp exon đã gộp), ty_le_exon_chong."""
    khoi = _gop(exon[["chrom", "start", "end", "gene_id"]], "gene_id")      # exon đã gộp trong từng gen
    doan, cap = [], []
    for _, g in khoi.groupby("chrom"):
        g = g.sort_values("start")
        s, e, gid = g.start.values, g.end.values, g.gene_id.values
        j_max = np.searchsorted(s, e, side="left")                            # các khối có start < end[i]
        for i in range(len(s)):
            for j in range(i + 1, j_max[i]):
                if gid[j] != gid[i]:
                    a, b = s[j], min(e[i], e[j])
                    if b > a:
                        doan += [(gid[i], a, b), (gid[j], a, b)]
                        cap += [(gid[i], gid[j]), (gid[j], gid[i])]
    tong = khoi.assign(bp=khoi.end - khoi.start).groupby("gene_id").bp.sum().rename("bp_exon")
    kq = pd.DataFrame(index=tong.index)
    kq["so_gen_chong"] = 0
    kq["bp_chong"] = 0
    if doan:
        d = pd.DataFrame(doan, columns=["gene_id", "start", "end"]).assign(chrom="")
        bp = _gop(d, "gene_id").assign(bp=lambda x: x.end - x.start).groupby("gene_id").bp.sum()
        so = pd.DataFrame(cap, columns=["gene_id", "gene_b"]).groupby("gene_id").gene_b.nunique()
        kq.loc[bp.index, "bp_chong"] = bp
        kq.loc[so.index, "so_gen_chong"] = so
    kq["bp_exon"] = tong
    kq["ty_le_exon_chong"] = (kq.bp_chong / kq.bp_exon).round(4)
    return kq.reset_index()
