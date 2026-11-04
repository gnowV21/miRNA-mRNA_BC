import gzip
from mirna_bc.chong_lan import doc_exon_gtf, gen_chong_lan

GTF = """##test
chr1\tT\texon\t101\t200\t.\t+\t.\tgene_id "A"; gene_name "A"; gene_type "protein_coding";
chr1\tT\texon\t151\t300\t.\t+\t.\tgene_id "A"; gene_name "A"; gene_type "protein_coding";
chr1\tT\texon\t251\t400\t.\t-\t.\tgene_id "B"; gene_name "B"; gene_type "lncRNA";
chr1\tT\texon\t1001\t1100\t.\t+\t.\tgene_id "C"; gene_name "C"; gene_type "protein_coding";
chr1\tT\texon\t1101\t1200\t.\t+\t.\tgene_id "D"; gene_name "D"; gene_type "protein_coding";
chr2\tT\texon\t251\t400\t.\t+\t.\tgene_id "E"; gene_name "E"; gene_type "protein_coding";
"""


def test_chong_lan(tmp_path):
    f = tmp_path / "t.gtf.gz"
    with gzip.open(f, "wt") as g:
        g.write(GTF)
    kq = gen_chong_lan(doc_exon_gtf(f)).set_index("gene_id")
    # A: exon gộp [100, 300) dài 200 bp; chồng với B (khác mạch) ở [250, 300) = 50 bp
    assert kq.loc["A", "bp_exon"] == 200 and kq.loc["A", "bp_chong"] == 50 and kq.loc["A", "so_gen_chong"] == 1
    assert kq.loc["B", "bp_chong"] == 50
    # C và D chỉ kề nhau (không chồng); E ở nhiễm sắc thể khác
    assert (kq.loc[["C", "D", "E"], "so_gen_chong"] == 0).all()
