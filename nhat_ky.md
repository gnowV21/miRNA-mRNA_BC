# Nhật ký dự án miRNA-mRNA_BC

Ghi mỗi quyết định và con số kèm ngày. Reviewer và hội đồng sẽ hỏi lý do của từng lựa chọn.

## Quyết định thiết kế đã chốt (10/2026)

Các dòng ghi "→ D…" đã được thay ở mốc kiểm tra 11/2026 (xem mục bên dưới).

- Đối tượng: TNBC theo IHC (ER−, PR−, HER2−); HER2 Equivocal chỉ giữ khi FISH âm. → D6 (thêm ca chỉ có ISH âm), D7 (cờ ER-low/LAR)
- Loại: bệnh nhân nam; đã điều trị tân bổ trợ; di căn xa lúc chẩn đoán (M1).
- Huấn luyện: TCGA-BRCA. Kiểm chứng ngoài: METABRIC. → D2 (nhánh mRNA kiểm chứng hai chiều)
- Endpoint chính: PFI (TCGA-CDR). Phụ: OS ở cả hai cohort; DFI khi cohort kiểm chứng chỉ có RFS/DFS. → D3
- METABRIC: định nghĩa IHC chặt; nếu thiếu IHC cho PR/HER2, chạy song song định nghĩa tổ hợp. (Đã xác minh: thiếu, chỉ dùng tổ hợp.)
- Mỗi bệnh nhân 1 mẫu: ưu tiên lọ A, rồi số đọc miRNA cao hơn.
- miRNA mức trưởng thành (MIMAT), RPM. mRNA: gen mã hóa protein, TPM. VST chỉ cho phân tích độ nhạy. → D8 (ngưỡng theo read)
- Nguồn cặp miRNA–mRNA: miRBind2, lọc tương quan âm trong từng vòng CV. → D4, D5
- So sánh: đặc trưng cặp vs BulkRNABert (bản không TCGA; thay thế scFoundation) vs BulkFormer vs biểu hiện gen thô.
- Mô hình: CoxNet (chính), RSF, Gradient Boosting Survival; nested CV; stability selection; quy tắc 1-SE. → D1
- Ngôn ngữ: Python; R chỉ cho immunedeconv và oncoPredict (Nội dung 5).
- Tài nguyên: Colab, Kaggle Notebooks; thuê GPU khi cần.

## Nội dung 1 – Thu thập và tiền xử lý

Số liệu chính thức là số trong `qc/` do notebook ND1 sinh ra (chạy 04/11/2026, GDC Data Release 46.0).

| Ngày | Bước | Nội dung ghi nhận | Con số |
| --- | --- | --- | --- |
| 04/11/2026 | Bước 0 | Phiên bản Python, pandas, pydeseq2, lifelines | Python 3.12.8; pandas 2.3.3; pydeseq2 0.5.4; lifelines 0.30.3; scipy 1.18.1; statsmodels 0.14.5. Đủ phiên bản ở `qc/pip_freeze.txt`; GDC ở `qc/phien_ban_gdc.json`; md5 mọi file: `qc/manifest_*.csv` (lệch 0) |
| 04/11/2026 | Bước 1 | TNBC sau từng tiêu chí loại; số ca thiếu IHC | 1097 → 972 (có ER, PR theo IHC và HER2 theo IHC hoặc ISH; loại 125 ca thiếu) → 202 (ER−, PR−) → 157 (HER2 âm; loại 37 IHC dương, 4 ISH dương, 3 mâu thuẫn, 5 Equivocal không có ISH âm) → 155 (loại M1). 17 ca HER2 chỉ có ISH. Ca mâu thuẫn: TCGA-AN-A0FL, TCGA-AN-A0FX (IHC 1+ nhưng ghi Positive), TCGA-LL-A73Y (IHC 3+ nhưng ghi Negative) |
| 04/11/2026 | Bước 2 | Số ca có PFI hợp lệ | 154. Loại sự kiện PFI (cohort cuối): 13 di căn xa, 5 u nguyên phát mới, 4 tái phát tại chỗ, 3 tử vong do u, 1 không rõ |
| 04/11/2026 | Bước 2b | Chú thích sinh học | PAM50 đủ cho cohort cuối (127 Basal, 12 Her2, 10 Normal, 4 LumA, 1 LumB); cBioPortal và TCGASubtype khớp 100% trên 139 ca. Purity ABSOLUTE và leukocyte fraction (Thorsson 2018) có cho phần lớn ca |
| 04/11/2026 | Bước 3 | Số miRNA trưởng thành; kiểm tra miR-21-5p | 2013 MIMAT (loại 56 MIMAT không có trong miRBase v22). miR-21-5p chiếm 24,8% tổng RPM (trung bình), top 10 chiếm 67%. **Cần tra miRBase:** MIMAT0017985 (trung vị 112 read, có ở 97% mẫu) và MIMAT0018073 (trung vị 13 read) bị loại vì không có trong v22 |
| 04/11/2026 | Bước 4 | Số gen mã hóa protein | 19.944. N_ambiguous trung vị 8,8% read mỗi mẫu. 9.154 gen có exon chồng lấn (GENCODE v36), 1.230 gen ≥50% exon bị chồng; CD68 (85,6%) và EIF4A1 được gắn cờ |
| 04/11/2026 | Bước 5 | Số bệnh nhân cuối có cả miRNA và mRNA | 152 (trước đây 136) |
| 04/11/2026 | Bước 5c | Đối chiếu STAR – RSEM (TOIL) | Chưa chạy: file TOIL 706 MB tải rất chậm (S3 ~100 KB/s, liên tục ngắt). Ghi kết quả khi có |
| 04/11/2026 | Bước 6–7 | 4 ngưỡng lọc + lý do; số miRNA, gen giữ lại | miRNA ≥10 read ở ≥50% mẫu; mRNA count ≥10 và TPM ≥1 ở ≥50% mẫu (D8). Lý do: RPM = 1 chỉ tương ứng 0,79–2,33 read. Giữ 333 miRNA; 13.179 gen qua ngưỡng, 12.825 gen sau giao với METABRIC |
| 04/11/2026 | Bước 7b | Chuẩn hóa miRNA thay thế | Spearman trung vị theo miRNA so với log2RPM gốc: bỏ miR-21 0,987; upper-quartile 0,944; DESeq2 0,946 |
| 04/11/2026 | Bước 8 | Sự kiện PFI / OS / DFI; theo dõi trung vị | PFI 26/152; OS 23/152; DFI 19/140; DSS 14/147. Theo dõi trung vị 2,7 năm. PFI: còn nguy cơ 54 ca ở 3 năm, 31 ca ở 5 năm; sự kiện muộn nhất 4,6 năm |
| 04/11/2026 | Bước 9 | METABRIC: TNBC theo 2 định nghĩa; cách đo PR, HER2; sự kiện RFS, OS | Không có PR_IHC, HER2_IHC → định nghĩa IHC chặt: 0 ca; tổ hợp: 276 ca. RFS 115 sự kiện (97 nếu cắt 5 năm); OS 146 (92); DSS 105 (85); theo dõi trung vị 14,6 năm. 18 ca ER bất đồng; 161 hóa trị; 64 điều trị nội tiết; PAM50: 151 Basal, 87 claudin-low (32%) |
| 04/11/2026 | Bước 10 | Kết quả kiểm tra Synapse / EGA | Chưa có dữ liệu miRNA METABRIC |
| 04/11/2026 | Bước 11 | QC thành phần mô | Cặp miRNA–mRNA âm có ý nghĩa (FDR < 0,05): 168.213 (thô) → 58.511 (hiệu chỉnh miễn dịch, mô đệm, biểu mô, purity) → 6.368 (hiệu chỉnh 5 PC); cặp dương thô 200.067. miR-200c-3p ~ ZEB1: −0,34 → −0,11; miR-21-5p ~ PDCD4: −0,22 → −0,22. 9 miRNA có \|ρ\| ≥ 0,5 với điểm miễn dịch (miR-150-5p, miR-142-5p, miR-4772-3p, miR-155-5p, miR-142-3p, miR-342-3p, miR-29c-3p, miR-146a-5p, miR-342-5p). PC2 (z theo gen) liên quan plate giải trình tự (p = 5·10⁻⁸) |

## Mốc kiểm tra 11/2026 – quyết định

Vương xác nhận ngày 04/11/2026 (chọn cả 8 phương án đề xuất).

| Mã | Quyết định | Lý do |
| --- | --- | --- |
| D1 | Kích hoạt phương án dự phòng: CoxNet là phân tích chính; RSF và GBS thăm dò (1 cấu hình bảo thủ). Đánh giá nội bộ chính: bootstrap hiệu chỉnh lạc quan (Harrell, 200 lần); repeated 5-fold CV là phụ. AUC tại 2 và 3 năm; Uno C cắt tại τ = 3 năm. Signature tối đa 5 cặp | 26 sự kiện PFI, mỗi fold chỉ 4–5 sự kiện |
| D2 | Kiểm chứng hai chiều cho nhánh chỉ dùng mRNA: huấn luyện METABRIC (RFS cắt 5 năm, 97 sự kiện) → kiểm tra TCGA (PFI) | METABRIC giàu thông tin gấp khoảng 4 lần |
| D3 | Huấn luyện: PFI. Kiểm chứng chính: RFS cắt ở 5 năm. Phụ: DSS ↔ "Died of Disease", OS. Báo cáo loại sự kiện PFI | TCGA theo dõi ngắn; TNBC tái phát sớm; 5/26 sự kiện PFI là u nguyên phát mới |
| D4 | Lọc cặp bằng Spearman từng phần (hiệu chỉnh purity, miễn dịch, mô đệm), chạy một lần trên toàn bộ u nguyên phát ER− của TCGA-BRCA (232 bệnh nhân có cả miRNA và mRNA), không dùng sống còn. Lọc trong từng vòng CV là độ nhạy | Lọc không giám sát không rò rỉ kết cục; tương quan thô chủ yếu phản ánh thành phần mô |
| D5 | f_rank (hiệu phân vị trong mẫu) là đặc trưng chính; f_raw và f_ref là độ nhạy | Chỉ f_rank bất biến với biến đổi đơn điệu khi chuyển RNA-seq → microarray |
| D6 | HER2 theo ASCO/CAP 2018: đưa lại 17 ca thiếu HER2 IHC nhưng ISH âm; loại ca HER2 mâu thuẫn | Đúng hướng dẫn; thêm 4 sự kiện |
| D7 | Giữ ca ER/PR "<10%", LAR, ESR1 cao trong cohort chính; gắn cờ và làm độ nhạy. TCGA không tách được ER-low (1–9%) | Cohort chính khớp thực hành lâm sàng |
| D8 | miRNA ≥10 read ở ≥50% mẫu; mRNA count ≥10 và TPM ≥1 ở ≥50% mẫu; chỉ giữ phân tử đo được trên cả hai nền tảng | RPM = 1 chỉ 0,8–2,3 read, là nhiễu đếm |

Đề cương đã cập nhật theo các quyết định này: `docs/proposal_miRNA-mRNA_v2_theo_doi_thay_doi.docx` (tóm tắt ở `docs/thay_doi_de_cuong_v2.md`).

## Mốc kiểm tra cuối 11/2026

- Số sự kiện PFI: 26 (152 bệnh nhân; trước khi sửa định nghĩa HER2 là 22/136).
- Quyết định: áp phương án dự phòng (D1–D8).
- Lý do: dưới ngưỡng 40 sự kiện của đề cương; chỉ 31 ca còn nguy cơ ở 5 năm.

## Phát hiện từ phân tích mô tả (11/2026)

Chưa nhìn kết cục. Phân tích tạm ngoài notebook (04/11/2026, trên cohort cũ 136 ca) gợi ra 5 nhóm vấn đề; số chính thức là số trong `qc/` (bảng trên).

1. **Ít sự kiện, theo dõi ngắn:** cohort kiểm chứng METABRIC có nhiều sự kiện hơn cohort huấn luyện.
2. **Tương quan âm phản ánh thành phần mô:** phần lớn cặp âm có ý nghĩa biến mất sau khi hiệu chỉnh miễn dịch, mô đệm, biểu mô và purity; miRNA của tế bào miễn dịch (miR-142, miR-150, miR-155) tương quan mạnh với điểm miễn dịch.
3. **miRNA-seq mang tính thành phần:** miR-21-5p chiếm khoảng một phần tư tổng RPM.
4. **Lỗi định lượng gen chồng lấn trong STAR – Counts:** CD68 và EIF4A1 có số đếm thấp bất thường.
5. **Hai cohort không đồng nhất:** LAR, claudin-low, ER bất đồng, khác thời kỳ điều trị.
