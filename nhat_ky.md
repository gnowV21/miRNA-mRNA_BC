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
| 04/11/2026 | Bước 3 | Số miRNA trưởng thành; kiểm tra miR-21-5p | 2013 MIMAT (loại 56 MIMAT không có trong miRBase v22). miR-21-5p chiếm 24,8% tổng RPM (trung bình), top 10 chiếm 67%. MIMAT0017985 (trung vị 112 read, có ở 97% mẫu) và MIMAT0018073 (trung vị 13 read) bị loại vì không có trong v22. Tra 04/11/2026: MIMAT0017985 = **hsa-miR-3607-3p** (ACUGUAAACGCUUUCUGAUG), còn trên trang tra cứu miRBase nhưng không có trong `mature.fa` chính thức. **Quyết định (Vương, 04/11/2026):** loại khỏi phân tích chính theo quy tắc định trước; giữ số liệu và trình tự để chạy độ nhạy ở ND2 (`data/b5_mir_*_khong_v22.parquet`, `data/mirbase_bo_sung_do_nhay.csv`). MIMAT0018073 chưa tra |
| 04/11/2026 | Bước 4 | Số gen mã hóa protein | 19.944. N_ambiguous trung vị 8,8% read mỗi mẫu. 9.154 gen có exon chồng lấn (GENCODE v36), 1.230 gen ≥50% exon bị chồng; CD68 (85,6%) và EIF4A1 được gắn cờ |
| 04/11/2026 | Bước 5 | Số bệnh nhân cuối có cả miRNA và mRNA | 152 (trước đây 136) |
| 04/11/2026 | Bước 5c | Đối chiếu STAR – RSEM (TOIL) | 152/152 mẫu có trong TOIL. Trung vị ρ Spearman STAR–RSEM theo gen 0,964. 726 gen `lech_star` (TPM STAR/RSEM < 0,25, có exon chồng lấn, RSEM ≥ 1 TPM); 815 gen nếu bỏ điều kiện chồng lấn. Xác nhận lỗi định lượng: CD68 (ρ 0,43; trung vị 0,04 so với 98,8 TPM) và EIF4A1 (ρ 0,39; 0,9 so với 535 TPM); PTPRC, ZEB1, PDCD4 khớp tốt (ρ ≥ 0,98). Phần lớn gen `lech_star` là gen biểu hiện cao (protein ribosome, ty thể, ALDOA). Trong bộ 12.825 gen giữ lại: 300 gen `lech_star` (có TPM1, đích kinh điển của miR-21) và 223 gen khác có ρ < 0,8. 352/19.944 gen mã hóa protein không có trong TOIL (v23), trong đó 2 gen thuộc bộ giữ lại. Xử lý: xem **D9** (loại 314 gen) |
| 04/11/2026 | Bước 6–7 | 4 ngưỡng lọc + lý do; số miRNA, gen giữ lại | miRNA ≥10 read ở ≥50% mẫu; mRNA count ≥10 và TPM ≥1 ở ≥50% mẫu (D8). Lý do: RPM = 1 chỉ tương ứng 0,79–2,33 read. Giữ 333 miRNA; 13.179 gen qua ngưỡng, 12.825 gen sau giao với METABRIC; **12.511 gen phân tích chính** sau khi bỏ 314 gen D9. Bản RSEM cho độ nhạy D9: 12.614 gen (11.918 chung với bản chính; 696 chỉ có ở RSEM, trong đó 587 gen D9 hoặc bị STAR đếm thiếu nên không qua ngưỡng; 593 chỉ có ở bản chính) |
| 04/11/2026 | Bước 7b | Chuẩn hóa miRNA thay thế | Spearman trung vị theo miRNA so với log2RPM gốc: bỏ miR-21 0,987; upper-quartile 0,944; DESeq2 0,946 |
| 04/11/2026 | Bước 8 | Sự kiện PFI / OS / DFI; theo dõi trung vị | PFI 26/152; OS 23/152; DFI 19/140; DSS 14/147. Theo dõi trung vị 2,7 năm. PFI: còn nguy cơ 54 ca ở 3 năm, 31 ca ở 5 năm; sự kiện muộn nhất 4,6 năm |
| 04/11/2026 | Bước 9 | METABRIC: TNBC theo 2 định nghĩa; cách đo PR, HER2; sự kiện RFS, OS | Không có PR_IHC, HER2_IHC → định nghĩa IHC chặt: 0 ca; tổ hợp: 276 ca. RFS 115 sự kiện (97 nếu cắt 5 năm); OS 146 (92); DSS 105 (85); theo dõi trung vị 14,6 năm. 18 ca ER bất đồng; 161 hóa trị; 64 điều trị nội tiết; PAM50: 151 Basal, 87 claudin-low (32%) |
| 04/11/2026 | Bước 10 | Kết quả kiểm tra Synapse / EGA | Chưa có dữ liệu miRNA METABRIC |
| 04/11/2026 | Bước 11 | QC thành phần mô | Cặp miRNA–mRNA âm có ý nghĩa (FDR < 0,05): tính trên 333 miRNA × 12.511 gen phân tích chính (sau D9): 165.509 (thô) → 57.566 (hiệu chỉnh miễn dịch, mô đệm, biểu mô, purity) → 6.267 (hiệu chỉnh 5 PC); cặp dương 196.113 → 65.745 → 19.297. (Trước D9, trên 12.825 gen: 168.213 → 59.219 → 6.368.) Điểm thành phần tính trên 17.383 gen chung hai nền tảng, đã bỏ gen D9 ở cả hai cohort. miR-200c-3p ~ ZEB1: −0,34 → −0,11; miR-21-5p ~ PDCD4: −0,22 → −0,22 (sau 5 PC: −0,19 và −0,11). 9 miRNA có \|ρ\| ≥ 0,5 với điểm miễn dịch (miR-150-5p, miR-142-5p, miR-4772-3p, miR-155-5p, miR-142-3p, miR-342-3p, miR-29c-3p, miR-146a-5p, miR-342-5p). PC2 (z theo gen) liên quan plate giải trình tự (p = 5·10⁻⁸) |

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

### D9 – Gen bị STAR đếm thiếu (04/11/2026)

Vương chốt ngày 04/11/2026, trước khi tạo bất kỳ cặp nào ở ND2. Quyết định không dùng dữ liệu sống còn.

- **Tiêu chí:** loại gen `lech_star` (exon chồng lấn; TPM STAR thấp hơn RSEM của TOIL từ 4 lần trở lên; RSEM ≥ 1 TPM), cộng với gen có ρ(STAR, RSEM) < 0,8 **và** tỷ số STAR/RSEM ngoài khoảng [0,25; 4]. Số gen bị loại khỏi bộ 12.825: **314** (300 `lech_star` + 14), còn **12.511** gen.
  - 14 gen thêm: 3 paralog bị multimap không có exon chồng lấn (EIF3C, NOMO2, LSM8) và 11 gen STAR cao hơn RSEM từ 4 lần trở lên, nhiều khả năng do mô hình gen đổi từ GENCODE v23 sang v36 (ALOX12, C22orf23, CFAP126, CYP3A5, H4C9, HLA-G, INCA1, NOX1, PALM2AKAP2, POMK, TRPV1).
  - Đề xuất ban đầu là loại mọi gen có ρ < 0,8. Cách này thêm 223 gen, nhiều hơn mức "vài chục" đã đặt, nên đã báo lại trước khi chốt. 209/223 gen có mức biểu hiện ở STAR và RSEM gần nhau (tỷ số 0,25–4), chủ yếu là gen dài, nhiều isoform (UBR4, PRPF8, SRRM2, MDN1, HERC2, DROSHA, MDM2, SPTAN1). Với các gen này chưa rõ STAR hay RSEM sai, vì RSEM của TOIL chia read theo isoform bằng EM trên chú thích v23 và không phải chuẩn vàng. Chỉ 13 gen là ZNF, nên nhóm này không phải chủ yếu là paralog. Tỷ lệ gen có ρ < 0,8 gần như không đổi theo mức biểu hiện (2–5% ở mỗi nhóm ngũ phân vị), nên đây cũng không phải nhiễu do biểu hiện thấp. Các gen này chỉ gắn cờ `rho_thap` (209 gen) và được phân tích độ nhạy bằng RSEM bao phủ.
  - 2 gen không có trong TOIL không kiểm chứng được và được giữ lại.
- **Phân tích chính:** loại khỏi tập ứng viên tạo cặp, khỏi tập gen tính phân vị của f_rank (cả hai cohort: METABRIC lấy đúng các gen trong `data/b7_gen.csv` có `phan_tich_chinh`), khỏi baseline biểu hiện gen thô, và khỏi ứng viên gen tham chiếu cho f_ref. Điểm thành phần mô ở Bước 9b cũng bỏ các gen này ở cả hai cohort. Giữ trong đầu vào ND3 (`b7_mr_log2tpm_tat_ca.parquet`), ghi chú là đầu vào có gen bị STAR đếm thiếu.
- **Độ nhạy:** chạy lại bước chọn cặp ND2 với RSEM TPM cho toàn bộ gen (`data/b7_mr_log2tpm_rsem.parquet`: mã v23 ánh xạ sang v36 theo mã Ensembl gốc, đã chuyển log2(TPM + 0,001) về TPM, áp ngưỡng TPM ≥ 1 ở ≥50% mẫu; TOIL không có số đếm, nhưng với STAR điều kiện count chỉ loại thêm 4 gen). Báo cáo các cặp chỉ xuất hiện ở bản RSEM (ví dụ miR-21-5p–TPM1). Nếu nhiều cặp mạnh chỉ có ở bản RSEM thì cân nhắc lại việc dùng STAR làm nguồn chính.
- **Lý do:** STAR – Counts bỏ các read mơ hồ ở vùng chồng lấn. Kiểm chứng: CD68 có ρ = 0,43 (STAR 0,04 so với RSEM 98,8 TPM), EIF4A1 có ρ = 0,39. Các gen đối chứng đều có ρ ≥ 0,98. Không chuyển hẳn sang RSEM vì đề cương và pipeline dùng STAR của GDC (có md5, phiên bản dữ liệu), còn TOIL là bộ tính lại của bên thứ ba trên chú thích cũ hơn; loại 314 gen chỉ mất 2,4% số gen.
- **Ghi chú cho ND5:** immunedeconv chạy trên TPM của STAR sẽ ước lượng thiếu đại thực bào (CD68 có trong nhiều bộ gen đại thực bào). Chạy thêm trên RSEM TPM, hoặc chỉ kết luận khi hai nguồn nhất quán.

## Mốc kiểm tra cuối 11/2026

- Số sự kiện PFI: 26 (152 bệnh nhân; trước khi sửa định nghĩa HER2 là 22/136).
- Quyết định: áp phương án dự phòng (D1–D8).
- Lý do: dưới ngưỡng 40 sự kiện của đề cương; chỉ 31 ca còn nguy cơ ở 5 năm.

## Phát hiện từ phân tích mô tả (11/2026)

Chưa nhìn kết cục. Phân tích tạm ngoài notebook (04/11/2026, trên cohort cũ 136 ca) gợi ra 5 nhóm vấn đề; số chính thức là số trong `qc/` (bảng trên).

1. **Ít sự kiện, theo dõi ngắn:** cohort kiểm chứng METABRIC có nhiều sự kiện hơn cohort huấn luyện.
2. **Tương quan âm phản ánh thành phần mô:** phần lớn cặp âm có ý nghĩa biến mất sau khi hiệu chỉnh miễn dịch, mô đệm, biểu mô và purity; miRNA của tế bào miễn dịch (miR-142, miR-150, miR-155) tương quan mạnh với điểm miễn dịch.
3. **miRNA-seq mang tính thành phần:** miR-21-5p chiếm khoảng một phần tư tổng RPM.
4. **Lỗi định lượng gen chồng lấn trong STAR – Counts:** CD68 và EIF4A1 có số đếm thấp bất thường. Đối chiếu với RSEM (Bước 5c) xác nhận: 726 gen bị STAR đếm thấp hơn 4 lần trở lên. Xử lý ở D9.
5. **Hai cohort không đồng nhất:** LAR, claudin-low, ER bất đồng, khác thời kỳ điều trị.
