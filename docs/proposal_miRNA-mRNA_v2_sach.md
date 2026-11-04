# Xây dựng mô hình tiên lượng tái phát ung thư vú bộ ba âm tính dựa trên đặc trưng cặp tương tác miRNA–mRNA và so sánh với biểu diễn từ mô hình nền tảng phiên mã

Oct 17, 2026 · @Vuong

## Thông tin chung

| Mục | Nội dung |
|---|---|
| Loại đề cương | Đề cương khóa luận tốt nghiệp |
| Sinh viên thực hiện | Trần Nguyễn Tiến Vương |
| Mã số sinh viên | 23180211 |
| Lớp | 23CSH3 |
| Ngành | Công nghệ Sinh học |
| Khoa | Sinh học – Công nghệ Sinh học |
| Trường | Trường Đại học Khoa học Tự nhiên, ĐHQG-HCM (HCMUS) |
| Giảng viên hướng dẫn | TS. Huỳnh Hữu Luân; PGS. TS. Nguyễn Thị Huệ |
| Địa điểm thực hiện | HG Lab |
| Thời gian thực hiện | 11/2026 – 8/2027 (bảo vệ dự kiến khoảng 8/2027) |
| Loại hình nghiên cứu | Nghiên cứu tính toán (in silico), dữ liệu công khai |

## 1. Đặt vấn đề và tính cấp thiết

Đề tài xây dựng một mô hình tiên lượng tái phát cho ung thư vú bộ ba âm tính (TNBC), dùng đặc trưng cặp miRNA–mRNA tính được trên từng mẫu, và kiểm tra xem đặc trưng đơn giản này có bền vững hơn biểu diễn từ mô hình nền tảng hay không.

**Nhu cầu lâm sàng.**  TNBC chiếm khoảng 10–15% ca ung thư vú, không có thụ thể ER, PR và HER2. Đây là phân nhóm ác tính nhất, ít liệu pháp nhắm đích, và nguy cơ tái phát sớm cao. Bác sĩ hiện thiếu công cụ phân tầng nguy cơ đáng tin cậy riêng cho nhóm này.

**Vai trò của miRNA.**  miRNA điều hòa sau phiên mã hơn 60% gen mã hóa protein ở người. Tác động của một miRNA lên khối u không nằm ở mức biểu hiện của riêng nó, mà ở mức nó ức chế các mRNA đích. Vì vậy, cặp miRNA–mRNA phản ánh hoạt tính điều hòa sát hơn so với từng phân tử đơn lẻ.

**Ba vấn đề của lĩnh vực hiện nay:**

- **Signature tiên lượng khó tái lập.**  Phần lớn signature miRNA từ TCGA chỉ được đánh giá trên một nền tảng đo, và hiệu năng giảm mạnh khi chuyển sang cohort hoặc công nghệ khác (RNA-seq sang microarray).

- **Khủng hoảng đánh giá trong dự đoán đích miRNA.**  Khi đánh giá bằng mẫu âm tính xác thực CLIP-seq (bộ miRBench), nhiều mô hình công bố hiệu năng trên 0,90 chỉ còn AU-PRC 0,55–0,80 (Lin et al., 2026). Chọn cặp tương tác cần dựa trên công cụ đã qua đánh giá nghiêm ngặt.

- **Mô hình nền tảng chưa được kiểm chứng đủ.**  Các mô hình nền tảng phiên mã (BulkRNABert, BulkFormer) đang được dùng cho tiên lượng, nhưng nhiều đánh giá độc lập cho thấy mô hình phức tạp thường không vượt baseline tuyến tính đơn giản (Ahlmann-Eltze et al., 2025).

**Câu hỏi nghiên cứu.**  Ở bệnh nhân TNBC, một chỉ số đơn giản tính từ cặp miRNA–mRNA có dự đoán nguy cơ tái phát hoặc tiến triển bệnh (PFI) tốt bằng hoặc tốt hơn đặc trưng do mô hình nền tảng tạo ra không?

Chỉ số đơn giản ở đây là hiệu số (theo phân vị hoặc theo log) giữa mức một miRNA và mức một mRNA đích của nó, tính ngay trên mẫu của từng bệnh nhân. Đặc trưng chính là hiệu phân vị của miRNA và mRNA trong cùng mẫu; hiệu số log và hiệu số log chuẩn hóa theo phân tử tham chiếu (dạng qPCR) được báo cáo làm phân tích độ nhạy. Đặc trưng do mô hình nền tảng tạo ra là embedding từ BulkRNABert và BulkFormer. Câu hỏi được tách thành ba câu nhỏ:

- Trên cohort huấn luyện TCGA (đo bằng RNA-seq), đặc trưng cặp có phân biệt bệnh nhân nguy cơ cao và thấp tốt bằng embedding không, đo bằng C-index của Uno (cắt tại 3 năm) sau hiệu chỉnh lạc quan bằng bootstrap, và so với các mô hình đối chứng null (mục 6.5)?

- Khi áp nguyên mô hình đã huấn luyện sang METABRIC (đo bằng microarray), loại đặc trưng nào giữ được hiệu năng tốt hơn khi các đặc trưng được định nghĩa bất biến với thang đo của nền tảng (hiệu phân vị trong mẫu)?

- Các cặp miRNA–mRNA trong mô hình cuối liên quan đến cơ chế sinh học nào của TNBC sau khi tách ảnh hưởng của thành phần mô (độ tinh khiết khối u, thâm nhiễm miễn dịch, mô đệm)?

## 2. Tổng quan tài liệu

Các hướng tiếp cận liên quan đã có, nhưng chưa nghiên cứu nào kết hợp đặc trưng cặp đơn mẫu, kiểm chứng xuyên nền tảng ở TNBC, và so sánh trực tiếp với mô hình nền tảng phiên mã.

### 2.1. Dự đoán đích miRNA

Công cụ cổ điển như TargetScan dựa trên đặc trưng thủ công: loại seed match, độ bảo tồn, năng lượng duplex. Thế hệ mới học trực tiếp từ trình tự. **miRBind2**  (*Bioinformatics* , 2026) dùng biểu diễn cặp nucleotide với kiến trúc CNN, được tinh chỉnh trên dữ liệu nhiễu loạn miRNA, và vượt TargetScan trên 50.549 cặp miRNA–gen (miRBind2, MDC Berlin, 2026). Công cụ này dự đoán cả mức ức chế chức năng ở cấp gen, không chỉ gắn kết vật lý.

### 2.2. Signature tiên lượng và mạng đặc hiệu mẫu

Nhiều nhóm đã chuyển từ signature miRNA đơn lẻ sang tương tác miRNA–mRNA:

| Nghiên cứu | Cách tiếp cận | Giới hạn liên quan đến đề tài |
|---|---|---|
| spongEffects | Điểm hoạt tính module ceRNA cho từng bệnh nhân | Dùng cho phân loại, chưa kiểm chứng tiên lượng xuyên nền tảng |
| BONOBO (Genome Research, 2024) | Mạng đồng biểu hiện miRNA–mRNA đặc hiệu mẫu, liên hệ với sống còn theo phân nhóm ung thư vú | Mức pathway, không ra signature dùng được trên từng bệnh nhân mới |
| Mô hình MTI cho TNBC (IJMS, 2025) | Chọn 8 miRNA qua tương quan điều hòa và ML | Đặc trưng vẫn là mức biểu hiện miRNA, không phải đặc trưng cặp |
| Mạng ceRNA đặc hiệu bệnh nhân | Bộ ba lncRNA–miRNA–mRNA tiên lượng ở ung thư vú, gan, phổi | Không phân tích riêng TNBC, không so với mô hình nền tảng |

### 2.3. Mô hình nền tảng phiên mã trong tiên lượng ung thư

- **BulkRNABert**  tiền huấn luyện kiểu BERT trên RNA-seq bulk, gồm cả dữ liệu TCGA, dùng cho phân loại ung thư và dự đoán sống còn (Gélard et al.).

- **BulkFormer**  có 150 triệu tham số, tiền huấn luyện trên hơn 500.000 hồ sơ phiên mã bulk, báo cáo vượt các mô hình đơn bào ở tác vụ tiên lượng (bioRxiv).

- Ở một bài toán khác là dự đoán thay đổi biểu hiện gen sau nhiễu loạn bằng CRISPR, một đánh giá độc lập cho thấy cả năm mô hình nền tảng thử nghiệm đều không vượt các baseline đơn giản (Ahlmann-Eltze et al., *Nature Methods*  2025). Các baseline đó gồm: mô hình "không thay đổi" (dự đoán biểu hiện giống mẫu đối chứng), mô hình cộng (cộng hiệu ứng của từng gen khi bị nhiễu loạn riêng) và mô hình tuyến tính trên biểu diễn PCA. Baseline tương ứng trong đề tài này là biểu hiện gen thô đưa thẳng vào mô hình sống còn (mục 6.3).

### 2.4. Khoảng trống nghiên cứu

- Chưa có signature tiên lượng TNBC nào dựa trên đặc trưng cặp miRNA–mRNA tính trong từng mẫu, không cần chuẩn hóa cấp cohort và bất biến với thang đo của nền tảng.

- Chưa có so sánh trực tiếp giữa đặc trưng sinh học đơn giản và embedding mô hình nền tảng khi chuyển nền tảng đo (RNA-seq sang microarray).

- Việc chọn cặp tương tác trong các nghiên cứu tiên lượng chưa dựa trên công cụ dự đoán đích đã qua đánh giá với mẫu âm tính CLIP-seq.

**Chuẩn hóa cấp cohort**  là biến đổi giá trị của một bệnh nhân bằng các con số tính từ cả nhóm bệnh nhân, như trung bình và độ lệch chuẩn (z-score), phân bố chung (quantile normalization) hay ngưỡng trung vị. Khi đó, cùng một mẫu có thể được xếp là "biểu hiện cao" ở cohort này nhưng "biểu hiện thấp" ở cohort khác. Mô hình phụ thuộc vào chuẩn hóa cấp cohort vì vậy khó áp cho từng bệnh nhân đơn lẻ trên lâm sàng, và dễ lệch khi chuyển sang cohort đo bằng nền tảng khác.

## 3. Mục tiêu nghiên cứu

**Mục tiêu tổng quát.**  Xây dựng và kiểm chứng ngoài một mô hình tiên lượng thời gian không tiến triển (PFI) cho bệnh nhân TNBC dựa trên đặc trưng cặp miRNA–mRNA, đồng thời đánh giá giá trị của đặc trưng này so với biểu diễn từ mô hình nền tảng phiên mã.

**Mục tiêu cụ thể:**

- Xây dựng pipeline thu thập và tiền xử lý dữ liệu miRNA, mRNA và lâm sàng cho bệnh nhân TNBC từ TCGA-BRCA và METABRIC.

- Xác định tập cặp miRNA–mRNA ứng viên bằng miRBind2, có tương quan âm từng phần (sau hiệu chỉnh thành phần mô) trong mô u ER−, và tính đặc trưng cặp trong từng mẫu (hiệu phân vị; log-ratio làm độ nhạy).

- Huấn luyện mô hình CoxNet (phân tích chính) trên bốn bộ đặc trưng: cặp miRNA–mRNA, embedding BulkRNABert, embedding BulkFormer, biểu hiện gen thô; Random Survival Forest và Gradient Boosting Survival là phân tích thăm dò.

- Kiểm chứng ngoài mô hình tốt nhất trên METABRIC, và so sánh với yếu tố lâm sàng và các signature đã công bố; so sánh với các mô hình đối chứng null (cặp ngẫu nhiên, đặc trưng riêng lẻ, lâm sàng + miễn dịch); với nhánh chỉ dùng mRNA, kiểm chứng thêm theo chiều ngược METABRIC → TCGA.

- Diễn giải sinh học signature qua SHAP, phân tích làm giàu chức năng, vi môi trường miễn dịch và dự đoán đáp ứng thuốc.

## 4. Đối tượng và phạm vi nghiên cứu

**Đối tượng.**  Bệnh nhân ung thư vú nguyên phát thuộc phân nhóm bộ ba âm tính, xác định theo hóa mô miễn dịch (IHC): ER âm tính, PR âm tính, HER2 âm tính (HER2 theo IHC hoặc ISH, ASCO/CAP 2018).

| Thành phần | Cohort huấn luyện | Cohort kiểm chứng ngoài |
|---|---|---|
| Bộ dữ liệu | TCGA-BRCA | METABRIC |
| Nền tảng đo miRNA | miRNA-seq (giải trình tự) | Microarray |
| Nền tảng đo mRNA | RNA-seq | Microarray |
| Mẫu | Khối u nguyên phát (Primary Tumor) | Khối u nguyên phát |
| Endpoint | PFI theo TCGA-CDR (chính); DSS, OS (phụ) | RFS cắt tại 5 năm (chính); DSS (Died of Disease), OS (phụ) |
| Định nghĩa TNBC | IHC từ dữ liệu lâm sàng: ER−, PR− theo IHC; HER2 âm theo ASCO/CAP 2018 (Wolff et al., 2018): IHC 0/1+, hoặc IHC 2+ kèm ISH âm, hoặc không có IHC nhưng ISH âm. Loại ca có kết quả HER2 mâu thuẫn. TCGA ghi ER/PR "<10%" nên không tách được nhóm ER-low (1–9%; Allison et al., 2020); nhóm này được gắn cờ và phân tích độ nhạy | Định nghĩa chính (IHC): chỉ giữ ca có kết quả nhuộm IHC cho cả ER, PR, HER2 và âm tính cả ba. Định nghĩa tổ hợp (chạy song song nếu METABRIC thiếu IHC cho PR hoặc HER2) ghép ba yếu tố đo bằng các phương pháp khác nhau: (1) ER âm theo IHC (cột ER_IHC); (2) PR âm theo trường trạng thái PR có sẵn (cột PR_STATUS, nhiều khả năng suy từ biểu hiện mRNA); (3) HER2 âm theo trường trạng thái HER2 có sẵn (cột HER2_STATUS, dựa trên số bản sao gen ERBB2 đo bằng chip SNP6 hoặc biểu hiện). Nguồn đo thực tế của từng cột được xác minh ở Bước 9. Kết quả xác minh (11/2026): bản công khai không có IHC cho PR và HER2, nên định nghĩa IHC chặt cho 0 ca; dùng định nghĩa tổ hợp (276 ca TNBC; 18 ca ER_IHC âm nhưng ER_STATUS dương được gắn cờ). |

**Phạm vi khóa luận.**  Toàn bộ phần tính toán: thu thập dữ liệu, xây dựng đặc trưng, huấn luyện, kiểm chứng ngoài và diễn giải sinh học. Bản thảo bài báo khoa học và kiểm chứng thực nghiệm nằm ngoài phạm vi khóa luận.

**Giới hạn.**  Chỉ dùng dữ liệu công khai, không thu mẫu bệnh nhân mới. Tài nguyên tính toán: máy cá nhân cho Nội dung 1; Google Colab và Kaggle Notebooks miễn phí cho các bước cần GPU; thuê GPU theo giờ khi cần cho bước trích embedding. Dữ liệu METABRIC thuộc thời kỳ điều trị cũ (161/276 ca được hóa trị, 64 ca ER− có điều trị nội tiết), nên phân tích kiểm chứng được hiệu chỉnh theo hóa trị.

## 5. Nội dung nghiên cứu

Đề tài gồm năm nội dung, mỗi nội dung ứng với một mục tiêu cụ thể ở mục 3.

- **Nội dung 1 – Dữ liệu.**  Tải, lọc bệnh nhân TNBC, ghép dữ liệu miRNA, mRNA và lâm sàng, chú thích sinh học (PAM50, độ tinh khiết, thâm nhiễm miễn dịch), chuẩn hóa trong từng cohort. Sản phẩm: bộ dữ liệu sạch và báo cáo số ca, số sự kiện.

- **Nội dung 2 – Đặc trưng cặp miRNA–mRNA.**  Chạy miRBind2 cho các cặp ứng viên, lọc theo tương quan âm từng phần, tính đặc trưng cặp trong từng mẫu (hiệu phân vị; log-ratio làm độ nhạy). Sản phẩm: ma trận đặc trưng cặp cho cả hai cohort.

- **Nội dung 3 – Đặc trưng so sánh.**  Trích embedding BulkRNABert và BulkFormer từ dữ liệu mRNA; chuẩn bị biểu hiện gen thô làm baseline. Sản phẩm: bốn bộ đặc trưng cùng tập bệnh nhân.

- **Nội dung 4 – Mô hình hóa và kiểm chứng.**  Huấn luyện CoxNet (chính; RSF và GBS thăm dò) trên bốn bộ đặc trưng, đánh giá bằng bootstrap hiệu chỉnh lạc quan; kiểm chứng ngoài trên METABRIC. Sản phẩm: bảng so sánh hiệu năng và mô hình cuối.

- **Nội dung 5 – Diễn giải sinh học.**  Phân tích SHAP, làm giàu chức năng, vi môi trường miễn dịch và đáp ứng thuốc cho signature cuối. Sản phẩm: danh sách trục điều hòa ưu tiên kèm bằng chứng sinh học.

Sơ đồ quy trình tổng thể được trình bày ở mục 6.

### 5.1. Nội dung 1 – Thu thập và tiền xử lý dữ liệu

Mục tiêu: có hai bộ dữ liệu TNBC sạch, cùng định dạng, và biết chính xác số ca, số sự kiện trước cuối tháng 11/2026.

**Các bước thực hiện**

- Tải từ GDC qua API bằng Python: Isoform Expression Quantification (miRNA mức trưởng thành), Gene Expression Quantification (STAR - Counts), bảng lâm sàng dạng BCR Biotab. Tải bảng TCGA-CDR để lấy PFI và PFI.time. Ghi phiên bản dữ liệu GDC (data release) và md5 của mọi file.

- Lọc mẫu: chỉ giữ khối u nguyên phát (mã loại mẫu 01); bệnh nhân có nhiều mẫu thì giữ một mẫu theo quy tắc cố định, ghi trong phương pháp.

- Xác định TNBC theo IHC: ER âm, PR âm, HER2 âm. HER2 âm theo ASCO/CAP 2018 (Wolff et al., 2018): IHC 0/1+, hoặc IHC 2+ kèm ISH âm, hoặc không có IHC nhưng ISH âm; loại ca có kết quả HER2 mâu thuẫn (điểm IHC trái với trạng thái ghi nhận). Gắn cờ, không loại, các ca ER/PR "<10%", ESR1 cao, kiểu LAR. Ghi số ca bị loại ở từng bước.

- Tiền xử lý: lọc theo số read (miRNA ≥10 read; mRNA ≥10 count và TPM ≥1; ở ≥50% mẫu); chỉ giữ miRNA và gen đo được trên cả hai nền tảng; chuyển log2(x + 1) cho dữ liệu đếm (METABRIC đã ở thang log2). Ghép miRNA và mRNA theo mã mẫu.

- METABRIC: tải mRNA và lâm sàng (gồm RFS) từ cBioPortal; dữ liệu miRNA theo Dvinge et al. (2013). Chuẩn hóa tên miRNA về cùng phiên bản miRBase (file mature.fa của miRBase), ánh xạ probe microarray về tên gen.

- Chú thích sinh học và QC (không dùng sống còn): PAM50, độ tinh khiết ABSOLUTE (Carter et al., 2012), thâm nhiễm miễn dịch (Thorsson et al., 2018), điểm thành phần mô trong từng mẫu, cờ LAR và ESR1 cao, thuốc điều trị; cờ gen chồng lấn trong STAR – Counts và đối chiếu với RSEM của TOIL (Vivian et al., 2017); phân tích độ nhạy chuẩn hóa miRNA (bỏ miR-21-5p, upper-quartile, DESeq2).

- Báo cáo: sơ đồ chọn mẫu, số ca, số sự kiện, thời gian theo dõi trung vị (tính bằng phương pháp Kaplan–Meier đảo); số ca còn nguy cơ tại 1, 2, 3, 5 năm; thời điểm sự kiện; loại sự kiện PFI.

**Đơn vị kiến thức cần áp dụng**

- Cấu trúc dữ liệu TCGA/GDC: mã barcode mẫu và bệnh nhân, các loại mẫu, cách truy vấn API của GDC, định dạng Parquet.

- Phân nhóm ung thư vú: tiêu chí IHC cho ER, PR, HER2; cách chấm điểm HER2 (0, 1+, 2+, 3+) và vai trò FISH.

- Dữ liệu sống còn: thời gian theo dõi, kiểm duyệt (censoring), định nghĩa endpoint PFI và RFS.

- Tiền xử lý omics: lọc biểu hiện thấp, biến đổi log, khác biệt thang đo giữa RNA-seq và microarray.

- Danh pháp miRNA: dạng -5p và -3p, sự thay đổi tên giữa các phiên bản miRBase.

- Dữ liệu dạng thành phần (compositional) của miRNA-seq (Gloor et al., 2017); giới hạn của định lượng gen chồng lấn khi đếm read theo gen.

**Đầu ra và tiêu chí hoàn thành:**  hai đối tượng dữ liệu (TCGA, METABRIC) và sơ đồ chọn mẫu; số sự kiện PFI được báo cáo để đưa vào mốc kiểm tra cuối tháng 11/2026.

Hướng dẫn thực hiện chi tiết kèm kiến thức nền: Hướng dẫn Nội dung 1

### 5.2. Nội dung 2 – Đặc trưng cặp miRNA–mRNA

Mục tiêu: có danh sách cặp ứng viên có bằng chứng tính toán (miRBind2) và bằng chứng trong bệnh cảnh (tương quan âm từng phần, sau hiệu chỉnh thành phần mô, ở mô u ER−), cùng hàm tính đặc trưng cặp (hiệu phân vị, log-ratio, log-ratio chuẩn hóa tham chiếu) dùng lại được cho mọi cohort.

**Các bước thực hiện**

- Chuẩn bị trình tự: miRNA trưởng thành từ miRBase; trình tự 3'UTR từ Ensembl (qua BioMart), chọn một bản phiên mã đại diện cho mỗi gen theo quy tắc cố định.

- Chỉ chấm các miRNA và gen vượt ngưỡng biểu hiện ở Nội dung 1, để giảm số cặp cần tính.

- Chạy miRBind2 trên Colab theo lô; lưu điểm sau mỗi lô lên Google Drive để không mất kết quả khi phiên bị ngắt.

- Chọn ứng viên: với mỗi miRNA, giữ các gen thuộc nhóm điểm cao nhất theo ngưỡng phần trăm đặt trước khi xem dữ liệu sống còn.

- Tính tương quan Spearman từng phần miRNA–mRNA, hiệu chỉnh theo độ tinh khiết, điểm miễn dịch và mô đệm, trên toàn bộ u nguyên phát ER− của TCGA-BRCA (khoảng 232 bệnh nhân có cả miRNA và mRNA); hiệu chỉnh đa kiểm định Benjamini–Hochberg, giữ cặp có hệ số âm và FDR dưới ngưỡng. Bước này không dùng dữ liệu sống còn nên không rò rỉ kết cục, và được làm một lần trước khi mô hình hóa. Lọc trong từng vòng CV trên cohort TNBC là phân tích độ nhạy.

- Kiểm tra miRBind2: đích dự đoán có lệch về phía tương quan âm so với gen không phải đích cùng mức biểu hiện không (kiểm định KS một phía).

- Tính ba dạng đặc trưng cho từng bệnh nhân theo công thức ở mục 6.2: hiệu phân vị (chính), log-ratio, log-ratio chuẩn hóa theo phân tử tham chiếu (chọn không dùng kết cục).

- Kiểm tra: so sánh phân bố đặc trưng giữa TCGA và METABRIC; đo tương quan giữa các cặp dùng chung miRNA hoặc chung gen (đa cộng tuyến).

- Mở rộng nếu kịp: isomiR đầu 5′ (seed dịch chuyển nên tập đích khác; miRBind2 chấm trực tiếp trình tự isomiR; Telonis et al., 2015), tỷ lệ 5p/3p.

**Đơn vị kiến thức cần áp dụng**

- Cơ chế điều hòa của miRNA: sinh tổng hợp, phức hợp RISC, vùng seed, vị trí gắn trên 3'UTR, gắn kết không chính tắc.

- Dự đoán đích miRNA: khác biệt giữa gắn kết vật lý và ức chế chức năng; vì sao mẫu âm tính CLIP-seq làm thay đổi đánh giá.

- Deep learning ở mức vận hành: mã hóa trình tự làm đầu vào CNN, chạy suy luận theo lô, quản lý bộ nhớ GPU.

- Thống kê: tương quan hạng Spearman, đa kiểm định và FDR.

- Rò rỉ dữ liệu: vì sao mọi bước học từ dữ liệu phải nằm trong vòng cross-validation.

- Dữ liệu dạng tỷ lệ: ý nghĩa của log-ratio và vai trò của pseudo-count.

**Đầu ra và tiêu chí hoàn thành:**  bảng điểm miRBind2 cho toàn bộ cặp đã chấm; hàm Python tính đặc trưng từ một ma trận miRNA và mRNA bất kỳ; báo cáo số cặp sau từng bước lọc.

### 5.3. Nội dung 3 – Embedding mô hình nền tảng

Mục tiêu: có ma trận embedding BulkRNABert và BulkFormer cho cùng tập bệnh nhân ở cả hai cohort, đúng định dạng đầu vào mà mỗi mô hình yêu cầu.

**Các bước thực hiện**

- Đọc tài liệu từng mô hình: danh sách gen đầu vào, kiểu chuẩn hóa yêu cầu (ví dụ TPM, log), cách xử lý gen bị thiếu.

- Chuyển dữ liệu mRNA của TCGA và METABRIC sang định dạng đó: ánh xạ tên gen, sắp đúng thứ tự gen, xử lý gen thiếu theo hướng dẫn của tác giả.

- Chạy suy luận theo lô trên Colab (không huấn luyện lại); lưu embedding lên Google Drive.

- Kiểm tra chất lượng: vẽ PCA hoặc UMAP của embedding, tô màu theo cohort và theo phân nhóm PAM50. Nếu điểm tách theo cohort mạnh hơn theo sinh học, đó là dấu hiệu embedding nhạy với nền tảng đo.

- Giảm chiều bằng PCA bên trong từng vòng cross-validation, ở Nội dung 4.

- Nhánh chỉ dùng mRNA (embedding, gen thô) được kiểm chứng hai chiều: TCGA → METABRIC và METABRIC → TCGA.

- Tải phiên bản BulkRNABert tiền huấn luyện không có dữ liệu TCGA; ghi rõ phiên bản và nguồn trọng số trong phương pháp.

**Đơn vị kiến thức cần áp dụng**

- Mô hình nền tảng ở mức khái niệm: kiến trúc transformer, tiền huấn luyện tự giám sát kiểu che (masked modeling), embedding, học chuyển giao.

- Lập trình Python: thư viện transformers của Hugging Face, pandas, chạy suy luận trên GPU.

- Giảm chiều: PCA, phương sai giải thích, vì sao phải khớp PCA trên dữ liệu huấn luyện.

- Hiệu ứng lô và khác biệt nền tảng (RNA-seq và microarray).

- Rò rỉ thông tin qua dữ liệu tiền huấn luyện.

**Đầu ra và tiêu chí hoàn thành:**  bốn ma trận embedding (2 mô hình × 2 cohort); hình kiểm tra chất lượng; ghi chú các gen bị thiếu và cách đã xử lý.

### 5.4. Nội dung 4 – Mô hình hóa và kiểm chứng ngoài

Mục tiêu: so sánh CoxNet trên 4 bộ đặc trưng (phân tích chính; 12 tổ hợp mô hình × đặc trưng là phân tích thăm dò) trên TCGA, khóa một mô hình cuối, rồi đánh giá nó một lần duy nhất trên METABRIC.

**Các bước thực hiện**

- Viết kế hoạch phân tích trước khi chạy: chỉ số chính, tiêu chí chọn mô hình, ngưỡng nhóm nguy cơ. Lưu kèm ngày tháng để chứng minh không điều chỉnh sau khi thấy kết quả.

- Đánh giá nội bộ chính: bootstrap hiệu chỉnh lạc quan (200 lần, gồm toàn bộ pipeline); repeated stratified 5-fold CV (20 lần) làm phân tích phụ; chọn siêu tham số bằng CV lặp ở vòng trong theo quy tắc 1-SE.

- Đóng gói mọi bước học từ dữ liệu (lọc cặp, PCA, chuẩn hóa, chọn đặc trưng) thành một pipeline chạy lại bên trong mỗi lần bootstrap và mỗi phần của CV.

- Huấn luyện 3 mô hình trên 4 bộ đặc trưng; lưu dự đoán vòng ngoài của mọi tổ hợp.

- Tính chỉ số: C-index (Harrell; Uno cắt tại τ = 3 năm, Uno et al., 2011), AUC phụ thuộc thời gian tại 2 và 3 năm, Integrated Brier Score (đến 3 năm ở TCGA, 5 năm ở METABRIC), đường hiệu chuẩn. Kiểm tra giả định tỷ lệ nguy cơ của mô hình Cox bằng phần dư Schoenfeld. Nếu vi phạm, báo cáo HR theo khoảng thời gian (0–2, 2–5 năm).

- Đối chứng null: (i) đặc trưng riêng lẻ, (ii) cặp ngẫu nhiên cùng mức biểu hiện (200 lần, p thực nghiệm), (iii) lâm sàng + miễn dịch + purity, (iv) (iii) + đặc trưng cặp.

- Chọn mô hình cuối theo tiêu chí đã viết ở bước 1; khi chênh lệch nhỏ, ưu tiên mô hình đơn giản hơn.

- Huấn luyện lại mô hình cuối trên toàn bộ TCGA; khóa hệ số và ngưỡng chia nhóm nguy cơ.

- Áp mô hình đã khóa lên METABRIC với RFS cắt tại 5 năm: tính lại các chỉ số, vẽ Kaplan–Meier, Cox đa biến kèm hóa trị, giai đoạn, hạch, điểm miễn dịch, so sánh với signature đã công bố, decision curve analysis, so sánh bootstrap giữa các bộ đặc trưng. Phân tích độ nhạy: loại ca ER bất đồng; chỉ giữ ca có hóa trị.

**Đơn vị kiến thức cần áp dụng**

- Phân tích sống còn: hàm nguy cơ, kiểm duyệt, Kaplan–Meier, kiểm định log-rank, mô hình Cox và giả định tỷ lệ nguy cơ, tỷ số nguy cơ (HR).

- Học máy có phạt: LASSO, ridge, elastic-net; vì sao phạt giúp mô hình ổn định khi ít mẫu.

- Mô hình cây: random forest và gradient boosting cho dữ liệu sống còn; các siêu tham số chính và dấu hiệu overfitting.

- Thiết kế đánh giá: cross-validation và nested cross-validation, rò rỉ dữ liệu, số sự kiện trên mỗi biến (EPV).

- Đánh giá mô hình dự báo: phân biệt (C-index, AUC), hiệu chuẩn, lợi ích lâm sàng (DCA), bootstrap.

- Chuẩn báo cáo TRIPOD+AI.

**Đầu ra và tiêu chí hoàn thành:**  bảng so sánh CoxNet trên 4 bộ đặc trưng kèm đối chứng null (12 tổ hợp là thăm dò); mô hình cuối có công thức điểm nguy cơ; toàn bộ kết quả trên METABRIC; mã nguồn chạy lại ra cùng con số.

### 5.5. Nội dung 5 – Diễn giải sinh học

Mục tiêu: giải thích signature bằng cơ chế sinh học có kiểm chứng chéo giữa hai cohort, không chỉ báo cáo con số dự đoán. Thứ tự cố định: bước 1–2 bắt buộc; bước 3 rồi bước 4 làm nếu kịp tiến độ.

**Các bước thực hiện**

- SHAP trên mô hình cuối: xếp hạng tầm quan trọng toàn cục của từng cặp, và giải thích điểm nguy cơ cho từng bệnh nhân điển hình.

- Làm giàu chức năng theo hai cách: phân tích quá đại diện (ORA) cho các gen đích trong signature, với nền là các gen đã được chấm điểm; và GSEA trên toàn bộ gen xếp hạng theo tương quan với điểm nguy cơ. Dùng GO, KEGG, MSigDB Hallmark.

- Vi môi trường miễn dịch: ước lượng thành phần tế bào bằng nhiều phương pháp trong gói R immunedeconv (gọi qua %%R); so sánh nhóm nguy cơ cao và thấp bằng kiểm định Wilcoxon kèm FDR. Chỉ kết luận khi kết quả nhất quán giữa các phương pháp và giữa hai cohort. Kiểm tra mỗi cặp trong signature còn liên quan với kết cục sau khi hiệu chỉnh điểm miễn dịch và độ tinh khiết; cặp mất ý nghĩa được diễn giải là chỉ dấu thành phần mô, không phải trục điều hòa.

- Đáp ứng thuốc: gói R oncoPredict (gọi qua %%R) với dữ liệu GDSC2, so sánh IC50 dự đoán giữa hai nhóm nguy cơ; ưu tiên các nhóm thuốc dùng trong TNBC như taxane, platinum, ức chế PARP. Kết quả ghi rõ là giả thuyết.

- Tra cứu y văn và miRTarBase cho từng cặp trong signature để biết cặp nào đã có bằng chứng thực nghiệm, cặp nào là phát hiện mới.

**Đơn vị kiến thức cần áp dụng**

- Giá trị Shapley và giới hạn khi diễn giải (tầm quan trọng không đồng nghĩa với quan hệ nhân quả).

- Phân tích làm giàu: phân phối siêu bội trong ORA, nguyên lý GSEA, cấu trúc cơ sở dữ liệu GO, KEGG, Hallmark.

- Miễn dịch học khối u: tế bào lympho xâm nhập u (TILs), các phân nhóm miễn dịch của TNBC, nguyên lý tách thành phần tế bào (deconvolution) từ dữ liệu bulk.

- Dược lý ung thư: IC50, dữ liệu sàng lọc thuốc trên dòng tế bào, phác đồ điều trị TNBC hiện hành.

- Kiểm định phi tham số và đa kiểm định.

**Đầu ra và tiêu chí hoàn thành:**  bộ hình và bảng cho chương kết quả; danh sách trục miRNA–mRNA ưu tiên, mỗi trục kèm bằng chứng tính toán, bằng chứng y văn và mức độ nhất quán giữa hai cohort.

### 5.6. Tổng hợp đơn vị kiến thức

Thống kê sống còn và thiết kế đánh giá mô hình là hai mảng cần thành thạo nhất, vì chúng quyết định kết quả có đứng vững trước hội đồng và phản biện hay không. Nguồn học dưới đây là gợi ý, cần kiểm tra phiên bản mới nhất.

| Đơn vị kiến thức | Dùng ở | Mức cần đạt | Nguồn học gợi ý |
|---|---|---|---|
| Thống kê sống còn (kiểm duyệt, Kaplan–Meier, Cox, giả định tỷ lệ nguy cơ) | ND1, ND4 | Thành thạo | An Introduction to Statistical Learning with Applications in Python (ISLP), chương 11; tài liệu scikit-survival và lifelines |
| Thiết kế và đánh giá mô hình dự báo (nested CV, rò rỉ dữ liệu, EPV, hiệu chuẩn, DCA) | ND4 | Thành thạo | Steyerberg, Clinical Prediction Models; hướng dẫn TRIPOD+AI |
| Học máy có phạt và mô hình cây (elastic-net, random forest, boosting) | ND4 | Thành thạo | An Introduction to Statistical Learning, chương 6 và 8 |
| Kiểm định và đa kiểm định (Spearman, Wilcoxon, FDR) | ND2, ND5 | Thành thạo | Giáo trình thống kê sinh học của khoa |
| Sinh học miRNA (sinh tổng hợp, RISC, seed, 3'UTR) | ND2, ND5 | Hiểu sâu | Bartel, Cell 2018, "Metazoan MicroRNAs" |
| Ung thư vú và TNBC (IHC, HER2, PAM50, điều trị) | ND1, ND5 | Hiểu | Các bài tổng quan về TNBC trên Nature Reviews Clinical Oncology |
| Dữ liệu TCGA/GDC và METABRIC | ND1 | Thực hành | Tài liệu API của GDC; tài liệu cBioPortal |
| Tiền xử lý omics và hiệu ứng nền tảng | ND1, ND3 | Thực hành | Các workflow RNA-seq trên Bioconductor |
| Giảm chiều (PCA) | ND3 | Hiểu | An Introduction to Statistical Learning, chương 12 |
| Deep learning cơ bản (CNN, transformer, embedding, học chuyển giao) | ND2, ND3 | Hiểu khái niệm, chạy được suy luận | Khóa học miễn phí của Hugging Face |
| Diễn giải mô hình (SHAP) | ND5 | Hiểu | Molnar, Interpretable Machine Learning (bản trực tuyến miễn phí) |
| Phân tích làm giàu (ORA, GSEA) | ND5 | Thực hành | Tài liệu chính thức của gseapy |
| Miễn dịch khối u và deconvolution | ND5 | Hiểu | Sturm et al., Bioinformatics 2019 |
| Lập trình Python (pandas, scikit-learn), R cơ bản cho Nội dung 5, Git | Tất cả | Thực hành | Tài liệu chính thức của từng công cụ |

## 6. Dữ liệu, công cụ và phương pháp

Toàn bộ phân tích chạy bằng Python, trên máy cá nhân (Nội dung 1) hoặc trên Google Colab, Kaggle Notebooks (các bước cần GPU); riêng ước lượng vi môi trường miễn dịch và dự đoán đáp ứng thuốc dùng gói R, gọi qua ô %%R (thư viện rpy2) trong cùng notebook.

quy trình tổng thể · 2 cohort, 3 bộ đặc trưng, 5 nội dung

Đặc trưng cặp (khung tô màu) là giả thuyết chính; hai bộ còn lại là đối chứng, đi qua cùng cách chia dữ liệu và cùng các chỉ số đánh giá.

### 6.1. Thu thập và tiền xử lý dữ liệu

- **TCGA-BRCA.**  Tải qua API của GDC bằng Python: miRNA mức trưởng thành từ file Isoform Expression Quantification (đơn vị RPM), RNA-seq từ file STAR - Counts (đơn vị TPM), dữ liệu lâm sàng dạng BCR Biotab. Chỉ giữ mẫu khối u nguyên phát; ghép miRNA và mRNA theo mã mẫu (barcode mẫu), không chỉ mã bệnh nhân. Ghi phiên bản dữ liệu GDC (Data Release 46.0) và md5 của mọi file; lưu cả số read của miRNA.

- **Endpoint.**  PFI lấy từ bộ dữ liệu lâm sàng chuẩn hóa TCGA-CDR (Liu et al., *Cell*  2018), không dùng trường lâm sàng thô của GDC. Endpoint phụ: DSS và OS ở cả hai cohort (METABRIC: DSS theo "Died of Disease"); DFI khi cohort kiểm chứng chỉ có RFS hoặc DFS. Ở METABRIC, RFS (endpoint kiểm chứng chính) được cắt tại 5 năm cho khớp thời gian theo dõi của TCGA.

- **Lọc TNBC.**  Dùng các trường trạng thái ER, PR, HER2 đo bằng IHC (HER2 dùng ISH khi thiếu IHC, theo ASCO/CAP 2018) trong dữ liệu lâm sàng bổ sung của TCGA. Bệnh nhân thiếu bất kỳ trường nào bị loại; số ca loại được báo cáo. Ca có kết quả HER2 mâu thuẫn bị loại.

- **METABRIC.**  Dữ liệu mRNA và lâm sàng (gồm RFS) từ cBioPortal; dữ liệu miRNA từ nghiên cứu Dvinge et al. (*Nature*  2013). Dữ liệu miRNA có cho khoảng 1.283 mẫu; thử tải qua Synapse (tài khoản miễn phí) trong tháng 11/2026, không lấy được thì theo phương án ở mục 8. Bản công khai chỉ có IHC cho ER; PR và HER2 lấy từ trường trạng thái có sẵn (định nghĩa tổ hợp, 276 ca TNBC).

- **Lọc và biến đổi.**  Loại miRNA, gen biểu hiện thấp theo số read (miRNA ≥10 read; mRNA ≥10 count và TPM ≥1; ở ≥50% mẫu), chỉ giữ phân tử đo được trên cả hai nền tảng; chuyển về thang log2 trong từng cohort. METABRIC đã ở thang log2 nên không lấy log thêm.

### 6.2. Xây dựng đặc trưng cặp miRNA–mRNA

- **Chấm điểm cặp bằng miRBind2.**  Đầu vào là trình tự miRNA trưởng thành (miRBase) và trình tự 3'UTR của gen. Đầu ra là điểm dự đoán ức chế ở cấp gen. Chỉ chấm các miRNA và gen biểu hiện trong TNBC; giữ các cặp thuộc nhóm điểm cao nhất của mỗi miRNA (ngưỡng phần trăm cố định trước).

- **Lọc theo tương quan âm** **từng phần** **.**  Tính tương quan Spearman từng phần giữa miRNA và mRNA đích, hiệu chỉnh độ tinh khiết, điểm miễn dịch và mô đệm, trên toàn bộ u nguyên phát ER− của TCGA-BRCA. Giữ cặp có tương quan âm có ý nghĩa sau hiệu chỉnh FDR. Bước lọc không dùng dữ liệu sống còn nên không rò rỉ kết cục, và được làm một lần trước khi mô hình hóa; lọc bên trong từng vòng cross-validation trên cohort TNBC là phân tích độ nhạy.

- **Tính log-ratio trong từng mẫu** **(phân tích độ nhạy)** **.**  Với bệnh nhân *s* , cặp miRNA *i*  và mRNA *j* :

    f_{ij}(s) = \log_2\big(\text{miRNA}_i(s) + 1\big) - \log_2\big(\text{mRNA}_j(s) + 1\big)

**Đặc trưng chính**  là hiệu phân vị trong mẫu; P là phân vị của phân tử trong mẫu *s* , tính trên cùng tập miRNA và gen ở mọi cohort:

    f^{\mathrm{rank}}_{ij}(s) = P^{\mathrm{miR}}_i(s) - P^{\mathrm{mRNA}}_j(s)

**Dạng qPCR**  (chuẩn hóa theo phân tử tham chiếu trong cùng mẫu, tương đương ΔCt):

    f^{\mathrm{ref}}_{ij}(s) = \big[\log_2 \text{miR}_i(s) - \overline{\log_2 \text{miR}_{\mathrm{ref}}}(s)\big] - \big[\log_2 \text{mRNA}_j(s) - \overline{\log_2 \text{gen}_{\mathrm{ref}}}(s)\big]

Cả ba dạng đặc trưng tính được trên một bệnh nhân duy nhất, không cần chuẩn hóa theo cohort. Đây là tính chất cần cho một xét nghiệm lâm sàng.

Với mô hình tuyến tính, log-ratio tương đương việc đưa miRNA và mRNA thành hai đặc trưng có hệ số ràng buộc ngược dấu; lợi thế nếu có đến từ giảm số tham số, và được kiểm tra trực tiếp bằng đối chứng (i). Độ lệch cộng thêm theo nền tảng chỉ dịch điểm nguy cơ một hằng số; thứ làm hỏng khi chuyển nền tảng là độ dốc và ngưỡng tuyệt đối. Vì vậy đặc trưng chính dùng hiệu phân vị trong mẫu (cùng họ với relative expression ordering / k-TSP; Geman et al., 2004; Tan et al., 2005), bất biến với mọi biến đổi đơn điệu.

### 6.3. Đặc trưng so sánh

- **Embedding BulkRNABert.**  Mô hình transformer tiền huấn luyện trên RNA-seq bulk, có sẵn trên Hugging Face. Chạy suy luận (không huấn luyện) để lấy vector biểu diễn cho mỗi bệnh nhân. Dùng phiên bản tiền huấn luyện trên dữ liệu không ung thư (GTEx, ENCODE), không có TCGA, để tránh rò rỉ thông tin; cần xác minh trọng số của phiên bản này được công bố. Nếu không có, thay bằng embedding scFoundation, mô hình tiền huấn luyện trên dữ liệu đơn bào, không có TCGA.

- **Embedding BulkFormer.**  Mô hình 150 triệu tham số tiền huấn luyện trên hơn 500.000 mẫu bulk. Chạy suy luận tương tự.

- **Biểu hiện gen thô.**  Baseline không dùng mô hình nền tảng.

- Embedding có hàng trăm chiều; giảm chiều bằng PCA, khớp PCA chỉ trên phần dữ liệu huấn luyện của mỗi vòng.

- Với METABRIC (microarray), dữ liệu được ánh xạ tên gen và đưa về thang tương đương trước khi đưa vào mô hình nền tảng. Mức suy giảm hiệu năng khi chuyển nền tảng là một kết quả cần đo, không phải lỗi cần che.

### 6.4. Mô hình hóa

- **CoxNet**  (scikit-survival: CoxnetSurvivalAnalysis): mô hình Cox có phạt elastic-net, tự động loại đặc trưng không quan trọng. Mô hình tham chiếu chính. Phân tích chính: l1_ratio cố định 0,5; alpha chọn theo quy tắc 1-SE trên CV lặp ở vòng trong; tối đa 5 cặp trong signature (EPV khoảng 26/5 ≈ 5; Riley et al., 2020).

- **Random Survival Forest**  (scikit-survival: RandomSurvivalForest): rừng cây quyết định cho dữ liệu sống còn, bắt được quan hệ phi tuyến. Phân tích thăm dò (một cấu hình bảo thủ).

- **Gradient Boosting Survival**  (scikit-survival: GradientBoostingSurvivalAnalysis, hoặc xgboost với hàm mục tiêu Cox): cộng dồn nhiều cây nhỏ, mỗi cây sửa sai số của các cây trước. Phân tích thăm dò (một cấu hình bảo thủ).

- **Đánh giá nội bộ.**  Chỉ số chính ước lượng bằng bootstrap hiệu chỉnh lạc quan (200 lần, khớp lại toàn bộ pipeline mỗi lần); repeated stratified 5-fold CV (20 lần, phân tầng theo sự kiện) là phân tích phụ. Siêu tham số chọn bằng CV lặp ở vòng trong. Cài đặt bằng StratifiedKFold và GridSearchCV của scikit-learn; các bước học từ dữ liệu được đóng gói trong Pipeline. Mọi bước học từ dữ liệu (lọc cặp, PCA, chọn đặc trưng) nằm trong mỗi lần bootstrap và mỗi phần CV. Chọn đặc trưng bằng stability selection: lặp lại trên nhiều mẫu bootstrap, chỉ giữ cặp được chọn ổn định.

- CoxNet trên 4 bộ đặc trưng là phân tích chính; 12 tổ hợp 3 mô hình × 4 bộ đặc trưng là phân tích thăm dò, đánh giá trên cùng các mẫu bootstrap và cùng cách chia dữ liệu để so sánh công bằng. Mô hình cuối chọn theo quy tắc 1-SE; RSF và Gradient Boosting dùng không gian siêu tham số bảo thủ.

### 6.5. Đánh giá hiệu năng

- **C-index**  (Harrell và Uno; Uno cắt tại τ = 3 năm): khả năng xếp đúng thứ tự bệnh nhân nào tiến triển trước.

- **AUC phụ thuộc thời gian**  tại 2 và 3 năm (TCGA chỉ còn 31 ca được theo dõi đến 5 năm) (sksurv.metrics.cumulative_dynamic_auc).

- **Integrated Brier Score**  và **đường hiệu chuẩn**  (calibration): mức sai lệch giữa nguy cơ dự đoán và thực tế.

- **Decision curve analysis**  (tính lợi ích ròng theo Vickers và Elkin, 2006): lợi ích lâm sàng ròng so với chiến lược điều trị tất cả hoặc không ai.

- **Kaplan–Meier**  cho nhóm nguy cơ cao và thấp, ngưỡng cắt cố định từ cohort huấn luyện rồi áp nguyên cho METABRIC.

- **Giá trị cộng thêm:**  Cox đa biến kèm tuổi, giai đoạn, tình trạng hạch; so sánh với các signature đã công bố.

- **So sánh thống kê**  giữa các bộ đặc trưng bằng bootstrap chênh lệch C-index.

- **Đối chứng null:**  (i) cùng miRNA và gen đưa vào thành đặc trưng riêng; (ii) cặp ngẫu nhiên giữ miRNA, thay gen bằng gen không phải đích cùng thập phân vị biểu hiện (200 lần, p thực nghiệm); (iii) tuổi, giai đoạn, hạch + điểm miễn dịch + purity; (iv) (iii) + đặc trưng cặp.

- **Đánh giá nội bộ:**  bootstrap hiệu chỉnh lạc quan; ở METABRIC, theo dõi được cắt tại 5 năm (kiểm duyệt hành chính) để khớp với TCGA.

- Báo cáo theo hướng dẫn TRIPOD+AI.

### 6.6. Diễn giải sinh học

- **SHAP**  (shap, survshap): đóng góp của từng cặp vào điểm nguy cơ của từng bệnh nhân.

- **Làm giàu chức năng**  (gseapy, bộ gen GO, KEGG, MSigDB Hallmark): các mRNA đích trong signature tập trung vào con đường sinh học nào.

- **Vi môi trường miễn dịch**  (gói R immunedeconv gọi qua %%R, gồm các phương pháp như quanTIseq, xCell, MCP-counter): ước lượng thành phần tế bào miễn dịch, so sánh giữa nhóm nguy cơ.

- **Dự đoán đáp ứng thuốc**  (gói R oncoPredict gọi qua %%R, dữ liệu dòng tế bào GDSC2): ước tính nhóm nguy cơ cao nhạy với thuốc nào. Kết quả mang tính tạo giả thuyết.

### 6.7. Bảng tổng hợp công cụ

| Công cụ | Ngôn ngữ / nền tảng | Công dụng trong đề tài | Nội dung |
|---|---|---|---|
| API của GDC (requests) | Python | Tải dữ liệu TCGA | 1 |
| TCGA-CDR | Bảng dữ liệu | Endpoint PFI chuẩn hóa | 1 |
| cBioPortal | Gói dữ liệu tải về | Dữ liệu METABRIC | 1 |
| miRBase (mature.fa) | Cơ sở dữ liệu | Tên và trình tự miRNA trưởng thành | 1, 2 |
| PanCanAtlas (ABSOLUTE; Thorsson et al., 2018) | Bảng dữ liệu | Độ tinh khiết, thâm nhiễm miễn dịch | 1, 2, 5 |
| GENCODE v36 GTF (quét bằng pandas) | Python | Cờ gen chồng lấn trong STAR – Counts | 1 |
| UCSC Xena TOIL RSEM | Bảng dữ liệu | Độ nhạy định lượng mRNA | 1 |
| pandas, pyarrow | Python | Xử lý bảng, lưu Parquet | Tất cả |
| PyDESeq2 | Python | VST cho phân tích độ nhạy | 1 |
| scipy, statsmodels | Python | Tương quan Spearman, hiệu chỉnh FDR, tương quan từng phần | 1, 2 |
| miRBind2 | Python, Colab | Chấm điểm cặp miRNA–gen | 2 |
| BulkRNABert | Python, Hugging Face | Embedding mô hình nền tảng (thay thế: scFoundation nếu không có bản không TCGA) | 3 |
| BulkFormer | Python | Embedding mô hình nền tảng | 3 |
| scikit-learn | Python | Pipeline, nested CV, PCA | 3, 4 |
| scikit-survival | Python | CoxNet, RSF, Gradient Boosting Survival; C-index, AUC theo thời gian, Brier | 4 |
| lifelines | Python | Kaplan–Meier, log-rank, Cox đa biến, kiểm tra giả định tỷ lệ nguy cơ | 1, 4 |
| shap, survshap | Python | Giải thích mô hình | 5 |
| gseapy | Python | Làm giàu chức năng (ORA, GSEA) | 5 |
| immunedeconv | R, gọi qua %%R | Ước lượng tế bào miễn dịch | 5 |
| oncoPredict | R, gọi qua %%R | Dự đoán đáp ứng thuốc | 5 |

## 7. Kết quả dự kiến và sản phẩm

Khóa luận cho ra một câu trả lời định lượng cho câu hỏi nghiên cứu, dù đặc trưng cặp thắng hay thua embedding mô hình nền tảng; cả hai chiều đều là kết quả có giá trị công bố.

**Kết quả khoa học:**

- Bộ dữ liệu TNBC đã ghép và tiền xử lý từ TCGA-BRCA và METABRIC, kèm báo cáo số ca và số sự kiện.

- Danh sách cặp miRNA–mRNA ứng viên trong TNBC, có điểm miRBind2 và tương quan âm.

- Bảng so sánh hiệu năng 12 tổ hợp (3 mô hình × 4 bộ đặc trưng) (CoxNet là phân tích chính) đánh giá bằng bootstrap hiệu chỉnh lạc quan và trên METABRIC.

- Một signature tiên lượng cuối, gồm số ít cặp miRNA–mRNA, có công thức điểm nguy cơ cụ thể.

- Mô tả sinh học các trục điều hòa trong signature: con đường liên quan, khác biệt miễn dịch, giả thuyết về đáp ứng thuốc.

**Sản phẩm:**

- Khóa luận tốt nghiệp.

- Mã nguồn có thể tái lập (Python, kèm phần R cho Nội dung 5), lưu trên GitHub.

- Bộ kết quả đủ để viết bản thảo bài báo gửi tạp chí Q1–Q2 sau khi tốt nghiệp.

## 8. Rủi ro và phương án dự phòng

Rủi ro lớn nhất là cỡ mẫu TNBC nhỏ trong TCGA; nó quyết định mô hình chịu được bao nhiêu đặc trưng và phải được đo ngay tháng đầu.

| Rủi ro | Mức độ | Dấu hiệu nhận biết | Phương án dự phòng |
|---|---|---|---|
| Ít sự kiện PFI ở TNBC trong TCGA (ước tính 110–160 ca, vài chục sự kiện) | Cao | Dưới khoảng 40 sự kiện sau khi lọc. TRẠNG THÁI 11/2026: ĐÃ XẢY RA – 22 → 26 sự kiện PFI sau khi sửa định nghĩa HER2 (152 ca; 54 ca còn nguy cơ ở 3 năm, 31 ca ở 5 năm); phương án dự phòng đã kích hoạt 11/2026 | Stability selection: chỉ giữ cặp được chọn ở phần lớn các lần lấy mẫu bootstrap; báo cáo khoảng tin cậy bootstrap cho mọi chỉ số hiệu năng; bổ sung 11/2026: CoxNet là phân tích chính, bootstrap hiệu chỉnh lạc quan, AUC tại 2 và 3 năm, tối đa 5 cặp, kiểm chứng hai chiều nhánh mRNA, cắt METABRIC tại 5 năm |
| RSF và Gradient Boosting overfitting | Cao | Hiệu năng vòng trong cao, vòng ngoài thấp | Kết hợp bốn biện pháp: CoxNet là mô hình chính; RSF và Gradient Boosting dùng siêu tham số bảo thủ (nút lá lớn, cây nông, tốc độ học nhỏ, dừng sớm); chọn mô hình theo quy tắc 1-SE; báo cáo chênh lệch hiệu năng huấn luyện–kiểm tra |
| Dữ liệu miRNA METABRIC thuộc diện truy cập có kiểm soát | Cao | Không tải được từ nguồn công khai | Thử tải miRNA qua Synapse (tài khoản miễn phí) trước; luôn dùng phần công khai trên cBioPortal (mRNA, lâm sàng) để kiểm chứng embedding và biểu hiện gen thô. Không lấy được miRNA thì đặc trưng cặp kiểm chứng trên cohort GEO thay thế; không có bộ đơn lẻ phù hợp thì gộp nhiều cohort GEO nhỏ, lọc TNBC trong từng bộ và hiệu chỉnh khác biệt giữa các bộ. Không dùng bản sao dữ liệu không chính thức. Ứng viên GEO cụ thể: GSE22220/GSE22219 (Buffa et al., 2011; miRNA + mRNA, tái phát xa); SCAN-B GSE96058 (RNA-seq, chỉ kiểm chứng nhánh mRNA, cần kiểm tra endpoint có sẵn); cần kiểm tra số ca TNBC trước khi dùng |
| Định nghĩa TNBC ở METABRIC không hoàn toàn theo IHC | Trung bình | PR hoặc HER2 suy ra từ biểu hiện hoặc số bản sao | Chỉ giữ ca có đủ kết quả IHC cho cả ER, PR và HER2; báo cáo số ca bị loại. Nếu METABRIC thiếu IHC cho PR hoặc HER2, chạy song song hai định nghĩa: IHC chặt (nếu còn đủ ca) và tổ hợp có sẵn (ER theo IHC, HER2 theo số bản sao, PR theo trường có sẵn); báo cáo cả hai. Kết quả 11/2026: không có IHC cho PR và HER2; định nghĩa tổ hợp cho 276 ca; 18 ca ER bất đồng được gắn cờ |
| PFI (TCGA) và RFS (METABRIC) không tương đương | Trung bình | Định nghĩa sự kiện khác nhau | Thêm OS làm endpoint phụ ở cả hai cohort; khi cohort kiểm chứng dùng RFS hoặc DFS, huấn luyện thêm mô hình với DFI của TCGA-CDR để khớp định nghĩa (DFI có ít ca và ít sự kiện hơn PFI); bổ sung 11/2026: thêm DSS ở cả hai cohort, cắt RFS của METABRIC tại 5 năm, báo cáo loại sự kiện PFI (5/26 sự kiện là u nguyên phát mới) |
| Rò rỉ thông tin do BulkRNABert đã tiền huấn luyện trên TCGA | Trung bình | Hiệu năng trên TCGA cao bất thường so với METABRIC | Dùng phiên bản BulkRNABert tiền huấn luyện trên dữ liệu không ung thư (GTEx, ENCODE), không có TCGA; cần xác minh tác giả có công bố trọng số của phiên bản này; nếu không có, thay bằng embedding scFoundation (tiền huấn luyện trên dữ liệu đơn bào, không có TCGA) |
| Colab miễn phí bị ngắt phiên hoặc thiếu bộ nhớ | Trung bình | Lỗi hết bộ nhớ, phiên hết hạn | Chạy theo lô nhỏ, lưu kết quả trung gian lên Google Drive; dùng Kaggle Notebooks làm nguồn GPU miễn phí thứ hai; thuê GPU theo giờ cho bước trích embedding khi cần (chạy một lần, lưu lại) |
| Khối lượng công việc vượt thời gian (3 nhánh phân tích sinh học) | Trung bình | Trễ mốc tháng 5/2027 | Thứ tự ưu tiên cố định: SHAP và làm giàu chức năng bắt buộc; vi môi trường miễn dịch, rồi đáp ứng thuốc, làm nếu kịp tiến độ |
| Bộ lọc tương quan và signature phản ánh thành phần mô (TIL, độ tinh khiết) thay vì điều hòa | Cao | Số cặp dương có ý nghĩa ≥ số cặp âm; miRNA tương quan mạnh với điểm miễn dịch (11/2026: 196.113 cặp dương so với 165.509 cặp âm trên 12.511 gen; 9 miRNA có |ρ| ≥ 0.5 với điểm miễn dịch) | Tương quan từng phần; đối chứng (iii)–(iv); kiểm tra cặp sau hiệu chỉnh ở ND5 |
| miRNA-seq mang tính thành phần (miR-21-5p chiếm khoảng 24,8% RPM) | Trung bình | Đặc trưng thay đổi khi miR-21 thay đổi | Ngưỡng theo read; độ nhạy bỏ miR-21, upper-quartile, DESeq2; đặc trưng hiệu phân vị |
| Lỗi định lượng gen chồng lấn trong STAR – Counts | Trung bình | Gen biểu hiện cao nhưng TPM thấp bất thường (CD68, EIF4A1); đối chiếu RSEM (11/2026): 726 gen bị STAR đếm thiếu từ 4 lần trở lên | Loại 314 gen bị STAR định lượng sai (đối chiếu RSEM của TOIL) khỏi phân tích chính: ứng viên tạo cặp, tập gen tính phân vị f_rank và điểm thành phần ở cả hai cohort, baseline gen thô, ứng viên gen tham chiếu f_ref; giữ trong đầu vào ND3. Độ nhạy: chọn cặp lại bằng RSEM TPM (D9) |
| Hai cohort không đồng nhất (LAR, claudin-low, ER bất đồng, thời kỳ điều trị) | Trung bình | Phân bố PAM50, cờ ER, điều trị khác nhau (METABRIC: 32% claudin-low; 161/276 ca hóa trị) | Bảng 1; hiệu chỉnh hóa trị; phân tích độ nhạy loại ca bất đồng |
| Đặc trưng cặp không vượt embedding mô hình nền tảng | Không phải rủi ro | Kết quả ngược giả thuyết | Vẫn là kết quả có giá trị; báo cáo trung thực, phân tích nguyên nhân |

## 9. Kế hoạch thực hiện

Mốc kiểm tra cuối tháng 11/2026 quyết định phần còn lại của kế hoạch: số sự kiện PFI giới hạn số cặp tối đa trong signature, và quyền truy cập dữ liệu miRNA METABRIC quyết định cohort kiểm chứng.

**Kết quả mốc kiểm tra 11/2026**  (chốt 04/11/2026): 152 bệnh nhân TNBC, 26 sự kiện PFI, dưới ngưỡng 40, nên áp phương án dự phòng ở mục 8. Tám quyết định D1–D8 (đánh giá, kiểm chứng hai chiều, endpoint, bộ lọc tương quan, định nghĩa đặc trưng, HER2, ER-low/LAR, ngưỡng lọc) được ghi trong nhật ký dự án và tóm tắt ở docs/thay_doi_de_cuong_v2.md. Dữ liệu miRNA METABRIC chưa có.

kế hoạch thực hiện · 5 nội dung, 1 mốc kiểm tra, bảo vệ

Các nội dung chồng lấn có chủ đích: ND3 chạy song song cuối ND2 vì chỉ cần dữ liệu mRNA. ND4 chỉ bắt đầu khi ND2 có ma trận đặc trưng. Tết Nguyên đán (đầu tháng 2/2027) nằm trong ND2 và ND3, đã tính vào thời lượng.

## 10. Hướng phát triển

Bước tiếp theo sau khóa luận là kiểm chứng thực nghiệm các cặp miRNA–mRNA trong signature, chuyển kết quả tính toán thành một xét nghiệm đo được trong phòng thí nghiệm.

- **Định lượng miRNA bằng stem-loop RT-qPCR.**  Thiết kế mồi stem-loop RT, mồi xuôi đặc hiệu và mồi ngược phổ quát cho từng miRNA trong signature.

- **Định lượng mRNA đích bằng RT-qPCR.**  Đo trên cùng mẫu RNA để tính đặc trưng f^ref theo đúng định nghĩa của mô hình.

- **Mẫu thử nghiệm.**  Dòng tế bào TNBC hoặc mẫu mô u lưu trữ, tùy điều kiện phòng thí nghiệm và phê duyệt đạo đức.

- **Ưu thế của đặc trưng chuẩn hóa theo tham chiếu.**  qPCR không đo ra RPM hay TPM; đặc trưng chuyển giao là f^ref (chuẩn hóa theo miRNA và gen tham chiếu, tương đương ΔCt), nên vẫn đo được trên từng bệnh nhân mà không cần một cohort tham chiếu. Các phân tử tham chiếu được chọn ở ND2 theo độ ổn định, không dùng kết cục.

Hướng này phụ thuộc vào điều kiện phòng thí nghiệm, mẫu và phê duyệt đạo đức; khóa luận không cam kết thực hiện.

## Tài liệu tham khảo

Tài liệu có liên kết đã được mở và kiểm tra; tài liệu không có liên kết được trích từ hiểu biết sẵn có, cần đối chiếu số trang và năm trước khi nộp.

**Dữ liệu và endpoint**

- Liu J, et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. *Cell* . 2018;173(2):400–416.

- Curtis C, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. *Nature* . 2012;486:346–352.

- Dvinge H, et al. The shaping and functional consequences of the microRNA landscape in breast cancer. *Nature* . 2013;497:378–382.

- Pölsterl S. scikit-survival: a library for time-to-event analysis built on top of scikit-learn. *Journal of Machine Learning Research* . 2020;21(212):1–6.

- Wolff AC, et al. Human epidermal growth factor receptor 2 testing in breast cancer: American Society of Clinical Oncology/College of American Pathologists clinical practice guideline focused update. *Journal of Clinical Oncology* . 2018;36(20):2105–2122.

- Allison KH, et al. Estrogen and progesterone receptor testing in breast cancer: ASCO/CAP guideline update. *Journal of Clinical Oncology* . 2020;38(12):1346–1366.

- Lehmann BD, et al. Refinement of triple-negative breast cancer molecular subtypes: implications for neoadjuvant chemotherapy selection. *PLoS One* . 2016;11(6):e0157368.

- Thorsson V, et al. The immune landscape of cancer. *Immunity* . 2018;48(4):812–830.

- Carter SL, et al. Absolute quantification of somatic DNA alterations in human cancer. *Nature Biotechnology* . 2012;30(5):413–421.

- Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. *Nature Communications* . 2015;6:8971.

- Vivian J, et al. Toil enables reproducible, open source, big biomedical data analyses. *Nature Biotechnology* . 2017;35(4):314–316.

- Buffa FM, et al. microRNA-associated progression pathways and potential therapeutic targets identified by integrated mRNA and microRNA expression profiling in breast cancer. *Cancer Research* . 2011;71(17):5635–5645.

**miRNA và dự đoán đích**

- Friedman RC, et al. Most mammalian mRNAs are conserved targets of microRNAs. *Genome Research* . 2009;19:92–105.

- Kozomara A, Birgaoanu M, Griffiths-Jones S. miRBase: from microRNA sequences to function. *Nucleic Acids Research* . 2019;47:D155–D162.

- miRBind2 enables sequence-only prediction of miRNA binding and transcript repression. *Bioinformatics* . 2026;42:btag410.

- Sammut S, et al. miRBench: novel benchmark datasets for microRNA binding site prediction that mitigate against prevalent microRNA frequency class bias. *Bioinformatics* . 2025;41:i542–i551.

- Lin W-H, et al. DeepExoMir: a reproducible RNA language model framework for CLIP-seq-supported microRNA target-site prioritization. *International Journal of Molecular Sciences* . 2026;27:6184.

- Telonis AG, et al. Beyond the one-locus-one-miRNA paradigm: microRNA isoforms enable deeper insights into breast cancer heterogeneity. *Nucleic Acids Research* . 2015;43(19):9158–9175.

**Tương tác miRNA–mRNA và tiên lượng**

- spongEffects: ceRNA modules offer patient-specific insights into the miRNA regulatory landscape. *Bioinformatics* . 2023.

- Bayesian inference of sample-specific coexpression networks. *Genome Research* . 2024;34:1397–1410.

- Unveiling novel miRNA–mRNA interactions and their prognostic roles in triple-negative breast cancer. *International Journal of Molecular Sciences* . 2025;26(5):1916.

- Constructing lncRNA-miRNA-mRNA networks specific to individual cancer patients and finding prognostic biomarkers.

**Mô hình nền tảng**

- Gélard M, Richard G, Pierrot T, Cournède P-H. BulkRNABert: cancer prognosis from bulk RNA-seq based language models.

- A large-scale foundation model for bulk transcriptomes (BulkFormer). *bioRxiv* .

- Ahlmann-Eltze C, Huber W, Anders S. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nature Methods* . 2025;22(8):1657–1661.

**Phương pháp và báo cáo**

- Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. *Journal of Statistical Software* . 2010;33(1):1–22.

- Ishwaran H, et al. Random survival forests. *Annals of Applied Statistics* . 2008;2(3):841–860.

- Lundberg SM, Lee S-I. A unified approach to interpreting model predictions. *NeurIPS* . 2017.

- Fang Z, Liu X, Peltz G. GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics* . 2023;39(1):btac757.

- Sturm G, et al. Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. *Bioinformatics* . 2019;35(14):i436–i445.

- Maeser D, Gruener RF, Huang RS. oncoPredict: an R package for predicting in vivo or cancer patient drug response and biomarkers from cell line screening data. *Briefings in Bioinformatics* . 2021;22(6):bbab260.

- Collins GS, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ* . 2024;385:e078378.

- Davidson-Pilon C. lifelines: survival analysis in Python. *Journal of Open Source Software* . 2019;4(40):1317.

- Muzellec B, et al. PyDESeq2: a python package for bulk RNA-seq differential expression analysis. *Bioinformatics* . 2023;39(9):btad547.

- Vickers AJ, Elkin EB. Decision curve analysis: a novel method for evaluating prediction models. *Medical Decision Making* . 2006;26(6):565–574.

- Geman D, et al. Classifying gene expression profiles from pairwise mRNA comparisons. *Statistical Applications in Genetics and Molecular Biology* . 2004;3(1):Article 19.

- Tan AC, et al. Simple decision rules for classifying human cancers from gene expression profiles. *Bioinformatics* . 2005;21(20):3896–3904.

- Uno H, et al. On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data. *Statistics in Medicine* . 2011;30(10):1105–1117.

- Riley RD, et al. Calculating the sample size required for developing a clinical prediction model. *BMJ* . 2020;368:m441.

- Gloor GB, et al. Microbiome datasets are compositional: and this is not optional. *Frontiers in Microbiology* . 2017;8:2224.

