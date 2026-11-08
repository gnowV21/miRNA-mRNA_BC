# Thay đổi đề cương v2 (mốc kiểm tra 11/2026)

- **File:** `docs/proposal_miRNA-mRNA_v2_theo_doi_thay_doi.docx` (tracked changes, tác giả "Claude Code"). Bản gốc `docs/proposal_miRNA-mRNA.docx` không bị sửa.
- **Số thao tác:** 117 (đợt 04/11/2026) + 27 (đợt 08/11/2026: D5, D10, D11, thứ bậc phân tích); tổng 292 phần tử `w:ins`, 39 phần tử `w:del`. Đợt 08/11 đánh số lại `w:id` của mọi phần tử sửa cho duy nhất (bản 04/11 có 9 id trùng do sao chép hàng bảng).
- **Bản đọc soát (đã chấp nhận thay đổi):** `docs/proposal_miRNA-mRNA_v2_sach.md`. Máy chưa có pandoc nên bản này được xuất bằng script Python (bỏ phần xóa, giữ phần chèn). Chưa xuất PDF vì máy chưa có LibreOffice.
- **Số liệu** lấy tự động từ `qc/` (notebook ND1 chạy ngày 04/11/2026), không chép tay.
- **Mã quyết định** D1–D11: xem mục "Mốc kiểm tra 11/2026 – quyết định" (D1–D9) và "Quyết định 07–08/11/2026" (D5 bổ sung, D10, D11, thứ bậc) trong `nhat_ky.md`.

## Khác với yêu cầu

- **Bảng rủi ro (mục 8):** không thêm cột "trạng thái" (Word theo dõi thay đổi kiểu thêm cột rất kém). Trạng thái 11/2026 được chèn vào ô "Dấu hiệu nhận biết" của từng dòng.
- **Công thức (mục 6.2):** công thức gốc là đoạn văn bản LaTeX kiểu `Code`, không phải equation OMML. Hai công thức mới chèn thành đoạn `Code` cùng kiểu; công thức gốc giữ nguyên.
- **Bảng công cụ (mục 6.7):** dòng GENCODE ghi "quét bằng pandas" vì `pyranges` không cài được trên Windows.
- **Tài liệu tham khảo:** 14 tài liệu mới, metadata đối chiếu qua DOI (doi.org) ngày 04/11/2026. Danh mục của đề cương không gắn liên kết và tự đánh số theo danh sách Word, trích dẫn trong thân bài dạng tác giả–năm, nên không cần đánh số lại thủ công.

## Lưu ý

Đề cương gốc là bản xuất từ Claude Doc. Sửa file Word **không** cập nhật Claude Doc. Sau khi duyệt, cần chép thay đổi lên Claude Doc, hoặc dùng bản Word làm bản chính từ đây.

## Danh sách thay đổi

| Mục | Nội dung cũ | Nội dung mới | Lý do | Mã |
|---|---|---|---|---|
| 1. Câu hỏi nghiên cứu | là hiệu số log giữa | là hiệu số (theo phân vị hoặc theo log) giữa | Đồng bộ với đặc trưng chính | D5 |
| 1. Câu hỏi nghiên cứu | (sau «tính ngay trên mẫu của từng bệnh nhân.») |  Đặc trưng chính là hiệu phân vị của miRNA và mRNA trong cùng mẫu; hiệu số log và hiệu số log chuẩn hóa theo phân tử tham chiếu (dạng qPCR) được báo cáo làm phân tích độ nhạy. | Chỉ hiệu phân vị bất biến khi chuyển RNA-seq sang microarray | D5 |
| 1. Câu hỏi nghiên cứu | (sau «đo bằng C-index») |  của Uno (cắt tại 3 năm) sau hiệu chỉnh lạc quan bằng bootstrap, và so với các mô hình đối chứng null (mục 6.5) | Ít sự kiện; cần đối chứng null | D1 |
| 1. Câu hỏi nghiên cứu | (sau «giữ được hiệu năng tốt hơn») |  khi các đặc trưng được định nghĩa bất biến với thang đo của nền tảng (hiệu phân vị trong mẫu) | Bất biến thang đo | D5 |
| 1. Câu hỏi nghiên cứu | (sau «cơ chế sinh học nào của TNBC») |  sau khi tách ảnh hưởng của thành phần mô (độ tinh khiết khối u, thâm nhiễm miễn dịch, mô đệm) | Tương quan thô chủ yếu phản ánh thành phần mô | D4 |
| 2.4 Khoảng trống | (sau «không cần chuẩn hóa cấp cohort») |  và bất biến với thang đo của nền tảng | Bất biến thang đo | D5 |
| 3. Mục tiêu cụ thể | (sau «có tương quan âm») |  từng phần (sau hiệu chỉnh thành phần mô) | Hiệu chỉnh thành phần mô | D4 |
| 3. Mục tiêu cụ thể | trong mô u TNBC | trong mô u ER− | Lọc không giám sát trên toàn bộ u ER− của TCGA-BRCA | D4 |
| 3. Mục tiêu cụ thể | tính đặc trưng log-ratio trên từng mẫu | tính đặc trưng cặp trong từng mẫu (hiệu phân vị; log-ratio làm độ nhạy) | Đặc trưng chính | D5 |
| 3. Mục tiêu cụ thể | và so sánh ba mô hình sống còn (CoxNet, Random Survival Forest, Gradient Boosting Survival) | mô hình CoxNet (phân tích chính) | Ít sự kiện: CoxNet chính | D1 |
| 3. Mục tiêu cụ thể | (sau «biểu hiện gen thô») | ; Random Survival Forest và Gradient Boosting Survival là phân tích thăm dò | Ít sự kiện | D1 |
| 3. Mục tiêu cụ thể | (sau «các signature đã công bố») | ; so sánh với các mô hình đối chứng null (cặp ngẫu nhiên, đặc trưng riêng lẻ, lâm sàng + miễn dịch); với nhánh chỉ dùng mRNA, kiểm chứng thêm theo chiều ngược METABRIC → TCGA | Đối chứng null; METABRIC giàu sự kiện hơn | D1, D2 |
| 4. Đối tượng và phạm vi | (sau «HER2 âm tính») |  (HER2 theo IHC hoặc ISH, ASCO/CAP 2018) | Hướng dẫn ASCO/CAP 2018 | D6 |
| 4. Đối tượng và phạm vi, bảng – Endpoint | (sau «PFI theo TCGA-CDR») |  (chính); DSS, OS (phụ) | Thêm DSS | D3 |
| 4. Đối tượng và phạm vi, bảng – Endpoint | RFS (xấp xỉ PFI) | RFS cắt tại 5 năm (chính); DSS (Died of Disease), OS (phụ) | TCGA theo dõi ngắn; TNBC tái phát sớm | D3 |
| 4. Đối tượng và phạm vi, bảng – Định nghĩa TNBC | (sau «IHC từ dữ liệu lâm sàng») | : ER−, PR− theo IHC; HER2 âm theo ASCO/CAP 2018 (Wolff et al., 2018): IHC 0/1+, hoặc IHC 2+ kèm ISH âm, hoặc không có IHC nhưng ISH âm. Loại ca có kết quả HER2 mâu thuẫn. TCGA ghi ER/PR "<10%" nên không tách được nhóm ER-low (1–9%; Allison et al., 2020); nhóm này được gắn cờ và phân tích độ nhạy | Đúng hướng dẫn; ER-low là giới hạn | D6, D7 |
| 4. Đối tượng và phạm vi, bảng – Định nghĩa TNBC | (cuối đoạn) |  Kết quả xác minh (11/2026): bản công khai không có IHC cho PR và HER2, nên định nghĩa IHC chặt cho 0 ca; dùng định nghĩa tổ hợp (276 ca TNBC; 18 ca ER_IHC âm nhưng ER_STATUS dương được gắn cờ). | Đã xác minh ở Bước 9 | D7 |
| 4. Đối tượng và phạm vi – Giới hạn | (cuối đoạn) |  Dữ liệu METABRIC thuộc thời kỳ điều trị cũ (161/276 ca được hóa trị, 64 ca ER− có điều trị nội tiết), nên phân tích kiểm chứng được hiệu chỉnh theo hóa trị. | Khác biệt điều trị giữa hai cohort | D3 |
| 5. Nội dung nghiên cứu (tóm tắt) | (sau «ghép dữ liệu miRNA, mRNA và lâm sàng») | , chú thích sinh học (PAM50, độ tinh khiết, thâm nhiễm miễn dịch) | Hiệp biến cần cho D4 | D4 |
| 5. Nội dung nghiên cứu (tóm tắt) | (sau «lọc theo tương quan âm») |  từng phần | Hiệu chỉnh thành phần mô | D4 |
| 5. Nội dung nghiên cứu (tóm tắt) | tính log-ratio trong từng mẫu | tính đặc trưng cặp trong từng mẫu (hiệu phân vị; log-ratio làm độ nhạy) | Đặc trưng chính | D5 |
| 5. Nội dung nghiên cứu (tóm tắt) | Huấn luyện ba mô hình trên bốn bộ đặc trưng bằng nested cross-validation | Huấn luyện CoxNet (chính; RSF và GBS thăm dò) trên bốn bộ đặc trưng, đánh giá bằng bootstrap hiệu chỉnh lạc quan | Ít sự kiện | D1 |
| 5.1 Nội dung 1 – bước 1 | (cuối đoạn) |  Ghi phiên bản dữ liệu GDC (data release) và md5 của mọi file. | Tái lập |  |
| 5.1 Nội dung 1 – bước 3 | HER2 IHC 2+ (không xác định) dùng kết quả FISH nếu có; không có FISH thì loại. | HER2 âm theo ASCO/CAP 2018 (Wolff et al., 2018): IHC 0/1+, hoặc IHC 2+ kèm ISH âm, hoặc không có IHC nhưng ISH âm; loại ca có kết quả HER2 mâu thuẫn (điểm IHC trái với trạng thái ghi nhận). Gắn cờ, không loại, các ca ER/PR "<10%", ESR1 cao, kiểu LAR. | Đúng hướng dẫn; thêm 17 ca chỉ có ISH | D6, D7 |
| 5.1 Nội dung 1 – bước 4 | lọc miRNA và gen biểu hiện thấp theo ngưỡng định trước; chuyển log2(x + 1). | lọc theo số read (miRNA ≥10 read; mRNA ≥10 count và TPM ≥1; ở ≥50% mẫu); chỉ giữ miRNA và gen đo được trên cả hai nền tảng; chuyển log2(x + 1) cho dữ liệu đếm (METABRIC đã ở thang log2). | RPM = 1 chỉ tương ứng 0,8–2,3 read | D8 |
| 5.1 Nội dung 1 – bước mới | (đoạn mới) | Chú thích sinh học và QC (không dùng sống còn): PAM50, độ tinh khiết ABSOLUTE (Carter et al., 2012), thâm nhiễm miễn dịch (Thorsson et al., 2018), điểm thành phần mô trong từng mẫu, cờ LAR và ESR1 cao, thuốc điều trị; cờ gen chồng lấn trong STAR – Counts và đối chiếu với RSEM của TOIL (Vivian et al., 2017); phân tích độ nhạy chuẩn hóa miRNA (bỏ miR-21-5p, upper-quartile, DESeq2). | Hiệp biến thành phần mô; lỗi định lượng gen chồng lấn; miRNA-seq mang tính thành phần | D4 |
| 5.1 Nội dung 1 – bước 6 | (sau «(tính bằng phương pháp Kaplan–Meier đảo)») | ; số ca còn nguy cơ tại 1, 2, 3, 5 năm; thời điểm sự kiện; loại sự kiện PFI | Theo dõi ngắn; PFI gồm cả u nguyên phát mới | D3 |
| 5.1 Nội dung 1 – đơn vị kiến thức | (đoạn mới) | Dữ liệu dạng thành phần (compositional) của miRNA-seq (Gloor et al., 2017); giới hạn của định lượng gen chồng lấn khi đếm read theo gen. | Phát hiện từ dữ liệu |  |
| 5.2 Nội dung 2 – mục tiêu | (sau «(tương quan âm») |  từng phần, sau hiệu chỉnh thành phần mô, | Hiệu chỉnh thành phần mô | D4 |
| 5.2 Nội dung 2 – mục tiêu | ở TNBC | ở mô u ER− | Lọc trên toàn bộ u ER− | D4 |
| 5.2 Nội dung 2 – mục tiêu | hàm tính đặc trưng log-ratio | hàm tính đặc trưng cặp (hiệu phân vị, log-ratio, log-ratio chuẩn hóa tham chiếu) | Ba dạng đặc trưng | D5 |
| 5.2 Nội dung 2 – bước 5 | Bên trong mỗi vòng cross-validation (không phải trên toàn bộ dữ liệu): tính tương quan Spearman miRNA–mRNA, hiệu chỉnh đa kiểm định Benjamini–Hochberg, giữ cặp có hệ số âm và FDR dưới ngưỡng. | Tính tương quan Spearman từng phần miRNA–mRNA, hiệu chỉnh theo độ tinh khiết, điểm miễn dịch và mô đệm, trên toàn bộ u nguyên phát ER− của TCGA-BRCA (khoảng 232 bệnh nhân có cả miRNA và mRNA); hiệu chỉnh đa kiểm định Benjamini–Hochberg, giữ cặp có hệ số âm và FDR dưới ngưỡng. Bước này không dùng dữ liệu sống còn nên không rò rỉ kết cục, và được làm một lần trước khi mô hình hóa. Lọc trong từng vòng CV trên cohort TNBC là phân tích độ nhạy. | Tương quan thô chủ yếu phản ánh thành phần mô (165.509 → 57.566 cặp âm sau hiệu chỉnh, trên 12.511 gen) | D4 |
| 5.2 Nội dung 2 – bước mới | (đoạn mới) | Kiểm tra miRBind2: đích dự đoán có lệch về phía tương quan âm so với gen không phải đích cùng mức biểu hiện không (kiểm định KS một phía). | QC miRBind2 trong bệnh cảnh | D4 |
| 5.2 Nội dung 2 – bước 6 | Tính log-ratio cho | Tính ba dạng đặc trưng cho | Ba dạng đặc trưng | D5 |
| 5.2 Nội dung 2 – bước 6 | (sau «mục 6.2») | : hiệu phân vị (chính), log-ratio, log-ratio chuẩn hóa theo phân tử tham chiếu (chọn không dùng kết cục) | Ba dạng đặc trưng | D5 |
| 5.2 Nội dung 2 – mở rộng | (đoạn mới) | Mở rộng nếu kịp: isomiR đầu 5′ (seed dịch chuyển nên tập đích khác; miRBind2 chấm trực tiếp trình tự isomiR; Telonis et al., 2015), tỷ lệ 5p/3p. | Dữ liệu isoform đã có sẵn |  |
| 5.3 Nội dung 3 – bước mới | (đoạn mới) | Nhánh chỉ dùng mRNA (embedding, gen thô) được kiểm chứng hai chiều: TCGA → METABRIC và METABRIC → TCGA. | METABRIC có 97 sự kiện RFS trong 5 năm, TCGA có 26 sự kiện PFI | D2 |
| 5.4 Nội dung 4 – mục tiêu | so sánh công bằng 12 tổ hợp mô hình × đặc trưng | so sánh CoxNet trên 4 bộ đặc trưng (phân tích chính; 12 tổ hợp mô hình × đặc trưng là phân tích thăm dò) | Ít sự kiện | D1 |
| 5.4 Nội dung 4 – bước 2 | Thiết lập nested cross-validation: vòng ngoài 5 phần, lặp lại nhiều lần với hạt giống ngẫu nhiên cố định, phân tầng theo trạng thái sự kiện; vòng trong 5 phần để chọn siêu tham số. | Đánh giá nội bộ chính: bootstrap hiệu chỉnh lạc quan (200 lần, gồm toàn bộ pipeline); repeated stratified 5-fold CV (20 lần) làm phân tích phụ; chọn siêu tham số bằng CV lặp ở vòng trong theo quy tắc 1-SE. | Mỗi fold chỉ có 4–5 sự kiện (26 sự kiện PFI) | D1 |
| 5.4 Nội dung 4 – bước 3 | mỗi phần của vòng ngoài | mỗi lần bootstrap và mỗi phần của CV | Đồng bộ D1 | D1 |
| 5.4 Nội dung 4 – bước 5 | C-index (Harrell, Uno) | C-index (Harrell; Uno cắt tại τ = 3 năm, Uno et al., 2011) | Theo dõi ngắn | D1 |
| 5.4 Nội dung 4 – bước 5 | tại 3 và 5 năm | tại 2 và 3 năm | Chỉ 31 ca còn nguy cơ ở 5 năm | D1 |
| 5.4 Nội dung 4 – bước 5 | (sau «Integrated Brier Score») |  (đến 3 năm ở TCGA, 5 năm ở METABRIC) | Theo dõi ngắn | D1, D3 |
| 5.4 Nội dung 4 – bước 5 | (cuối đoạn) |  Nếu vi phạm, báo cáo HR theo khoảng thời gian (0–2, 2–5 năm). | TNBC tái phát sớm, nguy cơ không tỷ lệ | D3 |
| 5.4 Nội dung 4 – bước mới | (đoạn mới) | Đối chứng null: (i) đặc trưng riêng lẻ, (ii) cặp ngẫu nhiên cùng mức biểu hiện (200 lần, p thực nghiệm), (iii) lâm sàng + miễn dịch + purity, (iv) (iii) + đặc trưng cặp. | Chứng minh quan hệ đích mang thông tin; tách khỏi thành phần mô | D1, D4 |
| 5.4 Nội dung 4 – bước 8 | (sau «lên METABRIC») |  với RFS cắt tại 5 năm | Khớp thời gian theo dõi | D3 |
| 5.4 Nội dung 4 – bước 8 | Cox đa biến kèm yếu tố lâm sàng | Cox đa biến kèm hóa trị, giai đoạn, hạch, điểm miễn dịch | Khác biệt điều trị; thành phần mô | D3, D4 |
| 5.4 Nội dung 4 – bước 8 | (cuối đoạn) |  Phân tích độ nhạy: loại ca ER bất đồng; chỉ giữ ca có hóa trị. | Hai cohort không đồng nhất | D7 |
| 5.4 Nội dung 4 – đầu ra | bảng so sánh 12 tổ hợp | bảng so sánh CoxNet trên 4 bộ đặc trưng kèm đối chứng null (12 tổ hợp là thăm dò) | Đồng bộ D1 | D1 |
| 5.5 Nội dung 5 – bước 3 | (cuối đoạn) |  Kiểm tra mỗi cặp trong signature còn liên quan với kết cục sau khi hiệu chỉnh điểm miễn dịch và độ tinh khiết; cặp mất ý nghĩa được diễn giải là chỉ dấu thành phần mô, không phải trục điều hòa. | 9 miRNA có \|ρ\| ≥ 0,5 với điểm miễn dịch | D4 |
| 6.1 Thu thập và tiền xử lý – TCGA-BRCA | (cuối đoạn) |  Ghi phiên bản dữ liệu GDC (Data Release 46.0) và md5 của mọi file; lưu cả số read của miRNA. | Tái lập; ngưỡng theo read | D8 |
| 6.1 Thu thập và tiền xử lý – Endpoint | (sau «Endpoint phụ: ») | DSS và  | Thêm DSS | D3 |
| 6.1 Thu thập và tiền xử lý – Endpoint | (sau «ở cả hai cohort») |  (METABRIC: DSS theo "Died of Disease") | Thêm DSS | D3 |
| 6.1 Thu thập và tiền xử lý – Endpoint | (cuối đoạn) |  Ở METABRIC, RFS (endpoint kiểm chứng chính) được cắt tại 5 năm cho khớp thời gian theo dõi của TCGA. | Theo dõi ngắn ở TCGA | D3 |
| 6.1 Thu thập và tiền xử lý – Lọc TNBC | (sau «HER2 đo bằng IHC») |  (HER2 dùng ISH khi thiếu IHC, theo ASCO/CAP 2018) | Đúng hướng dẫn | D6 |
| 6.1 Thu thập và tiền xử lý – Lọc TNBC | (cuối đoạn) |  Ca có kết quả HER2 mâu thuẫn bị loại. | Dữ liệu mâu thuẫn (vd. TCGA-LL-A73Y) | D6 |
| 6.1 Thu thập và tiền xử lý – METABRIC | Chỉ giữ ca có đủ IHC cho ER, PR, HER2. | Bản công khai chỉ có IHC cho ER; PR và HER2 lấy từ trường trạng thái có sẵn (định nghĩa tổ hợp, 276 ca TNBC). | Đã xác minh ở Bước 9 | D7 |
| 6.1 Thu thập và tiền xử lý – Lọc và biến đổi | (sau «gen biểu hiện thấp») |  theo số read (miRNA ≥10 read; mRNA ≥10 count và TPM ≥1; ở ≥50% mẫu), chỉ giữ phân tử đo được trên cả hai nền tảng | RPM = 1 là nhiễu đếm | D8 |
| 6.1 Thu thập và tiền xử lý – Lọc và biến đổi | (cuối đoạn) |  METABRIC đã ở thang log2 nên không lấy log thêm. | Thang đo microarray | D8 |
| 6.2 Đặc trưng cặp – lọc | (sau «Lọc theo tương quan âm») |  từng phần | Hiệu chỉnh thành phần mô | D4 |
| 6.2 Đặc trưng cặp – lọc | (sau «Tính tương quan Spearman») |  từng phần | Hiệu chỉnh thành phần mô | D4 |
| 6.2 Đặc trưng cặp – lọc |  trong mô u. | , hiệu chỉnh độ tinh khiết, điểm miễn dịch và mô đệm, trên toàn bộ u nguyên phát ER− của TCGA-BRCA. | Lọc không giám sát | D4 |
| 6.2 Đặc trưng cặp – lọc | Bước lọc này chạy  |  | Đổi chỗ lọc | D4 |
| 6.2 Đặc trưng cặp – lọc | bên trong |  | Đổi chỗ lọc | D4 |
| 6.2 Đặc trưng cặp – lọc |  từng vòng cross-validation, không chạy trên toàn bộ dữ liệu, để tránh rò rỉ thông tin. | Bước lọc không dùng dữ liệu sống còn nên không rò rỉ kết cục, và được làm một lần trước khi mô hình hóa; lọc bên trong từng vòng cross-validation trên cohort TNBC là phân tích độ nhạy. | Lọc không giám sát ổn định hơn 108 mẫu/vòng | D4 |
| 6.2 Đặc trưng cặp – công thức | (sau «Tính log-ratio trong từng mẫu») |  (phân tích độ nhạy) | f_raw thành độ nhạy | D5 |
| 6.2 Đặc trưng cặp – công thức mới | (đoạn mới) | Đặc trưng chính là hiệu phân vị trong mẫu; P là phân vị của phân tử trong mẫu s, tính trên cùng tập miRNA và gen ở mọi cohort: | Bất biến thang đo | D5 |
| 6.2 Đặc trưng cặp – công thức mới | (đoạn mới) | f^{\mathrm{rank}}_{ij}(s) = P^{\mathrm{miR}}_i(s) - P^{\mathrm{mRNA}}_j(s) | f_rank | D5 |
| 6.2 Đặc trưng cặp – công thức mới | (đoạn mới) | Dạng qPCR (chuẩn hóa theo phân tử tham chiếu trong cùng mẫu, tương đương ΔCt): | Chuyển giao sang qPCR | D5 |
| 6.2 Đặc trưng cặp – công thức mới | (đoạn mới) | f^{\mathrm{ref}}_{ij}(s) = \big[\log_2 \text{miR}_i(s) - \overline{\log_2 \text{miR}_{\mathrm{ref}}}(s)\big] - \big[\log_2 \text{mRNA}_j(s) - \overline{\log_2 \text{gen}_{\mathrm{ref}}}(s)\big] | f_ref | D5 |
| 6.2 Đặc trưng cặp – lập luận | Đặc trưng này | Cả ba dạng đặc trưng | Áp cho cả ba công thức | D5 |
| 6.2 Đặc trưng cặp – lập luận | , và là giả thuyết chính về độ bền xuyên nền tảng. | . | Lập luận cũ chưa đúng | D5 |
| 6.2 Đặc trưng cặp – lập luận | (đoạn mới) | Với mô hình tuyến tính, log-ratio tương đương việc đưa miRNA và mRNA thành hai đặc trưng có hệ số ràng buộc ngược dấu; lợi thế nếu có đến từ giảm số tham số, và được kiểm tra trực tiếp bằng đối chứng (i). Độ lệch cộng thêm theo nền tảng chỉ dịch điểm nguy cơ một hằng số; thứ làm hỏng khi chuyển nền tảng là độ dốc và ngưỡng tuyệt đối. Vì vậy đặc trưng chính dùng hiệu phân vị trong mẫu (cùng họ với relative expression ordering / k-TSP; Geman et al., 2004; Tan et al., 2005), bất biến với mọi biến đổi đơn điệu. | Lập luận bất biến thang đo | D5 |
| 6.4 Mô hình hóa | (cuối đoạn) |  Phân tích chính: l1_ratio cố định 0,5; alpha chọn theo quy tắc 1-SE trên CV lặp ở vòng trong; tối đa 5 cặp trong signature (EPV khoảng 26/5 ≈ 5; Riley et al., 2020). | Ít sự kiện | D1 |
| 6.4 Mô hình hóa | (cuối đoạn) |  Phân tích thăm dò (một cấu hình bảo thủ). | Ít sự kiện | D1 |
| 6.4 Mô hình hóa | (cuối đoạn) |  Phân tích thăm dò (một cấu hình bảo thủ). | Ít sự kiện | D1 |
| 6.4 Mô hình hóa | Nested cross-validation. | Đánh giá nội bộ. | Đồng bộ D1 | D1 |
| 6.4 Mô hình hóa |  Vòng ngoài 5 phần, lặp lại nhiều lần, phân tầng theo trạng thái sự kiện, dùng để ước lượng hiệu năng. Vòng trong 5 phần dùng để chọn siêu tham số. |  Chỉ số chính ước lượng bằng bootstrap hiệu chỉnh lạc quan (200 lần, khớp lại toàn bộ pipeline mỗi lần); repeated stratified 5-fold CV (20 lần, phân tầng theo sự kiện) là phân tích phụ. Siêu tham số chọn bằng CV lặp ở vòng trong. | Mỗi fold chỉ có 4–5 sự kiện | D1 |
| 6.4 Mô hình hóa | nằm trong vòng ngoài | nằm trong mỗi lần bootstrap và mỗi phần CV | Đồng bộ D1 | D1 |
| 6.4 Mô hình hóa | Mỗi tổ hợp 3 mô hình × 4 bộ đặc trưng (12 tổ hợp) được đánh giá cùng một cách chia dữ liệu để so sánh công bằng. | CoxNet trên 4 bộ đặc trưng là phân tích chính; 12 tổ hợp 3 mô hình × 4 bộ đặc trưng là phân tích thăm dò, đánh giá trên cùng các mẫu bootstrap và cùng cách chia dữ liệu để so sánh công bằng. | Ít sự kiện | D1 |
| 6.5 Đánh giá hiệu năng | (sau «(Harrell và Uno») | ; Uno cắt tại τ = 3 năm | Theo dõi ngắn | D1 |
| 6.5 Đánh giá hiệu năng | tại 3 và 5 năm | tại 2 và 3 năm (TCGA chỉ còn 31 ca được theo dõi đến 5 năm) | Theo dõi ngắn | D1 |
| 6.5 Đánh giá hiệu năng | (đoạn mới) | Đối chứng null: (i) cùng miRNA và gen đưa vào thành đặc trưng riêng; (ii) cặp ngẫu nhiên giữ miRNA, thay gen bằng gen không phải đích cùng thập phân vị biểu hiện (200 lần, p thực nghiệm); (iii) tuổi, giai đoạn, hạch + điểm miễn dịch + purity; (iv) (iii) + đặc trưng cặp. | Đối chứng null | D1, D4 |
| 6.5 Đánh giá hiệu năng | (đoạn mới) | Đánh giá nội bộ: bootstrap hiệu chỉnh lạc quan; ở METABRIC, theo dõi được cắt tại 5 năm (kiểm duyệt hành chính) để khớp với TCGA. | Ít sự kiện; theo dõi ngắn | D1, D3 |
| 6.7 Bảng công cụ – scipy, statsmodels | (cuối đoạn) | , tương quan từng phần | Tương quan từng phần | D4 |
| 6.7 Bảng công cụ – scipy, statsmodels | 2 | 1, 2 | Dùng ở Bước 11 của ND1 |  |
| 6.7 Bảng công cụ | (hàng mới) | PanCanAtlas (ABSOLUTE; Thorsson et al., 2018) \| Bảng dữ liệu \| Độ tinh khiết, thâm nhiễm miễn dịch \| 1, 2, 5 | Hiệp biến thành phần mô | D4 |
| 6.7 Bảng công cụ | (hàng mới) | GENCODE v36 GTF (quét bằng pandas) \| Python \| Cờ gen chồng lấn trong STAR – Counts \| 1 | Lỗi định lượng gen chồng lấn |  |
| 6.7 Bảng công cụ | (hàng mới) | UCSC Xena TOIL RSEM \| Bảng dữ liệu \| Độ nhạy định lượng mRNA \| 1 | Đối chiếu STAR với RSEM |  |
| 7. Kết quả dự kiến | trong nested cross-validation | (CoxNet là phân tích chính) đánh giá bằng bootstrap hiệu chỉnh lạc quan | Đồng bộ D1 | D1 |
| 8. Rủi ro – Ít sự kiện PFI | (cuối đoạn) | . TRẠNG THÁI 11/2026: ĐÃ XẢY RA – 22 → 26 sự kiện PFI sau khi sửa định nghĩa HER2 (152 ca; 54 ca còn nguy cơ ở 3 năm, 31 ca ở 5 năm); phương án dự phòng đã kích hoạt 11/2026 | Kết quả mốc kiểm tra (thay cho cột trạng thái) | D1 |
| 8. Rủi ro – Ít sự kiện PFI | (cuối đoạn) | ; bổ sung 11/2026: CoxNet là phân tích chính, bootstrap hiệu chỉnh lạc quan, AUC tại 2 và 3 năm, tối đa 5 cặp, kiểm chứng hai chiều nhánh mRNA, cắt METABRIC tại 5 năm | Phương án bổ sung | D1, D2, D3 |
| 8. Rủi ro – miRNA METABRIC | (cuối đoạn) | . Ứng viên GEO cụ thể: GSE22220/GSE22219 (Buffa et al., 2011; miRNA + mRNA, tái phát xa); SCAN-B GSE96058 (RNA-seq, chỉ kiểm chứng nhánh mRNA, cần kiểm tra endpoint có sẵn); cần kiểm tra số ca TNBC trước khi dùng | Ứng viên cụ thể |  |
| 8. Rủi ro – Định nghĩa TNBC METABRIC | (cuối đoạn) | . Kết quả 11/2026: không có IHC cho PR và HER2; định nghĩa tổ hợp cho 276 ca; 18 ca ER bất đồng được gắn cờ | Đã xác minh | D7 |
| 8. Rủi ro – PFI và RFS | (cuối đoạn) | ; bổ sung 11/2026: thêm DSS ở cả hai cohort, cắt RFS của METABRIC tại 5 năm, báo cáo loại sự kiện PFI (5/26 sự kiện là u nguyên phát mới) | Khớp endpoint | D3 |
| 8. Rủi ro – rủi ro mới | (hàng mới) | Bộ lọc tương quan và signature phản ánh thành phần mô (TIL, độ tinh khiết) thay vì điều hòa \| Cao \| Số cặp dương có ý nghĩa ≥ số cặp âm; miRNA tương quan mạnh với điểm miễn dịch (11/2026: 196.113 cặp dương so với 165.509 cặp âm trên 12.511 gen; 9 miRNA có \|ρ\| ≥ 0.5 với điểm miễn dịch) \| Tương quan từng phần; đối chứng (iii)–(iv); kiểm tra cặp sau hiệu chỉnh ở ND5 | Phát hiện từ phân tích mô tả 11/2026 | D4, D7, D8 |
| 8. Rủi ro – rủi ro mới | (hàng mới) | miRNA-seq mang tính thành phần (miR-21-5p chiếm khoảng 24,8% RPM) \| Trung bình \| Đặc trưng thay đổi khi miR-21 thay đổi \| Ngưỡng theo read; độ nhạy bỏ miR-21, upper-quartile, DESeq2; đặc trưng hiệu phân vị | Phát hiện từ phân tích mô tả 11/2026 | D4, D7, D8 |
| 8. Rủi ro – rủi ro mới | (hàng mới) | Lỗi định lượng gen chồng lấn trong STAR – Counts \| Trung bình \| Gen biểu hiện cao nhưng TPM thấp bất thường (CD68, EIF4A1); đối chiếu RSEM (11/2026): 726 gen bị STAR đếm thiếu từ 4 lần trở lên \| Loại 314 gen bị STAR định lượng sai (đối chiếu RSEM của TOIL) khỏi phân tích chính: ứng viên tạo cặp, tập gen tính phân vị f_rank và điểm thành phần ở cả hai cohort, baseline gen thô, ứng viên gen tham chiếu f_ref; giữ trong đầu vào ND3. Độ nhạy: chọn cặp lại bằng RSEM TPM (D9) | Phát hiện từ phân tích mô tả 11/2026 | D4, D7, D8, D9 |
| 8. Rủi ro – rủi ro mới | (hàng mới) | Hai cohort không đồng nhất (LAR, claudin-low, ER bất đồng, thời kỳ điều trị) \| Trung bình \| Phân bố PAM50, cờ ER, điều trị khác nhau (METABRIC: 32% claudin-low; 161/276 ca hóa trị) \| Bảng 1; hiệu chỉnh hóa trị; phân tích độ nhạy loại ca bất đồng | Phát hiện từ phân tích mô tả 11/2026 | D4, D7, D8 |
| 9. Kế hoạch | (đoạn mới) | Kết quả mốc kiểm tra 11/2026 (chốt 04/11/2026): 152 bệnh nhân TNBC, 26 sự kiện PFI, dưới ngưỡng 40, nên áp phương án dự phòng ở mục 8. Tám quyết định D1–D8 (đánh giá, kiểm chứng hai chiều, endpoint, bộ lọc tương quan, định nghĩa đặc trưng, HER2, ER-low/LAR, ngưỡng lọc) được ghi trong nhật ký dự án và tóm tắt ở docs/thay_doi_de_cuong_v2.md. Dữ liệu miRNA METABRIC chưa có. | Ghi kết quả mốc | D1–D8 |
| 10. Hướng phát triển | log-ratio | đặc trưng f^ref | qPCR không đo RPM/TPM | D5 |
| 10. Hướng phát triển | Ưu thế của đặc trưng log-ratio. | Ưu thế của đặc trưng chuẩn hóa theo tham chiếu. | qPCR không đo RPM/TPM | D5 |
| 10. Hướng phát triển |  Đặc trưng tính trên từng mẫu nên có thể đo bằng qPCR trên từng bệnh nhân mà không cần một cohort tham chiếu. |  qPCR không đo ra RPM hay TPM; đặc trưng chuyển giao là f^ref (chuẩn hóa theo miRNA và gen tham chiếu, tương đương ΔCt), nên vẫn đo được trên từng bệnh nhân mà không cần một cohort tham chiếu. Các phân tử tham chiếu được chọn ở ND2 theo độ ổn định, không dùng kết cục. | qPCR không đo RPM/TPM | D5 |
| Tài liệu tham khảo | (đoạn mới) | Wolff AC, et al. Human epidermal growth factor receptor 2 testing in breast cancer: American Society of Clinical Oncology/College of American Pathologists clinical practice guideline focused update. Journal of Clinical Oncology. 2018;36(20):2105–2122. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Allison KH, et al. Estrogen and progesterone receptor testing in breast cancer: ASCO/CAP guideline update. Journal of Clinical Oncology. 2020;38(12):1346–1366. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Lehmann BD, et al. Refinement of triple-negative breast cancer molecular subtypes: implications for neoadjuvant chemotherapy selection. PLoS One. 2016;11(6):e0157368. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Thorsson V, et al. The immune landscape of cancer. Immunity. 2018;48(4):812–830. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Carter SL, et al. Absolute quantification of somatic DNA alterations in human cancer. Nature Biotechnology. 2012;30(5):413–421. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. Nature Communications. 2015;6:8971. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Vivian J, et al. Toil enables reproducible, open source, big biomedical data analyses. Nature Biotechnology. 2017;35(4):314–316. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Buffa FM, et al. microRNA-associated progression pathways and potential therapeutic targets identified by integrated mRNA and microRNA expression profiling in breast cancer. Cancer Research. 2011;71(17):5635–5645. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Telonis AG, et al. Beyond the one-locus-one-miRNA paradigm: microRNA isoforms enable deeper insights into breast cancer heterogeneity. Nucleic Acids Research. 2015;43(19):9158–9175. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Geman D, et al. Classifying gene expression profiles from pairwise mRNA comparisons. Statistical Applications in Genetics and Molecular Biology. 2004;3(1):Article 19. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Tan AC, et al. Simple decision rules for classifying human cancers from gene expression profiles. Bioinformatics. 2005;21(20):3896–3904. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Uno H, et al. On the C-statistics for evaluating overall adequacy of risk prediction procedures with censored survival data. Statistics in Medicine. 2011;30(10):1105–1117. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Riley RD, et al. Calculating the sample size required for developing a clinical prediction model. BMJ. 2020;368:m441. | Trích dẫn mới; metadata đối chiếu qua DOI |  |
| Tài liệu tham khảo | (đoạn mới) | Gloor GB, et al. Microbiome datasets are compositional: and this is not optional. Frontiers in Microbiology. 2017;8:2224. | Trích dẫn mới; metadata đối chiếu qua DOI |  |

## Đợt 08/11/2026 – D5 (ngoại lệ), D10, D11, thứ bậc phân tích

- **Mục 1:** đặc trưng chính f_rank (D5) đã có từ đợt 04/11/2026, không sửa thêm. Ngoại lệ của D5 (không có miRNA METABRIC thì f_ref làm chính) được ghi ở mục 6.2.
- **Mục 2.2, dòng Xu 2025:** mô tả cũ ("chọn 8 miRNA", "đặc trưng vẫn là mức biểu hiện miRNA") không khớp với tóm tắt bài (Crossref): bài dùng 4 miRNA, 8 gen đích và random forest trên 10 cặp miRNA–mRNA. Đã sửa theo yêu cầu.
- **Metadata** của Lee 2023, Hussen 2025 và Xu 2025 đối chiếu qua Crossref ngày 08/11/2026. Trang MDPI chặn truy cập tự động nên **chưa đọc được toàn văn**; các chi tiết phương pháp ghi trong bảng 2.2 (CPM, |ΔPCC|, 10 loại ung thư, lọc theo nhãn trước khi chia dữ liệu, cohort kiểm chứng 23 ca, AUC trên dữ liệu huấn luyện) lấy theo yêu cầu của Vương, cần đối chiếu toàn văn trước khi nộp.
- **Tài liệu số 12** (Xu 2025, IJMS 26:1916) đã có trong danh mục: chỉ bổ sung tác giả, tiêu đề đầy đủ và DOI, không thêm trùng.
- **Hussen 2025** là bài tổng quan (AI và dấu ấn RNA trong ung thư), được trích ở mục 2.4 cho ý kiểm chứng ngoài là rào cản chính.

| Mục | Nội dung cũ | Nội dung mới | Lý do | Mã |
|---|---|---|---|---|
| 2.2 Tổng quan – bảng – dòng Xu 2025 | (sau «(») | Xu et al.,  | Ghi tác giả |  |
| 2.2 Tổng quan – bảng – dòng Xu 2025 | Chọn 8 miRNA qua tương quan điều hòa và ML | Chọn miRNA theo biểu hiện khác biệt u–lành; random forest trên 10 cặp miRNA–mRNA, đặc trưng là tích miRNA × mRNA × ρ | Mô tả cũ chưa đúng: bài dùng đặc trưng cặp dạng tích |  |
| 2.2 Tổng quan – bảng – dòng Xu 2025 | Đặc trưng vẫn là mức biểu hiện miRNA, không phải đặc trưng cặp | Cohort kiểm chứng 23 ca; AUC tính trên dữ liệu huấn luyện; đặc trưng tích phụ thuộc thang đo nên khó chuyển nền tảng | Mở rộng giới hạn |  |
| 2.2 Tổng quan – bảng | (hàng mới) | Lee et al. (IJMS, 2023) \| ΔPCC – độ thay đổi tương quan miRNA–RNA khi thêm một mẫu u vào tập mô lành; 10 loại ung thư TCGA \| miRNA đo từ RNA-seq dài (không phải dạng trưởng thành), không dùng bằng chứng trình tự, lọc đặc trưng theo nhãn trước khi chia dữ liệu, tiên lượng không kiểm chứng ngoài | Đối chứng thăm dò ΔPCC | D11 |
| 2.4 Khoảng trống | (đoạn mới) | Các đặc trưng cặp đã công bố thuộc dạng tích (Xu et al., 2025) hoặc độ lệch so với mô lành (Lee et al., 2023); chưa có so sánh các họ đặc trưng này trong cùng một thiết kế có kiểm chứng ngoài, trong khi kiểm chứng ngoài vẫn là rào cản chính khi đưa dấu ấn RNA vào lâm sàng (Hussen et al., 2025). | Định vị so với Lee 2023 và Xu 2025 | D11 |
| 5.2 Nội dung 2 – các bước | (đoạn mới) | Tính điểm hoạt động miRNA (D10, bộ đặc trưng phụ): trung bình f^rank trên các đích đã qua bộ lọc của mỗi miRNA, chỉ cho miRNA có ≥ 5 đích; trọng số cố định, không học từ kết cục. | Điểm hoạt động miRNA | D10 |
| 5.2 Nội dung 2 – các bước | (đoạn mới) | Tính ΔPCC (D11, đối chứng thăm dò) trên cùng tập cặp, tham chiếu là 104 mẫu mô lành TCGA-BRCA có cả miRNA và mRNA (đếm ở Bước 12 của Nội dung 1), xử lý cùng cách và cùng tập phân tử với mẫu u. | ΔPCC dùng cùng tập cặp với f_rank | D11 |
| 5.4 Nội dung 4 – các bước – bước 1 | (cuối đoạn) |  Kế hoạch gồm thứ bậc phân tích đã chốt ngày 08/11/2026 (mục 6.4): một phép so sánh chính (CoxNet trên f^rank cặp so với biểu hiện gen thô và embedding), các phân tích phụ và thăm dò, quy tắc dự phòng D10 và thứ tự bỏ bớt khi thiếu thời gian. | Thứ bậc phân tích đăng ký trước | D10, D11 |
| 6.2 Xây dựng đặc trưng cặp | (đoạn mới) | Ngoại lệ (D5): nếu không có dữ liệu miRNA METABRIC, f^ref là đặc trưng chính, vì khi đó không kiểm chứng được tính bất biến xuyên nền tảng của f^rank và hướng chuyển giao còn lại là qPCR (mục 10). | Ngoại lệ của D5 | D5 |
| 6.2 Xây dựng đặc trưng cặp | (đoạn mới) | Điểm hoạt động miRNA (D10, bộ đặc trưng phụ). Với miRNA i có tập đích Tᵢ đã qua bộ lọc, chỉ tính khi \|Tᵢ\| ≥ 5: | Định nghĩa | D10 |
| 6.2 Xây dựng đặc trưng cặp | (đoạn mới) | A_i(s) = \frac{1}{\|T_i\|} \sum_{j \in T_i} f^{\mathrm{rank}}_{ij}(s) | Công thức | D10 |
| 6.2 Xây dựng đặc trưng cặp | (đoạn mới) | Trọng số cố định (trung bình), không học từ kết cục; A bất biến với thang đo như f^rank và ít phụ thuộc vào việc chọn đúng từng cặp. Quy tắc dự phòng: nếu trong bootstrap không có cặp nào được chọn ở ≥ 50% số lần, điểm hoạt động miRNA thay f^rank cặp làm đặc trưng chính. | Quy tắc dự phòng | D10 |
| 6.2 Xây dựng đặc trưng cặp | (đoạn mới) | ΔPCC (D11, đối chứng thăm dò; Lee et al., 2023). Độ thay đổi hệ số tương quan Pearson của cặp khi thêm mẫu u s vào tập N gồm n mẫu mô lành: | Định nghĩa | D11 |
| 6.2 Xây dựng đặc trưng cặp | (đoạn mới) | \Delta\mathrm{PCC}_{ij}(s) = r_{ij}\big(\mathcal{N} \cup \{s\}\big) - r_{ij}\big(\mathcal{N}\big) | Công thức gốc Lee 2023 | D11 |
| 6.2 Xây dựng đặc trưng cặp | (đoạn mới) | Lee et al. (2023) dùng \|ΔPCC\| trên CPM. Đề tài tính trên thang log2 (cùng đầu vào với f^raw), vì Pearson trên CPM thô bị vài giá trị cực đoan chi phối; ΔPCC có dấu là biến thể chính, \|ΔPCC\| là biến thể phụ; tập cặp là tập đã qua bộ lọc ở trên, để so sánh công bằng với f^rank. Theo hàm ảnh hưởng của hệ số tương quan, với z chuẩn hóa theo trung bình và độ lệch chuẩn của mô lành và r = rᵢⱼ(N): | Quy ước D11 | D11 |
| 6.2 Xây dựng đặc trưng cặp | (đoạn mới) | \Delta\mathrm{PCC}_{ij}(s) \approx \frac{z_i(s)\, z_j(s) - r\,\big(z_i(s)^2 + z_j(s)^2\big)/2}{n + 1} | Dạng xấp xỉ bằng hàm ảnh hưởng | D11 |
| 6.2 Xây dựng đặc trưng cặp | (đoạn mới) | Như vậy ΔPCC gần với tích hai giá trị z so với mô lành, tức cùng họ với đặc trưng tích miRNA × mRNA × ρ của Xu et al. (2025), nhưng phụ thuộc vào trung bình và độ lệch chuẩn của một cohort tham chiếu. Đặc trưng tích là đối chứng thăm dò thứ hai. | Liên hệ ΔPCC với đặc trưng tích | D11 |
| 6.4 Mô hình hóa | (đoạn mới) | Thứ bậc phân tích (đăng ký trước, 08/11/2026). (1) Chính, một phép so sánh: CoxNet trên f^rank cặp so với biểu hiện gen thô và embedding. (2) Phụ: điểm hoạt động miRNA; đối chứng null (i)–(iv) ở mục 6.5. (3) Thăm dò: ΔPCC, đặc trưng tích (Xu et al., 2025), TabPFN, RSF và GBS; f^raw và f^ref là phân tích độ nhạy. | Thứ bậc phân tích | D5, D10, D11 |
| 6.4 Mô hình hóa | (đoạn mới) | Quy tắc dự phòng và thứ tự bỏ bớt. Nếu trong bootstrap không có cặp nào được chọn ở ≥ 50% số lần, điểm hoạt động miRNA thay f^rank cặp làm đặc trưng chính (D10). Thiếu thời gian thì bỏ theo thứ tự: ΔPCC → TabPFN → đặc trưng tích; không bỏ điểm hoạt động miRNA. | Quy tắc dự phòng D10 | D10 |
| 6.5 Đánh giá hiệu năng | (đoạn mới) | Đối chứng thăm dò: ΔPCC (Lee et al., 2023; tham chiếu mô lành TCGA-BRCA, bỏ nếu < 50 mẫu) và đặc trưng tích (Xu et al., 2025), tính trên cùng tập cặp và đánh giá cùng cách với f^rank; không dùng để chọn mô hình cuối. METABRIC không có mẫu mô lành trong gói công khai nên ΔPCC chỉ đánh giá trong TCGA. | Thêm đối chứng thăm dò | D11 |
| 8. Rủi ro – rủi ro mới | (hàng mới) | Signature cặp không ổn định do ít sự kiện \| Cao \| Không có cặp nào được chọn ở ≥ 50% số lần bootstrap \| Quy tắc dự phòng D10: điểm hoạt động miRNA thay f^rank cặp làm đặc trưng chính | Quyết định 08/11/2026 | D10 |
| 8. Rủi ro – rủi ro mới | (hàng mới) | Không đủ mẫu mô lành làm tham chiếu cho ΔPCC \| Thấp \| < 50 mẫu lành TCGA-BRCA có cả miRNA và mRNA (11/2026: có 104 mẫu) \| Bỏ ΔPCC (D11). METABRIC không có mẫu lành trong gói cBioPortal nên ΔPCC không kiểm chứng ngoài được | Quyết định 08/11/2026 | D11 |
| Tài liệu tham khảo | (đầu đoạn) | Xu J, et al.  | Bổ sung tác giả (Crossref) |  |
| Tài liệu tham khảo | (sau «triple-negative breast cancer») | : insights into miR-210, miR-183, miR-21, and miR-181b | Tiêu đề đầy đủ (Crossref) |  |
| Tài liệu tham khảo | (cuối đoạn) |  doi:10.3390/ijms26051916 | Bổ sung DOI; đã có trong danh mục, không thêm trùng |  |
| Tài liệu tham khảo | (đoạn mới) | Lee S, Cho M, Park B, Han K. Finding miRNA–RNA network biomarkers for predicting metastasis and prognosis in cancer. International Journal of Molecular Sciences. 2023;24(5):5052. doi:10.3390/ijms24055052 | Trích dẫn mới; metadata đối chiếu Crossref | D11 |
| Tài liệu tham khảo | (đoạn mới) | Hussen BM, Abdullah SR, Hidayat HJ, Samsami M, Taheri M. Integrating AI and RNA biomarkers in cancer: advances in diagnostics and targeted therapies. Cell Communication and Signaling. 2025;23:430. doi:10.1186/s12964-025-02434-2 | Trích dẫn mới; metadata đối chiếu Crossref | D11 |
