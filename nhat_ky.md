# Nhật ký dự án miRNA-mRNA_BC

Ghi mỗi quyết định và con số kèm ngày. Reviewer và hội đồng sẽ hỏi lý do của từng lựa chọn.

## Quyết định thiết kế đã chốt (10/2026)

- Đối tượng: TNBC theo IHC (ER−, PR−, HER2−); HER2 Equivocal chỉ giữ khi FISH âm.
- Loại: bệnh nhân nam; đã điều trị tân bổ trợ; di căn xa lúc chẩn đoán (M1).
- Huấn luyện: TCGA-BRCA. Kiểm chứng ngoài: METABRIC.
- Endpoint chính: PFI (TCGA-CDR). Phụ: OS ở cả hai cohort; DFI khi cohort kiểm chứng chỉ có RFS/DFS.
- METABRIC: định nghĩa IHC chặt; nếu thiếu IHC cho PR/HER2, chạy song song định nghĩa tổ hợp.
- Mỗi bệnh nhân 1 mẫu: ưu tiên lọ A, rồi số đọc miRNA cao hơn.
- miRNA mức trưởng thành (MIMAT), RPM. mRNA: gen mã hóa protein, TPM. VST chỉ cho phân tích độ nhạy.
- Nguồn cặp miRNA–mRNA: miRBind2, lọc tương quan âm trong từng vòng CV.
- So sánh: đặc trưng cặp vs BulkRNABert (bản không TCGA; thay thế scFoundation) vs BulkFormer vs biểu hiện gen thô.
- Mô hình: CoxNet (chính), RSF, Gradient Boosting Survival; nested CV; stability selection; quy tắc 1-SE.
- Ngôn ngữ: Python; R chỉ cho immunedeconv và oncoPredict (Nội dung 5).
- Tài nguyên: Colab, Kaggle Notebooks; thuê GPU khi cần.

## Nội dung 1 – Thu thập và tiền xử lý

| Ngày | Bước | Nội dung ghi nhận | Con số |
| --- | --- | --- | --- |
|  | Bước 0 | Phiên bản Python, pandas, pydeseq2, lifelines |  |
|  | Bước 1 | TNBC sau từng tiêu chí loại; số ca thiếu IHC |  |
|  | Bước 2 | Số ca có PFI hợp lệ |  |
|  | Bước 3 | Số miRNA trưởng thành; kiểm tra miR-21-5p |  |
|  | Bước 4 | Số gen mã hóa protein |  |
|  | Bước 5 | Số bệnh nhân cuối có cả miRNA và mRNA |  |
|  | Bước 6–7 | 4 ngưỡng lọc + lý do; số miRNA, gen giữ lại |  |
|  | Bước 8 | Sự kiện PFI / OS / DFI; theo dõi trung vị |  |
|  | Bước 9 | METABRIC: TNBC theo 2 định nghĩa; cách đo PR, HER2; sự kiện RFS, OS |  |
|  | Bước 10 | Kết quả kiểm tra Synapse / EGA |  |

## Mốc kiểm tra cuối 11/2026

- Số sự kiện PFI:
- Quyết định: giữ thiết kế / áp phương án dự phòng
- Lý do:
