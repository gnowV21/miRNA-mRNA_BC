# miRNA-mRNA_BC

Khóa luận tốt nghiệp: **Xây dựng mô hình tiên lượng tái phát ung thư vú bộ ba âm tính dựa trên đặc trưng cặp tương tác miRNA–mRNA và so sánh với biểu diễn từ mô hình nền tảng phiên mã**

 Trần Nguyễn Tiến Vương  * Khoa Sinh học – Công nghệ Sinh học, HCMUS * 

## Cấu trúc thư mục

| Thư mục / file | Nội dung |
| --- | --- |
| `notebooks/` | Notebook Python, mỗi Nội dung một file. Hiện có `ND1_tien_xu_ly_du_lieu.ipynb`. |
| `docs/` | Bản xuất của đề cương và tab hướng dẫn (Word hoặc PDF), tải từ Claude Doc. |
| `data/` | Dữ liệu đã xử lý (Parquet, CSV) do notebook lưu vào. |
| `qc/` | Hình kiểm tra, sơ đồ chọn mẫu, bảng báo cáo số sự kiện. |
| `nhat_ky.md` | Nhật ký quyết định và con số của dự án. Ghi mỗi khi chốt một quyết định. |

## Thư mục làm việc

Thư mục này là thư mục làm việc chính: `C:\Users\TIEN VUONG\Desktop\TienVuong\miRNA-mRNA_BC`.
Notebook tự nhận biết đang chạy ở đâu:

| Chạy ở đâu | Thư mục làm việc | Dùng cho |
| --- | --- | --- |
| Máy cá nhân (Jupyter, VS Code, Anaconda) | Thư mục này | Nội dung 1 (không cần GPU) |
| Google Colab | `MyDrive/miRNA-mRNA_BC` trên Google Drive | Các bước cần GPU (Nội dung 2–3) |

File thô tải từ GDC khi chạy trên máy nằm ở `_tam\` (có thể xóa sau khi xong Nội dung 1, đã loại khỏi Git).
Khi cần chạy trên Colab, chép thư mục `data\` lên `MyDrive/miRNA-mRNA_BC/data` trên Drive.

## Cách chạy Nội dung 1 trên máy cá nhân

1. Cài Python (khuyến nghị Anaconda) hoặc VS Code kèm tiện ích Jupyter.
2. Mở `notebooks\ND1_tien_xu_ly_du_lieu.ipynb`.
3. Chạy lần lượt từ Bước 0 đến Bước 10. Bước 0 tự cài các gói còn thiếu.
4. Làm theo tab **Hướng dẫn Nội dung 1** trong đề cương: phần kiểm tra và ghi nhật ký của từng bước.

## Lưu ý về dữ liệu

- Không đưa dữ liệu bệnh nhân (thư mục `data/`) lên GitHub công khai. File `.gitignore` đã loại sẵn.
- miRNA METABRIC là dữ liệu truy cập có kiểm soát: chỉ dùng bản tải hợp lệ (Synapse hoặc EGA, bộ EGAD00010000438), không dùng bản sao không chính thức.
