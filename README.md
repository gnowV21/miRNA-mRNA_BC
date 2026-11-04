# miRNA-mRNA_BC

Khóa luận tốt nghiệp: **Xây dựng mô hình tiên lượng tái phát ung thư vú bộ ba âm tính dựa trên đặc trưng cặp tương tác miRNA–mRNA và so sánh với biểu diễn từ mô hình nền tảng phiên mã**

 Trần Nguyễn Tiến Vương  * Khoa Sinh học – Công nghệ Sinh học, HCMUS * 

## Cấu trúc thư mục

| Thư mục / file | Nội dung |
| --- | --- |
| `notebooks/` | Notebook Python, mỗi Nội dung một file. Hiện có `ND1_tien_xu_ly_du_lieu.ipynb`. Bản trước khi sửa 11/2026 ở `notebooks/_luu/`. |
| `src/mirna_bc/` | Hàm dùng chung cho các notebook: điểm thành phần mô (`thanh_phan.py`), đặc trưng cặp và tương quan từng phần (`dac_trung.py`), cờ gen chồng lấn (`chong_lan.py`), khung đánh giá mô hình cho ND4 (`danh_gia.py`). |
| `tests/` | Test cho `src/mirna_bc` trên dữ liệu giả. Chạy: `python -m pytest tests`. |
| `docs/` | Bản xuất của đề cương và tab hướng dẫn (Word hoặc PDF), tải từ Claude Doc. Đề cương v2 (tracked changes) và bản đọc soát; `docs/thay_doi_de_cuong_v2.md` tóm tắt các thay đổi. |
| `data/` | Dữ liệu đã xử lý (Parquet, CSV) do notebook lưu vào. |
| `qc/` | Hình kiểm tra, sơ đồ chọn mẫu, bảng báo cáo số sự kiện, phiên bản dữ liệu và manifest md5. |
| `nhat_ky.md` | Nhật ký quyết định và con số của dự án. Ghi mỗi khi chốt một quyết định. |

## Thư mục làm việc

Thư mục này là thư mục làm việc chính: `D:\Users\TIEN VUONG\Desktop\TienVuong\miRNA-mRNA_BC`.
Notebook tự dò thư mục dự án (đi ngược từ thư mục đang mở tới khi gặp `nhat_ky.md`); có thể ghi đè bằng biến môi trường `MIRNA_BC_DIR`.
Notebook tự nhận biết đang chạy ở đâu:

| Chạy ở đâu | Thư mục làm việc | Dùng cho |
| --- | --- | --- |
| Máy cá nhân (Jupyter, VS Code, Anaconda) | Thư mục này | Nội dung 1 (không cần GPU) |
| Google Colab | `MyDrive/miRNA-mRNA_BC` trên Google Drive | Các bước cần GPU (Nội dung 2–3) |

File thô tải từ GDC và các nguồn khác khi chạy trên máy nằm ở `_tam\` (đã loại khỏi Git). **Không xóa** các file isoform miRNA trong `_tam\`: phần isomiR ở Nội dung 2 cần dùng.
Khi cần chạy trên Colab, chép thư mục `data\` lên `MyDrive/miRNA-mRNA_BC/data` trên Drive.

## Cách chạy Nội dung 1 trên máy cá nhân

1. Cài Python (khuyến nghị Anaconda) hoặc VS Code kèm tiện ích Jupyter.
2. Mở `notebooks\ND1_tien_xu_ly_du_lieu.ipynb`.
3. Chạy lần lượt từ trên xuống (Run All được). Bước 0 tự cài các gói còn thiếu. Thứ tự các bước (đã đổi 11/2026, các ô trong notebook đã xếp đúng thứ tự này):
   0 → 1 → 2 → 2b → 3 → 4 → 5 → 5c → 9 → 9b → 6 → 7 → 7b → 8 → 10 → 11.
   Bước 9 (METABRIC) chạy trước Bước 6–7 vì Bước 7 cần danh sách gen đo được trên microarray.
4. Làm theo tab **Hướng dẫn Nội dung 1** trong đề cương: phần kiểm tra và ghi nhật ký của từng bước.

## Lưu ý về dữ liệu

- Không đưa dữ liệu bệnh nhân (thư mục `data/`) lên GitHub công khai. File `.gitignore` đã loại sẵn.
- miRNA METABRIC là dữ liệu truy cập có kiểm soát: chỉ dùng bản tải hợp lệ (Synapse hoặc EGA, bộ EGAD00010000438), không dùng bản sao không chính thức.
