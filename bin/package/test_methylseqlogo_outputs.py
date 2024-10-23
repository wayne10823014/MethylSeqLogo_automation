import unittest
import os
import glob
import numpy as np
import pandas as pd

class TestMethylSeqLogoOutputs(unittest.TestCase):
    def setUp(self):
        # 設置生成文件和驗證文件的目錄路徑
        self.generated_dir = '/home/wayne/MethylSeqLogo_automation/Output1'
        self.validation_dir = '/home/wayne/MethylSeqLogo_automation/Output validation data'

        # 獲取生成目錄中的所有文件
        self.generated_files = glob.glob(os.path.join(self.generated_dir, '*'))

    def test_compare_all_files(self):
        # 用來儲存不一致的檔案名稱
        inconsistent_files = []

        # 對生成目錄中的所有文件進行比較
        for generated_file in self.generated_files:
            # 獲取文件名
            file_name = os.path.basename(generated_file)
            validation_file = os.path.join(self.validation_dir, file_name)

            # 檢查驗證文件是否存在
            if not os.path.exists(validation_file):
                inconsistent_files.append(f"驗證文件不存在：{file_name}")
                continue

            # 判斷文件類型，根據文件類型選擇比較方法
            try:
                if file_name.endswith(('.csv', '.txt', '.tsv')):
                    # 對於文本文件，讀取並比較內容
                    self.compare_text_files(generated_file, validation_file)
                elif file_name.endswith(('.png', '.jpg')):
                    # 對於圖像文件，使用圖像處理庫進行比較
                    self.compare_image_files(generated_file, validation_file)
                else:
                    # 其他類型的文件，按二進制文件比較
                    self.compare_binary_files(generated_file, validation_file)
            except AssertionError:
                # 如果比較失敗，加入不一致的檔案名稱
                inconsistent_files.append(file_name)

        # 如果有不一致的文件，顯示它們
        if inconsistent_files:
            self.fail(f"以下文件不一致：{', '.join(inconsistent_files)}")

    def compare_text_files(self, file1, file2):
        # 讀取文件內容
        try:
            # 嘗試讀取為DataFrame（適用於CSV、TSV等表格文件）
            df1 = pd.read_csv(file1, sep=None, engine='python')
            df2 = pd.read_csv(file2, sep=None, engine='python')
            # 比較DataFrame
            pd.testing.assert_frame_equal(df1, df2, check_dtype=False, check_exact=False, rtol=1e-5, atol=1e-8)
        except Exception:
            # 如果不是表格文件，按文本比較
            with open(file1, 'r') as f1, open(file2, 'r') as f2:
                content1 = f1.readlines()
                content2 = f2.readlines()
                if content1 != content2:
                    raise AssertionError

    def compare_binary_files(self, file1, file2):
        # 比較二進制文件
        with open(file1, 'rb') as f1, open(file2, 'rb') as f2:
            content1 = f1.read()
            content2 = f2.read()
            if content1 != content2:
                raise AssertionError

    def compare_image_files(self, file1, file2):
        # 使用PIL庫比較圖像文件
        from PIL import Image

        img1 = Image.open(file1).convert('RGB')
        img2 = Image.open(file2).convert('RGB')
        arr1 = np.array(img1)
        arr2 = np.array(img2)
        if not np.array_equal(arr1, arr2):
            raise AssertionError

if __name__ == '__main__':
    unittest.main()
