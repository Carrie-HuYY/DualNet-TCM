import requests
import re
import os
import time
from concurrent.futures import ThreadPoolExecutor

# 访问网址并获取原链接
headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/96.0.4664.45 Safari/537.36'
}

url = 'https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/'
res = requests.get(url, headers=headers)
res.encoding = 'utf-8'
page_content = res.text

pattern = r'<a href="([^"]+\.xml\.gz)"'
matches = re.findall(pattern, page_content)

matches = matches[1160:len(matches)]

print(matches)

base_url = 'https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/'

# 下载进指定文件夹
download_folder = 'D:/NBSD/Auxiliary_code/pubmed'
extract_folder = 'D:/NBSD/Auxiliary_code/pubmed_extracted'

if not os.path.exists(download_folder):
    os.makedirs(download_folder)

if not os.path.exists(extract_folder):
    os.makedirs(extract_folder)


def download_file(file):
    full_url = base_url + file
    response = requests.get(full_url, headers=headers, stream=True)
    # 检查请求是否成功
    if response.status_code == 200:
        # 获取文件大小
        total_size = int(response.headers.get('content-length', 0))
        downloaded_size = 0
        start_time = time.time()

        # 构建完整的文件路径
        file_path = os.path.join(download_folder, file.split('/')[-1])
        with open(file_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded_size += len(chunk)
                    speed = downloaded_size / (time.time() - start_time)
                    print(
                        f'\rDownloading {file}: {downloaded_size / (1024 * 1024):.2f} MB / {total_size / (1024 * 1024):.2f} MB, Speed: {speed / (1024 * 1024):.2f} MB/s',
                        end='')
        print(f'\nDownloaded {file} to {download_folder}')

        if file.endswith('.gz'):
            extract_path = os.path.join(extract_folder, file.split('/')[-1].replace('.gz', ''))
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(extract_path)
            print(f'Extracted {file} to {extract_path}')
        else:
            print(f'Skipped extraction for {file} as it is not a .gz file')
    else:
        print(f'Failed to download {file}')


# 使用ThreadPoolExecutor来并发下载文件
with ThreadPoolExecutor(max_workers=10) as executor:
    executor.map(download_file, matches)
