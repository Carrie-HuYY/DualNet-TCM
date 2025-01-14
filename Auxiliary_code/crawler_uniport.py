import requests
import re
import os
from concurrent.futures import ThreadPoolExecutor

# 访问网址并获取原链接
headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/96.0.4664.45 Safari/537.36'}

url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/pan_proteomes/'
res = requests.get(url, headers=headers)
res.encoding = 'utf-8'
page_content = res.text

pattern = r'<a href="([^"]+\.fasta\.gz)"'
matches = re.findall(pattern, page_content)

base_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/pan_proteomes/"

# 下载进指定文件夹
download_folder = 'D:/github_planting/MGM/data/Uniport_data'

if not os.path.exists(download_folder):
    os.makedirs(download_folder)

def download_file(file):
    full_url = base_url + file
    response = requests.get(full_url, headers=headers)
    # 检查请求是否成功
    if response.status_code == 200:
        # 构建完整的文件路径
        file_path = os.path.join(download_folder, file.split('/')[-1])
        with open(file_path, 'wb') as f:
            f.write(response.content)
        print(f'Downloaded {file} to {download_folder}')
    else:
        print(f'Failed to download {file}')

# 使用ThreadPoolExecutor来并发下载文件
with ThreadPoolExecutor(max_workers=10) as executor:
    executor.map(download_file, matches)