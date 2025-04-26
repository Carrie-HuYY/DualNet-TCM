# import urllib.request
# import os
# import time

# url = "https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/"
# path = "F:/Project/OTTM/OTTM/data/ABS2/"

# if not os.path.exists(path):
#     os.makedirs(path)

# start_time = time.time()
# for i in range(1, 1167):
#     file_name = "pubmed23n" + str(i).zfill(4) + ".xml.gz"
#     file_url = url + file_name
#     file_path = path + file_name
#     if not os.path.exists(file_path):
#         urllib.request.urlretrieve(file_url, file_path)
#         print("Downloaded", file_name)
#     else:
#         print(file_name, "already exists")
#     if time.time() - start_time > 300:
#         os.system("python F:/Project/OTTM/OTTM/data/ABS/download.py")

# 重启程序


import os
import urllib.request
import time
import concurrent.futures

def download_file(file_url, file_path):
    if not os.path.exists(file_path):
        urllib.request.urlretrieve(file_url, file_path)
        print("Downloaded", file_path)
    else:
        print(file_path, "already exists")


# 执行程序
# 这里range内的数字是文件名的数字，从1开始，需要手动数次每次下载100到200个
# 不能一次下载完成 中间会出现不下载的情况
if __name__ == "__main__":


    url = "https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/"
    path = "F:/Project/OTTM/OTTM/data/ABS2/"


    start_time = time.time()
# 一共1166个 序列到1167
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        for i in range(1100, 1167):
            file_name = "pubmed23n" + str(i).zfill(4) + ".xml.gz"
            file_url = url + file_name
            file_path = path + file_name
            executor.submit(download_file, file_url, file_path)

