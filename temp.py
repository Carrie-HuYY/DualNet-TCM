import pandas as pd
from elasticsearch import Elasticsearch
import os
import time


# 配置Elasticsearch连接
es = Elasticsearch(
    ['http://localhost:9200/']
)
# 定义索引名称
index_name = "abstract22"

folder_path = "D:/Snoopy/result"
base_filename = "pubmed25n"
start_num = 1016
end_num = 1275

print("start parsing")
for num in range(start_num, end_num + 1):
    start_time = time.time()
    num_str = f"{num:04d}"
    file_path = os.path.join(folder_path, f"{base_filename}{num_str}.xlsx")
    df = pd.read_excel(file_path)
    data_to_insert = df[["PMID", "Abstract"]].rename(columns={"PMID": "pubMedId", "Abstract": "abstract"}).to_dict(orient="records")

    for idx, item in enumerate(data_to_insert, start=1):
        try:
            response = es.index(index=index_name, body=item)
            print(f"第 {idx} 条数据插入成功")
        except Exception as e:
            print(f"第 {idx} 条数据插入失败：", e)

    # 验证数据
    query = {
        "query": {
            "match_all": {}
        }
    }

    response = es.search(index=index_name, body=query)
    print("查询结果：", response)

    print(f"第 {num} 个xml文件插入成功")