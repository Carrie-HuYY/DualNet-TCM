import os
import json
import time
from elasticsearch import Elasticsearch
from elasticsearch import helpers

# 获取所有摘要的路径
def get_abstract_path():
    abstract_path = []
    for root, dirs, files in os.walk(os.getcwd()):
        for file in files:
            if 'json' in file:
                abstract_path.append(os.path.join(root,file))
    return abstract_path

# 读取摘要文件并保存为字典
def get_abstract_dict(abstract_path):
    abstract_dict = {}
    with open (abstract_path, 'r', encoding='utf-8') as f:
        abstract_dict = json.load(f)
    return abstract_dict 

#将摘要信息保存成可以被elasticsearch识别的json格式
def normalize_data(abstract_dict):
    data = []
    for key,value in abstract_dict.items():
        data.append({
            # 在这里更改保存的index
            '_index':'abstract22',
            '_id':key,
            '_source':{
                'pubMedId':key,
                'abstract':value
                }
            }
        )
    return data

def data_to_elasticsearch():
    abstractPath = get_abstract_path()# 获取摘要的路径
    for path in abstractPath:
        start = time.perf_counter()
        abstract_dict = get_abstract_dict(path)
        data = normalize_data(abstract_dict)
        end = time.perf_counter()
        print('获取标准化数据的时间')
        print(end - start)

        start = time.perf_counter()
        # 将数据保存到elasticsearch中
        helpers.bulk(es, data)
        end = time.perf_counter()
        print('获取数据存入elasticsearch的时间')
        print(end - start)

if __name__ == '__main__':
    es = Elasticsearch(timeout=30, max_retries=10, retry_on_timeout=True)
    es.indices.delete(index='abstract22', ignore=[400, 404])
    #es.indices.create(index='abstract', ignore=400) #创建索引
    es.indices.create(index='abstract22', ignore=400)
    data_to_elasticsearch()

    abstractPath = get_abstract_path()# 获取摘要的路径
    for path in abstractPath:
        print(path)
        start = time.perf_counter()
        abstract_dict = get_abstract_dict(path)
        data = normalize_data(abstract_dict)
        end = time.perf_counter()
        print('获取标准化数据的时间')
        print(end - start)

        start = time.perf_counter()
        # 将数据保存到elasticsearch中
        helpers.bulk(es, data)
        end = time.perf_counter()
        print('获取数据存入elasticsearch的时间')
        print(end - start)



import concurrent.futures

es = Elasticsearch(timeout=30, max_retries=10, retry_on_timeout=True)
es.indices.delete(index='abstract22', ignore=[400, 404])
es.indices.create(index='abstract22', ignore=400)

abstractPath = get_abstract_path()

def process_data(path):
    abstract_dict = get_abstract_dict(path)
    data = normalize_data(abstract_dict)
    helpers.bulk(es, data)
    print( path + 'done')

with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = []
    for path in abstractPath:
        print(path)
        futures.append(executor.submit(process_data, path))
    for future in concurrent.futures.as_completed(futures):
        pass






