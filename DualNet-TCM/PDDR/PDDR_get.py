import os
import json
from elasticsearch import Elasticsearch


def get_targetNum_dict(symbol_list, interaction_num, PPI_DICT):
    """
    将组学数据中与每个差异表达蛋白相互作用的蛋白进行汇总，统计每个蛋白出现的次数
    根据interaction_num的值，筛选出与差异表达蛋白相互作用次数大于等于interaction_num的蛋白
    """
    ALL_PPI_PROTEIN = []

    for i in symbol_list:
        if i in PPI_DICT.keys():
            ALL_PPI_PROTEIN.extend(PPI_DICT[i])

    PPI_NUMBER = {}
    TARGET_PPI = {}

    for i in ALL_PPI_PROTEIN:
        if i not in PPI_NUMBER.keys():
            PPI_NUMBER[i] = 1
        else:
            PPI_NUMBER[i] += 1
            if PPI_NUMBER[i] >= interaction_num and i not in TARGET_PPI.keys():
                TARGET_PPI[i] = PPI_NUMBER[i]
            elif PPI_NUMBER[i] >= interaction_num and i in TARGET_PPI.keys():
                TARGET_PPI[i] = PPI_NUMBER[i]

    return TARGET_PPI



def get_data(protein_list_path, interaction_num):

    file_name = protein_list_path
    Symbol_list = get_Symbol(file_name)

    print('The number of proteins in the file is: ', len(Symbol_list))

    # step_2.2 获取symbol对应的PPI蛋白 在这里只要跟差异表达蛋白拥有相互作用就进行保留
    # interaction_num = 0 已经存在相互作用的蛋白 最少应该和几个差异表达蛋白列表中的蛋白相互作用,可以设置为参数

    Symbol_PPI_list = get_PPI_Symbol_List(Symbol_list, interaction_num)

    print('The number of PPI proteins list is: ', len(Symbol_PPI_list))
    print('Please wait for a while, the program is running...')

    with open ('data/Drug/Symbol_To_Target.json', 'r') as f:
        Symbol_To_Target_wm = json.load(f)

    with open ('data/ID_Transformed/Symbol_To_Fullname.json', 'r') as f:
        Symbol_To_Fullname = json.load(f)

    return Symbol_PPI_list, Symbol_To_Target_wm, Symbol_To_Fullname, Symbol_list


def get_txt():
    '''
    获取当前地址文件夹中所有txt的文件名
    :return: 返回文件名
    '''
    for file_name in os.listdir():
        if file_name.endswith('.txt'):
            return file_name

#
def get_Symbol(file_name):
    """
    返回
    :param file_name:
    :return:
    """
    with open(file_name, 'r', encoding='utf-8') as f:
        f = f.readlines()
    Symbol = []

    for i in f:
        i = i.replace('\t', ' ').replace('\n', '').split(' ')
        Symbol.append(i[0])

    return Symbol

# 'hepatocellular carcinoma'
# 查询药物关于疾病的报道信息
def get_drug_report_info(drug_ap, drug_cl, disease, input_num, es):
    drug_ap_not_report, drug_ap_report, drug_cl_not_report, drug_cl_report = [], [], [], []
    for drug_name in drug_ap:
        # 这个名字需要进行处理
        # 会出现特殊字符无法处理的情况[Avastin+/-Tarceva]
        drug_name = drug_name.replace(
            '+/-', ' ').replace(
            '/', ' ').replace(
            '[', '').replace(
            ']', '').replace(
            '-', ' ')
        query = {
            'query': {
                'bool': {
                    'must': [
                        {
                            "match": {
                                "abstract": drug_name
                            }
                        },
                        {
                            "match_phrase": {
                                "abstract": disease
                            }
                        },
                    ]
                }
            }
        }
        res = es.search(index='abstract22', body=query, scroll='5m')
        reported_number = res['hits']['total']['value']
        es.clear_scroll(scroll_id=res['_scroll_id'])
        if reported_number > input_num:
            drug_ap_report.append(drug_name)
        else:
            drug_ap_not_report.append(drug_name)
    for drug_name in drug_cl:
        drug_name = drug_name.replace(
            '+/-', ' ').replace(
            '/', ' ').replace(
            '[', '').replace(
            ']', '').replace(
            '-', ' ')
        query = {
            'query': {
                'bool': {
                    'must': [
                        {
                            "match": {
                                "abstract": drug_name
                            }
                        },
                        {
                            "match_phrase": {
                                "abstract": disease
                            }
                        },
                    ]
                }
            }
        }
        res = es.search(index='abstract22', body=query, scroll='5m')
        reported_number = res['hits']['total']['value']
        es.clear_scroll(scroll_id=res['_scroll_id'])
        if reported_number > input_num:
            drug_cl_report.append(drug_name)
        else:
            drug_cl_not_report.append(drug_name)
    return drug_ap_not_report, drug_ap_report, drug_cl_not_report, drug_cl_report


# 从药物列表获取药物的频率
def get_drug_frequency(drug_not_report, drug_report, es):
    drug_frequency = []
    if drug_not_report != []:
        for drug in drug_not_report:
            # 这个名字需要进行处理
            # 会出现特殊字符无法处理的情况[Avastin+/-Tarceva]
            drug_name = drug.replace(
                '+/-', ' ').replace(
                '/', ' ').replace(
                '[', '').replace(
                ']', '').replace(
                '-', ' ')
            query = {
                'query': {
                    "match": {
                        "abstract": drug_name
                    }
                }
            }
            res = es.search(index='abstract22', body=query, scroll='5m')
            reported_number = res['hits']['total']['value']
            drug_frequency.append(reported_number)
            es.clear_scroll(scroll_id=res['_scroll_id'])
    else:
        for drug in drug_report:
            drug_name = drug.replace(
                '+/-', ' ').replace(
                '/', ' ').replace(
                '[', '').replace(
                ']', '').replace(
                '-', ' ')
            query = {
                'query': {
                    "match": {
                        "abstract": drug_name
                    }
                }
            }
            res = es.search(index='abstract22', body=query, scroll='5m')
            reported_number = res['hits']['total']['value']
            drug_frequency.append(reported_number)
            es.clear_scroll(scroll_id=res['_scroll_id'])
    return drug_frequency


def get_PPI_Symbol_List(symbol_list, interaction_num):
    """
    获取PPI列表
    """
    with open ('data/PPI/PPI.json', 'r') as f:
        PPI_DICT = json.load(f)

    TARGET_PPI = get_targetNum_dict(symbol_list, interaction_num, PPI_DICT)
    TARGET_PPI_LIST = [i for i in TARGET_PPI.keys() if i not in symbol_list]

    return TARGET_PPI_LIST


if __name__ == '__main__':
    es = Elasticsearch(['https://localhost:9200'])