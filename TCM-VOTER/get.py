from collections import Counter
import pandas as pd
import os
import json
from elasticsearch import Elasticsearch

# TODO: 为各函数增加抛出异常功能，若无法查询到相关信息，则抛出异常。
# TODO: 修改函数解释以及其中的变量名，使之与函数作用统一

def get_formula(by, items) -> pd.DataFrame:
    """
        读取HerbiV_formula数据集，返回items中复方的信息。
        Read the HerbiV_formula dataset and return the formula(s) information in items.

        Args:
            by (str): 数据集中与items相匹配的列的列名。Column name of the column in the dataset that matches items.
            items (collections.abc.Iterable): 要查询的复方。Formula(s) to be queried.

        Returns:
            formula: items中复方的信息。Formula(s) information in items.
    """

    formula_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) + r'/Data/Formula.xlsx')
    formula = formula_all.loc[formula_all[by].isin(items)].copy()
    formula.index = range(formula.shape[0])

    return formula


def get_formula_tcm_links(by, items) -> pd.DataFrame:
    """
        读取HerbiV_formula_tcm_links数据集，返回items中复方/中药的复方-中药连接信息。
        Read the HerbiV_formula_tcm_links dataset
        and return the formula(s)-TCM connection information of formula(s)/TCM in items.

        Args:
            by (str):数据集中与items相匹配的列的列名。Column name of the column in the dataset that matches items.
            items (collections.abc.Iterable): 要查询的复方/中药。Formula(s)/TCM to be queried.

        Returns:
            formula_tcm_links(pandas.DataFrame): items中复方/中药的复方-中药连接信息。
            Formula(s)-TCM connection information of formula(s)/TCM in items.
    """

    # 读取HerbiV_formula_tcm_links数据集
    formula_tcm_links_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) +
                                        r'/Data/Formula_TCM_Links.xlsx')

    # 在数据集中获取items中复方/中药的复方-中药连接信息
    formula_tcm_links = formula_tcm_links_all.loc[formula_tcm_links_all[by].isin(items)].copy()

    # 重新设置索引
    formula_tcm_links.index = range(formula_tcm_links.shape[0])

    return formula_tcm_links


def get_tcm(by, items) -> pd.DataFrame:
    """
        读取HerbiV_tcm数据集，返回items中中药的信息。
        Read the HerbiV_tcm dataset and return the TCM information in items.

        Args:
            by (str): 数据集中与items相匹配的列的列名。Column name of the column in the dataset that matches items.
            items (collections.abc.Iterable): 要查询的中药。TCM to be queried.

        Returns:
            pandas.DataFrame: items中中药的信息。TCM information in items.
    """

    # 读取HerbiV_tcm数据集
    tcm_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) + r'/Data/TCM.xlsx')

    # 在数据集中获取items中中药的信息
    tcm = tcm_all.loc[tcm_all[by].isin(items)].copy()

    # 重新设置索引
    tcm.index = range(tcm.shape[0])

    return tcm


def get_tcm_chem_links(by, items) -> pd.DataFrame:
    """
        读取HerbiV_tcm_chemical_links数据集，返回items中中药/化合物的中药-成分（化合物）连接信息。
        Read the HerbiV_tcm_chemical_links dataset and
        return the TCM-ingredient(s)(chemical(s)) information of TCM/chemical(s) in items.

        Args:
            by (str): 数据集中与items相匹配的列的列名。Column name of the column in the dataset that matches items.
            items (collections.abc.Iterable): 要查询的中药/化合物。TCM/chemical(s) to be queried.

        Returns:
            pandas.DataFrame: items中中药/化合物的中药-成分连接信息。TCM-ingredient(s) information of TCM/chemical(s) in items.
    """

    # 读取HerbiV_tcm_chemical_links数据集
    current_directory = os.path.dirname(os.path.abspath(__file__))
    tcm_chem_links_all = pd.read_excel(current_directory + r'/Data/TCM_Chemical_Links.xlsx')

    # 在数据集中获取items中中药/化合物的中药-成分连接信息
    tcm_chem_links = tcm_chem_links_all.loc[tcm_chem_links_all[by].isin(items)].copy()

    # 重新设置索引
    tcm_chem_links.index = range(tcm_chem_links.shape[0])

    return tcm_chem_links


def get_chemicals(by, items) -> pd.DataFrame:
    """
        读取HerbiV_chemicals数据集，返回items中化合物的信息。
        Read the HerbiV_chemicals dataset and return the chemical(s) information in items.

        Args:
            by (str): 数据集中与items相匹配的列的列名。Column name of the column in the dataset that matches items.
            items (collections.abc.Iterable): 要查询的化合物。Chemical(s) to be queried.

        Returns:
            pandas.DataFrame: items中化合物的信息。Chemical(s) information in items.
    """

    # 读取HerbiV_chemical_protein_links数据集
    current_directory = os.path.dirname(os.path.abspath(__file__))
    chem_all = pd.read_excel(current_directory + r'/Data/Chemical.xlsx')

    # 在数据集中获取items中化合物的信息
    chem = chem_all.loc[chem_all[by].isin(items)].copy()

    # 重新设置索引
    chem.index = range(chem.shape[0])

    return chem


def get_chem_protein_links(by, items, score=900) -> pd.DataFrame:
    """
        读取HerbiV_chemical_protein_links数据集，
        返回items中化合物/蛋白的化合物-靶点（蛋白）连接的combined_score(s)大于等于score的连接信息。
        Read the HerbiV_chemical_protein_links dataset and
        return chemical(s)-target(s)(protein(s)) connection information
        for which the combined_score of the chemical(s)/protein(s) in items is no less than the score.

        Args:
            by (str): 数据集中与items相匹配的列的列名。Column name of the column in the dataset that matches items.
            items (collections.abc.Iterable): 要查询的化合物/蛋白。Chemical(s)/protein(s) to be queried.
            score (int): 仅combined_score大于等于score的记录会被筛选出，默认为900，最大为1000，最小为0。
            Record(s) with combined_score no less than score will be filtered out, 900 by default.

        Returns:
            pandas.DataFrame: items中化合物/蛋白的化合物-靶点（蛋白）连接的combined_score大于等于score的连接信息。
            Chemical(s)-target(s)(protein(s)) connection information for which
            the combined_score of the chemical(s)/protein(s) is no less than the score in items.
    """

    # 读取HerbiV_chemical_protein_links数据集
    current_directory = os.path.dirname(os.path.abspath(__file__))
    chem_protein_links_all = pd.read_excel(current_directory + r'/Data/Chemical_Protein_Links.xlsx')

    # 在数据集中获取items中化合物/蛋白的化合物-靶点（蛋白）连接的combined_score大于等于score的连接信息
    chem_protein_links = chem_protein_links_all.loc[
        (chem_protein_links_all[by].isin(items)) &
        (chem_protein_links_all['Combined_score'] >= score)].copy()

    # 将Combined_score变换为0-1的浮点数
    chem_protein_links['Combined_score'] = chem_protein_links['Combined_score'].astype(float)
    chem_protein_links.loc[:, 'Combined_score'] = chem_protein_links.loc[:, 'Combined_score'].apply(
        lambda x: x / 1000)

    # 重新设置索引
    chem_protein_links.index = range(chem_protein_links.shape[0])

    return chem_protein_links


def get_proteins(by, items) -> pd.DataFrame:
    """
        读取HerbiV_proteins数据集，返回items中蛋白的信息。
        Read the HerbiV_proteins dataset and return the protein(s) information in items.

        Args:
            by (str): 数据集中与items相匹配的列的列名。Column name of the column in the dataset that matches items.
            items (collections.abc.Iterable): 要查询的蛋白。Protein(s) to be queried.

        Returns:
            pandas.DataFrame: items中蛋白的信息。Protein information in items.
    """

    # 读取HerbiV_proteins数据集
    current_directory = os.path.dirname(os.path.abspath(__file__))
    proteins_all = pd.read_excel(current_directory + r'/Data/Protein.xlsx')

    # 在数据集中获取items中蛋白的信息
    proteins = proteins_all.loc[proteins_all[by].isin(items)].drop_duplicates(subset=['Ensembl_ID'])

    # 重置索引
    proteins.index = range(proteins.shape[0])

    return proteins



def get_SD(by, items) -> pd.DataFrame:
    """
        读取HerbiV_proteins数据集，返回items中蛋白的信息。
        Read the HerbiV_proteins dataset and return the protein(s) information in items.

        Args:
            by (str): 数据集中与items相匹配的列的列名。Column name of the column in the dataset that matches items.
            items (collections.abc.Iterable): 要查询的蛋白。Protein(s) to be queried.

        Returns:
            pandas.DataFrame: items中蛋白的信息。Protein information in items.
    """

    # 读取HerbiV_proteins数据集
    current_directory = os.path.dirname(os.path.abspath(__file__))
    proteins_all = pd.read_excel(current_directory + r'/Data/SD.xlsx')
    proteins = proteins_all.loc[proteins_all[by].isin(items)].copy()
    proteins.index = range(proteins.shape[0])

    return proteins


def get_SD_Formula_links(by, items) -> pd.DataFrame:
    """
        读取HerbiV_formula_tcm_links数据集，返回items中复方/中药的复方-中药连接信息。
        Read the HerbiV_formula_tcm_links dataset
        and return the formula(s)-TCM connection information of formula(s)/TCM in items.

        Args:
            by (str):数据集中与items相匹配的列的列名。Column name of the column in the dataset that matches items.
            items (collections.abc.Iterable): 要查询的复方/中药。Formula(s)/TCM to be queried.

        Returns:
            formula_tcm_links(pandas.DataFrame): items中复方/中药的复方-中药连接信息。
            Formula(s)-TCM connection information of formula(s)/TCM in items.
    """

    formula_tcm_links_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) + r'/Data/SD_Formula_Links.xlsx')
    formula_tcm_links = formula_tcm_links_all.loc[formula_tcm_links_all[by].isin(items)].copy()
    formula_tcm_links.index = range(formula_tcm_links.shape[0])

    return formula_tcm_links


def get_targetNum_dict(symbol_list, interaction_num, PPI_DICT):
    ALL_PPI_PROTEIN = [
        p for symbol in symbol_list
        if symbol in PPI_DICT
        for p in PPI_DICT[symbol]
    ]
    PPI_NUMBER = Counter(ALL_PPI_PROTEIN)
    return {k: v for k, v in PPI_NUMBER.items() if v >= interaction_num}



def get_data(protein_list_path, interaction_num):

    file_name = protein_list_path
    Symbol_list = get_Symbol(file_name)

    print('The number of proteins in the file is: ', len(Symbol_list))

    # step_2.2 获取symbol对应的PPI蛋白 在这里只要跟差异表达蛋白拥有相互作用就进行保留
    # interaction_num = 0 已经存在相互作用的蛋白 最少应该和几个差异表达蛋白列表中的蛋白相互作用,可以设置为参数

    Symbol_PPI_list = get_PPI_Symbol_List(Symbol_list, interaction_num)

    print('The number of PPI proteins list is: ', len(Symbol_PPI_list))
    print('Please wait for a while, the program is running...')

    with open ('Data/Drug/Symbol_To_Target.json', 'r') as f:
        Symbol_To_Target_wm = json.load(f)

    with open ('Data/ID_Transformed/Symbol_To_Fullname.json', 'r') as f:
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
    Symbol = pd.read_excel(file_name)

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
    with open ('Data/PPI/PPI.json', 'r') as f:
        PPI_DICT = json.load(f)

    symbol_list = symbol_list['gene_name'].tolist()

    TARGET_PPI = get_targetNum_dict(symbol_list, interaction_num, PPI_DICT)
    TARGET_PPI_LIST = [i for i in TARGET_PPI.keys() if i not in symbol_list]

    return TARGET_PPI_LIST


if __name__ == '__main__':
    with open ('Data/PPI/PPI.json', 'r') as f:
        PPI_DICT = json.load(f)

    symbol_list = ["ARF5"]
    interaction_num = 0

    TARGET_PPI = get_targetNum_dict(symbol_list, interaction_num, PPI_DICT)
    TARGET_PPI_LIST = [i for i in TARGET_PPI.keys() if i not in symbol_list]

    print(TARGET_PPI_LIST)