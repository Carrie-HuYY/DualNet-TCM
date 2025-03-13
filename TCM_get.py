import pandas as pd
import os


# TODO: 为各函数增加抛出异常功能，若无法查询到相关信息，则抛出异常。

def get_formula(by, items) -> pd.DataFrame:
    """
        读取HerbiV_formula数据集，返回items中复方的信息。
        Read the HerbiV_formula dataset and return the formula(s) information in items.

        Args:
            by (str): 数据集中与items相匹配的列的列名。Column name of the column in the dataset that matches items.
            items (collections.abc.Iterable): 要查询的复方。Formula(s) to be queried.

        Returns:
            formula: items中复方的信息。Formula(s) information in items.

        Examples:
            >>> get_formula('HVPID', ['HVP1625'])# 获取HVPID为HVP1625的复方（小柴胡汤）的信息
                 HVPID  ... Source Document
            0  HVP1625  ...   shang han lun
            [1 rows x 6 columns]
    """
    # 读取HerbiV_formula数据集
    formula_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) + r'/Data/TCM/Formula.xlsx')

    # 在数据集中获取items中复方的信息
    formula = formula_all.loc[formula_all[by].isin(items)].copy()

    # 重新设置索引
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

        Examples:
            >>> get_formula_tcm_links('HVPID', ['HVP1625'])# 获取HVPID为HVP1625的复方（小柴胡汤）的复方-中药连接信息
                 HVPID    HVMID
            0  HVP1625  HVM0367
            1  HVP1625  HVM0735
            2  HVP1625  HVM0766
            3  HVP1625  HVM1695
            4  HVP1625  HVM3203
            5  HVP1625  HVM4463
    """

    # 读取HerbiV_formula_tcm_links数据集
    formula_tcm_links_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) +
                                        r'/Data/TCM/Formula_Herb_Links.xlsx')

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

        Examples:
            >>> get_tcm('cn_name', ['柴胡', '黄芩'])# 获取cn_name（中文名）为柴胡和黄芩的中药的信息（不建议使用中文名检索）
                 HVMID cn_name pinyin_name  ... TCM_ID_id SymMap_id TCMSP_id
            0  HVM0367      柴胡     CHAI HU  ...    3396.0      58.0     80.0
            1  HVM1695      黄芩   HUANG QIN  ...    6700.0     188.0    371.0
            [2 rows x 19 columns]
    """

    # 读取HerbiV_tcm数据集
    tcm_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) + r'/Data/TCM/Herb.xlsx')

    # 在数据集中获取items中中药的信息
    tcm = tcm_all.loc[tcm_all[by].isin(items)].copy()

    # 重新设置索引
    tcm.index = range(tcm.shape[0])

    return tcm


def get_SD(by, items) -> pd.DataFrame:
    """
        读取HerbiV_proteins数据集，返回items中蛋白的信息。
        Read the HerbiV_proteins dataset and return the protein(s) information in items.

        Args:
            by (str): 数据集中与items相匹配的列的列名。Column name of the column in the dataset that matches items.
            items (collections.abc.Iterable): 要查询的蛋白。Protein(s) to be queried.

        Returns:
            pandas.DataFrame: items中蛋白的信息。Protein information in items.

        Examples:
            >>> get_proteins('protein_name', ['PDCD1 PD1'])# 获取gene_name（基因名）为PDCD1 PD1的蛋白的信息（不建议使用名称检索）
                    Ensembl_ID  ...  gene_name
            0  ENSP00000335062  ...  PDCD1 PD1
    """

    # 读取HerbiV_proteins数据集
    current_directory = os.path.dirname(os.path.abspath(__file__))
    proteins_all = pd.read_excel(current_directory + r'/Data/TCM/SD.xlsx')

    # 在数据集中获取items中蛋白的信息
    proteins = proteins_all.loc[proteins_all[by].isin(items)].copy()

    # 重置索引
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

        Examples:
            >>> get_formula_tcm_links('HVPID', ['HVP1625'])# 获取HVPID为HVP1625的复方（小柴胡汤）的复方-中药连接信息
                 HVPID    HVMID
            0  HVP1625  HVM0367
            1  HVP1625  HVM0735
            2  HVP1625  HVM0766
            3  HVP1625  HVM1695
            4  HVP1625  HVM3203
            5  HVP1625  HVM4463
    """

    # 读取HerbiV_formula_tcm_links数据集
    formula_tcm_links_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) +
                                        r'/Data/TCM/SD_Formula_Links.xlsx')

    # 在数据集中获取items中复方/中药的复方-中药连接信息
    formula_tcm_links = formula_tcm_links_all.loc[formula_tcm_links_all[by].isin(items)].copy()

    # 重新设置索引
    formula_tcm_links.index = range(formula_tcm_links.shape[0])

    return formula_tcm_links

