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
    """
    # 读取HerbiV_formula数据集
    formula_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) + r'/Data/WM/Formula.xlsx')

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
    """

    # 读取HerbiV_formula_tcm_links数据集
    formula_tcm_links_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) +
                                        r'/Data/WM/Formula_Herb_Links.xlsx')

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
    tcm_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) + r'/Data/WM/Herb.xlsx')

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
    tcm_chem_links_all = pd.read_excel(current_directory + r'/Data/WM/Herb_Chemical_Links.xlsx')

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
    chem_all = pd.read_excel(current_directory + r'/Data/WM/Chemical.xlsx')

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
    chem_protein_links_all = pd.read_excel(current_directory + r'/Data/WM/Chemical_Protein_Links.xlsx')

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
    proteins_all = pd.read_excel(current_directory + r'/Data/WM/Protein.xlsx')

    # 在数据集中获取items中蛋白的信息
    proteins = proteins_all.loc[proteins_all[by].isin(items)].drop_duplicates(subset=['Ensembl_ID'])

    # 重置索引
    proteins.index = range(proteins.shape[0])

    return proteins

def get_tcm_formula(by, items) -> pd.DataFrame:
    """
        读取HerbiV_formula数据集，返回items中复方的信息。
        Read the HerbiV_formula dataset and return the formula(s) information in items.

        Args:
            by (str): 数据集中与items相匹配的列的列名。Column name of the column in the dataset that matches items.
            items (collections.abc.Iterable): 要查询的复方。Formula(s) to be queried.

        Returns:
            formula: items中复方的信息。Formula(s) information in items.
    """
    # 读取HerbiV_formula数据集
    formula_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) + r'/Data/TCM/Formula.xlsx')

    # 在数据集中获取items中复方的信息
    formula = formula_all.loc[formula_all[by].isin(items)].copy()

    # 重新设置索引
    formula.index = range(formula.shape[0])

    return formula


def get_tcm_formula_tcm_links(by, items) -> pd.DataFrame:
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
                                        r'/Data/TCM/Formula_Herb_Links.xlsx')

    # 在数据集中获取items中复方/中药的复方-中药连接信息
    formula_tcm_links = formula_tcm_links_all.loc[formula_tcm_links_all[by].isin(items)].copy()

    # 重新设置索引
    formula_tcm_links.index = range(formula_tcm_links.shape[0])

    return formula_tcm_links


def get_tcm_tcm(by, items) -> pd.DataFrame:
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
    tcm_all = pd.read_excel(os.path.dirname(os.path.abspath(__file__)) + r'/Data/TCM/Herb.xlsx')

    # 在数据集中获取items中中药的信息
    tcm = tcm_all.loc[tcm_all[by].isin(items)].copy()

    # 重新设置索引
    tcm.index = range(tcm.shape[0])

    return tcm


def get_tcm_SD(by, items) -> pd.DataFrame:
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
    proteins_all = pd.read_excel(current_directory + r'/Data/TCM/SD.xlsx')

    # 在数据集中获取items中蛋白的信息
    proteins = proteins_all.loc[proteins_all[by].isin(items)].copy()

    # 重置索引
    proteins.index = range(proteins.shape[0])

    return proteins


def get_tcm_SD_Formula_links(by, items) -> pd.DataFrame:
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
                                        r'/Data/TCM/SD_Formula_Links.xlsx')

    # 在数据集中获取items中复方/中药的复方-中药连接信息
    formula_tcm_links = formula_tcm_links_all.loc[formula_tcm_links_all[by].isin(items)].copy()

    # 重新设置索引
    formula_tcm_links.index = range(formula_tcm_links.shape[0])

    return formula_tcm_links