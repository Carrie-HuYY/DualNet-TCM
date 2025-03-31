from typing import Union
from tqdm import tqdm
import warnings
import numpy as np
import pandas as pd
from math import ceil
import random
from typing import List, Tuple

warnings.filterwarnings("ignore")
np.seterr(all="ignore")


def calculate_herbiv_score(data: pd.DataFrame, id_col: str, link_data: pd.DataFrame, score_col: str) -> pd.Series:
    """
    计算 HerbiV Score。

    Args:
        data: 包含 ID 列的数据。
        id_col: ID 列的名称。
        link_data: 连接数据，包含 ID 和目标列。
        score_col: 目标列的名称。

    Returns:
        计算后的 HerbiV Score。
    """
    # 使用 groupby 和 apply 优化计算
    def _calculate_score(group):
        return 1 - (1 - group[score_col]).prod()

    scores = link_data.groupby(id_col).apply(_calculate_score)
    return data[id_col].map(scores).fillna(0)  # 如果没有匹配的分数，默认为 0


def score(
    tcm: pd.DataFrame,
    tcm_chem_links: pd.DataFrame,
    chem: pd.DataFrame,
    chem_protein_links: pd.DataFrame,
    formula: Union[pd.DataFrame, None] = None,
    formula_tcm_links: Union[pd.DataFrame, None] = None,
    weights: Union[dict, None] = None
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    计算复方、中药和化合物的 HerbiV Score。

    Args:
        tcm: 中药信息。
        tcm_chem_links: 中药-化合物连接信息。
        chem: 化合物信息。
        chem_protein_links: 化合物-蛋白质连接信息。
        formula: 复方信息。默认为 None。
        formula_tcm_links: 复方-中药连接信息。默认为 None。
        weights: 各靶点的权重，各权重的和应为 1。默认为 None。

    Returns:
        tcm_and_score: 中药信息及 HerbiV Score。
        chem_and_score: 化合物信息及 HerbiV Score。
        formula_and_score: 复方信息及 HerbiV Score。
    """
    # 创建副本以避免修改原始数据
    formula_and_score = formula.copy() if formula is not None else None
    tcm_and_score = tcm.copy()
    chem_and_score = chem.copy()

    # 获取所有蛋白质 ID
    proteins_id = chem_protein_links['Ensembl_ID'].unique()

    # 预计算化合物-蛋白质连接信息
    chem_protein_dict = chem_protein_links.groupby('Ensembl_ID').apply(
        lambda x: x.set_index('DNCID')['Combined_score'].to_dict()
    ).to_dict()

    # 预计算中药-化合物连接信息
    tcm_chem_dict = tcm_chem_links.groupby('DNHID')['DNCID'].apply(list).to_dict()

    # 预计算复方-中药连接信息
    if formula is not None:
        formula_tcm_dict = formula_tcm_links.groupby('DNFID')['DNHID'].apply(list).to_dict()

    # 使用 tqdm 显示蛋白质计算的进度
    for protein in tqdm(proteins_id, desc="Calculating Chemical Scores"):
        # 计算每一个化合物的 HerbiV Score
        chem_scores = chem_protein_dict.get(protein, {})
        chem_and_score[f'{protein} HerbiV Score'] = chem['DNCID'].map(chem_scores).fillna(0).apply(
            lambda x: 1 - (1 - x) if x != 0 else 0
        )

        # 计算各中药的 HerbiV Score
        tcm_scores = tcm['DNHID'].apply(
            lambda x: 1 - (1 - chem_and_score.loc[chem_and_score['DNCID'].isin(tcm_chem_dict.get(x, [])), f'{protein} HerbiV Score']).prod()
        )
        tcm_and_score[f'{protein} HerbiV Score'] = tcm_scores

        # 若传入了复方相关信息，则还需计算各复方的 HerbiV Score
        if formula is not None:
            formula_scores = formula['DNFID'].apply(
                lambda x: 1 - (1 - tcm_and_score.loc[tcm_and_score['DNHID'].isin(formula_tcm_dict.get(x, [])), f'{protein} HerbiV Score']).prod()
            )
            formula_and_score[f'{protein} HerbiV Score'] = formula_scores

    # 设置默认权重
    if weights is None:
        weights = {f'{protein} HerbiV Score': 1 for protein in proteins_id}

    # 加权计算各复方、中药、成分（化合物）的 Importance Score
    def calculate_weighted_score(data: pd.DataFrame) -> pd.Series:
        return (data[list(weights.keys())] * pd.Series(weights)).mean(axis=1)

    if formula is not None:
        formula_and_score['Importance Score'] = calculate_weighted_score(formula_and_score)
    tcm_and_score['Importance Score'] = calculate_weighted_score(tcm_and_score)
    chem_and_score['Importance Score'] = calculate_weighted_score(chem_and_score)

    # 根据 Importance Score 降序排序
    if formula is not None:
        formula_and_score = formula_and_score.sort_values(by='Importance Score', ascending=False).reset_index(drop=True)
    tcm_and_score = tcm_and_score.sort_values(by='Importance Score', ascending=False).reset_index(drop=True)
    chem_and_score = chem_and_score.sort_values(by='Importance Score', ascending=False).reset_index(drop=True)

    return tcm_and_score, chem_and_score, formula_and_score

def component(items_and_score: pd.DataFrame, random_state=None, num=1000, c=10) -> pd.DataFrame:
    """
    生成复方/中药的组合，并计算其重要性评分。

    Args:
        items_and_score: 包含复方/中药信息及重要性评分的 DataFrame。
        random_state: 随机种子。
        num: 需要生成的组合数。
        c: 背包容量。

    Returns:
        包含组合及其重要性评分的 DataFrame。
    """
    if 'DNFID' in items_and_score.columns:
        by = 'DNFID'
        name = 'name'
    else:
        by = 'DNHID'
        name = 'cn_name'

    dps = []
    items_ls = []
    n = len(items_and_score)
    weights = np.ones(n, dtype=int)  # 使用 NumPy 数组替代列表
    names = items_and_score[by].values  # 使用 NumPy 数组替代列表
    values = items_and_score['Importance Score'].values  # 使用 NumPy 数组替代列表
    n = ceil(len(weights) / 10)

    if random_state is not None:
        random.seed(random_state)

    for _ in tqdm(range(num), desc="Generating Components"):
        random_indices = random.sample(range(len(weights)), n)
        sampled_weights = weights[random_indices]
        sampled_names = names[random_indices]
        sampled_values = values[random_indices]

        # 不能再得出之前的解
        dp, items = knapsack(sampled_weights, n, items_ls, sampled_names, sampled_values, c)
        dps.append(dp)
        items_ls.append(items)

    # 用 pd.DataFrame 存储结果
    components = pd.DataFrame({'Importance Score': dps, 'items': items_ls})
    components['items'] = components['items'].apply(lambda x: ';'.join(x))

    # 计算 Score 的提升量
    components['Boost'] = components.apply(boost, axis=1, args=(items_and_score, by))

    # 根据 Boost 降序排序
    components = components.sort_values(by='Boost', ascending=False).reset_index(drop=True)

    return components


def boost(row: pd.Series, items_and_score: pd.DataFrame, by: str) -> float:
    """
    计算组合的重要性评分提升量。

    Args:
        row: 包含组合及其重要性评分的行。
        items_and_score: 包含复方/中药信息及重要性评分的 DataFrame。
        by: 用于匹配的列名。

    Returns:
        重要性评分的提升量。如果无法计算（没有匹配项），返回0。
    """
    ls = row['items'].split(';')
    scores = items_and_score.loc[items_and_score[by].isin(ls), 'Importance Score'].values

    # Handle empty case
    if len(scores) == 0:
        return 0.0

    max_score = max(scores)

    # Avoid division by zero
    if max_score == 0:
        return 0.0

    return (row['Importance Score'] - max_score) / max_score


def knapsack(weights: np.ndarray, n: int, forbidden_combinations: List[List[str]], names: np.ndarray, values: np.ndarray, c: int = 10) -> Tuple[float, List[str]]:
    """
    使用动态规划解决背包问题。

    Args:
        weights: 物品的重量。
        n: 物品数量。
        forbidden_combinations: 禁止的组合列表。
        names: 物品的名称。
        values: 物品的价值。
        c: 背包容量。

    Returns:
        最大价值及选择的物品列表。
    """
    # 创建一个二维数组用于存储计算结果
    dp = np.zeros((n + 1, c + 1), dtype=float)
    # 创建一个二维列表用于记录选择的物品
    items = [[[] for _ in range(c + 1)] for _ in range(n + 1)]

    for i in range(1, n + 1):
        for j in range(1, c + 1):
            if weights[i - 1] <= j:
                # 检查当前是否与禁止组合冲突
                conflict = False
                for combination in forbidden_combinations:
                    if names[i - 1] in combination and any(
                            item in items[i - 1][j - weights[i - 1]] for item in combination):
                        conflict = True
                        break
                if conflict:
                    dp[i][j] = dp[i - 1][j]
                    items[i][j] = items[i - 1][j]
                else:
                    new_value = 1 - (1 - values[i - 1]) * (1 - dp[i - 1][j - weights[i - 1]])
                    if new_value > dp[i - 1][j]:
                        dp[i][j] = new_value
                        items[i][j] = [names[i - 1]] + items[i - 1][j - weights[i - 1]]
                    else:
                        dp[i][j] = dp[i - 1][j]
                        items[i][j] = items[i - 1][j]
            else:
                dp[i][j] = dp[i - 1][j]
                items[i][j] = items[i - 1][j]

    # 计算累计 Score 比例
    score_ratio = np.cumsum(dp[-1]) / np.sum(dp[-1])

    # 最大似然估计
    mle_estimates = (score_ratio - 1 / len(score_ratio)) / np.sqrt(2 / len(score_ratio))

    # 应选择的复方/中药数（即索引）但不能只选择一个
    num_components = np.argmin(mle_estimates) + 1
    num_components = 2 if num_components <= 1 else num_components

    return dp[-1][num_components], items[-1][num_components]