import get
import compute
import output
from tqdm import tqdm
import pandas as pd


# TODO: 将文档修改为get中的格式。
def from_tcm_or_formula(tcm_or_formula_id,
                        proteins_id=None,
                        score=990,
                        out_for_cytoscape=True,
                        out_graph=True,
                        re=True,
                        path='results'):
    """
        进行经典的正向网络药理学分析

        Args:
            tcm_or_formula_id: 任何可以使用in判断一个元素是否在其中的组合数据类型，拟分析的中药或复方的ID。
            proteins_id: None 或任何可以使用in判断一个元素是否在其中的组合数据类型，存储拟分析蛋白（靶点）在STITCH中的Ensembl_ID。
                        默认为None
            score (int): HerbiV_chemical_protein_links数据集中仅combined_score大于等于score的记录会被筛选出，默认为990。
            out_for_cytoscape (bool): 是否输出用于Cytoscape绘图的文件，默认为True。
            out_graph (bool): 是否输出基于ECharts的html格式的网络可视化图，默认为True。
            re (bool): 是否返回原始分析结果（中药、化合物（中药成分）、蛋白（靶点）及其连接信息）。
            path (str): 存放结果的目录。


        Returns:
            formula: 复方信息（仅在输入的tcm_or_formula为HVPID时返回）。
            formula_tcm_links: 复方-中药连接信息（仅在输入的tcm_or_formula为HVPID时返回）。
            tcm: 中药信息。
            tcm_chem_links: 中药-化合物（中药成分）连接信息。
            chem: 化合物（中药成分）信息。
            chem_protein_links: 化合物（中药成分）-蛋白（靶点）连接信息。
            proteins: 蛋白（靶点）信息。
    """

    if tcm_or_formula_id[0][2] == 'F':  # 判断输入是否为复方的HVPID
        formula = get.get_formula('DNFID', tcm_or_formula_id)  # 获取该复方的信息
        formula_tcm_links = get.get_formula_tcm_links('DNFID', formula['DNFID'])
        tcm = get.get_tcm('DNHID', formula_tcm_links['DNHID'])
    else:
        formula = None
        formula_tcm_links = None
        tcm = get.get_tcm('DNHID', tcm_or_formula_id)

    tcm_chem_links = get.get_tcm_chem_links('DNHID', tcm['DNHID'])
    chem = get.get_chemicals('DNCID', tcm_chem_links['DNCID'])
    chem_protein_links = get.get_chem_protein_links('DNCID', chem['DNCID'], score)

    if proteins_id is None:
        proteins = get.get_proteins('Ensembl_ID', chem_protein_links['Ensembl_ID'])
    else:
        proteins = get.get_proteins('Ensembl_ID', proteins_id)

    # 深度优先搜索筛选有效节点
    formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links, proteins = dfs_filter(
        formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links, proteins)

    tcm, chem, formula = compute.score(tcm, tcm_chem_links, chem, chem_protein_links, formula, formula_tcm_links)

    if out_for_cytoscape:
        output.out_for_cyto(tcm, tcm_chem_links, chem, chem_protein_links, proteins, path)

    if out_graph:
        output.vis(tcm, tcm_chem_links, chem, chem_protein_links, proteins, path)

    if re:
        if tcm_or_formula_id[0][2] == 'F':
            return formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links, proteins
        else:
            return tcm, tcm_chem_links, chem, chem_protein_links, proteins


def from_proteins(proteins,
                  score=0,
                  random_state=None,
                  num=1000,
                  tcm_component=True,
                  formula_component=True,
                  out_for_cytoscape=True,
                  re=True,
                  path='result'):
    """
    进行逆向网络药理学分析
    """
    # 初始化进度条
    with tqdm(total=10, desc="Overall Progress") as pbar:
        # 获取蛋白质信息
        proteins = get.get_proteins('Ensembl_ID', proteins)
        pbar.update(1)

        # 获取化合物-蛋白质连接信息
        chem_protein_links = get.get_chem_protein_links('Ensembl_ID', proteins['Ensembl_ID'], score)
        pbar.update(1)

        # 异常处理：如果没有找到化合物-蛋白质连接
        if chem_protein_links.empty:
            raise ValueError(
                f"No compound-protein links found based on the set score value (score={score}). "
                f"Please try lowering the score to obtain more results."
            )

        # 获取化合物信息
        chem = get.get_chemicals('DNCID', chem_protein_links['DNCID'])
        pbar.update(1)

        # 获取中药-化合物连接信息
        tcm_chem_links = get.get_tcm_chem_links('DNCID', chem['DNCID'])
        pbar.update(1)

        # 获取中药信息
        tcm = get.get_tcm('DNHID', tcm_chem_links['DNHID'])
        pbar.update(1)

        # 获取复方-中药连接信息
        formula_tcm_links = get.get_formula_tcm_links('DNHID', tcm['DNHID'])
        pbar.update(1)

        # 获取复方信息
        formula = get.get_formula('DNFID', formula_tcm_links['DNFID'])
        pbar.update(1)

        # 深度优先搜索筛选有效节点
        formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links, proteins = dfs_filter(
            formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links, proteins
        )
        pbar.update(1)

        # 计算 Score
        tcm, chem, formula = compute.score(tcm, tcm_chem_links, chem, chem_protein_links, formula, formula_tcm_links)
        pbar.update(1)

        # 调用优化模型
        if tcm_component:
            tcms = compute.component(tcm.loc[tcm['Importance Score'] != 1.0], random_state, num)
            pbar.update(1)
        else:
            tcms = None

        if formula_component:
            formulas = compute.component(formula.loc[formula['Importance Score'] != 1.0], random_state, num)
            pbar.update(1)
        else:
            formulas = None

        # 输出用于 Cytoscape 绘图的文件
        if out_for_cytoscape:
            output.out_for_cyto(tcm, tcm_chem_links, chem, chem_protein_links, proteins, path)
            pbar.update(1)

        # 返回结果
        if re:
            return formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links, proteins, tcms, formulas


def from_tcm_formula(formula_id, out_graph=True, re=True, path='results'):
    """
    进行经典的正向网络药理学分析，并将结果保存到 Excel 中，同时使用 pyecharts 进行可视化。

    Args:
        formula_id: 复方的 ID。
        out_graph (bool): 是否输出基于 pyecharts 的网络可视化图，默认为 True。
        re (bool): 是否返回原始分析结果。
        path (str): 存放结果的目录。
    """
    import os
    if not os.path.exists(path):
        os.makedirs(path)

    # 总进度条，分为 6 个阶段
    with tqdm(total=8, desc="（中医维度）从复方进行分析") as pbar:
        # 阶段 1：获取复方信息
        formula = get.get_tcm_formula('DNFID', formula_id)  # 获取该复方的信息
        pbar.update(1)

        # 阶段 2：获取复方-中药连接信息
        formula_tcm_links = get.get_tcm_formula_tcm_links('DNFID', formula['DNFID'])
        pbar.update(1)

        # 阶段 3：获取中药信息
        tcm = get.get_tcm_tcm('DNHID', formula_tcm_links['DNHID'])
        pbar.update(1)

        # 阶段 4：获取复方-辩证连接信息
        formula_SD_links = get.get_tcm_SD_Formula_links('DNFID', formula['DNFID'])
        pbar.update(1)

        # 阶段 5：获取辩证信息
        SD = get.get_tcm_SD('DNSID', formula_SD_links['DNSID'])
        pbar.update(1)

        Formula_illness_links = get.get_tcm_Formula_illness_links('DNFID', formula['DNFID'])
        pbar.update(1)

        illness = get.get_tcm_illness('DNIID', Formula_illness_links['DNIID'])
        pbar.update(1)

        # 阶段 6：保存结果到 Excel
        formula_df = pd.DataFrame(formula)
        formula_tcm_links_df = pd.DataFrame(formula_tcm_links)
        tcm_df = pd.DataFrame(tcm)
        formula_SD_links_df = pd.DataFrame(formula_SD_links)
        SD_df = pd.DataFrame(SD)
        illness_df = pd.DataFrame(illness)

        if out_graph:
            output.tcm_vis(formula_df, formula_tcm_links_df, formula_SD_links_df, path)

        if re:
            with pd.ExcelWriter(f"{path}/results.xlsx") as writer:
                formula_df.to_excel(writer, sheet_name="复方信息", index=False)
                formula_tcm_links_df.to_excel(writer, sheet_name="复方-中药连接信息", index=False)
                tcm_df.to_excel(writer, sheet_name="中药信息", index=False)
                formula_SD_links_df.to_excel(writer, sheet_name="复方-辩证连接信息", index=False)
                SD_df.to_excel(writer, sheet_name="辩证信息", index=False)
                illness_df.to_excel(writer, sheet_name="疾病信息", index=False)
            pbar.update(1)

    return None

def from_tcm_SD(SD_id, out_graph=True, re=True, path='results'):
    """
    进行经典的正向网络药理学分析，并将结果保存到 Excel 中，同时使用 pyecharts 进行可视化。

    Args:
        SD_id: 复方的 ID。
        out_graph (bool): 是否输出基于 pyecharts 的网络可视化图，默认为 True。
        re (bool): 是否返回原始分析结果。
        path (str): 存放结果的目录。
    """
    import os
    if not os.path.exists(path):
        os.makedirs(path)

    with tqdm(total=8, desc="（中医维度）从辩证进行分析") as pbar:
        SD = get.get_tcm_SD('DNSID', SD_id)
        pbar.update(1)

        formula_SD_links = get.get_tcm_SD_Formula_links('DNSID', SD['DNSID'])
        pbar.update(1)

        formula = get.get_tcm_formula('DNFID', formula_SD_links['DNFID'])  # 获取该复方的信息
        pbar.update(1)

        formula_tcm_links = get.get_tcm_formula_tcm_links('DNFID', formula['DNFID'])
        pbar.update(1)

        SD_illness_links = get.get_tcm_SD_illness_links('DNIID', SD['DNSID'])
        pbar.update(1)

        tcm = get.get_tcm_tcm('DNHID', formula_tcm_links['DNHID'])
        pbar.update(1)

        illness = get.get_tcm_illness('DNIID', SD_illness_links['DNIID'])
        pbar.update(1)

        formula_df = pd.DataFrame(formula)
        formula_tcm_links_df = pd.DataFrame(formula_tcm_links)
        tcm_df = pd.DataFrame(tcm)
        formula_SD_links_df = pd.DataFrame(formula_SD_links)
        SD_df = pd.DataFrame(SD)
        illness_df = pd.DataFrame(illness)

        if out_graph:
                output.tcm_vis(formula_df, formula_tcm_links_df, formula_SD_links_df, path)
                pbar.update(1)

        if re:
            with pd.ExcelWriter(f"{path}/results.xlsx") as writer:
                formula_df.to_excel(writer, sheet_name="复方信息", index=False)
                formula_tcm_links_df.to_excel(writer, sheet_name="复方-中药连接信息", index=False)
                tcm_df.to_excel(writer, sheet_name="中药信息", index=False)
                formula_SD_links_df.to_excel(writer, sheet_name="复方-辩证连接信息", index=False)
                SD_df.to_excel(writer, sheet_name="辩证信息", index=False)
                illness_df.to_excel(writer, sheet_name="疾病信息", index=False)
            return formula_df, formula_tcm_links_df, tcm_df, formula_SD_links_df, SD_df
            pbar.update(1)


def dfs_filter(formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links, proteins):
    """
    深度优先搜索筛选有效节点（优化版）
    """
    formula_id = set()
    tcm_id = set()
    chem_id = set()
    proteins_id = set()

    # ----------------------------------------------
    # 1. 预先生成映射关系，加速查询
    # ----------------------------------------------
    # 确保 proteins 是 DataFrame
    if isinstance(proteins, set):
        raise TypeError("proteins 参数必须是一个 Pandas DataFrame，而不是集合。")

    # 目标蛋白质的集合（快速查找）
    target_proteins = set(proteins['Ensembl_ID'])

    # 复方 -> 中药的映射（如果存在复方-中药连接）
    formula_to_tcms = {}
    if formula_tcm_links is not None:
        for f in formula['DNFID']:
            formula_to_tcms[f] = set(formula_tcm_links.loc[formula_tcm_links['DNFID'] == f, 'DNHID'])

    # 中药 -> 化合物的映射
    tcm_to_chems = tcm_chem_links.groupby('DNHID')['DNCID'].agg(set).to_dict()

    # 化合物 -> 蛋白质的映射
    chem_to_proteins = chem_protein_links.groupby('DNCID')['Ensembl_ID'].agg(set).to_dict()

    # ----------------------------------------------
    # 2. 优化遍历逻辑，避免重复查询 DataFrame
    # ----------------------------------------------
    # 遍历复方（如果存在）或直接遍历所有中药
    if formula_tcm_links is not None:
        # 存在复方-中药连接：遍历每个复方对应的中药
        for f in formula['DNFID']:
            for m in formula_to_tcms.get(f, set()):
                chems = tcm_to_chems.get(m, set())
                for c in chems:
                    proteins_in_chem = chem_to_proteins.get(c, set())
                    valid_proteins = proteins_in_chem & target_proteins  # 直接取交集
                    if valid_proteins:
                        formula_id.add(f)
                        tcm_id.add(m)
                        chem_id.add(c)
                        proteins_id.update(valid_proteins)
    else:
        # 不存在复方-中药连接：直接遍历所有中药
        for m in tcm['DNHID']:
            chems = tcm_to_chems.get(m, set())
            for c in chems:
                proteins_in_chem = chem_to_proteins.get(c, set())
                valid_proteins = proteins_in_chem & target_proteins
                if valid_proteins:
                    tcm_id.add(m)
                    chem_id.add(c)
                    proteins_id.update(valid_proteins)

    # ----------------------------------------------
    # 3. 根据有效节点更新数据
    # ----------------------------------------------
    formula = None if formula is None else formula.loc[formula['DNFID'].isin(formula_id)]
    tcm = tcm.loc[tcm['DNHID'].isin(tcm_id)]
    chem = chem.loc[chem['DNCID'].isin(chem_id)]
    proteins = proteins.loc[proteins['Ensembl_ID'].isin(proteins_id)]

    if formula_tcm_links is not None:
        formula_tcm_links = formula_tcm_links.loc[
            formula_tcm_links['DNFID'].isin(formula_id) &
            formula_tcm_links['DNHID'].isin(tcm_id)
        ]
    tcm_chem_links = tcm_chem_links.loc[
        tcm_chem_links['DNHID'].isin(tcm_id) &
        tcm_chem_links['DNCID'].isin(chem_id)
    ]
    chem_protein_links = chem_protein_links.loc[
        chem_protein_links['DNCID'].isin(chem_id) &
        chem_protein_links['Ensembl_ID'].isin(proteins_id)
    ]

    # 重新索引
    tcm_chem_links.reset_index(drop=True, inplace=True)
    chem_protein_links.reset_index(drop=True, inplace=True)
    proteins.reset_index(drop=True, inplace=True)

    return formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links, proteins

if __name__ == '__main__':
    #from_tcm_formula(['DNF045905'])
    #from_tcm_SD(['DNS001'])
    #from_tcm_or_formula(['DNF109423'], score=990)
    from_proteins(['ENSP00000381588', 'ENSP00000252519'])