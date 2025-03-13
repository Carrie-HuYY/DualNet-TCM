import TCM_get
from tqdm import tqdm

# TODO: 将文档修改为get中的格式。
def from_formula(tcm_or_formula_id,
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

    # 获取复方信息并显示进度
    with tqdm(total=1, desc="获取复方信息") as pbar:
        formula = TCM_get.get_formula('DNFID', tcm_or_formula_id)  # 获取该复方的信息
        pbar.update(1)

    # 获取复方-中药连接信息并显示进度
    with tqdm(total=1, desc="获取复方-中药连接信息") as pbar:
        formula_tcm_links = TCM_get.get_formula_tcm_links('DNFID', formula['DNFID'])
        pbar.update(1)

    # 获取中药信息并显示进度
    with tqdm(total=1, desc="获取中药信息") as pbar:
        tcm = TCM_get.get_tcm('DNHID', formula_tcm_links['DNHID'])
        pbar.update(1)

    # 获取复方-SD连接信息并显示进度
    with tqdm(total=1, desc="获取复方-SD连接信息") as pbar:
        formula_SD_links = TCM_get.get_SD_Formula_links('DNFID', formula['DNFID'])
        pbar.update(1)

    # 获取SD信息并显示进度
    with tqdm(total=1, desc="获取SD信息") as pbar:
        SD = TCM_get.get_SD('DNSID', formula_SD_links['DNSID'])
        pbar.update(1)

    print(formula, tcm, SD)

if __name__ == '__main__':
    from_formula(['DNF110222'])