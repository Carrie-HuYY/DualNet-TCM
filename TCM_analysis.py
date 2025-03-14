import TCM_get
import pandas as pd
from tqdm import tqdm
from TCM_output import vis


def from_tcm_formula(formula_id, out_graph=True, re=True, path='results'):
    """
    进行经典的正向网络药理学分析，并将结果保存到 Excel 中，同时使用 pyecharts 进行可视化。

    Args:
        formula_id: 复方的 ID。
        out_graph (bool): 是否输出基于 pyecharts 的网络可视化图，默认为 True。
        re (bool): 是否返回原始分析结果。
        path (str): 存放结果的目录。
    """
    # 创建结果目录
    import os
    if not os.path.exists(path):
        os.makedirs(path)

    # 获取复方信息
    with tqdm(total=1, desc="获取复方信息") as pbar:
        formula = TCM_get.get_formula('DNFID', formula_id)  # 获取该复方的信息
        pbar.update(1)

    # 获取复方-中药连接信息
    with tqdm(total=1, desc="获取复方-中药连接信息") as pbar:
        formula_tcm_links = TCM_get.get_formula_tcm_links('DNFID', formula['DNFID'])
        pbar.update(1)

    # 获取中药信息
    with tqdm(total=1, desc="获取中药信息") as pbar:
        tcm = TCM_get.get_tcm('DNHID', formula_tcm_links['DNHID'])
        pbar.update(1)

    # 获取复方-SD 连接信息
    with tqdm(total=1, desc="获取复方-SD 连接信息") as pbar:
        formula_SD_links = TCM_get.get_SD_Formula_links('DNFID', formula['DNFID'])
        pbar.update(1)

    # 获取 SD 信息
    with tqdm(total=1, desc="获取 SD 信息") as pbar:
        SD = TCM_get.get_SD('DNSID', formula_SD_links['DNSID'])
        pbar.update(1)

    with tqdm(total=1, desc="保存结果到 Excel") as pbar:
        # 将结果整理为 DataFrame
        formula_df = pd.DataFrame(formula)
        formula_tcm_links_df = pd.DataFrame(formula_tcm_links)
        tcm_df = pd.DataFrame(tcm)
        formula_SD_links_df = pd.DataFrame(formula_SD_links)
        SD_df = pd.DataFrame(SD)

        # 保存到 Excel 文件
        with pd.ExcelWriter(f"{path}/results.xlsx") as writer:
            formula_df.to_excel(writer, sheet_name="复方信息", index=False)
            formula_tcm_links_df.to_excel(writer, sheet_name="复方-中药连接信息", index=False)
            tcm_df.to_excel(writer, sheet_name="中药信息", index=False)
            formula_SD_links_df.to_excel(writer, sheet_name="复方-SD 连接信息", index=False)
            SD_df.to_excel(writer, sheet_name="SD 信息", index=False)
        pbar.update(1)

    if out_graph:
        vis(formula_df, formula_tcm_links_df, formula_SD_links_df, path)

if __name__ == '__main__':
    from_tcm_formula(['DNF102324'])