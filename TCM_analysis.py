import TCM_get
import pandas as pd
from tqdm import tqdm
from pyecharts import options as opts
from pyecharts.charts import Graph

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

    with tqdm(total=1, desc="生成可视化图表") as pbar:
        # 构建节点和边
        nodes = []
        links = []

        for index, row in formula_df.iterrows():
            formula_id = row['DNFID']
            nodes.append({'name': str(formula_id), 'symbol_size': 40, 'category': 0})

        for tcm_id in formula_tcm_links_df['DNHID'].tolist():
            nodes.append({'name': str(tcm_id), 'symbol_size': 20, 'category': 1})

        for sd_id in formula_SD_links_df['DNSID'].tolist():
            nodes.append(opts.GraphNode(name=str(sd_id), symbol_size=20, category=2))

        for index, row in formula_tcm_links_df.iloc[0:].iterrows():
            formula = row[0]
            tcm = row[1]
            links.append({'source': str(formula), 'target': str(tcm)})

        for index, row in formula_SD_links_df.iloc[0:].iterrows():
            SD = row[0]
            formula = row[1]
            links.append({'source': str(formula), 'target': str(SD)})

        unique_list = list(set(tuple(item.items()) for item in nodes))
        nodes = [dict(item) for item in unique_list]

        graph = (
            Graph(init_opts=opts.InitOpts(width="1000px", height="600px"))
            .add(
                series_name="",
                nodes=nodes,
                links=links,
                categories=[
                    {"name": "复方"},  # 对应category 0
                    {"name": "中药"},  # 对应category 1
                    {"name": "SD"}  # 对应category 2
                ],
                repulsion=8000,  # 节点之间的斥力
                layout="force",  # 使用力导向布局
                linestyle_opts=opts.LineStyleOpts(),  # 边的样式
                label_opts=opts.LabelOpts(is_show=True, position="inside"),  # 节点标签
            )
            .set_global_opts(
                title_opts=opts.TitleOpts(title="网络图示例"),  # 图表标题
                tooltip_opts=opts.TooltipOpts(trigger="item", formatter="{b}")  # 鼠标悬停提示
            )
        )

        # 保存为 HTML 文件
        graph.render(f"{path}/graph.html")
        pbar.update(1)

if __name__ == '__main__':
    from_tcm_formula(['DNF102324'])