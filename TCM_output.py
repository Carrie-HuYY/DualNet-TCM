from tqdm import tqdm
import pandas as pd
from pyecharts import options as opts
from pyecharts.charts import Graph

def vis(formula_df, formula_tcm_links_df, formula_SD_links_df, path):
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