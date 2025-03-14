import os
import pandas as pd
from pyecharts import options as opts
from pyecharts.charts import Graph


def re_name(tcm, tcm_chem_links, chem, chem_protein_links, protein):
    """
    :param tcm: pd.DataFrame类型，中药信息
    :param tcm_chem_links: pd.DataFrame类型，中药-化合物（中药成分）连接信息
    :param chem: pd.DataFrame类型，化合物（中药成分）信息
    :param chem_protein_links: pd.DataFrame类型，化合物（中药成分）-蛋白质（靶点）连接信息
    :param protein: pd.DataFrame类型，靶点信息
    :return: 返回清洗过的数据
    """
    tcm_c = tcm.copy()
    tcm_chem_links_c = tcm_chem_links.copy()
    chem_c = chem.copy()
    chem_protein_links_c = chem_protein_links.copy()
    protein_c = protein.copy()

    out_chem_protein_links = chem_protein_links_c.iloc[:, 0:2]
    out_chem_protein_links.columns = ['SourceNode', 'TargetNode']

    out_chem_protein_links.loc[:, 'SourceNode'] = out_chem_protein_links.loc[:, 'SourceNode'].apply(
        lambda x: chem_c.loc[chem_c['HVCID'] == x]['Name'].iloc[0] if len(
            chem_c.loc[chem_c['HVCID'] == x]['Name']) > 0 else None)

    out_chem_protein_links.dropna(subset=['SourceNode'], inplace=True)

    out_chem_protein_links.loc[:, 'TargetNode'] = out_chem_protein_links.loc[:, 'TargetNode'].apply(
        lambda x: protein_c.loc[protein_c['Ensembl_ID'] == x]['gene_name'].iloc[0] if len(
            protein_c.loc[protein_c['Ensembl_ID'] == x]['gene_name']) > 0 else None)

    out_chem_protein_links.dropna(subset=['TargetNode'], inplace=True)

    out_tcm_chem = tcm_chem_links_c.iloc[:, 0:2]
    out_tcm_chem.columns = ['SourceNode', 'TargetNode']

    out_tcm_chem.loc[:, 'SourceNode'] = out_tcm_chem.loc[:, 'SourceNode'].apply(
        lambda x: tcm_c.loc[tcm_c['HVMID'] == x]['cn_name'].iloc[0] if len(
            tcm_c.loc[tcm_c['HVMID'] == x]['cn_name']) > 0 else None)

    out_tcm_chem.dropna(subset=['SourceNode'], inplace=True)

    out_tcm_chem.loc[:, 'TargetNode'] = out_tcm_chem.loc[:, 'TargetNode'].apply(
        lambda x: chem_c.loc[chem_c['HVCID'] == x]['Name'].iloc[0] if len(
            chem_c.loc[chem_c['HVCID'] == x]['Name']) > 0 else None)

    out_tcm_chem.dropna(subset=['TargetNode'], inplace=True)

    out_chem = chem_c.loc[:, ['Name']]
    out_chem.columns = ['Key']
    out_chem['Attribute'] = 'Chemicals'

    out_tcm = tcm_c.loc[:, ['cn_name']]
    out_tcm.columns = ['Key']
    out_tcm['Attribute'] = 'TCM'

    out_gene = protein_c.loc[:, ['gene_name']]
    out_gene.columns = ['Key']
    out_gene['Attribute'] = 'Proteins'

    return out_tcm, out_tcm_chem, out_chem, out_chem_protein_links, out_gene


def out_for_cyto(tcm, tcm_chem_links, chem, chem_protein_links, protein, path='result'):
    """
    输出Cytoscape用于作图的网络文件和属性文件
    :param protein:
    :param tcm: pd.DataFrame类型，中药信息
    :param tcm_chem_links: pd.DataFrame类型，中药-化合物（中药成分）连接信息
    :param chem: pd.DataFrame类型，化合物（中药成分）信息
    :param chem_protein_links: pd.DataFrame类型，化合物（中药成分）-蛋白质（靶点）连接信息
    :param path: 字符串类型，存放结果的目录
    """
    # 若无path目录，先创建该目录
    if not os.path.exists(path):
        os.mkdir(path)

    tcm, tcm_chem_links, chem, chem_protein_links, protein = \
        re_name(tcm, tcm_chem_links, chem, chem_protein_links, protein)

    # 输出Network文件
    pd.concat([chem_protein_links, tcm_chem_links]).to_csv(os.path.join(path, 'Network.csv'), index=False)

    # 输出Type文件
    pd.concat([tcm, chem, protein]).to_csv(os.path.join(path, "Type.csv"), index=False)


def vis(tcm, tcm_chem_links, chem, chem_protein_links, protein, path='result'):
    """
    使用pyecharts可视化分析结果
    :param tcm: pd.DataFrame类型，中药信息
    :param tcm_chem_links: pd.DataFrame类型，中药-化合物（中药成分）连接信息
    :param chem: pd.DataFrame类型，化合物（中药成分）信息
    :param chem_protein_links: pd.DataFrame类型，化合物（中药成分）-蛋白质（靶点）连接信息
    :param protein: pd.DataFrame类型，蛋白质（靶点）连接信息
    :param path: 字符串类型，存放结果的目录
    """
    # 若无path目录，先创建该目录
    if not os.path.exists(path):
        os.mkdir(path)

    tcm, tcm_chem_links, chem, chem_protein_links, protein = \
        re_name(tcm, tcm_chem_links, chem, chem_protein_links, protein)

    nodes = []
    edges = []

    categories = [
        {"name": "中药", "color": "#61a0a8"},
        {"name": "化学成分", "color": "#f47920"},
        {"name": "靶点", "color": "#ca8622"},
    ]

    for index, row in tcm_chem_links.iloc[0:].iterrows():
        chinese_medicine = row[0]
        chemical_component = row[1]
        nodes.append({'name': chinese_medicine, "symbolSize": 20, 'category': 0, "color": "#1FA9E9"})
        nodes.append({'name': chemical_component, "symbolSize": 20, 'category': 1, "color": "#FFFF00"})
        edges.append({'source': chinese_medicine, 'target': chemical_component})

    for index, row in chem_protein_links.iloc[0:].iterrows():
        chemical_component = row[0]
        target = row[1]
        nodes.append({'name': chemical_component, "symbolSize": 20, 'category': 1, "color": "#FFFF00"})
        nodes.append({'name': target, "symbolSize": 20, 'category': 2, "color": "#000000"})
        edges.append({'source': chemical_component, 'target': target})

    unique_list = list(set(tuple(item.items()) for item in nodes))
    nodes = [dict(item) for item in unique_list]

    unique_list = list(set(tuple(item.items()) for item in nodes))
    nodes = [dict(item) for item in unique_list]

    Graph(init_opts=opts.InitOpts(width="2400px", height="1200px")) \
        .add(
        '',
        nodes=nodes,
        links=edges,
        categories=categories,
        repulsion=8000,
        layout="circular",
        is_rotate_label=True,
        linestyle_opts=opts.LineStyleOpts(color="source", curve=0.3),
        label_opts=opts.LabelOpts(position="right")
    ) \
        .set_global_opts(
        title_opts=opts.TitleOpts(title=''),
        legend_opts=opts.LegendOpts(orient="vertical", pos_left="2%", pos_top="20%")
    ) \
        .render(path=os.path.join(path, "Graph.html"))


def tcm_vis(formula_df, formula_tcm_links_df, formula_SD_links_df, path):
    nodes = []
    links = []

    # 添加复方节点
    for index, row in formula_df.iterrows():
        formula_id = row['DNFID']
        nodes.append({'name': str(formula_id), 'symbol_size': 40, 'category': 0})

    # 添加中药节点
    for tcm_id in formula_tcm_links_df['DNHID'].tolist():
        nodes.append({'name': str(tcm_id), 'symbol_size': 20, 'category': 1})

    # 添加SD节点
    for sd_id in formula_SD_links_df['DNSID'].tolist():
        nodes.append({'name': str(sd_id), 'symbol_size': 40, 'category': 2})

    # 添加复方-中药边
    for index, row in formula_tcm_links_df.iloc[0:].iterrows():
        formula = row.iloc[0]
        tcm = row.iloc[1]
        links.append({'source': str(formula), 'target': str(tcm)})

    for index, row in formula_SD_links_df.iloc[0:].iterrows():
        SD = row.iloc[0]
        formula = row.iloc[1]
        links.append({'source': str(formula), 'target': str(SD)})

    # 去重节点
    unique_list = list(set(tuple(item.items()) for item in nodes))
    nodes = [dict(item) for item in unique_list]

    # 创建图表
    graph_network = (
        Graph(init_opts=opts.InitOpts(width="2400px", height="1200px"))  # 使用百分比设置宽高
        .add(
            series_name="",
            nodes=nodes,
            links=links,
            categories=[
                {"name": "复方"},  # 对应category 0
                {"name": "中药"},  # 对应category 1
                {"name": "辩证"}  # 对应category 2
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

    graph_circle = (
        Graph(init_opts=opts.InitOpts(width="2400px", height="1200px"))  # 使用百分比设置宽高
        .add(
            series_name="",
            nodes=nodes,
            links=links,
            repulsion=8000,
            layout="circular",
            is_rotate_label=True,
            linestyle_opts=opts.LineStyleOpts(color="source", curve=0.3),
            label_opts=opts.LabelOpts(position="right"),
            categories=[
                {"name": "复方"},  # 对应category 0
                {"name": "中药"},  # 对应category 1
                {"name": "辩证"}  # 对应category 2
            ],
        )
        .set_global_opts(
            title_opts=opts.TitleOpts(title="网络图示例"),  # 图表标题
            tooltip_opts=opts.TooltipOpts(trigger="item", formatter="{b}")  # 鼠标悬停提示
        )
    )

    graph_network.render(f"{path}/graph_network.html")
    graph_circle.render(f"{path}/graph_circle.html")
