import json
import os
import get
import textwrap
import pandas as pd
from pyecharts import options as opts
from pyecharts.charts import Sunburst ,Tree ,Bar ,Page, Graph
from elasticsearch import Elasticsearch


def all_targets_tree(fda_unre, fda_re, clinical_unre, clinical_re, dir_name):

    with open('config.json', 'r') as f:
        config = json.load(f)

    disease_name = config['disease_name']
    reported_number = config['reported_number']

    fda_all_nmb = len(fda_unre) + len(fda_re)
    cli_all_nmb = len(clinical_unre) + len(clinical_re)
    fda_unre_dict = [{'name': i} for i in fda_unre]
    fda_re_dict = [{'name': i} for i in fda_re]
    clinical_unre_dict = [{'name': i} for i in clinical_unre]
    clinical_re_dict = [{'name': i} for i in clinical_re]
    data = [
        {
            "children": [

                {
                    "children": [
                        {
                            "children": fda_unre_dict,
                            "name": "Not Reported",
                            'value': len(fda_unre)
                        },
                        {
                            "children": fda_re_dict,
                            "name": "Reported",
                            'value': len(fda_re)
                        }],
                    "name": "FDA Approved",
                    'value': fda_all_nmb,
                },
                {
                    "children": [
                        {
                            "children": clinical_unre_dict,
                            "name": "Not Reported",
                            'value': len(clinical_unre)
                        },
                        {
                            "children": clinical_re_dict,
                            "name": "Reported",
                            'value': len(clinical_re)
                        }
                    ],
                    "name": "Clinical",
                    'value': cli_all_nmb
                },
            ],
            "name": "Target",
            'value': fda_all_nmb + cli_all_nmb
        }
    ]

    c = (
        Tree(init_opts=opts.InitOpts(width="1200px", height="800px", renderer="svg"))
        .add("", data, collapse_interval=2,
             symbol_size=10,
             symbol="emptyCircle",
             #leaves_label_opts=opts.LabelOpts(position="right"),
             itemstyle_opts=opts.ItemStyleOpts(border_width=1, border_color="#48466d"),
             )
        .set_global_opts(title_opts=opts.TitleOpts(title="All Targets"))
        .set_series_opts(label_opts=opts.LabelOpts(
            font_size=20,
            font_weight='bold',
            color="#48466d"
        ),
        )

        .render('results/' + disease_name + ' reported_number_' + str(
            reported_number) + '/' + dir_name + "/Targets_tree.html")
    )



# 通过phase将药物进行分类
def drug_classify(target_name):
    with open('data\Drug\Target_To_Drug.json', 'r') as f:
        Target_To_Drug = json.load(f)

    drug_phase = {'Approved': [], 'Clinical_trial': [], 'Others': []}
    drug_ap_cl, drug_ap, drug_cl = [], [], []
    for drug in Target_To_Drug[target_name]:
        value = [*drug.values()][0]
        key = [*drug.keys()][0]
        if value == 'Approved':
            drug_phase['Approved'].append({
                'name': key
            })
            drug_ap_cl.append(key)
            drug_ap.append(key)
        elif value.startswith('Phase') or value.startswith('Clinical'):
            drug_phase['Clinical_trial'].append({
                'name': key
            })
            drug_ap_cl.append(key)
            drug_cl.append(key)
        else:
            drug_phase['Others'].append({
                'name': key
            })
    return drug_phase, drug_ap_cl, drug_ap, drug_cl



# 将字符串处理，过长的字符串换行
def wrap_text(text, max_length=20):
    if len(text) > max_length:
        return textwrap.fill(text, max_length)
    return text


# 将药物转化为tree图需要的数据类型
def drug_treetype_data(drugs):
    drug = []
    for i in drugs:
        name = wrap_text(i)
        drug.append({'name': name})
    return drug


# 生成每个靶标对应药物的树状图以及柱状图
def target_tree_bar(dir_name, symbol, drug_frequency,
                    drug_ap_not_report, drug_ap_report,
                    drug_cl_not_report, drug_cl_report,
                    disease_name, reported_number):
    drug_ap_cl = drug_ap_not_report + drug_ap_report + drug_cl_not_report + drug_cl_report
    drug_not_report = drug_ap_not_report + drug_cl_not_report
    drug_report = drug_ap_report + drug_cl_report
    drug_ap_not_report = drug_treetype_data(drug_ap_not_report)
    drug_ap_report = drug_treetype_data(drug_ap_report)
    drug_cl_not_report = drug_treetype_data(drug_cl_not_report)
    drug_cl_report = drug_treetype_data(drug_cl_report)
    tree_data = [
        {
            "children": [

                {
                    "children": [
                        {
                            "children": drug_ap_not_report,
                            "name": "Not Reported",
                            'value': len(drug_ap_not_report)
                        },
                        {
                            "children": drug_ap_report,
                            "name": "Reported",
                            'value': len(drug_ap_report)
                        }],
                    "name": "FDA Approved",
                    'value': len(drug_ap_not_report + drug_ap_report),
                },
                {
                    "children": [
                        {
                            "children": drug_cl_not_report,
                            "name": "Not Reported",
                            'value': len(drug_cl_not_report)
                        },
                        {
                            "children": drug_cl_report,
                            "name": "Reported",
                            'value': len(drug_cl_report)
                        }
                    ],
                    "name": "Clinical",
                    'value': len(drug_cl_report + drug_cl_not_report)
                },
            ],
            "name": symbol,
            'value': len(drug_ap_cl)
        }
    ]

    tree = (
        Tree(init_opts=opts.InitOpts(width="1600px", height="800px", renderer="svg"))
        .add("", tree_data, collapse_interval=2,
             symbol_size=10,
             symbol="emptyCircle",
             #leaves_label_opts=opts.LabelOpts(position="right"),
             itemstyle_opts=opts.ItemStyleOpts(border_width=1, border_color="#48466d"),
             edge_fork_position="100%",
             )
        .set_global_opts(title_opts=opts.TitleOpts(title=symbol + ' corresponding drugs'))
        .set_series_opts(label_opts=opts.LabelOpts(
            font_size=15,
            font_weight='bold',
            color="#48466d"
        ),
        )
    )
    if drug_not_report == []:
        drug_data = drug_report
    else:
        drug_data = drug_not_report

    bar = (
        Bar(init_opts=opts.InitOpts(width="2000px", height="800px", renderer="svg"))
        .add_xaxis(drug_data)
        .add_yaxis("Drug Frequency", drug_frequency, color='#617bdb')
        .set_global_opts(
            xaxis_opts=opts.AxisOpts(
                is_show=False,
                axislabel_opts=opts.LabelOpts(font_size=15,
                                              font_weight='bold',
                                              color="#48466d",
                                              ),
            ),
            yaxis_opts=opts.AxisOpts(

                axislabel_opts=opts.LabelOpts(font_size=15,
                                              font_weight='bold',
                                              color="#48466d"),

            ),
            legend_opts=opts.LegendOpts(is_show=False),

        )
        .set_series_opts(
            label_opts=opts.LabelOpts(position="right",
                                      color="#48466d",
                                      font_size=15,
                                      font_weight='bold', ),
        )
        .reversal_axis()
    )

    (
        Page()
        .add(tree, bar)
    ).render(
        'results/' + disease_name + ' reported_number_' + str(
            reported_number) + '/' + dir_name + '/' + symbol + '/' + symbol + '.html'
    )


# 制作sunburst图
def get_sunburst(un_relevant_targets_recommend_drug, fa, dir_name):

    with open('config.json', 'r') as f:
        config = json.load(f)

    disease_name = config['disease_name']
    reported_number = config['reported_number']

    data2 = [
        {
            'name': 'FDA approve',
            "itemStyle": {"color": '#fac858'},
            'children': []

        },
        {
            'name': 'Clinical trial',
            "itemStyle": {"color": '#73c0de'},
            'children': []

        },
    ]

    for key, value in un_relevant_targets_recommend_drug.items():
        if key in fa:

            children = {
                "name": key,
                'value': 1,
                "itemStyle": {"color": '#fac858'},
                'children': [
                    {'name': value,
                     'value': 1,
                     "itemStyle": {"color": '#fac858'},
                     }
                ]
            }
            data2[0]['children'].append(children)
        else:
            children = {
                "name": key,
                'value': 1,
                "itemStyle": {"color": '#73c0de'},
                'children': [
                    {'name': value,
                     'value': 1,
                     "itemStyle": {"color": '#73c0de'}
                     }
                ]
            }
            data2[1]['children'].append(children)

    c = (
        Sunburst(init_opts=opts.InitOpts(width="1200px", height="1200px", renderer="svg"))
        .add(
            "",
            data_pair=data2,
            highlight_policy="ancestor",
            sort_="null",
            radius=[0, "95%"],
            # center=["55%", "55%"],
            # 居中

            levels=[
                {},
                {
                    "r0": "15%",
                    "r": "35%",
                    "itemStyle": {"borderWidth": 2},
                    "label": {"rotate": "tangential", },
                },
                {"r0": "35%", "r": "60%", "label": {"align": "right"}},
                {
                    "r0": "60%",
                    "r": "62%",
                    "label": {"position": "outside", "padding": 3, "silent": False},
                    "itemStyle": {"borderWidth": 3},
                },
            ],
        )
        .set_global_opts(title_opts=opts.TitleOpts(title="Suggestions for drug targets",
                                                   pos_left='center',
                                                   ))
        .set_series_opts(label_opts=opts.LabelOpts(formatter="{b}",
                                                   color='black',
                                                   font_weight='bold',
                                                   font_size=15,
                                                   font_family='Microsoft YaHei',
                                                   ))

        .render('results/' + disease_name + ' reported_number_' + str(
            reported_number) + '/' + dir_name + "/drug_suggestion.html")
    )


# 生成excel表格
def get_excel(un_relevant_targets_recommend_drug, dir_name, disease_name, reported_number):
    # 生成excel表格
    df = pd.DataFrame(un_relevant_targets_recommend_drug.items(), columns=['Target', 'Drug'])
    df.to_excel('results/' + disease_name + ' reported_number_' + str(
        reported_number) + '/' + dir_name + "/drug_suggestion.html", index=False)


# 生成靶标对应药物的sunburst图和每个靶标对应的药物信息
def get_sunburst_tree_bar(dir_name, fda_no_review, ct_no_review, fa, disease, input, Symbol_To_Target, es):

    with open('config.json', 'r') as f:
        config = json.load(f)

    disease_name = config['disease_name']
    reported_number = config['reported_number']

    target_not_report = fda_no_review + ct_no_review
    un_relevant_targets_recommend_drug = {}
    for symbol in target_not_report:
        target = [*Symbol_To_Target[symbol].keys()][0]
        # print(symbol)
        # 一个靶标对应的药物信息
        drug_phase, drug_ap_cl, drug_ap, drug_cl = drug_classify(target)

        (drug_ap_not_report,
         drug_ap_report,
         drug_cl_not_report,  # 需要提供两个关键词 1.疾病的名字 2.命中的数量
         drug_cl_report) = get.get_drug_report_info(drug_ap, drug_cl, disease, input, es)

        # 药物热度频率
        drug_not_report = drug_ap_not_report + drug_cl_not_report
        drug_report = drug_ap_report + drug_cl_report
        drug_frequency = get.get_drug_frequency(drug_not_report, drug_report, es)
        if drug_frequency:
            os.makedirs(
                'results/' + disease_name + ' reported_number_' + str(reported_number) + '/' + dir_name + '/' + symbol,
                exist_ok=True)
            target_tree_bar(dir_name, symbol, drug_frequency,
                            drug_ap_not_report, drug_ap_report,
                            drug_cl_not_report, drug_cl_report,
                            disease_name, reported_number)

            number_index = drug_frequency.index(max(drug_frequency))
            if drug_not_report == []:
                suggest_drug = drug_report[number_index]
            else:
                suggest_drug = drug_not_report[number_index]
            un_relevant_targets_recommend_drug[symbol] = suggest_drug

    # 输出为excel文件
    df = pd.DataFrame(un_relevant_targets_recommend_drug.items(), columns=['Target', 'Recommend Drug'])
    df.to_excel('results/' + disease_name + ' reported_number_' + str(
        reported_number) + '/' + dir_name + "/drug_suggestion.xlsx", index=False)
    get_sunburst(un_relevant_targets_recommend_drug, fa, dir_name)


# 对推荐靶标数量进行控制
def sort_targets(no_review, target_max_number, es):
    # 将靶标和对应文献数量做成列表
    sort_list = []
    for target in no_review:
        res = es.search(index="abstract22", body={"query": {"match": {"abstract": target}}}, scroll='5m')
        target_hot = res['hits']['total']['value']
        es.clear_scroll(scroll_id=res['_scroll_id'])
        sort_list.append([target, target_hot])

    # 使用sorted函数对列表进行排序
    sort_list = sorted(sort_list, key=lambda x: x[1], reverse=True)
    sort_list = [x[0] for x in sort_list]

    # 判断推荐的靶标数量是否大于靶标推荐最大值
    if len(sort_list) > target_max_number:
        sort_list = sort_list[:target_max_number]
    else:
        pass
    return sort_list


# 生成新的靶标列表
def new_targets_list(list, sort_list):
    new_list = [x for x in list if x in sort_list]
    return new_list


def re_name(SD, SD_Formula_Links, formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links, protein):
    """
    清洗和重命名数据。

    :param tcm: pd.DataFrame, 中药信息
    :param tcm_chem_links: pd.DataFrame, 中药-化合物连接信息
    :param chem: pd.DataFrame, 化合物信息
    :param chem_protein_links: pd.DataFrame, 化合物-蛋白质连接信息
    :param protein: pd.DataFrame, 蛋白质信息
    :return: 返回清洗后的数据
    """
    # 创建副本以避免修改原始数据
    sd_c = SD.copy()
    sd_formula_links_c = SD_Formula_Links.copy()
    formula_C = formula.copy()
    formula_tcm_links_c = formula_tcm_links.copy()
    tcm_c = tcm.copy()
    tcm_chem_links_c = tcm_chem_links.copy()
    chem_c = chem.copy()
    chem_protein_links_c = chem_protein_links.copy()
    protein_c = protein.copy()

    # 清洗辩证-方剂连接信息
    out_sd_formula_links = sd_formula_links_c.iloc[:, 0:2].rename(
        columns={sd_formula_links_c.columns[0]: 'SourceNode', sd_formula_links_c.columns[1]: 'TargetNode'}
    )

    out_sd_formula_links['SourceNode'] = out_sd_formula_links['SourceNode'].apply(
        lambda x: sd_c.loc[sd_c['DNSID'] == x, '证候'].iloc[0] if not sd_c.loc[sd_c['DNSID'] == x, '证候'].empty else None
    )
    out_sd_formula_links.dropna(subset=['SourceNode'], inplace=True)

    out_sd_formula_links['TargetNode'] = out_sd_formula_links['TargetNode'].apply(
        lambda x:formula_C.loc[formula_C['DNFID'] == x, 'name'].iloc[0] if not formula_C.loc[formula_C['DNFID'] == x, 'name'].empty else None
    )

    # 清洗方剂-中药连接信息
    out_formula_tcm_links = formula_tcm_links_c.iloc[:, 0:2].rename(
        columns={formula_tcm_links.columns[0]: 'SourceNode', formula_tcm_links.columns[1]: 'TargetNode'}
    )
    out_formula_tcm_links['SourceNode'] = out_formula_tcm_links['SourceNode'].apply(
        lambda x: formula_C.loc[formula_C['DNFID'] == x, 'name'].iloc[0] if not formula_C.loc[formula_C['DNFID'] == x, 'name'].empty else None
    )
    out_formula_tcm_links.dropna(subset=['SourceNode'], inplace=True)

    out_formula_tcm_links['TargetNode'] = out_formula_tcm_links['TargetNode'].apply(
        lambda x: tcm_c.loc[tcm_c['DNHID'] == x, 'cn_name'].iloc[0] if not tcm_c.loc[tcm_c['DNHID'] == x, 'cn_name'].empty else None
    )
    out_formula_tcm_links.dropna(subset=['TargetNode'], inplace=True)

    # 清洗中药-化合物连接信息
    out_tcm_chem = tcm_chem_links_c.iloc[:, 0:2].rename(
        columns={tcm_chem_links_c.columns[0]: 'SourceNode', tcm_chem_links_c.columns[1]: 'TargetNode'}
    )
    out_tcm_chem['SourceNode'] = out_tcm_chem['SourceNode'].apply(
        lambda x: tcm_c.loc[tcm_c['DNHID'] == x, 'cn_name'].iloc[0] if not tcm_c.loc[tcm_c['DNHID'] == x, 'cn_name'].empty else None
    )
    out_tcm_chem.dropna(subset=['SourceNode'], inplace=True)

    out_tcm_chem['TargetNode'] = out_tcm_chem['TargetNode'].apply(
        lambda x: chem_c.loc[chem_c['DNCID'] == x, 'Name'].iloc[0] if not chem_c.loc[chem_c['DNCID'] == x, 'Name'].empty else None
    )
    out_tcm_chem.dropna(subset=['TargetNode'], inplace=True)

    # 清洗化合物-蛋白质连接信息
    out_chem_protein_links = chem_protein_links_c.iloc[:, 0:2].rename(
        columns={chem_protein_links_c.columns[0]: 'SourceNode', chem_protein_links_c.columns[1]: 'TargetNode'}
    )
    out_chem_protein_links['SourceNode'] = out_chem_protein_links['SourceNode'].apply(
        lambda x: chem_c.loc[chem_c['DNCID'] == x, 'Name'].iloc[0] if not chem_c.loc[chem_c['DNCID'] == x, 'Name'].empty else None
    )
    out_chem_protein_links.dropna(subset=['SourceNode'], inplace=True)

    out_chem_protein_links['TargetNode'] = out_chem_protein_links['TargetNode'].apply(
        lambda x: protein_c.loc[protein_c['Ensembl_ID'] == x, 'gene_name'].iloc[0] if not protein_c.loc[protein_c['Ensembl_ID'] == x, 'gene_name'].empty else None
    )
    out_chem_protein_links.dropna(subset=['TargetNode'], inplace=True)

    # 清洗辩证信息
    out_sd = sd_c[['证候']].rename(columns={'证候': 'Key'})
    out_sd['Attribute'] = 'SD'

    # 清洗方剂信息
    out_formula = formula_C[['name']].rename(columns={'name': 'Key'})
    out_formula['Attribute'] = 'Formula'

    # 清洗化合物信息
    out_chem = chem_c[['Name']].rename(columns={'Name': 'Key'})
    out_chem['Attribute'] = 'Chemicals'

    # 清洗中药信息
    out_tcm = tcm_c[['cn_name']].rename(columns={'cn_name': 'Key'})
    out_tcm['Attribute'] = 'TCM'

    # 清洗蛋白质信息
    out_gene = protein_c[['gene_name']].rename(columns={'gene_name': 'Key'})
    out_gene['Attribute'] = 'Proteins'

    return out_sd, out_sd_formula_links, out_formula, out_formula_tcm_links, out_tcm, out_tcm_chem, out_chem, out_chem_protein_links, out_gene


def out_for_cyto(SD,
                 SD_Formula_Links,
                 formula,
                 formula_tcm_links,
                 tcm,
                 tcm_chem_links,
                 chem,
                 chem_protein_links,
                 protein,
                 path='results'):
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

    SD, SD_Formula_Links, formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links, protein = \
        re_name(SD, SD_Formula_Links, formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links, protein)

    # 输出Network文件
    pd.concat([SD_Formula_Links, formula_tcm_links, tcm_chem_links, chem_protein_links, tcm_chem_links]).to_csv(os.path.join(path, 'Network.csv'), index=False)

    # 输出Type文件
    pd.concat([SD, formula, tcm, chem, protein]).to_csv(os.path.join(path, "Type.csv"), index=False)


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
