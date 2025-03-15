import json
import pandas as pd
from pyecharts.charts import Sunburst ,Tree ,Bar ,Page
from pyecharts import options as opts
from elasticsearch import Elasticsearch
import os
import textwrap
import PDDR_get


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

        .render('output/' + disease_name + ' reported_number_' + str(
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
        'output/' + disease_name + ' reported_number_' + str(
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

        .render('output/' + disease_name + ' reported_number_' + str(
            reported_number) + '/' + dir_name + "/drug_suggestion.html")
    )


# 生成excel表格
def get_excel(un_relevant_targets_recommend_drug, dir_name, disease_name, reported_number):
    # 生成excel表格
    df = pd.DataFrame(un_relevant_targets_recommend_drug.items(), columns=['Target', 'Drug'])
    df.to_excel('output/' + disease_name + ' reported_number_' + str(
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
         drug_cl_report) = PDDR_get.get_drug_report_info(drug_ap, drug_cl, disease, input, es)

        # 药物热度频率
        drug_not_report = drug_ap_not_report + drug_cl_not_report
        drug_report = drug_ap_report + drug_cl_report
        drug_frequency = PDDR_get.get_drug_frequency(drug_not_report, drug_report, es)
        if drug_frequency:
            os.makedirs(
                'output/' + disease_name + ' reported_number_' + str(reported_number) + '/' + dir_name + '/' + symbol,
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
    df.to_excel('output/' + disease_name + ' reported_number_' + str(
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

if __name__ == '__main__':
    es = Elasticsearch(['https://localhost:9200'])