import output
import get
import json
from bs4 import BeautifulSoup
import os
import time
from elasticsearch import Elasticsearch


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


def classify_targets_wm(Symbol_To_Target, Symbol_list):
    """
    classify_targets: 根据给定的 Symbol_To_Target 字典和 Symbol_list 列表，对靶标（targets）进行西药分类。

    :param Symbol_To_Target:
    :param Symbol_list:
    :return:
        有药物的靶标（target_have_drug）：这些靶标在 Symbol_To_Target 字典中有对应的条目，即存在与之相关的药物信息。
        无药物的靶标（target_no_drug）：这些靶标在 Symbol_To_Target 字典中没有对应的条目，即没有与之相关的药物信息。
        FDA批准的靶标（target_FDA_approved）：这些靶标不仅有药物，而且这些药物已经成功通过FDA审批。
        临床试验中的靶标（target_clinical_trial）：这些靶标有药物，且这些药物目前处于临床试验阶段。
        其他靶标（target_others）：这些靶标有药物，但药物既不是FDA批准的，也不在临床试验中。
    """
    target_have_drug, target_no_drug = [], []
    target_FDA_approved, target_clinical_trial, target_others = [], [], []

    for symbol in Symbol_list:
        if symbol in Symbol_To_Target.keys():
            target_have_drug.append(symbol)
        else:
            target_no_drug.append(symbol)

    for symbol in target_have_drug:
        target_phage = [*Symbol_To_Target[symbol].values()][0]
        target_name = [*Symbol_To_Target[symbol].keys()][0]
        drug_phase, drug_ap_cl, drug_ap, drug_cl = output.drug_classify(target_name)
        # 增加一个判断，如果symbol所在阶段确实存在对应的药物，则将其加入到对应的列表中
        if target_phage == 'Successful target' and drug_ap != []:
            target_FDA_approved.append(symbol)
        elif target_phage == 'Clinical Trial target' and drug_cl != []:
            target_clinical_trial.append(symbol)
        else:
            target_others.append(symbol)

    return target_have_drug, target_no_drug, target_FDA_approved, target_clinical_trial, target_others


def classify_targets_html(target_have_drug, target_no_drug, target_FDA_approved,
                          target_clinical_trial, target_others):
    """
    将html中的数据进行更改后输出
    :param target_have_drug:
    :param target_no_drug:
    :param target_FDA_approved:
    :param target_clinical_trial:
    :param target_others:
    :return:
    """
    with open('config.json', 'r') as f:
        config = json.load(f)

    disease_name = config['disease_name']
    reported_number = config['reported_number']

    text_html = open(r'Template/target_pie_template.html',
                     'r', encoding='utf-8').read()
    text_html = text_html.replace(
        'Compound data', str(len(target_have_drug))).replace(
        'No-drug data', str(len(target_no_drug))).replace(
        'FDA Approved data', str(len(target_FDA_approved))).replace(
        'Others data', str(len(target_others))).replace(
        'Clinical data', str(len(target_clinical_trial)))

    soup = BeautifulSoup(text_html, 'html.parser')
    with open('results/' + disease_name + '/Targets_pie_chart.html', 'w', encoding='utf-8') as fp:
        fp.write(str(soup))


def query_target(symbol, Symbol_To_PubMedID, Symbol_To_UniprotID, Symbol_To_Fullname, es, keywords):
    pubMedId = Symbol_To_PubMedID[symbol]

    sql1 = {
        'query': {
            'bool': {
                'must': [
                    {
                        'terms': {
                            'pubMedId': pubMedId
                        }
                    },
                    {
                        "match_phrase": {
                            "abstract": keywords
                        }
                    },
                    {
                        "match_phrase": {  # abstract中还要存在另一个关键词
                            "abstract": symbol
                        }
                    },
                ]
            }
        }
    }
    res = es.search(index='abstract22', body=sql1, scroll='5m')
    reported_number_1 = res['hits']['total']['value']
    es.clear_scroll(scroll_id=res['_scroll_id'])

    if reported_number_1 == 0:
        if symbol in Symbol_To_Fullname.keys():
            fullName = Symbol_To_Fullname[symbol]
            sql2 = {
                "query": {
                    "bool": {
                        'must': [
                            {
                                "match_phrase": {
                                    "abstract": keywords
                                }
                            },
                            {
                                "bool": {
                                    "should": [
                                        {
                                            "match_phrase": {
                                                "abstract": fullName
                                            }
                                        },
                                        {
                                            "match_phrase": {
                                                "abstract": symbol
                                            }
                                        }
                                    ]
                                }
                            }
                        ]
                    }
                }
            }
            res = es.search(index='abstract22', body=sql2, scroll='5m')
            reported_number_2 = res['hits']['total']['value']
            es.clear_scroll(scroll_id=res['_scroll_id'])
        else:
            reported_number_2 = 0
    else:
        reported_number_2 = 0

    return reported_number_1 + reported_number_2


# 通过摘要中的关键词进行查询，将靶标分为对于该疾病报道过的靶标和没有报道过的靶标
def report_info(fa, ct, es, keywords, input_num):
    with open('Data/ID_Transformed/Symbol_To_PubMedID.json', 'r') as f:
        Symbol_To_PubMedID = json.load(f)
    with open('Data/ID_Transformed/Symbol_To_UniprotID.json', 'r') as f:
        Symbol_To_UniprotID = json.load(f)
    with open('Data/ID_Transformed/Symbol_To_Fullname.json', 'r') as f:
        Symbol_To_Fullname = json.load(f)
    fda_no_review, fda_review, ct_no_review, ct_review = [], [], [], []
    for symbol in fa:
        # 存在有的symbol没有对应的uniprotID或者pubMedID 对于这样的symbol进行剔除
        if symbol in Symbol_To_PubMedID.keys() and symbol in Symbol_To_UniprotID.keys():
            query_num = query_target(symbol, Symbol_To_PubMedID, Symbol_To_UniprotID, Symbol_To_Fullname, es, keywords)
            if query_num > input_num:
                fda_review.append(symbol)
            else:
                fda_no_review.append(symbol)
    for symbol in ct:
        if symbol in Symbol_To_PubMedID.keys() and symbol in Symbol_To_UniprotID.keys():
            query_num = query_target(symbol, Symbol_To_PubMedID, Symbol_To_UniprotID, Symbol_To_Fullname, es, keywords)
            if query_num > input_num:
                ct_review.append(symbol)
            else:
                ct_no_review.append(symbol)
    return fda_no_review, fda_review, ct_no_review, ct_review


def set_config_auto():
    with open('config.json', 'r') as f:
        config = json.load(f)

    disease_name = config['disease_name']
    reported_number = config['reported_number']
    interaction_num = config['interaction_num']
    target_max_number = config['target_max_number']  # 靶标推荐最大值,根据文献数量进行排序

    os.makedirs('results/' + disease_name, exist_ok=True)

    return disease_name, reported_number, interaction_num, target_max_number


def research_status_test(protein_list_path: str) -> None:
    """分析蛋白质靶标研究状态并生成可视化报告。

    该函数执行完整的靶标分析流程，包括：
    1. 从输入文件加载蛋白质列表
    2. 通过PPI网络扩展靶标集合
    3. 对靶标进行分类和研究状态验证
    4. 生成交互式可视化报告

    Args:
        protein_list_path (str): 蛋白质列表文件路径，文件应为每行一个基因符号的文本文件

    Returns:
        None: 无直接返回值，但会生成以下输出文件：
            - HTML报告: Target_classification.html, PPI_Target_classification.html
            - JSON可视化数据: Target_tree.json, PPI_Target_sunburst.json 等
            - 控制台输出统计信息

    Raises:
        FileNotFoundError: 如果输入的蛋白质列表文件不存在
        ConnectionError: 如果无法连接本地Elasticsearch服务

    Notes:
        1. 需要预先配置本地Elasticsearch服务(默认端口9200)
        2. 依赖以下辅助模块:
           - get: 数据获取模块
           - output: 结果输出模块
        3. 分类标准:
           - h_dr: 已报道的疾病相关靶标
           - fa: FDA批准药物靶标
           - ct: 临床试验阶段靶标
           - ot: 其他靶标
    """
    disease_name, reported_number, interaction_num, target_max_number = set_config_auto()

    Symbol_PPI_list, Symbol_To_Target_wm, Symbol_To_Fullname, Symbol_list = (
        get.get_data(protein_list_path, interaction_num))

    p_h_dr, p_no_dr, p_fa, p_ct, p_ot = classify_targets_wm(Symbol_To_Target_wm, Symbol_PPI_list)
    h_dr, no_dr, fa, ct, ot = classify_targets_wm(Symbol_To_Target_wm, Symbol_list)

    classify_targets_html(h_dr, no_dr, fa, ct, ot)
    classify_targets_html(p_h_dr, p_no_dr, p_fa, p_ct, p_ot)

    es = Elasticsearch(
        ['http://localhost:9200/']
    )

    t0 = time.time()
    print('Program start time:', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

    # symbol对应的靶标
    fda_no_review, fda_review, ct_no_review, ct_review = report_info(fa, ct, es, disease_name, reported_number)
    p_fda_no_review, p_fda_review, p_ct_no_review, p_ct_review = report_info(p_fa, p_ct, es, disease_name,
                                                                             reported_number)

    # 生成靶标信息的Tree图
    output.all_targets_tree(fda_no_review, fda_review, ct_no_review, ct_review)
    output.all_targets_tree(p_fda_no_review, p_fda_review, p_ct_no_review, p_ct_review)

    t1 = time.time()
    print('Program end time:', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print("Running Time：%.6fs" % (t1 - t0))

    print('The number of recommend targets is: ', len(fda_no_review + ct_no_review))
    print('The number of recommend PPI targets is: ', len(p_fda_no_review + p_ct_no_review))

    # 生成靶标列表
    no_review = fda_no_review + ct_no_review
    p_no_review = p_fda_no_review + p_ct_no_review

    fda_no_review = output.new_targets_list(fda_no_review, output.sort_targets(no_review, target_max_number, es))
    ct_no_review = output.new_targets_list(ct_no_review, output.sort_targets(no_review, target_max_number, es))
    p_fda_no_review = output.new_targets_list(p_fda_no_review, output.sort_targets(p_no_review, target_max_number, es))
    p_ct_no_review = output.new_targets_list(p_ct_no_review, output.sort_targets(p_no_review, target_max_number, es))

    # 获得全部靶标药物推荐的旭日图，以及每个靶标对应的药物信息和药物热度
    output.get_sunburst_tree_bar(fda_no_review, ct_no_review, fa, disease_name, reported_number,
                                 Symbol_To_Target_wm, es)
    output.get_sunburst_tree_bar(p_fda_no_review, p_ct_no_review, p_fa, disease_name,
                                 reported_number, Symbol_To_Target_wm, es)


def update_config(disease_name, target_max_number, reported_number, interaction_num, config_path="config.json"):
    """读取用户输入并更新配置文件"""

    updates = {"disease_name": disease_name, "target_max_number": int(target_max_number),
               "interaction_num": int(interaction_num), "reported_number": int(reported_number)}

    with open(config_path, "r+", encoding="utf-8") as f:
        config = json.load(f)
        config.update(updates)  # 合并更新

        # 写回文件
        f.seek(0)
        json.dump(config, f, indent=4)
        f.truncate()

    print("配置文件已更新！")
