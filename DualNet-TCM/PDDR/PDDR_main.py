import time
import os
import json
from elasticsearch import Elasticsearch
import PDDR_output, PDDR_analysis, PDDR_get


def set_config_auto():
    with open('config.json', 'r') as f:
        config = json.load(f)

    disease_name = config['disease_name']
    reported_number = config['reported_number']
    interaction_num = config['interaction_num']
    target_max_number = config['target_max_number']  # 靶标推荐最大值,根据文献数量进行排序

    os.makedirs('output/' + disease_name + ' reported_number_' + str(reported_number), exist_ok=True)
    os.makedirs('output/' + disease_name + ' reported_number_' + str(reported_number) + '/Target', exist_ok=True)
    os.makedirs('output/' + disease_name + ' reported_number_' + str(reported_number) + '/PPI_Target', exist_ok=True)

    return disease_name, reported_number, interaction_num, target_max_number



def PDDR(protein_list_path):
    """
    整体Snoopy代码运行流
    :return: 返回result文件夹中的可视化方案及文件（可选）
    """
    disease_name, reported_number, interaction_num, target_max_number = set_config_auto()

    Symbol_PPI_list, Symbol_To_Target_wm, Symbol_To_Fullname, Symbol_list = PDDR_get.get_data(protein_list_path, interaction_num)

    p_h_dr,p_no_dr,p_fa,p_ct,p_ot = PDDR_analysis.classify_targets_wm(Symbol_To_Target_wm, Symbol_PPI_list)
    h_dr,no_dr,fa,ct,ot = PDDR_analysis.classify_targets_wm(Symbol_To_Target_wm, Symbol_list)

    PDDR_analysis.classify_targets_html(h_dr,no_dr,fa,ct,ot,'Target')
    PDDR_analysis.classify_targets_html(p_h_dr,p_no_dr,p_fa,p_ct,p_ot,'PPI_Target')

    #es = Elasticsearch(timeout=30, max_retries=10, retry_on_timeout=True)

    es = Elasticsearch(
        ['http://localhost:9200/']
    )

    t0 = time.time()
    print('Program start time:', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))

    # symbol对应的靶标
    fda_no_review, fda_review, ct_no_review, ct_review = PDDR_analysis.report_info(fa, ct, es, disease_name, reported_number)
    p_fda_no_review, p_fda_review, p_ct_no_review, p_ct_review = PDDR_analysis.report_info(p_fa, p_ct, es, disease_name, reported_number)

    # 生成靶标信息的Tree图
    PDDR_output.all_targets_tree(fda_no_review, fda_review, ct_no_review, ct_review, 'Target')
    PDDR_output.all_targets_tree(p_fda_no_review, p_fda_review, p_ct_no_review, p_ct_review, 'PPI_Target')

    t1 = time.time()
    print('Program end time:', time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
    print("Running Time：%.6fs" % (t1 - t0))

    print('The number of recommend targets is: ', len(fda_no_review + ct_no_review))
    print('The number of recommend PPI targets is: ', len(p_fda_no_review + p_ct_no_review))

    # 生成靶标列表
    no_review = fda_no_review + ct_no_review
    p_no_review = p_fda_no_review + p_ct_no_review

    fda_no_review = PDDR_output.new_targets_list(fda_no_review, PDDR_output.sort_targets(no_review, target_max_number, es))
    ct_no_review = PDDR_output.new_targets_list(ct_no_review, PDDR_output.sort_targets(no_review, target_max_number, es))
    p_fda_no_review = PDDR_output.new_targets_list(p_fda_no_review, PDDR_output.sort_targets(p_no_review, target_max_number, es))
    p_ct_no_review = PDDR_output.new_targets_list(p_ct_no_review, PDDR_output.sort_targets(p_no_review, target_max_number, es))

    # 获得全部靶标药物推荐的旭日图，以及每个靶标对应的药物信息和药物热度
    PDDR_output.get_sunburst_tree_bar('Target', fda_no_review, ct_no_review, fa, disease_name, reported_number, Symbol_To_Target_wm, es)
    PDDR_output.get_sunburst_tree_bar('PPI_Target', p_fda_no_review, p_ct_no_review, p_fa, disease_name, reported_number, Symbol_To_Target_wm, es)

    print('The program is finished!')
    print('-------------------------------------------------------------------')
    input('Please enter any key to exit')


if __name__ == "__main__":
    PDDR(protein_list_path = "Example_protein_list.txt")