import analysis, report, output, get, compute
import os
from tqdm import tqdm
import pandas as pd


def from_tcm_or_formula(tcm_or_formula_id,
                        proteins_id=None,
                        score=990,
                        out_for_excel=True,
                        out_for_cytoscape=True,
                        research_status_test=True,
                        safety_research=False,
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
            :param research_status_test: 
            :param out_for_excel: 
    """


    if tcm_or_formula_id[0][2] == 'F':  # 判断输入是否为复方的
        formula = get.get_formula('DNFID', tcm_or_formula_id)  # 获取该复方的信息
        formula_tcm_links = get.get_formula_tcm_links('DNFID', formula['DNFID'])
        tcm = get.get_tcm('DNHID', formula_tcm_links['DNHID'])
        sd_formula_links = get.get_SD_Formula_links('DNFID', formula['DNFID'])
        sd = get.get_SD('DNSID', sd_formula_links['DNSID'])
    else:
        tcm = get.get_tcm('DNHID', tcm_or_formula_id)
        formula_tcm_links = get.get_formula_tcm_links('DNHID', tcm['DNHID'])
        formula = get.get_formula('DNFID', formula_tcm_links['DNFID'])
        sd_formula_links = get.get_SD_Formula_links('DNFID', formula['DNFID'])
        sd = get.get_SD('DNSID', sd_formula_links['DNSID'])

    tcm_chem_links = get.get_tcm_chem_links('DNHID', tcm['DNHID'])
    chem = get.get_chemicals('DNCID', tcm_chem_links['DNCID'])
    chem_protein_links = get.get_chem_protein_links('DNCID', chem['DNCID'], score)

    if proteins_id is None:
        proteins = get.get_proteins('Ensembl_ID', chem_protein_links['Ensembl_ID'])
    else:
        proteins = get.get_proteins('Ensembl_ID', proteins_id)

    SD_df = pd.DataFrame(sd)
    SD_Formula_Links_df = pd.DataFrame(sd_formula_links)
    formula_df = pd.DataFrame(formula)
    formula_tcm_links_df = pd.DataFrame(formula_tcm_links)
    tcm_df = pd.DataFrame(tcm)
    tcm_chem_links_df = pd.DataFrame(tcm_chem_links)
    chem_df = pd.DataFrame(chem)
    chem_protein_links_df = pd.DataFrame(chem_protein_links)
    protein_df = pd.DataFrame(proteins)

    if out_graph:
        output.vis(SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                   tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df, path)

    if out_for_cytoscape:
        output.out_for_cyto(SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                           tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df, path)

    if out_for_excel:
        with pd.ExcelWriter(f"{path}/results.xlsx") as writer:
            SD_df.to_excel(writer, sheet_name="辩证信息", index=False)
            SD_Formula_Links_df.to_excel(writer, sheet_name="辩证-复方信息", index=False)
            formula_df.to_excel(writer, sheet_name="复方信息", index=False)
            formula_tcm_links_df.to_excel(writer, sheet_name="复方-中药连接", index=False)
            tcm_df.to_excel(writer, sheet_name="中药信息", index=False)
            tcm_chem_links_df.to_excel(writer, sheet_name="中药-化合物连接", index=False)
            chem_df.to_excel(writer, sheet_name="化合物信息", index=False)
            chem_protein_links_df.to_excel(writer, sheet_name="化合物-靶点连接", index=False)
            protein_df.to_excel(writer, sheet_name="靶点信息", index=False)


    if research_status_test:
        analysis.update_config()

        protein_research_test = protein_df['gene_name']
        protein_research_test = pd.DataFrame(protein_research_test)

        df_expanded = protein_research_test.copy()
        df_expanded["split_values"] = df_expanded["gene_name"].str.split(r"[\s/]+")
        df_expanded = df_expanded.explode("split_values")

        df_expanded = pd.DataFrame(df_expanded['split_values'])
        df_expanded.rename(columns={"split_values": "gene_name"}, inplace=True)
        df_expanded.to_excel("Protein_List.xlsx", index=False)

        analysis.research_status_test(f"Protein_List.xlsx")

    if safety_research:
        return 0


    if re:
        return (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df)
    else:
        return 0


import time
from datetime import datetime

def from_proteins(proteins,
                  score=990,
                  out_graph=True,
                  out_for_cytoscape=True,
                  out_for_excel=True,
                  research_status_test=True,
                  safety_research=False,
                  random_state=None,
                  num=1000,
                  tcm_component=False,
                  formula_component=False,
                  re=True,
                  path='results/'):

    # 初始化计时和日志
    start_time = time.time()
    def log_step(step_name):
        elapsed = time.time() - start_time
        print(f"[{datetime.now().strftime('%H:%M:%S')}] [Step: {step_name}] | Time elapsed: {elapsed:.2f}s")

    # 创建输出目录
    log_step("Creating output directory")
    if not os.path.exists(path):
        os.makedirs(path)

    # 逐步执行并记录时间
    log_step("Fetching proteins")
    proteins = get.get_proteins('Ensembl_ID', proteins)

    log_step("Fetching chem-protein links")
    chem_protein_links = get.get_chem_protein_links('Ensembl_ID', proteins['Ensembl_ID'], score)
    if chem_protein_links.empty:
        raise ValueError(f"未找到化合物-蛋白质连接(score={score})，请降低score值")

    log_step("Fetching chemicals")
    chem = get.get_chemicals('DNCID', chem_protein_links['DNCID'])

    log_step("Fetching TCM-chem links")
    tcm_chem_links = get.get_tcm_chem_links('DNCID', chem['DNCID'])

    log_step("Fetching TCMs")
    tcm = get.get_tcm('DNHID', tcm_chem_links['DNHID'])

    log_step("Fetching formula-TCM links")
    formula_tcm_links = get.get_formula_tcm_links('DNHID', tcm['DNHID'])

    log_step("Fetching formulas")
    formula = get.get_formula('DNFID', formula_tcm_links['DNFID'])

    log_step("Fetching SD-formula links")
    sd_formula_links = get.get_SD_Formula_links('DNFID', formula['DNFID'])

    log_step("Fetching SDs")
    sd = get.get_SD('DNSID', sd_formula_links['DNSID'])

    log_step("Computing scores")
    tcm, chem, formula = compute.score(tcm, tcm_chem_links, chem,
                                      chem_protein_links, formula, formula_tcm_links)

    tcms, formulas = None, None
    if tcm_component:
        log_step("Computing TCM components")
        tcms = compute.component(tcm.loc[tcm['Importance Score'] != 1.0], random_state, num)
    if formula_component:
        log_step("Computing formula components")
        formulas = compute.component(formula.loc[formula['Importance Score'] != 1.0], random_state, num)

    # 转换为 DataFrame
    log_step("Converting to DataFrames")
    SD_df = pd.DataFrame(sd)
    SD_Formula_Links_df = pd.DataFrame(sd_formula_links)
    formula_df = pd.DataFrame(formula)
    formula_tcm_links_df = pd.DataFrame(formula_tcm_links)
    tcm_df = pd.DataFrame(tcm)
    tcm_chem_links_df = pd.DataFrame(tcm_chem_links)
    chem_df = pd.DataFrame(chem)
    chem_protein_links_df = pd.DataFrame(chem_protein_links)
    protein_df = pd.DataFrame(proteins)

    # 输出结果
    if out_for_cytoscape:
        log_step("Exporting for Cytoscape")
        output.out_for_cyto(SD_df, SD_Formula_Links_df, formula_df,
                           formula_tcm_links_df, tcm_df,
                           tcm_chem_links_df, chem_df,
                           chem_protein_links_df, protein_df, path)

    if out_graph:
        log_step("Generating graphs")
        output.vis(SD_df, SD_Formula_Links_df, formula_df,
                  formula_tcm_links_df, tcm_df,
                  tcm_chem_links_df, chem_df,
                  chem_protein_links_df, protein_df, path)

    if out_for_excel:
        log_step("Exporting to Excel")
        with pd.ExcelWriter(f"{path}/results.xlsx") as writer:
            SD_df.to_excel(writer, sheet_name="辩证信息", index=False)
            SD_Formula_Links_df.to_excel(writer, sheet_name="辩证-复方信息", index=False)
            formula_df.to_excel(writer, sheet_name="复方信息", index=False)
            formula_tcm_links_df.to_excel(writer, sheet_name="复方-中药连接", index=False)
            tcm_df.to_excel(writer, sheet_name="中药信息", index=False)
            tcm_chem_links_df.to_excel(writer, sheet_name="中药-化合物连接", index=False)
            chem_df.to_excel(writer, sheet_name="化合物信息", index=False)
            chem_protein_links_df.to_excel(writer, sheet_name="化合物-靶点连接", index=False)
            protein_df.to_excel(writer, sheet_name="靶点信息", index=False)

    if research_status_test:
        log_step("Running research status test")
        analysis.update_config()
        protein_research_test = protein_df['gene_name'].to_frame()
        df_expanded = protein_research_test.assign(
            split_values=protein_research_test['gene_name'].str.split(r"[\s/]+")
        ).explode("split_values")[['split_values']].rename(columns={"split_values": "gene_name"})
        df_expanded.to_excel("Protein_List.xlsx", index=False)
        analysis.research_status_test("Protein_List.xlsx")

    if safety_research:
        log_step("Safety research (placeholder)")
        return 0

    log_step("Completed")
    if re:
        return (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df,
                tcms, formulas)


def from_SD(SD_ID,
            score=990,
            out_graph=True,
            out_for_cytoscape=True,
            out_for_excel=True,
            research_status_test=True,
            safety_research=False,
            re=True,
            path='results/'
            ):
    """
    从中医辨证开始，构建中药整合药理学

    Args: (原有参数说明保持不变)
    """
    # 初始化计时和日志
    start_time = time.time()

    def log_step(step_name):
        elapsed = time.time() - start_time
        print(f"[{datetime.now().strftime('%H:%M:%S')}] [Step: {step_name}] | Time elapsed: {elapsed:.2f}s")

    log_step("Creating output directory")
    if not os.path.exists(path):
        os.makedirs(path)

    log_step("Fetching SD data")
    SD = get.get_SD('DNSID', SD_ID)

    log_step("Fetching SD-Formula links")
    SD_Formula_Links = get.get_SD_Formula_links('DNSID', SD['DNSID'])

    log_step("Fetching formulas")
    formula = get.get_formula('DNFID', SD_Formula_Links['DNFID'])

    log_step("Fetching formula-TCM links")
    formula_tcm_links = get.get_formula_tcm_links('DNFID', formula['DNFID'])

    log_step("Fetching TCMs")
    tcm = get.get_tcm('DNHID', formula_tcm_links['DNHID'])

    log_step("Fetching TCM-chem links")
    tcm_chem_links = get.get_tcm_chem_links('DNHID', tcm['DNHID'])

    log_step("Fetching chemicals")
    chem = get.get_chemicals('DNCID', tcm_chem_links['DNCID'])

    log_step("Fetching chem-protein links")
    chem_protein_links = get.get_chem_protein_links('DNCID', chem['DNCID'], score=score)

    log_step("Fetching proteins")
    protein = get.get_proteins('Ensembl_ID', chem_protein_links['Ensembl_ID'])

    log_step("Converting to DataFrames")
    # 转换为DataFrame
    SD_df = pd.DataFrame(SD)
    SD_Formula_Links_df = pd.DataFrame(SD_Formula_Links)
    formula_df = pd.DataFrame(formula)
    formula_tcm_links_df = pd.DataFrame(formula_tcm_links)
    tcm_df = pd.DataFrame(tcm)
    tcm_chem_links_df = pd.DataFrame(tcm_chem_links)
    chem_df = pd.DataFrame(chem)
    chem_protein_links_df = pd.DataFrame(chem_protein_links)
    protein_df = pd.DataFrame(protein)

    if out_graph:
        log_step("Generating visualization graphs")
        output.vis(SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                   tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df, path)

    if out_for_cytoscape:
        log_step("Exporting for Cytoscape")
        output.out_for_cyto(SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                            tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df, path)

    if out_for_excel:
        log_step("Exporting to Excel")
        with pd.ExcelWriter(f"{path}/results.xlsx") as writer:
            SD_df.to_excel(writer, sheet_name="辩证信息", index=False)
            SD_Formula_Links_df.to_excel(writer, sheet_name="辩证-复方信息", index=False)
            formula_df.to_excel(writer, sheet_name="复方信息", index=False)
            formula_tcm_links_df.to_excel(writer, sheet_name="复方-中药连接", index=False)
            tcm_df.to_excel(writer, sheet_name="中药信息", index=False)
            tcm_chem_links_df.to_excel(writer, sheet_name="中药-化合物连接", index=False)
            chem_df.to_excel(writer, sheet_name="化合物信息", index=False)
            chem_protein_links_df.to_excel(writer, sheet_name="化合物-靶点连接", index=False)
            protein_df.to_excel(writer, sheet_name="靶点信息", index=False)

    if research_status_test:
        log_step("Running research status test")
        analysis.update_config()

        protein_research_test = protein_df['gene_name']
        protein_research_test = pd.DataFrame(protein_research_test)

        df_expanded = protein_research_test.copy()
        df_expanded["split_values"] = df_expanded["gene_name"].str.split(r"[\s/]+")
        df_expanded = df_expanded.explode("split_values")

        df_expanded = pd.DataFrame(df_expanded['split_values'])
        df_expanded.rename(columns={"split_values": "gene_name"}, inplace=True)
        df_expanded.to_excel("Protein_List.xlsx", index=False)

        analysis.research_status_test(f"Protein_List.xlsx")

    if safety_research:
        log_step("Running safety research")
        # 这里添加安全性研究的代码
        return 0

    log_step("Completed all operations")
    if re:
        return (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df)
    else:
        return 0

if __name__ == '__main__':
    from_SD(['DNS0001'])
    from_proteins(['ENSP00000000233'], score=660, out_for_cytoscape=False, out_graph=False,tcm_component=False, formula_component=False)