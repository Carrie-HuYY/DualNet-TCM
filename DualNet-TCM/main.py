import get
import compute
import os
import output
from tqdm import tqdm
import pandas as pd
import analysis


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

    total_steps = 14
    pbar = tqdm(total=total_steps, desc="中药网络药理学分析进度")

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


    # 转换为DataFrame
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
        pbar.set_postfix_str("当前步骤: 生成可视化图表")
        output.vis(SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                   tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df, path)
        pbar.update(1)

    if out_for_cytoscape:
        pbar.set_postfix_str("当前步骤: 生成Cytoscape文件")
        output.out_for_cyto(SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                           tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df, path)
        pbar.update(1)

    if out_for_excel:
        pbar.set_postfix_str("当前步骤: 导出Excel结果")
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
        pbar.update(1)

    if research_status_test:
        pbar.set_postfix_str("当前步骤: 研究价值评定")
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
        pbar.update(1)

    if safety_research:
        pbar.set_postfix_str("当前步骤: 安全性研究")
        # 这里添加安全性研究的代码
        pbar.update(1)

    pbar.set_postfix_str("当前步骤: 完成")
    pbar.close()

    if re:
        return (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df)
    else:
        return 0


def from_proteins(proteins,
                  score=990,
                  out_graph=True,
                  out_for_cytoscape=True,
                  out_for_excel=True,
                  research_status_test=False,
                  safety_research=False,
                  random_state=None,
                  num=1000,
                  tcm_component=True,
                  formula_component=True,
                  re=True,
                  path='results/'):
    """
    进行逆向网络药理学分析（带详细进度条）
    """
    # 初始化主进度条
    main_steps = 13 + sum([out_graph, out_for_cytoscape, out_for_excel,
                           research_status_test, safety_research])
    main_pbar = tqdm(total=main_steps, desc="总体进度", position=0)

    try:
        # 创建结果目录
        if not os.path.exists(path):
            os.makedirs(path)

        # 阶段1：数据获取
        with tqdm(total=8, desc="数据获取阶段", position=1, leave=False) as data_pbar:
            main_pbar.set_postfix_str("获取蛋白质数据")
            proteins = get.get_proteins('Ensembl_ID', proteins)
            data_pbar.update(1)
            main_pbar.update(1)

            main_pbar.set_postfix_str("获取化合物-蛋白质连接")
            chem_protein_links = get.get_chem_protein_links('Ensembl_ID', proteins['Ensembl_ID'], score)
            if chem_protein_links.empty:
                raise ValueError(f"未找到化合物-蛋白质连接(score={score})，请降低score值")
            data_pbar.update(1)
            main_pbar.update(1)

            main_pbar.set_postfix_str("获取化合物数据")
            chem = get.get_chemicals('DNCID', chem_protein_links['DNCID'])
            data_pbar.update(1)
            main_pbar.update(1)

            main_pbar.set_postfix_str("获取中药-化合物连接")
            tcm_chem_links = get.get_tcm_chem_links('DNCID', chem['DNCID'])
            data_pbar.update(1)
            main_pbar.update(1)

            main_pbar.set_postfix_str("获取中药数据")
            tcm = get.get_tcm('DNHID', tcm_chem_links['DNHID'])
            data_pbar.update(1)
            main_pbar.update(1)

            main_pbar.set_postfix_str("获取复方-中药连接")
            formula_tcm_links = get.get_formula_tcm_links('DNHID', tcm['DNHID'])
            data_pbar.update(1)
            main_pbar.update(1)

            main_pbar.set_postfix_str("获取复方数据")
            formula = get.get_formula('DNFID', formula_tcm_links['DNFID'])
            data_pbar.update(1)
            main_pbar.update(1)

            main_pbar.set_postfix_str("获取辩证-复方连接")
            sd_formula_links = get.get_SD_Formula_links('DNFID', formula['DNFID'])
            data_pbar.update(1)
            main_pbar.update(1)

            main_pbar.set_postfix_str("获取辩证数据")
            sd = get.get_SD('DNSID', sd_formula_links['DNSID'])
            data_pbar.update(1)
            main_pbar.update(1)

        # 阶段2：数据处理
        with tqdm(total=2, desc="数据处理阶段", position=1, leave=False) as process_pbar:
            main_pbar.set_postfix_str("数据过滤")
            formula, formula_tcm_links, tcm, tcm_chem_links, chem, chem_protein_links, proteins = \
                analysis.dfs_filter(formula, formula_tcm_links, tcm, tcm_chem_links,
                                    chem, chem_protein_links, proteins)
            process_pbar.update(1)
            main_pbar.update(1)

            main_pbar.set_postfix_str("计算评分")
            tcm, chem, formula = compute.score(tcm, tcm_chem_links, chem,
                                               chem_protein_links, formula, formula_tcm_links)
            process_pbar.update(1)
            main_pbar.update(1)

        # 阶段3：优化模型
        with tqdm(total=2, desc="模型优化阶段", position=1, leave=False) as model_pbar:
            tcms, formulas = None, None

            if tcm_component:
                main_pbar.set_postfix_str("中药组成优化")
                tcms = compute.component(tcm.loc[tcm['Importance Score'] != 1.0], random_state, num)
                model_pbar.update(1)
                main_pbar.update(1)

            if formula_component:
                main_pbar.set_postfix_str("方剂组合优化")
                formulas = compute.component(formula.loc[formula['Importance Score'] != 1.0], random_state, num)
                model_pbar.update(1)
                main_pbar.update(1)

        # 阶段4：结果输出
        with tqdm(total=5, desc="结果输出阶段", position=1, leave=False) as output_pbar:
            # 转换为DataFrame
            main_pbar.set_postfix_str("数据转换")
            SD_df = pd.DataFrame(sd)
            SD_Formula_Links_df = pd.DataFrame(sd_formula_links)
            formula_df = pd.DataFrame(formula)
            formula_tcm_links_df = pd.DataFrame(formula_tcm_links)
            tcm_df = pd.DataFrame(tcm)
            tcm_chem_links_df = pd.DataFrame(tcm_chem_links)
            chem_df = pd.DataFrame(chem)
            chem_protein_links_df = pd.DataFrame(chem_protein_links)
            protein_df = pd.DataFrame(proteins)
            output_pbar.update(1)
            main_pbar.update(1)

            if out_for_cytoscape:
                main_pbar.set_postfix_str("生成Cytoscape文件")
                output.out_for_cyto(SD_df, SD_Formula_Links_df, formula_df,
                                    formula_tcm_links_df, tcm_df,
                                    tcm_chem_links_df, chem_df,
                                    chem_protein_links_df, protein_df, path)
                output_pbar.update(1)
                main_pbar.update(1)

            if out_graph:
                main_pbar.set_postfix_str("生成可视化图表")
                output.vis(SD_df, SD_Formula_Links_df, formula_df,
                           formula_tcm_links_df, tcm_df,
                           tcm_chem_links_df, chem_df,
                           chem_protein_links_df, protein_df, path)
                output_pbar.update(1)
                main_pbar.update(1)

            if out_for_excel:
                main_pbar.set_postfix_str("导出Excel结果")
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
                output_pbar.update(1)
                main_pbar.update(1)

            if research_status_test:
                main_pbar.set_postfix_str("研究价值评定")
                analysis.update_config()
                protein_research_test = protein_df['gene_name'].to_frame()
                df_expanded = protein_research_test.assign(
                    split_values=protein_research_test['gene_name'].str.split(r"[\s/]+")
                ).explode("split_values")[['split_values']].rename(columns={"split_values": "gene_name"})
                df_expanded.to_excel("Protein_List.xlsx", index=False)
                analysis.research_status_test("Protein_List.xlsx")
                output_pbar.update(1)
                main_pbar.update(1)

            if safety_research:
                main_pbar.set_postfix_str("安全性研究")
                # 添加安全性研究代码
                output_pbar.update(1)
                main_pbar.update(1)

        main_pbar.set_postfix_str("完成")
        if re:
            return (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                    tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df)

    except Exception as e:
        main_pbar.set_postfix_str(f"错误: {str(e)}")
        raise
    finally:
        main_pbar.close()


def from_SD(SD_ID,
            score=990,
            out_graph=True,
            out_for_cytoscape=True,
            out_for_excel=True,
            research_status_test=False,
            safety_research=False,
            re=True,
            path='results/'
            ):
    """
    从中医辨证开始，构建中药整合药理学

    Args: (原有参数说明保持不变)
    """
    total_steps = 10
    if out_graph:
        total_steps += 1
    if out_for_cytoscape:
        total_steps += 1
    if out_for_excel:
        total_steps += 1
    if research_status_test:
        total_steps += 1
    if safety_research:
        total_steps += 1

    # 初始化进度条
    pbar = tqdm(total=total_steps, desc="中药网络药理学分析进度")

    # 创建目录
    if not os.path.exists(path):
        os.makedirs(path)
    pbar.update(1)
    pbar.set_postfix_str("当前步骤: 获取辩证信息")
    SD = get.get_SD('DNSID', SD_ID)
    pbar.update(1)

    pbar.set_postfix_str("当前步骤: 获取辩证-复方连接")
    SD_Formula_Links = get.get_SD_Formula_links('DNSID', SD['DNSID'])
    pbar.update(1)

    pbar.set_postfix_str("当前步骤: 获取复方信息")
    formula = get.get_formula('DNFID', SD_Formula_Links['DNFID'])
    pbar.update(1)

    pbar.set_postfix_str("当前步骤: 获取复方-中药连接")
    formula_tcm_links = get.get_formula_tcm_links('DNFID', formula['DNFID'])
    pbar.update(1)

    pbar.set_postfix_str("当前步骤: 获取中药信息")
    tcm = get.get_tcm('DNHID', formula_tcm_links['DNHID'])
    pbar.update(1)

    pbar.set_postfix_str("当前步骤: 获取中药-化合物连接")
    tcm_chem_links = get.get_tcm_chem_links('DNHID', tcm['DNHID'])
    pbar.update(1)

    pbar.set_postfix_str("当前步骤: 获取化合物信息")
    chem = get.get_chemicals('DNCID', tcm_chem_links['DNCID'])
    pbar.update(1)

    pbar.set_postfix_str("当前步骤: 获取化合物-靶点连接")
    chem_protein_links = get.get_chem_protein_links('DNCID', chem['DNCID'], score=score)
    pbar.update(1)

    pbar.set_postfix_str("当前步骤: 获取靶点信息")
    protein = get.get_proteins('Ensembl_ID', chem_protein_links['Ensembl_ID'])
    pbar.update(1)

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
        pbar.set_postfix_str("当前步骤: 生成可视化图表")
        output.vis(SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                   tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df, path)
        pbar.update(1)

    if out_for_cytoscape:
        pbar.set_postfix_str("当前步骤: 生成Cytoscape文件")
        output.out_for_cyto(SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                           tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df, path)
        pbar.update(1)

    if out_for_excel:
        pbar.set_postfix_str("当前步骤: 导出Excel结果")
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
        pbar.update(1)

    if research_status_test:
        pbar.set_postfix_str("当前步骤: 研究价值评定")
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
        pbar.update(1)

    if safety_research:
        pbar.set_postfix_str("当前步骤: 安全性研究")
        # 这里添加安全性研究的代码
        pbar.update(1)

    pbar.set_postfix_str("当前步骤: 完成")
    pbar.close()

    if re:
        return (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df)
    else:
        return 0


if __name__ == '__main__':
    from_tcm_or_formula(['DNH0011'])