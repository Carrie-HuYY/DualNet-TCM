import output
import get
import Assist
import compute
import os
import time
from datetime import datetime
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('network_pharmacology_analysis.log')
    ]
)
logger = logging.getLogger(__name__)


def from_chemical(chemical_id,
                  score=990,
                  DiseaseName="cough",
                  target_max_number=70,
                  report_number=0,
                  interaction_number=0,
                  out_for_excel=True,
                  out_for_cytoscape=True,
                  research_status_test=True,
                  safety_research=True,
                  out_graph=True,
                  re=True,
                  path='results'):
    """
    执行基于化学物的正向网络药理学分析

    该函数实现完整的网络药理学分析流程，包括：
    1. 化学物、蛋白质、中药成分的数据获取
    2. 网络构建与可视化
    3. 靶点蛋白研究状态分析
    4. 成分安全性评估
    5. 多种格式结果输出(Excel、Cytoscape、可视化图形)

    参数:
        chemical_id (str): 待分析化学物的ID(DNCID格式)
        score (int, 可选): 蛋白质相互作用的综合分数阈值，默认为990
        DiseaseName (str, 可选): 用于研究状态分析的疾病名称，默认为"cough"(咳嗽)
        target_max_number (int, 可选): 最大靶点数量，默认为70
        report_number (int, 可选): 生成报告数量，默认为0
        interaction_number (int, 可选): 考虑的相互作用数量，默认为0
        out_for_excel (bool, 可选): 是否输出Excel文件，默认为True
        out_for_cytoscape (bool, 可选): 是否准备Cytoscape文件，默认为True
        research_status_test (bool, 可选): 是否执行研究状态测试，默认为True
        safety_research (bool, 可选): 是否进行安全性研究，默认为True
        out_graph (bool, 可选): 是否生成可视化图形，默认为True
        re (bool, 可选): 是否返回原始结果，默认为True
        path (str, 可选): 输出目录路径，默认为'results'

    返回:
        tuple 或 int: 如果re=True返回包含分析结果的多个DataFrame(按顺序):
            - SD_df: 辨证信息数据
            - SD_Formula_Links_df: 辨证-复方关联数据
            - formula_df: 复方信息数据
            - formula_tcm_links_df: 复方-中药关联数据
            - tcm_df: 中药信息数据
            - tcm_chem_links_df: 中药-化学物关联数据
            - chem_df: 化学物信息数据
            - chem_protein_links_df: 化学物-靶点关联数据
            - protein_df: 靶点信息数据
        如果re=False则返回0

    异常:
        Exception: 执行过程中发生错误时会记录日志并重新抛出异常
    """

    # 记录开始时间
    start_time = time.time()
    logger.info("Starting network pharmacology analysis...")

    try:
        # 数据获取阶段
        logger.info("Fetching chemical data...")
        chem = get.get_chemicals('DNCID', chemical_id)

        logger.info("Fetching protein links data...")
        chem_protein_links = get.get_chem_protein_links('DNCID', chemical_id, score=score)

        logger.info("Fetching protein data...")
        proteins = get.get_proteins('Ensembl_ID', chem_protein_links['Ensembl_ID'])

        logger.info("Fetching TCM-chemical links...")
        tcm_chem_links = get.get_tcm_chem_links('DNCID', chem['DNCID'])

        logger.info("Fetching TCM data...")
        tcm = get.get_tcm('DNHID', tcm_chem_links['DNHID'])

        logger.info("Fetching formula-TCM links...")
        formula_tcm_links = get.get_formula_tcm_links('DNHID', tcm['DNHID'])

        logger.info("Fetching formula data...")
        formula = get.get_formula('DNFID', formula_tcm_links['DNFID'])

        logger.info("Fetching SD-formula links...")
        sd_formula_links = get.get_SD_Formula_links('DNFID', formula['DNFID'])

        logger.info("Fetching SD data...")
        sd = get.get_SD('DNSID', sd_formula_links['DNSID'])

        # 转换为DataFrame
        logger.info("Converting data to DataFrames...")

        (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
         tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df) = Assist.create_dataframes(
            sd, sd_formula_links, formula, formula_tcm_links,
            tcm, tcm_chem_links, chem, chem_protein_links, proteins
        )

        # 可视化输出
        if out_graph:
            logger.info("Generating network visualizations...")
            output.vis(SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                       tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df, path)

        # Cytoscape输出
        if out_for_cytoscape:
            logger.info("Preparing files for Cytoscape...")
            output.out_for_cyto(SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                                tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df, path)

        # Excel输出
        if out_for_excel:
            logger.info("Writing results to Excel file...")
            # 假设您已经有了所有需要的DataFrame
            Assist.save_results_to_excel(
                path="results/results.xlsx",
                SD_df=SD_df,
                SD_Formula_Links_df=SD_Formula_Links_df,
                formula_df=formula_df,
                formula_tcm_links_df=formula_tcm_links_df,
                tcm_df=tcm_df,
                tcm_chem_links_df=tcm_chem_links_df,
                chem_df=chem_df,
                chem_protein_links_df=chem_protein_links_df,
                protein_df=protein_df
            )

        # 研究状态测试
        if research_status_test:
            logger.info("Performing research status test...")

            Assist.analyze_proteins(DiseaseName, target_max_number, report_number, interaction_number, protein_df)

        # 安全性研究
        if safety_research:
            logger.info("Performing safety research analysis...")
            Assist.generate_toxicity_report(protein_df, chem_df, formula_df, tcm_df)

        # 计算并记录总时间
        elapsed_time = time.time() - start_time
        logger.info(f"Analysis completed successfully in {elapsed_time:.2f} seconds")

        if re:
            return (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df, tcm_df,
                    tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df)
        else:
            return 0

    except Exception as e:
        logger.error(f"Error occurred during analysis: {str(e)}", exc_info=True)
        raise


def from_tcm_or_formula(tcm_or_formula_id,
                        proteins_id=None,
                        score=990,
                        DiseaseName="cough",
                        target_max_number=70,
                        report_number=0,
                        interaction_number=0,
                        out_for_excel=True,
                        out_for_cytoscape=True,
                        research_status_test=True,
                        safety_research=True,
                        out_graph=True,
                        re=True,
                        path='results'):
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    logger_TCM_VOTER = logging.getLogger(__name__)

    logger_TCM_VOTER.info("Starting network pharmacology analysis")

    # Step 1: Get formula/TCM information
    logger_TCM_VOTER.info("Fetching formula/TCM information")
    if tcm_or_formula_id[0][2] == 'F':  # Check if input is formula
        formula = get.get_formula('DNFID', tcm_or_formula_id)
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

    # Step 2: Get chemical information
    logger_TCM_VOTER.info("Fetching TCM-chemical links")
    tcm_chem_links = get.get_tcm_chem_links('DNHID', tcm['DNHID'])
    logger_TCM_VOTER.info("Fetching chemical information")
    chem = get.get_chemicals('DNCID', tcm_chem_links['DNCID'])

    # Step 3: Get protein information
    logger_TCM_VOTER.info(f"Fetching chem-protein links with score >= {score}")
    chem_protein_links = get.get_chem_protein_links('DNCID', chem['DNCID'], score)

    if proteins_id is None:
        logger_TCM_VOTER.info("Fetching protein information")
        proteins = get.get_proteins('Ensembl_ID', chem_protein_links['Ensembl_ID'])
    else:
        logger_TCM_VOTER.info("Fetching specified protein information")
        proteins = get.get_proteins('Ensembl_ID', proteins_id)

    # Convert to DataFrames
    logger_TCM_VOTER.info("Converting results to DataFrames")

    (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
     tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df) = Assist.create_dataframes(
        sd, sd_formula_links, formula, formula_tcm_links,
        tcm, tcm_chem_links, chem, chem_protein_links, proteins
    )

    # Output visualization
    if out_graph:
        logger_TCM_VOTER.info("Generating visualization")

        output.vis(SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                   tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df, path)

    # Output for Cytoscape
    if out_for_cytoscape:
        logger_TCM_VOTER.info("Generating Cytoscape output")
        output.out_for_cyto(SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                            tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df, path)

    # Output to Excel
    if out_for_excel:
        logger_TCM_VOTER.info("Writing results to Excel")
        # 假设您已经有了所有需要的DataFrame
        Assist.save_results_to_excel(
            path="results/",
            SD_df=SD_df,
            SD_Formula_Links_df=SD_Formula_Links_df,
            formula_df=formula_df,
            formula_tcm_links_df=formula_tcm_links_df,
            tcm_df=tcm_df,
            tcm_chem_links_df=tcm_chem_links_df,
            chem_df=chem_df,
            chem_protein_links_df=chem_protein_links_df,
            protein_df=protein_df
        )

    # Research status test
    if research_status_test:
        logger_TCM_VOTER.info("Performing research status test")

        Assist.analyze_proteins(DiseaseName, target_max_number, report_number, interaction_number, protein_df)

    # Safety research (placeholder)
    if safety_research:
        logger_TCM_VOTER.info("Performing safety research test")

        Assist.generate_toxicity_report(protein_df, chem_df, formula_df, tcm_df)

    logger_TCM_VOTER.info("Analysis completed")

    if re:
        return (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df, tcm_df,
                tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df)
    else:
        return 0


def from_proteins(proteins,
                  score=990,
                  DiseaseName="cough",
                  target_max_number=70,
                  report_number=0,
                  interaction_number=0,
                  out_graph=True,
                  out_for_cytoscape=True,
                  out_for_excel=True,
                  research_status_test=True,
                  safety_research=True,
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
    (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df, tcm_df, tcm_chem_links_df,
     chem_df, chem_protein_links_df, protein_df) = (
        Assist.create_dataframes(sd, sd_formula_links, formula, formula_tcm_links, tcm,
                                 tcm_chem_links, chem, chem_protein_links, proteins))

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
        Assist.save_results_to_excel(path="results/results.xlsx",
                                     SD_df=SD_df,
                                     SD_Formula_Links_df=SD_Formula_Links_df,
                                     formula_df=formula_df,
                                     formula_tcm_links_df=formula_tcm_links_df,
                                     tcm_df=tcm_df,
                                     tcm_chem_links_df=tcm_chem_links_df,
                                     chem_df=chem_df,
                                     chem_protein_links_df=chem_protein_links_df,
                                     protein_df=protein_df)

    if research_status_test:
        log_step("Running research status test")

        Assist.analyze_proteins(DiseaseName, target_max_number, report_number, interaction_number, protein_df)

    if safety_research:
        log_step("Safety research (placeholder)")

        Assist.generate_toxicity_report(protein_df, chem_df, formula_df, tcm_df)

    log_step("Completed")
    if re:
        return (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df,
                tcms, formulas)


def from_SD(SD_ID,
            score=990,
            DiseaseName="cough",
            target_max_number=70,
            report_number=0,
            interaction_number=0,
            out_graph=True,
            out_for_cytoscape=True,
            out_for_excel=True,
            research_status_test=True,
            safety_research=True,
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

    (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df, tcm_df, tcm_chem_links_df,
     chem_df, chem_protein_links_df, protein_df) = (
        Assist.create_dataframes(SD, SD_Formula_Links, formula, formula_tcm_links, tcm,
                                 tcm_chem_links, chem, chem_protein_links, protein))

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
        Assist.save_results_to_excel("results/results.xlsx",
                                     SD_df,
                                     SD_Formula_Links_df,
                                     formula_df,
                                     formula_tcm_links_df,
                                     tcm_df,
                                     tcm_chem_links_df,
                                     chem_df,
                                     chem_protein_links_df,
                                     protein_df)

    if research_status_test:
        log_step("Running research status test")

        Assist.analyze_proteins(DiseaseName, target_max_number, report_number, interaction_number, protein_df)

    if safety_research:
        log_step("Running safety research")

        Assist.generate_toxicity_report(protein_df, chem_df, formula_df, tcm_df)

    log_step("Completed all operations")
    if re:
        return (SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df,
                tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df)
    else:
        return 0


def TCM_VOTER(SearchType,
              SearchName,
              DiseaseName="cough",
              target_max_number=70,
              report_number=0,
              interaction_number=0,
              score=990,
              out_graph=True,
              out_for_cytoscape=True,
              out_for_excel=True,
              research_status_test=True,
              safety_research=True,
              re=True,
              path='results/'
              ):
    """
    中药网络药理学分析入口函数（TCM_VOTER）

    该函数是网络药理学分析的统一入口，根据不同的搜索类型调用相应的分析函数

    参数:
        SearchType (int): 搜索类型，可选值:
            0 - 辨证(证候)
            1 - 方剂
            2 - 中药
            3 - 化学成分
            4 - 靶点蛋白
        SearchName (str): 搜索名称，根据SearchType不同而不同:
            证候名/方剂名/中药名/成分名/基因名
        DiseaseName (str, 可选): 目标疾病名称，用于研究状态分析，默认为"cough"
        target_max_number (int, 可选): 最大靶点数量限制，默认为70
        report_number (int, 可选): 报告生成数量，默认为0
        interaction_number (int, 可选): 相互作用数量限制，默认为0
        score (int, 可选): 蛋白质相互作用分数阈值，默认为990
        out_graph (bool, 可选): 是否输出可视化图形，默认为True
        out_for_cytoscape (bool, 可选): 是否生成Cytoscape格式文件，默认为True
        out_for_excel (bool, 可选): 是否输出Excel结果文件，默认为True
        research_status_test (bool, 可选): 是否执行研究状态测试，默认为True
        safety_research (bool, 可选): 是否进行安全性研究，默认为True
        re (bool, 可选): 是否返回结果数据，默认为True
        path (str, 可选): 结果输出路径，默认为'results/'

    返回:
        int: 固定返回0

    示例:
        >>> # 辨证分析
        >>> TCM_VOTER(SearchType=0, SearchName="气虚证")
        >>>
        >>> # 方剂分析
        >>> TCM_VOTER(SearchType=1, SearchName="四物汤")
        >>>
        >>> # 中药分析
        >>> TCM_VOTER(
        ...     SearchType=2,
        ...     SearchName="当归",
        ...     DiseaseName="anemia",
        ...     out_graph=False
        ... )
        >>>
        >>> # 化学成分分析
        >>> TCM_VOTER(SearchType=3, SearchName="槲皮素")
        >>>
        >>> # 靶点蛋白分析
        >>> TCM_VOTER(
        ...     SearchType=4,
        ...     SearchName="TNF",
        ...     path="custom_results/"
        ... )

    工作流程:
        1. 根据SearchType确定搜索类型
        2. 调用get模块获取对应的ID
        3. 根据类型调用相应的分析函数:
            - 辨证: from_SD()
            - 方剂/中药: from_tcm_or_formula()
            - 化学成分: from_chemical()
            - 靶点蛋白: from_proteins()
        4. 所有分析函数共享相同的参数配置
        5. 返回固定值0
    """

    if SearchType == 0:
        SearchID = get.get_SD('证候', SearchName)['DNSID']

        from_SD(SearchID,
                score=score,
                DiseaseName=DiseaseName,
                target_max_number=target_max_number,
                report_number=report_number,
                interaction_number=interaction_number,
                out_graph=out_graph,
                out_for_cytoscape=out_for_cytoscape,
                out_for_excel=out_for_excel,
                research_status_test=research_status_test,
                safety_research=safety_research,
                re=re,
                path=path
                )

    if SearchType == 1:
        SearchID = get.get_formula('name', SearchName)['DNFID']

        from_tcm_or_formula(SearchID, score=score,
                            DiseaseName=DiseaseName,
                            target_max_number=target_max_number,
                            report_number=report_number,
                            interaction_number=interaction_number,
                            out_graph=out_graph,
                            out_for_cytoscape=out_for_cytoscape,
                            out_for_excel=out_for_excel,
                            research_status_test=research_status_test,
                            safety_research=safety_research,
                            re=re,
                            path=path
                            )

    if SearchType == 2:
        SearchID = get.get_tcm('cn_name', SearchName)['DNHID']

        from_tcm_or_formula(SearchID, score=score,
                            DiseaseName=DiseaseName,
                            target_max_number=target_max_number,
                            report_number=report_number,
                            interaction_number=interaction_number,
                            out_graph=out_graph,
                            out_for_cytoscape=out_for_cytoscape,
                            out_for_excel=out_for_excel,
                            research_status_test=research_status_test,
                            safety_research=safety_research,
                            re=re,
                            path=path
                            )

    if SearchType == 3:
        SearchID = get.get_chemicals('Name', SearchName)['DNCID']

        from_chemical(SearchID, score=score,
                      DiseaseName=DiseaseName,
                      target_max_number=target_max_number,
                      report_number=report_number,
                      interaction_number=interaction_number,
                      out_graph=out_graph,
                      out_for_cytoscape=out_for_cytoscape,
                      out_for_excel=out_for_excel,
                      research_status_test=research_status_test,
                      safety_research=safety_research,
                      re=re,
                      path=path
                      )

    if SearchType == 4:
        SearchID = get.get_proteins('gene_name', SearchName)['Ensembl_ID']

        from_proteins(SearchID, score=score,
                      DiseaseName=DiseaseName,
                      target_max_number=target_max_number,
                      report_number=report_number,
                      interaction_number=interaction_number,
                      out_graph=out_graph,
                      out_for_cytoscape=out_for_cytoscape,
                      out_for_excel=out_for_excel,
                      research_status_test=research_status_test,
                      safety_research=safety_research,
                      re=re,
                      path=path
                      )

    return 0


if __name__ == '__main__':
    TCM_VOTER(SearchType=1,
              SearchName=['定喘汤'],
              score=900,
              DiseaseName="cough",
              target_max_number=70,
              report_number=0,
              interaction_number=0,
              out_graph=True,
              out_for_cytoscape=True,
              out_for_excel=True,
              research_status_test=True,
              safety_research=True,
              re=True,
              path='results/'
              )

# from_tcm_or_formula(['DNH0158'])
# from_chemical(['DNC0007'], score=0)
# from_SD(['DNS0001'])
# from_proteins(['ENSP00000000233'], score=660, out_for_cytoscape=True, out_graph=True,
# tcm_component=False, formula_component=False)
