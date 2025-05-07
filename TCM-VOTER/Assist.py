import pandas as pd
import analysis
import report


def create_dataframes(sd, sd_formula_links, formula, formula_tcm_links,
                     tcm, tcm_chem_links, chem, chem_protein_links, proteins):
    """返回多个DataFrame而不是字典"""
    return (
        pd.DataFrame(sd),                    # SD_df
        pd.DataFrame(sd_formula_links),      # SD_Formula_Links_df
        pd.DataFrame(formula),               # formula_df
        pd.DataFrame(formula_tcm_links),     # formula_tcm_links_df
        pd.DataFrame(tcm),                   # tcm_df
        pd.DataFrame(tcm_chem_links),       # tcm_chem_links_df
        pd.DataFrame(chem),                 # chem_df
        pd.DataFrame(chem_protein_links),    # chem_protein_links_df
        pd.DataFrame(proteins)               # protein_df
    )


def save_results_to_excel(path,
                          SD_df,
                          SD_Formula_Links_df,
                          formula_df,
                          formula_tcm_links_df,
                          tcm_df,
                          tcm_chem_links_df,
                          chem_df,
                          chem_protein_links_df,
                          protein_df):
    """
    将多个DataFrame保存到Excel文件的不同工作表中

    参数:
    path (str): 保存Excel文件的目录路径
    SD_df (DataFrame): 辩证信息数据
    SD_Formula_Links_df (DataFrame): 辩证-复方信息数据
    formula_df (DataFrame): 复方信息数据
    formula_tcm_links_df (DataFrame): 复方-中药连接数据
    tcm_df (DataFrame): 中药信息数据
    tcm_chem_links_df (DataFrame): 中药-化合物连接数据
    chem_df (DataFrame): 化合物信息数据
    chem_protein_links_df (DataFrame): 化合物-靶点连接数据
    protein_df (DataFrame): 靶点信息数据
    """
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


def analyze_proteins(DiseaseName, target_max_number, report_number, interaction_number, protein_df):
    """
    分析蛋白质数据并生成研究报告

    参数:
    analysis: 分析工具对象
    DiseaseName (str): 疾病名称
    target_max_number (int): 最大靶点数量
    report_number (int): 报告数量
    interaction_number (int): 交互数量
    protein_df (DataFrame): 包含蛋白质/基因数据的DataFrame

    返回:
    无 (结果会保存到文件并通过analysis对象处理)
    """
    # 更新分析配置
    analysis.update_config(DiseaseName, target_max_number, report_number, interaction_number)

    try:
        # 提取基因名称并处理
        protein_research_test = protein_df['gene_name'].copy()
        protein_research_test = pd.DataFrame(protein_research_test)

        # 分割可能包含多个基因的字段(处理用空格或斜杠分隔的基因名)
        df_expanded = protein_research_test.copy()
        df_expanded["split_values"] = df_expanded["gene_name"].str.split(r"[\s/]+")
        df_expanded = df_expanded.explode("split_values")

        # 清理并保存基因列表
        df_expanded = pd.DataFrame(df_expanded['split_values'])
        df_expanded.rename(columns={"split_values": "gene_name"}, inplace=True)
        df_expanded.drop_duplicates(inplace=True)  # 去除重复基因名
        df_expanded.dropna(inplace=True)  # 去除空值

        # 保存处理后的基因列表
        df_expanded.to_excel("Protein_List.xlsx", index=False)

        # 执行研究状态测试
        analysis.research_status_test("Protein_List.xlsx")

        print("蛋白质分析完成，结果已保存到 Protein_List.xlsx")

    except Exception as e:
        print(f"分析过程中发生错误: {str(e)}")
        raise


def generate_toxicity_report(protein_df, chem_df, formula_df, tcm_df):
    """
    生成毒性报告

    参数:
    report: 报告生成工具对象
    protein_df (DataFrame): 包含蛋白质/基因数据的DataFrame
    chem_df (DataFrame): 包含化合物数据的DataFrame
    formula_df (DataFrame): 包含复方数据的DataFrame
    tcm_df (DataFrame): 包含中药数据的DataFrame

    返回:
    无 (结果会通过report对象处理并生成报告)
    """
    try:
        # 提取所需列并转换为DataFrame
        protein_data = pd.DataFrame(protein_df['gene_name'])
        chem_data = pd.DataFrame(chem_df['Name'])
        formula_data = pd.DataFrame(formula_df['name'])
        tcm_data = pd.DataFrame(tcm_df['cn_name'])

        # 过滤毒性数据
        toxic_data = report.filter_toxic_data(
            protein_data, chem_data, formula_data, tcm_data
        )
        toxic_formula, toxic_herb, toxic_chemical, toxic_protein = toxic_data

        # 生成毒性报告
        report.generate_toxicity_report(
            toxic_formula,
            toxic_herb,
            toxic_chemical,
            toxic_protein
        )

        print("毒性报告生成完成")

    except KeyError as e:
        print(f"数据框中缺少必要的列: {str(e)}")
        raise
    except Exception as e:
        print(f"生成毒性报告时发生错误: {str(e)}")
        raise
