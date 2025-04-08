import pandas as pd
from datetime import datetime


def read_toxicity_data():
    """
    从Excel读取毒性数据

    参数:

    返回:
        四个DataFrame: (targets_df, chem_df, formula_df, herb_df)
    """

    try:
        # 读取各表格数据
        targets_df = pd.read_excel("Data/Toxicity/靶点.xlsx")
        chem_df = pd.read_excel("Data/Toxicity/成分.xlsx")
        formula_df = pd.read_excel("Data/Toxicity/方剂.xlsx")
        herb_df = pd.read_excel("Data/Toxicity/中药.xlsx")

        print("数据读取成功！")
        return targets_df, chem_df, formula_df, herb_df

    except Exception as e:
        print(f"读取数据时出错: {str(e)}")
        return None, None, None, None


def generate_toxicity_report(toxic_formula, toxic_herb, toxic_chemical, toxic_protein,
                             output_file="results/toxicity_report.txt"):
    """
    生成毒性报告

    参数:
        toxic_targets_df: 有毒靶点DataFrame
        toxic_chem_df: 有毒成分DataFrame
        toxic_formula_df: 有毒方剂DataFrame
        toxic_herb_df: 有毒中药DataFrame
        output_file: 输出报告文件名

    返回:
        生成文本报告文件
    """
    # 获取当前日期
    current_date = datetime.now().strftime("%Y-%m-%d")

    # 初始化报告内容
    report_content = f"中药毒性分析报告\n生成日期: {current_date}\n\n"
    report_content += "=" * 50 + "\n"

    # 1. 有毒方剂汇总
    report_content += "1. 有毒方剂汇总:\n"
    if not toxic_formula.empty:
        for _, row in toxic_formula.iterrows():
            report_content += f"- 方剂名称: {row['formula_name']} ({row['formula_name_pinyin']})\n"
            report_content += f"  剂型: {row['dosage_form']}\n"
            report_content += f"  毒性效应: {row['toxicity_effect']}\n"
            report_content += f"  英文毒性效应: {row['toxicity_effect_en']}\n\n"
    else:
        report_content += "未发现有毒方剂记录\n\n"

    # 2. 有毒中药汇总
    report_content += "2. 有毒中药汇总:\n"
    if not toxic_herb.empty:
        for _, row in toxic_herb.iterrows():
            report_content += f"- 中药名称: {row['herb_name']} ({row['herb_name_pinyin']})\n"
            report_content += f"  拉丁名: {row['herb_name_latin']}\n"
            report_content += f"  毒性程度: {row['toxicity_degree']} ({row['toxicity_degree_en']})\n"
            report_content += f"  功效: {row['action']} ({row['action_en']})\n"
            report_content += f"  毒性效应: {row['toxic_effect']}\n"
            report_content += f"  英文毒性效应: {row['toxic_effect_en']}\n\n"
    else:
        report_content += "未发现有毒中药记录\n\n"

    # 3. 有毒成分汇总
    report_content += "3. 有毒成分汇总:\n"
    if not toxic_chemical.empty:
        for _, row in toxic_chemical.iterrows():
            report_content += f"- 成分名称: {row['component_name']} ({row['component_name_en']})\n"
            report_content += f"  分类: {row['ingredient_classification']} ({row['ingredient_classification_en']})\n"
            report_content += f"  分子量: {row['molecular_weight']}\n"
            report_content += f"  分子式: {row['molecular_formula']}\n"
            report_content += f"  CAS号: {row['CAS']}\n\n"
    else:
        report_content += "未发现有毒成分记录\n\n"

    # 4. 有毒靶点汇总
    report_content += "4. 有毒靶点汇总:\n"
    if not toxic_protein.empty:
        for _, row in toxic_protein.iterrows():
            report_content += f"- 基因符号: {row['gene_symbol']}\n"
            report_content += f"  基因全名: {row['gene_full_name']}\n"
            report_content += f"  Uniprot ID: {row['UniprotID']}\n"
            report_content += f"  毒性类型: {row['TCMSTD_Target']} ({row['TCMSTD_Target_cn']})\n\n"
    else:
        report_content += "未发现有毒靶点记录\n\n"

    # 5. 毒性类型统计
    report_content += "5. 毒性类型统计:\n"

    # 从方剂中提取毒性类型
    formula_toxicity = []
    if not toxic_formula.empty:
        for effect in toxic_formula['toxicity_effect']:
            parts = effect.split("||")
            for part in parts:
                if ":" in part:
                    tox_type = part.split(":")[1].strip()
                    formula_toxicity.append(tox_type.split("[")[0].strip())

    # 从中药中提取毒性类型
    herb_toxicity = []
    if not toxic_herb.empty:
        for effect in toxic_herb['toxic_effect']:
            parts = effect.split("||")
            for part in parts:
                if ":" in part:
                    tox_type = part.split(":")[1].strip()
                    herb_toxicity.append(tox_type.split("[")[0].strip())

    # 从靶点中提取毒性类型
    target_toxicity = []
    if not toxic_protein.empty:
        target_toxicity = toxic_protein['TCMSTD_Target_cn'].tolist()

    # 合并所有毒性类型
    all_toxicity = formula_toxicity + herb_toxicity + target_toxicity

    # 统计毒性类型
    if all_toxicity:
        toxicity_counts = pd.Series(all_toxicity).value_counts()
        for tox_type, count in toxicity_counts.items():
            report_content += f"- {tox_type}: {count}次提及\n"
    else:
        report_content += "未发现毒性类型记录\n"

    # 写入报告文件
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(report_content)

    print(f"毒性报告已生成: {output_file}")


def filter_toxic_data(targets_df, chem_df, formula_df, herb_df):
    """
    筛选有毒数据

    参数:
        targets_df: 靶点DataFrame
        chem_df: 成分DataFrame
        formula_df: 方剂DataFrame
        herb_df: 中药DataFrame

    返回:
        四个DataFrame: (toxic_targets, toxic_chem, toxic_formula, toxic_herb)
    """
    targets_toxicity, chem_toxicity, formula_toxicity, herb_toxicity = read_toxicity_data()

    toxic_formula_mask = formula_toxicity['formula_name'].isin(formula_df['name'])
    toxic_formula_df = formula_toxicity[toxic_formula_mask]

    toxic_herb_mask = herb_toxicity['herb_name'].isin(herb_df['cn_name'])
    toxic_herb_df = herb_toxicity[toxic_herb_mask]

    toxic_protein_mask = targets_toxicity['gene_symbol'].isin(targets_df['gene_name'])
    toxic_protein_df = targets_toxicity[toxic_protein_mask]

    toxic_chemical_mask = chem_toxicity['component_name'].isin(chem_df['Name'])
    toxic_chemical_df = chem_toxicity[toxic_chemical_mask]

    return toxic_formula_df, toxic_herb_df, toxic_chemical_df, toxic_protein_df


# 示例使用
if __name__ == "__main__":
    import main

    SD_df, SD_Formula_Links_df, formula_df, formula_tcm_links_df, tcm_df, tcm_chem_links_df, chem_df, chem_protein_links_df, protein_df = (
        main.from_tcm_or_formula(['DNH0158'], research_status_test=False, safety_research=True))

    protein_df = pd.DataFrame(protein_df['gene_name'])
    chem_df = pd.DataFrame(chem_df['Name'])
    formula_df = pd.DataFrame(formula_df['name'])
    tcm_df = pd.DataFrame(tcm_df['cn_name'])

    toxic_formula, toxic_herb, toxic_chemical, toxic_protein = filter_toxic_data(
        protein_df, chem_df, formula_df, tcm_df
    )

    generate_toxicity_report(toxic_formula, toxic_herb, toxic_chemical, toxic_protein)