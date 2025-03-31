import pandas as pd
import re


def parse_toxicity_data(toxic_effect):
    """解析毒性数据字符串，返回结构化字典"""
    toxicity_data = {
        'human': [],
        'rat': [],
        'mice': [],
        'rabbit': []
    }

    species_entries = re.split(r';\s*(?=[A-Za-z]+:)', toxic_effect)

    for entry in species_entries:
        if not entry:
            continue

        match = re.match(r'([A-Za-z]+):\s*(.+)', entry)
        if not match:
            continue

        species = match.group(1).lower()
        effects = match.group(2)

        effect_list = re.split(r'\|\|', effects)

        for effect in effect_list:
            effect_match = re.match(r'([^\[]+)(\[.+\])', effect.strip())
            if effect_match:
                effect_type = effect_match.group(1).strip()
                refs = effect_match.group(2)
                ref_count = len(re.findall(r'\[\d+-\d+\]', refs))

                toxicity_data[species].append({
                    'type': effect_type,
                    'references': refs,
                    'count': ref_count
                })

    return toxicity_data


def generate_toxicity_report(component_name, toxic_effect):
    """生成单个成分的毒性报告"""
    toxicity_data = parse_toxicity_data(toxic_effect)

    report_lines = [
        f"成分安全评估报告：{component_name}",
        "=" * 50
    ]

    species_order = ['human', 'rat', 'mice', 'rabbit']

    for species in species_order:
        if not toxicity_data[species]:
            continue

        species_names = {
            'human': '人类',
            'rat': '大鼠',
            'mice': '小鼠',
            'rabbit': '兔'
        }

        report_lines.append(f"\n▶ {species_names[species]}实验数据：")

        sorted_effects = sorted(toxicity_data[species],
                                key=lambda x: x['count'],
                                reverse=True)

        for effect in sorted_effects:
            type_translation = {
                'Cardiotoxicity': '心脏毒性',
                'Hepatotoxicity': '肝毒性',
                'Nephrotoxicity': '肾毒性',
                'Neurotoxicity': '神经毒性',
                'Reproductive Toxicity': '生殖毒性',
                'Developmental Toxicity': '发育毒性',
                'Genotoxicity': '遗传毒性',
                'Respiratory Toxicity': '呼吸系统毒性',
                'Digestive System Toxicity': '消化系统毒性',
                'Other Side-effects': '其他毒副作用'
            }

            chinese_type = type_translation.get(effect['type'], effect['type'])
            report_lines.append(
                f"• 观察到{chinese_type}（文献支持：{effect['references']}，共{effect['count']}篇相关研究）")

    report_lines.append("\n⚠ 风险提示：")

    if toxicity_data['human']:
        human_effects = [e['type'] for e in toxicity_data['human']]
        if 'Cardiotoxicity' in human_effects:
            report_lines.append("- 该成分对人类显示心脏毒性，临床使用需特别谨慎！")
        if 'Hepatotoxicity' in human_effects:
            report_lines.append("- 存在明确的肝毒性风险，长期使用需监测肝功能指标")
    else:
        report_lines.append("- 缺乏人类直接毒性数据，建议参考动物实验数据评估风险")

    if toxicity_data['rat'] or toxicity_data['mice']:
        report_lines.append("- 动物实验显示多系统毒性，建议进行更全面的安全性评估")

    return "\n".join(report_lines)


def generate_all_reports(df, output_file=None):
    """为整个数据框生成报告"""
    full_report = []

    for _, row in df.iterrows():
        if pd.notna(row['toxic_effect']):
            report = generate_toxicity_report(
                row['component_name'],
                row['toxic_effect']
            )
            full_report.append(report)
            full_report.append("\n" + "=" * 80 + "\n")

    full_report_text = "\n".join(full_report)

    if output_file:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(full_report_text)

    return full_report_text


if __name__ == '__main__':
    df = pd.read_excel('Data/Toxicity/中药.xlsx')  # 加载数据
    report = generate_all_reports(df, 'toxicity_report.txt')