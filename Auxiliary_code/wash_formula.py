import pandas as pd
import time

# 读取中药名和方剂数据
Herb = pd.read_excel("D:/NBSD/Data/Tcm/Herb.xlsx")
Herb_set = set(Herb['Herb_cn_name'])  # 将中药名转换为集合，提高查找效率

Formula = pd.read_excel("D:/NBSD/Data/Tcm/Formula.xlsx")
Formula_name = Formula['方剂']
Formula_components = Formula['成分']

# 初始化结果列表
result_data = []

# 开始计时
start_time = time.time()

# 定义一个函数，用于提取文本中的中药名
def extract_herbs(text, herb_set):
    found_herbs = set()
    for herb in herb_set:
        if isinstance(herb, str) and isinstance(text, str):
            if herb in text:
                found_herbs.add(herb)
    return found_herbs

# 遍历所有方剂数据
for i in range(len(Formula_name)):  # 动态处理所有数据
    formula_name = str(Formula_name[i])
    text = str(Formula_components[i])
    found_herbs = extract_herbs(text, Herb_set)  # 调用函数提取中药名
    result_data.append({'Formula': formula_name, 'Found_Herbs': ', '.join(found_herbs)})

# 将结果转换为 DataFrame
result_df = pd.DataFrame(result_data)

# 展平成分列
expanded_df = result_df.assign(Found_Herbs=result_df['Found_Herbs'].str.split(", ")).explode('Found_Herbs')

# 结束计时
end_time = time.time()

# 保存结果到 Excel 文件
expanded_df.to_excel("D:/NBSD/Data/Tcm/Result.xlsx", index=False)

# 打印运行时间和保存信息
print(f"运行时间：{end_time - start_time:.2f} 秒")
print("结果已保存到 D:/NBSD/Data/Tcm/Result.xlsx")