import pandas as pd

# 读取原始数据
data = pd.read_excel("D:/NBSD/Data/TCM/SD.xlsx")

# 选择需要的列
data = data[['证候', '对应方剂']]

# 初始化结果列表
result = []

# 遍历数据框的每一行
for index, row in data.iterrows():
    zhengxing = row['证候']  # 证型
    tangzheng = row['对应方剂']  # 汤证部分

    # 检查汤证部分是否为空
    if pd.notnull(tangzheng):
        # 将汤证部分按逗号分割，去掉多余的空格
        tangzheng_list = [t.strip() for t in tangzheng.split(',')]

        # 将每个汤证与对应的证型组合并添加到结果中
        for tang in tangzheng_list:
            result.append([zhengxing, tang])

# 将结果转换为数据框
result_df = pd.DataFrame(result, columns=['证候', '对应方剂'])

# 保存结果到新的 Excel 文件
result_df.to_excel("D:/NBSD/Data/TCM/SD_Processed.xlsx", index=False)

print("处理完成，结果已保存到 SD_Processed.xlsx 文件中。")