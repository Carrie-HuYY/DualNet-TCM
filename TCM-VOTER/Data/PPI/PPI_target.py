import pandas as pd
import json
import re
'''
输出一个与差异表达蛋白存在相互作用的蛋白列表

'''


'''
1.处理STRING数据库中的PPI数据

'''
def get_PPI_dict():
    '''
    将ENSP转换为uniprotID
    '''
    def ENSPToUniprotId():
        with open (r'uniprot_sprot_human.dat', 'r') as f:
            f = f.read()
            # 通过// 将每个uniprotID的信息分开
            f = f.split('\n//')
        ENSP_to_uniprotId = {}
        for i in f:
            f_list = i.split('\n')
            flag = False
            for i in f_list:
                if i.startswith('ID   '):
                    value = i.split()[1].split('_')[0]
                elif i.startswith('DR   STRING'):
                    key = i.split(';')[1].replace(' ', '')
                    flag = True
                if flag:
                    ENSP_to_uniprotId[key] = value
                    break
        return ENSP_to_uniprotId
    # 将STRING数据库的PPI数据转换为字典保存
    ENSP_to_uniprotId = ENSPToUniprotId()       
    PPI_DF=pd.read_csv(r'9606.protein.physical.links.full.v11.5.txt', sep=' ')
    # 保留相互作用强度大于700的蛋白质对，因为大于这个数值就可以认为是直接相互作用的
    PART_PPI_DF = PPI_DF.loc[PPI_DF['combined_score'] > 700][['protein1','protein2']]
    protein1_list = PART_PPI_DF['protein1'].tolist()
    protein2_list = PART_PPI_DF['protein2'].tolist()

    PPI_DICT = {}
    for num in range(len(protein1_list)):
        '''
        将在uniprot中具有uniprotId的ENSP转换为uniprotId
        其余的ENSP不进行转换
        '''
        if protein1_list[num] in ENSP_to_uniprotId.keys() and protein2_list[num] in ENSP_to_uniprotId.keys():
            symbol_1 = ENSP_to_uniprotId[protein1_list[num]]
            symbol_2 = ENSP_to_uniprotId[protein2_list[num]]
            if symbol_1 not in PPI_DICT.keys():
                PPI_DICT[symbol_1] = [symbol_2]
            else:
                PPI_DICT[symbol_1].append(symbol_2)
    return PPI_DICT


# 将字典转化为json格式
# PPI_DICT = get_PPI_dict()
# with open("PPI.json", "w") as f:
#     f.write(json.dumps(PPI_DICT, ensure_ascii=False, indent=4, separators=(',', ':')))

# 把差异表达蛋白对应的每个PPi的蛋白列表进行组合,并统计其中每个蛋白出现的次数就可以获得
# 哪些蛋白与这些差异表达蛋白有较强的关联性,可以作为靶标进行分析

def get_symbol_list():
    with open (r'PDPM_Protein_List.txt') as f:
        symbol_list = []
        for line in f:
            line = line.replace('\t', ' ').replace('\n', '').split(' ')
            symbol_list.append(line[0])
    return symbol_list
        

# 将组学数据中与每个差异表达蛋白相互作用的蛋白进行汇总，统计每个蛋白出现的次数
def get_targetNum_dict(symbol_list, interaction_num, PPI_DICT):
    '''
    将组学数据中与每个差异表达蛋白相互作用的蛋白进行汇总，统计每个蛋白出现的次数
    根据interaction_num的值，筛选出与差异表达蛋白相互作用次数大于等于interaction_num的蛋白
    '''
    ALL_PPI_PROTEIN = []

    for i in symbol_list:
        if i in PPI_DICT.keys():
            ALL_PPI_PROTEIN.extend(PPI_DICT[i])

    PPI_NUMBER = {}
    TARGET_PPI = {}

    for i in ALL_PPI_PROTEIN:
        if i not in PPI_NUMBER.keys():
            PPI_NUMBER[i] = 1
        else:
            PPI_NUMBER[i] += 1
            if PPI_NUMBER[i] >= interaction_num and i not in TARGET_PPI.keys():
                TARGET_PPI[i] = PPI_NUMBER[i]
            elif PPI_NUMBER[i] >= interaction_num and i in TARGET_PPI.keys():
                TARGET_PPI[i] = PPI_NUMBER[i]

    return TARGET_PPI

# 将geneID转换为pubmedID
def geneIDToSymbol():
    GeneToPubMed = pd.read_csv('gene2pubmed', sep='\t')

    HuMan_DF = GeneToPubMed.loc[GeneToPubMed['#tax_id'] == 9606][['GeneID', 'PubMed_ID']]


    HuMan_geneId_list = HuMan_DF['GeneID'].tolist()
    HuMan_PubMedId_list = HuMan_DF['PubMed_ID'].tolist()

    GeneID_To_Symbol = {}
    for i in range(len(HuMan_geneId_list)):
        if HuMan_geneId_list[i] not in GeneID_To_Symbol.keys():
            GeneID_To_Symbol[HuMan_geneId_list[i]] = [HuMan_PubMedId_list[i]]
        else:
            GeneID_To_Symbol[HuMan_geneId_list[i]].append(HuMan_PubMedId_list[i])


    with open("GeneID_To_PubMedID.json", "w") as f:
        f.write(json.dumps(GeneID_To_Symbol, ensure_ascii=False, indent=4, separators=(',', ':')))

# 将symbol转换为geneID
def symbolToGeneID():
    ALL_DF = pd.read_csv(r'protein-coding_gene.txt', sep='\t')
    symbolToGeneID_DF = ALL_DF[['symbol', 'entrez_id']]
    symbol_list = symbolToGeneID_DF['symbol'].tolist()
    GeneID_list = symbolToGeneID_DF['entrez_id'].tolist()
    symbolToGeneID_dict = {}
    for i in range(len(symbol_list)):
        if symbol_list[i] not in symbolToGeneID_dict.keys():
            symbolToGeneID_dict[symbol_list[i]] = GeneID_list[i]
    with open("Symbol_To_GeneID.json", "w") as f:
        f.write(json.dumps(symbolToGeneID_dict, ensure_ascii=False, indent=4, separators=(',', ':')))



# 将每个symbol对应的uniprotID进行汇总
def symbolToUniProtID():
    with open (r'uniprot_sprot_human.dat', 'r') as f:
        f = f.read()
        # 通过// 将每个uniprotID的信息分开
        f = f.split('\n//')
    symbol_to_uniprotId = {}
    for i in f:
        f_list = i.split('\n')
        flag = False
        for i in f_list:
            if i.startswith('ID   '):
                value = i.split()[1].split('_')[0]
            elif i.startswith('GN   Name=') and i.endswith(';'):
                key = i.split(';')[0].split('=')[1].split()[0]
                flag = True
            elif i.startswith('GN   Name=') and not i.endswith(';'):
                key = i.split('=')[1].split()[0]
            if flag:
                symbol_to_uniprotId[key] = value
                break
    with open("Symbol_To_UniProtID.json", "w") as f:
        f.write(json.dumps(symbol_to_uniprotId, ensure_ascii=False, indent=4, separators=(',', ':')))

