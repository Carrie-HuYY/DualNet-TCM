import json

'''
1.将symbol转化为target

2.将target和其研究阶段进行对应

3.将target和药物及其研究阶段进行对应

'''

# 将uniprotID转化为target
# 输入 P2-01-TTD_uniprot_all.txt的路径
# 输出 UniprotID_To_Target.json
def uniprotID_To_Target(TTD_uniprot_all_txt):
    with open ('P2-01-TTD_uniprot_all.txt','r') as f:
        target_info = f.read().split('\n\n')[2].split('\n\t')
    UniprotID_To_Target = {}
    for i in target_info:
        '''
        通过flag来控制是否将uniprotID和target对应起来
        '''
        flag = 0
        ta_1 = i.split('\n')
        for i in ta_1:
            if 'TARGETID' in i:
                targetid = i.split('\t')[2]
                # 保证uniprotID是来自于人类的
            elif 'UNIPROID' in i and 'HUMAN' in i:
                uniproid = i.split('\t')[2].split('_')[0]
                flag += 1
                # 将mRNA去掉
            elif 'TARGNAME' in i and 'mRNA' not in i:
                flag += 1
            elif 'TARGTYPE' in i :
                targtype = i.split('\t')[2]
        if flag == 2:                                                                                          
            UniprotID_To_Target[uniproid] = {targetid : targtype}

    with open("UniprotID_To_Target.json", "w", encoding='utf-8') as f:
        f.write(json.dumps(UniprotID_To_Target, ensure_ascii=False, indent=4, separators=(',', ':')))

# 将symbol转化为target
def symbol_to_target():
    with open ('../ID_Transformed/Symbol_To_UniprotID.json', 'r') as f:
        Symbol_To_UniprotID = json.load(f)
        
    with open ('UniprotID_To_Target.json', 'r') as f:
        UniprotID_To_Target = json.load(f)


    Symbol_To_Target = {}
    uniprotID_list = []
    for symbol in Symbol_To_UniprotID.keys():
        '''
        先将symbol转化成uniprotID,再将uniprotID转化成target
        '''
        uniprotID = Symbol_To_UniprotID[symbol]
        if uniprotID in UniprotID_To_Target.keys():
            Symbol_To_Target[symbol] = UniprotID_To_Target[uniprotID]
            # 查看有多少个uniprotID是对应target的
            uniprotID_list.append(uniprotID)

    with open("Symbol_To_Target.json", "w", encoding='utf-8') as f:
        f.write(json.dumps(Symbol_To_Target, ensure_ascii=False, indent=4, separators=(',', ':')))

# 存在部分靶标没有对应的uniprotID 可能不是人源蛋白
def tar_uni_sym():
    with open ('P1-01-TTD_target_download.txt','r') as f:
        drug_info = f.read().split('\n\n')[2].split('\n')

    ta_info = {}
    for line in drug_info:
        # 其中有一行都是\t需要删去
        if line != '\t\t\t\t':
            if line[:6] not in ta_info.keys():
                ta_info[line[:6]] = []
                ta_info[line[:6]].append(line)
            else:
                ta_info[line[:6]].append(line)

    # {'T47101':{'uniprotid':'','symbol':''}}
    # 有的靶标没有uniprotid,有的靶标没有symbol
    ta_uni_sym = {}
    for key,value in ta_info.items():
        ta_uni_sym[key] = {'uniprotid':'','symbol':''}
        uni_num = 0
        sym_num = 0
        for i in value:
            if 'UNIPROID' in i:
                uniproid = i.split('\t')[2]
                print(uniproid)
                uni_num = 1
            elif 'GENENAME' in i:
                symbol = i.split('\t')[2]
                sym_num = 1
            if uni_num == 1 and sym_num == 1:
                ta_uni_sym[key]['uniprotid'] = uniproid
                ta_uni_sym[key]['symbol'] = symbol
                break
            if uni_num == 1 and sym_num == 0:
                ta_uni_sym[key]['uniprotid'] = uniproid
                ta_uni_sym[key]['symbol'] = ''
            if uni_num == 0 and sym_num == 1:
                ta_uni_sym[key]['uniprotid'] = ''
                ta_uni_sym[key]['symbol'] = symbol
    with open("Target_To_UNIandSYM.json", "w", encoding='utf-8') as f:
        f.write(json.dumps(ta_uni_sym, ensure_ascii=False, indent=4, separators=(',', ':')))


'''
有的target是没有药的,需要先把所有的targetid找出来
再通过targetid找到对应的药物的信息

'''
# 将target和药物及其研究阶段进行对应
def get_target_drug():
    with open ('P1-01-TTD_target_download.txt','r') as f:
        drug_info = f.read().split('\n\n')[2].split('\n')
    UniprotID_To_Target2 = {}
    Target_To_Drug = {}
    for line in drug_info:
        if 'TARGETID' in line:
            targetid = line.split('\t')[2]
            print(targetid)
            if targetid not in Target_To_Drug.keys():
                Target_To_Drug[targetid] = []
            flag = 0
        elif 'DRUGINFO' in line:
            drug_name = line.split('\t')[-2]
            drug_phase = line.split('\t')[-1]
            Target_To_Drug[targetid].append({drug_name : drug_phase})
    with open("Target_To_Drug.json", "w", encoding='utf-8') as f:
        f.write(json.dumps(Target_To_Drug, ensure_ascii=False, indent=4, separators=(',', ':')))

        






        
        
        
        
        

        
        
        
        
