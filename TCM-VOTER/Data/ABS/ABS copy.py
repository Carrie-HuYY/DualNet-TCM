import xml.etree.cElementTree as ET
import json
import threading

# 处理一个xml文件需要30s或者更多，所以处理1500个文件需要较长的时间
# 需要开多个线程进行处理
# 若单一线程则需要将近10个小时才能完成（可能程序有点问题）
# 多线程处理100文件需要30min左右



# 将摘要处理成json格式   
def abstract_to_json():
    # 获取所有摘要的路径
    ABS_PH_list = []
    for i in range(1,1552,1):
        num = str(i).zfill(4)
        ABS_PH = 'pubmed22n' + num 
        ABS_PH_list.append(ABS_PH)

    for path in ABS_PH_list[100:200]:
        print(path)
        tree = ET.parse('G:/ABS/'+ path + '.xml') 
        PubmedArticle = tree.getroot().findall('PubmedArticle')
        pubmedid_to_abstract = {}
        for i in PubmedArticle:
            flag = False
            MedlineCitation = i.find('MedlineCitation')
            for i in MedlineCitation.find('Article'):
                if 'Abstract' in i.tag:
                    flag = True
            if flag == True:     
                Abstract = ''   
                PMID = MedlineCitation.find('PMID').text
                for i in MedlineCitation.find('Article').find('Abstract').findall('AbstractText'):
                    if i.text is not None:
                        # 若存在其他标签则将标签中的内容也加入到摘要中
                        for subelem in i.iter():
                            if subelem.text:
                                Abstract += subelem.text.strip()
                            if subelem.tail:
                                Abstract += subelem.tail.strip()
                pubmedid_to_abstract[PMID] = Abstract
            elif flag == False:
                PMID = MedlineCitation.find('PMID').text
                Abstract = ''
                pubmedid_to_abstract[PMID] = Abstract

        with open('G:/ABS/ABS_json/' + path + ".json", "w", encoding='utf-8') as f:
            f.write(json.dumps(pubmedid_to_abstract, ensure_ascii=False, indent=4, separators=(',', ':')))
        f.close

# 将摘要保存到elasticsearch中
abstract_to_json()



    
    
    
    
    
    
    
    
    














import os

# Define the directory path where the txt files are located
directory_path = "G:/ABS/ABS_json/"
file_list = []
# Loop through all files in the directory
for filename in os.listdir(directory_path):
    if filename.endswith(".json"):
        # Open the file and count the number of lines
        with open(os.path.join(directory_path, filename), "r") as file:
            num_lines = sum(1 for line in file)
        # Check if the number of lines is less than 30000 and print the filename if it is
        if num_lines < 20000:
            print(filename)
            file_list.append(filename)

file_list = [name[:-5] for name in file_list ]
ABS_PH_list = []
for i in range(1,1552,1):
    num = str(i).zfill(4)
    ABS_PH = 'pubmed22n' + num 
    ABS_PH_list.append(ABS_PH)

for path in file_list:
    print(path)
    tree = ET.parse('G:/ABS/'+ path + '.xml') 
    PubmedArticle = tree.getroot().findall('PubmedArticle')
    pubmedid_to_abstract = {}
    for i in PubmedArticle:
        flag = False
        MedlineCitation = i.find('MedlineCitation')
        for i in MedlineCitation.find('Article'):
            if 'Abstract' in i.tag:
                flag = True
        if flag == True:     
            Abstract = ''   
            PMID = MedlineCitation.find('PMID').text
            for i in MedlineCitation.find('Article').find('Abstract').findall('AbstractText'):
                if i.text is not None:
                    # 若存在其他标签则将标签中的内容也加入到摘要中
                    for subelem in i.iter():
                        if subelem.text:
                            Abstract += subelem.text.strip()
                        if subelem.tail:
                            Abstract += subelem.tail.strip()
            pubmedid_to_abstract[PMID] = Abstract
        elif flag == False:
            PMID = MedlineCitation.find('PMID').text
            Abstract = ''
            pubmedid_to_abstract[PMID] = Abstract

    with open('G:/ABS/ABS_json/' + path + ".json", "w", encoding='utf-8') as f:
        f.write(json.dumps(pubmedid_to_abstract, ensure_ascii=False, indent=4, separators=(',', ':')))
    f.close


