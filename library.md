## NBSD

多模态药物靶标寻找工具

- 文本：中医古籍、《中草药》、《南中医学报》、pubmed
- 基因：基因表达谱
- 中药资源分布：
- String数据库：PPI网络图

### 2024.11.17 

| Dataset                | Download Address                                                                          | Update Time |
|------------------------|-------------------------------------------------------------------------------------------|-------------|
| String PPI information | https://stringdb-downloads.org/download/protein.physical.links.v12.0.txt.gz               | 2024.11.17  |        
| TCM Classics           | https://github.com/xiaopangxia/TCM-Ancient-Books.git                                      | 2021.5.6    |
| TTD targets            | https://db.idrblab.net/ttd/sites/default/files/ttd_database/P1-01-TTD_target_download.txt | 2024.01.10  |
| TTD drug               | https://db.idrblab.net/ttd/sites/default/files/ttd_database/P1-02-TTD_drug_download.txt   | 2024.01.10  |
| Human uniprot          |https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/pan_proteomes/| 2024.10.2   |
| TCM formuals dataset   |https://www.nmpa.gov.cn/| 2024.10.2   |

### 2024.11.19

编写爬虫代码爬取ETCM上所有的方剂，主要是通过正则表达式的方式

### 2024.11.23

成功爬取所有方剂及其组合，一共3959条，保存在Data_ETCM.xlsx中


### 2024.11.26

听澳门大学欧阳德方教授的讲座，思考能不能加人体微分方程模拟药物代谢

### 2024.11.27

获取包含54152条中医辨证数据集，涵盖148个综合证，训练一个辩证预测的样本

### 2024.11.29

我将中医辩证数据集的分类定义成多标签分类任务，格式如下

| 标题                |描述| 具体样例           |
|-------------------|---|----------------|
| User_id           |病人ID| 479372         |
| ICD_ID & ICD_name |国际疾病分类| BNP120吐血病      |
| Norm_syndrome     |病人证型| 气虚不摄证          |
| chief_complaint   |主诉| 呕血2小时          |
| description       |病人病史描述| 患者2小时前无明显诱因出现吐鲜血，色红，反酸、嗳气、烧心偶作，偶有恶心干呕， |
| dection           |中医四诊信息|中医四诊摘要：神志清晰，|

### 2024.11.30

首先把模型理解成按照主诉、病人病史描述和中医四诊信息，来预测国际疾病分类和病人证型的流程。
首先把lcd_name、description、chief_complaint和detection作为问题输入，syndrome和国际疾病分类作为标签，整理成数据库格式