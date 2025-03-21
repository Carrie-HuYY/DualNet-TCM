# DualNet-TCM(中药整合药理学数据库)

***DualNet-TCM***(TCM (中药) + Dual (双重网络) + Net (网络药理学)):
从网络药理学和中医理论双重视角出发，构建"方剂-中药-成分-靶点-病症"和"方剂-对症-症状-病证"双重关联网络。
同时，结合基因研究和文献计量学方法，开展转化医学研究：通过药物可用性评估和文献挖掘技术，对大量候选蛋白或基因进行系统筛选和分类，
从而精准缩小需要实验验证的目标范围，为中药现代化研究提供新的思路和方法。



- [1. 简介](#简介)

- [2. 安装](#安装)

- [3. 使用](#使用)
  - [3.1 from_SD：从辩证出发](#from_sd)
  - [3.2 from_TCM_or_Formula：从中药/方剂出发](#from_tcm_or_formula)
  - [3.3 from_proteins：从靶点出发](#from_protein)

- [4. 结果展示](#结果展示)
  - [4.1 基于pyecharts的可交互图](#out_graph)
  - [4.2 可用于cytoscape绘图的文件](#out_for_cytoscape)
  - [4.3 保存结果的excel表格](#out_for_excel)
  - [4.4 靶点研究进展的查询结果](#research_status_test)
  - [4.5 安全性评估结果](#safety_research)

- [5. 注意](#注意)
  - [5.1 关于Python版本问题](#1关于python版本)
  - [5.2 关于数据下载问题](#2关于数据下载)

## 简介

## 安装

可以使用pip安装DualNet-TCM

```cmd
pip install DualNet-TCM
```

requirement包括(Python版本需要≥3.9，其余没有强制要求)

```cmd
pandas~=1.3.5
rdkit~=2023.3.2
tqdm~=4.67.1
pyecharts~=2.0.7
numpy~=1.21.6
```

## 使用

### `from_SD`
从中医辩证开始，构建中药整合药理学

```python
from DualNet-TCM import analysis
analysis.from_SD(
  SD, 
  score=990, 
  out_graph=False, 
  out_for_cytoscape=False, 
  out_for_excel=True,
  research_status_test=False,
  safety_research=False,
  re=True, 
  path='results/'
)
```

`from_SD`需要一个必需形参SD：任何可以使用in判断一个元素是否在辩证数据集中的组合数据类型，
存储拟分析的辩证的ID， 如['DNS001'],具体ID编号及SD信息可见于[SD辩证数据集](/DualNet-TCM/Data/SD.xlsx)

`from_SD`的可选形参有：
- `score`: int类型，[Chemical_Protein_Links数据集](/DualNet-TCM/Data/Chemical_Protein_Links.xlsx)
中仅combined_score大于等于score的记录会被筛选出，默认为990；
- `out_graph`: boolean类型，是否输出基于ECharts的html格式的网络可视化图，默认为`False`；
- `out_for_cytoscape`: boolean类型，是否输出用于Cytoscape绘图的文件，默认为`False`；
- `our_for_excel`: boolean类型，是否将结果输出到excel中，默认为`True`；
- `research_status_test`: boolean类型，是否将得到的靶点进行研究价值评定，默认为`False`；
- `safety_research`: boolean类型，是否对得到的靶点进行安全性研究，默认为`False`；
- `re`: boolean类型，是否返回原始分析结果（辩证、复方、中药、化合物（中药成分）、蛋白（靶点）及其连接信息），
默认为True。若re为True，则函数将返回运行结果sd、sd_formula_links.xlsx, formula、formula_tcm_links、tcm、tcm_chem_links、chem、
chem_protein_links和proteins，它们均为pd.DataFrame类型，分别存储了辩证信息、辩证-复方连接信息，
复方信息、复方-中药连接信息、中药信息、中药-化合物（中药成分）连接信息、化合物（中药成分）信息、 化合物（中药成分）-蛋白（靶点）连接信息和蛋白（靶点）信息；
- `path`: str类型，存放结果的路径，默认为`results/`。若无此路径，将自动建立相应的目录。

### `from_TCM_or_Formula`
从中药或者方剂开始，构建中药整合药理学网络

```python
from DualNet-TCM import analysis
analysis.from_TCM_or_Formula(
  tcm_or_formula, 
  score=990, 
  out_graph=False, 
  out_for_cytoscape=False, 
  out_for_excel=True,
  Research_status_test=False,
  safety_research=False,
  re=True, 
  path='results/'
)
```

`from_TCM_or_Formula`需要一个必需形参tcm_or_formula：任何可以使用in判断一个元素是否在其中的组合数据类型，
存储拟分析的中药或复方的ID， 如['DNH0367', 'DNH1695']，具体ID编号及中药/方剂信息可见于[方剂数据集](/DualNet-TCM/Data/Formula.xlsx)
和[中药数据集](/DualNet-TCM/Data/Herb.xlsx)

`from_TCM_or_Formula`的可选形参有：
- `score`: int类型，[Chemical_Protein_Links数据集](/DualNet-TCM/Data/Chemical_Protein_Links.xlsx)
中仅combined_score大于等于score的记录会被筛选出，默认为990；
- `out_graph`: boolean类型，是否输出基于ECharts的html格式的网络可视化图，默认为`False`；
- `out_for_cytoscape`: boolean类型，是否输出用于Cytoscape绘图的文件，默认为`False`；
- `our_for_excel`: boolean类型，是否将结果输出到excel中，默认为`True`；
- `research_status_test`: boolean类型，是否将得到的靶点进行研究价值评定，默认为`False`；
- `safety_research`: boolean类型，是否对得到的靶点进行安全性研究，默认为`False`；
- `re`: boolean类型，是否返回原始分析结果（辩证、复方、中药、化合物（中药成分）、蛋白（靶点）及其连接信息），
默认为True。若re为True，则函数将返回运行结果sd、sd_formula_links.xlsx, formula、formula_tcm_links、tcm、tcm_chem_links、chem、
chem_protein_links和proteins，它们均为pd.DataFrame类型，分别存储了辩证信息、辩证-复方连接信息，
复方信息、复方-中药连接信息、中药信息、中药-化合物（中药成分）连接信息、化合物（中药成分）信息、 化合物（中药成分）-蛋白（靶点）连接信息和蛋白（靶点）信息；
- `path`: str类型，存放结果的路径，默认为`results/`。若无此路径，将自动建立相应的目录。

### `from_Protein`
从靶点开始，构建中药整合药理学

```python
from DualNet-TCM import analysis
analysis.from_Protein(
  protein,
  score=0,
  random_state=None,
  num=1000, 
  tcm_component=True, 
  formula_component=True,
  out_graph=False, 
  out_for_cytoscape=False, 
  out_for_excel=True,
  Research_status_test=False,
  safety_research=False,
  re=True, 
  path='results/'
)
```

`from_Protein`需要一个必需形参protein：任何可以使用in判断一个元素是否在其中的组合数据类型，
存储拟分析的中药或复方的ID， 如['ENSP00000381588', 'ENSP00000252519']，具体ID编号及靶点信息可见于[靶点数据集](/DualNet-TCM/Data/Protein.xlsx)

`from_Protein`的可选形参有：
- `score`: int类型，[Chemical_Protein_Links数据集](/DualNet-TCM/Data/Chemical_Protein_Links.xlsx)
中仅combined_score大于等于score的记录会被筛选出，默认为`0`；
- `random_state`: int类型，指定优化模型使用的随机数种子，默认为`None`，即不指定随机数种子；
- `num`: int类型，指定优化时需生成的解的组数，默认为`1000`；
- `tcm_component`: boolean类型，是否进行中药组合优化，默认为`True`；
- `formula_component`: boolean类型，是否进行复方组合优化，默认为`True`；
- `out_graph`: boolean类型，是否输出基于ECharts的html格式的网络可视化图，默认为`False`；
- `out_for_cytoscape`: boolean类型，是否输出用于Cytoscape绘图的文件，默认为`False`；
- `our_for_excel`: boolean类型，是否将结果输出到excel中，默认为`True`；
- `research_status_test`: boolean类型，是否将得到的靶点进行研究价值评定，默认为`False`；
- `safety_research`: boolean类型，是否对得到的靶点进行安全性研究，默认为`False`；
- `re`: boolean类型，是否返回原始分析结果（辩证、复方、中药、化合物（中药成分）、蛋白（靶点）及其连接信息），
默认为True。若re为True，则函数将返回运行结果sd、sd_formula_links.xlsx, formula、formula_tcm_links、tcm、tcm_chem_links、chem、
chem_protein_links和proteins，它们均为pd.DataFrame类型，分别存储了辩证信息、辩证-复方连接信息， 复方信息、复方-中药连接信息、中药信息、
- 中药-化合物（中药成分）连接信息、化合物（中药成分）信息、 化合物（中药成分）-蛋白（靶点）连接信息和蛋白（靶点）信息；
- `path`: str类型，存放结果的路径，默认为`results/`。若无此路径，将自动建立相应的目录。

## 结果展示

### out_graph

### out_for_cytoscape

### out_for_excel

### research_status_test

### safety_research

## Star History

[![Star History Chart](https://api.star-history.com/svg?repos=Carrie-HuYY/DualNet-TCM&type=Date)](https://star-history.com/#Carrie-HuYY/DualNet-TCM&Date)

## 注意

### 1.关于Python版本

由于在 Python 3.9 之前的版本中，`tuple[...]` 和 `list[...]` 
这样的类型注解语法不被支持。
Python 3.9 引入了原生的类型注解支持（PEP 585），
但在早期版本中，需要使用 typing 模块中的 Tuple、List 等类型. 所以需要Python≥3.9，如果＜3.9的话，
可将`compute.py`函数中修改如下：

```python
from typing import Tuple, Union

def score(weights: Union[dict, None] = None) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # 函数逻辑
    pass
```

### 2.关于数据下载
整体大小为9个G，由于百度网盘限制，所以拆分成三个压缩包，解压后放data/文件夹即可

[下载链接1](https://pan.baidu.com/s/1zIlTjstJMscKdZnP30wc1g?pwd=2n2t) 

[下载链接2](https://pan.baidu.com/s/1tg8WQtJiJi70A8HIRYG_PA?pwd=9bvh) 

[下载链接3](https://pan.baidu.com/s/1tg8WQtJiJi70A8HIRYG_PA?pwd=9bvh)
