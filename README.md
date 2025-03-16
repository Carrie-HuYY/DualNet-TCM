## DualNet-TCM

***DualNet-TCM***(TCM (中药) + Dual (双重网络) + Net (网络药理学)):
从网络药理学和中医理论双重视角出发，构建"方剂-中药-成分-靶点-病症"和"方剂-对症-症状-病证"双重关联网络。
同时，结合基因研究和文献计量学方法，开展转化医学研究：通过药物可用性评估和文献挖掘技术，对大量候选蛋白或基因进行系统筛选和分类，
从而精准缩小需要实验验证的目标范围，为中药现代化研究提供新的思路和方法。

- [安装](#安装)

- [使用](#使用)

- [环境要求](#requirements)

- [致谢](#致谢)

### 安装

### 注意

#### 1.关于Python版本

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

#### 2.数据下载 
整体大小为9个G，由于百度网盘限制，所以拆分成三个压缩包，解压后放data/文件夹即可

[下载链接1](https://pan.baidu.com/s/1zIlTjstJMscKdZnP30wc1g?pwd=2n2t) 

[下载链接2](https://pan.baidu.com/s/1tg8WQtJiJi70A8HIRYG_PA?pwd=9bvh) 

[下载链接3](https://pan.baidu.com/s/1tg8WQtJiJi70A8HIRYG_PA?pwd=9bvh)

### 使用

### Requirements

- pandas
- pyecharts
- numpy
- tqdm
- requests
- os

### Star History

[![Star History Chart](https://api.star-history.com/svg?repos=Carrie-HuYY/DualNet-TCM&type=Date)](https://star-history.com/#Carrie-HuYY/DualNet-TCM&Date)
