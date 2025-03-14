import os
import re
import sys
import time

import pandas as pd
from pubchempy import get_compounds
from rdkit import Chem
from rdkit import RDLogger
from tqdm import tqdm
from retrying import retry

sys.path.append(os.path.abspath('../src'))

# 读取数据
tcmbank_ingredient = pd.read_csv(
    'ingredient_all.CSV',
    encoding='ISO-8859-1',
    dtype=object,
    usecols=[
        'TCMBank_ID', 'name', 'Alias', 'Molecular_Formula', 'Smiles', 'Molecular_Weight',
        'Molecular_Volume', 'Ingredient_id', 'CAS_id', 'SymMap_id', 'TCMID_id', 'TCMSP_id',
        'TCM-ID_id', 'PubChem_id', 'DrugBank_id'
    ]
)
tcmbank_ingredient.rename(columns={'name': 'Ingredient_name', 'Smiles': 'Ingredient_Smile'}, inplace=True)

# 打印列名
print(tcmbank_ingredient.columns)

# 禁用 RDKit 日志
RDLogger.DisableLog('rdApp.*')


def smiles_standardization(smiles):
    """
    标准化 SMILES 字符串。
    """
    if pd.isna(smiles) or not isinstance(smiles, str):
        return None  # 如果输入为空或非字符串，返回 None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None  # 如果 SMILES 无效，返回 None

    # 标准化分子
    mol = Chem.RemoveHs(mol)  # 去除氢原子
    mol = Chem.AddHs(mol)  # 添加显式氢原子
    standardized_smiles = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
    return standardized_smiles


def is_smiles(smiles):
    """
    判断输入是否为有效的 SMILES 字符串。
    """
    if pd.isna(smiles) or not isinstance(smiles, str):
        return False  # 如果输入为空或非字符串，返回 False

    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False


@retry(
    stop_max_attempt_number=5,  # 最大重试次数
    wait_exponential_multiplier=100,  # 初始等待时间（毫秒）
    wait_exponential_max=1000  # 最大等待时间（毫秒）
)
def get_from_pubchem(compound):
    """
    根据 CAS 号、PubChem_id、名称、别名等信息查找 PubChem 数据库中对应的化合物信息。
    """
    # 先根据 CAS 号查找
    cas_id = compound['CAS_id']
    if isinstance(cas_id, str) and cas_id.strip():
        com_obj_ls = get_compounds(cas_id, 'name')
        if len(com_obj_ls) == 1:
            return com_obj_ls[0]

    # 没有则根据 PubChem_id 查找
    pubchem_id = compound['PubChem_id']
    if isinstance(pubchem_id, str) and pubchem_id.strip():
        try:
            com_obj_ls = get_compounds(int(pubchem_id), 'cid')
            if len(com_obj_ls) == 1:
                return com_obj_ls[0]
        except ValueError:
            pass  # 如果 PubChem_id 不是数字，跳过

    # 没有则根据名称查找
    name = compound['Ingredient_name']
    if isinstance(name, str) and name.strip():
        com_obj_ls = get_compounds(name, 'name')
        if len(com_obj_ls) == 1:
            return com_obj_ls[0]

    # 还可以根据别名查找
    alias = compound['Alias']
    if isinstance(alias, str) and alias.strip():
        for alias_name in alias.split(';'):
            com_obj_ls = get_compounds(alias_name.strip(), 'name')
            if len(com_obj_ls) == 1:
                return com_obj_ls[0]
    return None


# 标准化 SMILES 列
tcmbank_ingredient['Ingredient_Smile'] = tcmbank_ingredient['Ingredient_Smile'].apply(smiles_standardization)

# 保存标准化后的数据
tcmbank_ingredient.to_csv('tcmbank_ingredient.csv', index=False)

# 初始化统计变量
origin_smile_count, smile_check = 0, 0
unmatched_smiles_rows = []

# 遍历数据框，进行 PubChem 查询和 SMILES 匹配
for index, row in tqdm(tcmbank_ingredient.iterrows(), total=tcmbank_ingredient.shape[0]):
    try:
        compound_obj = get_from_pubchem(row)
        if compound_obj is not None:
            origin_smile = row['Ingredient_Smile']
            if is_smiles(origin_smile):
                origin_smile_count += 1

                # 如果是有立体异构信息的 SMILES，则和 isomeric_smiles 比较，否则和 canonical_smiles 比较
                if '@' in origin_smile or '\\' in origin_smile or '/' in origin_smile:
                    if smiles_standardization(compound_obj.isomeric_smiles) == origin_smile:
                        smile_check += 1
                elif smiles_standardization(compound_obj.canonical_smiles) == origin_smile:
                    smile_check += 1
                else:
                    # 匹配不上则记录
                    row['pubchempy_smiles'] = smiles_standardization(compound_obj.canonical_smiles)
                    row['pubchempy_isomeric_smiles'] = smiles_standardization(compound_obj.isomeric_smiles)
                    unmatched_smiles_rows.append(row)
                    continue

            # 更新 PubChem_id 和 SMILES
            tcmbank_ingredient.loc[index, 'PubChem_id'] = compound_obj.cid
            tcmbank_ingredient.loc[index, 'Ingredient_Smile'] = smiles_standardization(compound_obj.canonical_smiles)
            tcmbank_ingredient.loc[index, 'Isoneric_Smiles'] = smiles_standardization(compound_obj.isomeric_smiles)
    except Exception as e:
        print(f"Error processing row {index}: {e}")
        unmatched_smiles_rows.append(row)  # 记录错误行
    time.sleep(3)  # 避免频繁请求 PubChem

# 保存未匹配的 SMILES 数据
unmatched_smiles_df = pd.DataFrame(unmatched_smiles_rows)
unmatched_smiles_df.to_csv('unmatched_smiles_tcmbank.csv', index=False)

# 启用 RDKit 日志
RDLogger.EnableLog('rdApp.*')