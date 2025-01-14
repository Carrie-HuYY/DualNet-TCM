import requests
import re
import os


def step1_visit_website(target_url):
    """
    访问网站链接并获取网页源码
    :param target_url:网址链接
    :return: 网页源代码的txt格式
    """
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36'}

    proxies = {"http": None, "https": None}
    res = requests.get(url, headers=headers, proxies=proxies)
    res.encoding = 'utf-8'
    print(res.text)
    return res.text


def step2_get_infor(page_content):
    """
    获取各类方剂及其对应信息
    :param page_content: 网页源码
    :return:
    """
    # 获取方剂名称信息
    formula_name_pattern = r'<a.*?>(.*?)target="_blank">(.*?)<\/a>'
    formula_name_matches = re.findall(formula_name_pattern, page_content)

    # 获取方剂剂型
    Dosage_Form_pattern = r'<a.*?>(.*?)style = "max-height: 200px;overflow: auto" > (.*?)<\/a>'
    Dosage_Form_matches = re.findall(formula_name_pattern, page_content)

    return formula_name_matches, Dosage_Form_matches


if __name__ == '__main__':
    base_url = 'http://www.tcmip.cn/ETCM/index.php/Home/Index/fj_details.html?id='
    for i in range(1):
        extra_url = str(i + 1)
        url = base_url + extra_url
        print(url)
        url_content = step1_visit_website(url)
        formula_name, Dosage_Form = step2_get_infor(url_content)
