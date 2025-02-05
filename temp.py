from sklearn.feature_extraction.text import CountVectorizer, TfidfVectorizer
import spacy

nlp = spacy.load("en_core_web_sm")

documents = ('In March 1977, a large volume of the industrial chemical hexachlorocyclopentadiene (HCCPD) was dumped into '
             'a municipal sewage system in Kentucky. We evaluated the health effects of exposure to HCCPD in '
             '145 sewage treatment plant workers. We found that 85 (59%) had noted eye irritation, 65 (45%) had headaches, '
             'and 39 (27%) had throat irritation. Symptoms occurred throughout the plant; however, highest attack rates '
             'occurred in primary sewage treatment areas. Medical examination of 41 employees three days after the plant'
             ' was closed showed proteinuria and elevation of serum lactic dehydrogenase levels; these findings were not '
             'present three weeks later. This episode demonstrates the toxicity of HCCPD and emphasizes the vulnerability '
             'of sewage workers to chemical toxins in wastewater systems.')
doc = nlp(documents)

genes = []
diseases = []

for ent in doc.ents:
    if ent.label_ == "CHEMICAL":  # 假设病名被识别为地理位置（GPE）
        diseases.append(ent.text)
    elif ent.label_ == "ORG":  # 假设基因名被识别为组织（ORG）
        genes.append(ent.text)

print("Genes:", genes)
print("Diseases:", diseases)