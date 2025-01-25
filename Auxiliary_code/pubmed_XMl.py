from lxml import etree

def parse_article_meta(tree):
    """
    Parse PMID, PMC and DOI from given article tree
    """
    article_meta = tree.find(".//article-meta")
    if article_meta is not None:
        pmid_node = article_meta.find('article-id[@pub-id-type="pmid"]')
        pmc_node = article_meta.find('article-id[@pub-id-type="pmc"]')
        pub_id_node = article_meta.find('article-id[@pub-id-type="publisher-id"]')
        doi_node = article_meta.find('article-id[@pub-id-type="doi"]')
    else:
        pmid_node = None
        pmc_node = None
        pub_id_node = None
        doi_node = None

    pmid = pmid_node.text if pmid_node is not None else ""
    pmc = pmc_node.text if pmc_node is not None else ""
    pub_id = pub_id_node.text if pub_id_node is not None else ""
    doi = doi_node.text if doi_node is not None else ""

    dict_article_meta = {"pmid": pmid, "pmc": pmc, "doi": doi, "publisher_id": pub_id}

    return dict_article_meta


def parse_date(tree, date_type):
    """Parse publication dates based on the provided date type."""

    def get_text(node):
        return node.text if node is not None else None

    pub_date_path = f".//pub-date[@pub-type='{date_type}' or @date-type='{date_type}']"
    date_node = tree.xpath(pub_date_path)

    if not date_node:
        return {}

    date_dict = {}
    for part in ["year", "month", "day"]:
        text = get_text(date_node[0].find(part))
        if text is not None:
            date_dict[part] = text

    return date_dict

def read_xml(path, nxml=False):
    """
    Parse tree from given XML path
    """
    try:
        tree = etree.parse(path)
        if ".nxml" in path or nxml:
            remove_namespace(tree)  # strip namespace when reading an XML file
    except:
        try:
            tree = etree.fromstring(path)
        except Exception:
            print(
                "Error: it was not able to read a path, a file-like object, or a string as an XML"
            )
            raise
    return tree


def parse_pubmed_xml(path, include_path=False, nxml=False):
    """
    Given an input XML path to PubMed XML file, extract information and metadata
    from a given XML file and return parsed XML file in dictionary format.
    You can check ``ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/`` to list of available files to download

    Parameters
    ----------
    path: str
        A path to a given PumMed XML file
    include_path: bool
        if True, include a key 'path_to_file' in an output dictionary
        default: False
    nxml: bool
        if True, this will strip a namespace of an XML after reading a file
        see https://stackoverflow.com/questions/18159221/remove-namespace-and-prefix-from-xml-in-python-using-lxml to
        default: False

    Return
    ------
    dict_out: dict
        A dictionary contains a following keys from a parsed XML path
        'full_title', 'abstract', 'journal', 'pmid', 'pmc', 'doi',
        'publisher_id', 'author_list', 'affiliation_list', 'publication_year',
        'publication_date', 'epublication_date', 'subjects'
    """
    tree = read_xml(path, nxml)

    tree_title = tree.find(".//title-group/article-title")
    if tree_title is not None:
        title = [t for t in tree_title.itertext()]
        sub_title = tree.xpath(".//title-group/subtitle/text()")
        title.extend(sub_title)
        title = [t.replace("\n", " ").replace("\t", " ") for t in title]
        full_title = " ".join(title)
    else:
        full_title = ""

    try:
        abstracts = list()
        abstract_tree = tree.findall(".//abstract")
        for a in abstract_tree:
            for t in a.itertext():
                text = t.replace("\n", " ").replace("\t", " ").strip()
                abstracts.append(text)
        abstract = " ".join(abstracts)
    except BaseException:
        abstract = ""

    journal_node = tree.findall(".//journal-title")
    if journal_node is not None:
        journal = " ".join(["".join(node.itertext()) for node in journal_node])
    else:
        journal = ""

    dict_article_meta = parse_article_meta(tree)

    pub_date_dict = parse_date(tree, "ppub")
    if "year" not in pub_date_dict:
        pub_date_dict = parse_date(tree, "collection")
    pub_date = format_date(pub_date_dict)

    try:
        pub_year = int(pub_date_dict.get("year"))
    except TypeError:
        pub_year = None

    epub_date = format_date(parse_date(tree, "epub"))

    subjects_node = tree.findall(".//article-categories//subj-group/subject")
    subjects = list()
    if subjects_node is not None:
        for s in subjects_node:
            subject = " ".join([s_.strip() for s_ in s.itertext()]).strip()
            subjects.append(subject)
        subjects = "; ".join(subjects)
    else:
        subjects = ""

    # create affiliation dictionary
    affil_id = tree.xpath(".//aff[@id]/@id")
    if len(affil_id) > 0:
        affil_id = list(map(str, affil_id))
    else:
        affil_id = [""]  # replace id with empty list

    affil_name = tree.xpath(".//aff[@id]")
    affil_name_list = list()
    for e in affil_name:
        name = stringify_affiliation_rec(e)
        name = name.strip().replace("\n", " ")
        affil_name_list.append(name)
    affiliation_list = [[idx, name] for idx, name in zip(affil_id, affil_name_list)]

    tree_author = tree.xpath('.//contrib-group/contrib[@contrib-type="author"]')
    author_list = list()
    for author in tree_author:
        author_aff = author.findall('xref[@ref-type="aff"]')
        try:
            ref_id_list = [str(a.attrib["rid"]) for a in author_aff]
        except BaseException:
            ref_id_list = ""
        try:
            author_list.append(
                [
                    author.find("name/surname").text,
                    author.find("name/given-names").text,
                    ref_id_list,
                ]
            )
        except BaseException:
            author_list.append(["", "", ref_id_list])
    author_list = flatten_zip_author(author_list)

    coi_statement = '\n'.join(parse_coi_statements(tree))

    dict_out = {
        "full_title": full_title.strip(),
        "abstract": abstract,
        "journal": journal,
        "pmid": dict_article_meta["pmid"],
        "pmc": dict_article_meta["pmc"],
        "doi": dict_article_meta["doi"],
        "publisher_id": dict_article_meta["publisher_id"],
        "author_list": author_list,
        "affiliation_list": affiliation_list,
        "publication_year": pub_year,
        "publication_date": pub_date,
        "epublication_date": epub_date,
        "subjects": subjects,
        "coi_statement": coi_statement,
    }
    if include_path:
        dict_out["path_to_file"] = path
    return dict_out


dict_out = parse_pubmed_xml("pubmed25n0001.xml/pubmed25n0001.xml")
print(dict_out)