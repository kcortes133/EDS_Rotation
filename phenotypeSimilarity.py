import requests
from itertools import chain
from io import StringIO
from tqdm import tqdm

import numpy as np

from scipy.sparse import coo_matrix
import seaborn as sns
import pandas as pd


def getWikiIDs(dataframe):

    api_url = "http://garbanzo.sulab.org/"
    c = ' '.join([x.upper() for x in dataframe.gene])
    endpoint = "exactmatches/"
    params = {'c': c}
    print('url',api_url+endpoint)
    r = requests.get(api_url + endpoint, params=params)
    #qids = [x for x in r.json() if "wd" in x]

    return

def getPhenoTerms(url, genes):
    url = "https://api.monarchinitiative.org/api/mart/gene/phenotype/NCBITaxon:9606"
    d = requests.get(url).json()
    phenos = {x['subject']: set(x['objects']) for x in d}
    return phenos

def getPhenotypes(gene_id):
    """Query Monarch to determine the phenotypes associated with a NCBI gene id."""
    url = "https://api.monarchinitiative.org/api/bioentity/gene/{}/phenotypes/".format(gene_id)
    res = requests.get(url)
    return set(res.json()["objects"])
