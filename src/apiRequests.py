# Author: Katherina Cortes
# Date: July 26, 2022
# Purpose: api requests

import requests


def getOrthologs(url):
    d = requests.get(url).json()
    orthologs = {x['subject']: set(x['objects']) for x in d}
    return orthologs


def getGeneOrthologs(allOrthologs, genes):
    geneOrthologs = []
    for g in genes:
        geneOrthologs.extend(list(allOrthologs[g]))
    return geneOrthologs


def getPhenotypes(url, gene):
    url = url.format(gene)
    d = requests.get(url).json()
    phenotypes = []
    for i in d['associations']:
        phenotypes.append(i['object']['id'])
    return phenotypes



#def getPhenotypes(gene_id):
#    """Query Monarch to determine the phenotypes associated with a NCBI gene id."""
#    url = "https://api.monarchinitiative.org/api/bioentity/gene/{}/phenotypes/".format(gene_id)
#    res = requests.get(url)
#    return set(res.json()["objects"])

def getDiseasePhenotypes(url):
    res = requests.get(url).json()
    phenotypes = []
    for i in res['associations']:
        phenotypes.append(i['object']['id'])
    return phenotypes

# sim/search
# ['matches']['id']
def getSimSearch(url):
    res = requests.get(url).json()

    matches = []
    for i in res['matches']:
        matches.append(i['id'])
    return matches
