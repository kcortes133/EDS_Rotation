# TODO: ontobio doesnt work

import pandas as pd
import requests
from pprint import pprint
#import ontobio
from ontobio.ontol_factory import OntologyFactory

from ontobio.io.gafparser import GafParser
from ontobio.assoc_factory import AssociationSetFactory

import BioLinkWrapper

def getOrthologs(genes):
    genes_orthologs = list()
    wrapper = BioLinkWrapper.BioLinkWrapper()
    for index, gene in genes.iterrows():
        orthologs = wrapper.get_orthologs(gene[0])
        for orth in orthologs['associations']:
            orth_dict = {
                'gene': gene[1],
                'symbol': gene[0],
                'ortholog_name': orth['object']['label'],
                'ortholog_curie': orth['object']['id'],
                'orth_tax_label': orth['object']['taxon']['label'],
                'orth_tax_id': orth['object']['taxon']['id']
            }
            genes_orthologs.append(orth_dict)

    return genes_orthologs



def query_mygene(gene_curie):
    gene_curie = gene_curie.replace('NCBIGene:', '')
    url = 'https://mygene.info/v3/query?q={}&fields=all'.format(gene_curie)
    hit = requests.get(url)
    hit = hit.json()
    uniprot = hit['hits'][0]['uniprot']
    if 'Swiss-Prot' in uniprot.keys():
        if isinstance(uniprot['Swiss-Prot'], list):
            return uniprot['Swiss-Prot']
        elif isinstance(uniprot['Swiss-Prot'], str):
            return [uniprot['Swiss-Prot']]
        else:
            return None


def make_assocs(group, parse=False):
    p = GafParser()
    afactory = AssociationSetFactory()
    ofactory = OntologyFactory()
    ont = ofactory.create('pato')

    url = "http://geneontology.org/gene-associations/gene_association.{}.gz".format(group)
    if group == 'human':
        url = "http://geneontology.org/gene-associations/goa_human.gaf.gz"
    assocs = p.parse(url)
    print(type(assocs))
    if parse:
        return assocs
    else:
        return afactory.create_from_assocs(assocs, ontology=ont)


def parse_gene_functions(curie):
    BLW1 = BioLinkWrapper()
    function_list =list()
    functions = BLW1.get_gene_function(gene_curie=curie)
    if 'associations' in functions.keys():
        for assoc in functions['associations']:
            function_list.append(assoc['object']['label'])
    function_set = set(function_list)
    return ", ".join(function_set)


def funcSimGenes(asoc_hsap, genes):
    human_sims = list()
    for index, row in genes.iterrows():
        fa_gene_curie = row[0]
        fa_gene_name = row[1]
        ukb = query_mygene(fa_gene_curie)
        for hgene in list(asoc_hsap.subject_label_map.keys()):
            amScore = asoc_hsap.jaccard_similarity('UniProtKB:{}'.format("".join(ukb)), hgene)
            if amScore > .7 and amScore < 1:
                human_sims.append({
                    'fa_gene_name': fa_gene_name,
                    'fa_gene_curie': fa_gene_curie,
                    'sim_gene_name': asoc_hsap.label(hgene),
                    'sim_hit_curie': hgene,
                    'sim_score': amScore,
                })
    human_simDF = pd.DataFrame(data=human_sims)
    return human_simDF