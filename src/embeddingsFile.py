# Author: Katherina Cortes
# Date: July 16, 2022
# Purpose: get genes from embedded monarch graph using cosine sim

import pandas as pd
import tarfile
from scipy.spatial.distance import cosine


def extractTarGz():
    file = tarfile.open('katherina (1).tar.gz')
    file.extractall('./embeddings')
    file.close()


def loadGraph():
    embeddedGraph = pd.read_csv('./embeddings/katherina/Monarch_Embedding', sep=',')
    humanGenes = []

    rows = embeddedGraph.iloc[:,0]
    for r in range(len(rows)):
        if 'https://' in rows[r]:
            newR = rows[r].split('/')[-1]
            rows[r] = newR
        if 'HGNC' in rows[r]:
            humanGenes.append(rows[r])
    embeddedGraph.index = rows
    embeddedGraph = embeddedGraph.drop('Unnamed: 0',1)
    # problematic - graph has duplicates that do not match
    embeddedGraph = embeddedGraph[~embeddedGraph.index.duplicated(keep='first')]
    print(embeddedGraph)

    return embeddedGraph, humanGenes


def geneMappings():
    mappingsF = 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/quarterly/tsv/hgnc_complete_set_2022-04-01.txt'
    mappingDB = pd.read_csv(mappingsF, sep='\t')
    print(mappingDB)
    return mappingDB[['hgnc_id', 'symbol']]


def cosineSimGenes(graph, genes, all_human_genes):
    # ehlers danlos
    eds_curie = "MONDO:0020066"
    # ehlers danlos hypermobility type
    eds4_curie = "MONDO:0007523"
    cosine_sims_eds_vs_genes = {}
    for this_gene_curie in all_human_genes:
        cosine_sims_eds_vs_genes[this_gene_curie] = cosine(graph.loc[eds4_curie], graph.loc[this_gene_curie])
    return cosine_sims_eds_vs_genes


mappingsDB = geneMappings()
geneMappingsDict = dict(zip(mappingsDB['hgnc_id'], mappingsDB['symbol']))
graph, humanGenes = loadGraph()
cosine_sims_eds_vs_genes = cosineSimGenes(graph, [], humanGenes)

rankedGenes = sorted(cosine_sims_eds_vs_genes.items(), key=lambda item: item[1], reverse=True)
for gene in rankedGenes[-20:]:
    if 'HGNC' in gene[0]:
        print(gene)
        #print(geneMappingsDict[gene[0]])


