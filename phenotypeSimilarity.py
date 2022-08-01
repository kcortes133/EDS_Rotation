import pandas as pd
import apiRequests


def getGenes2Pheno(g2pfile):
    df = pd.read_csv(g2pfile, sep='\t')
    return {k: g['HPO-Term-ID'].tolist() for k,g in df.groupby('entrez-gene-symbol')}


# TODO: can get ancestors as well
# additively build
def jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection

    return float(intersection)/union


def simSearchJ(genes, g2p):
    similarity = {}

    for gene in genes:
        similarity[gene] = {}
        for g in g2p:
            simScore = jaccard(genes[gene], g2p[g])
            similarity[gene][g] = simScore

    return similarity



def simSearch(genes, simScoreURL):
    simMatches = {}

    for g in genes:
        phenos = genes[g]
        # remove EFO - experimental factor ontology ids
        phenos = [x for x in phenos if 'EFO:' not in x]
        m = apiRequests.getSimSearch(simScoreURL + '&id='.join(phenos))
        simMatches[g] = m
    return simMatches


def writeGenestoFile(genes, file):
    with open(file, 'w') as f:
        for g in genes:
            f.write(g+'\n')
    return


def getGeneSyms(geneMappings, revOrthologs, geneIDS):
    syms = []
    for gene in geneIDS:
        if 'HGNC' in gene:
            syms.append(geneMappings[gene])
        elif 'MGI' in gene:
            if gene in revOrthologs:
                syms.extend(list(revOrthologs[gene]))
        else: syms.append(gene)
    print(syms)
    return syms
