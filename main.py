# Author: Katherina Cortes
# Date: March 24, 2022
# Purpose: Gene

'''
Genes associated with EDS
1. EDS Core gene set (provided by Ellen)
2. All non-human orthologs of EDS core gene set (api.monarchinitiative.org - works)
    https://api.monarchinitiative.org/api/mart/ortholog/NCBITaxon%3A9606/NCBITaxon%3A10090
3. Get phenotypes for EDS and orthologs set (api.monarchinitiative.org - works (gene/phenotype call)
4. Do phenotype similarity to return ranked phenotypically similar genes (api.monarch returns 500 for sim search)
    we need all human g2p associations (api.monarch mart call times out)
    we can get these from HPO downloads from a gene2phenotype tab file with bad headers
5. Return candidates above some sim index threshold (python)
6. Get human orthologs of candidate subset. (api.monarchinitiative.org - already have these from the above call.)

We also want to do Phenotypes associated with EDS
'''

# install jupyter ntbk
# https://github.com/NCATS-Tangerine/cq-notebooks/tree/master/Orange_Demonstrator_1_CQs/OrangeQ1.1_PPI_Network

import pandas as pd
from src import apiRequests, PPI, phenotypeSimilarity, embeddingsFile
import argparse

parser = argparse.ArgumentParser(description='Ehlers Danlos Exploratory Work')
parser.add_argument('--ppi', metavar='Protien-Protien-Interactors', type=bool, default=False, help='get interactors')
parser.add_argument('--pheno', metavar='EDS Phenotype', type=bool, default=False, help='get genes by EDS phenotypes')
parser.add_argument('--ortho', metavar='EDS Mouse Orthologs', type=bool, default=False, help='get genes by orthologs to '
                                                                                             'Ellen genes')
args = parser.parse_args()

def main():
    core_set = './genes.txt'
    columns = ['gene', 'interactor_id', 'interactor_symbol', 'qualifier', 'inferred_gene']
    dataframe = pd.read_csv(core_set, sep="   ", names=['NCBI_ID', 'HGNC_ID', 'gene'])  # sep on 3 spaces

    g2pFile = './genes_to_phenotype.txt'
    g2p = phenotypeSimilarity.getGenes2Pheno(g2pFile)

    monarchURL = 'https://api.monarchinitiative.org/api/'
    orthologURL = monarchURL + 'mart/ortholog/NCBITaxon:9606/NCBITaxon:10090'
    revOrthologURL = monarchURL + 'mart/ortholog/NCBITaxon:10090/NCBITaxon:9606'
    phenotypeURL = monarchURL + 'bioentity/gene/{}/phenotypes'
    # EDS general
    #disPhenotypeURL = monarchURL + 'bioentity/disease/MONDO:0020066/phenotypes'
    # EDS hypermobility type
    disPhenotypeURL = monarchURL + 'bioentity/disease/MONDO:0007523/phenotypes'
    simScoreURL = 'https://api-gcp.monarchinitiative.org/api/sim/search?id='
    mappingsDB = embeddingsFile.geneMappings()
    geneMappingsDict = dict(zip(mappingsDB['hgnc_id'], mappingsDB['symbol']))
    revOrthologs = apiRequests.getOrthologs(revOrthologURL)

    # Protien Protien Interactions
    if args.ppi:
        solr_url = 'https://solr.monarchinitiative.org/solr/golr/select'

        interaction_params = {
            'wt': 'json',
            'rows': 100,
            'start': 0,
            'q': '*:*',
            'fl': 'subject, subject_label, subject_closure, \
                   object, object_label, object_taxon',
            'fq': ['relation_closure: "RO:0002434"']
        }


        interactor_table = PPI.makeInteractorTable(columns, dataframe, solr_url, interaction_params)

        topOverallInteractors = PPI.getTopInteractors(interactor_table)

    if args.pheno:
        # get phenotypes associated with EDS
        disPhenotypes = apiRequests.getDiseasePhenotypes(disPhenotypeURL)

        # jaccard sim between all genes and EDS genes
        #disSim = phenotypeSimilarity.simSearch({'EDS': disPhenotypes}, simScoreURL)
        disSim = phenotypeSimilarity.simSearchJ({'EDS': disPhenotypes}, g2p)
        rankedGenes = sorted(disSim['EDS'].items(), key=lambda item: item[1], reverse=True)
        print(rankedGenes[:20])
        #phenotypeSimilarity.writeGenestoFile(rankedGenes[:20], './outputs/EDSGenes_PhenoSim.txt')
        gPS = list(map(lambda x: x[0], rankedGenes[:20]))
        #gPS = phenotypeSimilarity.getGeneSyms(geneMappingsDict,revOrthologs, disSim['EDS'])


    if args.ortho:
        orthologs = apiRequests.getOrthologs(orthologURL)
        geneOrthologs = apiRequests.getGeneOrthologs(orthologs, dataframe['HGNC_ID'])
        allGenesAndOrthos = geneOrthologs + list(dataframe['HGNC_ID'])

        # Phenotypes of disease
        # //TODO: sim search top 20
        allPhenotypes = {}
        for gene in allGenesAndOrthos:
            phenotypes = apiRequests.getPhenotypes(phenotypeURL, gene)
            allPhenotypes[gene] = phenotypes



        for p in allPhenotypes:
            print(p)
        similarity = phenotypeSimilarity.simSearchJ(allPhenotypes, g2p)
        #similarity = phenotypeSimilarity.simSearch(allPhenotypes, simScoreURL)

        allGenes = []
        rankedGenes = {}
        for gene in similarity:
            ranked = sorted(similarity[gene].items(), key=lambda item: item[1], reverse=True)[:20]
            #phenotypeSimilarity.writeGenestoFile(ranked[:20], './outputs/EDSGenes_'+ str(gene).replace(':', '_')+'_OrthologSim.txt')
            #ranked = phenotypeSimilarity.getGeneSyms(geneMappingsDict, revOrthologs, similarity[gene])
            for r in ranked:
                if r[0] in rankedGenes:
                    rankedGenes[r[0]] += 1
                else: rankedGenes[r[0]] = 1
            allGenes.extend(ranked)
        #gOS = set(allGenes)
        gOS = set(map(lambda x: x[0], allGenes))
        #phenotypeSimilarity.writeGenestoFile(gOS, './outputs/EDSAllGenes_OrthologSim.txt')
        print(sorted(rankedGenes.items(), key=lambda item: item[1], reverse=True)[:20])

    if args.ortho and args.pheno:
        print(gPS)
        print(gOS.intersection(gPS))
        #phenotypeSimilarity.writeGenestoFile(gOS.intersection(gPS), './outputs/EDS3Genes_Ortho&PhenoSim.txt')


main()


