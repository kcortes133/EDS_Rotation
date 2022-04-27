# Author: Katherina Cortes
# Date: March 24, 2022
# Purpose: rotation intro

# install jupyter ntbk
# https://github.com/NCATS-Tangerine/cq-notebooks/tree/master/Orange_Demonstrator_1_CQs/OrangeQ1.1_PPI_Network

import pandas as pd
import PPI, phenotypeSimilarity
import argparse

import functSim

parser = argparse.ArgumentParser(description='Ehlers Danlos Exploratory Work')
parser.add_argument('--ppi', metavar='Protien-Protien-Interactors', type=bool, default=False, help='get interactors')

args = parser.parse_args()

def main():
    core_set = 'https://raw.githubusercontent.com/kcortes133/EDS_Rotation/master/genes.txt'

    columns = ['gene', 'interactor_id', 'interactor_symbol', 'qualifier', 'inferred_gene']
    dataframe = pd.read_csv(core_set, sep="   ", names=['gene', 'symbol'])  # sep on 3 spaces

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



    # phenotype Similarity
    #url = "https://api.monarchinitiative.org/api/mart/gene/phenotype/NCBITaxon:9606"
    #phenos = phenotypeSimilarity.getPhenotypes(url)

    orthologs = functSim.getOrthologs(dataframe)
    orth_columns = ['gene_name', 'gene_curie', 'ortholog_name', 'ortholog_curie', 'orth_tax_label', 'orth_tax_id']
    orth_df = pd.DataFrame(data=orthologs, columns=orth_columns)
    print('{} Orthologs found'.format(len(orth_df.index)))

    asoc_hsap = functSim.make_assocs('human')
    human_simDF = functSim.funcSimGenes(asoc_hsap, dataframe)
    print(human_simDF)

main()


