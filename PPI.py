import pandas as pd
import copy
from query import *


# @param solr: url to query
# @param params: parameters for solr request
def get_solr_results(solr, params):
    resultCount = params['rows']
    while params['start'] < resultCount:
        solr_request = requests.get(solr, params=params)
        response = solr_request.json()
        resultCount = response['response']['numFound']
        params['start'] += params['rows']
        for doc in response['response']['docs']:
            yield doc
    return


def makeInteractorTable(columns, dataframe, solr_url, interaction_params):
    # Make new dataframe for results
    interact_table = pd.DataFrame(columns=columns)

    # Get interactions, both direct and inferred
    for index, row in dataframe.iterrows():
        params = copy.deepcopy(interaction_params)
        params['fq'].append('subject_closure: "{0}" \
                            OR subject_ortholog_closure: "{0}"'
                            .format(row['gene']))
        for doc in get_solr_results(solr_url, params):
            result = {}
            result['gene'] = row['symbol']
            result['interactor_id'] = doc['object']
            result['interactor_symbol'] = doc['object_label']
            if row['gene'] in doc['subject_closure']:
                result['qualifier'] = "direct"
            else:
                result['qualifier'] = "homology"
            interact_table = interact_table.append(result, ignore_index=True)
    return interact_table



# @param interactorsTable: table of genes and interactors
def getTopInteractors(interact_table):
    interactors = interact_table[['gene', 'interactor_symbol']].drop_duplicates()
    all_interactors = interactors.value_counts('interactor_symbol').reset_index()
    multi_interactors = all_interactors[all_interactors[0] > 1]

    print(multi_interactors)



def getPairWiseInteractors(geneList, interact_table, minInteractions):
    interactors = {}
    topInteractors = {}
    for gene in geneList:
        geneInteract = list(interact_table[interact_table['gene'].isin([gene])]['interactor_symbol'])
        for interactor in geneInteract:
            if interactor in interactors:
                interactors[interactor].append(gene)
            else:
                interactors[interactor] = [gene]


    for i in interactors:
        if len(interactor[i]) >= minInteractions:
            topInteractors[i] = copy(interactors[i])

    return topInteractors
