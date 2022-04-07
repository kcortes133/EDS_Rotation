import requests
import pandas as pd
import copy, json
from query import *
import graphviz


solr_url = 'https://solr.monarchinitiative.org/solr/golr/select'
url = 'https://api.monarchinitiative.org/api/bioentity/gene/NCBIGene:1289/interactions'
params ={
    'wt': 'json',
    'rows': 1
}
request = requests.get(url,params=params)
json_obj = json.loads(request.text)
json_formatt = json.dumps(json_obj, indent=2)

print(json_formatt)
with open('response_data.json', 'w') as outF:
    outF.write(json_formatt)


#core_set = 'https://raw.githubusercontent.com/ncats-tangerine/cq-notebooks/master/fa_gene_sets/fa_1_core_complex.txt'
#core_set = './genes.txt'

#columns = ['gene', 'interactor_id', 'interactor_symbol', 'qualifier', 'inferred_gene']
#dataframe = pd.read_csv(core_set, sep='\t', names=['gene', 'symbol'])
#interaction_params = {
#    'wt': 'json',
#    'rows': 100,
#    'start': 0,
#    'q': '*:*',
#    'fl': 'subject, subject_label, subject_closure, \
#           object, object_label, object_taxon',
#    'fq': ['relation_closure: "RO:0002434"']
#}

#interactTable = pd.DataFrame(columns=columns)

#for index, row in dataframe.iterrows():
#    params = copy.deepcopy(interaction_params)
#    params['fq'].append('subject_closure: "{0}" \
#                        OR subject_ortholog_closure: "{0}"'
#                        .format(row['gene']))
#    for doc in get_solr_results(solr_url, params):
#        result = {}
#        result['gene'] = row['symbol']
#        result['interactor_id'] = doc['object']
#        result['interactor_symbol'] = doc['object_label']
#        if row['gene'] in doc['subject_closure']:
#            result['qualifier'] = "direct"
#        else:
#            result['qualifier'] = "homology"
#        interact_table = interact_table.append(result, ignore_index=True)

#interact_table.head(10)