import requests

def get_solr_results(solr, params):
    resultCount = params['rows']
    while params['start'] < resultCount:
        solr_request = requests.get(solr, params=params)
        response = solr_request.json()
        resultCount = response['response']['numFound']
        params['start'] += params['rows']
        for doc in response['response']['docs']:
            yield doc
