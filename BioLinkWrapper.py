import requests


class BioLinkWrapper(object):
    def __init__(self):
        self.endpoint = 'https://api.monarchinitiative.org/api/'
        self.params = {
            'fetch_objects': 'true',
            'rows': '100'
        }

    def get_gene(self, gene_curie):
        params = {}
        url = '{}bioentity/gene/{}'.format(self.endpoint, gene_curie)
        response = requests.get(url, params)
        return response.json()

    def get_orthologs(self, gene_curie, orth_taxon=None):
        params = {}
        if orth_taxon:
            params['homolog_taxon'] = orth_taxon
        url = '{}bioentity/gene/{}/homologs/'.format(self.endpoint, gene_curie)
        response = requests.get(url, params)
        return response.json()

    def get_phenotypes(self, gene_curie):
        url = '{}bioentity/gene/{}/phenotypes/'.format(self.endpoint, gene_curie)
        response = requests.get(url)
        return response.json()

    def get_diseases(self, gene_curie):
        url = '{}bioentity/gene/{}/diseases/'.format(self.endpoint, gene_curie)
        response = requests.get(url)
        return response.json()

    def get_interactions(self, gene_curie):
        url = '{}bioentity/gene/{}/interactions/'.format(self.endpoint, gene_curie)
        response = requests.get(url)
        return response.json()

    def get_functions(self, gene_curie):
        url = '{}bioentity/gene/{}/function/'.format(self.endpoint, gene_curie)
        response = requests.get(url)
        return response.json()

    def get_disease_models(self, disease_curie):
        url = '{}/bioentity/disease/{}/models/'.format(self.endpoint, disease_curie)
        response = requests.get(url)
        return response.json()

    def get_all_phenotypes_for_taxon(self, taxon_curie):
        # get phenotypes associated with taxid
        url = "mart/gene/phenotype/{}".format(self.endpoint, taxon_curie)
        response = requests.get(url)
        return response.json()

    def get_gene_function(self, gene_curie):
        # get function associated with gene
        url = "{}bioentity/gene/{}/function/".format(self.endpoint, gene_curie)
        response = requests.get(url, params=self.params)
        return response.json()