import tensorflow
from ensmallen.datasets.monarchinitiative import Monarch
from embiggen.pipelines import compute_node_embedding

g = Monarch()
node_embedding_method_name = 'SkipGram'

embedding, _ = compute_node_embedding(
    g,
    node_embedding_method_name=node_embedding_method_name
)
