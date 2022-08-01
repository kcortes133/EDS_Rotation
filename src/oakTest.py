# Author: Katherina Cortes
# Date: July 27, 2022
# Purpose: See if I can load Monarch grah into OAK

# TODO: cant get graph to load
import oaklib
from oaklib.interfaces.obograph_interface import OboGraphInterface
from oaklib.implementations.pronto.pronto_implementation import ProntoImplementation
from oaklib.resource import OntologyResource
from oaklib.utilities.obograph_utils import graph_to_image, default_stylemap_path

g = OboGraphInterface.load_graph(graph='embeddings/katherina/graphs/monarchinitiative/Monarch/latest/monarch-kg.tar.gz', replace=False)