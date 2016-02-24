from indra.biopax import biopax_api
from indra.trips import trips_api
from indra_to_kami import nodes_to_kami, IndraKamiConverter
import indra.statements
import json

# Get a biopax processor from a biopax query
#bp = biopax_api.process_pc_pathsfromto(['BRAF'], ['MAP2K1'])
#bp.get_phosphorylation()

tp = trips_api.process_text('MEK2 phosphorylates ERK1 at Thr-202 and Tyr-204')

ikc = IndraKamiConverter()
nodes = set([])
# Collect the nodes to generate from the INDRA statements
for stmt in tp.statements:
    if isinstance(stmt, indra.statements.Phosphorylation):
        new_nodes = ikc.phosphorylation(stmt)
        nodes.update(new_nodes)

# Create the JSON output
output = nodes_to_kami(nodes)

json_str = json.dumps(output, indent=2)
with open('indra_to_kami_example1.json', 'w') as f:
    f.write(json_str)


