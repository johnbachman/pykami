import pickle
import sys
import belpy.statements
from core import *

flag_names = {
        'Phosphorylation': 'phos',
        'PhosphorylationTyrosine': 'Yphos',
        'PhosphorylationSerine': 'Sphos',
        'PhosphorylationThreonine': 'Tphos',
        }

residue_name = {
        'Phosphorylation': '',
        'PhosphorylationTyrosine': 'Y',
        'PhosphorylationSerine': 'S',
        'PhosphorylationThreonine': 'T',
        }

class BelpyKamiConverter(object):

    def __init__(self):
        self.nodes_dict = {}

    def get_create_agent(self, name):
        if name in self.nodes_dict:
            return self.nodes_dict[name]
        else:
            agent = Agent(name)
            self.nodes_dict[name] = agent
            return agent

    def phosphorylation(self, bp_stmt):
        # Get the names of enzyme and substrate
        enz_agent = self.get_create_agent(bp_stmt.enz_name)
        sub_agent = self.get_create_agent(bp_stmt.sub_name)
        nodes = [enz_agent, sub_agent]
        # Get the flag name from the modification
        flag_name = flag_names[bp_stmt.mod]
        # Do we know the site? If so, create a site for the agent
        if bp_stmt.mod_pos:
            # Does this agent already have a site for this residue?
            kr_name = '%s%s' % (residue_name[bp_stmt.mod], bp_stmt.mod_pos)
            kr = sub_agent.get_create_key_residue(kr_name)
            # Add the flag to the site
            flag = kr.get_create_flag(flag_name)
            nodes += [kr, flag]
            # Add edge between enzyme and substrate site
            phos = Phosphorylation(enz_agent, kr)
            relationships = [phos]
            if flag.formula is None:
                flag.formula = str(phos.id)
            else:
                flag.formula += ' or %s' % phos.id
        else:
            # Add edge between enzyme and substrate agent
            phos = Phosphorylation(enz_agent, sub_agent)
            relationships = [phos]
            # Add the flag to the substrate agent, dependent on the specific
            # phosphorylation relationship
            flag = sub_agent.get_create_flag(flag_name)
            if flag.formula is None:
                flag.formula = str(phos.id)
            else:
                flag.formula += ' or %s' % phos.id

            nodes += [flag]
        # Return the 3-4 nodes we've obtained/created
        return (nodes, relationships)

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        bps = pickle.load(f)

    bkc = BelpyKamiConverter()

    nodes = set([])
    relationships = set([])
    for stmt in bps:
        if isinstance(stmt, belpy.statements.Phosphorylation):
            (new_nodes, new_relationships) = bkc.phosphorylation(stmt)
            nodes.update(new_nodes)
            relationships.update(new_relationships)

    kg = Graph('example_graph', agents=nodes, relationships=relationships)
    kg.render()
    kg.write()
