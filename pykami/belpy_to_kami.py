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
            # FIXME This should be canonicalized in some way
            kr_name = '%s%s' % (residue_name[bp_stmt.mod], bp_stmt.mod_pos)
            kr = sub_agent.get_create_key_residue(kr_name)
            # Add edge between enzyme and substrate site
            phos = Phosphorylation(enz_agent, kr)
            # Add the flag to the site
            flag = kr.get_create_flag(flag_name)
            if flag.formula is None:
                flag.formula = str(phos.id)
            else:
                flag.formula += ' or %s' % phos.id
            nodes += [kr, flag, phos]
        else:
            # Add edge between enzyme and substrate agent
            phos = Phosphorylation(enz_agent, sub_agent)
            # Add the flag to the substrate agent, dependent on the specific
            # phosphorylation relationship
            flag = sub_agent.get_create_flag(flag_name)
            if flag.formula is None:
                flag.formula = str(phos.id)
            else:
                flag.formula += ' or %s' % phos.id
            nodes += [flag, phos]
        # Return the 3-4 nodes we've obtained/created
        return nodes

    def activity_modification(self, bp_stmt):
        # Get the agent
        agent = self.get_create_agent(bp_stmt.monomer_name)
        nodes = [agent]
        # Add an "activity" attribute to the agent
        activity_name = '%s_active' % bp_stmt.activity
        active_attr = agent.get_create_attribute(activity_name)
        # FIXME This should be canonicalized in some way
        kr_name = '.'.join(['%s%s' %(residue_name[m],p) for m,p in zip(bp_stmt.mod,bp_stmt.mod_pos)])
        flag_name = '.'.join(flag_names[m] for m in bp_stmt.mod)
        # Get the right statement for the agent/site condition
        
        if bp_stmt.mod_pos is None:
            condition = '%s.%s' % (agent.name, flag_name)
        else:
            condition = '%s.%s.%s' % (agent.name, kr_name, flag_name)
        # Get the right formula for increases vs. decreases
        if bp_stmt.relationship == 'directlyDecreases':
            qualifier = 'not '
        else:
            qualifier = ''
        # Build up the formula
        if active_attr.formula is None:
            active_attr.formula = '%s%s' % (qualifier, condition)
        else:
            active_attr.formula += ',\\n%s%s' % (qualifier, condition)
        nodes += [active_attr]
        return nodes

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        bps = pickle.load(f)

    bkc = BelpyKamiConverter()

    nodes = set([])
    for stmt in bps:
        if isinstance(stmt, belpy.statements.Phosphorylation):
            new_nodes = bkc.phosphorylation(stmt)
            nodes.update(new_nodes)
        elif isinstance(stmt, belpy.statements.ActivityModification):
            new_nodes = bkc.activity_modification(stmt)
            nodes.update(new_nodes)

    kg = Graph('example_graph', nodes=nodes)
    kg.render()
    kg.write()
