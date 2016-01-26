import pickle
import sys
import indra.statements
from core import *
import json

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

kami_types = {
        'Agent': 'agent',
        'Site': 'region',
        'KeyResidue': 'key_r',
        'Attribute': 'attr',
        'Flag': 'flag'
        }

class IndraKamiConverter(object):

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
        enz_agent = self.get_create_agent(bp_stmt.enz.name)
        sub_agent = self.get_create_agent(bp_stmt.sub.name)
        nodes = [enz_agent, sub_agent]
        # Get the flag name from the modification
        flag_name = flag_names[bp_stmt.mod]
        # Do we know the site? If so, create a site for the agent
        if bp_stmt.mod_pos:
            # Does this agent already have a site for this residue?
            # FIXME This should be canonicalized in some way
            site_name = '%s%s' % (residue_name[bp_stmt.mod], bp_stmt.mod_pos)
            site = sub_agent.get_create_site(site_name)
            # Add the flag to the site
            flag = site.get_create_flag(flag_name)
            #if flag.formula is None:
            #    flag.formula = str(phos.id)
            #else:
            #    flag.formula += ' or %s' % phos.id

            # Add edge between enzyme and substrate flag
            phos = Phosphorylation(enz_agent, flag)
            nodes += [site, flag, phos]
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
        agent = self.get_create_agent(bp_stmt.monomer.name)
        nodes = [agent]
        # Add an "activity" attribute to the agent
        activity_name = '%s_active' % bp_stmt.activity
        active_attr = agent.get_create_attribute(activity_name)

        # FIXME This should be canonicalized in some way
        site_name = '.'.join(['%s%s' %(residue_name[m],p) \
                    for m,p in zip(bp_stmt.mod,bp_stmt.mod_pos)])
        flag_name = '.'.join(flag_names[m] for m in bp_stmt.mod)
        #site_name = '%s%s' % (residue_name[bp_stmt.mod], bp_stmt.mod_pos)
        #''.join(['%s%s' %(residue_name[m],p) for m,p in zip(bp_stmt.mod,bp_stmt.mod_pos)])
        #flag_name = flag_names[bp_stmt.mod]

        # Get the right statement for the agent/site condition
        if bp_stmt.mod_pos is None:
            condition = '%s.%s' % (agent.name, flag_name)
        else:
            condition = '%s.%s.%s' % (agent.name, site_name, flag_name)
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

def get_path(source):
    sub_source = source
    source_path = []
    while True:
        source_path.append(sub_source.id)
        if sub_source.parent:
            sub_source = sub_source.parent
        else:
            break
    source_path.reverse()
    return source_path

def nodes_to_kami(nodes):
    output = {}

    # Compile list of agents
    agent_list = []
    site_list = []
    key_rs_list = []
    attributes_list = []
    flags_list = []
    actions_list = []
    actions_binder_list = []
    edges_list = []

    def add_flags_attributes(node, path):
        node_class = kami_types[node.__class__.__name__]
        for flag_name, flag_obj in node.flags.iteritems():
            flags_list.append({'class': ['node', 'flag'],
                               'name': flag_obj.id,
                               'label': flag_name,
                               'dest_class': ['node', node_class],
                               'dest_path': path,
                               'values': [], # FIXME
                               'v_equiv': None})
        for attr_name, attr_obj in node.attributes.iteritems():
            attributes_list.append({'class': ['node', 'attr'],
                               'name': attr_obj.id,
                               'label': attr_name,
                               'dest_class': ['node', node_class],
                               'dest_path': path,
                               'values': [], # FIXME
                               'v_equiv': None})


    for node in nodes:
        # AGENTS
        if isinstance(node, Agent):
            print "agent", node.name
            agent_list.append({'class':['node', 'agent'],
                               'name':node.id,
                               'label': node.name,
                               'cx':None,
                               'cy':None,
                               'family':None,
                               'abstract':False})
            add_flags_attributes(node, get_path(node))
        # SITES
        if isinstance(node, Site):
            site_list.append({'class': ['node', 'region'],
                              'name': node.id,
                              'label': node.name,
                              'ag_name': node.parent.id,
                              'color': None})
            add_flags_attributes(node, get_path(node))
        # KEY RESIDUES
        if isinstance(node, KeyResidue):
            # If this node belongs to a region, 
            if isinstance(node.parent, Site):
                assert isinstance(node.parent.parent, Agent)
                region_name = node.parent.id
                ag_name = node.parent.parent.id
            # Otherwise, we should belong to an Agent
            else:
                assert isinstance(node.parent, Agent)
                ag_name = node.parent.id
                region_name = None

            key_rs_list.append({'class': ['node', 'key_r'],
                                'name': node.id,
                                'label': node.name,
                                'ag_name': ag_name,
                                'region_name': region_name,
                                'angle': None})
            add_flags_attributes(node, get_path(node))

            # Get flags/attributes for agent
            # Iterate over sites
            # For each site, get flags/attributes
            # For each site, get key residues
            # For each keyr, get flags/attributes
            # Iterate over keyr
            # For each keyr, get flags/attributes
        # ACTIONS
        if isinstance(node, Phosphorylation):
            print "phosphorylation", node.source, node.target
            phos_action_binder = [{'class': ['node', 'binder'],
                                   'name': 'left',
                                   'act_name': node.id},
                                  {'class': ['node', 'binder'],
                                   'name': 'right',
                                   'act_name': node.id}]
            actions_binder_list += phos_action_binder

            # EDGES LIST
            edges_list.append({
                        'class': ['edge'],
                        'in_class': ['node', 'binder'],
                        'in_path': [node.id, 'right'],
                        'out_class': ['node',
                                      kami_types[node.target.__class__.__name__]],
                        'out_path': get_path(node.target)})

            # Treat the enzyme and substrate the same
            context_list = []
            for source in [node.source, node.target]:

                # Build context for this node and all of its parents
                while True:
                    source_type = kami_types[source.__class__.__name__]
                    # Get the path to this node
                    source_path = get_path(source)
                    context = {'el_cl': ['node', source_type],
                               'el_path': source_path}

                    if source_type == 'flag' or \
                       source_type == 'attr':
                        context['el_value'] = ['unphos']
                    context_list.append(context)

                    # Go up the hierarchy for the next node
                    if source.parent:
                        source = source.parent
                    else:
                        break

            actions_list.append({'class': ['node', 'action', 'mod'],
                                'name': node.id,
                                'label': node.id,
                                'context': context_list})

    output['infos'] = [{'scale':1, 'center':'A'}]
    output['agents'] = agent_list
    output['regions'] = site_list
    output['key_rs'] = []
    output['attributes'] = attributes_list
    output['flags'] = flags_list
    output['actions'] = actions_list
    output['actions_binder'] = actions_binder_list
    output['edges'] = edges_list

    json_str = json.dumps(output, indent=2)
    with open('example.json', 'w') as f:
        f.write(json_str)

    return output

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        bps = pickle.load(f)

    bkc = IndraKamiConverter()

    nodes = set([])
    for stmt in bps:
        if isinstance(stmt, indra.statements.Phosphorylation):
            new_nodes = bkc.phosphorylation(stmt)
            nodes.update(new_nodes)
        elif isinstance(stmt, indra.statements.ActivityModification):
            new_nodes = bkc.activity_modification(stmt)
            nodes.update(new_nodes)

    output = nodes_to_kami(nodes)

    #kg = Graph('example_graph', nodes=nodes)
    #kg.render()
    #kg.write()
