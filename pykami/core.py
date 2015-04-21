from pygraphviz import *
edge_style = {'fontname': 'arial', 'fontsize': 9}

id_counter = 0

def get_id():
    global id_counter
    id_counter += 1
    return id_counter

class Graph(object):
    def __init__(self, name, agents, relationships):
        self.name = name
        self.agents = agents
        self.relationships = relationships
        self.g = AGraph(name=name, directed=True)

    def render(self):
        # Iterate over all of the agents...
        for agent in self.agents:
            if isinstance(agent, Agent):
                agent.render(self.g)
        # For all of the relationships...
        for rel in self.relationships:
            # Render the relationship
            rel.render(self.g)

    def write(self):
        self.g.write('%s.dot' % self.name)

class Component(object):
    def __init__(self, name=None, is_abstract=False, flags=None, attributes=None,
                 annotations=None):
        self.id = get_id()
        self.name = name
        self.is_abstract = is_abstract
        # Flags
        self.flags = {}
        if flags is None:
            flags = []
        for flag in flags:
            self.add_flag(flag)
        # Attributes
        self.attributes = {}
        if attributes is None:
            attributes = []
        for attribute in attributes:
            self.add_attribute(attribute)
        # Any additional information
        self.annotations = annotations

    def add_flag(self, flag):
        self.flags[flag.name] = flag

    def get_create_flag(self, flag_name):
        if flag_name in self.flags:
            return self.flags[flag_name]
        else:
            flag = Flag(flag_name, None)
            self.flags[flag_name] = flag
            return flag

    def add_attribute(self, attribute):
        self.attributes[attribute.name] = attribute

    def get_create_attribute(self, attribute_name):
        if attribute_name in self.attributes:
            return self.attributes[attribute_name]
        else:
            attribute = Attribute(attribute_name, None)
            self.attributes[attribute_name] = attribute
            return attribute

    def __str__(self):
        return self.name

class Agent(Component):
    def __init__(self, name, sites=None, key_residues=None, is_abstract=False,
                 flags=None, attributes=None, annotations=None):
        super(Agent, self).__init__(name, is_abstract=is_abstract, flags=flags,
                                    attributes=attributes,
                                    annotations=annotations)
        # Sites
        self.sites = {}
        if sites is None:
            sites = []
        for site in sites:
            self.add_site(site)
        # Key residues
        self.key_residues = {}
        if key_residues is None:
            key_residues = []
        for kr in self.key_residues:
            self.add_key_residue(kr)
        # Annotations
        if annotations is None:
            annotations = []
        self.annotations = annotations

    def add_site(self, site):
        self.sites[site.name] = site

    def get_create_site(self, site_name):
        if site_name in self.sites:
            return self.sites[site_name]
        else:
            site = Site(site_name)
            self.sites[site_name] = site
            return site

    def add_key_residue(self, kr):
        self.key_residues[kr.name] = kr

    def get_create_key_residue(self, kr_name):
        if kr_name in self.key_residues:
            return self.key_residues[kr_name]
        else:
            kr = KeyResidue(kr_name)
            self.key_residues[kr_name] = kr
            return kr

    def render(self, g):
        # Style parameters for agent nodes
        agent_style = {'color': 'lightgrey', 'style': 'filled',
                       'fontname': 'arial'}
        # Add the agent to the graph
        g.add_node(self.id, label=self.name, **agent_style)
        # Keep a list of all of of the sites and key residues contained
        # within this agent, so that we can create a subgraph around the
        # agent when we're done
        agent_nodes = [self.id]
        # Iterate over all of the sites directly attached to the agent...
        for site_name, site in self.sites.iteritems():
            site_nodes = site.render(g)
            # Add the added nodes to the nodes for the subgraph
            agent_nodes += site_nodes
            g.add_edge(self.id, site.id, label='site', **edge_style)
        # Iterate over key residues associated directly with the agent and
        # add them
        for kr_name, kr in self.key_residues.iteritems():
            kr_nodes = kr.render(g)
            agent_nodes += kr_nodes
            g.add_edge(self.id, kr.id, label='kr', **edge_style)
        # Iterate over any flags associated with this agent
        for flag_name, flag in self.flags.iteritems():
            flag_nodes = flag.render(g)
            agent_nodes += flag_nodes
            g.add_edge(self.id, flag.id, label='flag', **edge_style)
        # Create a subgraph for the agent, its sites and key residues
        g.add_subgraph(agent_nodes, 'cluster_%s' % self.id)
        return agent_nodes

class Site(Component):
    def __init__(self, name, key_residues=None, is_abstract=False, flags=None,
                 attributes=None, annotations=None):
        super(Site, self).__init__(name, is_abstract=is_abstract, flags=flags,
                                   attributes=attributes, annotations=annotations)
        # Key residues
        self.key_residues = {}
        if key_residues is None:
            key_residues = []
        for kr in self.key_residues:
            self.add_key_residue(kr)

    def render(self, g):
        site_style = {'color': 'red', 'style': 'filled', 'fontname': 'arial'}
        # Add this site to the graph and the list of nodes
        g.add_node(self.id, label=self.name, **site_style)
        site_nodes = [self.id]
        # Iterate over any key residues associated with this site and add
        for kr_name, kr in self.key_residues.iteritems():
            kr.render(g)
            site_nodes.append(kr.id)
            g.add_edge(self.id, kr.id, 'kr', **edge_style)
        # Iterate over any flags associated with this agent
        for flag_name, flag in self.flags.iteritems():
            flag.render(g)
            site_nodes.append(flag.id)
            g.add_edge(self.id, flag.id, label='flag', **edge_style)
        return site_nodes

class KeyResidue(Component):
    def __init(self, name):
        self.name = name

    def render(self, g):
        kr_style = {'color': 'green', 'style': 'filled', 'fontname': 'arial'}
        g.add_node(self.id, label=self.name, **kr_style)
        kr_nodes = [self.id]
        # Iterate over any flags associated with this key residue
        for flag_name, flag in self.flags.iteritems():
            flag.render(g)
            kr_nodes.append(flag.id)
            g.add_edge(self.id, flag.id, label='flag', **edge_style)
        return kr_nodes

class Flag(object):
    def __init__(self, name, formula=None):
        self.id = get_id()
        self.name = name
        self.formula = formula
        self.style = {'color':'pink', 'style': 'filled', 'shape':'component',
                      'fontname': 'arial', 'fontsize': 10, 'size': 15}

    def render(self, g):
        g.add_node(self.id, label='%s: %s' % (self.name, self.formula),
                   **self.style)
        flag_nodes = [self.id]
        return flag_nodes

class Attribute(Flag):
    def __init__(self, name, formula):
        super(Attribute, self).__init__(name, formula)
        self.style = {'color':'green', 'style': 'filled', 'shape':'component',
                      'fontname': 'arial', 'fontsize': 10, 'size': 15}

# Relationships ===============================================================

class Relationship(object):
    def __init__(self):
        self.style = {'shape': 'square', 'fontname': 'arial',
                'color': 'lightblue', 'style': 'filled', 'fontsize': 10}
        self.edge_style = {'fontname': 'arial', 'fontsize': 9,
                           'style': 'dotted'}
        # Get an ID for this node
        self.id = get_id()

class DirectedBinary(Relationship):
    def __init__(self, source, target):
        super(DirectedBinary, self).__init__()
        self.source = source
        self.target = target

    def render(self, g):
        # Add this node to the graph
        #import ipdb; ipdb.set_trace()
        g.add_node(self.id, label='%s' % self.id, **self.style)
        # Add directed edges from source to this node, and from this node
        # to the target
        g.add_edge(self.source.id, self.id, **self.edge_style)
        g.add_edge(self.id, self.target.id, **self.edge_style)

class Phosphorylation(DirectedBinary):
    def __init__(self, source, target):
        super(Phosphorylation, self).__init__(source, target)
        self.short_name = 'p'
        self.id = '%s%s' % (self.short_name, self.id)

class UndirectedNAry(Relationship):
    def __init__(self):
        super(UndirectedNAry, self).__init__()

    def render(self, g):
        # Add this node to the graph
        g.add_node(self.id, label=self.id, **self.style)
        # Add directed edges into this bind node from its two participating
        # nodes
        for node in self.node_list:
            g.add_edge(node.id, self.id, **self.edge_style)

class Bind(UndirectedNAry):
    def __init__(self, node1, node2):
        super(Bind, self).__init__()
        self.node_list = [node1, node2]
        self.short_name = 'b'
        self.id = '%s%d' % (self.short_name, self.id)

class Complex(UndirectedNAry):
    def __init__(self, node_list):
        super(Complex, self).__init__()
        self.node_list = node_list
        self.short_name = 'cplx'
        self.id = '%s%d' % (self.short_name, self.id)

if __name__ == '__main__':
    agents = []
    rel = []

    ITGA3 = Agent('ITGA3')
    agents.append(ITGA3)

    ITGAV = Agent('ITGAV')
    agents.append(ITGAV)

    c = Complex([ITGA3, ITGAV])
    rel.append(c)

    SHC1 = Agent('SHC1')
    agents.append(SHC1)

    GRB2 = Agent('GRB2')
    agents.append(GRB2)

    SRC = Agent('SRC')
    agents.append(SRC)

    RAF1 = Agent('RAF1')
    agents.append(RAF1)
    RAF1_Y341 = KeyResidue('Y341')
    RAF1.add_key_residue(RAF1_Y341)
    RAF1_S338 = KeyResidue('S338')
    RAF1.add_key_residue(RAF1_S338)
    RAF1_active = Flag('active', 'example_formula')
    RAF1.add_flag(RAF1_active)

    # Src phosphorylates RAF1 on Y341
    p = Phosphorylation(SRC, RAF1_Y341)
    rel.append(p)

    p = Phosphorylation(RAF1, RAF1_S338)
    rel.append(p)

    b = Bind(SHC1, GRB2)
    rel.append(b)

    SOS1 = Agent('SOS1')
    agents.append(SOS1)
    b = Bind(SHC1, SOS1)
    rel.append(b)

    kg = Graph('example_graph', agents=agents, relationships=rel)
    kg.render()
    kg.write()
