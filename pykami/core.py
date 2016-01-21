from pygraphviz import *

edge_style = {'fontname': 'arial', 'fontsize': 9}

id_counter = 0

def get_id():
    """Get a unique identifier for each component as it is created."""
    global id_counter
    id_counter += 1
    return id_counter


class Graph(object):
    """A container for the nodes and relationships in the graph.

    Contains the top-level method for rendering the graph in GraphViz.

    Parameters
    ----------
    name : string
        The name of the graph.
    nodes : list of nodes
        The list should contain all relevant Agents and Relationships.  Though
        Sites and KeyResidues can be included in this list, they are ignored
        in the top-level rendering step because they are (currently)
        presumed to be contained with Agents.

    Attributes
    ----------
    g : pygraphviz.AGraph
        Instance of the PyGraphviz graph to which the nodes are rendered.
    """
    def __init__(self, name, nodes):
        self.name = name
        self.nodes = nodes
        self.g = AGraph(name=name, directed=True)

    def render(self):
        """Recursively render nodes, sub-nodes and edges.

        Builds up the pygraphviz graph data structure but does not write the
        graph to a file.

        The agents and their sub-nodes (sites, key residues) are rendered
        first, and the relationships (with edges connecting agents) are
        rendered afterwards, because the node labels and properties do not
        come out correctly if the nodes are referenced by an edge before
        they are explicitly created.
        """
        relationships = []
        # Iterate over all of the nodes, rendering only the agents (but
        # collecting the relationships for the next round of rendering)
        for node in self.nodes:
            if isinstance(node, Agent):
                node.render(self.g)
            elif isinstance(node, Relationship):
                relationships.append(node)
        # Iterate again, rendering the relationships this time
        for rel in relationships:
            rel.render(self.g)

    def write(self):
        """Write the graph to a file after rendering."""
        self.g.write('%s.dot' % self.name)


class Component(object):
    """Parent class for Agents, Sites, and KeyResidues.

    Agents, Sites and KeyResidues all shared functionality, in particular they
    can contain flags, attributes, and annotations.

    Attributes
    ----------
    flags : dict
        Dictionary mapping the name of a flag to the Flag object.
    attributes : dict
        Dictionary mapping the name of an attribute to the Attribute object.
    """
    def __init__(self, name=None, parent=None, is_abstract=False, flags=None,
                 attributes=None, annotations=None):
        # Get a globally unique, numeric identifier for the component
        self.id = get_id()
        self.name = name
        self.parent = parent
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
        if annotations is None:
            annotations = []
        self.annotations = annotations

    def add_flag(self, flag):
        """Add a flag to the flags dict."""
        self.flags[flag.name] = flag

    def get_create_flag(self, flag_name):
        """Return the flag with the given name if present, or create it."""
        if flag_name in self.flags:
            return self.flags[flag_name]
        else:
            flag = Flag(flag_name, self, None)
            self.flags[flag_name] = flag
            return flag

    def add_attribute(self, attribute):
        """Add an attribute to the attributes dict."""
        self.attributes[attribute.name] = attribute

    def get_create_attribute(self, attribute_name):
        """Return the attribute with the given name if present, or create it."""
        if attribute_name in self.attributes:
            return self.attributes[attribute_name]
        else:
            attribute = Attribute(attribute_name, self, None)
            self.attributes[attribute_name] = attribute
            return attribute

    def render(self, g):
        """Render flags and attributes for the component."""
        component_nodes = []
        # Iterate over any flags associated with this component
        for flag_name, flag in self.flags.iteritems():
            flag_nodes = flag.render(g)
            component_nodes += flag_nodes
            g.add_edge(self.id, flag.id, label='flag', **edge_style)
        # Iterate over any attributes associated with this component
        for attribute_name, attribute in self.attributes.iteritems():
            attribute_nodes = attribute.render(g)
            component_nodes += attribute_nodes
            g.add_edge(self.id, attribute.id, label='attr', **edge_style)
        return component_nodes

    def __str__(self):
        return ("%s(%s)" % (type(self).__name__, self.name))


class Agent(Component):
    """Top-level nodes containing sites and key residues (e.g., proteins)"""
    def __init__(self, name, sites=None, key_residues=None,
                 is_abstract=False, flags=None, attributes=None,
                 annotations=None):
        # Parent constructor
        super(Agent, self).__init__(name, parent=None, is_abstract=is_abstract,
                                    flags=flags, attributes=attributes,
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

    def add_site(self, site):
        """Add a site to the sites dict."""
        self.sites[site.name] = site

    def get_create_site(self, site_name):
        """Return the site with the given name if present, or create it."""
        if site_name in self.sites:
            return self.sites[site_name]
        else:
            site = Site(site_name, self)
            self.sites[site_name] = site
            return site

    def add_key_residue(self, kr):
        """Add a key residue ite to the key residues dict."""
        self.key_residues[kr.name] = kr

    def get_create_key_residue(self, kr_name):
        """Return the residue with the given name if present, or create it."""
        if kr_name in self.key_residues:
            return self.key_residues[kr_name]
        else:
            kr = KeyResidue(kr_name, self)
            self.key_residues[kr_name] = kr
            return kr

    def render(self, g):
        """Build the graph for the agent and its subnodes.

        Subnodes include sites, key residues, flags, and attributes.
        """
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
        # Call the Component render method to render flags and attributes
        flags_attrs = super(Agent, self).render(g)
        agent_nodes += flags_attrs
        # Create a subgraph for the agent, its sites and key residues
        g.add_subgraph(agent_nodes, 'cluster_%s' % self.id)
        return agent_nodes


class Site(Component):
    """Nodes representing logical or physical states of proteins.

    Sites can contain key residues but not agents or other sites.
    """
    def __init__(self, name, parent, key_residues=None, is_abstract=False,
                 flags=None, attributes=None, annotations=None):
        # Parent constructor
        super(Site, self).__init__(name, parent=parent, is_abstract=is_abstract,
                                   flags=flags, attributes=attributes,
                                   annotations=annotations)
        # Key residues
        self.key_residues = {}
        if key_residues is None:
            key_residues = []
        for kr in self.key_residues:
            self.add_key_residue(kr)

    def render(self, g):
        """Build the graph for the site and its subnodes.

        Subnodes include key residues, flags, and attributes.
        """
        site_style = {'color': 'red', 'style': 'filled', 'fontname': 'arial'}
        # Add this site to the graph and the list of nodes
        g.add_node(self.id, label=self.name, **site_style)
        site_nodes = [self.id]
        # Iterate over any key residues associated with this site and add
        for kr_name, kr in self.key_residues.iteritems():
            kr.render(g)
            site_nodes.append(kr.id)
            g.add_edge(self.id, kr.id, 'kr', **edge_style)
        # Call the Component render method to render flags and attributes
        flags_attrs = super(Site, self).render(g)
        site_nodes += flags_attrs
        return site_nodes

class KeyResidue(Component):
    """Amino acid residues defining agent/site functionality.

    Because they only have the functionality of the base Component class,
    the constructor only calls the parent class.
    """
    def __init__(self, name, parent, is_abstract=False,
                 flags=None, attributes=None, annotations=None):
        # Parent constructor: makes sure 'parent' field is filled out
        super(KeyResidue, self).__init__(name, parent=parent,
                                         is_abstract=is_abstract, flags=flags,
                                         attributes=attributes,
                                         annotations=annotations)

    def render(self, g):
        """Build the graph for the key residue and its subnodes.

        Subnodes include flags and attributes.
        """
        kr_style = {'color': 'green', 'style': 'filled', 'fontname': 'arial'}
        g.add_node(self.id, label=self.name, **kr_style)
        kr_nodes = [self.id]
        # Call the Component render method to render flags and attributes
        flags_attrs = super(KeyResidue, self).render(g)
        kr_nodes += flags_attrs
        return kr_nodes

class Flag(object):
    def __init__(self, name, parent, formula=None):
        self.id = get_id()
        #self.name = '%s%s' % (name, self.id)
        self.name = name
        self.parent = parent
        self.formula = formula
        self.style = {'color':'pink', 'style': 'filled', 'shape':'component',
                      'fontname': 'arial', 'fontsize': 10, 'size': 15}

    def render(self, g):
        g.add_node(self.id, label='%s: %s' % (self.name, self.formula),
                   **self.style)
        flag_nodes = [self.id]
        return flag_nodes

class Attribute(Flag):
    def __init__(self, name, parent, formula):
        super(Attribute, self).__init__(name, parent, formula)
        self.style = {'color':'sandybrown', 'style': 'filled',
                      'shape':'component', 'fontname': 'arial',
                      'fontsize': 10, 'size': 15}

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

if __name__ == '__main__':
    nodes = []

    ITGA3 = Agent('ITGA3')
    nodes.append(ITGA3)

    ITGAV = Agent('ITGAV')
    nodes.append(ITGAV)

    c = Bind(ITGA3, ITGAV)
    nodes.append(c)

    SHC1 = Agent('SHC1')
    nodes.append(SHC1)

    GRB2 = Agent('GRB2')
    nodes.append(GRB2)

    SRC = Agent('SRC')
    nodes.append(SRC)

    RAF1 = Agent('RAF1')
    nodes.append(RAF1)
    RAF1_Y341 = KeyResidue('Y341')
    RAF1.add_key_residue(RAF1_Y341)
    RAF1_S338 = KeyResidue('S338')
    RAF1.add_key_residue(RAF1_S338)
    RAF1_active = Flag('active', 'example_formula')
    RAF1.add_flag(RAF1_active)

    # Src phosphorylates RAF1 on Y341
    p = Phosphorylation(SRC, RAF1_Y341)
    nodes.append(p)

    p = Phosphorylation(RAF1, RAF1_S338)
    nodes.append(p)

    b = Bind(SHC1, GRB2)
    nodes.append(b)

    SOS1 = Agent('SOS1')
    nodes.append(SOS1)
    b = Bind(SHC1, SOS1)
    nodes.append(b)

    kg = Graph('example_graph', nodes)
    kg.render()
    kg.write()
