#!/usr/bin/env python
# coding: utf-8

# In[5]:


import networkx as nx
import itertools
from random import choice
import pandas as pd
import matplotlib.pyplot as plt
from statistics import median
from sys import stderr
from mpl_toolkits.mplot3d import axes3d
from random import choice
from math import inf


# In[2]:


def vertex_ref(g, i):
    neighbors = {n: all_neighbours(g, n, i-1) for n in g.nodes()}
    degs = {}
    for k,v in neighbors.items():
        if i == 0:
            degs[k] = [0]
        else:
            degs[k] = sorted([g.degree(n) for n in v])

    print('  Vertex Ref[{} i] {:.3%}\t\r'.format(i, (list(neighbors.keys()).index(k) + 1)/len(neighbors.items())), file=stderr, end='')
    print()
    return degs


def all_neighbours(g, start, k):
    nbrs = set([start])
    for i in range(k):
        nbrs = set((nbr for n in nbrs for nbr in g[n]))
    return nbrs

def edge_facts_subgraph(g, g_pert, edge_facts):
    degs = {}
    nodes_number = len(g.nodes) * len(g_pert.nodes)

    i = 0
    for start in g.nodes:
        fact = all_subgraphs(g, edge_facts, start)
        degs[start] = []
        for node in g_pert.nodes:
            subgraph = all_subgraphs(g_pert, edge_facts, node)
            if nx.is_isomorphic(fact, subgraph):
                    degs[start].append(node)          
            i += 1
    print('Edge Facts[{} edges] {:.3%}\t\r'.format(edge_facts, i/nodes_number), file=stderr, end='')
    
    print()
    return degs


def all_subgraphs(g, edge_facts, start):
    graph = nx.Graph()
    graph.add_node(start)

    for node, edges in nx.bfs_successors(g, start):
        for edge in edges:
            if edge_facts == 0:
                break

            graph.add_edge(node, edge)
            edge_facts -= 1

    return graph

def read_graph(path, name):
    g = nx.read_edgelist(path)
    g.name = name
    return g

def create_scale_free_graph(n):
    g = nx.scale_free_graph(n)
    g.name = 'scale-free-{}'.format(n)
    return nx.Graph(g)

def draw_graph(g, pert, layout):
    p = int(pert*100)
    nx.draw_networkx(g, pos=layout, node_size=5, font_size=2, font_color='b', arrowsize=3)
    plt.draw()
    plt.savefig('img/pert_{}.png'.format(p), dpi=500)
    plt.close()

def giant_comp(g):
    return max((g.subgraph(c) for c in nx.connected_components(g)), key=len)
def get_measurements(g):
    infos = pd.Series()
    infos['nodes'] = len(g)
    infos['edges'] = len(g.edges())
    infos['diameter'] = nx.diameter(giant_comp(g))
    infos['betweenness'] = median(nx.betweenness_centrality(g).values())
    infos['clustering'] = median(nx.clustering(g).values())
    infos['components'] = nx.number_connected_components(g)
    infos['closeness'] = median(nx.closeness_centrality(g).values())
    return infos

def perturbation(graph, p):
    g = graph.copy()
    edges_to_remove = int(len(g.edges()) * p)
    removed_edges = []
    for i in range(edges_to_remove):
        random_edge = choice(list(g.edges()))
        g.remove_edges_from([random_edge])
        removed_edges.append(random_edge)
    while(edges_to_remove > 0):
        first_node = choice(list(g.nodes()))
        second_node = choice(list(g.nodes()))
        if(second_node == first_node):
            continue
        if g.has_edge(first_node, second_node) or (first_node, second_node) in removed_edges or (second_node, first_node) in removed_edges:
            continue
        else:
            g.add_edge(first_node, second_node)
            edges_to_remove -= 1
    return g


def deanonymize_vertexref(g, pert_g, i):
    vertexref = vertex_ref(g, i)
    vertexref_pert = vertex_ref(pert_g, i)
    eq = eq_class(vertexref)
    eq_pert = eq_class(vertexref_pert)
    result_eq = {}
    for index in range(0, max(len(eq), len(eq_pert))):
        result = []
        if index < len(eq_pert):
            for value in eq_pert[index]:
                if index < len(eq): 
                    if value in eq[index]:
                        result.append(value)
            result_eq[index] = result
    return deanonymize(result_eq.values(), len(pert_g.nodes))


def deanonymize_edgefacts(g, pert_g, edge_facts):
    edgefacts = edge_facts_subgraph(g, pert_g, edge_facts)
    eq = eq_class(edgefacts)

    return deanonymize(eq, len(pert_g.nodes))
def deanonymize(eq, nodes_number):
    out = {}    
    out['20-inf'] = get_k_anonimity_range(eq,21,inf) / nodes_number
    out['11-20']  = get_k_anonimity_range(eq,11,20)  / nodes_number
    out['5-10']   = get_k_anonimity_range(eq,5,10)   / nodes_number
    out['2-4']    = get_k_anonimity_range(eq,2,4)    / nodes_number
    out['1']      = get_k_anonimity_range(eq,1,1)    / nodes_number
    return out
def eq_class(facts: dict):
    eq_class = {}
    for key, degrees in facts.items():
        k = tuple(sorted(degrees))
        
        if k not in eq_class:
            eq_class[k] = []
        eq_class[k].append(key)
    return list(eq_class.values())

def get_k_anonimity_range(vals, minv, maxv):
    result = 0
    if len(vals) > 0:
        for value in vals:
            if len(value) >= minv and len(value) <= maxv:
                result += len(value)
    return result


# In[4]:


g =  create_scale_free_graph(100)
#g = read_graph('enron.txt','eron')
def draw_edge_factor(ef, pert):
    for i in range (0, 5):
        x, y = zip(*ef[i].items())        
        plt.plot(x, y, label="edge_facts=" + str((i+1) * 10))
    plt.title("Perturbation: {:.0%}".format(pert))
    plt.legend()
    plt.savefig('img/edge_facts_pert_{}.png'.format(int(pert * 100)), dpi=500)        
    plt.close()
def draw_vrtx_rf(vr, pert):
    for i in range (0, 5):
        x, y = zip(*vr[i].items())        
        plt.plot(x, y,label="i=" + str(i))
    plt.title("Perturbation: {:.0%}".format(pert))
    plt.legend()
    plt.savefig('img/vertex_ref_pert_{}.png'.format(int(pert * 100)), dpi=500)        
    plt.close() 
vr_range = [1, 2, 3, 4, 5]
ef_range = [10, 20, 30, 40, 50]
for pert in [0, 0.02, 0.05, 0.1]:
    print('\nPerturbation ({:.0%})'.format(pert))
    pert_g = perturbation(g, pert)
    vr = [deanonymize_vertexref(g, pert_g, i) for i in vr_range]
    ef = [deanonymize_edgefacts(g, pert_g, n) for n in ef_range]
    draw_vrtx_rf(vr, pert)
    draw_edge_factor(ef, pert)
    draw_graph(pert_g, pert, nx.spring_layout(g))
    

print("Executed successfully")


# In[ ]:





# In[ ]:




