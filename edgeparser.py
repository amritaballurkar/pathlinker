import networkx as nx
import os
import matplotlib.pyplot as plt

folder = 'NetPath-pathways'

egfr1 = nx.DiGraph()
tgfbeta = nx.DiGraph()
tnfalpha = nx.DiGraph()
wnt = nx.DiGraph()
graphlist = [egfr1, tgfbeta, tnfalpha, wnt]

edges = []
node = False
for file in os.listdir(folder):
    f = os.path.join(folder, file)
    if os.path.isfile(f):
        
        with open(f) as text:
            lines = text.readlines()
            if node:
                node = False
            else:
                edges.append(lines)
                node = True

counter = 0

for edgeset in edges:
    i = 1
    for i in range(len(edgeset)):
        edge = edgeset[i].split()
        graphlist[counter].add_edge(edge[0], edge[1])
        if(edge[5] == 'physical'):
            graphlist[counter].add_edge(edge[1], edge[0])
    counter +=1


