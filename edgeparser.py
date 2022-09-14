import networkx as nx
import os
import gt
import matplotlib.pyplot as plt

folder = 'NetPath-pathways'

egfr1 = nx.DiGraph()
tgfbeta = nx.DiGraph()
tnfalpha = nx.DiGraph()
wnt = nx.DiGraph()
graphlist = [egfr1, tgfbeta, tnfalpha, wnt]

nodes = []
edges = []
node = False

for file in os.listdir(folder):
    f = os.path.join(folder, file)
    if os.path.isfile(f):
        with open(f) as text:
            lines = text.readlines()
            if node:
                listofnodes = dict()
                for n in lines:
                    n = n.split()
                    if(n[0] != '#node'):
                        if(n[1] == 'target'):
                            listofnodes[n[0]] = 'tf'
                        else:
                            listofnodes[n[0]] = n[1]
                nodes.append(listofnodes)
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
        #if(edge[5] == 'physical'):
        #    graphlist[counter].add_edge(edge[1], edge[0])
    counter +=1

whpin = nx.Graph()
if os.path.isfile('WHPIN.txt'):
    with open('WHPIN.txt') as text:
        lines = text.readlines()
        i = 1
        for i in range(len(lines)):
            edge = lines[i].split()
            whpin.add_edge(edge[0], edge[1], weight=edge[2])
            

is_subgraph = lambda G, H: all(h_edge in G.edges() for h_edge in H.edges())

print(is_subgraph(whpin, graphlist[0]))
print(is_subgraph(whpin, graphlist[1]))
print(is_subgraph(whpin, graphlist[2]))
print(is_subgraph(whpin, graphlist[3]))


for graph in graphlist:
    counter = 0
    totalcounter = 0
    for ed in graphlist[3].edges():
        if(ed not in whpin.edges()):
            counter+=1
            #print(ed)
        totalcounter+=1
    print(counter == 0)
