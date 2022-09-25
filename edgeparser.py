import networkx as nx
import os
from collections import OrderedDict
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
        if(edge[0] != '-' and edge[1] != '-'):
            graphlist[counter].add_edge(edge[0], edge[1])
            if(edge[5] == 'physical'):
                graphlist[counter].add_edge(edge[1], edge[0])
    counter +=1

whpin = nx.Graph()
if os.path.isfile('WHPIN.txt'):
    with open('WHPIN.txt') as text:
        lines = text.readlines()
        i = 1
        for i in range(len(lines)):
            edge = lines[i].split()
            whpin.add_edge(edge[0], edge[1], weight=edge[2])
            
# sanity check
for graph in graphlist:
    counter = 0
    totalcounter = 0
    for ed in graph.edges():
        if(ed not in whpin.edges()):
            counter+=1
            print(ed)
        totalcounter+=1
    print(counter == 0)

pegfr1 = OrderedDict()
ptgfbeta = OrderedDict()
ptnfalpha = OrderedDict()
pwnt = OrderedDict()

predicts = [pegfr1, ptgfbeta, ptnfalpha, pwnt]

paths = True
predictededges = []
for file in os.listdir("PathLinker-results"):
    f = os.path.join("PathLinker-results", file)
    if os.path.isfile(f):
        with open(f) as text:
            if not paths:
                lines = text.readlines()
                predictededges.append(lines)
                paths = True

            else:
                paths = False
counter = 0
for edgeset in predictededges:
    edgenum = 1
    for i in range(1,len(edgeset)):
        edge = edgeset[i].split()
        
        if(edge[0] != '-' and edge[1] != '-'):
            if int(edge[2]) in predicts[counter].keys():
                predicts[counter][int(edge[2])].add_edge(edge[0], edge[1])
            else:
                predicts[counter][int(edge[2])] = nx.DiGraph()
                predicts[counter][int(edge[2])].add_edge(edge[0], edge[1])
            
            edgenum = int(edge[2]) + 1
    counter +=1

results = [[], [], [], [], [], [], [], []]

for i in range(len(predicts)):
    p = len(graphlist[i].edges)
    counter = 1
    positives = 0
    for graph in predicts:
        for key in graph.keys():
            ed = graph[key].edges
            for edge in ed:
                if edge in graphlist[i].edges:
                    positives+=1
                counter+=1
            results[i * 2].append(positives/counter)
            results[(i * 2) + 1].append(positives / p)

plt.plot(results[1], results[0], label="EGFR1")
plt.plot(results[3], results[2], label="TGF_beta")
plt.plot(results[5], results[4], label="TGF_alpha")
plt.plot(results[7], results[6], label="Wnt")
plt.legend()
plt.title("Evaluation of Reconstructed Pathways")
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.show()