from socket import create_server
import networkx as nx
from networkx import dijkstra_path
import scipy
from scipy import sparse
from scipy.sparse import csr_array
import numpy
import os
from collections import OrderedDict
import matplotlib.pyplot as plt

# some code that does not need to be run every time
# may be commented out to speed up the program

# --------------------- Part 1 --------------------------
folder = 'NetPath-pathways'

egfr1 = nx.DiGraph()
tgfbeta = nx.DiGraph()
tnfalpha = nx.DiGraph()
wnt = nx.DiGraph()
graphlist = [egfr1, tgfbeta, tnfalpha, wnt]

nodes = []
edges = []
node = False

# reads files from NetPath-pathways folder and sorts them 
# based on whether they are describing nodes or edges
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

# adds edges to the digraphs 
for edgeset in edges:
    i = 1
    for i in range(len(edgeset)):
        edge = edgeset[i].split()
        if(edge[0] != '-' and edge[1] != '-'):
            graphlist[counter].add_edge(edge[0], edge[1])
            if(edge[5] == 'physical'):
                graphlist[counter].add_edge(edge[1], edge[0])
    counter +=1

# reads from the whpin file and adds corresponding edges to the whpin graph
whpin = nx.DiGraph()
if os.path.isfile('WHPIN.txt'):
    with open('WHPIN.txt') as text:
        lines = text.readlines()
        i = 1
        for i in range(len(lines)):
            edge = lines[i].split()
            whpin.add_edge(edge[0], edge[1])
            
# sanity check (all results are true so it has been commented out
# to speed up the program)
"""
for graph in graphlist:
    counter = 0
    totalcounter = 0
    for ed in graph.edges():
        if(ed not in whpin.edges()):
            counter+=1
            print(ed)
        totalcounter+=1
    print(counter == 0)
"""
# --------------------- Part 2 --------------------------
# structures to hold the paths constructed by PathLinker
pegfr1 = OrderedDict()
ptgfbeta = OrderedDict()
ptnfalpha = OrderedDict()
pwnt = OrderedDict()

predicts = [pegfr1, ptgfbeta, ptnfalpha, pwnt]

# parses the PathLinker data sets
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

# ------------------------------ Part 3 -------------------------------------
# computes precision and recall for PathLinker for the 4 predicted graphs

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
            # precision
            results[i * 2].append(positives/counter)
            # recall
            results[(i * 2) + 1].append(positives / p)


# the code below was used to plot the precision recall curves

plt.plot(results[1], results[0], label="EGFR1")
plt.plot(results[3], results[2], label="TGF_beta")
plt.plot(results[5], results[4], label="TGF_alpha")
plt.plot(results[7], results[6], label="Wnt")
plt.legend()
plt.title("Evaluation of Reconstructed Pathways (PathLinker)")
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.savefig("precision_recall_PathLinker.png")


# ------------------------------ Part 4 -------------------------------------
# collecting tfs and receptors
tfs = [[], [], [], []]
sources = [[], [], [], []]
counter = 0
for graph in nodes:
    for node in graph.keys():
        if(graph[node] == 'tf'):
            tfs[counter].append(node)
        elif(graph[node] == 'receptor'):
            sources[counter].append(node)
    counter += 1

# using dijkstra's algorithm to find paths
paths = [nx.DiGraph(),nx.DiGraph(),nx.DiGraph(),nx.DiGraph()]
counter = 0
for i in range(len(tfs)):
    for tf in tfs[i]:
        for source in sources[i]:
            try:
                path = (dijkstra_path(whpin, source, tf))
                for j in range(len(path) - 1):
                    paths[i].add_edge(path[j], path[j + 1])
            except:
                counter+=1

# calculate and graph the precision recall curve

results = [[], [], [], [], [], [], [], []]

for i in range(len(paths)):
    p = len(graphlist[i].edges)
    counter = 1
    positives = 0
    for edge in paths[i].edges():
        if edge in graphlist[i].edges:
            positives+=1
        counter+=1
    # precision
    results[i * 2].append(positives/counter)
    # recall
    results[(i * 2) + 1].append(positives / p)

# plot precision recall points for the shortest path algorithm

plt.clf()
plt.plot(results[1], results[0], label="EGFR1", marker="o", markersize=10)
plt.plot(results[3], results[2], label="TGF_beta", marker="o", markersize=10)
plt.plot(results[5], results[4], label="TGF_alpha", marker="o", markersize=10)
plt.plot(results[7], results[6], label="Wnt", marker="o", markersize=10)
plt.legend()
plt.title("Evaluation of Reconstructed Pathways (Shortest Paths)")
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.savefig("precision_recall_ShortestPaths.png")

# ------------------------------ Part 5 -------------------------------------
# implementing RWR
matrix = nx.to_numpy_array(whpin)

dmatrix = numpy.tril(matrix)
dmatrix = numpy.triu(dmatrix)

# unfortunately I was not able to get the below line to work
# so I was not able to complete the RWR graph
# whenever I tried to take the inverse of the diagonalized
# matrix it outputted it as a singular matrix
# numpy.transpose(numpy.matmul(numpy.inv(dmatrix), matrix))