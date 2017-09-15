import networkx as nx   
import matplotlib.pyplot as plt  
import csv
import sys  
import random

reload(sys)  
sys.setdefaultencoding('utf8')

raw_top_nodes = {}

G = nx.DiGraph()
pos = {}


# f2 = open('../../Output/citation_data_mining/case_study_M/final_ifCS=0_k=10_alpha=0.9_N=100_M=10_ifSub=1_ifDelta=1','rb')
f2 = open('../../Output/venue/case_study_M/final_ifCS=0_k=20_alpha=0.9_N=20_M=10_ifSub=1_ifDelta=1','rb')
line = f2.readline()
line = f2.readline()
while line :
	row = line.split(' ')
	line = f2.readline()
	raw_top_nodes[int(row[0])] = 'top_list'
f2.close()

f1 = open('../../emb/venue/venue_indexNetworkMin10.txt','rb')
line = f1.readline()
while line:
	row = line.split(' ')
	line = f1.readline()
	pos[int(row[0])] = (float(row[1]), float(row[2]))
f1.close()

f1 = open('../../Dataset/venue/indexNetwork.txt','rb')
line = f1.readline()
while line:
	row = line.split(' ')
	if pos.has_key(int(row[0])) and pos.has_key(int(row[1])) and int(row[2]) >= 16:
		G.add_edge(int(row[0]), int(row[1]))
	line = f1.readline()
f1.close()

labels = {}
f1 = open('../../Dataset/venue/indexToName.txt','rb')
line = f1.readline()
while line:
	row = line.split(' ')
	labels[int(row[0])] = row[1].split('\t')[1]
	line = f1.readline()
f1.close()


plt.figure(figsize = (10,10), dpi = 500)
nx.draw(G, pos, node_color = 'grey', edge_color = 'grey', arrows=False, node_size = 1, alpha = 0.3)

G2 = nx.DiGraph()

top_labels = {}
for i in raw_top_nodes:
	G2.add_node(i)
	top_labels[i] = labels[i]

print(len(G2.nodes()))
nx.draw(G2, pos, labels=top_labels, node_color = 'red', node_size = 3, alpha = 1, font_size=10)

# plt.savefig('backbone_citation-DM')
plt.savefig('backbone_venue')