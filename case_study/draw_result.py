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
f2 = open('../../Output/citation_data_mining/baseline/page_rank_final_size=100_alpha=0.9_threshold=1e-15','rb')
line = f2.readline()
line = f2.readline()
while line :
	row = line.split(' ')
	line = f2.readline()
	raw_top_nodes[int(row[0])] = 'top_list'
f2.close()
print('finish reading positions')

f1 = open('../../emb/citation_data_mining/data_mining_indexNetwork.txt','rb')
line = f1.readline()
while line:
	row = line.split(' ')
	line = f1.readline()

r = 100

f1 = open('../../emb/citation_data_mining/data_mining_indexNetwork.txt','rb')
line = f1.readline()
while line:
	row = line.split(' ')
	line = f1.readline()
	if random.uniform(0,1) < 1:
		if abs(float(row[1]) + 5) < r and abs(float(row[2]) + 2) < r - 1:
			pos[int(row[0])] = (float(row[1]), float(row[2]))
			G.add_node(int(row[0]))
	else:
		if raw_top_nodes.get(int(row[0])) == 'top_list':
			pos[int(row[0])] = (float(row[1]), float(row[2]))
			G.add_node(int(row[0]))
f1.close()

print('finish reading nodes')


plt.figure(figsize = (10,10), dpi = 100)
nx.draw(G, pos, node_color = 'grey', node_size = 5, alpha = 0.5)

G2 = nx.DiGraph()

for i in raw_top_nodes:
	G2.add_node(i)

print(len(G2.nodes()))
nx.draw(G2, pos, node_color = 'red', node_size = 20, alpha = 1)

# plt.savefig('backbone_citation-DM')
plt.savefig('page_rank_citation-DM')