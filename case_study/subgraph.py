import networkx as nx 
import matplotlib.pyplot as plt

core_id = 343
node_num = 44161
f = open("../../emb/citation_data_mining/data_mining_indexNetwork.txt", "r")
line = f.readline() 
emb = [None] * node_num
line = f.readline()
while line:
	args = line.split(" ")
	emb[int(args[0])] = (float(args[1]), float(args[2]))
	line = f.readline()
f.close()
print "Read emb successfully."


G = nx.DiGraph()
pos = {}
fw = open("../../Output/citation_data_mining/subgraph/sub_index_network.txt", "w")
r1 = 2
r2 = 3
if_sub = [False] * node_num
for i in range(node_num):
	if emb[i] is not None:
		if abs(emb[i][0] - emb[core_id][0]) < r1 and abs(emb[i][1] - (emb[core_id][1] - 2)) < r2:
			if_sub[i] = True
			# fw.write(str(i) + " " + str(emb[i][0]) + " " + str(emb[i][1]) + "\n")
			G.add_node(i)
			pos[i] = emb[i]

print len(G.nodes())
print len(pos.keys())
f = open("../../Dataset/citation_data_mining/indexNetwork.txt", "r")
line = f.readline()
while line:
	args = line.split(" ")
	if if_sub[int(args[0])] and if_sub[int(args[1])]:
		fw.write(line)
		G.add_edge(int(args[0]), int(args[1]), weight=float(args[2]))
	line = f.readline()

plt.figure(figsize = (30,30), dpi = 30)
nx.draw(G, pos=pos, node_size=40, edge_color='grey', node_color='grey', edge_cmap=plt.cm.Blues, alpha = 0.5)

# G2 = nx.DiGraph()
# G2.add_node(1333)
# G2.add_node(2243)
# G2.add_node(509)
# G2.add_node(343)
# labels = {}
# labels[1333] = "A"
# labels[2243] = "B"
# labels[509] = "C"
# labels[343] = "D"
G2 = nx.DiGraph()
G2.add_node(1430)
G2.add_node(5856)
G2.add_node(21762)
G2.add_node(343)
labels = {}
labels[1430] = "A"
labels[5856] = "B"
labels[21762] = "C"
# labels[343] = "D"
nx.draw(G2, pos, node_color='b', node_size=400, alpha = 1, labels=labels, font_color = 'purple', font_size = 100)
#plt.savefig('backbone_subgraph')
plt.savefig('diversified_rank_subgraph')