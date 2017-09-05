selected_node = (1333, 2243, 509, 343)
selected_node = (1430, 5856, 21762, 343)
for i in range(4):
	num = 0
	f = open("../../Dataset/citation_data_mining/indexNetwork.txt", "r")
	line = f.readline()
	while line:
		args = line.split(" ")
		if int(args[1]) == selected_node[i]:
			num = num + 1
		line = f.readline()
	f.close()
	print selected_node[i], num