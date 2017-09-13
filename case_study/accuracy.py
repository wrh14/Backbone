import matplotlib.pyplot as plt

# f1 = open('../../Output/venue/case_study_M/final_ifCS=1_k=20_alpha=0.9_N=20_M=1_ifSub=1_ifDelta=1','rb')
# f2 = open('../../Output/venue/case_study_M/final_ifCS=1_k=20_alpha=0.9_N=20_M=10_ifSub=1_ifDelta=1','rb')
# f3 = open('../../Output/venue/case_study_M/final_ifCS=1_k=20_alpha=0.9_N=20_M=100_ifSub=1_ifDelta=1','rb')
# f4 = open('../../Output/venue/case_study_M/final_ifCS=0_k=20_alpha=0.9_N=20_M=10_ifSub=1_ifDelta=1','rb')

f1 = open('../../Output/citation_data_mining/case_study_M/final_ifCS=1_k=10_alpha=0.9_N=100_M=1_ifSub=1_ifDelta=1','rb')
f2 = open('../../Output/citation_data_mining/case_study_M/final_ifCS=1_k=10_alpha=0.9_N=100_M=10_ifSub=1_ifDelta=1','rb')
f3 = open('../../Output/citation_data_mining/case_study_M/final_ifCS=1_k=10_alpha=0.9_N=100_M=100_ifSub=1_ifDelta=1','rb')
f4 = open('../../Output/citation_data_mining/case_study_M/final_ifCS=1_k=10_alpha=0.9_N=100_M=1000_ifSub=1_ifDelta=1','rb')
f5 = open('../../Output/citation_data_mining/case_study_M/final_ifCS=0_k=10_alpha=0.9_N=100_M=10_ifSub=1_ifDelta=1','rb')

x1 = []
x2 = []
x3 = []
x4 = []
x5 = []
y1 = []
y2 = []
y3 = []
y4 = []
y5 = []

i = 1

line = f5.readline()
line = f5.readline()
while line:
	line= line.strip('\n')
	row = line.split(' ')
	line = f5.readline()
	x5.append(i)
	i = i + 1
	y5.append(float(row[1]))

f5.close()

i = 1

line = f1.readline()
line = f1.readline()
while line:
	line= line.strip('\n')
	row = line.split(' ')
	line = f1.readline()
	x1.append(i)
	y1.append(float(row[1]))
	i = i + 1

f1.close()

i = 1

line = f2.readline()
line = f2.readline()
while line:
	line= line.strip('\n')
	row = line.split(' ')
	line = f2.readline()
	x2.append(i)
	y2.append(float(row[1]))
	i = i + 1

f2.close()

i = 1

line = f3.readline()
line = f3.readline()
while line:
	line= line.strip('\n')
	row = line.split(' ')
	line = f3.readline()
	x3.append(i)
	y3.append(float(row[1]))
	i = i + 1

f3.close()

i = 1

line = f4.readline()
line = f4.readline()
while line:
	line= line.strip('\n')
	row = line.split(' ')
	line = f4.readline()
	x4.append(i)
	y4.append(float(row[1]))
	i = i + 1

f4.close()

plt.figure()  
plt.plot(x1,y1, lw=0.2, label = 'greedy-CS(1)', marker = 'o', ms = 1)
plt.plot(x2,y2, lw=0.2, label = 'greedy-CS(10)', marker = '<', ms = 1)
plt.plot(x3,y3, lw=0.2, label = 'greedy-CS(100)', marker = '>', ms = 1)
plt.plot(x4,y4, lw=0.2, label = 'greedy', marker = '+', ms = 1) 
plt.plot(x5,y5, lw=0.2, label = 'greedy-CS(1000)', marker = '+', ms = 1) 
plt.legend() 
plt.savefig("cs_accuracy_Citation_DM.jpg",dpi=100) 