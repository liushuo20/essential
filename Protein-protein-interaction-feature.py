import networkx as nx
import time


#----------------------define the time format--------------------------------#
timeFormat = '%Y-%m-%d %X'
print (time.strftime(timeFormat, time.localtime()))
#-----------------read data from txt and write it into a graph file---------------------#
ppi_in = open('F:\\project3supp\\basicBAS\\PPI\\224308.protein.links.v10.5.txt')
ppi_out = open('F:\\project3supp\\basicBAS\\PPI\\topo.txt', 'w')

PPI = nx.Graph()
ppi_in.readline()
for line in ppi_in.readlines():
	l = line.strip().split()
	if l[2] >= 400:
		PPI.add_node(l[0])
		if l[0] != l[1]:
			PPI.add_edge(l[0], l[1])
ppi_in.close()

NuOfNodes = PPI.number_of_nodes()
print (NuOfNodes)
print (time.strftime(timeFormat, time.localtime()))
gene = PPI.nodes()

#-------------------------calculate betweenness, closeness, cluster-------------------#
betweenness = nx.betweenness_centrality(PPI)
print ('betweenness is ok')
print (time.strftime(timeFormat, time.localtime()))
closeness = nx.closeness_centrality(PPI)
print ('closeness is ok')
print (time.strftime(timeFormat, time.localtime()))
cluster = nx.clustering(PPI)
print ('clustering is ok')
print (time.strftime(timeFormat, time.localtime()))

#------------------------calculate 2-5 step of coverage-------------------------------# 
stepInter = {}
for node1 in PPI.nodes():
	nodePath = nx.shortest_path(PPI, node1)
	du = PPI.degree(node1)
	s2 = s3 = s4 = s5 = 0
	for node2 in nodePath.keys():
		if len(nodePath[node2]) == 3:
			s2 += 1
		elif len(nodePath[node2]) == 4:
			s3 += 1
		elif len(nodePath[node2]) == 5:
			s4 += 1
		elif len(nodePath[node2]) == 6:
			s5 += 1
	stepInter[node1] = [du/float(NuOfNodes), (du + s2)/float(NuOfNodes), (du + s2 + s3)/float(NuOfNodes), (du + s2 + s3 + s4)/float(NuOfNodes), (du + s2 + s3 + s4 + s5)/float(NuOfNodes)]
print ('all step is ok')
print (time.strftime(timeFormat, time.localtime()))

#-------------------write into file------------------------------#
ppi_out.write('geneName	Betweenness	Closeness clusteringCoefficient	Degree step2	step3	step4	step5' + '\n')
for nod in PPI.nodes():
	str_out = nod + '	' + str(betweenness[nod]) + '	' + str(closeness[nod]) + '	' + str(cluster[nod]) + '	' + \
			str(stepInter[nod][0]) + '	' + str(stepInter[nod][1]) + '	' + str(stepInter[nod][2]) + '	' + str(stepInter[nod][3]) + '	' + str(stepInter[nod][4])+ '\n'
	ppi_out.write(str_out)
print (time.strftime(timeFormat, time.localtime()))
ppi_out.close()