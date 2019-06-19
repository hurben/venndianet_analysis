#18.04.11
#"visualization_for_cytoscape.py"
#
#not for web-visualization
#
#

def File_to_Readlines(file_name):

	file_open = open(file_name,'r')
	file_readlines = file_open.readlines()

	return file_readlines

def Create_Topology(adj_matrix_file):
#input : adj.matrix
	topology_dict = {}

	adj_matrix_readlines = File_to_Readlines(adj_matrix_file)

	#make gene list
	read = adj_matrix_readlines[0].replace('\n','')
	gene_list = read.split('\t')
	gene_list.remove('')
	print "Total # of genes in adjacency matrix :", len(gene_list)

	#make dict

	for i in range(1, len(adj_matrix_readlines)):
		read = adj_matrix_readlines[i]
		read = read.replace('\n','')

		token = read.split('\t')

		node_gene_1 = token[0]

		for j in range(1, len(token)):
			edge_existance = int(token[j])

			if edge_existance == 1:
				GENE_LIST_INDEX = j-1
				node_gene_2 = gene_list[GENE_LIST_INDEX]

				try : topology_dict[node_gene_1].append(node_gene_2)
				except KeyError: topology_dict[node_gene_1] = [node_gene_2]

	return topology_dict

def Create_P0_list(p0_file):
	p0_list = []
	p0_file_readlines = File_to_Readlines(p0_file)

	for i in range(len(p0_file_readlines)):
		read = p0_file_readlines[i]
		read = read.replace('\n','')
		token = read.split('\t')

		gene = token[0]
		value = int(token[1])
		if value == 1:
			p0_list.append(gene)

	return p0_list


def Create_RWR_Probability_Vector(rwr_file):
#input rwr
	rwr_dict = {}
	rwr_summary_file = str(rwr_file) +'.summary'

	rwr_file_readlines = File_to_Readlines(rwr_file)
	rwr_summary_file_readlines = File_to_Readlines(rwr_summary_file)

	for i in range(1, len(rwr_file_readlines)):
		read = rwr_file_readlines[i]
		read = read.replace('\n','')
		token = read.split()

		gene = token[0]
		value = token[1]

		if gene in rwr_dict.keys():
			"CRITICAL ERROR #1"
			quit()

		if gene == '2010111I01Rik':
			print gene
		rwr_dict[gene] = [value, 0]

	for i in range(len(rwr_summary_file_readlines)):
		read = rwr_summary_file_readlines[i]
		read = read.replace('\n','')
		gene = read.split('_')[0]
		rwr_dict[gene][1] = 1

	return rwr_dict


def Create_Result_Table(rwr_dict,topology_dict, p0_list, common_file_name):

	topology_txt = open(str(common_file_name) + '.topology.txt','w')
	topology_txt.write('gene\ttarget\n')
	topology_mini_txt = open(str(common_file_name) + '.topology.rwr.cut.txt','w')
	topology_mini_txt.write('gene\ttarget\n')

	for node1 in topology_dict.keys():

		for node2 in topology_dict[node1]:
			topology_txt.write(str(node1) + '\t' + str(node2) +'\n')

			if node1 in rwr_dict.keys():
				if node2 in rwr_dict.keys():
					if rwr_dict[node1][1] == 1:
						if rwr_dict[node2][1] == 1:
							topology_mini_txt.write(str(node1) + '\t' + str(node2) +'\n')
		
	p0_seeds_txt = open(str(common_file_name) + '.p0.seeds.txt','w')
	p0_seeds_txt.write('gene\tseeds\n')
	for gene in p0_list:
		p0_seeds_txt.write(str(gene) +'\t'+'1\n')



	rwr_result_txt = open(str(common_file_name) + '.rwr.result.cutoff.txt','w')
	rwr_result_txt.write('gene\trwr_score\tcut\n')
	for gene in rwr_dict.keys():
		score = rwr_dict[gene][0]
		cut = rwr_dict[gene][1]
		rwr_result_txt.write(str(gene) + '\t' + str(score) + '\t' + str(cut) +'\n')


	rwr_result_txt = open(str(common_file_name) + '.rwr.result.cyto.txt','w')
	rwr_result_txt.write('gene\ttarget\trwr_score\tseed\tcut\t\n')
	for node1 in topology_dict.keys():
		score = rwr_dict[gene][0]
		cut = rwr_dict[gene][1]
		seed = 0
		if node1 in p0_list:
			seed = 1
		for node2 in topology_dict[node1]:
			rwr_result_txt.write(str(node1) + '\t' + str(node2) + '\t' + str(score) + '\t' + str(seed) +'\t' + str(cut) +'\n')

			


import sys

mlv_file = sys.argv[1]
common_file_name = sys.argv[2]

adj_matrix_file = str(common_file_name) + '.adj.matrix'
p0_vector_file = str(common_file_name) + '.p0.vector' #a.k.a : intersection genes
rwr_vector_file = str(common_file_name) + '.rwr'
rwr_summary_vector_file = str(rwr_vector_file) + '.summary'

p0_list = Create_P0_list(p0_vector_file)
rwr_dict = Create_RWR_Probability_Vector(rwr_vector_file)
topology_dict = Create_Topology(adj_matrix_file)

Create_Result_Table(rwr_dict, topology_dict, p0_list, common_file_name)
