#Purpose : result representation after summary
#just to visualize 
import sys


def READLINES_TO_TOKEN(readlines, i, delimit):
	read = readlines[i]
	read = read.replace('\n','')

	if delimit == '':
		token = read.split()
	else:
		token = read.split(delimit)

	return read, token



def ppi_to_dict():

	ppi_dict = {}

	ppi_profile_file = open('//var/www/htdocs/MLV/src/ppi/STRING_PPI_Mus_musculus_Symbol.txt','r')
	ppi_profile_file_readlines = ppi_profile_file.readlines()

	for i in range(len(ppi_profile_file_readlines)):
		read, token = READLINES_TO_TOKEN(ppi_profile_file_readlines, i, '\t')
		node_1 = token[0]
		node_2 = token[1]

		try : ppi_dict[node_1].append(node_2)
		except KeyError : ppi_dict[node_1] = [node_2]

	return ppi_dict


def organize_mlv_file(mlv_file):

	mlv_open = open(mlv_file)
	mlv_readlines = mlv_open.readlines()

	mlv_dict = {}

	mlv_intersection_dict = {}
	mlv_complementary_dict = {}

	for i in range(len(mlv_readlines)):

		read, token = READLINES_TO_TOKEN(mlv_readlines, i, '\t')

		condition = token[0]

		if '__' not in condition:

			for j in range(1, len(token)):
				gene = token[j]
				genesymbol = gene.split('_')[0]
				try : mlv_dict[condition].append(gene)
				except KeyError : mlv_dict[condition] = [gene]

		if '__' in condition:

			for j in range(1, len(token)):
				gene = token[j]
				genesymbol = gene.split('_')[0]
				try : mlv_intersection_dict[condition].append(gene)
				except KeyError : mlv_intersection_dict[condition] = [gene]

	for condition in mlv_dict.keys():
		gene_list = mlv_dict[condition]
		new_gene_list = gene_list[:]

		for gene in gene_list:
			for intersection_condition in mlv_intersection_dict.keys():

				intersection_gene_list = mlv_intersection_dict[intersection_condition]

				if gene in intersection_gene_list:
					new_gene_list.remove(gene)
					break

		mlv_complementary_dict[condition] = new_gene_list

	for condition in mlv_intersection_dict.keys():
		mlv_dict[condition] = mlv_intersection_dict[condition]

	return mlv_dict, mlv_complementary_dict, mlv_intersection_dict


def orgnanize_rwr_files(rwr_file_list, ppi_dict):

	rwr_info_dict = {}

	rwr_ppi_dict = {}
	rwr_ppi_info_dict = {}
	rwr_ppi_cut_dict = {}

	open_rwr_file_list = open(rwr_file_list,'r')
	rwr_file_list_readlines = open_rwr_file_list.readlines()
	
	for i in range(len(rwr_file_list_readlines)):

		read, token = READLINES_TO_TOKEN(rwr_file_list_readlines, i, '/')
		
		rwr_file = read
		rwr_condition = token[len(token) - 1]

		rwr_file_open = open(rwr_file,'r')
		rwr_file_readlines = rwr_file_open.readlines()

		rwr_summary_file_open = open(str(rwr_file) + '.summary', 'r')
		rwr_summary_readlines = rwr_summary_file_open.readlines()
		
		rwr_info_dict[rwr_condition] = [len(rwr_file_readlines) - 1, len(rwr_summary_readlines)]


		p0_file = read.replace('matrix.rwr','matrix.p0.vector')
		p0_open = open(p0_file,'r')
		p0_readlines = p0_open.readlines()
		
		#p0 information : PART 1
		for j in range(1, len(rwr_file_readlines)):

			read_j, token_j = READLINES_TO_TOKEN(rwr_file_readlines, j, '')

			gene = token_j[0]
			prob = token_j[1]

			try :
				rwr_ppi_dict[rwr_condition].append(gene)
			except KeyError:
				rwr_ppi_dict[rwr_condition] = [gene]
				
			rwr_ppi_info_dict[rwr_condition, gene] = {
			'prob' : prob,
			'intersection' : 0
			}

		#p0 information : PART2 
		for j in range(len(p0_readlines)):

			read_j, token_j = READLINES_TO_TOKEN(p0_readlines, j, '\t')

			gene = token_j[0]
			intersection = token_j[1]

			rwr_ppi_info_dict[rwr_condition, gene]['intersection'] = intersection	


		#rwr cut info
		for j in range(len(rwr_summary_readlines)):

			read_j, token_j = READLINES_TO_TOKEN(rwr_summary_readlines, j, '_')
			gene = token_j[0]

			try : rwr_ppi_cut_dict[rwr_condition].append(gene)
			except KeyError : rwr_ppi_cut_dict[rwr_condition] = [gene]

	return rwr_info_dict, rwr_ppi_dict, rwr_ppi_info_dict, rwr_ppi_cut_dict
		

#############################################################

mlv_file = sys.argv[1]
rwr_file_list = sys.argv[2]

ppi_dict = ppi_to_dict()
mlv_dict, mlv_complementary_dict, mlv_intersection_dict = organize_mlv_file(mlv_file)
rwr_info_dict, rwr_ppi_dict, rwr_ppi_info_dict, rwr_ppi_cut_dict =  orgnanize_rwr_files(rwr_file_list, ppi_dict)


#[1] MLV file for Venndiagram level interpret
	
#[1-1] count number of genes in complementary, intersection

output = open(str(mlv_file) + '.stat','w')
output.write('#CONDITION\tENTRY_NUMBERS\n')

for condition in mlv_complementary_dict.keys():

	gene_list = mlv_complementary_dict[condition]
	output.write(str(condition) + '\t' + str(len(gene_list)) + '\n')
	
for condition in mlv_intersection_dict.keys():
	gene_list = mlv_intersection_dict[condition]
	output.write(str(condition) + '\t' + str(len(gene_list)) + '\n')

output.write('\n\n#CONDITION_INCLUDING_INTERSECTION\tENTRY_NUMBERS\n')

for condition in mlv_dict.keys():
	if '__' not in condition:
		gene_list = mlv_dict[condition]
		output.write(str(condition) + '\t' + str(len(gene_list)) + '\n')
output.close()



#[] RWR info
for condition in rwr_ppi_dict.keys():
	gene_list = rwr_ppi_dict[condition]

	output = open(str(condition) + '.ppi.info.cyto.txt','w')
	output.write('gene\tRWR_probability\tintersection\n')

	for gene in gene_list:

		intersection = rwr_ppi_info_dict[condition, gene]['intersection']
		prob = rwr_ppi_info_dict[condition, gene]['prob']

		output.write(str(gene) + '\t' + str(prob) + '\t' + str(intersection) + '\n')
	
	output.close()
		
#[] create RWR topology
for condition in rwr_ppi_dict.keys():

	edge_count = 0
	output = open(str(condition) + '.ppi.topology.cyto.txt','w')
	gene_list = rwr_ppi_dict[condition]

	for i in range(len(gene_list)):
		for j in range(i + 1, len(gene_list)):

			gene_i = gene_list[i]
			gene_j = gene_list[j]

			if gene_i != gene_j:

				gene_i_list = ppi_dict[gene_i]

				if gene_j in gene_i_list:
					output.write(str(gene_i) + '\t' + str(gene_j) +'\n')
					edge_count += 1
			else :
				print 'ERROR'
				print i,j

	rwr_info_dict[condition].append(edge_count)
	output.close()

#[] create RWR cut topology
for condition in rwr_ppi_cut_dict.keys():

	edge_count = 0
	output = open(str(condition) + '.ppi.cut.topology.cyto.txt','w')
	gene_list = rwr_ppi_cut_dict[condition]

	for i in range(len(gene_list)):
		for j in range(i + 1, len(gene_list)):

			gene_i = gene_list[i]
			gene_j = gene_list[j]

			if gene_i != gene_j:

				gene_i_list = ppi_dict[gene_i]

				if gene_j in gene_i_list:
					output.write(str(gene_i) + '\t' + str(gene_j) +'\n')
					edge_count += 1
			else :
				print 'ERROR'
				print i,j

	output.close()
	rwr_info_dict[condition].append(edge_count)

# RWR summary
output = open(str(mlv_file) + '.rwr.sum','w')
output.write('#CONDITION\tRWR_BEFORE\tRWR_AFTER\tRWR_BEFORE_EDGE\tRWR_AFTER_EDGE\n')
for condition in rwr_info_dict.keys():
	rwr_count = rwr_info_dict[condition][0]
	rwr_summ_count = rwr_info_dict[condition][1]
	rwr_edge_count = rwr_info_dict[condition][2]
	rwr_summ_edge_count = rwr_info_dict[condition][3]
	output.write(str(condition) + '\t' + str(rwr_count) + '\t' + str(rwr_summ_count) + '\t' + str(rwr_edge_count) + '\t' +str(rwr_summ_edge_count)+ '\n')

output.close()










