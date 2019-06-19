def Open_N_Readlines(single_file):

	open_file = open(single_file,'r')
	file_readlines = open_file.readlines()

	return file_readlines

def MLV_to_Dict(mlv_file):

	condition_id_dict = {}
	gene_condition_dict = {}
	condition_gene_dict = {}

	mlv_readlines = Open_N_Readlines(mlv_file)

	check_point = 0

	for i in range(1, len(mlv_readlines)):

		if check_point == 1:

			read = mlv_readlines[i]
			read = read.replace('\n','')
			token = read.split('\t')
			conditionID = token[0]

			for j in range(1,len(token)):
				gene = token[j]
				gene = gene.split('_')[0]

				try : condition_gene_dict[conditionID].append(gene)
				except KeyError: condition_gene_dict[conditionID] = [gene]

				gene_condition_dict[gene] = conditionID

		if check_point == 0:

			if mlv_readlines[i][0:3] != '###':

				read = mlv_readlines[i]
				read = read.replace('\n','')
				token = read.split('\t')

				print token
				conditionID = token[0]
				condition_tag = token[1]
				condition_id_dict[conditionID] = condition_tag

			else:
				check_point = 1

	return condition_id_dict, gene_condition_dict, condition_gene_dict

def Summary_to_Dict(summary_file, condition_of_interest_list):

	import operator

 	summary_readlines = Open_N_Readlines(summary_file)
	check = 0
	prob_dict = {}
	condition_dict = {}

	for i in range(len(summary_readlines)):
		read = summary_readlines[i]
		read = read.replace('\n','')
		token = read.split('\t')

		if check == 1:
			gene = token[0]
			prob = float(token[1])
			condition = token[2]


			if condition in condition_of_interest_list:

				prob_dict[gene] = prob
				condition_dict[gene] = condition


		if token[0] == 'Gene':
			check = 1
		

	sorted_prob_gene_list = sorted(prob_dict.items(), key=operator.itemgetter(1))


	return prob_dict, condition_dict, sorted_prob_gene_list


def TP_file_to_list(tp_file):
	tp_list = []

	tp_file_readlines = Open_N_Readlines(tp_file)
	for i in range(len(tp_file_readlines)):
		read = tp_file_readlines[i]
		read = read.replace('\n','')
		gene = read
		tp_list.append(gene)

	return tp_list
	


if __name__ == '__main__':

	import sys

	mlv_file = sys.argv[1]
	summary_file = sys.argv[2]
	tp_file = sys.argv[3]
	condition_of_interest = sys.argv[4]
	condition_of_interest_list = condition_of_interest.split('+')

	print condition_of_interest_list

	result_txt = open(str(summary_file) + '.a1','w')

	condition_id_dict, gene_condition_dict, condition_gene_dict = MLV_to_Dict(mlv_file)
	prob_dict, condition_dict, sorted_prob_gene_list = Summary_to_Dict(summary_file, condition_of_interest_list)
	tp_list = TP_file_to_list(tp_file)


	sorted_gene_list = []

	for gene_prob in sorted_prob_gene_list:
		gene = gene_prob[0]
		sorted_gene_list.append(gene)


	result_txt.write('INDEPENDANT\t' + str(len(prob_dict.keys())) + '\n')
	result_txt.write('TP_GENE\tRANK\tCONDITION\n')

	for tp_gene in tp_list:

		try : rank = sorted_gene_list.index(tp_gene) + 1
		except ValueError: rank = 'N/A'

		try : condition = condition_dict[tp_gene]
		except KeyError: condition = 'N/A'

		result_txt.write(str(tp_gene) + '\t' + str(rank) + '\t' + str(condition) + '\n')
	
	result_txt.write('DEPENDANT\n')
	result_txt.write('TP_GENE\tRANK\tCONDITION\n')


	reverse_sorted_gene_list = sorted_gene_list[::-1]

	for tp_gene in tp_list:

		try : rank = reverse_sorted_gene_list.index(tp_gene) + 1
		except ValueError: rank = 'N/A'

		try : condition = condition_dict[tp_gene]
		except KeyError: condition = 'N/A'

		result_txt.write(str(tp_gene) + '\t' + str(rank) + '\t' + str(condition) + '\n')
	

	print sorted_gene_list
	print reverse_sorted_gene_list






