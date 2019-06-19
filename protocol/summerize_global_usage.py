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

def TP_file_to_list(tp_file):
	tp_list = []

	tp_file_readlines = Open_N_Readlines(tp_file)
	for i in range(len(tp_file_readlines)):
		read = tp_file_readlines[i]
		read = read.replace('\n','')
		gene = read
		tp_list.append(gene)

	return tp_list
	

def RWR_Result_to_Dict(rwr_file):

	rwr_dict = {}
	
	RWR_readlines = Open_N_Readlines(rwr_file)

	for i in range(1, len(RWR_readlines)):
		read = RWR_readlines[i]
		read = read.replace('\n','')
		token = read.split()

		gene = token[0]
		prob = token[1]

		rwr_dict[gene] = prob

	return rwr_dict



if __name__ == '__main__':

	import sys

	mlv_file = sys.argv[1]
	rwr_file = sys.argv[2]
	tp_file = sys.argv[3]

	condition_id_dict, gene_condition_dict, condition_gene_dict = MLV_to_Dict(mlv_file)
	rwr_dict = RWR_Result_to_Dict(rwr_file)
	tp_list = TP_file_to_list(tp_file)

	output_txt = open(str(rwr_file)+ '.summary', 'w')

	#Header
	for conditionID in condition_gene_dict.keys():
		output_txt.write('ConditionID\tConditionName\n')
		output_txt.write(str(conditionID) + '\t' + str(condition_id_dict[conditionID]) + '\n')
	output_txt.write('Gene\tProbablity\tCondition\tTP\n')

	for gene in rwr_dict.keys():

		prob = rwr_dict[gene]
		output_txt.write(str(gene) +'\t' + str(prob))

		condition = gene_condition_dict[gene]
		output_txt.write('\t' + str(condition))

		if gene in tp_list:
			output_txt.write('\t1')
		else:
			output_txt.write('\t0')
		
		output_txt.write('\n')



	output_txt.close()



