def text_to_dict(file_dir):

	if '/' in file_dir:
		condition_name = file_dir.split('/')
		slash_n = len(condition_name)
		condition_name = condition_name[slash_n - 1]

	condition_name = condition_name.replace('.sig.gene.matrix','')

	print '--------------------'
	print '| ', condition_name
	print '--------------------'

	open_file = open(file_dir,'r')
	file_readlines = open_file.readlines()
	gene_list = []

	result_txt = open('%s.%s' % (condition_name,'UP'),'w')
	result_txt.write('TS_ID\tGENESYMBOL\tREAL_FC\n')

	for i in range(1, len(file_readlines)):
		read = file_readlines[i]
		read = read.replace('\n','')
		read = read.replace('\r','')
		token = read.split('\t')
		gene = token[1]
		fc = float(token[2])

		if fc > 1:
			result_txt.write('%s\t%s\t%s\n' % (token[0], gene, fc))
	result_txt.close()

if __name__ == '__main__':

	import sys
	input_list_file = sys.argv[1]

	open_input_dir = open(input_list_file,'r')
	input_list_readlines = open_input_dir.readlines()

	for i in range(len(input_list_readlines)):

		read = input_list_readlines[i]
		read = read.replace('\n','')
		condition_file = read
		text_to_dict(condition_file)
