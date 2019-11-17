# build a dict with keys b1,b2,...,bk 
# later calls revise dict without recalc by concat certain strings
# target and source genes should really be their indices

def build_pr_dict(data, all_nodes, included_nodes):
	# data should be a matrix: data[cells][genes]
	# returns dict of pr(target_gene, source_genes) and pr(source_genes)
	# ORDER OF SOURCE GENES MATTERS SINCE INDEXED IN ORDER
	# poss add pr(target_gene), but I think maxzn is invariant to it (and would form symm )

	pr_dict = {}
	num_cells = len(data)
	for i in range(num_cells):
		key = ''

		for node in all_nodes:
			if node in included_nodes:
				if data[i][node] == 1:
					key+='1'
				elif data[i][node] == 0: key+='0'
				else: 
					assert(False) #only accepting binary atm
			else: key+='0'

		if key not in pr_dict.keys(): 
			#poss faster way? such as init all poss combos = 0
			#or an array struct...but i think dict concat easier
			pr_dict[key] = 1
		else:
			pr_dict[key] += 1

	for k in pr_dict.keys():
		pr_dict[k] /= float(num_cells) 

	return pr_dict



def build_sub_pr_dict(pr_dict, excluded_source_gene):
	sub_pr_dict = {}
	for k in pr_dict.keys():
		new_key = k[:excluded_source_gene] + str(0) + k[excluded_source_gene+1:]
		#want to keep same str length, but now excluded gene has p(0)=1

		if new_key not in sub_pr_dict:
			sub_pr_dict[new_key] = pr_dict[k]

		else:
			sub_pr_dict[new_key] += pr_dict[k]

	return sub_pr_dict





