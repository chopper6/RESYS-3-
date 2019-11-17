import numpy as np, os, random as rd

debug = True

# TODO: import the file and use as local path
input_file = 'C:/Users/Crbn/Documents/MPRI M2/ReSys/project/data/wholecells_binary.csv'
#input_file = 'C:/Users/Cs/Documents/MPRI M2/ReSys/project/wholecells_binary.csv'
lizard_file = 'C:/Users/Crbn/Documents/MPRI M2/ReSys/project/data/lizards_bin.csv'
coronary_file = 'C:/Users/Crbn/Documents/MPRI M2/ReSys/project/data/coronary_bin.csv'
asia_file = 'C:/Users/Crbn/Documents/MPRI M2/ReSys/project/data/asia_bin.csv'


def gen_data(set='hema',include_stages=False, cutoff=None, hema_file=input_file): #'3NAND_AND_2OR'
	# data should be a matrix: data[cells][genes]
	if set == 'hema':
		data, genes = import_data(hema_file, cutoff,include_stages=include_stages)
		return data,genes

	elif set == 'lizards':
		data, genes = import_benchmark_data(lizard_file, cutoff,include_stages=False)
		return data,genes
	elif set == 'asia':
		data, genes = import_benchmark_data(asia_file, cutoff,include_stages=False)
		return data,genes
	elif set == 'coronary':
		data, genes = import_benchmark_data(coronary_file, cutoff,include_stages=False)
		return data,genes

	elif set == 'AND':
		AND = [
		[1,1,0,0],
		[1,0,1,0],
		[1,0,0,0]]

		genes = ['x1','x2','y'] #for import this would just be the first row
		data=np.array(AND).T
		return data, genes

	elif set == 'XOR':
		XOR = [
		[1,1,0,0],
		[1,0,1,0],
		[0,1,1,0]]

		genes = ['x1','x2','y'] 
		data=np.array(XOR).T
		return data, genes

	elif set == '3OR':
		data = [
		[1,1,1,1,0,0,0,0],
		[1,1,0,0,1,1,0,0],
		[1,0,1,0,1,0,1,0],
		[1,1,1,1,1,1,1,0]]

		genes = ['x1','x2','x3','y'] 
		data=np.array(data).T
		return data, genes

	elif set == '4AND':
		data = [
		[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
		[1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0],
		[1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0],
		[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],
		[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

		genes = ['x1','x2','x3','x4','y'] 
		data=np.array(data).T
		return data, genes

	elif set == '5AND':
		data = [
		[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
		[1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0],
		[1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0],
		[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],
		[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

		genes = ['x1','x2','x3','x4','x5','y'] 
		data=np.array(data).T
		return data, genes

	elif set == '5ANDx4':
		data = [
		[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
		[1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0],
		[1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0],
		[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],
		[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

		genes = ['x1','x2','x3','x4','x5','y'] 
		data=np.array(data).T
		return data, genes

	elif set == '5ANDx4noisy':
		data = [
		[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
		[1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0],
		[1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0],
		[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],
		[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]

		for i in range(len(data)):
			for j in range(len(data[0])):
				if rd.random()>.9: data[i][j]=(data[i][j]+1) %2
		genes = ['x1','x2','x3','x4','x5','y'] 
		data=np.array(data).T
		return data, genes

	elif set == '2OR_1USELESS':
		data = [
		[1,1,1,1,0,0,0,0],
		[1,1,0,0,1,1,0,0],
		[1,0,1,0,1,0,1,0],
		[1,1,1,1,1,1,0,0]]

		genes = ['x1','x2','u1','y'] 
		data=np.array(data).T
		return data, genes

	elif set == '3empty':
		data = [
		[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
		[1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0],
		[1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0],
		[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],
		]
		genes = ['x1','x2','x3'] 
		data=np.array(data).T
		return data, genes

	elif set == '4empty':
		data = [
		[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
		[1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0],
		[1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0],
		[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],
		]
		genes = ['x1','x2','x3','x4'] 
		data=np.array(data).T
		return data, genes

	elif set == 'redundant':
		data = [
		[1,1,1,0],
		[1,1,1,0],
		[1,1,1,0],
		[1,1,1,0]]
		genes = ['x1','x2','x3','x4'] 
		data=np.array(data).T
		return data, genes

	elif set == '3NAND_AND_2OR':
		data = [
		[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
		[1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0],
		[1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0],
		[1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0],
		[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0],
		[0,0,0,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0]]

		genes = ['x1','x2','x3','x4','x5','y'] 
		data=np.array(data).T
		return data, genes

	else: assert(False) #unknown data set






def import_benchmark_data(input_file, cutoff, include_stages=False):
	# WARNING: specific to this project (adds cols for cell states)
	assert(cutoff is None) #not implemented
	assert(os.path.isfile(input_file)) #check input_file path at top of data.py

	with open(input_file, 'r') as input:
		lines = input.readlines()

		gene_titles = lines[0].split(",")
		titles = lines[0].split(",")
		piece = titles[-1].split("\n")
		titles[-1] = piece[0]
		num_features = len(titles)
		num_instances = len(lines)-1    # rows, not counting title row
		data = np.empty((num_instances, num_features),dtype=int) 
		# will transform first col into 5 possible cell stages 

		for i in range(0,num_instances):
			row_orig = lines[i+1].split(",", num_features+1)
			row = [int(row_orig[r]) for r in range(len(row_orig))]
			#piece = row[-1].split("\n")
			#row[-1] = piece[0]
			data[i] = row
	#print_data(data, titles)

	return data, titles


def import_data(input_file, cutoff, include_stages=False):
	# WARNING: specific to this project (adds cols for cell states)
	if cutoff is not None: print("\nWARNING: cutting off genes > " + str(cutoff))
	assert(os.path.isfile(input_file)) #check input_file path at top of data.py

	with open(input_file, 'r') as input:
		lines = input.readlines()

		# add 'genes' to represent cells states
		titles = ['PS','NP','HF','4SG','4SFG']
		gene_titles = lines[0].split(",")
		titles += gene_titles[1:] #first title is just "samples"
		piece = titles[-1].split("\n")
		titles[-1] = piece[0]
		orig_num_features = len(titles)
		if cutoff is not None: titles = titles[:cutoff+4]

		if cutoff is not None: num_features = min(len(titles)-4,cutoff)     # columns
		else: num_features = len(titles)-4
		num_instances = len(lines)-1    # rows, not counting title row
		data = np.empty((num_instances, num_features+4)) 
		# will transform first col into 5 possible cell stages 

		for i in range(0,num_instances):
			row_orig = lines[i+1].split(",", orig_num_features)
			row_orig = [row_orig[0]] + [int(row_orig[r]) for r in range(1,len(row_orig))]
			# messy but hey, r[0] is str, others should be converted to int

			row = [0,0,0,0] + row_orig[:num_features]
			#piece = row[-1].split("\n")
			#row[-1] = piece[0]

			# add 'genes' representing spc cell states
			pieces = row[4].split('_')
			cell_stage = pieces[0]#[:-2] #last two corresp to organism letter and sample # jP
			#print("cell stage = " + str(cell_stage))
			if 'PS' in cell_stage:
				row[0], row[4] = 1,0
			elif 'NP' in cell_stage:
				row[1], row[4] = 1,0
			elif 'HF' in cell_stage:
				row[2], row[4] = 1,0
			elif '4SG' in cell_stage:
				row[3], row[4] = 1,0
			elif '4SFG' in cell_stage:
				row[4] = 1
			else: assert(False) #unknown cell stage

			data[i] = row

	if not include_stages:
		titles = titles[5:]
		data = data[:,5:]
	#if debug: print_data(data, titles)
	return data, titles


def print_data(data, titles):
	# just for debugging
	print("Imported titles = " + str(titles))
	print("First row of data = " + str(data[0]))
	print("Last row of data = " + str(data[-1,:]))
	# data[-1] would also work, just showing similarity to below notation

	# first 5 reversed for cell types, each cell should be only one of them
	if False:
		for row in data:
			assert(row[0]==1 or row[1]==1 or row[2]==1 or row[3]==1 or row[4]==1)
			for i in range(5):
				for j in range(5):
					if i!=j: assert(not(row[i]==1 and row[j]==1))

	print("First column of data = " + str(data[:,0]))
	print("Last column of data = " + str(data[:,-1]))
	# this is called slicing and can only be done with numpy lists, data[:] = data
