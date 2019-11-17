import networkx as nx
import pr, misc
from info_fns import *
from math import log

debug = False
# pass binary gene data where data[cells][genes] and list of corresp gene names
ROUND = 8

def build(data, gene_names, G=None, targets=None, thresh_mult=1):

	if G==None:
		G = nx.empty_graph(create_using=nx.DiGraph())
		for i in range(len(gene_names)):
			G.add_node(i, gene=gene_names[i], values=data.T[i]) #node number i is its index in data[cell][node]
			# curr values are not used (seperate data array instead)
	pr_dict_xy = pr.build_pr_dict(data, G.nodes(), G.nodes())

	if targets is not None:
		for node in G.nodes():
			if node not in targets: # no in edges from targets
				for node2 in targets:
					G.add_edge(node,node2) # only edges from non-targets to targets

		for node in targets:
			for node2 in G.nodes():
				if node2 not in targets:
					G.add_edge(node,node2) #and edges from targets to non-targets

	else:
		for node in G.nodes():
			for node2 in G.nodes():
				if node is not node2:
					G.add_edge(node,node2) # fully connected graph, no self loops

	#G = greedy_trim(G, data)
	#print("\nAFTER GREEDY TRIM:")
	#misc.postprocess(G)

	thresh = generate_thresholds(G,thresh_mult,trimmed=True)

	g=0
	if targets==None: targets = list(G.nodes())

	for node in G.nodes():
		if debug: print('\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
		print(node,"Determining edges for target ",G.nodes[node]['gene'])
		G=determine_in_edges(G, data, thresh, node, gene_names) #don't think need to return G anyhow
		g+=1

	return G


def greedy_trim(G, data):
	# NOT USED 
	#preliminary phase to reduce to log2(n) nodes where n is amount of data
	n=len(G.nodes[0]['values'])
	for node in G.nodes():
		print("Trimming for node", node)
		while(len(G.in_edges(node)) > log(n,2)+2):
			in_nodes = [edge[0] for edge in G.in_edges(node)]
			pr_dict_x = pr.build_sub_pr_dict(pr_dict_xy, node)
			chng_H_min, minNode = 10000, None
			for edge in G.in_edges(node):
				excluded_source_node=edge[0]
				sub_pr_dict_xy = pr.build_sub_pr_dict(pr_dict_xy, excluded_source_node)
				sub_pr_dict_x = pr.build_sub_pr_dict(pr_dict_x, excluded_source_node) 
				change_in_cond_entropy = H(sub_pr_dict_xy) - H(sub_pr_dict_x) - H(pr_dict_xy) + H(pr_dict_x)
				change_in_cond_entropy = round(change_in_cond_entropy,ROUND)

				#print(change_in_cond_entropy,'\t\t',H(sub_pr_dict_xy), H(sub_pr_dict_x), H(pr_dict_xy), H(pr_dict_x))
				if change_in_cond_entropy < chng_H_min:
					chng_H_min=change_in_cond_entropy
					minNode=excluded_source_node
					pr_dict_xy = sub_pr_dict_xy
					pr_dict_x = sub_pr_dict_x
			assert(minNode is not None)
			#print(chng_H_min)
			G.remove_edge(minNode,node)
	return G


def determine_in_edges(G, data, thresh, target_node, gene_names, pr_dict_xy=None, pr_dict_x=None, pr_dict_y=None):
	# pr_dict is used to evaluate entropy with all edges
	# sub_pr_dict is same but one input edge excluded
	# need H(x,y) and H(x), hence pr_dict_xy and _x 
	# thresh is an array such that thresh[i]=threshold for set of size i

	if pr_dict_xy == None: #then just started with this node, build a complete dict
		assert(pr_dict_x == None)
		source_nodes = []
		for edge in G.in_edges(target_node): #def a more concise way
			node = edge[0]
			assert(node is not target_node) 
			source_nodes += [node]
		pr_dict_xy = pr.build_pr_dict(data, G.nodes(), [target_node]+source_nodes)
		pr_dict_x = pr.build_sub_pr_dict(pr_dict_xy, target_node)
		pr_dict_y = pr.build_pr_dict(data, G.nodes(), [target_node])

	# TEST REVERSE ORDER:
	reverse=False
	edges = [] 
	for edge in G.in_edges(target_node):
		edges += [edge]
	for e in range(len(edges)):
		if reverse: edge = edges[-e]
		else: edge=edges[e]
	#for edge in G.in_edges(target_node):
		excluded_source_node = edge[0]
		k = len(G.in_edges(target_node)) #then use E(change k+1 -> k since include Y)
		sub_pr_dict_xy = pr.build_sub_pr_dict(pr_dict_xy, excluded_source_node)
		sub_pr_dict_x = pr.build_sub_pr_dict(pr_dict_x, excluded_source_node)  

		if debug: 
			print("\nTesting source: " + str(gene_names[excluded_source_node]))

		rm_edge = rm_evaluation(thresh, k, pr_dict_xy, pr_dict_x, pr_dict_y, sub_pr_dict_xy, sub_pr_dict_x) 

		if rm_edge:
			G.remove_edge(excluded_source_node, target_node)
			if len(G.in_edges(target_node)) > 0:

				#copies may not be nec, just being careful
				sub_pr_dict_xy_copy,sub_pr_dict_x_copy, pr_dict_y_copy = sub_pr_dict_xy.copy(), sub_pr_dict_x.copy(), pr_dict_y.copy(),
				determine_in_edges(G, data, thresh, target_node, gene_names, pr_dict_xy=sub_pr_dict_xy_copy, pr_dict_x=sub_pr_dict_x_copy, pr_dict_y=pr_dict_y_copy) #run algo again with remaining subset
			return G # other edges have been recursively checked
	
		# else keep trying the other edges

	return G




def generate_thresholds(G, mult,trimmed=False):
 	# mult = how many * the rd distrib should change in entropy be < to be significant?
	# tresholds for rm'g edges, one per size of node (#in_edges) 
	# thresh = 1 -> no edges remain; thresh = 0 -> all edges remain
	n = len(G.nodes[0]['values'])
	if trimmed: m=min(len(G.nodes()),int(log(n,2))+4)
	else: m=len(G.nodes())

	#print("\nEmpirically calculating Ee with m,n=",m,n)
	finite_m=False
	if finite_m:
		Ee = []
		print('m,n=',m,n)
		for i in range(m+1):
			Ee += [expected_entropy(n,i)]
			print(i,Ee)
		#Ee = Ee_empirical(n=n+1, mMax=m+1, reps=100)

		thresh = [0] + [2*Ee[k]-Ee[k-1]-Ee[k+1] for k in range(1,m)]

	else:
	# if m->inf ?
		#thresh = [0 for i in range(len(G.nodes()))]
		thresh = [0,1]+[-1*log((pow(k,2)-1)/pow(k,2),2) for k in range(2,len(G.nodes()))]
	
	#thresh = [0 for k in range(len(G.nodes))]
	for t in range(len(thresh)): thresh[t] *= mult

	#print("\nThresh pre max(.00001) = ", thresh)
	#thresh = [max(round(thresh[t],ROUND),.005) for t in range(len(thresh))]
	thresh = [round(thresh[t],ROUND) for t in range(len(thresh))]
	#if debug:
	print("\nUsing thresh = ", thresh)
	return thresh



def rm_evaluation(thresh, k, pr_dict_xy, pr_dict_x, pr_dict_y, sub_pr_dict_xy, sub_pr_dict_x, metric='conditional'):
	# k = starting # in edges, eval k -> k-1 in edges some sets of size k+1 -> k

	if metric=='conditional' or metric =='info':
		# rm edge if H(Y|{X}\x_j) - H(Y|{X}) > thresh[k]
		change_in_cond_entropy = H(sub_pr_dict_xy) - H(sub_pr_dict_x) - H(pr_dict_xy) + H(pr_dict_x)
		if metric == 'info': change_in_cond_entropy /= H(pr_dict_y)
		change_in_cond_entropy = round(change_in_cond_entropy,ROUND)
		assert(change_in_cond_entropy >= 0 and change_in_cond_entropy <= 1)

		if debug:
			#print('\npr_dict_xy',pr_dict_xy)
			#print('pr_dict_x',pr_dict_x)
			#print('H(XY),H(X),H(xY),H(x)=',H(pr_dict_xy), H(pr_dict_x), H(sub_pr_dict_xy), H(sub_pr_dict_x))
			#print("H(Y|X) = ",H(pr_dict_xy) - H(pr_dict_x))
			#print("H(Y|X -x) = ",H(sub_pr_dict_xy) - H(sub_pr_dict_x) )
			print("Change H(Y|X) = " + str(change_in_cond_entropy))
			print("vs thresh = ",thresh[k])

		if change_in_cond_entropy > thresh[k]: #thresh[k]: 
			if debug: print("Edge kept.")
			return False #keep it!
		else: 
			if debug: print("Edge rm'd")
			return True #rm it!

	else: assert(False) #unknown metric



############ TESTING ################

if __name__ == "__main__":
	print("TESTING thresh:")
	thresh = [0,1]+[-1*log((pow(k,2)-1)/pow(k,2),2) for k in range(2,20)]
	for t in thresh:
		print(t)