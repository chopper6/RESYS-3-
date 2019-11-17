import data, pr, info_fns, G_builder, misc
import networkx as nx, numpy as np
from matplotlib import pyplot as plt


def drawG(G, output_path):
	# draw network, change to diff file or fn later
	plt.figure(figsize=(20,20))
	node_size, nalpha, ealpha = 800, .3, .4
	pos = nx.circular_layout(G)  # positions for all nodes
	elist = sorted(list(G.edges()))
	nx.draw_networkx_edges(G, pos, arrows=True, edgelist=elist, alpha=ealpha)

	nodes = list(G.nodes())

	colors = [G.nodes[nodes[n]]['color'] for n in range(len(nodes))]
	#print('colors in misc.drawG():',colors)
	nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=colors, node_size=node_size, alpha=nalpha)
	labels = {n: G.nodes[n]['gene'] for n in G.nodes()}
	nx.draw_networkx_labels(G, pos, labels, font_size=8, font_color='blue')

	#nx.draw(G)
	#plt.show()
	plt.savefig(output_path + 'inferred_G.png')


def preprocess(data,gene_names):
	#TODO: other forms of normalization for entropy of indiv genes? poss in thresh?
	# just see if v low entropy genes (HoxD8 + HoxB2 are overrepresented in solution)

	print('\nBefore Preproc have',len(gene_names),'nodes:\n')
	Hmin = .0001
	print("\nEntropy of individual genes: ")
	all_nodes = [g for g in range(len(gene_names))]
	rm = []
	for g in range(len(gene_names)):
		pr_dict_x = pr.build_pr_dict(data, all_nodes, [g])
		Hx = info_fns.H(pr_dict_x)
		print(gene_names[g], Hx)
		if Hx < Hmin: rm += [g]
	i=0
	for r in rm:
		print("removing", gene_names[r-i],', since H<',Hmin)
		gene_names.remove(gene_names[r-i])
		np.delete(data, r-i, axis=1)

		i+=1

	print('\nAfter Preproc have',len(gene_names),'nodes:\n', gene_names)
	return data,gene_names


def postprocess(G):
	# write some output with hub nodes (in-deg, out-deg, centrality)
	rm = []
	for node in G.nodes():
		if len(G.in_edges(node))+len(G.out_edges(node))==0:
			rm += [node]

	for r in rm:
		G.remove_node(r)

	#sort by degs
	in_deg_nodes = list(G.nodes())
	in_degs = [G.in_degree(n) for n in in_deg_nodes]
	in_deg_nodes = [x for _,x in sorted(zip(in_degs,in_deg_nodes))] 
	print("\nGene : in-degree")
	for node in in_deg_nodes:
		print(G.nodes[node]['gene'],'\t',G.in_degree(node))

	out_deg_nodes = list(G.nodes())
	out_degs = [G.out_degree(n) for n in out_deg_nodes]
	out_deg_nodes = [x for _,x in sorted(zip(out_degs,out_deg_nodes))] 
	print("\nGene : out-degree")
	for node in out_deg_nodes:
		print(G.nodes[node]['gene'],'\t',G.out_degree(node))

	btwn_nodes = list(G.nodes())
	btwn_centrality_dict = nx.betweenness_centrality(G)
	btwn_centrality = [btwn_centrality_dict[node] for node in btwn_nodes]
	btwn_nodes = [x for _,x in sorted(zip(btwn_centrality,btwn_nodes))] 
	print("\nGene : betweenness centrality")
	#print(len(btwn_centrality),len(G.nodes()), len(btwn_nodes))
	for node in btwn_nodes:
		print(G.nodes[node]['gene'], btwn_centrality_dict[node])

	print('\nfinished with',len(G.nodes()),'nodes:\n',[G.nodes[node]['gene'] for node in G.nodes()])
	if len(in_degs)>0: print('\nAverage in degree = ',sum(in_degs)/len(in_degs))
	if len(out_degs)>0:print('\nAverage out degree = ',sum(out_degs)/len(out_degs))
	#print("\nfinished with edges: ")
	#for edge in G.edges():
	#	print(G.nodes[edge[0]]['gene'],' --> ', G.nodes[edge[1]]['gene'])

	return G

def assign_stages(G, rm_stages=True):
	#poss need to normalize colors in [0,1]
	stages = ['PS','NP','HF','4SG','4SFG']
	#stages_of_interest = ['PS','4SG','4SFG']
	#colors = [[0,0,255],[255,0,0],[0,255,0]]
	stages_of_interest = ['PS','4SG','4SFG']
	colors = [[0,0,255],[255,0,0],[0,255,0]]
	#colors = [[0,0,255],[128,0,255],[255,128,0],[0,255,0],[255,0,0]] #RGB values one for each of these
	for node in G.nodes():
		if G.nodes[node]['gene'] not in stages:
			num,rgb=0,[0,0,0]
			in_edge_check, out_edge_check = False, True
			if in_edge_check == True:
				for in_edge in G.in_edges(node):
					in_node = in_edge[0]
					for s in range(len(stages_of_interest)):
						if G.nodes[in_node]['gene'] == stages_of_interest[s]:
							for c in range(len(rgb)):
								rgb[c] += colors[s][c]
							num+=1
			if out_edge_check == True:
				for out_edge in G.out_edges(node):
					out_node = out_edge[1]
					for s in range(len(stages_of_interest)):
						if G.nodes[out_node]['gene'] == stages_of_interest[s]:
							for c in range(len(rgb)):
								rgb[c] += colors[s][c]
							num+=1
			if num>1:
				for c in range(len(rgb)):
					rgb[c] /= num

			for c in range(len(rgb)):
				rgb[c] /= float(255)
			G.nodes[node]['color'] = rgb
			if rgb == [0,0,0]: G.nodes[node]['color'] = [160/255,160/255,160/255]

	# rm stages nodes
	for node in G.nodes():
		if G.nodes[node]['gene'] in stages:
			G.nodes[node]['color'] = [1,0,0]
	if rm_stages:
		rm = []
		for node in G.nodes():
			if G.nodes[node]['gene'] in stages:
				rm += [node]

		for r in rm:
			G.remove_node(r)	

	else:print("Warning: not rm'g stages")

	return G