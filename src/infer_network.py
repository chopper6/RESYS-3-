import data, pr, info_fns, G_builder, misc

########## LOCAL PATHS ############

output_path = 'C:/Users/Crbn/Documents/MPRI M2/ReSys/project/output/'
hema_file = 'C:/Users/Crbn/Documents/MPRI M2/ReSys/project/data/wholecells_binary.csv'

########## MAIN PARAMETERS ###########

dataset='hema' #can also try others in data.py, such as '3NAND_AND_2OR'
cutoff=None
G=None
thresh_mult = 1

############# MAIN ################

dataa, gene_names = data.gen_data(set=dataset, include_stages=False, cutoff=cutoff, hema_file=hema_file)
dataa,gene_names = misc.preprocess(dataa, gene_names)
G = G_builder.build(dataa, gene_names, G=G, thresh_mult=thresh_mult)
G = misc.postprocess(G)
G = misc.assign_stages(G)

#print(G.nodes[edge[0]]['gene'], '->', G.nodes[edge[1]]['gene'])

misc.drawG(G,output_path) 