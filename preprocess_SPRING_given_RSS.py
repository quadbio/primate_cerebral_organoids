#!/usr/bin/env python

#import sompy
import sys
import numpy as np
import pandas as pd
import argparse
import os
import json
import matplotlib.pyplot as plt
#from preprocessing_python import *

DEFAULT = {
	'grouping' : None,
	'coarse_grain' : 1,
	'cg_sep' : None,
	'cg_method' : "nntplus",
	'cg_idx': "pseudo-idx.txt",
	'cg_input' : "pseudo-input.txt",
	'cg_grouping' : "pseudo-grouping.txt",
	'cg_dist' : None,
	'cg_thres' : 0,
	'dist_input': None,
	'dist_output' : None,
	'distMethod' : "cos0",
	'dist_weight' : "1",
	'numpc' : 50,
	'k' : 20,
	'num_fd_iter' : 100,
	'output' : 'SPRING'
}

# main function
#========================================================================#
def main():
	parser = argparse.ArgumentParser()
	
	# general input
	parser.add_argument("input", help="Path to the expression matrix CSV file (with comma as delimiters, no header; columns as samples). Use comma to separate file names if multiple input matrix are used.")
	parser.add_argument("-g", "--grouping", type = str, default = DEFAULT['grouping'], help="Path to the cell grouping CSV file (with comma as delimiters, no header; columns as samples).")
	
	# group input given grouping index
	parser.add_argument("--do-grouping-only", dest = "do_grouping", action = "store_true", help="Only do one thing, which is to group the input matrix into average input matrix of pseudo-cells given the pseudo-cell indexing table.")
	
	# coarse grain related
	parser.add_argument("-c", "--coarse-grain", dest = "coarse_grain", type = float, default = DEFAULT['coarse_grain'], help="Coarse-grain parameter (smaller than 1 for coarse grain, default: 1).")
	parser.add_argument("--cg-index-only", dest = "cg_idx_only", action = "store_true", help="Only apply coarse-grain and output the indices table of pseudo-cells with their grouping information.")
	parser.add_argument("--cg-sep", dest = "cg_sep", type = str, default = DEFAULT['cg_sep'], help="Apply coarse-grain procedure to cells in different groups separately, based on this grouping label involved in GROUPING (default: None)")
	parser.add_argument("--cg-method", dest = "cg_method", type = str, default = DEFAULT['cg_method'], help="The coarse grain method to generate pseudocells (default: nntplus [nntplus/downsample/knn/simple/nnterritory])")
	parser.add_argument("--cg-index", dest = "cg_idx", type = str, default = DEFAULT['cg_idx'], help="If coarse grain is applied, output the pseudo-cell aggregation index to this file (default: pseudo-idx.txt)")
	parser.add_argument("--cg-input", dest = "cg_input", type = str, default = DEFAULT['cg_input'], help="The output pseudo-cell expression matrix file (default: pseudo-input.txt)")
	parser.add_argument("--cg-grouping", dest = "cg_grouping", type = str, default = DEFAULT['cg_grouping'], help="The output pseudo-cell grouping labels file (default: pseudo-grouping.txt)")
	parser.add_argument("--cg-dist", dest = "cg_dist", type = str, default = DEFAULT['cg_dist'], help="The output pseudo-cell distance matrix file (default: None for no output)")
	parser.add_argument("--cg-thres", dest = "cg_thres", type = int, default = DEFAULT['cg_thres'], help="The threshold of cell number represented by each pseudo-cell after applying coarse gain (default: 3).")
	
	# distance matrix
	parser.add_argument("--dist-input", dest = "dist_input", type = str, default = DEFAULT['dist_input'], help="Path to the distance matrix (optional, with comma as delimiters, no header nor index)")
	parser.add_argument("--dist-output", dest = "dist_output", type = str, default = DEFAULT['dist_output'], help = "Output the distance matrix to this file (default: None for no output)")
	parser.add_argument("--dist-method", dest = "dist_method", type = str, default = DEFAULT['distMethod'], help = "The method used to calculate distance (cos0 / cos / pca / cor / eucl)")
	parser.add_argument("--dist-weight", dest = "dist_weight", type = str, default = DEFAULT['dist_weight'], help = "The weights to merge multiple distance matrix when multiple input matrices are given (default: 1). Use comma to separate different weights.")
	parser.add_argument("--numpc", type = int, default = DEFAULT['numpc'], help = "If 'pca' is used as the distance method, it determines the number of PCs used (default: 50)")
	
	# k & mutual neighbors
	parser.add_argument("-k", type = int, default=DEFAULT['k'], help="The numbers of nearest neighbors (default: 20).")
	parser.add_argument("--only-mutual", dest = "only_mutual", action = "store_true", help="Only connect mutual nearest neighbors")
	
	# SPRING version, parameter and output
	parser.add_argument("--use-old-spring", dest = "use_old_spring", action = "store_true", help="Instead of using ForceAtlas2 to generate kNN-network layout (SPRING v1.6), make the SPRING folder for SPRING Javascript running (SPRING v1.5)")
	parser.add_argument("--num-fd-iter", dest = "num_fd_iter", type = int, default = DEFAULT['num_fd_iter'], help = "The number of iterations to generate ForceAtlas2 force-directed kNN-network layout (only for SPRING v1.6).")
	parser.add_argument("--output", type = str, default = DEFAULT['output'], help="Output path of the SPRING folder for SPRING v1.5, or file prefix of the SPRING coordinates (SPRING v1.6) (default: SPRING).")
	
	args = parser.parse_args()
	
	# only do grouping if required
	if args.do_grouping:
		if args.cg_idx is None:
			sys.exit("CRITICAL ERROR: the pseudo-cell grouping indexing file is not given!")
		idx = pd.read_csv(args.cg_idx)
		print("PROGRESS: done reading the pseudo-cell grouping indices")
		input = pd.read_csv(args.input, header=None, index_col=0)
		print("PROGRESS: done reading the input file")
		inputPooled = groupInput(idx, input)
		print("PROGRESS: outputing the grouped input")
		inputPooled.to_csv(args.output + ".csv", header=False, index=True)
		print("DONE grouping input!")
		sys.exit()
	
	# read input and grouping labels
	files = args.input.split(",")
	input = [ pd.read_csv(file, index_col=0, header=None) for file in files ]
	inputMerged = None
	weights = np.repeat(1, len(input))
	if len(args.dist_weight.split(",")) > 1:
		weights = np.array([ float(x) for x in args.dist_weight.split(",") ])
	print("PROGRESS: done reading inputs and weights")
	
	cell_groupings_dict = dict()
	if args.grouping is not None:
		cell_groupings = pd.read_csv(args.grouping, header = None, index_col = 0)
		for i in range(cell_groupings.shape[0]):
			cell_groupings_dict[cell_groupings.index[i]] = list(cell_groupings.iloc[i,])
	print("PROGRESS: done reading grouping labels")
	
	# calculate distances
	dist = None
	if args.dist_input is not None:
		dist = pd.read_csv(args.dist_input, header = None, index_col = None)
		print("PROGRESS: done reading distance matrix")
	
	# do coarse grain to generate pseudo-cells
	if args.coarse_grain < 1:
		cellGroups = np.ones(input[0].shape[1])
		if args.cg_sep is not None and args.cg_sep in cell_groupings_dict:
			cellGroups = np.array(cell_groupings_dict[args.cg_sep])
		
		pseudoIdx = pd.DataFrame({ 'initial' : range(len(cellGroups)), 'pooled' : np.repeat(-1, len(cellGroups))})
		for group in np.unique(cellGroups):
			print("PROGRESS: start coarse grain to generate pseudo-cells for group-" + str(group))
			idxGroup = np.where(cellGroups == group)[0]
			
			inputGroup, distGroup = None, None
			if dist is not None:
				distGroup = dist.iloc[idxGroup, idxGroup]
			elif (len(idxGroup) > 1 and len(idxGroup) < 10000) or (len(input) == 1 and len(idxGroup) < 20000):
				inputGroup = [ x.iloc[:,idxGroup] for x in input ]
				distGroup = calAndIntDist(inputGroup, weights, args.dist_method, args = { 'numpc' : args.numpc })
				print("PROGRESS: got the distance matrix for group-" + str(group))
			elif len(idxGroup) == 1:
				inputGroup = [ x.iloc[:,idxGroup] for x in input ]
				distGroup = pd.DataFrame([[1]])
			else:
				inputGroup = [ x.iloc[:,idxGroup] for x in input ]
			
			pseudoIdxGroup = coarseGrainIdx(inputGroup, weights, distGroup, args.k, args.coarse_grain, method = args.cg_method, dist_args = { 'dist_method' : args.dist_method, 'args' : { 'numpc' : args.numpc } })
			pseudoIdx.iloc[idxGroup,1] = np.array(pseudoIdxGroup.iloc[:,1] + pseudoIdx.iloc[:,1].max() + 1)
			print("PROGRESS: done initial coarse grain for group-" + str(group) + ": " + str(pseudoIdxGroup.shape[0]) + " -> " + str(len(np.unique(pseudoIdxGroup.iloc[:,1]))))
		
		if args.cg_method != "downsample" and args.cg_thres > 1:
			pseudoLabels, pseudoCounts = np.unique(pseudoIdx.iloc[:,1], return_counts=True)
			pseudoPassed = pseudoLabels[pseudoCounts >= args.cg_thres]
			pseudoIdx.loc[np.isin(pseudoIdx.iloc[:,1], pseudoPassed, invert = True), "pooled"] = -1
			print("PROGRESS: done filtering of coarse grain: " + str(len(pseudoLabels)) + " -> " + str(len(pseudoPassed)))
		
		pseudoIdx.to_csv(args.cg_idx, header = True, index = False)
		cell_groupings_dict = groupGroupings(pseudoIdx, cell_groupings_dict)
		pd.DataFrame(cell_groupings_dict).transpose().to_csv(args.cg_grouping, header = False, index = True)
		
		if not args.cg_idx_only:
			input = [ groupInput(pseudoIdx, x) for x in input ]
			inputMerged = input[0]
			if len(input) > 1:
				for i in range(len(input) - 1):
					inputMerged = inputMerged.append(input[i+1])
			inputMerged.to_csv(args.cg_input, header = False, index = True)
			print("PROGRESS: generated grouped input")
			
			dist = calDistMultipleInputs(input, weights, args.dist_method, args = { 'numpc' : args.numpc })
			print("PROGRESS: calculated grouped distance")
			if args.cg_dist is not None:
				dist.to_csv(args.cg_dist, header = False, index = False)
		else:
			print("DONE")
	
	# prepare SPRING input folder
	if not args.cg_idx_only:
		if dist is None:
			dist = calDistMultipleInputs(input, weights, args.dist_method, args = { 'numpc' : args.numpc })
		
		if args.use_old_spring:
			if inputMerged is None:
				inputMerged = input[0]
				if len(input) > 1:
					for i in range(len(input) - 1):
						inputMerged = inputMerged.append(input[i+1])
			gene_list = list(inputMerged.index)
			
			E_np = np.array(inputMerged.transpose())
			D_np = np.array(dist)
			
			print("PROGRESS: start outputing")
			saveSpringDir(E_np, D_np, args.k, gene_list, args.output, cell_groupings = cell_groupings_dict, only_mutual_neighbors = args.only_mutual)
			print("DONE")
		else:
			links = None
			if args.only_mutual:
				links = getMutualEdges(getKNNFromDist(dist, args.k))
			else:
				links = getEdges(getKNNFromDist(dist, args.k))
			print("PROGRESS: got edges, start layout calculation")
			coord = get_force_layout(list(links), dist.shape[0], n_iter=args.num_fd_iter, edgeWeightInfluence=1, barnesHutTheta=2, scalingRatio=1, gravity=0.05, jitterTolerance=1, verbose=False)
			coord = coord / 5.0
			coord = coord - np.min(coord, axis = 0) - np.ptp(coord, axis = 0) / 2.0
			coord[:,0] = coord[:,0]  + 750
			coord[:,1] = coord[:,1]  + 250
			print("PROGRESS: got layout, outputing")
			pd.DataFrame(coord).to_csv(args.output + ".txt", index = True, header = False)
			print("DONE")

# for SPRING v1.5 - to generate the folder for the javascript run
# modified from save_spring_dir in preprocessing_python.py
#========================================================================#
def saveSpringDir(E, D, k, gene_list, project_directory, cell_groupings={}, custom_colors={}, use_genes=[], only_mutual_neighbors=False):
	if not os.path.exists(project_directory):
		os.makedirs(project_directory)
	if not project_directory[-1] == '/': project_directory += '/'
	
	# Build graph
	print('Building graph')
	if only_mutual_neighbors:
		edges = list(getMutualEdges(getKNNFromDist(pd.DataFrame(D), k)))
	else:
		edges = list(getEdges(getKNNFromDist(pd.DataFrame(D), k)))
		print('  Built a graph with ' + str(len(edges)) + " edges")
		# edges = get_knn_edges(D, k)
	
	# save genesets
	print('Saving gene sets')
	custom_colors['Uniform'] = np.zeros(E.shape[0])
	write_color_tracks(custom_colors, project_directory+'color_data_gene_sets.csv')
	all = []
	
	# save gene colortracks
	print('Saving coloring tracks')
	if not os.path.exists(project_directory+'gene_colors'):
		os.makedirs(project_directory+'gene_colors')
	II = len(gene_list) / 50 + 1
	for j in range(50):	
		fname = project_directory+'gene_colors/color_data_all_genes-'+repr(j)+'.csv'
		if len(use_genes) > 0: all_gene_colors = {g : E[:,i+II*j] for i,g in enumerate(gene_list[II*j:II*(j+1)]) if g in use_genes}
		else: all_gene_colors = {g : E[:,i+II*j] for i,g in enumerate(gene_list[II*j:II*(j+1)])}
		write_color_tracks(all_gene_colors, fname)
		all += all_gene_colors.keys()
	
	# Create and save a dictionary of color profiles to be used by the visualizer
	print('Color stats')
	color_stats = {}
	for i in range(E.shape[1]):
		mean = float(np.mean(E[:,i]))
		std = float(np.std(E[:,i]))
		max = float(np.max(E[:,i]))
		centile = float(np.percentile(E[:,i],99.6))
		color_stats[gene_list[i]] = (mean,std,0,max,centile)
	for k,v in custom_colors.items():
		color_stats[k] = (0,1,np.min(v),np.max(v)+.01,np.percentile(v,99))
	json.dump(color_stats,open(project_directory+'/color_stats.json','w'),indent=4, sort_keys=True)
	
	
	# save cell labels
	print('Saving categorical color data')
	categorical_coloring_data = {}
	for k,labels in cell_groupings.items():
		labels = [(l if type(l)==str else repr(l)) for l in labels]
		label_colors = {l:frac_to_hex(float(i)/len(set(labels))) for i,l in enumerate(list(set(labels)))}
		categorical_coloring_data[k] = {'label_colors':label_colors, 'label_list':labels}
	json.dump(categorical_coloring_data,open(project_directory+'categorical_coloring_data.json','w'),indent=4)
	
	
	print('Writing graph')
	nodes = [{'name':i,'number':i} for i in range(E.shape[0])]
	edges = [{'source':i, 'target':j, 'distance':0} for i,j in edges]
	out = {'nodes':nodes,'links':edges}
	open(project_directory+'graph_data.json','w').write(json.dumps(out,indent=4, separators=(',', ': ')))

# from processing_python.py
#========================================================================#
def write_color_tracks(ctracks, fname):
	out = []
	for name,score in ctracks.items():
		line = ','.join([name]+[repr(round(x,1)) for x in score])
		out += [line]
	out = sorted(out,key=lambda x: x.split(',')[0])
	open(fname,'w').write('\n'.join(out))

# from processing_python.py
#========================================================================#
def frac_to_hex(frac):
	rgb = tuple(np.array(np.array(plt.cm.jet(frac)[:3])*255,dtype=int))
	return '#%02x%02x%02x' % rgb

# from processing_python.py
#========================================================================#
def get_knn_edges(dmat, k):
	edge_dict = {}
	for i in range(dmat.shape[0]):
		for j in np.nonzero(dmat[i,:] <= sorted(dmat[i,:])[k])[0]:
			if i != j:
				ii,jj = tuple(sorted([i,j]))
				edge_dict[(ii,jj)] = dmat[i,j]
	return edge_dict.keys()

# for SPRING v1.6 - to use ForceAtlas2 to generate force-directed layout
# from spring_helper.py
#========================================================================#
def get_force_layout(links, n_cells, n_iter=100, edgeWeightInfluence=1, barnesHutTheta=2, scalingRatio=1, gravity=0.05, jitterTolerance=1, verbose=False):
	from fa2 import ForceAtlas2
	import networkx as nx
	
	G = nx.Graph()
	G.add_nodes_from(range(n_cells))
	G.add_edges_from(list(links))
	
	forceatlas2 = ForceAtlas2(
		# Behavior alternatives
		outboundAttractionDistribution=False,  # Dissuade hubs
		linLogMode=False,  # NOT IMPLEMENTED
		adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
		edgeWeightInfluence=edgeWeightInfluence,
		
		# Performance
		jitterTolerance=jitterTolerance,  # Tolerance
		barnesHutOptimize=True,
		barnesHutTheta=barnesHutTheta,
		multiThreaded=False,  # NOT IMPLEMENTED
		
		# Tuning
		scalingRatio=scalingRatio,
		strongGravityMode=False,
		gravity=gravity,
		# Log
		verbose=verbose)
	
	positions = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=n_iter)
	positions = np.array([positions[i] for i in sorted(positions.keys())])
	return positions


# functions for preprocessing and pseudo-cell generations
#========================================================================#
def calAndIntDist(inputs, weights, method, args = { 'numpc' : 50 }):
	dist = pd.DataFrame(np.zeros((inputs[0].shape[1], inputs[0].shape[1])))
	for i in range(len(inputs)):
		dist = dist + weights[i] * calDist(inputs[i], method, args)
	return dist

#========================================================================#
def calDist(input, method, args = { 'numpc' : 50 }):
	dist = None
	if method == "cor":
		dist = 1 - input.corr()
	elif method == "pca":
		from sklearn.decomposition import PCA
		from sklearn.metrics.pairwise import euclidean_distances
		
		pca = PCA(n_components = args['numpc'])
		input = pca.fit_transform(input.transpose())
		dist = pd.DataFrame(euclidean_distances(input))
	elif method == "cos":
		from sklearn.metrics.pairwise import cosine_distances
		dist = pd.DataFrame(cosine_distances(input.transpose()))
	elif method == "cos0":
		from sklearn.metrics.pairwise import cosine_distances
		input = input.apply(lambda x: x - np.mean(x), 0).transpose()
		dist = pd.DataFrame(cosine_distances(input))
	elif method == "eucl":
		from sklearn.metrics.pairwise import euclidean_distances
		dist = pd.DataFrame(euclidean_distances(input.transpose()))
	
	return dist

#========================================================================#
def calDistMultipleInputs(input, weights, method, args = { 'numpc' : 50 }):
	dist = pd.DataFrame(np.zeros((input[0].shape[1], input[0].shape[1])))
	for i in range(len(input)):
		dist = dist + weights[i] * calDist(input[i], method, args)
	return dist

#========================================================================#
def coarseGrain(input, weights, dist, k, coarse_grain_X, cell_groupings_dict, method = "nntplus", dist_args = { 'dist_method' : 'cos0', 'args' : { 'numpc' : 50 } }):
	idx = coarseGrainIdx(input, weights, dist, k, coarse_grain_X, method, dist_args)
	inputPooled, groupingPooled = groupInputAndGroup(idx, input, cell_groupings_dict)
	
	return (idx, inputPooled, groupingPooled)

#========================================================================#
def coarseGrainIdx(input, weights, dist, k, coarse_grain_X, method = "nntplus", dist_args = { 'dist_method' : 'cos0', 'args' : { 'numpc' : 50 } }):
	n = dist.shape[0] if dist is not None else input.shape[1]
	if k > n:
		k = n
	
	if method == "downsampled":
		knn = pd.DataFrame(index = range(n), columns = range(k))
	else:
		if dist is not None:
			knn = getKNNFromDist(dist, k)
		else:
			knn = getKNNFromInput(input, weights, k, dist_args['dist_method'], dist_args['args'])
	
	node_index = generatePseudoCells(knn, coarse_grain_X, method = method)
	
	idx = pd.DataFrame([ np.array(range(knn.shape[0])), np.array(node_index) ]).transpose()
	idx.columns = [ 'initial', 'pooled' ]
	
	return idx

#========================================================================#
def groupInputAndGroup(idx, input, cell_groupings_dict):
	inputPooled = groupInput(idx, input)
	groupingPooled = groupGroupings(idx, cell_groupings_dict)
	
	return (inputPooled, groupingPooled)

#========================================================================#
def groupInput(idx, input):
	pooledIdx = np.sort(np.unique(idx.loc[idx.loc[:,'pooled'] != -1, 'pooled']))
	inputPooled = pd.DataFrame([ input.iloc[:, np.where(idx.loc[:,'pooled'] == x)[0]].mean(1) for x in pooledIdx ]).transpose()
	return inputPooled

#========================================================================#
def groupGroupings(idx, cell_groupings_dict):
	pooledIdx = np.sort(np.unique(idx.loc[idx.loc[:,'pooled'] != -1, 'pooled']))
	groupingPooled = {k:[] for k in cell_groupings_dict}
	for i in pooledIdx:
		for k in groupingPooled:
			groupingPooled[k].append(mostCommon([cell_groupings_dict[k][ii] for ii in np.nonzero(idx.iloc[:,1]==i)[0]]))
	
	return groupingPooled

#========================================================================#
def getKNNFromDist(dist, k):
	edges = list()
	for i in range(dist.shape[0]):
		row = np.ma.array(dist.iloc[i,:], mask = False)
		row.mask[i] = True
		idx = row.argsort()[:k]
		if np.isnan(row).sum() > len(row) - 1 - k:
			idx[(len(row) - 1 - np.isnan(row).sum()):] = -1
		edges.append(list(idx))
	edgesDf = pd.DataFrame(edges)
	return edgesDf

#========================================================================#
def getKNNFromInput(input, weights, k, method, args = { 'numpc' : 50}):
	knn = None
	if method == "pca":
		dist = pd.DataFrame(np.zeros((input[0].shape[1], input[0].shape[1])))
		for i in range(len(input)):
			dist = dist + weights[i] * calDist(input[i], "pca", args)
		knn = getKNNFromDist(dist, k)
	elif method == "cos0" or method == "cor" or method == "cos":
		if method == "cos0" or method == "cor":
			input = [ inp.apply(lambda x: x - np.mean(x), 0).transpose() for inp in input ]
		else:
			input = [ inp.transpose() for inp in input ]
		
		knn = pd.DataFrame([ getCosKnnForOne(input, weights, i, k) for i in range(len(input[0])) ])
	elif method == "eucl":
		input = [ inp.transpose() for inp in input ]
		knn = pd.DataFrame([ getEuclKnnForOne(input, weights, i, k) for i in range(len(input[0])) ])
	
	return knn

#========================================================================#
def getCosKnnForOne(input, weights, idx, k):
	from sklearn.metrics.pairwise import cosine_distances
	dist = np.repeat(0, input[0].shape[0])
	for i in range(len(input)):
		dist = dist + weights[i] * cosine_distances(input[i].iloc[idx,:].values.reshape(1, -1), input[i])[0]
	idxNN = np.argsort(dist)[:(k+1)]
	if np.sum(idxNN != idx) == k:
		idxNN = idxNN[idxNN != idx]
	
	return idxNN

#========================================================================#
def getEuclKnnForOne(input, weights, idx, k):
	from sklearn.metrics.pairwise import euclidean_distances
	dist = np.repeat(0, input[0].shape[0])
	for i in range(len(input)):
		dist = dist + weights[i] * euclidean_distances(input[i].iloc[idx,:].values.reshape(1, -1), input[i])[0]
	idxNN = np.argsort(dist)[:(k+1)]
	if np.sum(idxNN != idx) == k:
		idxNN = idxNN[idxNN != idx]
	
	return idxNN

#========================================================================#
def getEdges(knn):
	edges = set()
	for i in range(knn.shape[0]):
		for j in range(knn.shape[1]):
			if knn.iloc[i,j] != -1:
				edges.add(tuple(sorted([ i, knn.iloc[i,j] ])))
	return edges

#========================================================================#
def getEdgesList(knn):
	edges = list()
	for i in range(knn.shape[0]):
		for j in range(knn.shape[1]):
			if knn.iloc[i,j] != -1:
				edges.append((i, knn.iloc[i,j]))
	return edges

#========================================================================#
def getMutualEdges(knn):
	edgesSet = dict()
	for i in range(knn.shape[0]):
		for j in range(knn.shape[1]):
			if knn.iloc[i,j] != -1:
				edgesSet[tuple(sorted([ i, knn.iloc[i,j] ]))] = edgesSet.get(tuple(sorted([ i, knn.iloc[i,j] ])), 0) + 1
	edges = set()
	for edge in edgesSet.keys():
		if edgesSet[edge] == 2:
			edges.add(edge)
	return edges

#========================================================================#
def mostCommon(labels):
	counts = { x : 0 for x in np.unique([ str(x) for x in labels ]) }
	for x in labels: counts[str(x)] += 1
	counts['nan'] = -1
	maxLabel = pd.DataFrame([counts]).iloc[0,:].idxmax()
	return maxLabel

#========================================================================#
def generatePseudoCells(knn, coarse_grain_X, method = "simple", args = { 'initdist' : 1, 'maxdist' : 2, 'mergedist' : 2, 'maxmergingsize' : 5 }):
	N = knn.shape[0]
	node_index = np.ones(N) * -1
	base_filter = np.random.uniform(0,1,len(node_index)) < coarse_grain_X
	if base_filter.sum() == 0:
		base_filter = np.random.uniform(0,1,len(node_index))
		base_filter = base_filter == base_filter.min()
	
	if method == "downsample":
		node_index[base_filter] = np.arange(np.sum(base_filter))
		return node_index
	
	if method == "simple":
		node_index[base_filter] = np.arange(np.sum(base_filter))
		edges = getEdges(knn)
		for count,(i,j) in enumerate(edges):
			if node_index[i] != -1 and node_index[j] == -1:
				node_index[j] = node_index[i]
			if node_index[j] != -1 and node_index[i] == -1:
				node_index[i] = node_index[j]
	
	elif method == "knn":
		selected = pd.DataFrame([ np.isin(range(knn.shape[0]), knn.iloc[i,:]) for i in np.where(base_filter)[0] ])
		node_index = list(selected.apply(lambda x: np.where(x)[0][0], axis=0))
		
	elif method == "nnterritory" or method == "nntplus":
		import igraph
		edges = getEdges(knn)
		graph = igraph.Graph(directed = False)
		graph.add_vertices(knn.shape[0])
		graph.add_edges(edges)
		
		shortestDistance = pd.DataFrame(graph.shortest_paths(source = np.where(base_filter)[0]))
		selected = shortestDistance <= args['initdist']
		for i in range(np.sum(base_filter)): selected.iloc[:, np.where(base_filter)[0][i]] = np.repeat([False, True, False], [i, 1, np.sum(base_filter)-i-1])
		
		belongToMulti = selected.apply(lambda x: np.sum(x)>1, axis = 0)
		selected.loc[:, belongToMulti] = selected.loc[:, belongToMulti].apply(lambda x: np.array(range(len(x))) == np.where(x)[0][np.random.rand(np.sum(x)).argmax()], axis=0)
		
		if args['maxdist'] > args['initdist']:
			belongToNowhereButNotFar = np.where(selected.mean(0) == 0 & shortestDistance.apply(lambda x: np.min(x[x!=0]) <= args['maxdist'], axis=0))[0]
			candidates = [ np.where(shortestDistance.iloc[:,i] == np.min(shortestDistance.iloc[:,i]))[0] for i in belongToNowhereButNotFar ]
			candidates = [ x[np.random.rand(len(x)).argmax()] for x in candidates ]
			newSel = pd.DataFrame([ np.repeat([False, True, False], [i, 1, np.sum(base_filter)-i-1]) for i in candidates ]).transpose()
			newSel.columns = belongToNowhereButNotFar
			selected.update(newSel)
		
		if method == "nntplus":
			shortestDistanceMerged = shortestDistance.copy()
			distPassed = shortestDistanceMerged.iloc[:,np.where(base_filter)[0]] <= args['mergedist']
			sizePassed = selected.sum(1) <= args['maxmergingsize']
			sizePassed = np.dot(pd.DataFrame(sizePassed), pd.DataFrame(sizePassed).transpose()) == 1
			mergePassed = np.array(distPassed & sizePassed)
			np.fill_diagonal(mergePassed, False)
			
			numIter = 0
			while np.sum(mergePassed.sum()) > 0:
				pairs = np.where(mergePassed)
				pairs = [ (pairs[0][i], pairs[1][i]) for i in range(len(pairs[0])) ]
				selPair = pairs[np.random.rand(len(pairs)).argmax()]
				
				selected.iloc[selPair[0],:] = selected.iloc[selPair[0],:] | selected.iloc[selPair[1],:]
				selected = selected.drop(index = selPair[1])
				selected.index = range(selected.shape[0])
				shortestDistanceMerged = shortestDistanceMerged.drop(index = selPair[1])
				shortestDistanceMerged.index = range(shortestDistanceMerged.shape[0])
				base_filter[np.where(base_filter)[0][selPair[1]]] = False
				
				distPassed = shortestDistanceMerged.iloc[:,np.where(base_filter)[0]] <= args['mergedist']
				sizePassed = selected.sum(1) <= args['maxmergingsize']
				sizePassed = np.dot(pd.DataFrame(sizePassed), pd.DataFrame(sizePassed).transpose()) == 1
				mergePassed = np.array(distPassed & sizePassed)
				np.fill_diagonal(mergePassed, False)
				
				numIter += 1
				#print "done iteration " + str(numIter) + ". # of mergable pairs: " + str(np.sum(mergePassed.sum()))
		
		node_index = list(selected.apply(lambda x: np.where(x)[0][0], axis=0))
	
	elif method == "som":
		pass
	
	return node_index



# run
#========================================================================#
if __name__ == "__main__":
	main()

