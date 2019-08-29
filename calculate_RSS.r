#' Load the Allen Brain Atlas reference set
#'
#' This function loads the Allen Brain Atlas reference data matrix. Different
#' pre-defined gene lists or sample lists allow prefiltering of the loaded data.
#'
#' @param sampleSet 	The pre-defined sample list for output.
#' @param geneSet 	The pre-defined gene list for output.
#' @param geneList 	The customized gene list for output. When it is NULL, the pre-defined gene list will be used.
#' @param excludeG2M 	When it is TRUE, G2 and M phrase marker genes will be excluded from the gene list.
#' @param logTransform 	When it is TRUE, log-transformation with pseudo-count will be applied to the expression matrix (in RPKM).
#' @return A list of the loaded data, in which there are three components: 'expr' is the data matrix, 'rowMeta' is the row (gene) information, and 'colMeta' is the column (reference sample) information.
#' @export
retrieveABARef <- function(sampleSet = c("full","fetal","fetal_before10pcw","fetal_after10pcw"),
	geneSet = c("full", "highvar", "fetal_highvar", "region_fetalB10", "region_fetalA10", "region_fetalB10-2", "regionGroup_fetalB10", "regionGroup_fetalA10"),
	geneList = NULL,
	excludeG2M = TRUE, logTransform = TRUE){
	
	sampleSet <- match.arg(sampleSet)
	geneSet <- match.arg(geneSet)
	load("data/rows_ABA.rda")
	load("data/columns_ABA.rda")
	load("data/expr_ABA.rda")
	
	if (sampleSet == "full"){
		idxCol <- 1:ncol(expr_ABA)
		idxRow <- which(rows_ABA$expressed_full & rows_ABA$gene_type == "protein_coding")
	} else if (sampleSet == "fetal"){
		idxCol <- which(columns_ABA$fetal_before10pcw | columns_ABA$fetal_after10pcw)
		idxRow <- which(rows_ABA$expressed_fetal & rows_ABA$gene_type == "protein_coding")
	} else{
		idxCol <- which(columns_ABA[,sampleSet])
		idxRow <- which(rows_ABA$expressed_fetal & rows_ABA$gene_type == "protein_coding")
	}
	colDat <- columns_ABA[idxCol,]
	rownames(colDat) <- sprintf("Ref%05d", 1:nrow(colDat))
	
	if (is.null(geneList)) {
		if (geneSet == "highvar"){
			idxRow <- which(p.adjust(rows_ABA$p_highvar, method="BH")<0.1)
		} else if (geneSet == "fetal_highvar"){
			idxRow <- which(p.adjust(rows_ABA$p_highvar_fetal, method="BH")<0.1)
		} else if (geneSet == "region_fetalB10"){
			idxRow <- which(p.adjust(rows_ABA$p_region_b10pcw, method="BH")<0.01 & rows_ABA$log2fc_region_b10pcw>1)
		} else if (geneSet == "region_fetalA10"){
			idxRow <- which(p.adjust(rows_ABA$p_region_a10pcw)<0.01 & rows_ABA$log2fc_region_a10pcw>1)
		} else if (geneSet == "region_fetalB10-2"){
			idxRow <- which(p.adjust(rows_ABA$p_region_b10pcw)<0.05 & rows_ABA$log2fc_region_b10pcw>1)
		} else if (geneSet == "regionGroup_fetalB10"){
			idxRow <- which(p.adjust(rows_ABA$p_regiongroup_b10pcw)<0.01 & rows_ABA$log2fc_regiongroup_b10pcw>1)
		} else if (geneSet == "regionGroup_fetalA10"){
			idxRow <- which(p.adjust(rows_ABA$p_regiongroup_a10pcw)<0.01 & rows_ABA$log2fc_regiongroup_a10pcw>1)
		}
	} else{
		numMatch <- apply(rows_ABA[,3:5], 2, function(x){ sum(x %in% geneList) })
		if (max(numMatch) < length(geneList) / 2)
			stop("Genes matching the ABA reference data set are too few. Make sure that they are in ENSEMBL IDs, gene symbols or Entrez IDs.")
		idxRow <- which(rows_ABA[,which.max(numMatch) + 2] %in% geneList)
	}
	if (excludeG2M) idxRow <- setdiff(idxRow, which(rows_ABA$G2M))
	rowDat <- rows_ABA[idxRow,]
	
	dat <- expr_ABA[idxRow, idxCol]
	colnames(dat) <- rownames(colDat)
	if (logTransform)
		dat <- log(dat + 1)
	
	res <- list(expr = dat, rowMeta = rowDat, colMeta = colDat)
	return(res)
}

#' Calculate RSS matrix
#'
#' This function calculates the Reference Similarity Spectrum (RSS) matrix,
#' based on the given input matrix and the reference data set matrix.
#'
#' @param input 	The input expression matrix, assuming its rows as genes and columns as samples. Its rows are expected to be named by the same types of gene representative (e.g. ENSEMBL IDs) as ref.
#' @param ref 	The reference expression matrix, assuming its rows as genes and columns as samples.
#' @param method 	The type of correlation to calculate. It allows any method used in the cor function.
#' @param scale 	When it is TRUE, similarities of one sample (cell) in the input matrix across the reference samples are central-scaled transformed.
#' @param threads 	The number of threads to calculate correlations.
#' @return The RSS matrix, with rows for the input samples and columns for the reference samples.
#' @export
calculateRSS <- function(input, ref, method = "pearson", scale = TRUE, threads = 1){
	candidates <- intersect(rownames(input), rownames(ref))

	corr <- NA
	corr <- cor(input[candidates,], ref[candidates,], method = method)
	
	if (scale){
		corr <- t(scale(t(corr)))
	}
	return(corr)
}


#' Plot feature's gradients/classes across all samples (cells) given the plotting coordinates
#'
#' This function plot the samples based on the given coordinates, coloring
#' each dot based on the value of its corresponding feature value/class.
#'
#' @param coord 	The plotting coordinates, expected to be of two columns. Each row represents one dot for one sample (cell).
#' @param value 	The values to be shown, with each value representing one sample (cell) with the same order as coord.
#' @param emphasize 	The Indices of samples (cells) to be emphasized. When it is set, the colors of all the other cells are set to be #bdbdbd30.
#' @param col 	The customized colors of dots. Either col or value must be given.
#' @param ... 	Other arguments passing to the plot function.
#' @export
plotFeature <- function(coord, values = NULL, emphasize = NULL, col = NULL, colorPal = NULL, edges = NULL, main = NA, xlab = "Dim-1", ylab = "Dim-2", cex.main = 1, cex.lab = 1, cex.axis = 1, lwd.edges = 0.2, col.edges = "#bdbdbd50", cex = 1, ...){
	if (is.null(col)){
		if (is.numeric(values)){
			if(is.null(colorPal)){
				colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
			}
			cellColor <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(values, breaks=30, right=F,include.lowest=T))]
			if(min(values, na.rm=T) == 0) cellColor[values == 0] <- "#bdbdbd30"
		} else{
			cols <- setNames(scales::hue_pal()(length(unique(values))), unique(values))
			cellColor <- cols[values]
		}
	} else{
		if (length(col) == 1) col <- rep(col, nrow(coord))
		cellColor <- col
	}
	
	plot(coord, type="n", main = main, xlab = xlab, ylab = ylab, cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis, ...)
	if (! is.null(edges) & (is.matrix(edges) | is.data.frame(edges))){
		for(i in 1:nrow(edges)) lines(coord[as.numeric(edges[i,]),1], coord[as.numeric(edges[i,]),2], lwd = lwd.edges, col = col.edges)
	}
	if (is.null(emphasize)){
		points(coord, col = cellColor, pch=16, cex = cex)
	} else{
		points(coord, col = "#efefef30", pch=16, cex = cex)
		points(coord[emphasize, ], col = cellColor[emphasize], pch=16, cex = cex)
	}
}


