# function to subsample cells/nuclei by control subtype heterogeneity and cell number differences among species
# used in adult expression comparison between human and chimp
sample_cells_ctrl_hetero <- function(species, cl, cl_annot, cts_to_sample = unique(cl_annot), num_cells_per_ct = 200){
	res <- setNames(lapply(cts_to_sample, function(ct){
		cl_ct <- names(cl_annot)[which(cl_annot == ct)]
		num_subtypes <- setNames(table(ceiling(runif(num_cells_per_ct) * length(cl_ct))), cl_ct)
		idx_sel_cl <- lapply(cl_ct, function(cl_sel){ sapply(unique(species), function(species_sel){
			idx_candidates <- which(cl == cl_sel & species == species_sel)
			idx_sel <- sample(idx_candidates, num_subtypes[cl_sel], replace = TRUE)
			return(idx_sel)
		}) })
		idx_sel_ct <- setNames(data.frame(do.call(rbind, idx_sel_cl)), unique(species))
		return(idx_sel_ct)
	}), cts_to_sample)
	return(res)
}

# function to test for difference of gene detection rates between two groups
# used in organoid chromatin accessibility comparison and adult expression comparison between human and chimp
differential_detection_rate_test <- function(expr1, expr2, genes2test = intersect(rownames(expr1), rownames(expr2)), num_threads = 1){
	require(doParallel)
	expr1 <- expr1[intersect(rownames(expr1), rownames(expr2)),]
	expr2 <- expr2[rownames(expr1),]
	log2fc <- log2(apply(exp(expr1) - 1, 1, mean) / apply(exp(expr2) - 1, 1, mean))

	expr1 <- Matrix(expr1)
	bin1 <- summary(expr1)
	bin1 <- sparseMatrix(i = bin1[,1], j = bin1[,2], x = 1, dims = dim(expr1), dimnames = dimnames(expr1))
	expr2 <- Matrix(expr2)
	bin2 <- summary(expr2)
	bin2 <- sparseMatrix(i = bin2[,1], j = bin2[,2], x = 1, dims = dim(expr2), dimnames = dimnames(expr2))
	
	bg <- mean(apply(bin1, 2, sum) / nrow(bin1)) / mean(apply(bin2, 2, sum) / nrow(bin2))
	g <- as.factor(rep(c("g1","g2"), c(ncol(bin1), ncol(bin2))))
	registerDoParallel(num_threads)
	res <- setNames(data.frame(log2fc = log2fc[intersect(genes2test, rownames(expr1))], foreach(gene = intersect(genes2test, rownames(expr1)), .combine = rbind) %dopar%{
		b1 <- as.logical(bin1[gene,])
		b2 <- as.logical(bin2[gene,])
		pt1 <- sum(b1) / length(b1)
		pt2 <- sum(b2) / length(b2)

		y <- c(b1, b2)
		m0 <- try(glm(y ~ I(ifelse(g=="g1", bg, 1)) + 0, family=binomial('identity')), silent=T)
		m1 <- glm(y ~ g, family = "binomial")
		p <- ifelse(sum("try-error" %in% class(m0))>0, NA, anova(m0, m1, test="Chisq")$Pr[2])
		res <- c(pt1 = pt1, pt2 = pt2, p = p)
		return(p)
	}, row.names = intersect(genes2test, rownames(expr1))), c("log2FC", "Pt1", "Pt2", "P"))
	stopImplicitCluster()

	return(res)
}

