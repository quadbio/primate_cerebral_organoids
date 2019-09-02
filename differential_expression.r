# function to get cells indices by random sampling with line number restricted and cell distribution along pseudotime controlled
# used in human-chimp organoid DE robustness estimation
sampling_ctrl_lines_pt <- function(lines, pt, pt_ref, num_lines, line_candidates = unique(lines), num_breaks_pt = 10){
	num_cells_int <- table(ceiling(pt_ref / max(c(pt,pt_ref)) * num_breaks_pt))
	lines_to_use <- sample(lines, min(c(num_lines, length(unique(lines)))))
	
	sel_idx <- sort(do.call(c, lapply(1:length(num_cells_int), function(idx_int){
		cell_pool_idx <- intersect(which(ceiling(pt/max(c(pt,pt_ref)) * num_breaks_pt) == idx_int), which(lines %in% lines_to_use))
		cells_sel <- sample(cell_pool_idx, num_cells_int[idx_int], replace=T)
		return(cells_sel)
	})))
	sel_pt <- pt[sel_idx]
	res <- list(sel_lines = lines_to_use, sel_idx = sel_idx, sel_pt = sel_pt)
	return(res)
}

# function to do DE along the aligned pseudotime trajectory
# used in human-chimp organoid DE analysis
DE_ftest_pt <- function(expr1, pt1, expr2, pt2, degree = 6, num_breaks_gdiff = 10, num_threads = 1){
	require(doParallel)
	require(splines)
	registerDoParallel(num_threads)

	expr1 <- expr1[intersect(rownames(expr1), rownames(expr2)),]
	expr2 <- expr2[rownames(expr1),]
	pt <- c(pt1, pt2)
	g <- as.factor(rep(c("g1","g2"), c(ncol(expr1), ncol(expr2))))
	res <- foreach(i = 1:nrow(expr1), .combine = rbind) %dopar%{
		e1 <- as.numeric(expr1[i,])
		e2 <- as.numeric(expr2[i,])
		avg1 <- sapply(1:num_breaks_gdiff, function(idx) mean(e1[which(ceiling(pt1 * num_breaks_gdiff) == idx)], na.rm=T))
		avg2 <- sapply(1:num_breaks_gdiff, function(idx) mean(e2[which(ceiling(pt2 * num_breaks_gdiff) == idx)], na.rm=T))
		gdiff <- mean(avg2 - avg1, na.rm=T)
		
		e <- c(e1, e2)
		m0 <- lm(e ~ ns(pt, df = degree))
		m1 <- lm(e ~ g * ns(pt, degree = degree))
		f <- anova(m1)$"Mean Sq"[4] / anova(m0)$"Mean Sq"[2]
		p <- pf(f, m1$df, m0$df, lower.tail=T)

		res <- c(g_diff = gdiff, f = f, p = p)
		return(res)
	}
	stopImplicitCluster()
	
	res <- setNames(data.frame(res, row.names = rownames(expr1)), c("generalizedDiff", "F", "P"))
	return(res)
}


