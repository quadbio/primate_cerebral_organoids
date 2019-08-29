# function: Dynamic-time-warping based alignment
align_diff <- function(diff, ref_name="ref", query_name="query", mode=c("global","fixed_start","fixed_end","partial"), rev=FALSE){
  mode <- match.arg(mode)
  diff[is.na(diff)] <- 0
  if(rev) diff <- t(diff)
  
  acc_dist <- matrix(NA, nrow=nrow(diff), ncol=ncol(diff))
  if(mode == "fixed_end" | mode == "partial"){
    acc_dist[1,] <- diff[1,]
  } else{
    acc_dist[1,1] <- diff[1,1]
    for(i in 2:ncol(diff)) acc_dist[1,i] <- diff[1,1] + sum(diff[1,2:(i-1)])*0.5
  }
  for(i in 2:nrow(acc_dist)) acc_dist[i,1] <- diff[1,1] + sum(diff[2:(i-1),1])*0.5
  
  direction <- matrix(NA, nrow=nrow(diff), ncol=ncol(diff))
  direction[1,-1] <- "left"
  direction[-1,1] <- "up"
  direction[1,1] <- "start"
  
  for(i in 2:nrow(acc_dist)){
    for(j in 2:ncol(acc_dist)){
      source <- c(acc_dist[i-1,j-1], acc_dist[i,j-1], acc_dist[i-1,j])
      idx <- which.min(source)
      direction[i,j] <- c("diag","left","up")[idx]
      acc_dist[i,j] <- min(source)+diff[i,j]*ifelse(idx == 1, 1, 0.5)
    }}
  
  links <- data.frame(c(),c())
  i <- nrow(acc_dist); j <- ifelse(mode == "partial" | mode == "fixed_start", which.min(acc_dist[nrow(acc_dist),]), ncol(acc_dist))
  while((mode %in% c("partial","fixed_end") & i>1) | (mode %in% c("global","fixed_start") & (i>1 | j>1))){
    links <- rbind(links, c(i,j))
    if(direction[i,j]=="up" | direction[i,j]=="diag") i <- i-1
    if(direction[i,j]=="left" | direction[i,j]=="diag") j <- j-1
  }
  links <- rbind(links, c(i,j))
  
  if(rev){
    links <- links[,2:1]
    diff <- t(diff)
    aligned_query_pseudotime <- sapply(sort(unique(links[,1])), function(i) mean(links[links[,1]==i, 2]) )
    if(length(aligned_query_pseudotime) < nrow(diff))
      aligned_query_pseudotime <- c(aligned_query_pseudotime, seq(max(links[,2])+1, max(links[,2])+nrow(diff)-length(aligned_query_pseudotime), 1))
  } else{
    aligned_query_pseudotime <- sapply(1:nrow(acc_dist), function(i) mean(links[links[,1]==i, 2]) ) # alignment trajectory
  }
  colnames(links) <- c(query_name,ref_name)
  
  res <- list(diff_mat = diff, alignment = links, aligned_query = aligned_query_pseudotime)
  return(res)
}

# function: to determine start/end of reference trajectory for truncated pseudotime alignment
get_truncate_points <- function(input_ref, pt_ref, input_query, pt_query, trunc_at = c("end","start","both"), method = c("svr","cor"), ..., nknots_cobs = 20, num_breaks_ref = 50, num_breaks_query = num_breaks_ref, tol = 0.05){
	trunc_at <- match.arg(trunc_at)
	method <- match.arg(method)

	if (method == "svr"){
		require(e1071)
		require(cobs)
		
		model_proj <-  svm(x = t(input_ref), y = pt_ref, ...)
		proj_pt_query <- predict(model_proj, t(input_query))
		pointwise_constraint <- NULL
		if (trunc_at == "end"){
			pointwise_constraint <- rbind(c(0,0,0))
		} else if (trunc_at == "start"){
			pointwise_constraint <- rbind(c(0,1,1))
		}
		model_cobs <- cobs(pt_query,proj_pt_query, constraint = "increase", pointwise=pointwise_constraint, nknots=nknots_cobs)
		pred_pt_query <- predict(model_cobs, pt_query)[,2]
		pt_trunc <- c(start = ifelse(trunc_at == "end", min(pt_ref), min(pred_pt_query)), end = ifelse(trunc_at == "start", max(pt_ref), max(pred_pt_query)))
	} else if (method == "cor"){
		average_ref <- t(apply(input_ref, 1, function(e) tapply(e, ceiling(pt_ref * num_breaks_ref), mean)))
		average_query <- t(apply(input_query, 1, function(e) tapply(e, ceiling(pt_query * num_breaks_query), mean)))
		cor_ref_query <- cor(average_ref, average_query)
		trunc_start <- (min(which(corr_ref_query[,1]>=max(corr_ref_query[,1])-(max(corr_ref_query[,1])-min(corr_ref_query[,1]))*tol))-1)/num_breaks_ref
		trunc_end <- max(which(corr_ref_query[,num_breaks_query]>=max(corr_ref_query[,num_breaks_query])-(max(corr_ref_query[,num_breaks_query])-min(corr_ref_query[,num_breaks_query]))*tol))/num_breaks_ref
		pt_trunc <- c(start = ifelse(trunc_at == "end", min(pt_ref), trunc_start), end = ifelse(trunc_at == "start", max(pt_ref), trunc_end))
	}

	return(pt_trunc)
}

# wrapper: alignment given expression matrix and pseudotimes of cells in reference trajectory and query trajectories
align_pt_traj <- function(expr_ref, pt_ref, expr_query, pt_query, # data input
		num_breaks_ref = 50, num_breaks_query = num_breaks_ref, dist_method = c("cor", "eucl"), degree = 1, # dtw alignment parameters
		ref_name = "ref", query_name = "query", mode = c("global", "fixed_start", "fixed_end", "partial"), rev = FALSE, # dtw alignment parameters
		nknots_cobs = 20){ # cobs parameters
	require(pdist)
	require(cobs)
	dist_method <- match.arg(dist_method)
	mode <- match.arg(mode)
	
	average_ref <- t(apply(expr_ref, 1, function(e) tapply(e, ceiling(pt_ref * num_breaks_ref), mean)))
	average_query <- t(apply(expr_query, 1, function(e) tapply(e, ceiling(pt_query * num_breaks_query), mean)))

	average_diff <- NULL
	if (dist_method == "cor"){
		average_diff <- 1-cor(average_query, average_ref)
	} else if (dist_method == "eucl"){
		average_diff <- as.matrix(pdist::pdist(t(average_query), t(average_ref)))
	}

	if (degree > 1)
		average_diff <- average_diff ^ degree
	
	alignment <- align_diff(average_diff, ref_name = ref_name, query_name = query_name, mode = mode, rev = rev)
	pointwise_constrain <- NULL
	if (mode == "global"){
		pointwise_constrain <- rbind(c(0,0,0), c(0,1,1))
	} else if (mode == "fixed_start"){
		pointwise_constrain <- rbind(c(0,0,0))
	} else if (mode == "fixed_end"){
		pointwise_constrain <- rbind(c(0,1,1))
	}
	model <- cobs(x = c(0,alignment$alignment[,1]) / num_breaks_query,
				  y = c(0,alignment$alignment[,2]) / num_breaks_ref,
				  constraint= "increase", pointwise=pointwise_constrain, nknots=nknots_cobs)
	
	pt_query_aligned <- predict(model, pt_query)[,2]
	res <- list(alignment = alignment, model = model, pt_query_aligned = pt_query_aligned)
	return(res)
}

# wrapper: truncated alignment given expression matrix, alternative representation matrix and pseudotimes of cells in reference and query trajectories
trunc_align_pt_traj <- function(expr_ref, alt_ref = expr_ref, pt_ref, expr_query, alt_query = expr_query, pt_query, # data input
		which_trunc = c("query", "ref"), trunc_at = "end", trunc_method = "svr", ..., nknots_cobs = 20, num_breaks_ref = 50, num_breaks_query = num_breaks_ref, tol = 0.05, # truncate pt point estimation parameters
		dist_method = "cor", degree = 1, ref_name = "ref", query_name = "query"){ # dtw alignment parameters
	which_trunc <- match.arg(which_trunc)
	
	trunc_pt <- c(start = NA, end = NA)
	expr_ref_trunc <- expr_ref
	pt_ref_trunc <- pt_ref
	expr_query_trunc <- expr_query
	pt_query_trunc <- pt_query
	if (which_trunc == "ref"){
		trunc_pt <- get_truncate_points(alt_ref, pt_ref, alt_query, pt_query, trunc_at = trunc_at, method = trunc_method, ..., nknots_cobs = nknots_cobs, num_breaks_ref = num_breaks_ref, num_breaks_query = num_breaks_query, tol = tol)
		idx_ref <- which(pt_ref >= trunc_pt["start"] & pt_ref <= trunc_pt["end"])
		expr_ref_trunc <- expr_ref[,idx_ref]
		pt_ref_trunc <- rank(pt_ref[idx_ref])/length(idx_ref)
	} else if (which_trunc == "query"){
		trunc_pt <- get_truncate_points(alt_query, pt_query, alt_ref, pt_ref, trunc_at = trunc_at, method = trunc_method, ..., nknots_cobs = nknots_cobs, num_breaks_ref = num_breaks_query, num_breaks_query = num_breaks_ref, tol = tol)
		idx_query <- which(pt_query >= trunc_pt["start"] & pt_query <= trunc_pt["end"])
		expr_query_trunc <- expr_query[,idx_query]
		pt_query_trunc <- rank(pt_query[idx_query])/length(idx_query)
	}

	average_ref <- t(apply(expr_ref_trunc, 1, function(e) tapply(e, ceiling(pt_ref_trunc * num_breaks_ref), mean)))
	average_query <- t(apply(expr_query_trunc, 1, function(e) tapply(e, ceiling(pt_query_trunc * num_breaks_query), mean)))
	average_diff <- NULL
	if (dist_method == "cor"){
		average_diff <- 1-cor(average_query, average_ref)
	} else if (dist_method == "eucl"){
		average_diff <- as.matrix(pdist::pdist(t(average_query), t(average_ref)))
	}
	alignment_trunc <- align_diff(average_diff, ref_name = ref_name, query_name = query_name, mode = "global", rev = F)
	model <- NULL
	if (which_trunc == "ref"){
		model <- cobs(x = c(0,alignment_trunc$alignment[,1]) / num_breaks_query,
					  y = c(0,alignment_trunc$alignment[,2]) / num_breaks_ref * (trunc_pt["end"] - trunc_pt["start"]) + trunc_pt["start"],
					  constraint = "increase", pointwise = rbind(c(0,0,trunc_pt["start"]), c(0,1,trunc_pt["end"])))
	} else if (which_trunc == "query"){
		model <- cobs(x = c(0,alignment_trunc$alignment[,1]) / num_breaks_query * (trunc_pt["end"] - trunc_pt["start"]) + trunc_pt["start"],
					  y = c(0,alignment_trunc$alignment[,2]) / num_breaks_ref,
					  constraint = "increase", pointwise = rbind(c(0,trunc_pt["start"],0), c(0,trunc_pt["end"],1)))
	}

	pt_query_aligned <- predict(model, pt_query)[,2]
	res <- list(trunc_traj = which_trunc, trunc_pt = trunc_pt, model = model, pt_query_aligned = pt_query_aligned)
	return(res)
}


