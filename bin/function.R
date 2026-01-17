
#get gabriel graph data
build_nb_gabriel <- function(loc) {
  xy <- as.matrix(loc[, c("x","y")])
  adegenet::chooseCN(xy = xy, result.type = "nb", plot.nb = FALSE, type = 2)
}

#Convert the Monmonier's  data into boudary data used for drawing
make_arrow_segments_all <- function(bound_df) {
  bound_df %>%
    arrange(run_id, dir, step) %>%
    group_by(run_id, dir) %>%
    mutate(xend = dplyr::lead(bound_x),
           yend = dplyr::lead(bound_y)) %>%
    ungroup() %>%
    filter(!is.na(xend) & !is.na(yend))
}

#Convert the adjacent nb into edge data
nb_to_edges_df <- function(nb, xy, samples = NULL) {
  elist <- lapply(seq_along(nb), function(i) {
    if (!length(nb[[i]])) return(NULL)
    data.frame(
      x    = rep(xy[i,1], length(nb[[i]])),
      y    = rep(xy[i,2], length(nb[[i]])),
      xend = xy[nb[[i]],1],
      yend = xy[nb[[i]],2],
      from = if (is.null(samples)) i else rep(samples[i], length(nb[[i]])),
      to   = if (is.null(samples)) nb[[i]] else samples[nb[[i]]],
      row.names = NULL, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, elist)
}

#Plot phylogenetic tree with tip clusters and sample labels
get_tree <- function(rpar,patient_SH,group){
  rpar1 <- drop.tip(rpar, "GL") 
  annotation <- data.frame(label=names(group),group=as.character(group),sample=gsub(paste0(patient,"_"),"",names(group)))
  trs <- full_join(rpar1,annotation,by='label')
  p1 <- ggtree(trs) + 
    geom_tippoint(aes(color=group),shape=16,size=20/.pt)+
    geom_treescale(x = -50, y = 0 , width = 50) + #控制tree的标尺位置
    geom_text2(aes(label=sample),size=8/.pt,color="white")+
    scale_color_manual(values = tree_color, name="class")+
    geom_rootedge(rootedge = rpar$edge.length[length(rpar$edge.length)]) +
    ggtitle(paste0(patient_SH))+
    geom_rootpoint()+
    theme(plot.title = element_text(size=10),
          legend.text = element_text(size=8),
          legend.title = element_text(size=10)
    ) +
    plot_layout(widths = unit(c(6), "cm"), heights = unit(c(10), "cm"))
  p1
  return(p1)
}

# Build a Gabriel neighborhood graph with optional tolerance.
# If tol <= 0 or the dataset is large, fall back to the fast spdep
# implementation. Otherwise, compute a tolerance-adjusted Gabriel graph
# by checking whether any third point falls inside the contracted circle
# between each pair of points.
build_gabriel_nb_tol <- function(coords, tol = 0, fallback_n = 400) {
  n <- nrow(coords)
  # Large datasets or tol <= 0: use spdep's Gabriel graph (no tolerance)
  if (n <= 1) return(vector("list", n))
  if (tol <= 0 || n > fallback_n) {
    gab <- try(spdep::gabrielneigh(coords), silent = TRUE)
    if (inherits(gab, "try-error")) {
      return(spdep::knn2nb(spdep::knearneigh(coords, k = pmin(2, n - 1))))
    } else {
      return(spdep::graph2nb(gab, sym = TRUE))
    }
  }
  
  # Small dataset + tol > 0: compute tolerance-adjusted Gabriel graph
  nb <- vector("list", n)
   # Check each pair (i, j)
  for (i in 1:(n - 1)) {
    xi <- coords[i, ]
    for (j in (i + 1):n) {
      xj <- coords[j, ]
      m  <- (xi + xj) / 2
      r  <- sqrt(sum((xi - xj)^2)) / 2
      r_ <- max(r - tol, 0)  # tolerance-shrunk radius
      
      inside <- FALSE
      if (r_ == 0) {
        for (k in 1:n) if (k != i && k != j) {
          if (all(abs(coords[k, ] - m) < 1e-12)) { inside <- TRUE; break }
        }
      } else {
        # General case: check if any point falls inside contracted circle
        for (k in 1:n) if (k != i && k != j) {
          if (sqrt(sum((coords[k, ] - m)^2)) <= r_ + 1e-12) { inside <- TRUE; break }
        }
      }
      
      if (!inside) {
        nb[[i]] <- c(nb[[i]], j)
        nb[[j]] <- c(nb[[j]], i)
      }
    }
  }
  nb
}


#Compute pairwise intratumor heterogeneity (ITH) between two tips
###get ITH
cal_ITH=function(dist_gene,two_tip){
  dist_gene=as.matrix(dist_gene)
  n_length1=dist_gene[two_tip[1],"GL"]
  n_length2=dist_gene[two_tip[2],"GL"]
  tlength=dist_gene[two_tip[1],two_tip[2]]
  ITH= 2*tlength/(n_length2+n_length1+tlength) #The proportion of non-shared mutations among all mutations
  return(ITH)
}
get_ITH_matrix <- function(binary_mat,select_info,r){
  #Pairwise physical distances between samples (as vector)
  pos_dist <- unlist(dist(select_info)) #sample distance
  diameter = quantile(pos_dist,probs = seq(0, 1, 0.25))[3] #sample distance median
  dist_gene=dist.hamming(binary_mat) ## Hamming distance between samples (genetic distance)
  #Define grid bounding box
  x1=min(select_info$X)-diameter;x2=max(select_info$X)+diameter
  y1=min(select_info$Y)-diameter;y2=max(select_info$Y)+diameter
  # Create regular grid with step size based on tumor radius r (largest diameter)
  net_x=seq(x1,x2,by=r/40) 
  net_y=seq(y1,y2,by=r/40) 
  # Create all grid points
  all_point=as.data.frame(x=c(),col.names = c("x","y"))
  for (x in net_x){
    all_point1=as.data.frame(x=c())
    for (y in net_y){
      all_point11=cbind.data.frame(x,y)
      all_point1=rbind(all_point1,all_point11)
    }
    all_point=rbind(all_point,all_point1)
  }
  #Compute local ITH for each grid location
  f=c()
  ITH=c()
  for (r in 1:nrow(all_point)){
    distance=apply(select_info[,1:2],1,function(x) sqrt((all_point[r,1]-x[1])**2+(all_point[r,2]-x[2])**2)) #Distance from grid point r to all samples
    count=length(which(distance<=diameter/2)) # Count samples within radius = diameter/2
    if (count>=2){
      f=c(f,r) 
      pair=combn(names(distance)[which(distance<=diameter/2)],2) #All pairs of samples within this local neighborhood
      ITH=c(ITH,mean(apply(pair,2,function(x) cal_ITH(as.matrix(dist_gene),x)))) #Average pairwise ITH within local region
    }
  }
  all_point_ITH <- cbind(all_point[f,],ITH)
  return(all_point_ITH)
}

#Compute Calinski–Harabasz (CH) index for a given clustering
compute_ch_index <- function(distance_matrix, group_labels) {
  sse <- 0
  N <- dim(distance_matrix)[1]
  x <- group_labels
  # Within-cluster sum of squares (SSE)
  for(i in names(table(x))){
    sample_location <- which(x==i)
    sse <- sse+sum(distance_matrix[sample_location,sample_location])/length(sample_location)
  }
  K <- length(unique(x)) # number of clusters
  ssb <- sum(distance_matrix)/N-sse # between-cluster sum of squares (SSB)
  ch <- (ssb*(N-K))/(sse*(K-1))
  return(ch)
}

#Check if a node is a tip (leaf) in a phylogenetic tree
is_leaf_node <- function(tree, node) {
  return(node <= length(tree$tip.label))
}

#Assign cluster labels to tree tips based on cut edges
assign_groups <- function(tree, cut_edge_ids, k) {
  # Initialize tip group labels as NA
  group_labels <- rep(NA, length(tree$tip.label))
  names(group_labels) <- tree$tip.label
  # Group counter
  group_counter <- 1
  # For each cut edge, assign all descendant tips to a new group
  for (edge_id in cut_edge_ids) {
    node <- tree$edge[edge_id, 2]
    # Extract clade rooted at this node
    clade <- extract.clade(tree, node)
    # Assign group label to all tips in this clade
    group_labels[clade$tip.label] <- group_counter
    group_counter <- group_counter + 1
  }
  # Assign remaining tips (not covered by any cut) to the last group
  remaining_tips <- names(group_labels)[is.na(group_labels)]
  if (length(remaining_tips) > 0) {
    group_labels[remaining_tips] <- group_counter
  }
  
  return(group_labels)
}

#Explore all combinations of internal-edge cuts for k=2..k_max and select clustering with maximum CH index for each k.
find_best_k_cutree_all_combinations <- function(tree, k_min = 2, k_max = 5, use_parallel = TRUE) {
  # Tree must be rooted
  if (!is.rooted(tree)) {
    stop("Tree must be rooted. Please root the tree using root().")
  }
  # Compute pairwise patristic distances on the tree
  dist_matrix <- cophenetic.phylo(tree)
  
  # Identify internal edges (those whose child node is not a leaf)
  internal_edges <- which(!sapply(tree$edge[,2], function(node) is_leaf_node(tree, node)))
  num_internal_edges <- length(internal_edges)
  
  # Check if enough internal edges exist for desired k
  if (num_internal_edges < (k_min - 1)) {
    stop("Not enough internal edges to perform the requested number of cuts.")
  }
  if(num_internal_edges < k_max) {
    k_max <- num_internal_edges  
  }
  
  # store best grouping for each k
  partition <- data.frame()
  # Loop over number of clusters k
  for (k in k_min:k_max) {
    cat("Processing k =", k, "\n")
    num_cuts <- k - 1
    cat("  Number of edges to cut:", num_cuts, "\n")
    cat("  Generating cut combinations...\n")
    
    # All combinations of internal edges to cut
    cut_combinations <- combn(internal_edges, num_cuts, simplify = FALSE)
    num_combinations <- length(cut_combinations)
    cat("  Total combinations:", num_combinations, "\n")
    
    # Track best CH index for this k
    best_ch_k <- -Inf
    best_grouping_k <- NULL
    best_combination_k <- NULL
    
    # Helper to process a single cut combination
    process_single_cut <- function(cut_edges, tree, dist_matrix, k) {
      group_labels <- assign_groups(tree, cut_edges, k)
      ch <- compute_ch_index(dist_matrix, group_labels)
      return(list(ch = ch, group = group_labels, cut = cut_edges))
    }
    
    # Evaluate all combinations (currently sequential)
    results <- lapply(cut_combinations, function(cut) process_single_cut(cut, tree, dist_matrix, k))
    
    # Select combination with highest CH index and correct number of groups
    for (res in results) {
      if(length(unique((res$group)))== k){
        if (res$ch > best_ch_k) {
          best_ch_k <- res$ch
          best_grouping_k <- res$group
          best_combination_k <- res$cut
        }
      }
    }
    
    # Store best grouping for this k in `partition`
    if (k>=2) {
      if (nrow(partition) == 0) {
        partition <- data.frame(matrix(ncol = 0, nrow = length(best_grouping_k)))
      }
      partition <- partition %>%
        mutate(!!sym(paste0(k, " ", "groups")) := best_grouping_k)
    }
    cat("  k =", k, "best CHindex value:", best_ch_k, "\n")
  }
  return(partition)
}

#For a given tree, find best clustering for k=2..4 using CH index and return both CH index values and best grouping.
get_CH_index <- function(rpar){
  tree <- root(rpar, outgroup = "GL", resolve.root = TRUE)
  tree <- drop.tip(tree, "GL") 
  dist_matrix <- cophenetic.phylo(tree) 
  
  # Search best cut for k=2..4
  result_simple <- find_best_k_cutree_all_combinations(tree, k_min = 2, k_max = 4,use_parallel = F)
  N <- dim(dist_matrix)[1] 
  
  # Compute CH index for each column in result_simple
  index <- apply(result_simple, 2, function(x){
    sse <- 0
    for(i in names(table(x))){
      sample_location <- which(x==i)
      sse <- sse+sum(dist_matrix[sample_location,sample_location])/length(sample_location)
    }
    K <- length(unique(x))
    ssb <- sum(dist_matrix)/N-sse
    ch <- (ssb*(N-K))/(sse*(K-1))
  })
  
  # Select the grouping with maximum CH index
  max_group <- which.max(index)
  group <- result_simple[,max_group]
  CHindex_result <- c()
  CHindex_result[["Index"]] <- index
  CHindex_result[["Genetic_group"]] <- group
  return(CHindex_result)
}

