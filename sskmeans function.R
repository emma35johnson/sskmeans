

library(tidyverse)
library(arrangements)

## `neg.eq` is BLOCK-LEVEL adjacency matrix (i.e must start at block 1, 2, ..., B)
## might want to change function input `k` to `centers` to match standard kmeans
## consider saving the sigma squared from the best `nstart` run

sskmeans <- function(data, k, pos.eq = NULL, neg.eq = NULL, max.iter = 50, nstart = 1) {
  
  ##########  If no `neg.eq` input  ############################################
  
  pos_only_mode <- is.null(neg.eq)  ## note that this also runs regular K-Means
  ## if no `pos.eq` is specified, either
  
  data_num <- select(data, where(is.numeric))   ## numeric data ONLY
  
  n <- nrow(data_num)    ## number of observations
  p <- ncol(data_num)    ## number of features (x-vars)
  
  ##########  Ensure `pos.eq` runs as intended  ################################
  
  if (is.null(pos.eq)) {
    
    pos.eq <- as.list(1:n)  ## imputes vals 1, 2, ..., n
    
  }
  
  # Table for obs in each Block:
  
  pos_lookup <- data.frame(       ## this is just for returned results
    block = seq_along(pos.eq),    ## only includes blocks specified by user
    IDs  = I(pos.eq)              ## I() removes list indices; turns into vector
  )
  
  # Impute missing block labels:
  
  all <- 1:n    ## all possible indices
  
  used <- unlist(pos.eq)    ## indices already used in `pos.eq`
  
  unused <- setdiff(all, used)    ## indices not used in `pos.eq`
  
  pos.eq <- c(pos.eq, as.list(unused))    ## appends unused indices as their own blocks
  
  # Initializations before looping:
  
  B_neg <- nrow(neg.eq)    ## original number of blocks specified (must start at 1)
  
  B_pos <- length(pos.eq)    ## total number of blocks after imputing unused obs
  ## Note: B_neg <= B_pos
  
  zb <- integer(B_pos)  ## cluster assignments per block
  
  # Total Sum of Squares:
  
  totss <- 0  ## total sum of squares
  
  for (i in 1:n) {    ## Euc Dist squared from overall mean for each obs
    
    totss <- totss + sum((data_num[i, ] - colMeans(data_num))^2)
    
  }
  
  min_totwss <- Inf    ## keeps minimal `tot.withinss`
  
  best_result <- NULL    ## result for best `nstart` run
  
  for (nstart_placeholder in 1:nstart) {    ## ensures `nstart` number restarts
    
    dist.sq <- matrix(nrow = n, ncol = k)   ## Euclidean distances squared per obs
    
    # Centroids / Random Centers Initialization:
    
    centers <- data_num[sample(1:n, k), ]  ## initializing `k` random centers
    
    # Vectors for individual cluster assignments
    
    cluster <- integer(n)       ## current iteration
    prev_cluster <- integer(n)  ## previous iteration
    
    # Iteration counter & stopper:
    
    iter <- 1
    stopper <- 0  ## stops when prev_cluster = cluster  -OR-  max.iter = iter
    
    
    ##############################################################################
    ##########  While Loop starts here  ##########################################
    ##############################################################################
    
    while (stopper == 0) {
      
      ############################################################################
      ##########  Positive K-M  ##################################################
      ############################################################################
      
      ##########  Euclidean Distances for Observations  ##########################
      
      for (i in 1:n) {
        
        distance <- data_num[i, ] %>%
          rbind(centers) %>%
          dist()    ## Euclidean Distance
        
        dist.sq[i, ] <- distance[1:k]^2   ## records the distance of the ith obs
        ## from the kth center
      }
      
      ##########  Block cluster assignments for positive constraints ONLY  #####
      
      for (b in 1:B_pos) {
        
        data.b <- dist.sq[unlist(pos.eq[b]), , drop = FALSE] ## data grouped by block
        
        if (nrow(data.b) == 1) {
          
          zb[b] <- which.min(data.b)  ## cannot perform colSums on singular row
          
        } else {
          
          zb[b] <- which.min(colSums(data.b))   ## minimal within-block SS
          
        }
        
        cluster[unlist(pos.eq[b])] <- zb[b]  ## assign cluster ID to all obs in block
        
      }
      
      ############################################################################
      ##########  Unbiased MLEs for mu & sigma   #################################
      ############################################################################
      
      if (!pos_only_mode) {
        
        mu_kj <- matrix(NA, nrow = k, ncol = p)   ## matrix of mean estimates
        
        ##########  Calculating each mu_kj #####################################
        
        for (j in 1:p) {
          
          for (l in 1:k) {
            
            mu_kj[l, j] <- mean(data_num[cluster == l, j])
            
          }
        }
        
        SSE <- 0  ## initializing SSE (within-cluster)
        
        X <- as.matrix(unname(data_num))  ## copy of data_num
        
        ##########  Calculating sigma squared ##################################
        
        for (j in 1:p) {  ## over x-vars
          
          for (l in 1:k) {  ## over centers
            
            for (i in which(cluster == l)) {  ## over clusters
              
              SSE <- SSE + (X[i, j] - mu_kj[l, j])^2  ## WCSS
              
            }
          }
        }
        
        sigma.sq <- SSE / (p * (n - k))  # unbiased
        
        ##########################################################################
        ##########  Negative K-M  ################################################
        ##########################################################################
        
        block_rels <- data.frame(    ## used for returned results & summation indices
          pos.eq   = 1:B_neg,
          neg.from = I(
            lapply(1:B_neg, function(b) {which(neg.eq[b, ] == 1)} )
          )
        )
        
        ##########  Euclidean Distances for Blocks  ##############################
        
        dist.sq_b <- matrix(nrow = B_neg, ncol = k)   ## initializing dist.sq_b for blocks
        
        for(b in 1:B_neg) {
          
          data.b <- as.matrix(data_num[unlist(pos.eq[b]), , drop = FALSE]) ## block `b` data
          
          dist.b_sum <- 0 ## initializing dist.b_sum
          
          for(i in 1:nrow(data.b)) {
            
            dist.b <- data.b[i,] %>%  ## similar to Euc Dist for single obs
              rbind(mu_kj) %>%
              dist()
            
            dist.b_sum <- dist.b_sum + dist.b[1:k]^2
            
          }
          
          dist.sq_b[b,] <- dist.b_sum   ## grouped by obs in block `b`
          
        }
        
        # Using each `dist.sq_b` to compute weights `w[b, l]` (to be added for z_b)
        
        w <- exp( -dist.sq_b / (2 * sigma.sq) ) ## w[b, l] = exp(-(dist.sq_b[b, ]) / (2 sigma^2))
        w[!is.finite(w)] <- 0    ## NaN as 0
        
        ##########  Block cluster assignments with negative constraints ########
        
        for (b in 1:B_neg) {
          
          neg_blocks <- as.integer(block_rels$neg.from[[b]])
          neg_blocks <- neg_blocks[!is.na(neg_blocks)]
          
          m <- length(neg_blocks)
          
          # Initializing objective for `z_b`:
          sum_k <- numeric(k)    ## max(sum_k) implies z_b assignment
          
          for (l in 1:k) {
            
            if (m == 0) {    ## no neg constraints for block `b`
              sum_k[l] <- w[b, l]
              next
            }
            
            ## clusters currently used by negatively related blocks
            neg_block_clusters <- unique(cluster[unlist(pos.eq[neg_blocks])])
            
            if (any(neg_block_clusters == l)) { ## if cluster `l` is taken by a neg block
              sum_k[l] <- 0
              next
            }
            
            if (w[b, l] == 0 || !is.finite(w[b, l])) { ## if weight is zero or NaN
              sum_k[l] <- 0
              next
            }
            
            # Calculating objective for `z_b`:
            minus_k <- setdiff(1:k, l)    ## all clusters besides focus `l`
            
            # Allow reuse of clusters if m > length(minus_k)
            use_replace <- m > length(minus_k)
            
            perms_mat <- arrangements::permutations(
              x = minus_k,
              k = m,
              replace = use_replace
            )
            
            # Defensive: if for some reason we got no permutations, treat as 0 contribution
            if (nrow(perms_mat) == 0) {
              sum_k[l] <- 0
              next
            }
            
            s <- 0  ## initializing summation for `z_b`
            
            for (i in 1:nrow(perms_mat)) {
              
              cols <- perms_mat[i, ]              ## cluster labels directly
              vals <- w[cbind(neg_blocks, cols)]   ## vector of weights for neg blocks
              
              s <- s + w[b, l] * prod(vals)        ## contribution to summation
            }
            
            sum_k[l] <- s
          }
          
          zb[b] <- which.max(sum_k)
          
          cluster[unlist(pos.eq[b])] <- zb[b]
        }
        
      }   ## end of `if (!pos_only_mode)`
      
      ############################################################################
      ##########  After each Iteration ###########################################
      ############################################################################
      
      if (all(cluster == prev_cluster)) {   ## ends on cluster convergence
        
        stopper <- 1
        
      }
      
      prev_cluster <- cluster   ## assigning last iteration to `prev_cluster`
      
      ##########  Recompute centers ##############################################
      
      centers <- data_num %>%
        mutate(cluster = cluster) %>%   ## necessary since `cluster` has changed
        group_by(cluster) %>%
        summarise(across(where(is.numeric), mean), .groups = "drop") %>%
        select(-cluster) %>%
        as.data.frame()
      
      ##########  Begin next Iteration  ##########################################
      
      iter <- iter + 1
      
      if (iter > max.iter) {  ## maximum iteration occurred before convergence
        
        message("Reached maximum iterations before convergence.")
        stopper <- 1
        
      }
      
      ############################################################################
    } ########## End `while` #####################################################
    ############################################################################
    
    ##########  Outputs ##########################################################
    
    sizes <- tabulate(cluster, nbins = k)    ## each n_b
    
    centers_out <- data_num %>%   ## final centers / centroids
      mutate(cluster = cluster) %>%
      group_by(cluster) %>%
      summarise(across(where(is.numeric), mean), .groups = "drop") %>%
      as.data.frame()
    
    ##########  Record best WCSS & result  #####################################
    
    WCSS <- numeric(k)  ## initializing WCSS
    
    X <- as.matrix(unname(data_num))  ## copy of data_num
    
    for (j in 1:p) {  ## over x-vars
      
      for (l in 1:k) {  ## over centers
        
        for (i in which(cluster == l)) {  ## over clusters
          
          WCSS[l] <- WCSS[l] + (X[i, j] - centers_out[l, j + 1])^2  ## WCSS
          
        }
      }
    }
    
    tot.withinss <- sum(WCSS)
    
    result <- list(
      "cluster" = cluster,
      "centers" = centers_out,
      "totss" = totss,
      "withinss" = WCSS,
      "tot.withinss" = tot.withinss,
      "betweenss"  = totss - tot.withinss,
      "size" = sizes,
      "iterations" = iter,
      "blocks" = pos_lookup
    )
    
    if (!pos_only_mode) {
      result$relationships = block_rels
    }
    
    if (tot.withinss < min_totwss) { ## keeps best result across `nstart`
      min_totwss <- tot.withinss
      best_result <- result
    }
    
  } ## end `for (nstart_placeholder in 1:nstart)`
  
  return(best_result)
}


################################################################################
##########  Example Usage ######################################################
################################################################################

## dummy data
df <- data.frame(
  x1 = c(2, 9, 8, 4.5, 3, 6, 8, 11),
  x2 = c(5, 7, 11, 4, 3, 9, 8, 8)
)

df %>%
  ggplot(aes(x1, x2)) +
  geom_point(size = 3) +
  geom_text(aes(label = 1:nrow(df)), vjust = -1) +
  theme_classic()

## positive constraints as a list:
## indices are blocks, elements are observations in that block
pos.eq <- list(c(1, 5), c(3, 6, 7), 2, 4)

## using adjacency matrix to represent the negative constraints
## i.e. each negative constraint in an edge coded `1`; no edges are represented by `0`

## could include a separate function in the package that creates this matrix:
B <- 4 # number of blocks INVOLVED (i.e. if points are in no constraints, do not include)
neg.eq <- matrix(data = rep(0, B^2), ## matrix of all zeros to start
                 nrow = B,
                 ncol = B)
neg.eq[1, 2] <- 1   ## 1 means a negative constraint exists between blocks
neg.eq[1, 3] <- 0   ## change to 1 to make cyclic
neg.eq[2, 3] <- 0   ## change to 1 to make cyclic
neg.eq[2, 4] <- 1
neg.eq[3, 4] <- 1
neg.eq[lower.tri(neg.eq)] <- t(neg.eq)[lower.tri(neg.eq)]

sskmeans(data = df, k = 3, pos.eq = pos.eq, neg.eq = neg.eq, nstart = 20)

## could use upper.tri() or lower.tri() in fxn to cut the searches in half

library(igraph)
graph_from_adjacency_matrix(neg.eq, mode = "undirected") %>% plot()
