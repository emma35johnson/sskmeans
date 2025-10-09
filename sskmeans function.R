
library(tidyverse)
library(arrangements)

sskmeans <- function(data, k, pos.eq = NULL, neg.eq = NULL, max.iter = 50) {
  
  ##########  Ensure `pos.eq` runs as intended  ################################
  
  if (is.null(pos.eq)) {
    
    data$pos.eq <- seq(nrow(data))  ## imputes vals 1, 2, ..., n
    
  } else if (!pos.eq %in% names(data)) {  ## must exist as a name in the data
    
    stop("`pos.eq` must be the name of a column in `data`")
    
  } else { ## good to go!
    
    names(data)[names(data) == pos.eq] <- "pos.eq"  ## ensures compatability when 
                                                    ## referencing `pos.eq`
  }
  
  # Impute missing block labels:
  
  data$pos.eq <- as.numeric(data$pos.eq)  ## convert to numeric for imputation
  
  pos.NAs <- which(is.na(data$pos.eq))  ## indexing the NA values
  
  if (length(pos.NAs) > 0) {  ## imputation step
    
    pos.input <- max(data$pos.eq, na.rm = TRUE)
    
    data$pos.eq[pos.NAs] <- pos.input + seq_along(pos.NAs)
    
  }
  
  data$pos.eq <- as.factor(as.integer(factor(data$pos.eq))) ## ensures no jumps; 
                                                            ## seq from 1 to B
  
  ##############################################################################
  
  ##########  Ensure `neg.eq` runs as intended  ################################
  
  if (is.null(neg.eq)) {
    
    data$neg.eq <- vector("list", nrow(data))   ## empty list --> runs pos_km
    
  } else if (!neg.eq %in% names(data)) {
    
    stop("`neg.eq` must be the name of a list column in `data`")
    
  } else { ## good to go!
    
    names(data)[names(data) == neg.eq] <- "neg.eq"  ## ensures compatability when 
                                                    ## referencing `neg.eq`
  }
  
  if (!is.list(data$neg.eq)) {
    
    data$neg.eq <- lapply(data$neg.eq, function(x) {   ## converts `neg.eq` to list
                                        if (is.na(x)) {NULL}
                                        else {x}
                                      })
  }
  
  ##############################################################################
  
  ##########  Initializing vaiables   ##########################################
  
  data_num <- select(data, where(is.numeric))   ## numeric data removes pos/neg
  
  n <- nrow(data)       ## number observations
  p <- ncol(data_num)   ## number features (x-vars)
  
  B <- max(as.integer(data$pos.eq)) ## total number of blocks
  zb <- integer(B)                  ## cluster assignments per block
  
  dist.sq <- matrix(nrow = n, ncol = k)   ## Euclidean distances squared per obs
  
  # Centroids
  centers <- data_num[sample(seq(n), k), ]  ## initializing `k` random centers
  
  # Vectors for individual cluster assignments
  cluster <- integer(n)       ## current iteration
  prev_cluster <- integer(n)  ## previous iteration
  
  # iteration counter & stopper
  iter <- 1
  stopper <- 0  ## stops when prev_cluster = cluster  -OR-  max.iter = iter
  
  ##############################################################################
  
  ##############################################################################
  ##########  While Loop starts here  ##########################################
  ##############################################################################

  while (stopper == 0) {
    
    ############################################################################
    ##########  Positive K-M  ##################################################
    ############################################################################

    ##########  Euclidean Distances for Observations  ##########################
    
    for (i in seq(n)) {
      
      distance <- data_num[i, ] %>%
        rbind(centers) %>%
        dist()  ## Euclidean Distance
      
      dist.sq[i, ] <- distance[1:k]^2   ## records the distance of the ith obs from
                                    ## the kth center
    }
    
    ########## Block cluster assignments for positive constraints ONLY`#########`
    
    for (b in 1:B) {
      
      data.b <- dist.sq[data$pos.eq == b, , drop = FALSE] ## data grouped by block
      
      if (nrow(data.b) == 1) {
        
        zb[b] <- which.min(data.b)  ## cannot perform colSums on singular row
        
      } else {
        
        zb[b] <- which.min(colSums(data.b))   ## minimal within-block SS
        
      }
      
      cluster[data$pos.eq == b] <- zb[b]  ## assign cluster ID to all obs in block
      
    }
    
    ############################################################################
    
    ############################################################################
    ##########  Unbiased MLEs for mu & sigma   #################################
    ############################################################################
    
    mu_kj <- matrix(NA, nrow = k, ncol = p)   ## matrix of mean estimates
    
    ##########  Calculating each mu_kj #########################################
    
    for (j in 1:p) {
      
      for (l in 1:k) {
        
        mu_kj[l, j] <- mean(data_num[cluster == l, j])
        
      }
    }
    
    SSE <- 0  ## initializing SSE (within-cluster)
    
    X <- as.matrix(unname(data_num))  ## copy of data_num

    ##########  Calculating sigma squared ######################################
    
    for (j in 1:p) {  ## over x-vars
      
      for (l in 1:k) {  ## over centers
        
        for (i in which(cluster == l)) {  ## over clusters
          
          SSE <- SSE + (X[i, j] - mu_kj[l, j])^2  ## WCSS
          
        }
      }
    }
    
    sigma.sq <- SSE / (p * (n - k))  # unbiased
    
    ############################################################################
    
    ############################################################################
    ##########  Negative K-M  ##################################################
    ############################################################################
    
    ##########  Euclidean Distances for Blocks  ################################
    
    dist.sq_b <- matrix(nrow = B, ncol = k)   ## initializing dist.sq_b for blocks
    
    for(b in 1:B) {
      
      data.b <- as.matrix(data_num[data$pos.eq == b, , drop = FALSE]) # block `b` data
      
      dist.b_sum <- 0 # initializing dist.b_sum
      
      for(i in 1:nrow(data.b)) {
        
        dist.b <- data.b[i,] %>%  ## similar to Euc Dist for single obs
          rbind(mu_kj) %>%
          dist()
        
        dist.b_sum <- dist.b_sum + dist.b[1:k]^2
        
      }
      
      dist.sq_b[b,] <- dist.b_sum   ## grouped by obs in block `b`
      
    }
    
    ##########  Table for obs in each Block ####################################
    
    pos_lookup <- data %>%
      mutate(ID = row_number()) %>%
      aggregate(ID ~ pos.eq, c)
    
    ##########  Table for Block relationships ##################################
    
    block_rels <- data %>%
      mutate(ID = row_number()) %>%             ## either neg is a list column -OR- IDs must be
      relocate(ID, .before = names(data)) %>%   ## repeated for each additional neg constraint
      unnest(neg.eq, keep_empty = TRUE) %>% 
      select(ID, pos.eq, neg.eq) %>%
      left_join(., unnest(pos_lookup, ID), 
                by = c("neg.eq" = "ID"), 
                relationship = "many-to-many") %>%
      select(pos.eq = pos.eq.x, neg.from = pos.eq.y) %>%
      aggregate(neg.from ~ pos.eq, data = ., c, na.action = na.pass) %>%
      mutate(neg.from = lapply(neg.from, function(x) x[!is.na(x)]))
    # NOTE: a block w/ no neg constraints is not NA, but length = 0
    
    ##########  Dealing with z_b for Negative Constraints ######################
    
    dist.sq_b[!is.finite(dist.sq_b)] <- Inf   ## precaution for large sums
    
    for (b in 1:B) {
      
      neg_blocks <- as.integer(block_rels$neg.from[[b]])  ## neg blocks for each
      neg_blocks <- neg_blocks[!is.na(neg_blocks)]        ## pos block
      
      sum_k <- numeric(k)   ## initialize each sum of the kth cluster
      
      for (l in 1:k) {
        
        if (!is.finite(dist.sq_b[b, l])) {  ## precaution for large sums
          
          sum_k[l] <- -Inf  ## will not matter since we take max k
          next
          
        }
        
        if (length(neg_blocks) == 0) {  ## for 'NA' negative constraints
          
          sum_k[l] <- exp(- dist.sq_b[b, l] / (2 * sigma.sq))   ## i.e. just sum
          next                                                  ## over kth cluster
          
        }
        
        ##########  Complex Negative Constraint Cases ##########################
        
        if (length(neg_blocks) > length(setdiff(1:k, l))) {   ## focus on kth cluster
          
          sum_k[l] <- 0   ## initializing sums for blocks w/ neg constraints
          next
          
        }
        
        if (length(neg_blocks) == 1) {
          
          perms_list <- lapply(setdiff(1:k, l), identity) ## returns (1:k)\l list
          
        } else {  ## NOTE this is necessary b/c must consider all cluster assignments
                  ## per each 'focus' cluster and each negative constraint
          
          combs <- combn(setdiff(1:k, l),    ## finds all combos of clusters
                         length(neg_blocks), ## minus the kth 'focus' cluster
                         simplify = FALSE)   ## `FALSE` indicates a list returns
          
          perms_list <- purrr::map(combs, ~ {   ## permutates all combos
            pm <- arrangements::permutations(.x)
            split(pm, seq_len(nrow(pm)))  ## writes elements as rows 
            }
            ) %>% 
            purrr::list_flatten()   ## coerces to one long list
        }
        
        ##########  Taking the sums ############################################
        
        s <- 0  ## sum[exp(-1/(2*sigma^2)*(sum dist.sq_b))]
        
        for (pl in perms_list) {
          
          pl_vec <- as.integer(pl)  ## coercing to numeric for compatibility
          
          tot <- dist.sq_b[b, l]  ## initializing sum dist.sq_b
          
          for (nb in seq_along(neg_blocks)) {   ## i.e. nb is specific neq block
            
            add <- dist.sq_b[neg_blocks[nb], pl_vec[nb]]
            
            if (!is.finite(add)) {  ## precaution for large sums
              
              tot <- Inf
              break 
              
              }
            
            tot <- tot + add  ## adds each neg Euc Dist (inner sums)
          }
          
          s <- s + exp(- tot / (2 * sigma.sq))  ## outer sums
        }
        
        sum_k[l] <- s
      }
      
      zb[b] <- which.max(sum_k)
      
      cluster[data$pos.eq == b] <- zb[b]  ## final cluster assignments
    }
    
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
  
  sizes <- data %>%   ## each n_b
    mutate(cluster = cluster) %>%   ## necessary since `cluster` has changed
    count(cluster) %>%
    pull(n)
  
  centers_out <- data_num %>%   ## final centers / centroids
    mutate(cluster = cluster) %>%
    group_by(cluster) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop") %>%
    as.data.frame()
  
  result <- list(
    "size" = sizes,
    "cluster" = cluster,
    "centers" = centers_out,
    "iterations" = iter,
    "blocks" = pos_lookup,
    "relationships" = block_rels
    
  )
  
  return(result)
}



################################################################################
##########  Example Usage ######################################################
################################################################################

## dummy data
df <- data.frame(
  x1 = c(2, 9, 8, 4.5, 3, 6, 8, 11),
  x2 = c(5, 7, 11, 4, 3, 9, 8, 8)
)
df$pos <- c(1, NA, 2, NA, 1, 2, 2, NA)
df$neg <- list(6, 4, NULL, c(2, 6), NULL, c(1, 4), NULL, NULL)
df

## with both pos and neg constraints
sskmeans(data = df, k = 3, pos.eq = 'pos', neg.eq = 'neg')

## with ONLY pos constraints
sskmeans(data = df[,1:3], k = 3, pos.eq = 'pos')

## with ONLY neg constraints
sskmeans(data = df[,c(1:2,4)], k = 3, neg.eq = 'neg')

## vs. NO constraints
sskmeans(data = df[,1:2], k = 3)
# compare w/
kmeans(df[,1:2], 3)

## if pos.eq is a column not indexed at 1
df2 <- df %>%
  mutate(pos = case_when(pos == 1 ~ 7,
                         pos!= 1 ~ pos))

sskmeans(data = df2, k = 3, pos.eq = 'pos', neg.eq = 'neg') ## works the same
## NOTE that the input block labels change in the algorithm (e.g. `7` becomes `1`)

## consider adding seeds
