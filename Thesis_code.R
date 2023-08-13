# Master Thesis project


#bayesian I optimality
#Information matrix




#Coordiante exchange + clustering



#

library(purrr)
library(ggplot2)
library(opdesmixr)
library(tidyverse)
library(here)

# To ensure exact reproducibility, make sure to install the same version of the opdesmixr package that was used to run this code and get the results of the paper:
devtools::install_github("mariobecerra/opdesmixr", ref = "777e67ffdd406c3b7a21952313e9aeb815cb4174")
#devtools::install_github("mariobecerra/opdesmixr", ref = "777e67ffdd406c3b7a21952313e9aeb815cb4174")
#designs_folder = here("out/cocktail_cornell_designs/")
#dir.create(designs_folder, showWarnings = F)
#n_cores = parallel::detectCores()


#1. Random design

mnl_create_random_initial_design = function(q, J, S, seed = NULL){
  X = array(rep(NA_real_, (q)*J*S), dim = c((q), J, S))

  if(!is.null(seed)) set.seed(seed)

  for(j in 1:J){
    for(s in 1:S){
      rands = runif(q)
      # mixture ingredients must sum up to 1
      X[1:q,j, s] = rands/sum(rands)
    }
  }

  return(X)
}

mnl_create_random_initial_design(3,2,7)

# choice probablities
MNL_compute_choice_probabilites = function(design, beta, order) {
  #
  q = dim(design)[1]
  J = dim(design)[2]
  S = dim(design)[3]


  p1 = q - 1
  p2 = q * (q - 1)/2
  p3 = q * (q - 1) * (q - 2)/6

  if (order == 1) {
    n_parameters = p1
  } else {
    if (order == 2) {
      n_parameters = p1 + p2
    } else {
      if (order == 3) {
        n_parameters = p1 + p2 + p3
      }
    }
  }



  # beta = rep(1,)
  beta_star = beta[1:p1]  # modify
  beta_2FI = beta[(p1 + 1):p2]
  beta_3FI = beta[(p2 + 1):p3]

  # stopifnot(length(beta) == n_parameters) stopifnot(order %in% 1:3)

  # define utility (matrix) of alternative j in choice set s for (i in 1:q){
  U = matrix(rep(NA_real_, J * S), ncol = S)
  for (j in 1:J) {
    for (s in 1:S) {

      # ORDER-1 X_JS

      x_js = design[, j, s]
      U_js_term1 = sum(x_js[1:p1] * beta_star)
      u_js = U_js_term1


      if (order >= 2) {
        # ORDER-2 X_IJS*X_KJS - 2FI terms
        x_second_order = rep(NA_real_, p2)
        counter = 0
        for (i in 1:(q - 1)) {
          for (k in (i + 1):q) {
            counter = counter + 1
            x_ijs_x_kjs = design[i, j, s] * design[k, j, s]
            x_second_order[counter] = x_ijs_x_kjs
          }
        }
        U_js_term2 = sum(x_second_order * beta_2FI)
        u_js = u_js + U_js_term2
      }
      # x_ijs_x_kjs = design[i,j,s]*design[k,j,s]
      if (order >= 3) {
        # 0RDER-3 X_IJS*X_KJS*XLJS - 3FI terms
        x_third_order = rep(NA_real_, p3)
        counter = 0
        # i = 1 k = i+1 l = k+1
        for (i in 1:(q - 2)) {
          for (k in (i + 1):(q - 1)) {
            for (l in (k + 1):q) {
              counter = counter + 1
              x_ijs_x_kjs_xljs = design[i, j, s] * design[k, j, s] * design[l,
                                                                            j, s]
              x_third_order[counter] = x_ijs_x_kjs_xljs
            }
          }
        }
        U_js_term3 = sum(x_third_order * beta_3FI)
        u_js = u_js + U_js_term3
      }
      # x_ijs_x_kjs_xljs = design[i,j,s]*design[k,j,s]*design[l,j,s]


      # U_js_term2 = sum(x_ijs_x_kjs[(p1+1):p2]*beta_2FI) U_js_term2 =
      # sum(x_second_order*beta_2FI) U_js_term3 = sum(x_ijs_x_kjs[(p2+1):p3]*beta_3FI)
      # U_js_term3 = sum(x_third_order*beta_3FI) u_js = U_js_term1 + U_js_term2 +
      # U_js_term3
      U[j, s] = u_js  # matrix of model expansion
    }
  }

  # define choice probability of alternative j in choice s
  P = matrix(rep(NA_real_, J * S), ncol = S)
  for (j in 1:J) {
    for (s in 1:S) {
      u_js = U[j, s]
      p_js = exp(u_js)/sum(exp(U[, s]))
      P[j, s] = p_js  # choice probabilities form matrix of model expansion
    }
  }
  return(P)
}


#model matrix
MNL_compute_model_matrix <- function(design, order) {
  q = dim(design)[1]
  J = dim(design)[2]
  S = dim(design)[3]


  p1 = q - 1
  p2 = q * (q - 1)/2
  p3 = q * (q - 1) * (q - 2)/6

  if (order == 1) {
    n_parameters = p1
  } else {
    if (order == 2) {
      n_parameters = p1 + p2
    } else {
      if (order == 3) {
        n_parameters = p1 + p2 + p3
      }
    }
  }

  model_array = array(rep(NA_real_, n_parameters * J * S), dim = c(n_parameters,
                                                                   J, S))

  if (order >= 1) {

    model_array[1:(q - 1), , ] = design[1:(q - 1), , ]
  }

  if (order >= 2) {
    counter = q - 1
    for (i in 1:(q - 1)) {
      for (j in (i + 1):q) {
        counter = counter + 1
        model_array[counter, , ] = design[i, , ] * design[j, , ]
      }
    }
  }
  if (order >= 3) {
    # counter = 0
    for (i in 1:(q - 2)) {
      for (j in (i + 1):(q - 1)) {
        for (k in (j + 1):q) {
          counter = counter + 1
          model_array[counter, , ] = design[i, , ] * design[j, , ] * design[k,
                                                                            , ]
        }
      }
    }
  }

  return(model_array)

}

#information matrix
MNL_compute_information_matrix = function(design, beta, order) {
  q = dim(design)[1]
  J = dim(design)[2]
  S = dim(design)[3]


  p1 = q - 1
  p2 = q * (q - 1)/2
  p3 = q * (q - 1) * (q - 2)/6

  if (order == 1) {
    n_parameters = p1
  } else {
    if (order == 2) {
      n_parameters = p1 + p2
    } else {
      if (order == 3) {
        n_parameters = p1 + p2 + p3
      }
    }
  }

  model_array = MNL_compute_model_matrix(design, order)
  P = MNL_compute_choice_probabilites(design, beta, order)

  information_matrix = matrix(rep(0, n_parameters * n_parameters), ncol = n_parameters)  # generalize this
  for (s in 1:S) {
    p_s = P[, s]
    I_s = model_array[1:n_parameters, , s] %*% (diag(p_s) - (t(t(p_s)) %*% t(p_s))) %*%
      t(model_array[1:n_parameters, , s])
    information_matrix = information_matrix + I_s
  }
  return(information_matrix)
}


#moment matrix
MNL_create_moment_matrix = function(q, order) {

  # q=ncol(design) # disregard, make it simplier define the number of parameters
  # for model order \in 1, 2, 3

  if (order == 1) {
    parameters = q - 1
  } else {
    if (order == 2) {
      parameters = q - 1 + q * (q - 1)/2
    } else {
      if (order == 3) {
        parameters = q - 1 + q * (q - 1)/2 + q * (q - 1) * (q - 2)/6
      }
    }
  }

  # Auxiliary matrix F is of size parameters with q columns; We fill it in with 1
  # if the parameter is assumed to be in the moment matrix

  auxiliary_matrix_f = matrix(rep(0, parameters * q), ncol = q)
  counter = 0
  for (i in 1:(q - 1)) {
    counter = counter + 1
    auxiliary_matrix_f[counter, i] = 1
  }

  # Indeces for order 2: \sum_{i}^{q-1}\sum_{j=i+1}^{q} j=i+1 and i=1

  if (order >= 2) {
    for (i in 1:(q - 1)) {
      for (j in (i + 1):q) {
        counter = counter + 1
        auxiliary_matrix_f[counter, i] = 1
        auxiliary_matrix_f[counter, j] = 1
      }
    }
  }

  # Indeces for order 3: \sum_{i}^{q-2}\sum_{j=i+1}^{q-1}\sum_{k=j+1}^{q} i=1,
  # j=i+1, k=j+1 each ith, jth, kth element are filled in with 1

  if (order >= 3) {
    for (i in 1:(q - 2)) {
      for (j in (i + 1):(q - 1)) {
        for (k in (j + 1):q) {
          counter = counter + 1
          auxiliary_matrix_f[counter, i] = 1
          auxiliary_matrix_f[counter, j] = 1
          auxiliary_matrix_f[counter, k] = 1
        }
      }
    }
  }

  MNL_moment_matrix = matrix(rep(NA_real_, parameters * parameters), ncol = parameters)

  for (i in 1:parameters) {
    for (j in 1:parameters) {

      # auxiliary_vector is the sum of ith and jth row of auxiliary matrix F

      auxiliary_vector = auxiliary_matrix_f[i, ] + auxiliary_matrix_f[j, ]
      numerator = prod(factorial(auxiliary_vector))
      denominator = factorial(q - 1 + sum(auxiliary_vector))
      MNL_moment_matrix[i, j] = numerator/denominator
    }
  }
  return(MNL_moment_matrix)
}


# I optimlaity
MNL_compute_i_optimality = function(design, order, beta) {
  q = dim(design)[1]
  printed_information_matrix = MNL_compute_information_matrix(design, beta, order)
  printed_moments_matrix = MNL_create_moment_matrix(q, order)
  # MNL_i_opt = inv(printed_information_matrix) %*% printed_moments_matrix
  MNL_i_opt = sum(diag(solve(printed_information_matrix, printed_moments_matrix)))  # trace is the sum of the diagonal elements
  return(MNL_i_opt)
}


#
#plot function
library(ggplot2)
library(ggtern)
plot_ternary_design <- function(design) {

  # Flatten design points
  ingredient1 <- as.vector(abs(design[1,,]))
  ingredient2 <- as.vector(abs(design[2,,]))
  ingredient3 <- as.vector(abs(design[3,,]))

  # Create data frame
  df <- data.frame(Ingredient1 = ingredient1, Ingredient2 = ingredient2, Ingredient3 = ingredient3)

  # Create plot object
  plot <- ggtern(data = df, aes(x = Ingredient1, y = Ingredient2, z = Ingredient3)) +
    geom_point(color = "blue", size = 2, shape = 16, alpha = 0.7) +
    theme_bw() +
    labs(title = "Design", x = "Ingredient 1", y = "Ingredient 2", z = "Ingredient 3")

  return(plot)
}



#clustering function
clustering <- function(random_design, k) {
  random_design_2d <- t(array(random_design, dim = c(3, 2 * 7)))

  # Compute distance
  distance <- dist(random_design_2d)

  # Perform hierarchical clustering
  clusters <- hclust(distance, method='average')

  # Extract the labels of each data point based on the number of clusters
  labels <- cutree(clusters, k)

  # Extract the coordinates of the data points for each cluster
  coords <- lapply(1:k, function(i) random_design_2d[labels == i, ])

  # Calculate the mean of each cluster
  means <- list()
  for (cluster in coords) {
    # Compute the mean of each dimension separately and append to means list
    if (length(cluster) == 3) {
      means <- c(means, list(cluster))
    } else {
      mean_x <- mean(cluster[, 1])
      mean_y <- mean(cluster[, 2])
      mean_z <- 1 - (mean_x + mean_y)
      means <- c(means, list(c(mean_x, mean_y, mean_z)))
    }
  }

  # Stack the mean coordinates of each cluster and transpose the result
  clustered_design = matrix(unlist(means), nrow = 3, ncol = k, byrow = FALSE)
  # Return the clustered design
  return(clustered_design)
}



#example

#initial design points
random_design=mnl_create_random_initial_design(3,2,7)
random_design
#clustering
clustered_design=clustering(random_design,k=12)
#reformat to dim(3,2,6)
clustered_design_reformat=array(clustered_design,dim=c(3,2,6))

beta=c(1,1,1,1,1,1,1)
plot_ternary_design(random_design)
plot_ternary_design(clustered_design_reformat)


clustered_design_reformat
#compute I optimality of the design
MNL_compute_i_optimality(array(random_design,dim=c(3,2,7)),3,beta)
MNL_compute_i_optimality(clustered_design_reformat,3,beta)


mnl_get_opt_crit_value(clustered_design_reformat,beta,3, "I")


#elbow curve

k=13
k%/%2
# Assume we have a range of k values to evaluate.
k_values <- 7:13

# Placeholder for the i-optimality values for each k.
i_opt_values <- numeric(length(k_values))

# For each k, perform the clustering and compute the I-optimality.
for (i in seq_along(k_values)) {
  k <- k_values[i]

  # Perform the clustering.
  clustered_design <- clustering(random_design, k)

  # Reformat the clustered_design
  clustered_design_reformat <- array(clustered_design, dim = c(3, 2, k%/%2))

  # Compute the I-optimality.
  i_opt_values[i] <- mnl_get_opt_crit_value(clustered_design_reformat,beta,3, "I")
}

# Combine the k_values and i_opt_values into a data frame for plotting.
elbow_data <- data.frame(k = k_values, i_optimality = i_opt_values)
# Use ggplot2 to create the elbow plot.
library(ggplot2)
ggplot(elbow_data, aes(x = k, y = i_optimality)) +
  geom_line() +
  geom_point() +
  xlab("Number of clusters (k)") +
  ylab("I-Optimality") +
  ggtitle("Elbow Curve for Clustering")










#VNS
# RcppArmadillo and arrangements packages are needed for this code
# install.packages("RcppArmadillo")
# install.packages("arrangements")

library(RcppArmadillo)
library(arrangements)

generate_simplex_lattice_design <- function(q, m) {
  # Ensure that q and m are positive
  if(q <= 0 || m <= 0){
    stop("values should be positive integer")
  }

  # Generate levels
  levels <- seq(0, 1, length.out = m + 1)

  # Generate points
  points <- as.matrix(expand.grid(rep(list(levels), q)))

  # Filter points that sum to 1
  design <- points[abs(rowSums(points) - 1) < 1e-9,]

  return(design)
}

###unique rows
unique_rows <- function(design) {
  # Get dimensions of the array
  q <- dim(design)[1]
  j <- dim(design)[2]
  s <- dim(design)[3]

  # Reshape the array and find unique rows
  arr <- array(t(design), c(j * s, q))
  unique_rows <- unique(arr)

  return(unique_rows)
}

#intial design and other points
lattice_design=generate_simplex_lattice_design(3,5)
lattice_design
n <- nrow(lattice_design)
n
k=10
j=2
s=7
lattice_indices <- sample(n, size=k, replace=FALSE)
lattice_points <- lattice_design[lattice_indices, ]
other_points <- lattice_design[-lattice_indices, ]
# Choose two random points from the lattice_design
set.seed(123)  # for reproducibility
random_indices <- sample(1:nrow(lattice_points), 2)
# Duplicate these points
duplicated_points <- lattice_points[random_indices,]
# Combine with original design to create new design
initial_design <- rbind(lattice_points, duplicated_points)
initial_design
#############

#VNS Functions
# Helper function: equivalent to numpy's unique rows for 3D arrays
unique_rows <- function(design) {
  q <- dim(design)[1]
  j <- dim(design)[2]
  s <- dim(design)[3]
  arr <- array(design, dim = c(q, j * s))
  return(unique(arr, MARGIN = 2))
}



# VNS First Neighborhood
initial_design_reformat=array(t(initial_design),dim=c(3,2,7))
initial_design=initial_design_reformat

alternatives <- dim(initial_design)[2]
choice <-dim(initial_design)[3]

i_opt_value <- mnl_get_opt_crit_value(initial_design_reformat, beta, order=3,"I")
rows <- unique_rows(initial_design_reformat)
i_opt_value
#######

choice
alternatives
nrow(rows)
initial_design
initial_design[,1,1]
rows[1,]
unique_points=rows[1,]
unique_points
candidate_design[, j, s]
######

improvement <- TRUE
while (improvement) {
  improvement <- FALSE
  for (s in 1:choice) {
    if (improvement) {
      break
    }
    for (j in 1:alternatives) {
      if (improvement) {
        break
      }
      candidate_design <- initial_design
      current_mix <- initial_design[, j, s]

      for (k in 1:nrow(rows)) {
        unique_points <- rows[,k]
        if (length(unique_points) != length(current_mix) || any(unique_points != current_mix)) {
          candidate_design <- initial_design
          candidate_design[, j, s] <- unique_points

          i_new_value <- NA
          try({
            i_new_value <- mnl_get_opt_crit_value(candidate_design, beta, order=3,"I")
          }, silent = TRUE)

          if (is.na(i_new_value)) {
            print("Singular matrix!")
            next
          }

          if (i_opt_value >= i_new_value && i_new_value > 0) {
            initial_design <- candidate_design
            i_opt_value <- i_new_value
            improvement <- TRUE
          }
        }
      }
    }
  }
}


######### VNS Second neigbhorhood
i_opt_value <- mnl_get_opt_crit_value(initial_design, beta, order=3,"I")
improvement <- FALSE

for (choice_idx in 1:choice) {
  if (improvement) {
    break
  }
  for (alternative_idx in 1:alternatives) {
    if (improvement) {
      break
    }
    for (other_choice_idx in 1:choice) {
      if (improvement) {
        break
      }
      for (other_alternative_idx in 1:alternatives) {
        if (choice_idx == other_choice_idx) {
          next
        }
        candidate_design <- initial_design
        candidate_design[, alternative_idx, choice_idx] <- initial_design[, other_alternative_idx, other_choice_idx]
        candidate_design[, other_alternative_idx, other_choice_idx] <- initial_design[, alternative_idx, choice_idx]

        i_new_value <- NA
        try({
          i_new_value <- mnl_get_opt_crit_value(candidate_design, beta, order=3,"I")
        }, silent = TRUE)

        if (is.na(i_new_value)) {
          print("Singular matrix!")
          next
        }

        # improvement with a minimum scale of 0.01
        if (i_opt_value >= (i_new_value + 0.01) && i_new_value > 0) {
          initial_design <- candidate_design
          i_opt_value <- i_new_value
          improvement <- TRUE
        }
      }
    }
  }
}

i_opt_value
####  VNS third neighborhood
i_opt_value <- mnl_get_opt_crit_value(initial_design, beta, order=3,"I")

design_points <- unique_rows(initial_design)
q <- dim(design_points)[2]

for (i in 1:nrow(design_points)) {
  improvement <- FALSE

  # Find indices where design points match
  indices <- which(apply(aperm(initial_design, c(2,1,3)), 1, function(x) all(x == design_points[i,])))

  # convert linear indices to array indices
  array_indices <- arrayInd(indices, .dim = dim(initial_design)[2:3])

  for (j in 1:length(other_points)) {
    candidate_design <- initial_design
    candidate_design[, array_indices[,1], array_indices[,2]] <- matrix(rep(other_points[j], each = length(indices)), nrow = q)

    i_new_value <- NA
    try({
      i_new_value <- mnl_get_opt_crit_value(candidate_design, beta, order=3,"I")
    }, silent = TRUE)

    if (is.na(i_new_value)) {
      print("Singular matrix!")
      next
    }

    if (i_opt_value >= i_new_value && i_new_value > 0) {
      initial_design <- candidate_design
      i_opt_value <- i_new_value
      other_points[j] <- design_points[i, ]
      improvement <- TRUE
      break
    }
  }

  if (improvement) {
    break
  }
}

i_opt_value






##########################Artificial Sweetener example

beta=c(1,0.86,0.21,3.07,2.34,3.24,-20.59)
#import I-optimal design
kappa0_5_Iopt <- read_csv("kappa0.5_Iopt.csv")
kappa0_5_Iopt=as.matrix(kappa0_5_Iopt)
kappa0_5_Iopt_reformat=array(kappa0_5_Iopt,dim=c(3,2,7))
kappa0_5_Iopt_reformat

z=c(0, 1.0, 0.0,0.61, 0.39, 0.00, 0.0, 1.0, 0.0, 0.24, 0.30, 0.47,
  0.0, 0.0, 1, 0, 0.67, 0.33, 1, 0, 0, 0.4, 0.6, 0.0, 0.22, 0.44,
  0.33, 0.59, 0, 0.41, 0.4, 0.0, 0.6, 1.0, 0, 0, 0, 0.33, 0.67,
  0.48, 0.25, 0.27)

kappa0_5_Iopt_reformat=array(z,dim=c(3,2,7))

#clustering
clustered_design=clustering(kappa0_5_Iopt_reformat,k=12)
#reformat to dim(3,2,6)
clustered_design_reformat=array(clustered_design,dim=c(3,2,6))

#plot
plot_ternary_design(kappa0_5_Iopt_reformat)
plot_ternary_design(clustered_design_reformat)
#I-optimality
mnl_get_opt_crit_value(clustered_design_reformat, beta, order=3,"I")


#elbow curve
k=13
k%/%2
# Assume we have a range of k values to evaluate.
k_values <- 7:13

# Placeholder for the i-optimality values for each k.
i_opt_values <- numeric(length(k_values))

# For each k, perform the clustering and compute the I-optimality.
for (i in seq_along(k_values)) {
  k <- k_values[i]

  # Perform the clustering.
  clustered_design <- clustering(kappa0_5_Iopt_reformat, k)

  # Reformat the clustered_design
  clustered_design_reformat <- array(clustered_design, dim = c(3, 2, k%/%2))

  # Compute the I-optimality.
  i_opt_values[i] <- mnl_get_opt_crit_value(clustered_design_reformat,beta,3, "I")
}

i_opt_values
# Combine the k_values and i_opt_values into a data frame for plotting.
elbow_data <- data.frame(k = k_values, i_optimality = i_opt_values)
# Use ggplot2 to create the elbow plot.
library(ggplot2)
ggplot(elbow_data, aes(x = k, y = i_optimality)) +
  geom_line() +
  geom_point() +
  xlab("Number of clusters (k)") +
  ylab("I-Optimality") +
  ggtitle("Elbow Curve for Clustering")



