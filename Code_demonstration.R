#get_beta_parameter
#intersect_term
#utility
#choice_probablity
#model_matrix


#design
original_design=array(c(
  0.00, 0.39, 0.00, 0.02, 1.0, 0.1, 0.00,
  0.52, 0.00, 0.39, 0.52, 0.17, 0.06, 0.0,
  0.0, 0.34, 1.00, 0.0, 0.0, 0.70, 0.00,
  0.00, 0.00, 0.34, 0.48, 0.15, 0.0, 0.51,
  1.0, 0.27, 0.00, 0.98, 0.0, 0.21, 1.00,
  0.48, 1.0, 0.27, 0.0, 0.68, 0.94, 0.49
), dim = c(7, 2, 3))
#design[1,1,] x_111,x_211,x311  for the first alternative  in the first choice set
#design[#numberofchoiceset,#numberofAlternative,#ingredient]
original_design[1,2,]
X_i11=original_design[,1,1]
#transform design
design <- aperm(original_design, c(3, 2, 1))

#beta
beta=c( 0.86      ,  -0.15762783,   2.46624476,   1.56402849,
        2.28214613, -21.61502097)


#################################################
#get_beta_parameter
get_parameters <- function(q, order) {
  # Calculate the number of parameters for linear, quadratic, and cubic effects
  p1 <- q - 1
  p2 <- q * (q - 1) %/% 2
  p3 <- q * (q - 1) * (q - 2) / 6

  # Return the appropriate number of parameters based on the order
  if (order == 1) {
    return(list(p1, 0, 0))
  } else if (order == 2) {
    return(list(p1, p2, 0))
  } else {
    return(list(p1, p2, p3))
  }
}
list=get_parameters(q=3,order=3)
p1=list[[1]]
p2=list[[2]]
p3=list[[3]]
########################################################
get_beta_coefficients <- function(beta, q, order) {
  params <- get_parameters(q, order)
  p1 <- params[[1]]
  p2 <- params[[2]]
  p3 <- params[[3]]

  if (length(beta) != p1 + p2 + p3) {
    stop("Number of beta coefficients does not match the number of parameters for the given order and number of ingredients.")
  }

  beta_star <- beta[1:p1]
  beta_2FI <- if (order >= 2) beta[(p1 + 1):(p1 + p2)] else numeric(0)
  beta_3FI <- if (order == 3) beta[(p1 + p2 + 1):(p1 + p2 + p3)] else numeric(0)

  return(list(beta_star, beta_2FI, beta_3FI))
}



list2=get_beta_coefficients(beta, 3, 3)
beta_star=list2[[1]]
beta_2FI=list2[[2]]
beta_3FI=list2[[3]]

######################################
#Choice probability and utility

#
design
design[1,1,]
design[1,2,]

#numbrt of beta parameters
compute_parameters <- function(q, order) {
  p1 <- q - 1
  p2 <- q * (q - 1) / 2
  p3 <- q * (q - 1) * (q - 2) / 6

  switch(order,
         '1' = p1,
         '2' = p1 + p2,
         '3' = p1 + p2 + p3,
         stop("Invalid order")
  )
}


# Compute choice probabilities using Multinomial Logit (MNL) model
# Function to compute utility
compute_utility <- function(x_js, beta_star, beta_2FI, beta_3FI, order) {
  U_js <- sum(x_js[1:length(beta_star)] * beta_star)

  if (order >= 2) {
    x_second_order <- colSums(outer(x_js, x_js) * lower.tri(matrix(1, length(x_js), length(x_js)), diag = FALSE))
    U_js <- U_js + sum(x_second_order * beta_2FI)
  }

  if (order == 3) {
    x_third_order <- apply(combn(x_js, 3, prod), 2, sum)
    U_js <- U_js + sum(x_third_order * beta_3FI)
  }

  return(U_js)
}

MNL_compute_choice_probabilities <- function(design, beta, order) {
  q <- dim(design)[1]
  J <- dim(design)[2]
  S <- dim(design)[3]

  n_parameters <- compute_parameters(q, order)

  beta_star <- beta[1:(q - 1)]
  beta_2FI <- beta[(q):(q + n_parameters - q)]
  beta_3FI <- beta[(q + n_parameters - q + 1):n_parameters]

  U <- matrix(0, J, S)
  for (j in 1:J) {
    for (s in 1:S) {
      x_js <- design[, j, s]
      U[j, s] <- compute_utility(x_js, beta_star, beta_2FI, beta_3FI, order)
    }
  }

  # Compute choice probabilities using matrix operations
  exp_U <- exp(U)
  P <- exp_U / rowSums(exp_U)

  return(P)
}



MNL_compute_choice_probabilites(design,beta,order=3)






#moment matrix



MNL_create_moment_matrix <- function(q, order) {

  parameters <- compute_parameters(q, order)

  # Initialize the auxiliary matrix with zeros
  auxiliary_matrix_f <- matrix(0, nrow = parameters, ncol = q)

  counter <- 0

  # First order terms
  for (i in 1:(q - 1)) {
    counter <- counter + 1
    auxiliary_matrix_f[counter, i] <- 1
  }

  # Second order terms
  if (order >= 2) {
    for (i in 1:(q - 1)) {
      for (j in (i + 1):q) {
        counter <- counter + 1
        auxiliary_matrix_f[counter, c(i, j)] <- 1
      }
    }
  }

  # Third order terms
  if (order >= 3) {
    for (i in 1:(q - 2)) {
      for (j in (i + 1):(q - 1)) {
        for (k in (j + 1):q) {
          counter <- counter + 1
          auxiliary_matrix_f[counter, c(i, j, k)] <- 1
        }
      }
    }
  }

  MNL_moment_matrix <- matrix(NA_real_, nrow = parameters, ncol = parameters)

  for (i in 1:parameters) {
    for (j in 1:parameters) {
      auxiliary_vector <- auxiliary_matrix_f[i, ] + auxiliary_matrix_f[j, ]
      numerator <- prod(factorial(auxiliary_vector))
      denominator <- factorial(q - 1 + sum(auxiliary_vector))
      MNL_moment_matrix[i, j] <- numerator / denominator
    }
  }

  return(MNL_moment_matrix)
}


MNL_create_moment_matrix(q=3,order=3)


#information matrix
MNL_compute_information_matrix <- function(design, beta, order) {
  q <- dim(design)[1]
  S <- dim(design)[3]

  n_parameters <- compute_parameters(q, order)

  model_array <- MNL_compute_model_matrix(design, order)
  P <- MNL_compute_choice_probabilites(design, beta, order)

  # Initialize the information matrix with zeros
  information_matrix <- matrix(0, nrow = n_parameters, ncol = n_parameters)

  for (s in 1:S) {
    p_s <- P[, s]
    I_s <- model_array[1:n_parameters, , s] %*% (diag(p_s) - tcrossprod(p_s)) %*% t(model_array[1:n_parameters, , s])
    information_matrix <- information_matrix + I_s
  }

  return(information_matrix)
}

MNL_compute_information_matrix(design,beta,order=3)





# I optimlaity
MNL_compute_i_optimality = function(design, order, beta) {
  q = dim(design)[1]
  printed_information_matrix = MNL_compute_information_matrix(design, beta, order)
  printed_moments_matrix = MNL_create_moment_matrix(q, order)
  # MNL_i_opt = inv(printed_information_matrix) %*% printed_moments_matrix
  MNL_i_opt = sum(diag(solve(printed_information_matrix, printed_moments_matrix)))  # trace is the sum of the diagonal elements
  return(MNL_i_opt)
}


MNL_compute_i_optimality(design,order=3,beta)



#Bayesian I optimality
beta_ma <- read.csv("beta_ma_keppa_05.csv")

list=list()
MNL_compute_bayesian_I_optiality= function(design,order,beta_ma){
  for(i in 1:128){


    beta_1=unlist(beta_ma[i,])
    I_optimality=MNL_compute_i_optimality(design, beta_1,order=3)

    list=append(list,I_optimality)

  }
  result=mean(unlist(list))
  return(result)
}
MNL_compute_bayesian_I_optiality(design,order=3,beta_ma)



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
    theme_nomask()+
    labs(x = "Ingredient 1", y = "Ingredient 2", z = "Ingredient 3")

  return(plot)
}








#keppa0_5_optimal_design
z0_5=c(0.00, 0.00, 1.00, 0.39, 0.34, 0.27, 0.00, 1.00, 0.00, 0.02, 0.00, 0.98, 1.00, 0.00, 0.00, 0.10, 0.70, 0.21, 0.00, 0.00, 1.00, 0.52, 0.00, 0.48, 0.00, 0.00, 1.00, 0.39, 0.34, 0.27, 0.52, 0.48, 0.00, 0.17, 0.15, 0.68, 0.06, 0.00, 0.94, 0.00, 0.51, 0.49)
kappa0_5_Iopt_reformat=array(z0_5,dim=c(3,2,7))
plot_ternary_design(kappa0_5_Iopt_reformat)
#keppa5_optimal_design
z5 =c(0.00, 0.00, 1.00, 0.47, 0.00, 0.53, 0.27, 0.43, 0.31, 0.08, 0.00, 0.92, 1.00, 0.00, 0.00, 0.39, 0.36, 0.25, 0.42, 0.27, 0.31, 0.00, 0.09, 0.91, 0.48, 0.52, 0.00, 0.18, 0.19, 0.62, 0.00, 1.00, 0.00, 0.33, 0.42, 0.24, 0.00, 0.00, 1.00, 0.00, 0.51, 0.49)
kappa5_Iopt_reformat=array(z5,dim=c(3,2,7))
plot_ternary_design(kappa5_Iopt_reformat)
#keppa10_optimal_design
z10=c(0.00, 0.00, 1.00, 0.39, 0.01, 0.60, 0.51, 0.25, 0.24, 0.98, 0.00, 0.02, 0.00, 0.10, 0.90, 0.38, 0.26, 0.36, 0.20, 0.28, 0.53, 0.44, 0.55, 0.01, 0.27, 0.41, 0.32, 0.12, 0.00, 0.88, 0.00, 0.00, 1.00, 0.01, 0.41, 0.58, 0.00, 0.96, 0.04, 0.27, 0.52, 0.21)
kappa10_Iopt_reformat=array(z10,dim=c(3,2,7))
plot_ternary_design(kappa10_Iopt_reformat)
#keppa30_optimal_design
z30=c(0.10, 0.59, 0.31, 0.00, 0.32, 0.68, 0.38, 0.00, 0.62, 0.66, 0.08, 0.26, 0.24, 0.54, 0.22, 0.00, 0.81, 0.19, 0.80, 0.00, 0.20, 0.56, 0.25, 0.19, 0.06, 0.26, 0.68, 0.00, 0.00, 1.00, 0.00, 0.00, 1.00, 0.26, 0.07, 0.67, 0.44, 0.50, 0.06, 0.29, 0.36, 0.35)
kappa30_Iopt_reformat=array(z30,dim=c(3,2,7))
plot_ternary_design(kappa30_Iopt_reformat)


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

clustered_design=clustering(design,k=12)
clustered_design_reformat=array(clustered_design,dim=c(3,2,6))
plot_ternary_design(design)
plot_ternary_design(clustered_design_reformat)

library(opdesmixr)
beta=c(0, 0.86      ,  -0.15762783,   2.46624476,   1.56402849,
       2.28214613, -21.61502097)
#elbow curve
elbow_plot <- function(random_design, beta, k_values = 7:13) {
  # Placeholder for the i-optimality values for each k.
  i_opt_values <- numeric(length(k_values))

  # For each k, perform the clustering and compute the I-optimality.
  for (i in seq_along(k_values)) {
    k <- k_values[i]

    # Perform the clustering.
    clustered_design <- clustering(random_design, k)

    # Reformat the clustered_design
    clustered_design_reformat <- array(clustered_design, dim = c(3, 2, k %/% 2))

    # Compute the I-optimality.
    i_opt_values[i] <- mnl_get_opt_crit_value(clustered_design_reformat, beta, 3, "I")
  }

  # Combine the k_values and i_opt_values into a data frame for plotting.
  elbow_data <- data.frame(k = k_values, i_optimality = i_opt_values)

  # Use ggplot2 to create the elbow plot.
  library(ggplot2)
  plot <- ggplot(elbow_data, aes(x = k, y = i_optimality)) +
    geom_line() +
    geom_point() +
    xlab("Number of clusters (k)") +
    ylab("I-Optimality") +
    ggtitle("Elbow Curve for Clustering")

  return(plot)
}

elbow_plot(design,beta)
#######################################
############V N S###############
#############################
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

#Final I optimality value
i_opt_value


#Save workspace
save.image("my_workspace11.RData")
load("my_workspace11.RData")

