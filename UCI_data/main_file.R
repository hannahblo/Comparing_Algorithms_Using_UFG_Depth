# This file discusses the UCI Data Set.

# For further information on the data see Appendix A.2 of Christoph Jansen, Malte
# Nalenz, Georg Schollmeyer and Thomas Augustin (2023): Statistical comparisons of
# classifiers by generalized stochastic dominance.  Journal of Machine Learning
# Research,  24: 1 - 37.


################################################################################
# R Session
################################################################################

### Information about packages ddandrda. This package is under developement on
# git. Installation can be done by:
# remove.packages("ddandrda")
# install.packages("devtools")
# devtools::install_github("hannahblo/ddandrda")
library(ddandrda)

### All the other R-packages are on CRAN packages (06.10.2023)
library(reshape)
library(dplyr)
library(data.table)
library(ggplot2)
library(forcats)
library(gridExtra)
library(stargazer)
library(prefmod)
library(reshape2)
library(utils)
library(hasseDiagram)
# For R versions >=3.5 do the following before installing hasseDiagramm:
# install.packages("BiocManager")
# BiocManager::install("Rgraphviz")
library(gurobi)
# This is a commercial solver that offers a free academic licenses which can be
# found here: https://www.gurobi.com/features/academic-named-user-license/ (accessed: 08.02.2023).
# To install this package, please follow the instructions there
# A documentation can be found here: https://www.gurobi.com/wp-content/plugins/hd_documentations/documentation/9.0/refman.pdf (page 643ff) (accessed: 08.02.2023).



setwd("UCI_data/")

# TODO
# warum testen ob wirklich poset ist


################################################################################
# Functions needed later
################################################################################
# Function which converts to to partial orders
convert_to_matrix <- function(single_data_eval) {
  list_learner <- colnames(single_data_eval)
  number_learner <- length(list_learner)
  graph_mat <- matrix(rep(0, number_learner * number_learner),
                      nrow = number_learner)
  rownames(graph_mat) <- colnames(graph_mat) <- list_learner
  diag(graph_mat) <- 1

  for (learner_i in list_learner) {
    learner_base <- single_data_eval[[learner_i]]
    for (learner_j in list_learner) {
      learner_comp <- single_data_eval[[learner_j]]
      if (all(learner_base < learner_comp)) {
        graph_mat[learner_i, learner_j] <- 1
      }
    }
  }
  return(graph_mat)
}




# this functions returns the gurobi model it is the basis of an MILP for
# computing the poset (or a poset) with the highest depth
define_gurobi_ufg <- function(depth_premises, constant_weight = 16){

  # Here, we define a gurobi model to compute the minimal and maximal depth values
  # More Information about the Gurobi Model can be found here:
  # https://www.gurobi.com/wp-content/plugins/hd_documentations/documentation/9.0/refman.pdf (page 643ff) (accessed: 08.02.2023).

  # general information
  n_items <- nrow((depth_premises$ufg_premises[[1]])[[1]]) # number of items

  ## Step 1: Defining the constraint matrix. The constraint matrix must capture
  # the following structure (included into the rows)
  # 1. if p is an element of an ufg set in S, than it must be a superset of the
  #     intersection
  # 2. if p is an element of an ufg set in S, then the union must be subset of
  #     the union od all elements in S
  # 3. if the intersection of an ufg set is a subset of p and if the union of the ufg set is a superset of p, then this ufg set MUST be counted for the depth
  # 4. p must be a partial order --> antisymmetric
  # 5. p must be a partial order --> transitive
  # 6. p must be a partial order --> reflexive


  # To achieve this, the columns of the constraint matrix are given by
  # col_1. The first seq(1, number_of_premises) columns represent each (set) element of S
  # col_2. the remaining rows represent the edges (or n_items^2 pairs of p)
  # In part col_1 an element is 1 if and only if  we are currently looking at
  # this element of Sscr



  # Constraints 1:
  # matrix for the constraints for the implications of the intersections:
  # If (a,b) is in the intersection of an ufg set of S, then for this premise
  # to be counted for a poset p, also p has to have the pair (a,b)
  A_intersection <- array(0,c(depth_premises$total_number_premises, depth_premises$total_number_premises + n_items^2))

  # Constraints 2:
  # matrix for the constraints for the implications of the (non-) union:
  # If a pair (a,b) is not in the union ('in the non-union") of the ufg premise,
  # then for this premise to be countet for a poset p, also p has to have NOT the pair (a,b)
  A_non_union <- array(0,c(depth_premises$total_number_premises, depth_premises$total_number_premises + n_items^2))

  # Contraints 3:
  # matrix for the constraints that model that if the intersection of an ufg set is a subset of p and p is a subset of the union of the ufg set, then the premise MUST be counted
  A_premise <- array(0,c(depth_premises$total_number_premises, depth_premises$total_number_premises + n_items^2))
  rhs_premise <- rep(0,depth_premises$total_number_premises)


   # objective of the MIL program: simply adding the set Sscr (weighted by the corresponding weight) which imply a poset p
  obj <- rep(0, depth_premises$total_number + n_items^2)

  # Now, we go through every set in $Sscr$ see Definition 2 of the corresponding paper.
  for (k in seq_len(depth_premises$total_number_premises)) {

    # If the partial order lies in the conclusion of the ufg premises, then we add
    # the following amount onto the ufg depth value.
    obj[k] <- (1/constant_weight)^length(depth_premises$ufg_premises[[k]]) / depth_premises$constant_cn

    # computing the intersection of all elements of the ufg premises
    intersection <- as.vector(Reduce("pmin", depth_premises$ufg_premises[[k]]))

    # computing the edges which never occur
    non_union <- as.vector(1 - Reduce("pmax", depth_premises$ufg_premises[[k]]))

    # check how large intersection and non existence of edges are
    # If both are empty, then every poset needs to lie in the conclusion. Else we
    # need a restriction from one of the two sides. Compare with Definition 1 in
    # the corresponding article.
    n_intersection <- sum(intersection)
    n_non_union <- sum(non_union)

    if (n_intersection >= 1) {
      # the following constraints model that if (a,b) in the intersection of a premise that is counted then also p has to have the pair (a,b).
      # premise counted => all pairs of the intersection in p: x_pairsofintersection/#intersection >= x_premiseiscounted
      A_intersection[k, k] <- -1 # this ufg premise S has an constraint which is given by the next line:
      A_intersection[k, depth_premises$total_number_premises +
                       which(intersection == 1)] <- 1/n_intersection
      # note that rhs part of the gurobi model will be set to 0. Thus, this constraint
      # implies that if this ufg premise wants to be used in the sum of the
      # depth, then the second part of the columns must be fulfilled.
    }

    if (n_non_union >= 1) {
      # the following constraints model that if (a,b) is in the non-union of a premise that is counted then also p has to have NOT the pair (a,b).
      # premise counted => all pairs of the non-union NOT in p: x_pairsofnonunion/#nonunion <= 1 - x_premiseiscounted
      A_non_union[k, k] <- -1
      A_non_union[k, depth_premises$total_number_premises +
                    which(non_union == 1)] <- -1/n_non_union
      # note that rhs part of the gurobi model will be set to -1. Thus, this constraint
      # implies that if this ufg premise wants to be used in the sum of the
      # depth, then the second part of the columns must be fulfilled.
    }

	if (n_non_union >= 1 | n_intersection >= 1) {
      # the following constraints model that if (a,b) is in the non-union of a premise that is counted then also p has to have NOT the pair (a,b).
      # premise counted => all pairs of the non-union NOT in p: x_pairsofnonunion/#nonunion <= 1 - x_premiseiscounted
      A_premise[k, k] <- 1
      A_premise[k, depth_premises$total_number_premises +
                    which(non_union == 1)] <- 1
	    A_premise[k, depth_premises$total_number_premises +
                    which(intersection == 1)] <- -1
	    rhs_premise[k] <- 1 - n_intersection
      # note that rhs part of the gurobi model will be set to -1. Thus, this constraint
      # implies that if this ufg premise wants to be used in the sum of the
      # depth, then the second part of the columns must be fulfilled.
    }

  }

  # Recall the comments in the beginning of this function. Part 5, the transitivity,
  # is now discussed
  # We go now through all pairs, and if (a,b) and (b,c) in the poset, then (a,c)
  # needs to be as well
  # Note that only the second part of the columns of the constraints are now of
  # interest. Since this has notehing to do with the ufg premises and conlcusions
  A_transitive <- array(0, c(n_items*(n_items - 1)*(n_items - 2), ncol(A_intersection)))
  # index[i,j] denotes column i and row j
  indexs <- seq_len(n_items^2)
  dim(indexs) <- c(n_items, n_items)
  t <- 1
  # transitivity means that if (i,j) and (j,k) in a poset, then also (i,k) must
  # be an element of the poset. Thus, x_ik >= x_ij +x_jk -1
  for (i in seq_len(n_items)) {
    for (j in seq_len(n_items)[-i]) {
      for (k in seq_len(n_items)[-c(i,j)]) {
        A_transitive[t, depth_premises$total_number_premises +
                       c(indexs[i,j],indexs[j,k])] <- -1
        A_transitive[t, depth_premises$total_number_premises +
                       indexs[i,k]] <- 1
        t <- t + 1

        # note that rhs part of the gurobi model will be set to -1. Thus, if [i,j]
        # and [j,k] in the poset, then [i,k] needs to be as well
      }
    }
  }


  # Recall the comments in the beginning of this function. Part 4, the antisymmetrie,
  # is now discussed
  # Note that only the second part of the columns of the constraints are now of
  # interest. Since this has notehing to do with the ufg premises and conlcusions
  A_antisym <- array(0,c(n_items*(n_items - 1),ncol(A_transitive)))
  t <- 1
  # we now model that if (i,j) in p, then (j,i) notin p
  # Thus, we ensure that x_ij +x_ji <=1
  for (i in seq_len(n_items - 1)) {
    for (j in seq((i + 1), n_items)) {

      A_antisym[t, depth_premises$total_number_premises +
                  c(indexs[i,j],indexs[j,i]) ] <- -1
      t <- t + 1

      # note that rhs part of the gurobi model will be set to -1. Thus, if [i,j]
      # and [j,i] cannot be in the poset
    }
  }


  # gurobi_model (without reflexivity part, see below)
  model <- list(modelsense = "max",
                obj = obj,
                lb = rep(0,depth_premises$total_number_premises + n_items^2),
                ub = rep(1,depth_premises$total_number_premises + n_items^2),
                vtype = rep("B",depth_premises$total_number_premises + n_items^2),
                A = rbind(A_intersection, A_non_union, A_transitive, A_antisym,A_premise),
                rhs = c(rep(0,depth_premises$total_number_premises),
                        rep(-1,depth_premises$total_number_premises),
                        rep(-1,nrow(A_transitive)),
                        rep(-1,nrow(A_antisym)),rhs_premise),
                sense = rep(">=", 3*depth_premises$total_number_premises +
                              nrow(A_transitive) + nrow(A_antisym)))

  # Part 6: reflexivity
  # Note that we need to ensure that the reflexive part is also given, thus
  # there we set the lowerbound to 1
  for (k in (1:n_items)) {
    model$lb[depth_premises$total_number_premises + indexs[k,k]] <- 1
  }
  return(model)
}


# Function which converts a list of posets as matrices (see convert_to_matrix
# function) to format needed for prefmod package (analogous to prefmod::design(prefmod::cemspc))
construct_design_bt <- function(list_graph) {

  heatmap <-  Reduce("+", list_graph)
  cl_names <- colnames(list_graph[[1]])
  df_design_bt <- data.frame(matrix(ncol = 2 + 3 + dim(heatmap)[1], nrow = 0))
  colnames(df_design_bt) <- c("cum_sum", "mu", "pref_a", "undecided", "pref_b", colnames(list_graph[[1]]))
  sum_graphs <- length(list_graph)


  mu <- 1
  for (i in seq_len(length(cl_names) - 1)) {
    for (j in seq(i + 1, length(cl_names))) {
      input <- rep(0, length(cl_names))

      input[c(i,j)] <- c(-1, 1)
      df_design_bt[nrow(df_design_bt) + 1, ] <- c(heatmap[i,j], mu, 1, 0, 0, input)

      input[c(i,j)] <- c(0, 0)
      df_design_bt[nrow(df_design_bt) + 1, ] <- c(sum_graphs - heatmap[i,j] - heatmap[j,i], mu, 0, 1, 0, input)

      input[c(i,j)] <- c(1, -1)
      df_design_bt[nrow(df_design_bt) + 1, ] <- c(heatmap[j,i], mu, 0, 0, 1, input)

      mu <- mu + 1
    }
  }


  return(df_design_bt)
}
################################################################################
# Upload Algorithm Performances
# The following code part as well as the performance evaluation was already done
# by Christoph Jansen, Malte Nalenz, Georg Schollmeyer and homas Augustin (2023):
# Statistical comparisons of classifiers by generalized tochastic dominance.
# Journal of Machine Learning Research,  24: 1 - 37.
################################################################################

res_files = list.files("performance_results/")
res_list = list()
for (i in 1:length(res_files)) {
  res_list[[i]] = read.csv(paste0("performance_results/", res_files[i]))
}


full_res  = do.call('rbind', res_list)

table_res_brier = full_res  %>% group_by(metric, dataset, method) %>%
  filter(metric == "brier_score") %>% summarise(mu = mean(val)) %>% ungroup()
output_brier  = reshape2::dcast(table_res_brier, dataset~method)
output_brier = cbind(rep("brier_score", 16), output_brier)
colnames(output_brier)[1] = "metric"


table_res_auc = full_res  %>% group_by(metric, dataset, method) %>%
  filter(metric == "auc") %>% summarise(mu = mean(val)) %>% ungroup()
output_auc    = reshape2::dcast(table_res_auc,dataset~method)
output_auc = cbind(rep("auc", 16), output_auc)
colnames(output_auc)[1] = "metric"

table_res_acc = full_res  %>% group_by(metric, dataset, method) %>%
  filter(metric == "acc") %>% summarise(mu = mean(val)) %>% ungroup()
output_acc    = reshape2::dcast(table_res_acc,dataset~method)
output_acc = cbind(rep("acc", 16), output_acc)
colnames(output_acc)[1] = "metric"

output = rbind(output_brier, output_acc, output_auc)
dim(output)


# Brier Score needs to be switched
index_brier <- which(output$metric == "brier_score")
output[index_brier, seq(3,10)] <- 1 - output[index_brier, seq(3, 10)]

# View(output)





################################################################################
# Construction of the Partial Orders
################################################################################
data_sets <- unique(full_res$dataset)
list_graph <- list()

for (data_name in data_sets) {
  single_data_eval <- filter(output, dataset == data_name)
  single_data_eval <- single_data_eval[, -c(1,2)]
  list_graph <- append(list_graph, list(convert_to_matrix(single_data_eval)))
}

################################################################################
#
# PART 1: FIRST IMPRESSION
#
################################################################################
setwd("depth_results/")


### Which edge exists
length(list_graph) # 80
length(unique(list_graph)) # 58
Reduce("|", list_graph)
Reduce("&", list_graph)
Reduce("+", list_graph)

duplicated(list_graph) # no duplications exist


pdf("all_observed.pdf", onefile = TRUE)
for (i in 1:length(list_graph)) {
  mat <- matrix(as.logical(list_graph[[i]]), ncol = 8)
  colnames(mat) <- rownames(mat) <- colnames(list_graph[[i]])
  hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()


# Heatmap
edges <- Reduce("+", list_graph)
colnames(edges) <- rownames(edges) <- c("BS", "CART", "EN", "GBM", "GLM", "LASSO", "RF", "RIDGE")
df_edge_exist <- melt(edges)
df_edge_exist <- df_edge_exist[df_edge_exist$value != 0, ]


jpeg(file = "heatmap_UCI.jpeg")
ggplot(df_edge_exist, aes(x = Var1, y = Var2)) +
  geom_raster(aes(fill = value)) +
  scale_fill_gradient(low = "lightcyan1", high = "darkcyan") +
  labs(x = "is below", y = "is above") +
  geom_text(aes(label = value)) +
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.3),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 13, vjust = -1),
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        plot.title = element_blank())
dev.off()
# note that the plotting is in some sense counterintutitve since it switches the
# x and y axis! Compare to edges




################################################################################
#
# PART 2: COMPUTATION OF UFG DEPHT + EXTENDED BRADLEY TERRY MODEL
#
################################################################################

################################################################################
# Computation of the UFG Depth value
# for further information, see the to the Repository corresponding article or
# https://isipta23.sipta.org/wp-content/uploads/2023/06/blocher23a.pdf
################################################################################
item_number <- 8
names_columns <- colnames(list_graph[[1]])

# Since we consider here only 16 partial orders going through all subsets of
# those 16 elements is faster than computing all posets with 8(!) items

start_time <- Sys.time()
depth_premises <- ddandrda::compute_ufg_depth_porder(list_graph,
                                            print_progress_text = TRUE,
                                            save_ufg_premises = TRUE)
total_time <- Sys.time() - start_time
total_time # Time difference of 1.978799 mins

# saveRDS(depth_premises, "depth_premises.rds")
# saveRDS(total_time, "total_time_depth_premises.rds")

max(depth_premises$depth_ufg) # [1] 0.31869

# observed deepest depth values
max_depth_index <- sort(depth_premises$depth_ufg, index.return = TRUE, decreasing = TRUE)$ix
pdf("plots_observed_from_highest_depth.pdf", onefile = TRUE)
for (i in max_depth_index) {
  current_poset <- matrix(as.logical(list_graph[[i]]), ncol = item_number)
  colnames(current_poset) <- rownames(current_poset) <- names_columns
  hasse(t(current_poset), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()



# Intersections (high to low) of the observed depth values
max_depth_index <- sort(depth_premises$depth_ufg, index.return = TRUE, decreasing = TRUE)$ix
pdf("plots_observed_intersect_from_highes.pdf", onefile = TRUE)
for (i in 1:max(max_depth_index)) {
  intersect <- matrix(rep(TRUE, item_number * item_number), ncol = item_number)
  for (j in seq(1, i)) {
    intersect <- intersect & matrix(as.logical(list_graph[[max_depth_index[j]]]), ncol = item_number)
  }
  colnames(intersect) <- rownames(intersect) <- names_columns
  hasse(t(intersect), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()



# Nevertheless, we want to obtain the poset with highest depth value based on
# ALL (not only observed posets). The number of all possible posets is
# 431723379 (see: https://oeis.org/A001035
# Therefore, we implemented an mixed integer linear programm and use gurobi.

# defining and optimize the model
start_time <- Sys.time()
gurobi_model <- define_gurobi_ufg(depth_premises)
total_time <- Sys.time() - start_time
total_time

N <- 20
start_time <- Sys.time()
highest_depth <- gurobi::gurobi(gurobi_model, params = list(OutputFlag = 1,
                                                            PoolSearchMode = 2,
                                                            PoolSolutions = N))
total_time <- Sys.time() - start_time
# saveRDS(total_time, paste0("total_time_", N, "_deepest_posets.rds"))


pdf(paste0("plots_all_", N, "_from_highest_depth.pdf"), onefile = TRUE)
for (k in seq(1,N)) {

  current_poset <- round(highest_depth$pool[[k]]$xn[-seq_len(depth_premises$total_number_premises)], 3) # TODOOO
  dim(current_poset) <- c(item_number, item_number)
  current_poset <- matrix(as.logical(current_poset), ncol = item_number)
  colnames(current_poset) <- rownames(current_poset) <- names_columns

  # TODO
  # FEHLT KOMMENTAR WARUM MAN POSET TESTEN MUSS
  if (ddandrda::test_if_porder(current_poset)) {
    print(paste0("Depth of ", k, " deepest poset is ", highest_depth$pool[[k]]$objval))
    hasse(t(current_poset), parameters = list(arrow = "backward", shape = "roundrect"))
  } else {
    print(paste0("Attention: ", k, " is not a poset!"))
  }
}
dev.off()






################################################################################
# Computation of the  Bradley-Terry Model with no convariables, but with ties
# for further informations see:
# - Hatzinger, Dittrich (2012): prefmod: An R Package for Modeling Preferences
#    Based on Paired Comparisons, Rankings, or Ratings
# - Hatzinger, Maier (2023): Package ‘prefmod’
# - Sinclaig (1982): GLIM for Preferences
# - Davidson (1970): On Extending the Bradley-Terry Model to Accommodate Ties in
#    Paired Comparison Experiments
################################################################################

# Constructing the design matrix
design_mat <- construct_design_bt(list_graph)

colnames(design_mat)[6] <- "Boosted_Stumps"
design_mat$mu <- as.factor(design_mat$mu)

result_UCI <- gnm(cum_sum ~ undecided + Boosted_Stumps + CART + ElasticNet + GBM +
                    glm + Lasso + rf + Ridge,
                  elim = mu, # see: https://cran.r-project.org/web/packages/gnm/gnm.pdf (page 33)
                  data = design_mat,
                  family = poisson)
summary(result_UCI)

## another approach to obtain the model is by the standard glm R function.
## Note that here mu is printed as well
# result_Uci_glm <- glm(cum_sum ~ -1 + mu + undecided + Boosted_Stumps + CART + ElasticNet + GBM +
#                         glm + Lasso + rf + Ridge,
#                       family = 'poisson', data = design_mat) # loglink is default
# summary(result_Uci_glm)


# lambdas
lambda <- result_UCI$coefficients[seq(2,9)]
lambda$Ridge <- 0

# worth
division <- sum(unlist(lapply(lambda, function(x){exp(2*x)})))
worth_UCI <- lapply(lambda, function(x){exp(2*x)/division})
sort_worth_UCI <- sort(unlist(worth_UCI))

graph_bt <- diag(TRUE, length(worth_UCI))
colnames(graph_bt) <- colnames(list_graph[[1]])
rownames(graph_bt) <- rownames(list_graph[[1]])
colnames(graph_bt)[1] <- "Boosted_Stumps"
rownames(graph_bt)[1] <- "Boosted_Stumps"

for (index in 1:(length(sort_worth_UCI) - 1)) {
  row <- names(sort_worth_UCI[index])
  higher <- which(sort_worth_UCI > sort_worth_UCI[row])
  columns <- names(higher)
  graph_bt[row, columns] <- TRUE
}

pdf("plot_bradley_terry.pdf", onefile = TRUE)
hasse(t(graph_bt), parameters = list(arrow = "backward", shape = "roundrect"))
dev.off()

################################################################################
# Classifier comparison by generalized stochastic dominance
# this was already computed, see Christoph Jansen, Malte Nalenz, Georg Schollmeyer
# and Thomas Augustin (2023): Statistical comparisons of classifiers by
# generalized stochastic dominance.  Journal of Machine Learning Research,24: 1 - 37.
################################################################################



