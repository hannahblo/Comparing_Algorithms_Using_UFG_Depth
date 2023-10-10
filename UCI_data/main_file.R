# This file discusses the UCI Data Set.
# For further infomration see Appendix A.2 of Christoph Jansen, Malte Nalenz,
# Georg Schollmeyer and Thomas Augustin (2023): Statistical comparisons of
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

### All the other R-packages are on CRAN packages (06.10.2023)

library(reshape)
library(dplyr)
library(data.table)
library(ggplot2)
library(forcats)
library(gridExtra)
library(stargazer)
library(prefmod)
library(ddandrda)
library(gurobi)
library(utils)
library(hasseDiagram)

setwd("UCI_data/")


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

# Function which converts a list of posets as matrices (see convert_to_matrix
# function) to format needed for prefmod package (analogouse to prefmod::cemspc)
convert_to_rowbt <- function(list_graph) {

  number_cols <- nrow(list_graph[[1]]) * (nrow(list_graph[[1]]) - 1) / 2
  model_comp <- rownames(list_graph[[1]])
  col_names <- list()
  for (i in 1:(length(model_comp) - 1)) {
    for (j in (i + 1):length(model_comp)) {
      col_names <- append(col_names, paste0(model_comp[i], ">", model_comp[j]))
    }
  }


  mat_bt <- matrix(ncol = number_cols, nrow = length(list_graph))
  colnames(mat_bt) <- col_names

  for (graph_index in 1:length(list_graph)) {
    col_index <- 1
    for (i in 1:(length(model_comp) - 1)) {
      for (j in (i + 1):length(model_comp)) {
        if (list_graph[[graph_index]][i,j] == 1) {
          mat_bt[graph_index, col_index] <- 2
        } else if (list_graph[[graph_index]][j,i] == 1) {
          mat_bt[graph_index, col_index] <- 0
        } else {
          mat_bt[graph_index, col_index] <- 1
        }
        col_index <- col_index + 1
      }
    }
  }

  return(mat_bt)
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
output_brier    = dcast(table_res_brier,dataset~method)
output_brier = cbind(rep("brier_score", 16), output_brier)
colnames(output_brier)[1] = "metric"
# brier score umdrehen dest kleiner desto besser


table_res_auc = full_res  %>% group_by(metric, dataset, method) %>%
  filter(metric == "auc") %>% summarise(mu = mean(val)) %>% ungroup()
output_auc    = dcast(table_res_auc,dataset~method)
output_auc = cbind(rep("auc", 16), output_auc)
colnames(output_auc)[1] = "metric"

table_res_acc = full_res  %>% group_by(metric, dataset, method) %>%
  filter(metric == "acc") %>% summarise(mu = mean(val)) %>% ungroup()
output_acc    = dcast(table_res_acc,dataset~method)
output_acc = cbind(rep("acc", 16), output_acc)
colnames(output_acc)[1] = "metric"

output = rbind(output_brier, output_acc, output_auc)
dim(output)


# Brier Score needs to be switched
index_brier <- which(output$metric == "brier_score")
output[, seq(3,10)] <- 1 - output[, seq(3, 10)]

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
# we convert as follows, the entry (i [row],j [column]) of the matrix denotes that
# i has a lower performance than j


mat_bt <- convert_to_rowbt(list_graph)

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

duplicated(list_graph)
index_duplicated <- which(duplicated(list_graph) == TRUE)

pdf("duplicates_in_observed.pdf", onefile = TRUE)
for (i in 1:length(index_duplicated)) {
  mat <- matrix(as.logical(list_graph[[i]]), ncol = 8)
  colnames(mat) <- rownames(mat) <- colnames(list_graph[[i]])
  hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()

pdf("all_observed.pdf", onefile = TRUE)
for (i in 1:length(list_graph)) {
  mat <- matrix(as.logical(list_graph[[i]]), ncol = 8)
  colnames(mat) <- rownames(mat) <- colnames(list_graph[[i]])
  hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()


# Heatmap
edges <- Reduce("+", list_graph)
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
# for further information, see the to the Repesotory corresponding article or
# https://isipta23.sipta.org/wp-content/uploads/2023/06/blocher23a.pdf
################################################################################
item_number <- 8
names_columns <- colnames(list_graph[[1]])

# Since we consider here only 16 partial orders going through all subsets of
# those 16 elements is faster than computing all posets with 8(!) items


depth_premises <- ddandrda::compute_ufg_depth_porder(list_graph,
                                            print_progress_text = FALSE,
                                            save_ufg_premises = FALSE)


## All partial orders with depth larger than 0
max_depth_index <- sort(depth_premises$depth_ufg, index.return = TRUE, decreasing = TRUE)$ix
pdf("plots_observed_from_highest_depth.pdf", onefile = TRUE)
for (i in max_depth_index) {
  mat <- matrix(as.logical(list_graph[[i]]), ncol = item_number)
  colnames(mat) <- rownames(mat) <- names_columns
  hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()



## Intersections (high to low)
max_depth_index <- sort(depth_premises$depth_ufg, index.return = TRUE, decreasing = TRUE)$ix
pdf("plots_observed_intersect_from_highes.pdf", onefile = TRUE)
for (i in 1:min(max_depth_index)) {
  intersect <- matrix(rep(TRUE, item_number * item_number), ncol = item_number)
  for (j in seq(1, i)) {
    intersect <- intersect & matrix(as.logical(list_graph[[max_depth_index[j]]]), ncol = item_number)
  }
  colnames(intersect) <- rownames(intersect) <- names_columns
  # print(mat * 1)
  # hasse(mat) plots the graph from top to bottom, with smallest value at bottom
  # -> change and set arrow to "backward"
  hasse(t(intersect), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()




# Nevertheless, we want to obtain the poset with highest depth value based on
# ALL (not only observed posets). The number of all possible posets is
# 431723379 (see: https://oeis.org/A001035)
# We use that we immediatly get that some posets need to be zero. More precise
# all posets with a depth larger than zero need to lie in the union of the edges
# of all observed posets. Since not all edges are observed
# (e.g. BoostedStumps < CART) we can reduced the posets to consider

union_edges <- Reduce("|", list_graph) * 1
edge_exists_index <- setdiff(which(union_edges == 1), c(1, 10, 19, 28, 37, 47, 55, 64))


length(edge_exists_index)
poset_list <- list()
for (cardinality in seq(1, length(edge_exists_index))) { # only needs a few seconds
  print(paste0("We are at cardinality ", cardinality))
  combinations <- utils::combn(x = edge_exists_index, m = cardinality,
                               simplify = FALSE)
  for (subset_index in combinations) {
    poset <- diag(x = 1, nrow = 8)
    poset[subset_index] <- 1
    if (ddandrda::test_if_porder(poset)) {
      poset_list <- append(poset_list, list(poset))
    }
  }
}


length(poset_list)


depth_poset_list <- ddandrda::compute_ufg_depth_porder(porder_observed = list_graph,
                                                       porder_depth = poset_list,
                                                       print_progress_text = FALSE) # needs a few seconds



## Some pictures and vlaues
print(paste0("The minimal value in list is ", min(depth_poset_list$depth_ufg), ". Note that the smalles value of ALL posets is 0."))
print(paste0("The maximal value is ", max(depth_poset_list$depth_ufg)))
#### mean, sd, distribution etc value are not meaningful as not all posets are considered!





## All partial orders with depth larger than 0
max_depth_index <- sort(depth_poset_list$depth_ufg, index.return = TRUE, decreasing = TRUE)$ix
stop_print <- 100
pdf("plots_all_from_highest_depth.pdf", onefile = TRUE)
for (i in max_depth_index[seq(1,stop_print)]) {
  mat <- matrix(as.logical(poset_list[[i]]), ncol = item_number)
  colnames(mat) <- rownames(mat) <- names_columns
  hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()



## Intersections (high to low)
max_depth_index <- sort(depth_poset_list$depth_ufg, index.return = TRUE, decreasing = TRUE)$ix
pdf("plots_all_intersect_from_highes.pdf", onefile = TRUE)
stop_print <- 100
for (i in 1:min(length(max_depth_index), stop_print)) {
  intersect <- matrix(rep(TRUE, item_number * item_number), ncol = item_number)
  for (j in seq(1, i)) {
    intersect <- intersect & matrix(as.logical(poset_list[[max_depth_index[j]]]), ncol = item_number)
  }
  colnames(intersect) <- rownames(intersect) <- names_columns
  # print(mat * 1)
  # hasse(mat) plots the graph from top to bottom, with smallest value at bottom
  # -> change and set arrow to "backward"
  hasse(t(intersect), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()


# minimal to maximal depth value does not make sense at this needs to contain all
# zeros which lead directly to an empty poset.





################################################################################
# Computation of the  Bradley-Terry Model with no convariables, but with ties
# for further informations see:
# - Hatzinger, Dittrich (2012): prefmod: An R Package for Modeling Preferences
#    Based on Paired Comparisons, Rankings, or Ratings
# - Hatzinger, Maier (2023): Package ‘prefmod’
################################################################################

# Constructing the design matrix
design_mat <- construct_design_bt(list_graph)

colnames(design_mat)[4] <- "Boosted_Stumps"
design_mat$mu <- as.factor(design_mat$mu)

result_UCI <- gnm(cum_sum ~ undecided + Boosted_Stumps + CART + ElasticNet + GBM +
                    glm + Lasso + rf + Ridge,
                  elim = mu, # see: https://cran.r-project.org/web/packages/gnm/gnm.pdf (page 33)
                  data = design_mat,
                  family = poisson)
summary(result_UCI)

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



