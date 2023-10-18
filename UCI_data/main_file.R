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
library(reshape2)
library(utils)
library(hasseDiagram) # For R versions >=3.5 do the following bevor installing hasseDiagramm:
# install.packages("BiocManager")
# BiocManager::install("Rgraphviz")



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


# this function is a help function for compute_upper_ufg_bound
sumup_prop_set <- function(set, weight_constant = (1/16)) {
  length_ufg <- lapply(X = set, FUN = function(x){length(x)})
  sum <- sum(unlist(lapply(X = length_ufg, FUN = function(x) {prod(rep(1/16, x))})))
  return(sum)
}



# fixed edges is a binary square matrix, needs to be rbind(c(row,col), c(row, col))
compute_upper_ufg_bound <- function(fixed_edges, fixed_none_edges, further_edge,
                                        depth_premises, constant_weight = (1/16)) {
  # Sorting ufg premises
  ufg_has_edge <- list()
  ufg_has_no_edge <- list()

  for (index_ufg in seq(1, depth_premises$total_number_premises)) {
    ufg <-  depth_premises$ufg_premises[[index_ufg]]


    if (all(Reduce("|", ufg)[rbind(fixed_edges, as.vector(further_edge))]) && (!(any(Reduce("&", ufg)[fixed_none_edges])))) {
      ufg_has_edge <- append(ufg_has_edge, list(ufg))
    }
    if (all(Reduce("|", ufg)[fixed_edges]) && (!(any(Reduce("&", ufg)[rbind(fixed_none_edges, as.vector(further_edge))])))) {
      ufg_has_no_edge <- append(ufg_has_no_edge, list(ufg))
    }

  }

  # Computing lower and upper bound
  upper_has_edge <- sumup_prop_set(ufg_has_edge) / depth_premises$constant_cn
  upper_has_not_edge <- sumup_prop_set(ufg_has_no_edge) / depth_premises$constant_cn

  return(list(fixed_edges = fixed_edges,
              fixed_none_edges = fixed_none_edges,
              further_edge = further_edge,
              upper_has_edge = upper_has_edge,
              upper_has_not_edge = upper_has_not_edge))

}



# Help Function which returns the cartesian product
help_fct_cartesian <- function(x) {
  result <- list()
  for (i in seq(1, 8)) {
    result <- append(result, list(c(x,i)))
  }
  return(result)
}

# This function computes the ufg of a poset where all ufg premises are already
# computed
compute_ufg_given_premises <- function(poset_interest, ufg, cn_constant,
                                       constant_weight = (1/16)) {

  depth <- 0
  for (premises in ufg) {
    if (all(ddandrda::test_porder_in_concl(premises, list(poset_interest)))) {
      depth <- depth + prod(rep(1/16, length(premises)))
    }
  }
  return(depth/cn_constant)
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
# brier score umdrehen dest kleiner desto besser


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
                                            print_progress_text = TRUE,
                                            save_ufg_premises = TRUE)

max(depth_premises$depth_ufg) # [1] 0.31869

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
for (i in 1:max(max_depth_index)) {
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
# We use that we immediately the Appendix of the corresponding article to get upper
# bounds for posets which contain or do not contain a certain edge. Since via the
# observed posets we already have a lower bound of the maximal ufg depth (max observed depth)
# we can use this to determine already some edges if they lie in the poset with
# maximal depth or not-

# IMPORTANT: Please read the following lines before proceeding.
# Line 366 to 440 gives the procedure to get all fixed edges (not edges). This code
# must be run twice.
# 1. In the first run, comment all code lines with comment: "comment out - Round 0"
#    and run the code from line 366 to 440 once
# 2. In the second step comment out all code lines with the comment "comment out - Round 1"
#    and run again through the entire code form line 366 to 440
#    Note that in Round 1 you also have to run the codes with te comment Round 0
# These two rounds give you the entire computation





###### Procedure: Check which edges are fixed by upper and lower bound discussions
## START
# Attention we overwrite the data frame
df_upper_bounds <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df_upper_bounds) <- c("fixed_edges", "fixed_none_edges", "further_edge", "upper_has_edge", "upper_has_not_edge")

all_edges_fixed <- list(c(1,1), c(2,2), c(3,3), c(4,4), c(5,5), c(6,6), c(7,7), c(8,8) #  Round 0: due to reflexivity these edges need to be set
                        # , c(2,1), c(2,4) # comment - Round 1
                        )
all_non_edges_fixed <- list(c(1,2), c(4,2), c(8,2) # Round 0: due to the heat map we can observe
                            # that these edge were never observed and therefore every poset which has
                            # this edge need to have depth zero
                            # , c(1,3), c(1,6), c(1,8), c(3,2), c(3,6), c(5,2), c(5,3),  # comment out - Round 1
                            # c(5,6), c(5,8), c(6,2), c(7,2)  # comment out - Round 1
                            )

row_col_interest <- unlist(lapply(seq_len(8), FUN = help_fct_cartesian), recursive = FALSE)
row_col_interest <- row_col_interest[-c(1, 10, 19, 28, 37, 46, 55, 64)] # Round 0
row_col_interest <- row_col_interest[-c(1, 23, 51)] # Round 0
row_col_interest <- row_col_interest[-c(7, 9)] # comment out - Round 1
row_col_interest <- row_col_interest[-c(1, 4, 6, 13, 16, 26, 27, 29, 31, 33, 40)] # comment out - Round 1
# length(row_col_interest)


# get lower bound of
# Round 0: we take the maximal observed ufg depth as starting point
max_depth_value <- max(depth_premises$depth_ufg)
# Round 1: check if the poset where we only add all edges fixed has higher depth then
# the one fixed at Round 0
# No, it is below the observed max --> stay at the same observed max as Round 0
# Note that indeed all_non_edges_fixed is true for observed max ufg depth
# check_poset <- diag(8)
# check_poset[rbind(c(2,1), c(2,4))] <- 1
# check_poset <- ddandrda::compute_transitive_hull(check_poset)
# check_poset # checken ob non edges fix passt
# depth_premises_check <- ddandrda::compute_ufg_depth_porder(porder_observed = list_graph,
#                                                            porder_depth = list(check_poset),
#                                                            print_progress_text = FALSE,
#                                                            save_ufg_premises = FALSE)
# depth_premises_check$depth_ufg



for (row_col_inner in row_col_interest) {
  bounds <- compute_upper_ufg_bound(fixed_edges = t(matrix(unlist(all_edges_fixed), nrow  = 2)),
                                        fixed_none_edges = t(matrix(unlist(all_non_edges_fixed), nrow  = 2)),
                                        further_edge = row_col_inner,
                                        depth_premises = depth_premises)
  df_upper_bounds[nrow(df_upper_bounds) + 1, ] <- list(list(bounds$fixed_edges), list(bounds$fixed_none_edges),
                                                       list(bounds$further_edge), bounds$upper_has_edge, bounds$upper_has_not_edge)
}


not_set_edge <- which(df_upper_bounds$upper_has_edge < max_depth_value)
set_edge <- which(df_upper_bounds$upper_has_not_edge < max_depth_value)
df_upper_bounds[not_set_edge, ]
# The following edges are not set in the maximal depth, as these edges lead to a lower depth then
# the current lower bound of max depth
# length(not_set_edge)
# for (i in not_set_edge) {print(df_upper_bounds[i, ]$further_edge)}
# Round 0: c(1,3), c(1,6), c(1,8), c(3,2), c(3,6), c(5,2), c(5,3), c(5,6), c(5,8), c(6,2), c(7,2)  # in total 11
# Round 1: empty
df_upper_bounds[set_edge, ]
# The following edges need to be set, else the depth value would be below the current max depth
# Note that in each row I just add those which we did not fix befor
# length(set_edge)
# for (i in set_edge) {print(df_upper_bounds[i, ]$further_edge)}
# Round 0: c(2,1), c(2,4) # in total 2
# Round 1: empty


## Thus, procedure ends after Round 1.

###### Procedure: Check which edges are fixed by upper and lower bound discussions
## END




# Now we deleted obtained all edges/none edges which are fixed by the upper
# and lower bound of the depth (see Appendix)

sum(unlist(lapply(X = seq(1,52), FUN = function(x) {choose(40,x)}))) # still too large?

max_depth_value <- max(depth_premises$depth_ufg)

poset_basic <- diag(8)
poset_basic[rbind(c(2,1), c(2,4))] <- 1
poset_basic <- ddandrda::compute_transitive_hull(poset_basic)

edges_unclear <- unlist(lapply(seq_len(8), FUN = help_fct_cartesian), recursive = FALSE)
edges_unclear <- edges_unclear[-c(1, 10, 19, 28, 37, 46, 55, 64)]
edges_unclear <- edges_unclear[-c(1, 23, 51)]
edges_unclear <- edges_unclear[-c(7, 9)]
edges_unclear <- edges_unclear[-c(1, 4, 6, 13, 16, 26, 27, 29, 31, 33, 40)]

poset_list <- list()
depth_values_list <- list()

# Note cardinality 1 to 7: no posets with higher depth than the observed one exists
start_time <- Sys.time()
for (card_sub in 1:40) {
  # iterating over all subsets of size k (for fixed n)
  # using Gosper's Hack
  # see: http://programmingforinsomniacs.blogspot.com/2018/03/gospers-hack-explained.html
  subset_binary <- c(rep(0, (length(edges_unclear) - card_sub)), rep(1, card_sub))
  print(paste0("now at cardinality ", card_sub))
  while (TRUE) {

    poset_inner <- poset_basic
    set_edges <- which(as.logical(subset_binary))
    for (edge in set_edges) {
      poset_inner[edges_unclear[[edge]][1], edges_unclear[[edge]][2]] <- 1
    }
    if (ddandrda::test_if_porder(poset_inner)) {
      depth_inner <- compute_ufg_given_premises(poset_inner, ufg = depth_premises$ufg_premises,
                                                cn_constant = depth_premises$constant_cn)
      if (depth_inner >= max_depth_value) {
        poset_list <- append(poset_list, list(poset_inner))
        depth_value_list <- append(depth_value_list, list(depth_inner))
      }
    }

    # switch to next subset or break while loop
    if (all(subset_binary[seq(1, card_sub)] == rep(1, card_sub))) {
      break
    }
    index_one <- which(subset_binary == 1)
    max_one <- max(index_one)
    max_zero_s_max_one <- max(
      which(subset_binary == 0)[which(subset_binary == 0) < max_one])
    subset_binary[c(max_zero_s_max_one, max_zero_s_max_one + 1)] <- c(1,0)
    ones_larger <- index_one[index_one > max_zero_s_max_one + 1]
    if (length(ones_larger) != 0) {
      subset_binary[seq(min(ones_larger), length(edges_unclear))] <- 0
      subset_binary[seq(length(edges_unclear) - length(ones_larger) + 1,
                        length(edges_unclear))] <- rep(1, length(ones_larger))
    }
  }
}
end_time <- Sys.time()
end_time - start_time



#
#
# ## Some pictures and vlaues
# print(paste0("The minimal value in list is ", min(depth_poset_list$depth_ufg), ". Note that the smalles value of ALL posets is 0."))
# print(paste0("The maximal value is ", max(depth_poset_list$depth_ufg)))
# #### mean, sd, distribution etc value are not meaningful as not all posets are considered!
#
#
#
#
#
# ## All partial orders with depth larger than 0
# max_depth_index <- sort(depth_poset_list$depth_ufg, index.return = TRUE, decreasing = TRUE)$ix
# stop_print <- 100
# pdf("plots_all_from_highest_depth.pdf", onefile = TRUE)
# for (i in max_depth_index[seq(1,stop_print)]) {
#   mat <- matrix(as.logical(poset_list[[i]]), ncol = item_number)
#   colnames(mat) <- rownames(mat) <- names_columns
#   hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))
# }
# dev.off()
#
#
#
# ## Intersections (high to low)
# max_depth_index <- sort(depth_poset_list$depth_ufg, index.return = TRUE, decreasing = TRUE)$ix
# pdf("plots_all_intersect_from_highes.pdf", onefile = TRUE)
# stop_print <- 100
# for (i in 1:min(length(max_depth_index), stop_print)) {
#   intersect <- matrix(rep(TRUE, item_number * item_number), ncol = item_number)
#   for (j in seq(1, i)) {
#     intersect <- intersect & matrix(as.logical(poset_list[[max_depth_index[j]]]), ncol = item_number)
#   }
#   colnames(intersect) <- rownames(intersect) <- names_columns
#   # print(mat * 1)
#   # hasse(mat) plots the graph from top to bottom, with smallest value at bottom
#   # -> change and set arrow to "backward"
#   hasse(t(intersect), parameters = list(arrow = "backward", shape = "roundrect"))
# }
# dev.off()
#
#
# # minimal to maximal depth value does not make sense at this needs to contain all
# # zeros which lead directly to an empty poset.





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



