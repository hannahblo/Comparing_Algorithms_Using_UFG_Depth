################################################################################
# Set Up
################################################################################

### Information about packages ddandrda. This package is under developement on
# git. Installation can be done by:
# remove.packages("ddandrda")
# install.packages("devtools")
# devtools::install_github("hannahblo/ddandrda")

### Information about packages oofos. This package is under developement on git.
# Note that to install this package, the R-package gurobi is needed. This is an
# remove.packages("oofos")
# install.packages("devtools")
# devtools::install_github("schollmeyer/oofos")

### All the other R-packages are on CRAN (07.02.2023)

library(ddandrda)
library(oofos)
library(gurobi)
library(OpenML)
library(dplyr)
library(farff)
library(ggplot2)
library(reshape2)
library(hasseDiagram)

setwd("comptutation_results/")



################################################################################
# Functions needed to perform the analysis
################################################################################

# Function which converts to to partial orders
convert_to_matrix <- function(single_data_eval) {
  list_learner <- single_data_eval$learner.name
  if (length(list_learner) != length(unique(list_learner))) {
    print(single_data_eval[ ,"data.name"])
    stop("Fehler im obigen Datensatz")
  }
  number_learner <- length(list_learner)
  graph_mat <- matrix(rep(0, number_learner * number_learner),
                      nrow = number_learner)
  rownames(graph_mat) <- colnames(graph_mat) <- list_learner
  diag(graph_mat) <- 1

  for (learner_i in list_learner) {
    learner_base <- single_data_eval[
      which(single_data_eval$learner.name == learner_i),
      -c(1, 2)]
    for (learner_j in list_learner) {
      learner_comp <- single_data_eval[
        which(single_data_eval$learner.name == learner_j),
        -c(1, 2)]
      if (all(learner_base < learner_comp)) {
        graph_mat[learner_i, learner_j] <- 1
      }
    }
  }
  return(graph_mat)
}


# Function to compute the weights of the fc
get_weighted_representation <- function(x, y = rep(1, dim(x)[1])) {
  ## computes weighted representation of a data matrix x with duplicated rows,
  ##  returns unique(x) together with counts: how often appears the column,
  # mean_y: mean of y in the set of the duplicated columns
  xd <- data.frame(cbind(x, y))
  names(xd)[1] <- "v1"
  v1 <- "v1"
  p <- dim(x)[2]
  result <- as.matrix(plyr::ddply(xd, names(xd[(1:p)]), dplyr::summarise, count = length(v1), mean.y = mean(y), sum.y = sum(y)))
  x_weighted <- result[, (1:p)]
  colnames(x_weighted) <- colnames(x)
  return(list(x_weighted = x_weighted, y_weighted = result[, p + 3], mean_y = result[, p + 2], counts = result[, p + 1]))
}


# Preparation of computing the ufg premises
prepare_ufg_premises <- function(list_mat_porders_ml,
                                 number_items) {

  fc_ml_porder <- ddandrda::compute_conceptual_scaling(input_porder = list_mat_porders_ml)
  porder_all <- ddandrda::compute_all_partial_orders(number_items, list = FALSE, complemented = TRUE)

  data_context <- get_weighted_representation(fc_ml_porder) # duplication
  n_row_context <- nrow(data_context$x_weighted)
  count_dup <- data_context$counts
  number_obs <- sum(data_context$counts)
  list_porder_premises <- ddandrda::convert_context_to_list(data_context$x_weighted[ ,(1:25)],  complemented = FALSE)

  whole_context <- rbind(data_context$x_weighted, porder_all) # context of all posets
  index <- which(!duplicated(whole_context))
  whole_context <- whole_context[index,]
  return(list(count_dup = count_dup,
              number_obs = number_obs,
              whole_context = whole_context,
              list_porder_premises = list_porder_premises,
              n_row_context = n_row_context))
}



# Computing the ufg depth based on already computed premises
compute_ufg_exist_premises <- function(poset_interest, ufg_premises,
                                       prep_ufg_premises) {

  emp_prob <- prep_ufg_premises$count_dup / prep_ufg_premises$number_obs
  depth_ufg <- rep(0, length(poset_interest))
  constant_c <- 0

  for (i in 1:length(ufg_premises)) {
    # print(paste0("Iteration ", i,  " of ", dim(ufg_premises)[1]))
    index_premise <- ufg_premises[[i]]
    prod_emp_ufg <- prod(emp_prob[index_premise])
    concl_ufg <- ddandrda::test_porder_in_concl(prep_ufg_premises$list_porder_premises[index_premise], poset_interest) * 1

    depth_ufg <- depth_ufg + concl_ufg * prod_emp_ufg
    constant_c <- constant_c + prod_emp_ufg
  }

  depth_value <- depth_ufg / constant_c

  return(depth_value)

}


# this plot returns summaries of the result
plot_result <- function(names_columns, item_number, depth_value,
                        list_mat_porders_ml,
                        max_plot_number = 100,
                        file_name_add) {

  max_plot_number <- min(max_plot_number, length(list_mat_porders_ml))

  sink(file = paste0("part", file_name_add, "_summery.txt"))
  print(paste0("The minimal value is ", min(depth_value)))
  print(paste0("The maximal value is ", max(depth_value)))
  print(paste0("The mean value is ", mean(depth_value)))
  print(paste0("The standard deviation is ", sd(depth_value)))
  print(paste0("The median is ", median(depth_value)))
  print(paste0("The number of depth value duplicates (reduced by duplicates given by the data) are ", length(depth_value) -
                 length(unique(depth_value))))
  sink(file = NULL)

  ### Distribution of Depth Values
  pdf(paste0("part", file_name_add,"_boxplot_depth.pdf"), onefile = TRUE)
  boxplot(depth_value, main = "Boxplot of the depth values")
  dev.off()



  ## partial orders
  max_depth_index <- sort(depth_value, index.return = TRUE, decreasing = TRUE)$ix
  pdf(paste0("part", file_name_add, "_plots_from_highest_to_lowest.pdf"), onefile = TRUE)
  for (i in max_depth_index[seq(1, max_plot_number)]) {
    mat <- matrix(as.logical(list_mat_porders_ml[[i]]), ncol = item_number)
    colnames(mat) <- rownames(mat) <- names_columns
    hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))
  }
  dev.off()



  ## Intersections (high to low)
  max_depth_index <- sort(depth_value, index.return = TRUE, decreasing = TRUE)$ix
  pdf(paste0("part", file_name_add, "_plots_intersect_from_highest_to_lowest.pdf"), onefile = TRUE)
  for (i in 1:max_plot_number) {
    intersect <- matrix(rep(TRUE, item_number * item_number), ncol = item_number)
    for (j in seq(1, i)) {
      intersect <- intersect & matrix(as.logical(list_mat_porders_ml[[max_depth_index[j]]]), ncol = item_number)
    }
    colnames(intersect) <- rownames(intersect) <- names_columns
    hasse(t(intersect), parameters = list(arrow = "backward", shape = "roundrect"))
  }
  dev.off()



  ## Intersections (low to high)
  min_depth_index <- sort(depth_value, index.return = TRUE, decreasing = FALSE)$ix
  pdf(paste0("part", file_name_add,"_plots_intersect_from_lowest_to_highest.pdf"), onefile = TRUE)
      for (i in 1:max_plot_number) {
        intersect <- matrix(rep(TRUE, item_number * item_number), ncol = item_number)
        for (j in seq(1, i)) {
          intersect <- intersect & matrix(as.logical(list_mat_porders_ml[[min_depth_index[j]]]), ncol = item_number)
        }
        colnames(intersect) <- rownames(intersect) <- names_columns
        hasse(t(intersect), parameters = list(arrow = "backward", shape = "roundrect"))
      }
    dev.off()

}

################################################################################
# ML Data Set Preparation
################################################################################
# Load the data set from OpenML (see https://www.openml.org/)
data_all <- OpenML::listOMLDataSets()
rel_datasets <- data_all %>% filter(status == "active" &
                                      number.of.classes == 2 &
                                      number.of.features < 100 &
                                      number.of.instances < 1000 &
                                      number.of.instances > 100 &
                                      number.of.instances.with.missing.values == 0 &
                                      max.nominal.att.distinct.values < 5
)

#### flows 2333 (rpart), 2330 (ranger), 2409 (knn), 2337 (xgboost), 2408 (glmnet), 2336(svm), 2317(logit), 2313 lda
# test = getOMLFlow(flow.id = 2333)
#test = getOMLTask(task.id = 3729)
#### 4689


flows_for_paper <- c(2333, 2330, 2317, 2337, 2408, 2336,  2409)
outls <- list()
for (i in 1:length(flows_for_paper)) {
  temp <- listOMLRunEvaluations(flow.id = flows_for_paper[i], limit = 10000)
  if (i == 1) {
    datasets = temp$data.name
  } else {
    datasets = intersect(datasets,temp$data.name)
  }
  outls[[i]] = temp
}
data_final <- do.call('rbind', outls)
data_final <-  data_final[data_final$data.name %in% datasets,]
data_final <- data_final %>% group_by(flow.id, data.name) %>% slice(n())



extractDataId <- function(taskid){
  print(taskid)
  if (length(tryCatch({res = getOMLTask(task.id = taskid)}, error = function(e){return(NA)})) <= 1) {
    return(NA)
  }
  return(res$input$data.set$desc$id)
}
data_final$data.id = sapply(data_final$task.id,function(x)extractDataId(taskid = x))


# data_final[which(is.na(data_final$data.id)), ]
# taskid = 3019: Data set has been deactivated. Thus delete this one
if (length(which(is.na(data_final$data.id))) > 0) {
  data_final <- data_final[-which(is.na(data_final$data.id)), ]
}

# We are interested in binary classification with sample size between 450 and
# 10000
data_final_filter = data_final %>%
  group_by(data.id) %>%
  dplyr::mutate(count = n()) %>%
  left_join(data_all, by = "data.id") %>%
  filter((count == length(flows_for_paper)) &
           (number.of.instances.x > 450) &
           (number.of.instances.x < 10000) & # 10000 500
           (number.of.classes == 2)
  )
data_final_filter <- data_final_filter[order(data_final_filter$data.name), ]

# saveRDS(data_final_filter, "data_final_filter.rds")
# data_final_filter <- readRDS("data_final_filter.rds")


################################################################################
#
# PART 1: ANALYSIS BASED ON ALL PERFORMANCE MEASURES AT ONCE
#
################################################################################

################################################################################
# Data Set
################################################################################


# We are only interested in the following performance measures and on the
# following ml classifiers

########
## ISIPTA Article:
data_set_eval <- data_final_filter[, c("data.name", "learner.name",
                                       "f.measure", "predictive.accuracy",
                                       "area.under.roc.curve", # Brier Score
                                       "root.mean.squared.error")]
# In contrast to the other performance measure, lower root.mean.squared.error
# is better.
data_set_eval[ ,"root.mean.squared.error"] <-
  1 - data_set_eval[ ,"root.mean.squared.error"]

# We are only interested in the following classifiers
data_set_eval <- data_set_eval[data_set_eval$learner.name %in%
                                 c("classif.ranger",
                                   "classif.rpart",
                                   "classif.multinom",
                                   "classif.kknn",
                                   "classif.glmnet"), ]
number_classifiers <- 5



########
## IJAR Article ueberlegungen
# data_set_eval <- data_final_filter[, c("data.name", "learner.name",
#                                        "f.measure", "predictive.accuracy",
#                                        "area.under.roc.curve", # Brier Score
#                                        "root.mean.squared.error",
#                                        "mean.absolute.error",
#                                        "kappa",
#                                        "kb.relative.information.score"
# )]
# # In contrast to the other performance measure, lower root.mean.squared.error
# # is better.
# data_set_eval[ ,"root.mean.squared.error"] <-
#   1 - data_set_eval[ ,"root.mean.squared.error"]
# data_set_eval[ ,"mean.absolute.error"] <-
#   1 - data_set_eval[ ,"mean.absolute.error"]
#
# # We are only interested in the following classifiers
# data_set_eval <- data_set_eval[data_set_eval$learner.name %in%
#                                  c("classif.ranger",
#                                    "classif.rpart",
#                                    "classif.multinom",
#                                    "classif.kknn",
#                                    "classif.glmnet"), ]
# number_classifiers <- 5



# Convert to porders
list_mat_porders_ml <- list()

for (i in seq(1, length(unique(data_final_filter$data.name)))) {
  list_mat_porders_ml[i] <- list(convert_to_matrix(data_set_eval[seq((i - 1) * number_classifiers + 1,
                                                                     i * number_classifiers), ]))
}


# saveRDS(list_mat_porders_ml, "part1_list_mat_porders_ml.rds")




################################################################################
# Descriptive analysis of existence of edges (Step 1: not with depth function)
################################################################################

### Which edge exists
length(list_mat_porders_ml) # 80
length(unique(list_mat_porders_ml)) # 58
Reduce("|", list_mat_porders_ml)
Reduce("&", list_mat_porders_ml)
Reduce("+", list_mat_porders_ml)


edges <- Reduce("+", list_mat_porders_ml)
colnames(edges) <- rownames(edges) <- c("LR", "RF", "CART", "LASSO", "KNN")
df_edge_exist <- melt(edges)
df_edge_exist <- df_edge_exist[df_edge_exist$value != 0, ]

pdf("part1_heightmap.pdf", onefile = TRUE)
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


################################################################################
# Computation of the ufg-depth
################################################################################



# Computation of S, see article (1)
prep_ufg_premises <- prepare_ufg_premises(list_mat_porders_ml, number_items = number_classifiers)
start_time <- Sys.time()
ufg_premises <- oofos::enumerate_ufg_premises(prep_ufg_premises$whole_context, prep_ufg_premises$n_row_context)
total_time <- Sys.time() - start_time

# saveRDS(total_time, "part1_total_time.rds")
# saveRDS(ufg_premises, "part1_ufg_premises.rds")
# length(ufg_premises)


# Computation of depth of observed posets
depth_value <- compute_ufg_exist_premises(poset_interest = prep_ufg_premises$list_porder_premises,
                                          ufg_premises,
                                          prep_ufg_premises)

depth_value_dupl <- compute_ufg_exist_premises(poset_interest = list_mat_porders_ml,
                                               ufg_premises,
                                               prep_ufg_premises)


# Computation of depth af all posets
list_porder_all <- ddandrda::compute_all_partial_orders(number_classifiers, list = TRUE, complemented = FALSE)
depth_value_all <- compute_ufg_exist_premises(poset_interest = list_porder_all,
                                              ufg_premises,
                                              prep_ufg_premises)


# saveRDS(depth_value, "part1_depth_value.rds")
# saveRDS(depth_value_dupl, "part1_depth_value_dupl.rds")
# saveRDS(depth_value_all, "part1_depth_value_all.rds")




################################################################################
# Descriptive Analysis: minimal, maximal, intersection etc
################################################################################
# observed data
plot_result(names_columns =  c("multinom", "ranger", "rpart", "glmnet", "kknn"),
            item_number = dim(prep_ufg_premises$list_porder_premises[[1]])[1],
            depth_value =  depth_value_dupl,
            list_mat_porders_ml = list_mat_porders_ml,
            file_name_add = "1_obs")

# all possible posets
plot_result(names_columns =  c("multinom", "ranger", "rpart", "glmnet", "kknn"),
            item_number = dim(prep_ufg_premises$list_porder_premises[[1]])[1],
            depth_value =  depth_value_all,
            list_mat_porders_ml = list_porder_all,
            max_plot_number = 50,
            file_name_add = "1_all")



################################################################################
# Descriptive Analysis: dispersion
################################################################################

# compute proportion
quantiles_ufg <- quantile(depth_value, probs = c(0.25, 0.5, 0.75))

# Note that the 25$ highest depth value corresponds to the 75 quantile
proportion_75 <- length(which(depth_value_all >= quantiles_ufg[1])) / length(depth_value_all)
proportion_50 <- length(which(depth_value_all >= quantiles_ufg[2])) / length(depth_value_all)
proportion_25 <- length(which(depth_value_all >= quantiles_ufg[3])) / length(depth_value_all)
proportion_25
# 75 ->  0.2562042, 50 -> 0.09974001, 25 -> 0.01890806








################################################################################
#
# PART 2: COMPARING POSETS BASED ON DIFFERENT PERFORMANCE MEASURES
#
################################################################################

################################################################################
# Data Set
################################################################################
# Now we produce two sets of posets based on the same data set and classifiers
# These two set of posets only differ in the used performance measures to obtain
# the posets

# The division of the performance measures is done based on  C. Ferri *,
# J. Hernández-Orallo, R. Modroiu (2008): An experimental comparison of
# performance measures for classification

data_set_eval_2 <- data_final_filter[, c("data.name", "learner.name",
                                       "f.measure", "predictive.accuracy",
                                       "area.under.roc.curve",
                                       "root.mean.squared.error"
)]
# In contrast to the other performance measure, lower root.mean.squared.error
# is better.
data_set_eval_2[ ,"root.mean.squared.error"] <-
  1 - data_set_eval_2[ ,"root.mean.squared.error"]


# based on Sensitivity and specificity
perform_measure_21 <- c("data.name", "learner.name",
                       "f.measure", "predictive.accuracy"
                       )

# based on L1 or L2 loss
perform_measure_22 <- c("data.name", "learner.name",
                        "area.under.roc.curve",
                        "root.mean.squared.error")  # Brier Score
# Note that kappa and kb.relative.information.score are ignored here

list_mat_porders_ml_21 <- list()
list_mat_porders_ml_22 <- list()

number_classifiers <- 5
for (i in seq(1, length(unique(data_final_filter$data.name)))) {
  list_mat_porders_ml_21[i] <- list(convert_to_matrix(data_set_eval_2[seq((i - 1) * number_classifiers + 1,
                                                                     i * number_classifiers), perform_measure_21]))
  list_mat_porders_ml_22[i] <- list(convert_to_matrix(data_set_eval_2[seq((i - 1) * number_classifiers + 1,
                                                                           i * number_classifiers), perform_measure_22]))
}

# saveRDS(list_mat_porders_ml_21, "part2_list_mat_porders_ml_21.rds")
# saveRDS(list_mat_porders_ml_22, "part2_list_mat_porders_ml_22.rds")

# list_mat_porders_ml_21 <- readRDS("part2_list_mat_porders_ml_21.rds")
# list_mat_porders_ml_22 <- readRDS("part2_list_mat_porders_ml_22.rds")

################################################################################
# Descriptive analysis of existence of edges and poset difference
################################################################################

### Which edge exists
current_interest <- list_mat_porders_ml_21 # list_mat_porders_ml_22

length(current_interest) # 80
length(unique(current_interest)) # 58
Reduce("|", current_interest)
Reduce("&", current_interest)
Reduce("+", current_interest)


edges <- Reduce("+", current_interest)
colnames(edges) <- rownames(edges) <- c("LR", "RF", "CART", "LASSO", "KNN")
df_edge_exist <- melt(edges)
df_edge_exist <- df_edge_exist[df_edge_exist$value != 0, ]

pdf("part2_heightmap_21.pdf", onefile = TRUE) # if necessary adjust to 22
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


#### Does any two poset equal?
same_posets <- lapply(seq(1, length(list_mat_porders_ml_21)), FUN = function(x) {all((list_mat_porders_ml_21[[x]] - list_mat_porders_ml_22[[x]]) == 0)})
which(unlist(same_posets))
length(which(unlist(same_posets)))

################################################################################
# Computation of the ufg-depth
################################################################################

# Computation of S, see article (1)
prep_ufg_premises_21 <- prepare_ufg_premises(list_mat_porders_ml_21, number_items = number_classifiers)
prep_ufg_premises_22 <- prepare_ufg_premises(list_mat_porders_ml_22, number_items = number_classifiers)
start_time <- Sys.time()
ufg_premises_21 <- oofos::enumerate_ufg_premises(prep_ufg_premises_21$whole_context, prep_ufg_premises_21$n_row_context)
ufg_premises_22 <- oofos::enumerate_ufg_premises(prep_ufg_premises_22$whole_context, prep_ufg_premises_22$n_row_context)
total_time <- Sys.time() - start_time

# saveRDS(total_time, "part2_total_time.rds")
# saveRDS(ufg_premises_21, "part2_ufg_premises_21.rds")
# saveRDS(ufg_premises_22, "part2_ufg_premises_22.rds")

# Computation of depth of observed posets
depth_value_21 <- compute_ufg_exist_premises(poset_interest = prep_ufg_premises_21$list_porder_premises,
                                             ufg_premises_21,
                                             prep_ufg_premises_21)

depth_value_dupl_21 <- compute_ufg_exist_premises(poset_interest = list_mat_porders_ml_21,
                                                  ufg_premises_21,
                                                  prep_ufg_premises_21)

depth_value_22 <- compute_ufg_exist_premises(poset_interest = prep_ufg_premises_22$list_porder_premises,
                                             ufg_premises_22,
                                             prep_ufg_premises_22)

depth_value_dupl_22 <- compute_ufg_exist_premises(poset_interest = list_mat_porders_ml_22,
                                                  ufg_premises_22,
                                                  prep_ufg_premises_22)


# Computation of depth of all posets
list_porder_all <- ddandrda::compute_all_partial_orders(number_classifiers, list = TRUE, complemented = FALSE)
depth_value_all_21 <- compute_ufg_exist_premises(poset_interest = list_porder_all,
                                                 ufg_premises_21,
                                                 prep_ufg_premises_21)
depth_value_all_22 <- compute_ufg_exist_premises(poset_interest = list_porder_all,
                                                 ufg_premises_22,
                                                 prep_ufg_premises_22)

# saveRDS(depth_value_21, "part2_depth_value_21.rds")
# saveRDS(depth_value_dupl_21, "part2_depth_value_dupl_21.rds")
# saveRDS(depth_value_all_21, "part2_depth_value_all_21.rds")
#
# saveRDS(depth_value_22, "part2_depth_value_22.rds")
# saveRDS(depth_value_dupl_22, "part2_depth_value_dupl_22.rds")
# saveRDS(depth_value_all_22, "part2_depth_value_all_22.rds")

# depth_value_21 <- readRDS("part2_depth_value_21.rds")
# depth_value_dupl_21 <- readRDS("part2_depth_value_dupl_21.rds")
# depth_value_all_21 <- readRDS("part2_depth_value_all_21.rds")
#
# depth_value_22 <- readRDS("part2_depth_value_22.rds")
# depth_value_dupl_22 <- readRDS("part2_depth_value_dupl_22.rds")
# depth_value_all_22 <- readRDS("part2_depth_value_all_22.rds")

################################################################################
# Analysis of their differences
################################################################################

# difference in the depth value
differences_all <- depth_value_all_21 - depth_value_all_22
pdf("part2_all_boxplot_differences.pdf", onefile = TRUE)
boxplot(differences_all)
dev.off()

differences_obs <- depth_value_dupl_21 - depth_value_dupl_22
pdf("part2_obs_boxplot_differences.pdf", onefile = TRUE)
boxplot(differences_obs)
dev.off()

# DD plot (Attention: iid assumption does not hold)
pdf("part2_all_DD_plot.pdf", onefile = TRUE)
plot(depth_value_all_21, depth_value_all_22)
dev.off()

pdf("part2_obs_DD_plot.pdf", onefile = TRUE)
plot(depth_value_dupl_21, depth_value_dupl_22)
dev.off()




# Rank structure
index_depth_21 <- match(depth_value_all_21, sort(unique(depth_value_all_21)))
index_depth_22 <- match(depth_value_all_22, sort(unique(depth_value_all_22)))
max(index_depth_21)
max(index_depth_22)
any(duplicated(depth_value_all_21))
any(duplicated(depth_value_all_22))
difference_rank <- index_depth_21 - index_depth_22
pdf("part2_all_difference_rank.pdf", onefile = TRUE)
boxplot(difference_rank)
dev.off()

summary(difference_rank)
# Min.     1st Qu.  Median    Mean 3rd Qu.    Max.
# -2134.0  -323.5    -7.0     0.0   323.0  2066.0
length(list_porder_all) # [1] 4231


# check the highest x % ranked posets
index_depth_21_quant <- which(index_depth_21 < quantile(seq(1, length(list_porder_all)), 0.10))
boxplot(difference_rank[index_depth_21_quant])

index_depth_22_quant <- which(index_depth_22 > quantile(seq(1, length(list_porder_all)), 0.95))
boxplot(difference_rank[index_depth_22_quant])


################################################################################
# Descriptive Analysis
################################################################################
# observed data
plot_result(names_columns =  c("multinom", "ranger", "rpart", "glmnet", "kknn"),
            item_number = number_classifiers,
            depth_value =  depth_value_dupl_21,
            list_mat_porders_ml = list_mat_porders_ml_21,
            file_name_add = "2_21_obs")

plot_result(names_columns =  c("multinom", "ranger", "rpart", "glmnet", "kknn"),
            item_number = number_classifiers,
            depth_value =  depth_value_dupl_22,
            list_mat_porders_ml = list_mat_porders_ml_22,
            file_name_add = "2_22_obs")

# all possible posets
plot_result(names_columns =  c("multinom", "ranger", "rpart", "glmnet", "kknn"),
            item_number = number_classifiers,
            depth_value =  depth_value_all_21,
            list_mat_porders_ml = list_porder_all,
            max_plot_number = 50,
            file_name_add = "2_21_all")

plot_result(names_columns =  c("multinom", "ranger", "rpart", "glmnet", "kknn"),
            item_number = number_classifiers,
            depth_value =  depth_value_all_22,
            list_mat_porders_ml = list_porder_all,
            max_plot_number = 50,
            file_name_add = "2_22_all")

################################################################################
#
# PART 3: Discussion on change of ufg depth structure adding performance measures
#
################################################################################

################################################################################
# Data Set
################################################################################
########
## ISIPTA Article:
data_set_eval <- data_final_filter[, c("data.name", "learner.name",
                                       "f.measure", "predictive.accuracy",
                                       "area.under.roc.curve", # Brier Score
                                       "root.mean.squared.error")]
# In contrast to the other performance measure, lower root.mean.squared.error
# is better.
data_set_eval[ ,"root.mean.squared.error"] <-
  1 - data_set_eval[ ,"root.mean.squared.error"]

# We are only interested in the following classifiers
data_set_eval <- data_set_eval[data_set_eval$learner.name %in%
                                 c("classif.ranger",
                                   "classif.rpart",
                                   "classif.multinom",
                                   "classif.kknn",
                                   "classif.glmnet"), ]
number_classifiers <- 5

# View(data_set_eval)
colnames(data_set_eval)
number_measures <- 4
name_measures <- c("f.measure", "predictive.accuracy", "area.under.roc.curve", "root.mean.squared.error")

# First, we compute the correlation between the performance measure
df_corr <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df_corr) <- c("measure_1", "measure_2", "correlation")
measure_corr <- diag(number_measures)
for (i in seq(1, number_measures - 1)) {
  for (j in seq(i + 1, number_measures)) {
    inner_cor <- cor(data_set_eval[[name_measures[[i]]]], data_set_eval[[name_measures[[j]]]], method = 'pearson')
    measure_corr[i,j] <- measure_corr[j,i] <- inner_cor
    df_corr[nrow(df_corr) + 1, ] <- c(name_measures[i], name_measures[j], inner_cor)
  }

}
colnames(measure_corr) <- rownames(measure_corr) <- name_measures
sort_index <- sort(df_corr$correlation, index.return = TRUE,  decreasing = FALSE)$ix

measure_corr

#### Case 1: Based on upper computed correlation
# We consider the following rounds (always uncomment higher rounds please)
# Round 0: two largest correlation vs three largest correlation difference
# Round 1: three larges correlation difference vs four largest correlation difference
perform_measure_31 <- c("data.name", "learner.name",
                        "area.under.roc.curve", "root.mean.squared.error" # Round 0
                        # , "predictive.accuracy" # Round 1
)
perform_measure_32 <- c("data.name", "learner.name",
                        "area.under.roc.curve", "root.mean.squared.error", "predictive.accuracy" # Round 0
                        # , "f.measure" # Round 1

)


# #### Case 2: Based on correlation computed by the Paper
# #### C. Ferri *, J. Hernández-Orallo, R. Modroiu (2008): An experimental
# #### comparison of performance measures for classification
# # We consider the following rounds (always uncomment higher rounds please)
# # Round 0: two largest correlation vs three largest correlation difference
# # Round 1: three larges correlation difference vs four largest correlation difference
# perform_measure_31 <- c("data.name", "learner.name",
#                         "area.under.roc.curve", "predictive.accuracy" # Round 0
#                         # ,  "root.mean.squared.error" # Round 1
# )
# perform_measure_32 <- c("data.name", "learner.name",
#                         "area.under.roc.curve", "predictive.accuracy","root.mean.squared.error"# Round 0
#                         # , "f.measure" # Round 1
#
# )

list_mat_porders_ml_31 <- list()
list_mat_porders_ml_32 <- list()

number_classifiers <- 5
for (i in seq(1, length(unique(data_final_filter$data.name)))) {
  list_mat_porders_ml_31[i] <- list(convert_to_matrix(data_set_eval[seq((i - 1) * number_classifiers + 1,
                                                                        i * number_classifiers), perform_measure_31]))
  list_mat_porders_ml_32[i] <- list(convert_to_matrix(data_set_eval[seq((i - 1) * number_classifiers + 1,
                                                                        i * number_classifiers), perform_measure_32]))
}




################################################################################
# Descriptive analysis of existence of edges and poset difference
################################################################################
#### Does any two poset which equal exist?
same_posets <- lapply(seq(1, length(list_mat_porders_ml_31)), FUN = function(x) {all((list_mat_porders_ml_31[[x]] - list_mat_porders_ml_32[[x]]) == 0)})
all(unlist(same_posets))
which(!unlist(same_posets))

# Observation Case 1:
# For Round 0: all observation equal
# For Round 1: only 4 posets differ: 11 14 21 68

## Observation Case 2:
# For Round 0: all observations equal
# For Roun 1: only 4 posets differ: 11 14 21 68

### Which edge exists
current_interest <- list_mat_porders_ml_31

length(current_interest) # 80
length(unique(current_interest)) # 58
Reduce("|", current_interest)
Reduce("&", current_interest)
Reduce("+", current_interest)


edges <- Reduce("+", current_interest)
colnames(edges) <- rownames(edges) <- c("LR", "RF", "CART", "LASSO", "KNN")
df_edge_exist <- melt(edges)
df_edge_exist <- df_edge_exist[df_edge_exist$value != 0, ]

pdf("part3_heatmap_31_Round0.pdf", onefile = TRUE) # if necessary adjust to 32 and/or Round1
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



