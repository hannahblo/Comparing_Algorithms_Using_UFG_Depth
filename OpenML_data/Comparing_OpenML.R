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
  ## computes weighthed representation of a data matrix x with duplicated rows,
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







################################################################################
# Data Set
################################################################################

setwd("comptutation_results/")


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

# learner_unique <- unique(data_final_filter$learner.name)
# data_set_unique <- unique(data_final_filter$data.name)

# View(data_final_filter)
# dim(data_final_filter)
# colnames(data_final_filter)
# length(learner_unique)
# learner_unique
# length(data_set_unique)

# We are only interested in the following performance measures

# @ALL: Hier alle performance maße einfügen die euch interessieren!
data_set_eval <- data_final_filter[, c("data.name", "learner.name",
                                       "f.measure", "predictive.accuracy",
                                       "area.under.roc.curve", # Brier Score
                                       "root.mean.squared.error",
                                       # die folgenden sind neu im vgl zu ISIPTA
                                       "mean.absolute.error",
                                       "kappa",
                                       "kb.relative.information.score"
                                       )]
# In contrast to the other performance measure, lower root.mean.squared.error
# is better.
data_set_eval[ ,"root.mean.squared.error"] <-
  1 - data_set_eval[ ,"root.mean.squared.error"]
data_set_eval[ ,"mean.absolute.error"] <-
  1 - data_set_eval[ ,"mean.absolute.error"]

# We are only interested in the following classifiers
data_set_eval <- data_set_eval[data_set_eval$learner.name %in%
                                 c("classif.ranger",
                                   "classif.rpart",
                                   "classif.multinom",
                                   "classif.kknn",
                                   "classif.glmnet"), ]


list_mat_porders_ml <- list()
number_classifiers <- 5
for (i in seq(1, length(unique(data_final_filter$data.name)))) {
  list_mat_porders_ml[i] <- list(convert_to_matrix(data_set_eval[seq((i - 1) * number_classifiers + 1,
                                                                     i * number_classifiers), ]))
}


# saveRDS(list_mat_porders_ml, "list_mat_porders_ml.rds")
# saveRDS(data_set_eval, "data_set_eval.rds")

## contingency table
cont_table <- Reduce("+", list_mat_porders_ml)


# @ALL: Hier können die Performance Measures in 2 Gruppen unterteilt werden

# Above are the posets based on all considered performance measures. Now we
# divide the performance measure and obtain for each set 80 posets

# based on Sensitivity and specificity
perform_measure_1 <- c("data.name", "learner.name",
                       "f.measure", "predictive.accuracy",
                       "area.under.roc.curve" # Brier Score
                       )
# based on L1 or L2 loss
perform_measure_2 <- c("data.name", "learner.name",
                       "root.mean.squared.error",
                       "mean.absolute.error")
# Note that kappa and kb.relative.information.score are ignored here

list_mat_porders_ml_divi1 <- list()
list_mat_porders_ml_divi2 <- list()

number_classifiers <- 5
for (i in seq(1, length(unique(data_final_filter$data.name)))) {
  list_mat_porders_ml_divi1[i] <- list(convert_to_matrix(data_set_eval[seq((i - 1) * number_classifiers + 1,
                                                                     i * number_classifiers), perform_measure_1]))
  list_mat_porders_ml_divi2[i] <- list(convert_to_matrix(data_set_eval[seq((i - 1) * number_classifiers + 1,
                                                                           i * number_classifiers), perform_measure_2]))
}


################################################################################
#
# PART 1: ANALYSIS BASED ON ALL PERFORMANCE MEASURES AT ONCE
#
################################################################################

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



################################################################################
# Computation of the ufg-depth
################################################################################

### Compute the VC dimension
# Formal context given by the partial orders in list_mat
fc_ml_porder <- ddandrda::compute_conceptual_scaling(input_porder = list_mat_porders_ml)
ml_porder_model <- oofos::compute_extent_vc_dimension(fc_ml_porder)
vc_fc_ml_porder <- gurobi::gurobi(ml_porder_model)
vc <- vc_fc_ml_porder$objval # 8



### Compute the ufg-depth
# Preparation of the computation, needed as input of
porder_all <- ddandrda::compute_all_partial_orders(5, list = FALSE, complemented = TRUE)
list_porder_all <- ddandrda::compute_all_partial_orders(5, list = TRUE, complemented = FALSE)

data_context <- get_weighted_representation(fc_ml_porder) # duplication
n_row_context <- nrow(data_context$x_weighted)
count_dup <- data_context$counts
number_obs <- sum(data_context$counts)
list_ml_porder_unique <- ddandrda::convert_context_to_list(data_context$x_weighted[ ,(1:25)],  complemented = FALSE)

whole_context <- rbind(data_context$x_weighted, porder_all) # context of all posets
index <- which(!duplicated(whole_context))
whole_context <- whole_context[index,]



# Computation of S, see article (1)
start_time <- Sys.time()
ufg_premises <- oofos::enumerate_ufg_premises(whole_context, n_row_context)
total_time <- Sys.time() - start_time

# saveRDS(total_time, "total_time.rds")
# saveRDS(ufg_premises, "ufg_premises.rds")
# length(ufg_premises)

# TODO @Hannah in ddandrda programmieren. Das folgende braucht man doch öfters als gedacht
# Insbesonder muss der ganze Code dazu angepasst werden

# ufg depth computation
emp_prob <- count_dup / number_obs
depth_ufg <- rep(0, length(list_ml_porder_unique))
constant_c <- 0

for (i in 1:length(ufg_premises)) {
  # print(paste0("Iteration ", i,  " of ", dim(ufg_premises)[1]))
  index_premise <- ufg_premises[[i]]
  if (length(index_premise) < 2) {
    print(paste0("cardinaltiy ufg_premise is ", length(index_premise)))
  }

  prod_emp_ufg <- prod(emp_prob[index_premise])
  concl_ufg <- ddandrda::test_porder_in_concl(list_ml_porder_unique[index_premise], list_ml_porder_unique) * 1

  depth_ufg <- depth_ufg + concl_ufg * prod_emp_ufg
  constant_c <- constant_c + prod_emp_ufg
}

depth_value <- depth_ufg / constant_c


# Adding duplicate values
depth_value_all <- c()
list_data_all <- vector("list", sum(count_dup))
saving <- 1
for (i in 1:length(depth_value)) {
  for (j in 1:count_dup[i]) {
    list_data_all[[saving]] <- list_ml_porder_unique[[i]]
    saving <- saving + 1
  }
  depth_value_all <- append(depth_value_all, rep(depth_value[i], count_dup[i]))

}


# saveRDS(constant_c, "constant_c.rds")
# saveRDS(depth_ufg, "ufg_depth.rds")
# saveRDS(vc, "vc.rds")
# saveRDS(depth_value_all, "depth_values.rds")





################################################################################
# Descriptive Analysis -> duplicates
################################################################################
names_columns <- c("multinom", "ranger", "rpart", "glmnet", "kknn")
item_number <- dim(list_ml_porder_unique[[1]])[1]

sort(count_dup, index.return = TRUE, decreasing = TRUE)
mat <- matrix(as.logical(list_ml_porder_unique[[27]]), ncol = item_number)
colnames(mat) <- rownames(mat) <- names_columns
hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))



################################################################################
# Descriptive Analysis: minimal, maximal ufg depth
################################################################################
print(paste0("The minimal value is ", min(depth_value_all)))
print(paste0("The maximal value is ", max(depth_value_all)))
print(paste0("The mean value is ", mean(depth_value_all)))
print(paste0("The standard deviation is ", sd(depth_value_all)))
print(paste0("The median is ", median(depth_value_all)))
print(paste0("The number of depth value duplicates (reduced by duplicates given by the data) are ", length(depth_value) -
               length(unique(depth_value))))

### Distribution of Depth Values
pdf("boxplot_depth.pdf", onefile = TRUE)
boxplot(depth_value_all, main = "Boxplot of the depth values")
dev.off()



## All partial orders
max_depth_index <- sort(depth_value_all, index.return = TRUE, decreasing = TRUE)$ix
pdf("plots_all_from_highest_to_lowest.pdf", onefile = TRUE)
for (i in max_depth_index) {
  mat <- matrix(as.logical(list_data_all[[i]]), ncol = item_number)
  colnames(mat) <- rownames(mat) <- names_columns
  # print(mat * 1)
  # hasse(mat) plots the graph from top to bottom, with smallest value at bottom
  # -> change and set arrow to "backward"
  hasse(t(mat), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()



## Intersections (high to low)
max_depth_index <- sort(depth_value_all, index.return = TRUE, decreasing = TRUE)$ix
pdf("plots_intersect_from_highest_to_lowest.pdf", onefile = TRUE)
for (i in 1:length(max_depth_index)) {
  intersect <- matrix(rep(TRUE, item_number * item_number), ncol = item_number)
  for (j in seq(1, i)) {
    intersect <- intersect & matrix(as.logical(list_data_all[[max_depth_index[j]]]), ncol = item_number)
  }
  colnames(intersect) <- rownames(intersect) <- names_columns
  # print(mat * 1)
  # hasse(mat) plots the graph from top to bottom, with smallest value at bottom
  # -> change and set arrow to "backward"
  hasse(t(intersect), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()



## Intersections (low to high)
min_depth_index <- sort(depth_value_all, index.return = TRUE, decreasing = FALSE)$ix
pdf("plots_intersect_from_lowest_to_highest.pdf", onefile = TRUE)
for (i in 1:length(min_depth_index)) {
  intersect <- matrix(rep(TRUE, item_number * item_number), ncol = item_number)
  for (j in seq(1, i)) {
    intersect <- intersect & matrix(as.logical(list_data_all[[min_depth_index[j]]]), ncol = item_number)
  }
  colnames(intersect) <- rownames(intersect) <- names_columns
  # print(mat * 1)
  # hasse(mat) plots the graph from top to bottom, with smallest value at bottom
  # -> change and set arrow to "backward"
  hasse(t(intersect), parameters = list(arrow = "backward", shape = "roundrect"))
}
dev.off()





################################################################################
# Descriptive Analysis: dispersion
################################################################################

# compute ufg for all depth values
emp_prob_all <- count_dup / number_obs
depth_ufg_all <- rep(0, length(list_porder_all))
constant_c_all <- 0

for (i in 1:length(ufg_premises)) {
  print(paste0("Iteration ", i,  " of ", length(ufg_premises)))
  index_premise_all <- ufg_premises[[i]]

  prod_emp_ufg_all <- prod(emp_prob_all[index_premise_all])
  concl_ufg_all <- ddandrda::test_porder_in_concl(list_ml_porder_unique[index_premise_all], list_porder_all) * 1

  depth_ufg_all <- depth_ufg_all + concl_ufg_all * prod_emp_ufg_all
  constant_c_all <- constant_c_all + prod_emp_ufg_all
}

depth_value_all <- depth_ufg_all / constant_c_all
unique(depth_value_all)
length(unique(depth_value_all))
length(list_porder_all)

# compute proportion
quantiles_ufg <- quantile(depth_value, probs = c(0.25, 0.5, 0.75))

# Note that the 25$ highest depth value corresponds to the 75 quantile
proportion_75 <- length(which(depth_value_all >= quantiles_ufg[1])) / length(depth_value_all)
proportion_50 <- length(which(depth_value_all >= quantiles_ufg[2])) / length(depth_value_all)
proportion_25 <- length(which(depth_value_all >= quantiles_ufg[3])) / length(depth_value_all)
proportion_75





################################################################################
#
# PART 2: COMPARING POSETS BASED ON DIFFERENT PERFORMANCE MEASURES
#
################################################################################

################################################################################
# Descriptive analysis of existence of edges (Step 1: not with depth function)
################################################################################

### Which edge exists
list_mat_interest <- list_mat_porders_ml_divi1 # list_mat_porders_ml_divi2

length(list_mat_interest)
length(unique(list_mat_interest))
Reduce("|", list_mat_interest)
Reduce("&", list_mat_interest)
Reduce("+", list_mat_interest)

edges <- Reduce("+", list_mat_interest)
colnames(edges) <- rownames(edges) <- c("LR", "RF", "CART", "LASSO", "KNN")
df_edge_exist <- melt(edges)
df_edge_exist <- df_edge_exist[df_edge_exist$value != 0, ]

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




################################################################################
# Computation of the ufg-depths
################################################################################
list_mat_interest <- list_mat_porders_ml_divi1 # list_mat_porders_ml_divi2

### Compute the VC dimension
# Formal context given by the partial orders in list_mat
fc_ml_porder_divi <- ddandrda::compute_conceptual_scaling(input_porder = list_mat_interest)
ml_porder_model_divi <- oofos::compute_extent_vc_dimension(fc_ml_porder_divi)
vc_fc_ml_porder_divi <- gurobi::gurobi(ml_porder_model_divi)
vc <- vc_fc_ml_porder_divi$objval # 8




### Compute the ufg-depth
# Preparation of the computation, needed as input of
porder_all <- ddandrda::compute_all_partial_orders(5, list = FALSE, complemented = TRUE)
list_porder_all <- ddandrda::compute_all_partial_orders(5, list = TRUE, complemented = FALSE)

data_context_divi <- get_weighted_representation(fc_ml_porder_divi) # duplication
n_row_context_divi <- nrow(data_context_divi$x_weighted)
count_dup_divi <- data_context_divi$counts
number_obs_divi <- sum(data_context_divi$counts)

list_ml_porder_unique_divi <- ddandrda::convert_context_to_list(data_context_divi$x_weighted[ ,(1:25)],  complemented = FALSE)

whole_context_divi <- rbind(data_context_divi$x_weighted, porder_all) # context of all posets
index_divi <- which(!duplicated(whole_context_divi))
whole_context_divi <- whole_context_divi[index_divi ,]



# Computation of S, see article (1)
start_time <- Sys.time()
ufg_premises_divi <- oofos::enumerate_ufg_premises(whole_context_divi, n_row_context_divi)
# ufg_premises_divi2 <- oofos::enumerate_ufg_premises(whole_context_divi, n_row_context_divi)
total_time <- Sys.time() - start_time


# ufg depth computation
emp_prob_divi <- count_dup_divi / number_obs_divi
depth_ufg_divi <- rep(0, length(unlist(list(list_mat_porders_ml_divi1, list_mat_porders_ml_divi2), recursive = FALSE)))
constant_c_divi <- 0

for (i in 1:length(ufg_premises_divi)) {
  # print(paste0("Iteration ", i,  " of ", dim(ufg_premises)[1]))
  index_premise <- ufg_premises_divi[[i]]
  if (length(index_premise) < 2) {
    print(paste0("cardinaltiy ufg_premise is ", length(index_premise)))
  }

  prod_emp_ufg <- prod(emp_prob_divi[index_premise])
  concl_ufg <- ddandrda::test_porder_in_concl(list_ml_porder_unique_divi[index_premise],
                                              obj_porder_obs = unlist(list(list_mat_porders_ml_divi1, list_mat_porders_ml_divi2), recursive = FALSE)) * 1

  depth_ufg_divi <- depth_ufg_divi + concl_ufg * prod_emp_ufg
  constant_c_divi <- constant_c_divi + prod_emp_ufg
}

depth_value_divi <- depth_ufg_divi / constant_c_divi


# # Adding duplicate values
# depth_value_all_divi <- c()
# list_data_all_divi <- vector("list", sum(count_dup_divi))
# saving <- 1
# for (i in 1:length(depth_value_divi)) {
#   for (j in 1:count_dup[i]) {
#     list_data_all_divi[[saving]] <- list_ml_porder_unique_divi[[i]]
#     saving <- saving + 1
#   }
#   depth_value_all_divi <- append(depth_value_all_divi, rep(depth_value_divi[i], count_dup_divi[i]))
#
# }


depth_value_all_divi1 <- depth_value_divi
# depth_value_all_divi2 <- depth_value_divi # oben anders setzen!!!

plot(depth_value_all_divi1, depth_value_all_divi2)

# saveRDS(constant_c, "constant_c.rds")
# saveRDS(depth_ufg, "ufg_depth.rds")
# saveRDS(vc, "vc.rds")
# saveRDS(depth_value_all, "depth_values.rds")













