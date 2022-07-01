#test parallel function
#tic()
# data("ecuador")
df <- ecuador

##Training set

nmod = 5
nfolds = 5
nrep = 5

required.packages <- c("sf","terra","data.table","tidyverse",
                       "mlr3","mlr3learners","mlr3spatiotempcv",
                       "mlr3viz","mlr3tuning","iml","future","furrr",
                       "purrr","xgboost","lattice","tictoc","scico")

purrr::walk(required.packages,library, character.only = T)

data = mlr3::as_data_backend(ecuador)
task = TaskClassifST$new("ecuador",
                         backend = data, target = "slides",
                         positive = "TRUE", extra_args = list(coordinate_names = c("x", "y"))
)

learner = mlr3::lrn("classif.xgboost", predict_type = "prob")
learner$fallback = lrn("classif.xgboost", predict_type = "prob")
learner$param_set$params$nthread

#set to use 4 CPUs
set_threads(learner, n = availableCores())

terminator = mlr3tuning::trm("evals", n_evals = nmod)
tuner = mlr3tuning::tnr("random_search")
perf.level = mlr3::rsmp("repeated_spcv_coords", folds = nfolds, repeats = nrep)
tune.level = mlr3::rsmp("spcv_coords", folds = nfolds)

##Parallization didn't work, mlr3 has a default nthread of 1 to make it 
#consistent with future package. So instead set nthread to 16 that will
#speed up the model training

search_space = paradox::ps(
  nrounds = paradox::p_int(lower = 100, upper = 300),
  max_depth = paradox::p_int(lower = 1, upper = 6),
  eta = paradox::p_dbl(lower = .1, upper = .5),
  min_child_weight = paradox::p_dbl(lower = 1, upper =10),
  colsample_bytree = paradox::p_dbl(lower = 0.5, upper =1),
  subsample = paradox::p_dbl(lower = 0.5, upper =1),
  lambda = paradox::p_dbl(lower = -1, upper = 0, trafo = function(x) 10^x))

autotune_xgboost <- mlr3tuning::AutoTuner$new(
  learner = learner,
  resampling = tune.level,
  measure = mlr3::msr("classif.auc"),
  search_space = search_space,
  terminator = terminator,
  tuner = tuner
)

#############Re-sampling #############

################################
tic()
resampled <- mlr3::resample(task = task,
                                #Inner resampling tuning level within this object
                                learner = autotune_xgboost,
                                # outer resampling (performance level)
                                resampling = perf.level,
                                store_models = TRUE,
                                encapsulate = "evaluate")

toc()
autoplot(resampled)

# #Extract the learners for future predictions
models <- mlr3misc::map(as.data.table(resampled)$learner, "learner")

# 
# #Get the learner
uncertainty <- models %>%
                map(pluck,"predict")

models <- mlr3misc::map(as.data.table(resampled)$learner, "learner")

l <- models %>% pluck(1)

m <- l$predict_newdata(df, task)$prob[,1] 

#####
myfunct <- function(x) {
  learners <- pluck(x)
  predict <- learners$predict_newdata(df, task)$prob[,1]
}

#Get the probability in a list
a <- models %>%
       map(myfunct)

#Convert to a dataframe
b <- t(sapply(a, c)) %>%
        as.data.frame()

p <- b %>%
    map_dfr(quantile,c(0.05, 0.95)) %>%
    set_names("low","upp")

p <- p %>%
  mutate(d = upp-low)

# 
# ##############################
# 
# # compute the AUROC as a data.table
auc.score <- resampled$score(
               measure = mlr3::msr("classif.auc")) %>%
               as.data.frame() %>%
               select(task_id, learner_id, resampling_id, classif.auc)
# 
# hist(auc.score$classif.auc)
# # 
mean(auc.score$classif.auc)
# 
# ##### Autotune using all of the data
autotuned = mlr3tuning::AutoTuner$new(
  learner = learner,
  resampling = mlr3::rsmp("spcv_coords", folds = nfolds), # spatial partitioning
  measure = mlr3::msr("classif.auc"), # performance measure
  search_space = search_space, # predefined hyper-parameter search space
  terminator = mlr3tuning::trm("evals", n_evals = nmod), # specify the number of iterations
  tuner = mlr3tuning::tnr("random_search") # specify random search
)
# 
# #Auto-tune a model
tic()
tuned.model = autotuned$train(task)
toc()
# 
m <- tuned.model$predict_newdata(df, task)
n <- m$prob

# # p <- pred$data
# # q <- do.call(cbind.data.frame, p)        

