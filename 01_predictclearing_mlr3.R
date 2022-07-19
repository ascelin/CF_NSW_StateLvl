##July 2022 using mlr3 modelling approach##

#Remove all data from the environment
rm(list = ls(all.names = TRUE))

#Project defaults: Albers Equal Projection
project.crs <- 'EPSG:3577'
project.res <- 100

#Load required libraries
required.packages <- c("sf","terra","data.table","tidyverse",
                       "mlr3","mlr3learners","mlr3spatiotempcv",
                       "mlr3viz","mlr3tuning","iml","future","furrr",
                       "purrr","xgboost","lattice","tictoc","scico","ggtext")

purrr::walk(required.packages,library, character.only = T)

#####User modified parameters
region <- c("state","Central West","Tablelands","Coastal","Western")
agent <- c("agri","fores","infra","combined")

#For testing
# region <- "state"
# agent <- "agri"

######### Modelling parameters #################

#n-samples <- 100 #Samples # all samples are selected in the code below use this later if there is a need
nfolds <- 5 #CV folds
nreps <- 2 #Number of times to repeat CV
nmod <- 50 #Hyper parameter search limit

##################################3 DON"T MODIFY ANYTHING BELOW THIS CODE ##########################

#Create a list of study area and bind them for loop

#Specify the data path based on the system
data.path <- case_when(
  Sys.info()["sysname"] == "Windows" ~ "./data/",
  #Sys.info()["sysname"] == "Darwin" ~ "mac",
  Sys.info()["sysname"] == "Linux" ~ "/home/ubuntu/data/"
)

#Warning message if it can't find the directory
if (dir.exists(data.path)){
  print("Directory exists - process will run")
} else {
  string.to.print <- paste("ERROR: can't find the directory: check data path", data.path)
  stop(string.to.print)
}

nsw <- st_read(str_c(data.path,"studyarea/state/NSW_STATE_POLYGON_shp_ex_islands_proj.shp"))
combined_bio <- st_read(str_c(data.path,"studyarea/Cfact_analysis_regions/Cfact_analysis_regions_prj.shp"))
bioregion <- st_read(str_c(data.path,"studyarea/ibra/IBRA_NSW_clipped.shp"))

studyarea <- rbind(nsw %>% transmute(name = "state"),
                   combined_bio %>% transmute(name = combined_bio$Cfact_Regi), 
                   bioregion %>% transmute(name = bioregion$REG_NAME_7))

###### SECTION 1: DATA PREPARATION ##############

# region <- "Coastal"
# agent <- "agri"

do_analysis <- function(region, agent) {

  # #Test
  # region <- "Coastal"
  # agent <- "agri"
  
  # Start the timer
  tic("Start analysis")
  
  results.path <- case_when(
    Sys.info()["sysname"] == "Windows" ~ str_c("./results/",region,"/"),
    #Sys.info()["sysname"] == "Darwin" ~ "mac",
    Sys.info()["sysname"] == "Linux" ~ str_c("/home/ubuntu/results/",region,"/")
  )

 if (!dir.exists(results.path)){dir.create(results.path)}

 print(str_c("Preparing data to run ", agent, " in ", region, " bioregion"))
  
 loss.path <- dir(str_c(data.path, "loss"), full.names=T, pattern = "majorityrule")
 cov.path <- dir(str_c(data.path, "covariates", sep=""), full.names = T, pattern = ".tif$")

 #Remove minimum lot size from the state level as there is no data in western states
 cov.path <- str_subset(cov.path, pattern = "minimumlotsize", negate = T)

 #Load BCT and WOODY DATA
 bct <- st_read(dir(str_c(data.path,"/bct/"), full.names = T, pattern = "*.shp$"))
 woody.private.land <- rast(dir(str_c(data.path, "woodyonprivateland", sep=""), full.names = T, pattern = ".tif$"))

 #Create a function to run analysis across the parameters
 roi <- studyarea %>% filter(name == region)
 
 ## 1.1 Get loss raster
 loss.file <- case_when(
  agent == "agri" ~ str_subset(loss.path,"agri"),
  agent == "fores" ~ str_subset(loss.path, "fores"),
  agent == "infra" ~ str_subset(loss.path,"infra"),
  agent == "afcombined" ~ c(str_subset(loss.path,"agri"),
                            str_subset(loss.path,"fores")))

 loss.file <- unique(loss.file) #Remove duplicate filenames

 loss.raster <- rast(loss.file)

 #Crop to the ROI only if it is not a state model
 ifelse(
   region == "state",
   loss.raster,
   loss.raster <- loss.raster %>% crop(roi) %>% mask(vect(roi))
 )
 
#loss.raster[loss.raster!=1] <- NA #Turn everything except 1 to NA
#This is not required in the new GEE state files - but if others then check this

#If there is agri and forestry then merge them into 1 file for combined loss
#if else(condition,true,false)

ifelse(
  nlyr(loss.raster) > 1, 
  loss.raster <- merge(loss.raster[[1]],loss.raster[[2]]),
  loss.raster
)

# 1.2 Mask out the losses outside of the woody in the private land

#Clip Loss Raster
all.losses <- rast(str_c(data.path,"loss/", "nsw_state_allagents_resampled100m_majorityrule.tif")) 

ifelse(
  region == "state",
  all.losses,
  all.losses <- all.losses %>% crop(roi) %>% mask(vect(roi))
)

woodyonprivate <- rast(str_c(data.path, "woodyonprivateland/","woodyonprivate.tif"))

ifelse(
  region == "state",
  woodyonprivate,
  woodyonprivate <- woodyonprivate %>% crop(roi) %>% mask(vect(roi))
)

##1.3 Now get the training points/pixels from the loss raster
# loss.raster <- loss.raster %>% 
#   crop(roi) %>%
#   mask(vect(roi))

loss.pts <- as.points(loss.raster) %>%
  geom() %>%
  as.data.frame() %>%
  transmute(x,y)

#Set n_samples to all loss points - change this if there is computational issue
nsamples <- as.numeric(nrow(loss.pts))

## 2. Collect background samples
bg.raster <- woodyonprivate %>%
  mask(., vect(bct), inverse = T) %>%
  mask(., all.losses, inverse = T)

###2.4 Samples to take from loss data
nsamples <- case_when(
  nrow(loss.pts) < nsamples ~ nrow(loss.pts),
  TRUE ~ as.integer(nsamples)
)

nsamples <- nsamples

loss.pts <- loss.pts %>%
  slice_sample(n=nsamples)

# 2.4 Take sample from the background same number
# as the loss points # Balanced data_set #

#Terra sample selected less points - change in code behavior so
# sample using tidyverse
#It this code had worked then we didn't have to mask and this would have
#saved some time - spatSample restrict to extent works but not to polygon
# check again at a later date if possible
# test.pts <- spatSample(bg.raster, nrow(loss.pts),
#                        method="random",
#                      #xy=T,
#                      as.points=T,
#                      ext = ext(roi),
#                      na.rm = T)

set.seed = 2022

bg.pts <- as.points(bg.raster) %>%
  geom() %>%
  as.data.frame() %>%
  transmute(x,y) %>%
  slice_sample(n = nrow(loss.pts))

## 3. Join loss and background points 1= loss and 0 = background
pts <- rbind(
  loss.pts %>% mutate(loss = as.factor(1)),
  bg.pts %>% mutate(loss = as.factor(0)))

##3.1 Convert to sf object
pts <- st_as_sf(pts, coords= c("x","y"))

model.name <- str_c("roi_",region,"_",
                    "agent_",agent,"_",
                    "samples_",nrow(loss.pts),"_",
                    "cv_",nfolds,"_",
                    "rep_",nreps,"_",
                    "mod_",nmod)

## 4.4 Covariates data

covariates <- rast(cov.path) 

ifelse(
  region == "state",
  covariates,
  covariates <- covariates %>% crop(roi) %>% mask(vect(roi))
)

print(str_c("No of covariates used for modelling: ", nlyr(covariates)))

#extract names from file_names

good.names <- str_extract(cov.path,pattern = "(?<=_)[^.]*(?=.)")
good.names
#Assign names now
names(covariates) <- good.names

## 5 Extract the covariate data into the training samples

df <- terra::extract(covariates,vect(pts), xy=T) %>%
  data.frame() %>%
  mutate(loss = as.factor(pts$loss)) %>%
  dplyr::select(-ID)%>% #At this step check which points have NA -
  drop_na() #mlr doesn't take na so remove if any column has NA

unique(df$forestryapprovals)

print(str_c("Data preparation complete for model:: ",model.name))

# #Plot the distribution of data
# df_long <- df %>% pivot_longer(
#   cols = - c(loss),
#   names_to = "variables") %>%
#     ggplot(aes(x = loss, y = value))+
#     geom_boxplot()+
#     facet_wrap(~variables, scales = "free_y")
#
# df_long
#
# boxplot_fun <- function(x,y){
#   ggplot(df, aes_string(x = x, y = y))+
#     geom_boxplot()
# }
#
# boxplot_fun (x = "elevation", y="loss")

###Get the coordinates for spatial partitioning
#coords <- df[ , c("lon","lat")]

################# STEP 2: MODELLING SECTION ###########################
print(str_c("Model training starting for ", agent, " in ", region, " bioregion"))

## Step 2.1 Create a classification task

task = mlr3spatiotempcv::TaskClassifST$new(
  id = "xgboost_model",
  backend = mlr3::as_data_backend(df),
  target = "loss",
  positive = "1",
  extra_args = list(
    coordinate_names = c("x", "y"), #Specify to use these columns as coordinates
    coords_as_features = FALSE,
    crs = "EPSG:3577") #Albers
)

levs = levels(task$truth())

## Step 2.2 Choose a learner and set predict output in probability
learner = mlr3::lrn("classif.xgboost",predict_type = "prob")

#booster = "gbtree",
#tree_method = "hist")

##Supposedly, specifying booster=gbtree and tree_method = "hist" in the learner was going to run xgboost faster
#tried on my laptop but it made the run slower, perhaps with GPU you can use "gpu_hist" sometime in the future

#To make sure the tuning doesn't stop if any model fails
learner$fallback = lrn("classif.xgboost", predict_type = "prob")

#set to use 4 CPUs
set_threads(learner, n = availableCores())

#Check the parameters you can set
# learner$param_set$ids()
# learner$help()

##Step 2.2 Performance estimation level using resampling technique

#Nested Spatial CV

#Create a outer level resampling strategy
perf.level = mlr3::rsmp("repeated_spcv_coords", folds = nfolds, repeats = nreps)

#Create an inner level resampling strategy
tune.level = mlr3::rsmp("spcv_coords", folds = nfolds)

cvplot = autoplot(tune.level, task, c(1, 2, 3, 4, 5),
                  crs = 3577, point_size = 3, axis_label_fontsize = 10,
                  plot3D = TRUE
)

#cvplot

ggsave(str_c(results.path,model.name,"_spatialcv.png"), cvplot, bg="white")

##specify the budget available for tuning - we use
##terminate after a given number of iterations of hyper-parameter combinations
evals = mlr3tuning::trm("evals", n_evals = nmod)

#Set the optimization of algorithm takes place - we choose random search
tuner = mlr3tuning::tnr("random_search")

# define the outer limits of the randomly selected hyper-parameters
search_space = paradox::ps(
  # The number of trees in the model (each one built sequentially)
  nrounds = paradox::p_int(lower = 100, upper = 300),
  # number of splits in each tree
  max_depth = paradox::p_int(lower = 1, upper = 6),
  # Learning rate - "shrinkage" - prevents overfitting
  eta = paradox::p_dbl(lower = .1, upper = .4),
  min_child_weight = paradox::p_dbl(lower = 1, upper =10),
  colsample_bytree = paradox::p_dbl(lower = 0.5, upper =1),
  subsample = paradox::p_dbl(lower = 0.5, upper =1),
  lambda = paradox::p_dbl(lower = -1, upper = 0, trafo = function(x) 10^x))

#Tuner for re-sampling only
at_resample <- mlr3tuning::AutoTuner$new(
  learner = learner,
  resampling = tune.level,
  measure = mlr3::msr("classif.auc"),
  search_space = search_space,
  terminator = evals,
  tuner = tuner
)

#Don't run the resampling for now##
#######################################################################
#######################################################################
#
# future::plan(list("sequential", "multisession"),
#              workers = floor(availableCores() - 2))


#
# progressr::with_progress(expr = {

resampled = mlr3::resample(task = task,
                           learner = at_resample,
                           # outer resampling (performance level)
                           resampling = perf.level,
                           store_models = F,
                           encapsulate = "evaluate")

# saveRDS(resampled, str_c(results.path,"resampled.Rds"))

# future:::ClusterRegistry("stop")
# # compute the AUROC as a data.table
auc <- resampled$score(
  measure = mlr3::msr("classif.auc")) %>%
  as.data.frame()

auc.mean <- mean(auc$classif.auc, na.rm=T)

auc_plot <- auc %>%
  ggplot() +
  geom_histogram(aes(x = classif.auc, y = (..count..)/sum(..count..)),
                 bins = 20, fill = "turquoise4", colour = "gray") +
  geom_vline(
    data = auc %>% summarise(m = mean(auc$classif.auc)),
    aes(xintercept = m), linetype="dotted", col = "black", size = 2)+
  theme_minimal() +
  xlab("AUC values") +
  ylab("Relative frequency")+
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = -0.5))

ggsave(str_c(results.path,model.name,"_auc.png"), auc_plot, bg="white")

#Get the models from the resampled data
#boot.models <- mlr3misc::map(as.data.table(resampled)$learner, "learner")

#Create a function to predict to a new dataframe
# boot.predict <- function(x) {
#   learners <- pluck(x)
#   predict <- learners$predict_newdata(new_df, task)$prob[,1]
# }

#########################################################################
##########################################################################

#Note;You should report this evaluation metric rather than
# from the final tuned model that we built next

#Until here is the model assessment, we don't use the hyper-parameters
# from these models, now we tune the hyper-parameters using
# all the data # We don't need the inner re sampling here
# With hyper-parameters combinations there will be 5 * 50 = 250 models here
# 
at_main = mlr3tuning::AutoTuner$new(
  learner = learner,
  resampling = mlr3::rsmp("spcv_coords", folds = nfolds), # spatial partitioning
  measure = mlr3::msr("classif.auc"), # performance measure
  search_space = search_space, # predefined hyper-parameter search space
  store_models = T,
  terminator = mlr3tuning::trm("evals", n_evals = nmod), # specify 50 iterations
  tuner = mlr3tuning::tnr("random_search") # specify random search
)
# 
# # hyper-parameter tuning
set.seed(2022)

#Auto-tune a model
tuned.model = at_main$train(task)
# 
# #VIP plot
vip <- tuned.model$learner$importance() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  set_names("variable","importance") %>%
  ggplot(aes(reorder(variable,importance),importance))+
  geom_bar(stat = "identity")+
  labs(y = "Variable importance", x = "")+
  coord_flip()+
  theme_minimal()
# 
# vip
# 
ggsave(str_c(results.path, model.name,"_vip.png"), vip, bg="white")

#Don't run model interpretation now
#################### Model interpretation ##############################
#######################################################################

#get a model object to use for iml
# model <- Predictor$new(tuned.model$learner, data = df %>% select(-loss), y = df$loss)
#
# effect <- Feature Imp$new(model, loss = "ce")
#
# #Get the feature effects
# num.features <- c("woody frag","slope","minimum lot size")
# feature.effect <- Feature Effects $new(model)
# plot(feature.effect, features = num.features)

######################################################################
######################################################################

#terra predict didn't have an option to predict to prob

# used a custom function from
# pfmlr = function(model, ...) {
#   if(model$predict_type == "prob") {
#     p = model$predict_newdata(...)$data$prob
#     if(length(levs) != ncol(p)) {
#       missing = setdiff(levs, colnames(p))
#       pm = matrix(0, ncol = length(missing), nrow = nrow(p), dimnames = list(NULL, missing))
#       p = cbind(p, pm)
#       p = p[, levs]
#     }
#     p[,1] #get the prob for 1st class
#   } else {
#     model$predict_newdata(...)$data$response
#   }
# }
# 
# #Make prediction
# pred <- terra::predict(covariates, tuned.model, fun = pfmlr, na.rm = TRUE)
# pred <- round(pred, 5)
# writeRaster(pred, str_c(results.path,model.name, "_predrisk",".tif"),overwrite=T)

toc(log = TRUE, quiet = TRUE)
log.txt <- unlist(tic.log(format = T))
tic.clearlog()

#Save all model details in a dataframe
model_details <- data.frame(
  roi = region,
  samples = nsamples,
  totalloss_pixels = nrow(loss.pts),
  cv = nfolds,
  mods = nmod,
  rep = nreps,
  agent = agent,
  auc_mean = auc.mean,
  auc_dist = I(list(auc$classif.auc)),
  time = log.txt)

write.csv(model_details, str_c(results.path, model.name,"_modeltime.csv"))

gc()

}

#Run the function over the list

## Get a list of all combinations to arrange by agent type

#Only state; 4 combined regions, and 2 bioregions

df.list <- crossing(studyarea$name,agent) %>%
            set_names("region", "agent") %>%
            filter(region %in% c("state","NSW North Coast","NSW South Western Slopes",
                                 combined_bio$Cfact_Regi)) %>%
            arrange(agent)

purrr::pwalk(list(
  region = df.list$region,
  agent = df.list$agent),
  .f = do_analysis)

#New predict data
# newdata = as.data.frame(as.matrix(covariates))
# ind = rowSums(is.na(newdata)) == 0
# 
# models <- mlr3misc::map(as.data.table(resampled)$learner, "learner")
# 
# #Create a function to predict to new data
# myfunct <- function(x) {
#   learners <- pluck(x)
#   predict <- learners$predict_newdata(newdata[ind, ], task)$prob[,1]
# }

#Get the probability in a list
# plan("multisession", workers = availableCores())
# options(future.globals.maxSize= 8912896000)

# pred_int <- models %>%
#   map(myfunct) %>% #apply the function
#   do.call(rbind, .) %>% #Bind the rows
#   data.frame() 
# 
# #Get the quarantines
# pred_int <- pred_int %>%
#   map_dfr(quantile,c(0.05, 0.95)) %>% #Get the percentiles
#   set_names("lower","upper")
# 
# future:::ClusterRegistry("stop")
# 
# newdata[ind, "lower"] <- pred_int$lower
# newdata[ind, "upper"] <- pred_int$upper
#Convert to a dataframe
# b <- t(sapply(a, c)) %>%
#   as.data.frame()
# 
# p <- b %>%
#   map_dfr(quantile,c(0.05, 0.95)) %>%
#   set_names("low","upp")
# 
# q <- a %>%
#   map_dfr(quantile,c(0.05, 0.95)) %>%
#   set_names("low","upp")
# 
# t <- do.call(bind_rows, a)

# pred_low = covariates$lscrocky
# pred_low[] = newdata$lower
# names(pred_low) = pred_low
# writeRaster(pred_low, str_c(results.path,"_","pred_lower",".tif"),overwrite=T)
# 
# pred_upper = covariates$lscrocky
# pred_upper[] = newdata$upper
# names(pred_upper) = "pred_upper"
# writeRaster(pred_upper, str_c(results.path,"_","pred_upper",".tif"),overwrite=T)
