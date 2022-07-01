
#Remove data from the environmentkkk
rm(list = ls(all.names = TRUE))

#Project defaults: Albers Equal Projection
project.crs <- 'EPSG:3577'
project.res <- 100

#Load required libraries

#Spatial
library(sf)
library(terra)

#library(raster)
library(data.table)

#Others
library(tidyverse)
library(mlr3)
library(mlr3learners)
library(mlr3spatiotempcv)
library(mlr3viz)
library(mlr3tuning)

#Model interpreation
library(iml)

#Parallel processing
library(future)
library(furrr)

#Data 
library(purrr)
library(xgboost)
library(lattice)
library(tictoc)
library(scico)

## User modified parameters
##----------------------------------------------------#
#Select which region to model

check <- tsk("ecuador")

## The region can be state; combined bio-region or single bio-region

#Select region
#regions.available <- c("state","Central West","Tablelands","Coastal","Western")
region <- "Coastal"

#agent.available <- c("agri","fores","infra","combined")
agent = "agri"

######### Modelling parameters #################

#How many samples to take, keep increasing until there is no accuracy gain

#Limit samples (no of loss and background points) for modelling
#Now model by increasing the number of samples - increase the sample to see
#if there will be increase in the accuracy

nsamples <- 5000
#No of cross-validation iterations and no.of models selection

nfolds = 5
nreps = 2
#Number of combinations from random grid search 
# for different hyper-parameter combination
nmod = 2

# Only relevant to windows and mac for now- parallel will be off for nectar
runparallel=F

##################################3 DON"T MODIFY ANYTHING BELOW THIS CODE ##########################

#Load datasets
nsw <- st_read("./data/studyarea/state/NSW_STATE_POLYGON_shp_ex_islands_proj.shp")
bioregion <- st_read("./data/studyarea/Cfact_analysis_regions/Cfact_analysis_regions_prj.shp")

ifelse(
  region == "Coastal",
  roi <- bioregion %>% filter(Cfact_Regi == region),
  roi <- nsw
)

print(str_c("Model starting for ", region))

#Specify the data path based on the system
system <- case_when(
  Sys.info()["sysname"] == "Windows" ~ "windows",
  Sys.info()["sysname"] == "Darwin" ~ "mac",
  Sys.info()["sysname"] == "Linux" ~ "nectar"
)

###### SECTION 1: DATA PREPARATION ##############

#Create a function that iterates over region, agent, and loss.pixel.type

nectar <- "/home/ubuntu/statelevelanalysis"

 
  data.path <- case_when(
      system == "windows" ~ str_c("./data/"),
      system == "nectar" ~ str_c(nectar,"/data/")
    )
  
  bct.path <- case_when(
    system == "windows" ~ str_c("./data/bct/"),
    system == "nectar" ~ str_c(nectar,"/data/bct/")
    )
  
  if (dir.exists(data.path)){
    print("Directory exists - process will run")
  } else {
    string.to.print <- paste("ERROR: can't find the directory: check data path", data.path)
    stop(string.to.print)
  }
  
  #Specify the export path
  results.path <- 
    case_when(
      system == "windows" ~ str_c("./results/",region,"/"),
      system == "nectar" ~ str_c(nectar,"/results/",region,"/")
    )
  
  if (!dir.exists(results.path)){dir.create(results.path)}

  loss.path <- dir(str_c(data.path, "loss"), full.names=T, pattern = "majorityrule")
  cov.path <- dir(str_c(data.path, "covariates", sep=""), full.names = T, pattern = ".tif$")
  
  #Load other datasets
  bct <- st_read(dir(bct.path, full.names = T, pattern = "*.shp$"))
  woody.private.land <- rast(dir(str_c(data.path, "woodyonprivateland", sep=""), full.names = T, pattern = ".tif$"))
  
    ## 1.1 Get loss raster
  loss.file <- case_when(
    agent == "agri" ~ str_subset(loss.path,"agri"),
    agent == "fores" ~ str_subset(loss.path, "fores"),
    agent == "infra" ~ str_subset(loss.path,"infra"),
    agent == "afcombined" ~ c(str_subset(loss.path,"agri"),
                              str_subset(loss.path,"fores")))
  
  loss.file <- unique(loss.file) #Remove duplicate filenames
  
  loss.raster <- rast(loss.file)
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
  #woody.on.private[woody.on.private==0] <- NA
  #This code is not required in the new import from GEE
  
  #Clip Loss Raster
  # loss.raster <- loss.raster %>% 
    # mask(.,woody.private.land) %>%
    # mask(., vect(roi))
  
  ##1.3 Now get the training points/pixels from the loss raster
  loss.pts <- as.points(loss.raster) %>%
               geom() %>%
               as.data.frame() %>%
               transmute(x,y)
  
  ### 2. Collect background samples
  
  ## Mask out all losses and BCT agreement sites from the background
  
  # 2.1 Load mask elements
  
  #All losses
  all.loss <- rast(str_subset(loss.path,"allagent"))
  #all.loss[all.loss==0] <- NA #remove the 0's 
  
  # 2.3 Get the background raster now by masking out
  bg.raster <- woody.private.land #%>%
    # mask(., vect(bct), inverse = T) %>%
    # mask(., all.loss, inverse = T) %>%
    # mask(., vect(roi)) #Clip to the roi boundary
  
  ##2.4 Samples to take from loss data
  no.to.sample <- case_when(
    nrow(loss.pts) < nsamples ~ nrow(loss.pts),
    TRUE ~ as.integer(nsamples)
  )
  
  loss.pts <- loss.pts %>%
    slice_sample(n=no.to.sample)
  
  used.samples <- nrow(loss.pts)
  
  # 2.4 Take sample from the background same number 
  # as the loss points # Balanced dataset #
  
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
  
  bg.pts <- as.points(bg.raster) %>%
    geom() %>%
    as.data.frame() %>%
    transmute(x,y) %>%
    slice_sample(n = nrow(loss.pts))
  
  ## 3. Join loss and background points
  pts <- rbind(
          loss.pts %>% mutate(loss = as.factor(1)), 
          bg.pts %>% mutate(loss = as.factor(0)))
  
  ##3.1 Convert to sf object
  pts <- st_as_sf(pts, coords= c("x","y"))
  
  #st_write(pts, "pts_test.shp")
  
  print(c("size of total points (loss and intact):", nrow(pts)))
  
  model.name <- str_c("roi_",region,"_",
                      "agent_",agent,"_",
                      "samples_",nrow(pts),"_",
                      "cv_",nfolds,"_",
                      "rep_",nreps)
  
  print(str_c("Starting model:: ",model.name))
  
  ## 4.4 Covariates data
  
  covariates <- rast(cov.path)
  
  #extract names from filenames
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
  
  
  #sf caused an issue so unload that before running mlr
  
  # devtools::unload("sf")
  # devtools::unload("terra")
  # devtools::unload("tidyverse")
  # 
  # devtools::unload("mlr3spatiotempcv")
  # library(mlr3spatiotempcv)
  
  #data_sf = sf::st_as_sf(df, coords = c("x", "y"))
  
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
  
  ## Step 2.1 Create a classification task
  
  # task = mlr3spatiotempcv::TaskClassifST$new(
  #   id = "xgboost_model",
  #   backend = df_sf,
  #   target = "loss",
  #   positive = "1")
  
  task = mlr3spatiotempcv::TaskClassifST$new(
    id = "xgboost_model",
    backend = mlr3::as_data_backend(df),
    target = "loss",
    positive = "1",
    extra_args = list(
      coordinate_names = c("x", "y"), #Specify to use these columns as coordinates
      coords_as_features = FALSE,
      crs = "EPSG:3577")
  )




