rm(list = ls())
setwd("C:/Users/Admin/Desktop/IZ/conefor_R/note_importance/result/t2t3_lulc_0.4rm") # Set working directory
#install.packages("usdm")
#devtools::install_github("rpatin/gbm")
#devtools::install_github("mrmaxent/maxnet")
#library(devtools)
#devtools::install_github("HemingNM/ENMwizard")
#devtools::install_github("biomodhub/biomod2", dependencies = TRUE) # Install the new biomod2 package
library(spThin)
library(car)
library(raster)
library(terra)
library(biomod2)
library(ENMwizard)
library(tidyverse)
library(usdm)

#folder_path <- "D:/7_Paper/hyh/var/c_v/a1999"
# Read all raster data in the folder
#rasters <- list.files(path = folder_path, pattern = ".tif", full.names = TRUE)
# Stack rasters
#raster_stack <- stack(rasters)
# Use the select_vars() function in the ENMwizard package to calculate correlations and select variables
#result <- ENMwizard::select_vars(
#  env = raster_stack,
#  cutoff = 0.7,
#  corr.mat = NULL,
#  names.only = T,
#  plot.dend = T,
#  rm.old = F,
#  sp.nm = "Trachypithecus francoisi")

# Export the calculated correlation matrix
#write.csv(result$corr.mat, file = "correlation_matrix.csv", row.names = TRUE)

# Extract selected variables
#selected_vars <- result$sel_vars
#myExpl <- raster_stack[[selected_vars]] # Stack raster files
#plot(myExpl)# Plotting


# Calculate multicollinearity (Variance Inflation Factor for variables)#### Use vifor function in usdm package
#myExpl<-as.data.frame(myExpl)
#vif(myExpl) 
#v1 <- vifcor(myExpl, th=0.9);v1   
#re1 <- exclude(myExpl,v1);re1    
#v2 <- vifstep(myExpl, th=5);v2  
#re2 <- exclude(myExpl, v2);re2    

##########################################################################################################################
# Filter species occurrence data
#occ <- read.csv("D:/7_Paper/hyh/occ/Point1987_1999.csv", header = TRUE)# Read species occurrence data
#head(occ)
#thin_occ<-thin(
#  loc.data=occ, ## Original species occurrence dataset
#  lat.col = "latitude",
#  long.col = "longitude",
#  spec.col = "sps",
#  thin.par=1,  # Thinning distance, unit: kilometers
#  reps=5,   # Number of thinning iterations
#  locs.thinned.list.return = FALSE,
#  write.files = TRUE,    # Whether to save csv files
#  max.files = 1,      # Number of output files
#  out.dir="C:\\Users\\10284\\Desktop\\QXN\\result", # Output file directory
#  out.base = "thinned_occ",     # Output file name
#  verbose = TRUE)
# Extract thinned species data
#thined_occ <- read.csv("C:\\Users\\10284\\Desktop\\QXN\\occ\\mh2km.csv")
#pre_occ<-thined_occ[,-1] # Remove the first column
#pre_occ# View species occurrence data
#plot(pre_occ)
# length(which(!is.na(values(subset(pre_var, 1))))) Calculate the number of non-missing values in a raster file
#pre_var# View predictor variables


######################## Moran's I autocorrelation test
#library(spdep)
# Read data from CSV file
#data <- read_csv("D:/7_Paper/hyh/occ/Point1987_1999.csv")

# Load necessary packages
#x<-data$species
#x
# Assume data contains 75 processed points with longitude and latitude in 'longitude' and 'latitude' columns
#coords <- cbind(data$longitude, data$latitude)

# Create adjacency matrix, here using k-nearest neighbors (k=4), can be adjusted according to actual data
#nb <- knn2nb(knearneigh(coords, k=4))

# Convert to adjacency weight matrix
#listw <- nb2listw(nb, style="W")

# Calculate Moran's I
#moran.test(data$species, listw)
#————————————————————————————————————————————————————————————————————————————————————————————————————————————————————
#——————————————————————————————————————————————————————————————————————————————————————————————————————————————
# Load species occurrence data
DataSpecies <- read.csv("C:/Users/Admin/Desktop/IZ/conefor_R/note_importance/occ/Point2000_2024_f.csv")
head(DataSpecies)
myRespName <- 'species' # Column name of the study species
myResp <- as.numeric(DataSpecies[, myRespName]) # Convert species column to continuous data 
N <- sum(myResp)
N

myRespXY <- DataSpecies[,c("longitude","latitude")]; head(myRespXY)
## Load environmental data
#folder_path <- "D:/7_Paper/hyh/var/study_var/lulc2019"
folder_path <- "C:/Users/Admin/Desktop/IZ/conefor_R/note_importance/tif/note_imp_t2t3_lulc_0.4rm.tif"
#rasters <- list.files(path = folder_path, pattern = ".tif", full.names = TRUE)
# Stack rasters
rasters <- raster(folder_path)
myExpl <- stack(rasters)
plot(myExpl )

### Convert data to the format required by the package and obtain pseudo-absence data
set.seed(1) # Set random seed
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, 
                                     expl.var = myExpl, 
                                     resp.xy = myRespXY,   
                                     resp.name = myRespName, 
                                     PA.nb.rep = 1, # Number of repetitions
                                     filter.raster = TRUE,
                                     PA.nb.absences = length(myResp),#10000,## Number of pseudo-absences
                                     PA.dist.min = 1000,                # Minimum distance from pseudo-absences to occurrence points
                                  #  PA.dist.max = 300,                # Maximum distance from pseudo-absences to occurrence points
                                     PA.strategy = 'random')
plot(myBiomodData)

################################
### Extract pseudo-absence data as a table
#patable <- NULL
#patable$species <- myBiomodData@data.species
#patable$coord <- myBiomodData@coord
#patable$env <- myBiomodData@data.env.var
#patable <- as.data.frame(patable)

#write.csv(patable,"D:/7_Paper/hyh/result/lulc2019/patableuva.csv")
###############################################

#--------------------------------------------------------------------------------- Automatic model parameter configuration
# k-fold selection
#cv.k <- bm_CrossValidation(bm.format = myBiomodData,
#                            strategy = "kfold",
#                            nb.rep = 5,  # Number of repetitions
#                            k =5)      # Number of cross-validation folds
#head(cv.k)

#allModels <- c("ANN","FDA","CTA","GBM", "MARS", "RF","MAXENT.Phillips.2","GAM")#"MARS","ANN","CTA", "FDA", "GBM", "GLM", 

# bigboss parameters
#myBiomodOption <- bm_ModelingOptions(data.type = 'binary',
#                             models = allModels,
#                             strategy = "bigboss",
#                            bm.format = myBiomodData
#                            ,calib.lines = cv.k
#                           )
# 
#
#myBiomodOption


###———————————————————————————————————————————————————————————— Custom model settings for each model

#user.ANN <- list('_allData_allRun' = list(size = 2),
#                 '_PA1_allRun' = list(size = 2))
#               #  '_PA1_RUN1' = list(size = 2),
#                 #'_PA1_RUN2' = list(size = 2))
#
#user.CTA <- list('_allData_allRun' = list(method = "class"),
#                 '_PA1_allRun' = list(method = "class"),
#                 '_PA1_RUN1' = list(method = "class"))
#                 #'_PA1_RUN2' = list(method = "class"))
#
#user.MARS <- list('_allData_allRun' = list(pmethod = "backward"),
#                  '_PA1_allRun' = list(pmethod = "backward"))
##
#user.GLM <- list('_allData_allRun' = list(type = 'quadratic', interaction.level = 0),
#                 '_PA1_allRun' = list(type = 'quadratic', interaction.level = 0))
#
#user.GBM <- list('_allData_allRun' = list(n.trees = 2000),
#                 '_PA1_allRun' = list(n.trees = 2000))
#
#user.GAM <- list('_allData_allRun' = list(algo = "GAM.gam.gam"))
#
#user.RF <- list('_allData_allRun' = list(ntree = 1000),
#                '_PA1_allRun' = list(ntree = 1000))
#
#user.MAXENT <- list('_allData_allRun' = list(path_to_maxent.jar = getwd(),
#                                             product = TRUE,
#                                             threshold = TRUE,
#                                             hinge = TRUE),
#                    '_PA1_allRun' = list(path_to_maxent.jar = getwd(),
#                                         product = TRUE,
#                                         threshold = TRUE,
#                                         hinge = TRUE)
#)
#

#
#user.val <- list(
#                  CTA.binary.rpart.rpart = user.CTA,
#                  ANN.binary.nnet.nnet = user.ANN,
#                  FDA.binary.mda.mda = user.FDA,
#                  GBM.binary.gbm.gbm = user.GBM,
#                  GLM.binary.stats.glm = user.GLM ,
#                  MARS.binary.earth.earth = user.MARS,
#                  MAXENT.binary.MAXENT.MAXENT = user.MAXENT,
#                  RF.binary.randomForest.randomForest = user.RF
#                  # GAM.binary.mgcv.gam = user.GAM
#                   )

#user.val <- list(
# ANN.binary.nnet.nnet = user.ANN,
# CTA.binary.rpart.rpart = user.CTA
# FDA.binary.mda.fda = user.FDA,
# GBM.binary.gbm.gbm = user.GBM,
# GLM.binary.stats.glm = user.GLM,
# MARS.binary.earth.earth = user.MARS,
# MAXENT.binary.MAXENT.MAXENT = user.MAXENT, 
# RF.binary.randomForest.randomForest = user.RF
# 
#)


#opt.u <- bm_ModelingOptions(data.type = 'binary',
#                            models = c("CTA", "ANN"),#,"FDA", "GBM", "GLM","MARS", "MAXENT", "RF"),
#                            strategy = 'user.defined',
#                            user.val = user.val,
#                            bm.format = myBiomodData,
#                            calib.lines  = NULL )

#opt.u@options$CTA.binary.rpart.rpart@args.values


################################## Custom model parameter optimization

#myBiomodOption<-bm_ModelingOptions(
#            data.type = "binary",
#             models = c("GLM","GAM", "GBM"),
#              strategy = 'user.defined',
#               user.val = user.val,
#                #user.base = "bigboss",
#                 bm.format = myBiomodData
#                
#)
#myBiomodOption




#---------------------------------------------------------------------------------------- Create different algorithm single models
#windowsFonts(Times = windowsFont("Times New Roman")) # Set plotting font

myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    #modeling.id = as.character(format(Sys.time(), "%s")),    # Set current time as unique model run ID
                                    models = c("CTA","GBM", "RF","MAXENT.Phillips.2"),
                                    CV.strategy = "kfold",
                                    CV.nb.rep = 5,                                      # Number of model repetitions
                                    CV.perc = 0.7,                                           # Proportion of training dataset
                                    CV.k=  5,                                                 # Number of K-fold cross-validation
                                    var.import = 3,                                          # Variable importance calculation method
                                    OPT.strategy = "bigboss",                                # Parameter optimization strategy                  
                                    #OPT.user = myBiomodOption,                             # Custom model parameters
                                    metric.eval = c('TSS','ROC'),                          # Model evaluation metrics
                                    seed.val = 123)                                      # Random seed
                                    


# nb.cpu = 1)                                       # Number of CPU cores for parallel running

myBiomodModelOut

######################## Evaluate single model scores and variable importance
# Model evaluation scores
#evaluation_data <-get_evaluations(myBiomodModelOut)       
#write.csv(evaluation_data, "D:/7_Paper/hyh/result/lulc2019/angle_models_evaluation_results.csv", row.names = FALSE)
# Model variable importance
#var_important_data<-get_variables_importance(myBiomodModelOut)   
#write.csv(var_important_data, "D:/7_Paper/hyh/result/lulc2019/angle_models_var_results.csv", row.names = FALSE)
# Model prediction results
#predictions_data<-get_predictions(myBiomodModelOut) 
#write.csv(predictions_data, "D:/7_Paper/hyh/result/lulc2019/angle_models_pridiction_results.csv", row.names = FALSE)


######################### Plot model results
# Comparison of model evaluation scores
bm_PlotEvalMean(bm.out= myBiomodModelOut,
                metric.eval = c("TSS","roc"),
                dataset = "calibration", 
                group.by = 'algo', 
                do.plot = TRUE
)


#bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'algo')) 

#### Variable importance
#bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c( 'algo', 'algo','expl.var')) # Group by variable

#bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c( 'expl.var','algo', 'algo')) # Group by model



##### Plot response curves
#myGLM1 <- BIOMOD_LoadModels(myBiomodModelOut, 
#                            PA = "PA1",         # Cross-validation partition
#                            run = "RUN1",      # Number of runs
#                            algo = "RF")         # Algorithm 

# myGLM <- BIOMOD_LoadModels(myBiomodModelOut, run = "RUN1") 
# myGLM <- BIOMOD_LoadModels(myBiomodModelOut, PA = "PA1") 
# myGLM <- BIOMOD_LoadModels(myBiomodModelOut, algo = "RF") 

# myGML<-get_built_models(myBiomodModelOut)[c(1:3, 12:14)]# Plot curves for specific models


#bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
#                      show.variables = c("bio2"),# Plot curves for specific variables
#                      models.chosen = "all",
#                      #fixed.var = 'min', # Plotting method
#                      fixed.var = 'median') 


#bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
#                      models.chosen = get_built_models(myBiomodModelOut)[3],
#                      fixed.var = 'median',
#                      do.bivariate = TRUE)



# -------------------------------------------------------------------------------—————————— Ensemble models and projections
##### Project single models
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,     # Single model output results
                                  proj.name = "Current",      # Result path name
                                  new.env = myExpl,          # Projection variable set
                                  models.chosen = 'all',    # Models to project
                            #     models.chosen = myGBM,
                                  metric.binary = 'all',    # Metrics used to convert to presence/absence
                                  metric.filter = 'all',      # Metrics used for model projection optimization
                                  overwrite=TRUE,              # Overwrite command
                                  build.clamping.mask = TRUE)
                            
#     ,on_0_1000 = FALSE)             # Whether to convert model prediction probabilities (0-1) to 0-1000 for space saving

#myBiomodProj


##### Select models of specific algorithms for ensembling or plotting response curves
#myMAXENT <- BIOMOD_LoadModels(myBiomodModelOut, algo = ("MAXNET")) 
#myRF <- BIOMOD_LoadModels(myBiomodModelOut, algo = ("RF")) 
#myCTA <- BIOMOD_LoadModels(myBiomodModelOut, algo = ("CTA")) 
#myFDA <- BIOMOD_LoadModels(myBiomodModelOut, algo = ("FDA"))
#myMARS <- BIOMOD_LoadModels(myBiomodModelOut, algo = ("MARS")) 
#myGBM <- BIOMOD_LoadModels(myBiomodModelOut, algo = ("GBM")) 


##### Ensemble models
myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      # models.chosen = myGBM,
                                      models.chosen ='all',                     # Models to include in the ensemble
                                      em.by = 'all',                            # Specify to ensemble all algorithms
                                      em.algo = c('EMwmean'),#em.algo = c('EMmean', 'EMcv', 'EMci', 'EMmedian', 'EMca', 'EMwmean'),  # Ensemble model algorithms
                                      #metric.select = c('ROC'),                  # Model metrics for selecting models to ensemble
                                      # metric.select.thresh = c(0.8),            # Threshold for model metrics when selecting models to ensemble
                                      metric.eval = c('TSS', 'ROC' ),          ## Evaluation metrics for ensemble models
                                      var.import = 4,                          # Variable importance calculation method
                                      EMci.alpha = 0.05,                      ## 
                                      EMwmean.decay = 'proportional')          ### Decay method for weighted average model ensembling

#myBiomodEM
#folder_path <- "D:/7_Paper/hyh/var/c_v/a2019"
#rasters <- list.files(path = folder_path, pattern = ".tif", full.names = TRUE)
# Stack rasters
#myExpl <- stack(rasters)
################################################### Project ensemble models (using single-model projection results and built ensemble models)
##### Use single-model projection results
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
                                             output.format = ".tif",           # Output format
                                             bm.proj = myBiomodProj,            # Single-model projection results
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
                                     




#       on_0_1000 = FALSE)                

#myBiomodEMProj
plot(myBiomodEMProj)


 ######################### Convert probability distribution map to landscape resistance surface using negative exponential function
library(raster)
library(rasterVis)

ensemble_raster <- raster("D:/7_Paper/hyh/Circuitscape/Assistance/HS_lulc_2009.tif")
plot(ensemble_raster)

# Negative exponential function conversion
resistance_raster <- 1000^((-1) * ensemble_raster)
plot(resistance_raster)

#writeRaster(resistance_raster, filename = "C:/Users/10284/Desktop/smaxent/circure/resistance_raster.tif", format = "GTiff")

# Calculate minimum and maximum values
min_resistance <- minValue(resistance_raster)
max_resistance <- maxValue(resistance_raster)

# Normalize landscape resistance values to 0-100
normalized_resistance <- (resistance_raster - min_resistance) / (max_resistance - min_resistance) * 99 + 1
plot(normalized_resistance)

# Export raster
writeRaster(normalized_resistance, filename="D:/7_Paper/hyh/Circuitscape/note_asc/resistance_lulc_2009.tif", format="GTiff")


#————————————————————————————————————————————————————————————————————————————————————————————————
#————————————————————————————————————————————————————————————————————————————————————————————————
library(raster)
library(terra)
rm(list = ls())

# Read model prediction probability distribution map (0-1)
raster_file <- "D:/7_Paper/hyh/Circuitscape/Assistance/binary_lulc_2019.tif"
pre <- raster(raster_file)

# Set coordinate system
proj4string(pre)
crs(pre) <- crs("+proj=utm +zone=48 +datum=WGS84")
threshold <-  0.5 # Set threshold

binary_map <- pre > threshold
binary_map <- rast(binary_map)  ## Convert to binary distribution map

plot(binary_map)

pach_raster<-patches(binary_map, directions=8, zeroAsNA=TRUE, allowGaps=TRUE)# Identify patches


pach_n <- zonal(cellSize(pach_raster, unit="km"), pach_raster, sum, as.raster=TRUE)

opt_raster <- ifel(pach_n <= 46.19, 0, pach_n)# Select patches - only select patches with value 1
#opt_raster <- ifel(pach_n <= 22.56, 0, 1)  # Select patches
opt_raster[is.na(pach_n)] <- 0  # Assign 0 to invalid pixels
plot(opt_raster)

writeRaster(opt_raster, filename=file.path( "D:/7_Paper/hyh/MAR_habitat/2019_lulc_MAR_habitat.tif"), overwrite=TRUE)








###################################################################################### Marginal response curve

library(ggplot2)
library(dplyr)
# Read CSV file
data <- read.csv("D:/7_Paper/hyh/result/response_curves_data/PCA3.csv")


# Ensure no NA values in the data frame
data <- na.omit(data)  # Remove rows with NA

# Ensure expl_val and pred_val are numeric
data$expl.val <- as.numeric(data$expl.val)
data$pred.val <- as.numeric(data$pred.val)


p <- ggplot(data = data) +
  # Plot scatter plot
  geom_smooth(mapping = aes(x = expl.val, y = pred.val, color = Time)) +
  
  
  
  # Custom colors
  scale_color_manual(values = c("red", "blue", "green")) +
  
  # Add ticks
  scale_y_continuous(n.breaks = 5) +
  scale_x_continuous(n.breaks = 7) +
  
  # Change tick label font
  theme(axis.text = element_text(family = "serif", colour = "black",size = 15)) +
  
  # Show ticks and add borders
  theme(axis.ticks = element_line(linewidth = 0.5, colour = "black")) +
  theme(axis.title.x = element_text(family = "serif", size = 30)) +
  # theme(axis.title.y = element_text(family = "serif", size = 20)) +
  
  # Add borders, remove internal fill color
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.background = element_rect(fill = NA),  # Remove panel background fill
        plot.background = element_rect(fill = NA)) +  # Remove overall plot background fill
  # Remove legend
  #theme(legend.position = "none") +
  # Label settings
  labs(x = "PCA3",y = "Presence probability", title = NULL) #,y = NULL

# Print the graph
print(p)


######################################################################################### Area change bar chart
library(ggplot2)
library(tidyr)
library(readr)

# Read CSV file
data <- read_csv("D:/7_Paper/hyh/result/csv/change_area.csv")

# Convert data to long format
data_long <- pivot_longer(data, cols = c("Stabled habitat", "Lost habitat", "Gained habitat"),
                          names_to = "Habitat_Type", values_to = "Value")

# Plot bar chart
p <- ggplot(data_long, aes(x = Scenario, y = Value, fill = Habitat_Type)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(), alpha = 0.7) +   
  facet_grid(~ Periods) +   
  theme_minimal() +  
  labs(x = "Models",
       y = expression("Area (Km"^2*")")) + 
  scale_fill_manual(values = c("Gained habitat" = "#B883D4", 
                               "Lost habitat" = "#FA7F6F", 
                               "Stabled habitat" = "#32B897")) +   # Custom fill colors
  geom_text(aes(label = Value), 
            position = position_dodge(0.75), 
            vjust = -0.5,  # Adjust text position, -0.5 means text above the bar
            size = 4) +  # Text size
  theme(legend.position = "none",
    text = element_text(size = 12),  # Set global font size
    axis.text = element_text(size = 15),  # Set axis tick font size
    axis.title = element_text(size = 16),  # Set axis title font size
    axis.title.x = element_text(size = 20),  # Set X-axis title font size
    axis.title.y = element_text(size = 20),  # Set Y-axis title font size
    #axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate X-axis labels
    #legend.position = "top",  # Legend position
    legend.text = element_text(size = 14),  # Set legend item font size
    legend.title = element_text(size = 16),  # Set legend title font size
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center title, increase font size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border
    #strip.background = element_rect(fill = "lightgrey", color = "black", size = 0.5),  # Set facet background
   # strip.text = element_text(face = "bold", size = 15),  # Set facet title style
    axis.ticks.length = unit(0.25, "cm"),  # Set axis tick line length
    axis.ticks = element_line(color = "black")  # Set tick line color to black
  )

print(p)




p <- ggplot(data_long, aes(x = Scenario, y = Value, fill = Habitat_Type)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(), alpha = 0.7) +   
  facet_wrap(~ Periods, scales = "free_y", nrow = 2) +  # Use facet_wrap to split into upper and lower panels
  theme_minimal() +  
  labs(x = "Models",
       y = expression("Area (Km"^2*")")) + 
  scale_fill_manual(values = c("Gained habitat" = "#B883D4", 
                               "Lost habitat" = "#FA7F6F", 
                               "Stabled habitat" = "#32B897")) +   # Custom fill colors
  geom_text(aes(label = Value), 
            position = position_dodge(0.75), 
            vjust = -0.5,  # Adjust text position, -0.5 means text above the bar
            size = 4.5) +  # Text size
  theme(
    text = element_text(size = 12),  # Set global font size
    axis.text = element_text(size = 13),  # Set axis tick font size
    axis.title = element_text(size = 13),  # Set axis title font size
    axis.title.x = element_text(size = 20),  # Set X-axis title font size
    axis.title.y = element_text(size = 20),  # Set Y-axis title font size
    axis.text.x = element_text( hjust = 1),  # Rotate X-axis labels
    legend.position = "none",  # Legend position
    legend.text = element_text(size = 14),  # Set legend item font size
    legend.title = element_text(size = 16),  # Set legend title font size
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center title, increase font size
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add border
    strip.background = element_rect(fill = "lightgrey", color = "black", size = 0.5),  # Set facet background
    strip.text = element_text(face = "bold", size = 15),  # Set facet title style
    axis.ticks.length = unit(0.25, "cm"),  # Set axis tick line length
    axis.ticks = element_line(color = "black")  # Set tick line color to black
  )


print(p)
ggsave("C:/Users/Admin/Desktop/面积变化3.pdf", plot = p, width = 5, height = 9)



# Plot bar chart
p <- ggplot(data_long, aes(x = Scenario, y = Value, fill = Habitat_Type)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(), alpha = 0.7) +   
  facet_grid(~ Periods) +   
  theme_minimal() +  
  labs(x = "Models",
       y = expression("Area (Km"^2*")")) + 
  #scale_fill_manual(values = c("Gained habitat" = "#B883D4", 
  ##                             "Lost habitat" = "#FA7F6F", 
   #                            "Stabled habitat" = "#32B897")) +   # Custom fill colors
  scale_fill_manual(values = c("Gained habitat" = rgb(197/255,0/255, 255/255), 
                               "Lost habitat" = rgb(255/255, 0/255, 0/255), 
                               "Stabled habitat" = rgb(85/255, 255/255, 0/255))) +  # Use correct RGB values
  geom_text(aes(label = Value), 
            position = position_dodge(0.75), 
            vjust = -0.5,  # Adjust text position, -0.5 means text above the bar
            size = 4) +  # Text size
  theme(legend.position = "none",
        text = element_text(size = 12),  # Set global font size
        axis.text = element_text(size = 15),  # Set axis tick font size
        axis.title = element_text(size = 16),  # Set axis title font size
        axis.title.x = element_text(size = 20),  # Set X-axis title font size
        axis.title.y = element_text(size = 20),  # Set Y-axis title font size
        #axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate X-axis labels
        #legend.position = "top",  # Legend position
        legend.text = element_text(size = 14),  # Set legend item font size
        legend.title = element_text(size = 16),  # Set legend title font size
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center title, increase font size
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        #panel.border = element_rect(color = "black", fill = NA, size = 1),  # Remove border
        axis.ticks.length = unit(0.25, "cm"),  # Set axis tick line length
        axis.ticks = element_line(color = "black"),  # Set tick line color to black
        axis.line = element_line(color = "black",size = 1)  # Add axis lines
  )

print(p)
ggsave("C:/Users/Admin/Desktop/面积变化3.pdf", plot = p, width = 8, height = 6)


########################################################################################################
########################################################################################################
##Spatiotemporal connectivity model 

# Load required libraries
library(sf)      # Read SHP files
library(raster)  # Process raster data
library(gdistance) # Calculate cost distance
library(sp)      # Spatial point operations
library(terra)

# 1. Read the SHP file (filename: t1t2.shp)
nodes_sf <- st_read("C:/Users/Admin/Desktop/IZ/conefor_R/MAR/t2t3/t2t3_0.8rm.shp")  # Replace with your file path (e.g., "data/t1t2.shp")

# 2. Check and set the projection (WGS84 UTM Zone 48N)
target_crs <- "+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs"
if (st_crs(nodes_sf) != target_crs) {
  nodes_sf <- st_transform(nodes_sf, crs = target_crs)
}

# 4. Extract node ID, area, and centroid coordinates
nodes <- data.frame(
  Node_ID = nodes_sf$Node_ID,  # Unique node identifier
  area_ha = nodes_sf$shape_Area  # Patch area (unit: hectares)
)
coordinates <- st_coordinates(nodes_sf)  # Extract centroid coordinates (X, Y)

# 5. Load and process the resistance surface (ass.tif, ensure consistent projection)
res <- raster("C:/Users/Admin/Desktop/IZ/conefor_R/resistance/resistance_t2t3.tif") 
# Convert the resistance surface projection to UTM Zone 48N

con <- 1 / res  # Convert to conductivity layer (reciprocal of resistance, higher values mean easier passage)
tr <- transition(con, transitionFunction = mean, directions = 8)  # 8-neighborhood average conductivity
tr <- geoCorrection(tr, multpl = FALSE, scl = FALSE)  # Geographical correction

# Calculate the median of the resistance surface
values_terra <- values(res)
global_median <- median(values_terra, na.rm = TRUE)
global_median

# 6. Set diffusion distance parameters (32km)
dist_medium_m <- 32000  # Convert 32km to meters (32×1000)
eff_dist <- dist_medium_m * 70  # Effective distance
probmin <- 0.001  # Minimum effective probability (retain connections with probability ≥0.1%)
k <- -log(0.5) / eff_dist  # Probability decay coefficient (distance-probability relationship)

# 7. Filter source nodes (loss=1 or stable=12) and target nodes (gain=2 or stable=12)
source_nodes <- nodes_sf[nodes_sf$change %in% c(1, 3), ]  # Source nodes: Disappeared/Stable
dest_nodes <- nodes_sf[nodes_sf$change %in% c(2, 3), ]    # Target nodes: Newly added/Stable

# 8. Convert to SpatialPoints objects (Conefor requires classic spatial point format)
xy_SP1 <- SpatialPoints(coords = cbind(source_nodes$x, source_nodes$y), proj4string = CRS("+proj=utm +zone=48 +datum=WGS84")) 
xy_SP2 <- SpatialPoints(coords = cbind(dest_nodes$x, dest_nodes$y), proj4string = CRS("+proj=utm +zone=48 +datum=WGS84"))

# 9. Calculate the minimum cost distance matrix (cumulative resistance for movement between nodes)
dist_mat <- costDistance(tr, xy_SP1, xy_SP2)
dist_table <- as.data.frame(dist_mat)  # Convert to data frame

# 10. Filter valid connections (distance less than maximum effective distance)
dist_max <- -log(probmin) / k  # Maximum distance corresponding to probmin
valid_links <- which(dist_table < dist_max, arr.ind = TRUE)  # Extract (i,j) indices

# 11. Construct the connection data frame (source, destination, distance)
dist <- data.frame(
  source_id = source_nodes$Node_ID[valid_links[, 1]],  # Source node ID
  dest_id = dest_nodes$Node_ID[valid_links[, 2]],      # Target node ID
  distance = as.numeric(dist_table[valid_links])       # Cost distance
)

# 12. Calculate spatial diffusion probability and integrate temporal probability
dist$spatial_prob <- exp(-k * dist$distance)  # Spatial diffusion probability (negative exponential function)
dist <- dist[dist$spatial_prob >= probmin, ]  # Filter low-probability connections

# Add node status (for temporal probability calculation)
dist$status_source <- nodes_sf$change[match(dist$source_id, nodes_sf$Node_ID)]
dist$status_dest <- nodes_sf$change[match(dist$dest_id, nodes_sf$Node_ID)]
dist$status_link <- paste(dist$status_source, dist$status_dest)  # Status combination

# Set temporal probability (based on node status combination)
dist$temporal_prob <- ifelse(
  dist$status_link %in% c("3 3", "3 2", "1 3"),  # Stable-Stable, Stable-Newly added, Disappeared-Stable
  1,  # Temporal probability is 1 (direct connection)
  ifelse(dist$status_link == "1 2", 0.5, 0)  # Disappeared-Newly added: Temporal probability 0.5, other invalid connections 0
)
dist$final_prob <- dist$spatial_prob * dist$temporal_prob  # Final effective probability

# 13. Prepare Conefor input files (nodes and connection probabilities)
# Node file format: Node_ID (integer)  area_ha (numeric)
write.table(
  nodes,
  file = "C:/Users/Admin/Desktop/IZ/conefor_R/MAR_result/nodes_t2t3_0.8rm.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = " "  # Conefor typically uses space-separated values
)

# Connection probability file format: source destination prob
write.table(
  dist[, c("source_id", "dest_id", "final_prob")],
  file = "C:/Users/Admin/Desktop/IZ/conefor_R/MAR_result/probs_t2t3_0.8rm.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = " "
)

# 14. Run Conefor (ensure the program path is correct)
setwd("C:/Users/Admin/Desktop/IZ/conefor_R/MAR_result")
shell("C:/Users/Admin/Desktop/IZ/conefor_R/MAR_result/conefor.exe -nodeFile nodes_t2t3_0.8rm.txt -conFile probs_t2t3_0.8rm.txt -t prob notall -PC -double -noout -* -wprobdir -wprobmax")







