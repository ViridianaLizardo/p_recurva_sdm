#' Species Distribution Modeling (SDM) Functions
#' 
#' This script contains functions to model species distributions using MaxEnt.
#' The models are trained in the species' native range and projected to America.

library(usdm)
library(RPushbullet)

# Load environmental raster data

wc_native <- rast('data/wc/native.tif')
wc_america <- rast('data/wc/america.tif')

# Apply Variance Inflation Factor (VIF) to remove collinearity
wc_america_vif <- rast('data/wc/america.tif') %>%
  vifstep(.,th = 10,
          keep = 'temp_mean_annual') %>%
  exclude(wc_america, .)

sp_name<- "Corymbia citriodora"

library(usdm)
library(RPushbullet)

# Function to model species distribution in native range and transfer to America
model_host<- function(sp_name){
      
      env <- rast(paste0('data/host_sdm/env/', sp_name, '.tif')) 
        
      # Load and process occurrence data
      
      p <- thinData(coords = read.csv(paste0('data/host_sdm/occ/', 
                                             sp_name, '.csv'))[,2:3],
                    env = env,x = 'X',y = 'Y',verbose = F)
      
      bg <- read.csv(paste0('data/host_sdm/bg/', sp_name, '.csv'))
      
      data <- prepareSWD(species = sp_name,
                         p = p,
                         a = bg,
                         env = env, verbose = F) 
      
      # Split data into training and testing
      data_split <- trainValTest(x = data,test = 0.15)
      
      train_default <- train(method = 'Maxent',data = data_split[[1]] )
      
      
      # Hyperparameter tuning
      
      h <- list(fc = c('lp', 'lq','pq', 
                       'lpq', 'lph', 'pqh', 'lpqh'),
                reg = seq(0.5,1.5,.5))
      
      optimized <- gridSearch(train_default, 
                              hypers = h, 
                              metric = "aicc",
                              env = env,
                              test = data_split[3],
                              save_models = F,
                              interactive = F)
      
      index <- which.min(optimized@results$AICc)
      
      # Train final model
      final_model <- train("Maxent", 
                           data = data_split[[1]], 
                           fc = optimized@results[index, 1],
                           reg = optimized@results[index, 2])
      
      # Prepare transfer layers 
      proj_native <- terra::subset(wc_native,
                                   names(env))
      
      proj_america <- terra::subset(wc_america,
                                    names(env))
      
      
      # Generate evaluation report
      eval_data <- data.frame(species = sp_name,
                              train = sum(data_split[[1]]@pa == 1),
                              test = sum(data_split[[2]]@pa == 1),
                         tss = tss(final_model, data_split[[2]]),
                         auc = auc(final_model, data_split[[2]]),
                         fc = optimized@results[index, 1],
                         reg = optimized@results[index, 2],
                         iter = optimized@results[index,3],
                         vars = paste0(names(env),
                                       collapse =', '),
                         trained = 'native')
      write.csv(eval_data,
                file = paste0('data/host_sdm/eval/',sp_name, '_native_tr.csv'),
                row.names = F)

      folder_native <- paste0('data/host_sdm/SDMtune report/native/', sp_name)
      folder_america <- paste0('data/host_sdm/SDMtune report/america/', sp_name)
      
      erase_and_rewrite <- function(folder_path) {
        if (dir.exists(folder_path)) {
          # Eliminar archivos dentro de la carpeta
          files <- list.files(folder_path, full.names = TRUE, recursive = TRUE)
          file.remove(files[file.info(files)$isdir == FALSE])
          
          # Eliminar subcarpetas después de los archivos
          unlink(folder_path, recursive = TRUE, force = TRUE)
          message('Previous report erased\n')
          
        }else{message('creating reports\n')}}
      
      folder_native <- paste0('data/host_sdm/SDMtune report/native/', sp_name)
      folder <- paste0('data/host_sdm/SDMtune report/america/', sp_name)
      
      
      erase_and_rewrite(folder_native)
      SDMtune::modelReport(model = final_model,
                           folder = folder_native,
                           test = data_split[[2]],type = 'cloglog', 
                           only_presence = T, 
                           jk = T, 
                           response_curves = T,clamp = F,
                           env = proj_native,verbose = F)
      
      erase_and_rewrite(folder)
      SDMtune::modelReport(model = final_model,
                           folder = folder,
                           test = data_split[[2]],type = 'cloglog',
                           only_presence = T, 
                           jk = T, 
                           response_curves = T,clamp = F,
                           env = proj_america,verbose = F)
      
      
      # Generate binary presence-absence maps
      ths <- thresholds(final_model, type = 'cloglog', test = data_split[[2]])
      
      ths_val <- ths %>% 
        filter(Threshold == 'Minimum training presence') %>% 
        pull('Cloglog value')
      write.csv(ths,file = paste0("data/host_sdm/thresholds/", sp_name, '_native_tr.csv'))
      
      ## america
      predict(final_model, proj_america,  type = "cloglog", clamp = F) %>% 
        plotPA(., colors = c('red', 'grey70'),ths_val, hr = F, overwrite = T, 
               filename = paste0("data/host_sdm/binary/america/", sp_name, '_native_tr.tif')
               )
      
      ## native
      predict(final_model, proj_native,  type = "cloglog", clamp = F) %>% 
        plotPA(., colors = c('red', 'grey70'),ths_val, hr = F, overwrite = T, 
               filename = paste0("data/host_sdm/binary/native/", sp_name, '_native_tr.tif')
               )

      message(sp_name)
      }

## Model in America
model_host_america<- function(sp_name){
  
  env <- wc_america_vif
  
  # Preparar datos
  

  p <- thinData(coords = read.csv(paste0('data/host_sdm/occ_america/', 
                             sp_name, '.csv'))[,2:3],
                env = env,x = 'X',y = 'Y',verbose = F)
  bg <- read.csv('data/host_sdm/america_bg.csv')
  
  data <- prepareSWD(species = sp_name,
                     p = p,
                     a = bg,
                     env = env, verbose = F) 
  
  ## Separar datos
  data_split <- trainValTest(x = data,test = 0.15)
  
  train_default <- train(method = 'Maxent',data = data_split[[1]] )
  
  
  # Fine tunning
  
  h <- list(fc = c('lp', 'lq','pq', 
                   'lpq', 'lph', 'pqh', 'lpqh'),
            reg = seq(0.5,1.5,.5))
  
  optimized <- gridSearch(train_default, 
                          hypers = h, 
                          metric = "aicc",
                          env = env,
                          test = data_split[3],
                          save_models = F,
                          interactive = F)
  
  index <- which.min(optimized@results$AICc)
  
  # Modelo final
  final_model <- train("Maxent", 
                       data = data_split[[1]], 
                       fc = optimized@results[index, 1],
                       reg = optimized@results[index, 2])
  
  # Proyectar en mapas 
  proj_native <- terra::subset(wc_native,
                               names(env))
  
  proj_america <- terra::subset(wc_america,
                                names(env))
  
  
  #Reporte
  eval_data <- data.frame(species = sp_name,
                          train = sum(data_split[[1]]@pa == 1),
                          test = sum(data_split[[2]]@pa == 1),
                          tss = tss(final_model, data_split[[2]]),
                          auc = auc(final_model, data_split[[2]]),
                          fc = optimized@results[index, 1],
                          reg = optimized@results[index, 2],
                          iter = optimized@results[index,3],
                          vars = paste0(names(env),
                                        collapse =', '),
                          trained = 'america')
  write.csv(eval_data,
          file = paste0('data/host_sdm/eval/',sp_name, '_america.csv'),
          row.names = F)
  
  erase_and_rewrite <- function(folder_path) {
    if (dir.exists(folder_path)) {
      # Eliminar archivos dentro de la carpeta
      files <- list.files(folder_path, full.names = TRUE, recursive = TRUE)
      file.remove(files[file.info(files)$isdir == FALSE])
      
      # Eliminar subcarpetas después de los archivos
      unlink(folder_path, recursive = TRUE, force = TRUE)
      message('Previous report erased\n')
      
    }else{message('creating reports\n')}}
  
  folder_native <- paste0('data/host_sdm/SDMtune report america/native/', sp_name)
  folder_america <- paste0('data/host_sdm/SDMtune report america/america/', sp_name)
  
  
  erase_and_rewrite(folder_native)
  SDMtune::modelReport(model = final_model,
                       folder = folder_native,
                       test = data_split[[2]],type = 'cloglog', 
                       only_presence = T, 
                       jk = T, 
                       response_curves = T,clamp = F,
                       env = proj_native,verbose = F)
  
  erase_and_rewrite(folder_america)
  SDMtune::modelReport(model = final_model,
                       folder = folder_america,
                       test = data_split[[2]],type = 'cloglog',
                       only_presence = T, 
                       jk = T, 
                       response_curves = T,clamp = F,
                       env = proj_america,verbose = F)
  
  
  ## Thresholds
  ths <- thresholds(final_model, type = 'cloglog', test = data_split[[2]])
  
  ths_val <- ths %>% 
    filter(Threshold == 'Minimum training presence') %>% 
    pull('Cloglog value')
  write.csv(ths,file = paste0("data/host_sdm/thresholds/", sp_name, '_america.csv'))
  #Binary
  
  ## america
  predict(final_model, proj_america,  type = "cloglog", clamp = F) %>% 
    plotPA(., colors = c('red', 'grey70'),ths_val, hr = F, overwrite = T, 
           filename = paste0("data/host_sdm/binary/america/", sp_name, '_america.tif')
    )
  
  ## native
  predict(final_model, proj_native,  type = "cloglog", clamp = F) %>% 
    plotPA(., colors = c('red', 'grey70'),ths_val, hr = F, overwrite = T, 
           filename = paste0("data/host_sdm/binary/native/", sp_name, '_america.tif')
    )
  
  message(sp_name)
}


# List of species to model
especies_names<- c("Corymbia citriodora",
                   "Eucalyptus camaldulensis",
                   "Eucalyptus cinerea",
                   "Eucalyptus globulus",
                   "Eucalyptus robusta")


lapply(especies_names, model_host)
lapply(especies_names, model_host_america)


# Send notification when finished
pbPost(type = 'note',title = 'acabó R', body = 'revisa los resultados')
