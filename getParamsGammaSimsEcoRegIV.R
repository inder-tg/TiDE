
# Written by Inder Tecuapetla on Dec 15, 2018
# Focused on 2004 (although EXP and GAMMA fits are also available for 2006)
# Get dist of rainfall in 4 ecoRegions within the SemiArid Southern Uplands class
# Get pixels by ecoRegion
# Get params of Exponential and Gamma fit on dist of rainfall

rm(list=ls())
library(raster)
library(ff)
source("myFunctions.R")

listFiles <- list.files(pattern = ".tif")

ecoRegionIV <- raster(listFiles[2])

index <- ff(vmode = "double", dim = c( nrow(ecoRegionIV), ncol(ecoRegionIV), 1 ),
            filename = paste0(getwd(), "/ecoRegionIV.ffdata"))

index[, , 1] <- ecoRegionIV[[1]][]
index_array <- index[, , 1]
save(index_array, file = paste0(getwd(), "/ecoRegionIV.RData"))


matEcoReg <- LoadToEnvironment(paste0(getwd(), "/ecoRegionIV.RData"))$index_array

# ecoLabels <- c("Great plains", "North American deserts", "California mediterranean",
#                "semi-arid Southern uplands", "temperate mountains", "dry tropical forest",
#                "tropical rain forest")

ecoLabels <- c("Madrense", "interiorPlains", "Mezquital", "plateauPlains")

mat_precipitation <- LoadToEnvironment(paste(getwd(), "/data/mat_precipitation_", 2004, ".RData", sep=""))
rain <- mat_precipitation$v

lengthMonths <- readRDS("lengthMonthsLeapYear")

if(dim(rain)[3] == 365){
  lengthMonths[2] <- 28
}

cumsumNumDays <- cumsum(lengthMonths) 

monthsLength <- c(1, cumsumNumDays)

# labelsEco <- c("GreatPlains", "NorthAmericanDeserts", "CalifMed",
#                "SemiSouthernUp", "TempMount", "DryTropicalForest",
#                "TropicalRainForest")

labelsEco <- ecoLabels

monthsName <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
                "Oct", "Nov", "Dec")

# =============================================================================
# let's get the rainfall amount by ecoRegion
# =============================================================================

for( k in 1:4 ){
  
  if( k == 1 ){
    ecoType <- 12111
  }
  
  if( k == 2 ){
    ecoType <- 12121
  }
  
  if( k ==3 ){
    ecoType <- 12211
  }
  
  if( k == 4 ){
    ecoType <- 12212
  }
  
  distMonthPREC <- vector("list", 12)
  
  for(i in 2:13){
    
    cat("Month:", monthsName[i-1], "\n")
    
    distMonth <- c()
    # test <- matrix(NA, nrow = nrow(ecoRegion), ncol = ncol(ecoRegion))
    for(j in 1:nrow(rain)){
      if(j %% 100 == 0){
        cat("row:", j, "\n")
      }
      
      prelimValidPixels <- sapply( 1:ncol(matEcoReg), function(s) ifelse(matEcoReg[j,s] == ecoType, 1, 0))
      
      validPixels <- (1:ncol(rain))[prelimValidPixels==1]
      
      # test[j,validPixels] <- 1
      
      for(t  in 1:length(validPixels)){
        distMonth <- c(distMonth,  rain[j, validPixels[t], monthsLength[(i-1)]:monthsLength[i]])
      }
    }
    
    # tempRaster <- matrixToRaster(test, ecoRegion)
    
    distMonthPREC[[i-1]] <- distMonth
  }
  
  save(distMonthPREC, file = paste0(getwd(), "/RData/PREC_", labelsEco[k], ".RData"))
}


# =============================================================================
# Let's get the PREC pixels in each ecoRegion
# =============================================================================

matEcoReg[is.na(matEcoReg)]<-0

for( k in 1:4){
  
  if( k == 1 ){
    ecoType <- 12111
  }
  
  if( k == 2 ){
    ecoType <- 12121
  }
  
  if( k ==3 ){
    ecoType <- 12211
  }
  
  if( k == 4 ){
    ecoType <- 12212
  }
  
  cat("EcoRegion:", labelsEco[k], "\n")
  tempValidPixels <- c()
  for(i in 1:nrow(rain)){
    prelimValidPixels <- sapply( 1:ncol(matEcoReg), 
                                 function(s) ifelse(matEcoReg[i,s] == ecoType, 1, 0))
    
    validPixels <- (1:ncol(rain))[prelimValidPixels==1]
    
    tempValidPixels <- c(tempValidPixels, validPixels)
  }
  
  pixelsTemp <- matrix(NA, nrow = length(tempValidPixels), ncol = dim(rain)[3])
  cont <- 0
  for(i in 1:nrow(rain)){
    
    if(i %% 100 == 0){
      cat("Row:", i, "\n")
    }
    
    prelimValidPixels <- sapply( 1:ncol(matEcoReg), 
                                 function(s) ifelse(matEcoReg[i,s] == ecoType, 1, 0) )
    
    if( sum(prelimValidPixels) > 0 ){
      validPixels <- (1:ncol(rain))[prelimValidPixels==1]
      
      for(s in 1:length(validPixels)){
        cont <- cont + 1
        pixelsTemp[cont,] <- rain[i, validPixels[s], ]
      }
    }
    
  }
  
  save(pixelsTemp, file = paste0(getwd(), "/RData/pixels_PREC_", labelsEco[k], ".RData"))
}


# =============================================================================
# Let's get the gamma fit on the monthly PREC data
# =============================================================================

library(nloptr)

GammaNLL <- function(pars, data){
  alpha <- pars[[1]]
  theta <- pars[[2]]
  return (-sum(dgamma(data, shape = alpha, scale = theta, log = TRUE)))
}

for(ecoType in 1:4){
  monthlyGammaParamTemp <- matrix(0, nrow = 12,  ncol = 2)
  monthlyExpParamTemp <- numeric(12)
  
  distMonthPREC <- LoadToEnvironment(paste0(getwd(), "/RData/PREC_", labelsEco[ecoType], ".RData"))$distMonthPREC
  
  cat("EcoRegion:", labelsEco[ecoType], "\n")
  
  for(k in 1:12){
    amountRain <- distMonthPREC[[k]]
    
    if( sum(amountRain, na.rm = T) != 0 ){
      amountRain <- amountRain[!is.na(amountRain)]
      
      png( paste0(getwd(), "/png_ecoRegionIV/exponential_gamma_", monthsName[k], "_", labelsEco[ecoType], ".png" ) )
      hist(amountRain, prob = T, main = paste("Hist. of PRECIPITATION in ", 
                                              monthsName[k], sep = ""))
      lambda <- 1 / mean(amountRain)
      
      curve(dexp(x, rate = lambda), col = 2, lty = 2, lwd = 4, add = TRUE)
      
      Scale <- var(amountRain) / mean(amountRain)
      Shape <- mean(amountRain) / Scale
      Fit <- nloptr(x0 = c(Shape, Scale), eval_f = GammaNLL, lb = c(0,0), data = amountRain,
                    opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1e5))
      
      curve(dgamma(x, shape = Fit$solution[1], scale = Fit$solution[2]), col = 3, 
            lty = 3, lwd = 4, add = TRUE)
      
      legend("topright", fill = c(2,3), legend = c("Exp", "Gamma"),
             bty = "n", xpd = NA, cex = 1.75)
      dev.off()
      
      monthlyExpParamTemp[k] <- lambda
      monthlyGammaParamTemp[k,] <- c(Fit$solution[1], Fit$solution[2])
    }
    
  }
  
  save(monthlyGammaParamTemp, monthlyExpParamTemp, 
       file = paste0(getwd(), "/RData/paramsRAINFALL_", labelsEco[ecoType], ".RData"))
}


