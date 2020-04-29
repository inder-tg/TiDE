
# how to use the occurence-amount model to simulate a precipitation process

# --- auxiliary functions

# --- functions to make simulation work
source(paste0(getwd(), "/simsCode/precipitationGEN.R"))
source(paste0(getwd(), "/simsCode/parameters.R"))

getPixel <- function(m){
  flag <- F
  
  while(flag == F){
    i <- sample(nrow(m),1)
    t <- m[i,]
    if(sum(is.na(t)) == 0){
      flag <- T
    }
  }
  t  
}

# --- brief explanation

# main function precipitationGenerator(data, yesterdayWeather, gammaParam, n) # found version for gamma distribution
# 
# data is going to be a randomly chosen pixel from some of the ecoregions of interest
# 
# yesterdayWeather is a Bernoulli trial (ifelse(rbinom(1,1,.5) == 1, "wet", "dry"))
# 
# gammaParam is going to be a 12 x 2 matrix with shape and scale parameters of gamma density.
# The values of these parameters were estimated from precipitation time series extracted from 
# each ecoregion of interest. See getParmsGammaSimsEcoRegIV.R for more details
# 
# n is sample size

# --- example

ecoLabels <- c("Madrense", "interiorPlains", "Mezquital", "plateauPlains")

params <- LoadToEnvironment(paste0(getwd(), "/data/RData/paramsRAINFALL_", ecoLabels[1], 
                                   ".RData"))
pixels <- LoadToEnvironment(paste0(getwd(), "/data/RData/pixels_PREC_", ecoLabels[1], 
                                   ".RData"))

gammaParams <- params$monthlyGammaParamTemp

pixelForSims <- pixels$pixelsTemp

pixel <- getPixel(pixelForSims)

typeWeather <- ifelse(rbinom(1,1,.5) == 1, "wet", "dry")

t2 <- precipitationGenerator(data = pixel,
                             yesterdayWeather = typeWeather,
                             gammaParam = gammaParams,
                             n = length(pixel))

plot(t2$fullSim, type = "h")

# --- overlapping of precipitation process on impulse signal

tau <- 10
Tau <- ceiling( tau * 366/100 )
stepTest1 <- createSignalPeak(0, 0, ceiling(3*366/10), ceiling(6*366/10), 1)
stepTest2 <- createSignalPeak(0, 0, ceiling(3*366/10) + Tau, ceiling(6*366/10) + Tau, 1)

n <- 366

ind <- seq(0, n, length = n)

signal1 <- stepTest1(ind)
signal2 <- stepTest2(ind)


copySignal <- signal1

int1 <- 1:ceiling(3*366/10) # 20
int2 <- (ceiling(3*366/10)+1):ceiling(6*366/10)
int3 <- (ceiling(6*366/10)+1):n

ind1 <- int1[t2$fullSim != 0]
ind2 <- int2[t2$fullSim == 0]
ind3 <- int3[t2$fullSim != 0]
indices <- c(ind1[!is.na(ind1)], ind2[!is.na(ind2)], ind3[!is.na(ind3)])

copySignal[indices] <- copySignal[indices] + t2$fullSim[indices]


data1 <- copySignal

b <- 0.0075
noise <- rnorm(n, sd = b)

data2 <- signal2 + noise

plot(data1, col = "blue", type = 'h')
par(new=T)
plot(data2, col = "red")

#----------------------------------------------------------------------

