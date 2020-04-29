
# lengthMonthsLeapYear required

LoadToEnvironment <- function(RData, env = new.env()) {
  load(RData, env)
  return(env) 
}

getEmpTransProb <- function(x){
  
  y <- x
  y[y != 0] <- 1
  
  matTransProb <- matrix(0, 2, 2)
  
  yesterday <- y[1:(length(y)-1)]
  today <- y[2:length(y)]
  
  for(i in 1:(length(y)-1)){
    if(yesterday[i] == 0 & yesterday[i] == today[i] ){
      matTransProb[1,1] <- matTransProb[1,1] + 1
    }
    
    if(yesterday[i] < today[i]){
      matTransProb[1,2] <- matTransProb[1,2] + 1
    }
    
    if(yesterday[i] > today[i] ){
      matTransProb[2,1] <- matTransProb[2,1] + 1
    }
    
    if(yesterday[i] == 1 & yesterday[i] == today[i] ){
      matTransProb[2,2] <- matTransProb[2,2] + 1
    }
  }
  
matTransProb
}

getMonthlyTransProb <- function(x, lengthMonths){
  
  cummulatedNumDays <- cumsum(lengthMonths)
  
  arrayTransProb <- array(0, dim = c(nrow = 2, ncol = 2, nmonths = 12))
  
  arrayTransProb[,,1] <- getEmpTransProb( x[ 1:lengthMonths[1] ] )
  
  arrayTransProb[,,2:12] <- sapply(2:12, function(s) 
    getEmpTransProb( x[ cummulatedNumDays[s-1] + 1:lengthMonths[s] ] ) )
  
arrayTransProb  
}

getMonth <- function(k, lengthMonths){
  
  cumsumNumDays <- cumsum(lengthMonths) 
  
  if(k > cumsumNumDays[length(cumsumNumDays)])
    stop("k is greater than length of time series")
  
  month <- numeric(1)
  
  monthsLength <- c(0, cumsumNumDays)
  
  month <- sum( monthsLength - k < 0 )

month
}

precipitationGenerator <- function(data, yesterdayWeather = c("dry", "wet"), gammaParam, n){
  
  # data = precipitation_17Years[[16]]
  # yesterdayWeather = "dry"
  # gammaParam = monthlyParam
  
  weatherYesterday <- match.arg(yesterdayWeather)

  if(missing(n)){
    lengthGen <- length(data)
  } else {
    if(n > length(data)){
      stop("n must be less or equal than length(data)")
    } else {
      lengthGen <- n
    }
  }
  
  lengthMonths <- readRDS(paste0(getwd(), "/simsCode/lengthMonthsLeapYear"))
  
  if(length(data) == 365){
    lengthMonths[2] <- 28
  }
  
  transProb <- getMonthlyTransProb(x = data, lengthMonths = lengthMonths)
  
  monthsProb <- array(0, dim = dim(transProb))
  
  for(i in 1:12){
    monthsProb[1,1,i] <- transProb[1,1,i] / (transProb[1,1,i] + transProb[1,2,i]) 
    monthsProb[1,2,i] <- 1 - monthsProb[1,1,i]
    monthsProb[2,1,i] <- transProb[2,1,i] / (transProb[2,1,i] + transProb[2,2,i]) 
    monthsProb[2,2,i] <- 1 - monthsProb[2,1,i]
  }
  
  simPrecipitation <- numeric(length(data))
  for(k in 1:length(data)){
    
    # k <- 1
    
    month <- getMonth(k = k, lengthMonths = lengthMonths)
    
    p00 <- monthsProb[1,1,month]
    p10 <- monthsProb[2,1,month]
    
    feasibility <- sum( is.nan(c(p00, p10)) )
    
    if(feasibility == 0){
      chancesRainToday <- runif(1)
      # tempParam <- gammaParam[[]]
      paramGamma <- gammaParam[month,]
      
      if(weatherYesterday == "dry"){
        if(chancesRainToday < p00){
          simPrecipitation[k] <- 0 
        } else {
          q <- runif(1)
          simPrecipitation[k] <- qgamma(p = q, shape = paramGamma[1], scale = paramGamma[2])
          weatherYesterday <- "wet"
        }
      }
      
      if(weatherYesterday == "wet"){
        if(chancesRainToday < p10){
          simPrecipitation[k] <- 0 
          weatherYesterday <- "dry"
        } else {
          q <- runif(1)
          simPrecipitation[k] <- qgamma(p = q, shape = paramGamma[1], scale = paramGamma[2])
        }
      }
    }
  }
  
  if(lengthGen != length(data)){
    sample <- sample(simPrecipitation, size = lengthGen)
  } else {
    sample <- simPrecipitation
  }
  
list(fullSim = simPrecipitation, sampleSim = sample)
}

# getSignal <- function(n, tau){
#   # stepF1Sym <- createSignalPeak(0, 0, 3, 6, 1)
#   # stepF2Sym <- createSignalPeak(0, 0, 4, 7, 1)
#   
#   tau <- 10
#   Tau <- tau * 365/100
#   stepTest1 <- createSignalPeak(0, 0, 3*365/10, 6*365/10, 1)
#   stepTest2 <- createSignalPeak(0, 0, 3*365/10 + Tau, 6*365/10 + Tau, 1)
#   
#   n <- 365
#   
#   ind <- seq(0, n, length = n)
#   
#   signal1 <- stepTest1(ind)
#   signal2 <- stepTest2(ind)
#   
#   
# }





