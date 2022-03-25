# Script for estimating logarithmised protein concentrations with optical density values from an ELISA test


# Loading necessary packages
library(minpack.lm)

# Set here the technical limit of the ELISA machine (for example 4 when using Infinite F50, Tecan Group AG)
technical_limit <- 4

# Reading in the data 
D <- read.table("TrialData.csv", sep=";", header=T)

# Create a new folder in the working directory for the visualisation of the serial dilution
dir.create("Graphs")


# Definition of function to estimate logarithmised protein concentrations from optical density values
estimate_log.vc <- function(par, y){
  I = as.numeric(par[1])
  S = as.numeric(par[2])
  A = as.numeric(par[3])
  log2(2^I / ((((technical_limit-median_buffercontrols)/((y-median_buffercontrols)))^(1/A)-1)^(1/S)))
}


# Data transformation for each ELISA plate separately 
result <- data.frame()
for(i in 1:11){
  
  D1 <- D[D$Plate==i,]
  
  # Calculate the median of the background noise via buffer controls (bc)
  buffercontrols <- D1[D1$sample=="bc",]
  median_buffercontrols <- median(buffercontrols$OD)
  
  # Extract serial dilution from the data and define the logarithmised protein concentration
  serialdilution <- D1[D1$sample=="serialdilution_0" | D1$sample=="serialdilution_1" | D1$sample=="serialdilution_2" |
                         D1$sample=="serialdilution_3" | D1$sample=="serialdilution_4" | D1$sample=="serialdilution_5" |
                         D1$sample=="serialdilution_6" | D1$sample=="serialdilution_7" | D1$sample=="serialdilution_8" |
                         D1$sample=="serialdilution_9" | D1$sample=="serialdilution_10" | D1$sample=="serialdilution_11",]
  
  serialdilution$log_c <- NA
  serialdilution[serialdilution$sample=="serialdilution_0",]$log_c <- 12
  serialdilution[serialdilution$sample=="serialdilution_1",]$log_c <- 11
  serialdilution[serialdilution$sample=="serialdilution_2",]$log_c <- 10
  serialdilution[serialdilution$sample=="serialdilution_3",]$log_c <- 9
  serialdilution[serialdilution$sample=="serialdilution_4",]$log_c <- 8
  serialdilution[serialdilution$sample=="serialdilution_5",]$log_c <- 7
  serialdilution[serialdilution$sample=="serialdilution_6",]$log_c <- 6
  serialdilution[serialdilution$sample=="serialdilution_7",]$log_c <- 5
  serialdilution[serialdilution$sample=="serialdilution_8",]$log_c <- 4
  serialdilution[serialdilution$sample=="serialdilution_9",]$log_c <- 3
  serialdilution[serialdilution$sample=="serialdilution_10",]$log_c <- 2
  serialdilution[serialdilution$sample=="serialdilution_11",]$log_c <- 1
  
  # Definition of the logistic regression model 
  non.lin.mod <- nlsLM(OD ~ median_buffercontrols+(technical_limit-median_buffercontrols)/(1+(2^I/2^log_c)^S)^A, data=serialdilution,
                       lower=c(I=0, S=0.0001, A=0.0001),
                       upper=c(I=Inf, S=Inf, A=Inf),
                       start=c(I=6, S=1, A=1),
                       weights=1/log_c)
  
  
  # Removing buffer controls and serial dilution from the data before data transformation
  D1 <- D1[D1$sample=="sample",]
  
  # Data transformation for optical density values that are greater than the background noise and smaller than the techincal limit of the machine 
  D1$log.est.vc <- sapply(D1$OD, function(x) {ifelse(median_buffercontrols<x & x<technical_limit, estimate_log.vc(coef(non.lin.mod), x), NA)})
  D1$Status <- NA
  D1$Status <- sapply(D1$OD, function(x) {ifelse(median_buffercontrols<x & x<technical_limit, "error-free transformation", ifelse(median_buffercontrols>=x, "OD too small", "OD at the technical limit")) })
  
  # Combining the results from the data transformation 
  result <- rbind(result, D1)
  
  # Visualisation of the serial dilution (black dots) and the three parameter logistic regression model (blue line)
  sequence <- seq(0,4,by=0.01)
  title <- paste("Plate", i, sep=" ")
  filename <- paste("Graphs/SerialDilution_Plate", i, ".jpg", sep="")
  jpeg(filename, width=1024, height=768)
  plot(sequence ~ estimate_log.vc(coef(non.lin.mod), sequence), type="l", xlim=c(1,12), 
       col="navyblue", xlab="log(protein concentration)", ylab="optical density",
       main=title)
  points(serialdilution$OD ~ serialdilution$log_c, pch=20)
  dev.off()
}


# Exporting the result as csv file
options(scipen=999)
write.table(result, "Result_DataTransformation.csv", sep=";", col.names=T, row.names=F, quote=F)