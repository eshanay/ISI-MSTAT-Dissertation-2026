library(Rcpp)
library(RcppArmadillo)

# Sourcing the C++ file
setwd("C:\\Users\\eshan\\Downloads\\ISI MSTAT references\\MSTAT Dissertation\\R codes")
sourceCpp("wigner_submatrix.cpp")

# Parameters
N <- 1000 #Wigner matrix dimension
p <- 4 #raised to the power p
iterations <- 1000 #number of replications

# Run simulation
#sqrtN <- write.csv(results,"C:\\Users\\eshan\\Downloads\\ISID MSTAT references\\MSTAT Dissertation\\outputs\\Fourth power sum outputs\\sqrtN.csv",row.names=F)
#Nby10 <- write.csv(results,"C:\\Users\\eshan\\Downloads\\ISID MSTAT references\\MSTAT Dissertation\\outputs\\Fourth power sum outputs\\Nby10.csv",row.names=F)

#submatrix dimension k's functional dependence on N is specified in the C code
results <- replicate(iterations,wigner_submatrix_sum_cpp(N, p))
results <- unname(unlist(results))
#standardizing the sum obtained from the entries of the submatrix
vec <- (results - mean(results))/sd(results)
shapiro.test(vec)$p.val

par(mar=c(4.5,4.5,1,1))
hist(vec, 
     main = "",
     prob=T,
     ylim=c(0,dnorm(0)),
     xlab = expression("k ~ sqrt(N)"),
     col = "lightblue",
     border = "black")
#overlaying the true standard normal density for graphical comparison
curve(dnorm(x),add=T,col="red",lwd=2)
