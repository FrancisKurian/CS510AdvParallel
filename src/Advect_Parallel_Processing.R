# Advect Parallel Processing
source("src/advect.R")
source("src/diffuse.R")

library(foreach)
library(doParallel)
library(parallel)
library(profvis)

# find no. of cores and register
cores <- detectCores()
cores
cluster <- makePSOCKcluster(cores)
registerDoParallel(cores=cores)

#### Defines parameters ####
D <- 1e-8        # Diffusion coefficient in m^2/s
delta.t <- 1e-3  # Time step size in s
end.time <- .10    # End simulation time in s
start.x <- 5.6
start.y <- -0.97

# Creates points to follow
num.dots <- 1000
dotsx <- rep(start.x, num.dots) # replicates start.x 1000 times
dotsy <- rep(start.y, num.dots) # replicates start.y 1000 times
dots.start <- matrix(c(dotsx, dotsy), num.dots, 2) # creates a matrix with num.dots rows and 2 columns filled with x & y

# Reads in position (x,y) and velocity (Ux,Uy) data
x <- as.matrix(read.table("data/x.csv", header = FALSE, sep = ","))
y <- as.matrix(read.table("data/y.csv", header = FALSE, sep = ","))
Ux <- as.matrix(read.table("data/Ux.csv", header = FALSE, sep = ","))
Uy <- as.matrix(read.table("data/Uy.csv", header = FALSE, sep = ","))
Ux[is.na(Ux)] <- 0
Uy[is.na(Uy)] <- 0

t <- 0
dots <- dots.start  # re-assigning the start matrix to dots

# define a sequence to use inside foreach function

ts=seq(delta.t,end.time,delta.t)

# Profile the original simulation using profvis

p1 <- profvis({
  trial1 <-  while(t <= end.time){
    dots <- advect(dots, x, y, Ux, Uy, 0.5*delta.t)
    dots <- diffuse(dots, D, delta.t)
    dots <- advect(dots, x, y, Ux, Uy, 0.5*delta.t)
    t <- t + delta.t
    print(paste("t =", t))
  }
  dots.end1 <- dots
  
  
})

#profile of modified parallel process

p2 <- profvis({
  dots.end2 <- foreach(t =ts,.combine ='cbind'  )%dopar%{
    dots <- advect(dots, x, y, Ux, Uy, 0.5*delta.t)
    dots <- diffuse(dots, D, delta.t)
    dots <- advect(dots, x, y, Ux, Uy, 0.5*delta.t)
    if(t==tail(ts, n=1)){
      return(dots)  
    }
  }
  
})

# dots.end1 
# dots.end2 

plot(dots.end1[, 1], dots.end1[, 2], col = "black", pch = 19)
points(dots.start[, 1], dots.start[, 2], col = "red", pch = 19)

plot(dots.end2[, 1], dots.end2[, 2], col = "black", pch = 19)
points(dots.start[, 1], dots.start[, 2], col = "red", pch = 19)



htmlwidgets::saveWidget(p1, "simulation_results.html")
htmlwidgets::saveWidget(p2, "simulation_results_parallel.html")


