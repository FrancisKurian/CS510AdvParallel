## Diffusion function
diffuse <- function(dots, D, delta.t){
  dist <- sqrt(4*D*delta.t)
  randos1 <- rnorm(n = nrow(dots), mean = dist, sd = 0.25*dist)
  randos2 <- runif(n = nrow(dots), min = 0, max = 2*pi)
  dots[, 1] <- dots[, 1] + randos1*cos(randos2)
  dots[, 2] <- dots[, 2] + randos1*sin(randos2)
  return(dots)
}