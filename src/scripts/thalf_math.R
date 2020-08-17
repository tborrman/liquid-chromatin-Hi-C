
exp_decay <- function(minutes, a, b, c) {
  # Sigmoid function with adjustable parameters
  return(a - (b*exp(-c*minutes)))
}

half_life_equation <- function(LOS_half, a, b, c) {
  # Get half-life using exponential decay function 
  # solved for minutes at half-Loss of structure
  hl <- -(log((a-LOS_half)/b, base=exp(1)) / c)
  return(hl)
}

get_half_life <- function(LOS, minutes, LOS_half) {
  # Fit LOS data with exponential decay curve and return 
  # half-life and sum of squared residuals
  fitM <- try(nls(LOS ~ exp_decay(minutes, a, b, c), start = list(a=0.75, b=0.5, c=0.01), 
                  trace=TRUE, control = nls.control(warnOnly = TRUE)))
  if(class(fitM) == "try-error") {
    return(c(NA, NA, NA))
  }
  else{
    params <- coef(fitM)
    ssr <- sum(as.numeric(residuals(fitM))^2)
    half_life <- half_life_equation(LOS_half, params["a"], params["b"], params["c"])
    tau <- 1/params["c"]
    return(c(half_life, ssr, tau))
  }
}

get_LOS_half <- function(row){
  LOS_start <- row[1]
  LOS_end <- row[length(row)]
  LOS_half <- LOS_start + ((LOS_end - LOS_start)/2)
  return(LOS_half)
}

min <- c(5, 60, 120, 180, 240, 960)

d <- c(70, 60, 30, 15,10, 5)

LOS1 <- (80-d)/ 80
LOS2 <- (40-d)/ 40
LOS3 <- (20-d)/ 20

LOShalf1 <- get_LOS_half(LOS1)
LOShalf2 <- get_LOS_half(LOS2)
LOShalf3 <- get_LOS_half(LOS3)



plot(min, LOS1, type="l", ylim=c(-2, 5))
lines(min, LOS2, type="l", ylim=c(-2, 5), col="blue")
lines(min, LOS3, type="l", ylim=c(-2, 5), col="red")

r1 <- get_half_life(LOS1, min, LOShalf1)
r2 <- get_half_life(LOS2, min, LOShalf2)
r3 <- get_half_life(LOS3, min, LOShalf3)

abline(v=r1[1], lty=2)

# Yep half-life is invariant to changing denominator of LOS

min <- c(5, 60, 120, 180, 240, 960)


d <- c(70, 60, 30, 15,10, 5)

LOS1 <- (80-d)/ 80
LOS2 <- (40-d)/ 80
LOS3 <- (20-d)/ 80

LOShalf1 <- get_LOS_half(LOS1)
LOShalf2 <- get_LOS_half(LOS2)
LOShalf3 <- get_LOS_half(LOS3)



plot(min, LOS1, type="l", ylim=c(-2, 5))
lines(min, LOS2, type="l", ylim=c(-2, 5), col="blue")
lines(min, LOS3, type="l", ylim=c(-2, 5), col="red")

r1 <- get_half_life(LOS1, min, LOShalf1)
r2 <- get_half_life(LOS2, min, LOShalf2)
r3 <- get_half_life(LOS3, min, LOShalf3)

abline(v=r1[1], lty=2)

# Yep half-life is invariant to changing denominator of controlcis in numerator of LOS

min <- c(5, 60, 120, 180, 240, 960)

d <- c(70, 60, 30, 15,10, 5)

LOS1 <- (80-d)/ 80
LOS2 <- (40-d)/ 100
LOS3 <- (20-d)/ 20

LOShalf1 <- get_LOS_half(LOS1)
LOShalf2 <- get_LOS_half(LOS2)
LOShalf3 <- get_LOS_half(LOS3)



plot(min, LOS1, type="l", ylim=c(-2, 5))
lines(min, LOS2, type="l", ylim=c(-2, 5), col="blue")
lines(min, LOS3, type="l", ylim=c(-2, 5), col="red")

r1 <- get_half_life(LOS1, min, LOShalf1)
r2 <- get_half_life(LOS2, min, LOShalf2)
r3 <- get_half_life(LOS3, min, LOShalf3)

abline(v=r1[1], lty=2)

# Yep half-life is invariant to changing control cis percent of LOS


