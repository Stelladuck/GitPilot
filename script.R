#### SECTION 0. Data preparation ####
# Example data of step function ￩
x       <- seq(from = 0, to = 4, by = 0.05)
y0      <- c(1, -1, 1, 2)
func    <- stepfun(1:3, y0, f = 0)
dy      <- runif(length(x), -0.1, 0.1)
y       <- func(x) + dy

#### SECTION 1. Objective function ####
fmodel  <- function(x, w) { w[1] + w[2]*x + w[3]*x**2 + w[4]*x**3 + w[5]*x**4 + w[6]*x**5}
fobj    <- function(w) { y_pred = fmodel(x, w); sqrt(sum((y - y_pred)**2) / length(y)) }

#### SECTION 2. Population ####
bound   <- c(-10, 10)
popsize <- 20
normed.population <- matrix(runif(120, 0, 1), nrow = popsize) # normalized
denorm.population <- min(bound) + normed.population*(max(bound) - min(bound))

#### SECTION 3. Evaluation ####
fitness <- rep(0, len = popsize)
for (i in seq(from = 1, to = popsize)){ w <- denorm.population[i,]; fitness[i] <- fobj(w) }
bestfit <- max(fitness)
bestarg <- normed.population[which.max(fitness),]

#### Differential evolution function ####
DE <- function(fobj, bound, mut = 0.8, crossp = 0.7, popsize = 20, iter = 1000) {
    for (i in seq(from = 1, to = iter)) {
        for (j in seq(from = 1, to = popsize)) {
            #### SECTION 4. Mutation ####
            choice  <- sample(seq(1:popsize)[c((1:popsize)!=j)], 3, replace = F)
            a <- normed.population[choice[1],]
            b <- normed.population[choice[2],]
            c <- normed.population[choice[3],]
            mutant  <- pmax(0, pmin(a + mut*(b-c), 1)) # normalized mutant
            cross   <- runif(ncol(normed.population), 0, 1) < crossp
            normed.trial    <- as.numeric(!cross)*normed.population[j,] + as.numeric(cross)*mutant
            denorm.trial    <- min(bound) + normed.trial*(max(bound) - min(bound))
            #### SECTION 5. Recombination ####
            f <- fobj(denorm.trial)
            if (f < fitness[j]) {
                fitness[j] <- f
                normed.population[j,] <- normed.trial
                denorm.population[j,] <- denorm.trial
            }
            if (f < bestfit) {
                bestfit <- f
                bestarg <- denorm.trial
            }
        }
    }
    return(c(bestarg, 1/bestfit))
}

#### Obtain best approximation ####
sol <- DE(fobj, bound)

#### Graphics ####
plot(x,y)
lines(x,fmodel(x,sol[1:length(sol)-1]))

