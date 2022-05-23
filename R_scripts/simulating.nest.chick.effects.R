#################################################
### SIMULATING JQ'S PROBLEM:
#################################################

# load lme4 package first:
if (!require("lme4")) install.packages("lme4")
library(lme4)

# First, input parameter values for the random effects:
nest.var= 1 # between-nest variance component 
chick.var= 1 # between-chicks (within nests) variance component  ## SET THIS TO EITHER 1 OR 0
resid.var<- 1 # residual variance component

### Now set up the factors for the model:
n.nests = 100 # no. nests
n.chicks = 10 # no. chicks per nest
n.repeats = 5 # no. repeated measures per chick

# Make up unique identifiers for nests:
nestID<- rep(1:n.nests,n.chicks*n.repeats)
nestID<- sort(nestID)
chickID<- 1:(n.chicks*n.nests)
chickID<- sort(rep(chickID, n.repeats))
chick.measure<- rep(1:n.repeats, n.chicks*n.nests)

# Put the above into a data frame:
sim.dat<- data.frame(nestID, chickID, chick.measure)

head(sim.dat)

# Now simulate a common nest effect:
nest.effect<- rnorm(n.nests, 0, nest.var)
temp.nest.IDs<- 1:n.nests

sim.dat$nest.effect <- nest.effect[match(sim.dat$nestID, temp.nest.IDs)]

# Now simulate a chick effect:
chick.effect<- rnorm(n.chicks*n.nests, 0, chick.var)
temp.chick.IDs<- 1:(n.chicks*n.nests)

sim.dat$chick.effect <- chick.effect[match(sim.dat$chickID, temp.chick.IDs)]

head(sim.dat)

# Now finally add random noise to each row in the data frame, corresponding to residual variance around the chick and nest effects:

n<- dim(sim.dat)[1]
sim.dat$resid.effect <- rnorm(n, 0, resid.var)

# finally, create the phenotype as a linear sum of each effect:
sim.dat$phen<- sim.dat$nest.effect + sim.dat$chick.effect + sim.dat$resid.effect

# have a look at the distribution of effects:
par(mfrow=c(2,2), mar=c(4,4,0,0))
hist(sim.dat$phen, main="", xlim=c(-8,8))
hist(sim.dat$nest.effect, main="", xlim=c(-8,8))
hist(sim.dat$chick.effect, main="", xlim=c(-8,8))
hist(sim.dat$resid.effect, main="", xlim=c(-8,8))

#  Now fit the mixed models:
sim.dat$nestID <- factor(sim.dat$nestID)
sim.dat$chickID <- factor(sim.dat$chickID)

# Model 1: just a nest effect:
mod1<- lmer(phen ~ 1 + (1|nestID), data=sim.dat)
summary(mod1)

# Model 2: just a chick effect:
mod2<- lmer(phen ~ 1 + (1|chickID), data=sim.dat)
summary(mod2)

# Model 3: both effects together, with chick nested within nest effect:
mod3<- lmer(phen ~ 1 + (1|nestID/chickID), data=sim.dat)
summary(mod3)



