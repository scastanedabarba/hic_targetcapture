library(tidyr)
library(MuMIn)
library(jtools)
library(MASS)

results = read.csv('read_simulations2.csv')

#convert library to factor and depth to numerical
results$Library = as.factor(results$Library)
results$Depth = as.numeric(results$Depth)

results$plated = c(rep(c(1.18e7),25), rep(c(1.37e6),25), rep(c(1.94e5),25), rep(c(1.44e4),25), rep(c(1.32e3),25))
#results$plated_optiprep = c(rep(c(6.83e5),25), rep(c(5.55e4),25), rep(c(2.66e4),25), rep(c(1.05e3),25), rep(c(1.61e2),25))


#########model for relating Hi-C reads to bacterial count##############################
model = glm(plated ~ log(pB10_pputida+.0000001) + log(Depth), 
              data=results, family='poisson')
summary(model)
r.squaredLR(model)
summ(model)


#make dataframe with predictions
pB10_pputida = sort(c(rep(c(1:250), 3)))
depth = c(rep(c(90000000),750))
desired_reads = c(rep(c(0, 10, 100),250))
predictions = data.frame(hic_reads = pB10_pputida,
                         Depth = depth,
                         desired_reads = desired_reads)
predictions['prediction'] = list(predict(model, predictions, type='response'))

#add probability of seeing more than x hi-c reads
predictions['probability'] = round(100-(ppois(predictions$desired_reads, predictions$hic_reads, lower.tail=TRUE)*100), digits = 2)
