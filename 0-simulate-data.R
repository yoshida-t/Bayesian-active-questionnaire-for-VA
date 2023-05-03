###########################
#
# Generate simulated data
#
###########################


source("functions.R")

###########################

n = 4000
train_size = 2000

N = 1
p = 50
C = 10

N2 = 100
st2 = p

###########################

set.seed(1)
trainid = sample(1:n, size=train_size, replace=FALSE)
smlted_misspecified = simulateXY(n, p, C, trainid=trainid, generate_option='Misspecified (K=3)')
smlted_correct = simulateXY(n, p, C, trainid=trainid, generate_option='Different theta')

save(n, train_size, st, N, p, C, N2, st2, trainid, smlted_misspecified, file = 'input_data/misspecified.RData')
save(n, train_size, st, N, p, C, N2, st2, trainid, smlted_correct, file = 'input_data/correct.RData')
