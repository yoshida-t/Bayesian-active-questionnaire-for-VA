## ---------------------------------------------------------  ##
## Prepare full PHMRC data
## ---------------------------------------------------------  ##
library(openVA)
PHMRC <- read.csv(getPHMRC_url("adult"))

# turn into binary data
binarydata <- ConvertData.phmrc(input = PHMRC, input.test = PHMRC, cause = "gs_text34")
causes <- as.character(unique(PHMRC$gs_text34))

tmp <- binarydata$output[, -c(1, 2)]
X0 <- matrix(0, dim(tmp)[1], dim(tmp)[2])
X0[tmp == "Y"] <- 1
X0[tmp == "."] <- NA
head(X0)
Y0 <- match(as.character(binarydata$output[,2]), causes)

sympnames <- colnames(binarydata$output)[-c(1, 2)]
# load better symptom names
symplist <- read.csv("data/PHMRC_symptoms.csv")
sympnames_new <- symplist[match(sympnames, symplist[,1]), 3]

colnames(X0) <- sympnames_new 

library(tidyverse)
library(caret)

C = length(unique(Y0))
n = dim(X0)[1]
p = dim(X0)[2]

set.seed(0)
folds = createFolds(Y0, k = 10, returnTrain = TRUE)

fit_inter <- fit_ins <- out.compare <- NULL
for (i in c(1:10)) {
	is.train <- which(1:dim(PHMRC)[1] %in% folds[[i]])
	train <- PHMRC[is.train, ]
	test <- PHMRC[-is.train, ]
	dim(test)
	dim(train)

	## ----phmrc-inter, messages=FALSE, results='hide'----------------------------------------
	fit_inter[[i]] <- codeVA(data = test, data.type = "PHMRC", model = "InterVA", 
	                     data.train = train, causes.train = "gs_text34", 
	                     phmrc.type = "adult", convert.type = "fixed")

	## ----phmrc-ins, messages=FALSE, results='hide'------------------------------------------
	fit_ins[[i]] <- codeVA(data = test, data.type = "PHMRC", model = "InSilicoVA",
	                    data.train = train, causes.train = "gs_text34", 
	                    phmrc.type = "adult", convert.type = "fixed", 
	                    Nsim=10000, auto.length = FALSE)

	pred.int <- getTopCOD(fit_inter[[i]], interVA.rule = FALSE)
	pred.int$truth <- test$gs_text34
	pred.int$type = "InterVA"
	pred.ins <- getTopCOD(fit_ins[[i]])
	pred.ins$truth <- test$gs_text34
	pred.ins$type = "InSilicoVA"
	tmp <- rbind(pred.int, pred.ins)
	tmp$fold <- i
	out.compare <- rbind(out.compare, tmp)
	save(out.compare, file = "output/PHMRC_10fold_compare.RData")
}

