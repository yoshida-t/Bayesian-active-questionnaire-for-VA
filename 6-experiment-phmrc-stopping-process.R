######################################################
#
# Post process the PHMRC stopping rule examples
#
######################################################


library(ggplot2)
source("functions.R")

K = 10
output <- NULL
for(fold in 1:K){
  file = paste0('output/fold', fold , '_PPmeans_mean_list.RData')
  load(file)
  file = paste0('output/fold', fold , '_PPmeans_post_0.7_list.RData')
  load(file)
  file = paste0('output/fold', fold , '_PPmeans_post_0.5_list.RData')
  load(file)
  PPmeans_list = list(PPmeans_mean_2, PPmeans_post_0.7_2, PPmeans_post_0.5_2)
  filename = paste0('output/stoppping_CV_table_fold', fold, '.csv')
  tmp <- createResultsTable(PPmeans_list, file=filename, type_names = c("mean", "0.7", "0.5"))
  output <- rbind(output, tmp)
}


library(openVA)
library(xtable)
PHMRC <- read.csv(getPHMRC_url("adult"))
tab <- sort(table(PHMRC$gs_text34), dec = TRUE)
xtable(tab)

binarydata <- ConvertData.phmrc(input = PHMRC, input.test = PHMRC, cause = "gs_text34")
causes <- as.character(unique(PHMRC$gs_text34))
causes.ordered <- rev(names(sort(table(PHMRC$gs_text34))))
causes[causes == "Other Non-communicable Diseases"] <- "Other NCD"
causes.ordered[causes.ordered == "Other Non-communicable Diseases"] <- "Other NCD"

data <- output
colnames(data)[2:7] = c("p1", "p2", "ID", "length", "true", "pred")
data$true <- causes[as.numeric(data$true)]
data$pred <- causes[as.numeric(data$pred)]
head(data)
data$name <- paste(data$p1, data$p2, data$Type)
# table(data$name)
data$correct <- data$true == data$pred
data$length <- as.numeric(data$length)

# data <- subset(data, p1 >= 0.8)
data$d <- -1
data$d[grepl("d= \\$ 0", data$p2)] <- 0
data$d[grepl("d= \\$ 0.75", data$p2)] <- 0.75
data$d[grepl("d= \\$ 0.5", data$p2)] <- 0.5
data$d[grepl("d= \\$ 0.25", data$p2)] <- 0.25

acc <- aggregate(correct ~ p1 + d + Type, data = data, FUN = mean)
len <- aggregate(length ~ p1 + d + Type, data = data, FUN = quantile, .5)
summary <- merge(acc, len)
g1 <- ggplot(subset(summary, d %in% c(0, 0.25, 0.5, 0.75))) + aes(x = length, y = correct, color = Type, shape = factor(d)) + geom_point(alpha = 0.5, size = 3)  + facet_wrap(~p1) + scale_color_brewer("Stopping rule", palette = "Set1") + theme_bw() + xlab("Median number of questions asked") + ylab("Classification accuracy") + scale_shape_discrete("Choice of d") + ylim(c(0.3, 0.5))


acc <- aggregate(correct ~ p1 + d + Type + true, data = data, FUN = mean)
len <- aggregate(length ~ p1 + d + Type + true, data = data, FUN = quantile, .5)
summary <- merge(acc, len)
summary$true <- factor(summary$true, levels = causes.ordered)
summary$Type[summary$Type == "mean"] <- "Point Estimate" 
summary$Type[summary$Type == "0.7"] <- "Predictive Prob > 0.7"
summary$Type[summary$Type == "0.5"] <- "Predictive Prob > 0.5"
g2 <- ggplot(subset(summary, d %in% c(0.5) & p1 < 0.95 & p1 > 0.6)) + 
      aes(x = length, y = correct, color = Type, shape = factor(p1)) + 
      geom_point(alpha = 0.5, size = 3) + 
      geom_line(aes(group = Type)) +  
      scale_color_brewer("Stopping Criterion", palette = "Set2") + 
      theme_bw() + 
      xlab("Median Number of Questions Asked") + 
      ylab("Probability of Correct Classification") + 
      scale_shape_discrete("Threshold (p_1st)") + 
      facet_wrap(~true ) + 
      theme(legend.position = "bottom") +
      ylim(0, 1)
ggsave(g2, file = "fig/6-experiment-phmrc-stop.pdf", width = 11, height = 12)


# add comparison
load("output/PHMRC_10fold_compare.RData")
out.compare$cause1[out.compare$cause1 == "Other Non-communicable Diseases"] <- "Other NCD"
out.compare$truth[out.compare$truth == "Other Non-communicable Diseases"] <- "Other NCD"
out.compare$correct <- out.compare$cause1 == out.compare$truth
summary.compare <- aggregate(correct ~ truth + type, data =  out.compare, FUN = mean)  
base1 <- subset(summary.compare, type == "InterVA")[, c(1, 3)]
base2 <- subset(summary.compare, type == "InSilicoVA")[, c(1, 3)]
colnames(base1) <- c("true", "InterVA")
colnames(base2) <- c("true", "InSilicoVA")
summary2 <- merge(summary, base1, all.x = TRUE)
summary2 <- merge(summary2, base2, all.x = TRUE)

label <- data.frame(sort(table(PHMRC$gs_text34)))
label[,1] <- as.character(label[,1])
label[label[,1] == "Other Non-communicable Diseases", 1] <- "Other NCD"
label$new <- paste0(label[,1], "\n", "(n = ", label[,2], ")")
summary2$label <- label[match(summary2$true, label[,1]), "new"]
label.ordered <- rev(label$new[order(label$Freq)])
summary2$label <- factor(summary2$label , levels = label.ordered)
g3 <- ggplot(subset(summary2, d %in% c(0.5) & p1 < 0.95 & p1 > 0.6)) + 
      aes(x = length, y = correct, color = Type, shape = factor(p1)) + 
      geom_hline(aes(yintercept = InterVA), color = "black", alpha = 0.4, linetype = "dashed") + 
      geom_hline(aes(yintercept = InSilicoVA), color = "red", alpha = 0.4, linetype = "dashed") + 
      geom_point(alpha = 0.5, size = 3) + 
      geom_line(aes(group = Type)) +  
      scale_color_brewer("Stopping Criterion", palette = "Set2") + 
      theme_bw() + 
      xlab("Median Number of Questions Asked") + 
      ylab("Probability of Correct Classification") + 
      scale_shape_discrete("Threshold (p_1st)") + 
      facet_wrap(~label ) + 
      theme(legend.position = "bottom") +
      ylim(0, 1)
ggsave(g3, file = "fig/6-experiment-phmrc-stop-with-comparison.pdf", width = 12, height = 12)

# for appendix
C = length(table(data$true))
nsplit = 5
causes.ordered.split = split(causes.ordered, ceiling(seq_along(causes.ordered)/ceiling(C/nsplit)))
g2_list = vector('list', nsplit)
for(ni in 1:nsplit){
  g2_list[[ni]] = ggplot(subset(summary, d %in% c(d_cand) & p1 < 0.95 & p1 > 0.6 & true %in% causes.ordered.split[[ni]])) +
    aes(x = length, y = correct, color = Type, shape = factor(p1)) + 
    geom_point(alpha = 0.5, size = 3) + 
    geom_line(aes(group = Type)) +  
    scale_color_brewer("Stopping Criterion", palette = "Set2") + 
    theme_bw() + 
    xlab("Median Number of Questions Asked") + 
    ylab("Probability of Correct Classification") + 
    scale_shape_discrete(bquote('Threshold'~p[1*st])) +
    facet_grid(true~d, labeller=labeller(.cols=label_both, .rows=label_value)) +
    theme(legend.position = "bottom") +
    ylim(0, 1)
}
ggsave(g2_list[[1]], file = "fig/6-experiment-phmrc-stop-full1.pdf", width = 11, height = 12)
ggsave(g2_list[[2]], file = "fig/6-experiment-phmrc-stop-full2.pdf", width = 11, height = 12)
ggsave(g2_list[[3]], file = "fig/6-experiment-phmrc-stop-full3.pdf", width = 11, height = 12)
ggsave(g2_list[[4]], file = "fig/6-experiment-phmrc-stop-full4.pdf", width = 11, height = 12)
ggsave(g2_list[[5]], file = "fig/6-experiment-phmrc-stop-full5.pdf", width = 11, height = 12)
