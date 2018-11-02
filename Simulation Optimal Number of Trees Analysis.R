

res.mmr <- res.const <- res.mmr2 <- res.const2 <- list()
for(i in 1:100){
  load(paste(getwd(), "/Gradient Boosting/Diff10/NTrees/simgbpen_", i, ".RData", sep = ""))
  res.mmr[[i]] <- res.gb$XGB.mmr$eval_best
  res.const[[i]] <- res.gb$XGB.const$eval_best
  res.mmr2[[i]] <- res.gb2$XGB.mmr$eval_best
  res.const2[[i]] <- res.gb2$XGB.const$eval_best
}

save(res.mmr, res.const, res.mmr2, res.const2, file = "sim_ntrees_perf.RData")

best.mmr <- setNames(data.frame(matrix(unlist(lapply(res.mmr, colMeans)), 100, 1, byrow = TRUE)),
                     "EO")
best.const <- setNames(data.frame(matrix(unlist(lapply(res.const, colMeans)), 100, 1, byrow = TRUE)), "EO")
best.mmr2 <- setNames(data.frame(matrix(unlist(lapply(res.mmr2, colMeans)), 100, 2, byrow = TRUE)),
                      c("AUC", "BS"))
best.const2 <- setNames(data.frame(matrix(unlist(lapply(res.const2, colMeans)), 100, 2, byrow = TRUE)),
                        c("AUC", "BS"))

summary(best.mmr)
summary(best.const)
summary(best.mmr2)
summary(best.const2)

library(ggplot2)
ggplot(data.frame(AUC = c(best.mmr2$AUC, best.const2$AUC), Type = rep(c("MMR", "Const"), each = 100)),
       aes(AUC)) +
  geom_histogram(aes(fill = as.factor(Type), color = as.factor(Type)),
                 binwidth = 5, alpha = 0.3,
                 linetype = "solid", position = "identity") +
  scale_fill_manual(values = c("red", "blue"),
                    name = "Initial Prediction", labels = c("Constant", "MMRpro")) +
  scale_color_manual(values = c("red", "blue"), guide = FALSE) +
  labs(x = "Optimal Number of Trees", y = "Count") +
  geom_vline(xintercept = mean(best.mmr2$AUC), col = "blue", size = 1.2) +
  geom_vline(xintercept = mean(best.const2$AUC), col = "red", size = 1.2)


plot(1, type = "n", xlim = c(0, 300), ylim = c(0.9, 1.1), ylab = "EO", xlab = "Iteration",
     main = "Starting with MMRpro")
for(i in 1:20){
  lines(res.gb$XGB.mmr$eval[[i]]$test_EO, lwd = 0.5)
  abline(v = which.min(abs(res.gb$XGB.mmr$eval[[i]]$test_EO - 1)), col = "red", lwd = 0.5)
}
abline(h = 1)

plot(1, type = "n", xlim = c(0, 300), ylim = c(0.3, 1.1), ylab = "EO", xlab = "Iteration",
     main = "Starting with constant")
for(i in 1:20){
  lines(res.gb$XGB.const$eval[[i]]$test_EO, lwd = 0.5)
  abline(v = which.min(abs(res.gb$XGB.const$eval[[i]]$test_EO - 1)), col = "red", lwd = 0.5)
}
abline(h = 1)



plot(1, type = "n", xlim = c(0, 300), ylim = c(0.9, 0.95), ylab = "AUC", xlab = "Iteration",
     main = "Starting with MMRpro")
for(i in 1:100){
  lines(res.gb2$XGB.mmr$eval[[i]]$test_auc, lwd = 0.5)
  abline(v = which.max(res.gb2$XGB.mmr$eval[[i]]$test_auc), col = "red", lwd = 0.5)
}

hist(res.gb$XGB.mmr$eval_best$AUC, breaks = 10, xlab = "Iteration that Maximizes AUC",
     main = "")


plot(1, type = "n", xlim = c(0, 300), ylim = c(0.9, 0.95), ylab = "AUC", xlab = "Iteration",
     main = "Starting with constant")
for(i in 1:100){
  lines(res.gb2$XGB.const$eval[[i]]$test_auc, lwd = 0.5)
  abline(v = which.max(res.gb2$XGB.const$eval[[i]]$test_auc), col = "red", lwd = 0.5)
}

hist(res.gb$XGB.const$eval_best$AUC, breaks = 10, xlab = "Iteration that Maximizes AUC",
     main = "")

plot(1, type = "n", xlim = c(0, 300), ylim = c(0.18, 0.22), ylab = "BS", xlab = "Iteration",
     main = "Starting with MMRpro")
for(i in 1:100){
  lines(res.gb2$XGB.mmr$eval[[i]]$test_rmse, lwd = 0.5)
  abline(v = which.min(res.gb2$XGB.mmr$eval[[i]]$test_rmse), col = "red", lwd = 0.5)
}

hist(res.gb$XGB.mmr$eval_best$BS, breaks = 10, xlab = "Iteration that Minimizes BS",
     main = "")


plot(1, type = "n", xlim = c(0, 300), ylim = c(0.18, 0.22), ylab = "BS", xlab = "Iteration",
     main = "Starting with constant")
for(i in 1:100){
  lines(res.gb2$XGB.const$eval[[i]]$test_rmse, lwd = 0.5)
  abline(v = which.min(res.gb2$XGB.const$eval[[i]]$test_rmse), col = "red", lwd = 0.5)
}

hist(res.gb$XGB.const$eval_best$BS, breaks = 10, xlab = "Iteration that Minimizes BS",
     main = "")


### EO

plot(1, type = "n", xlim = c(0, 300), ylim = c(0.85, 1.15), ylab = "EO", xlab = "Iteration",
     main = "XGBoost Starting with MMRpro")
for(i in 1:10){
  lines(res.gb$XGB.mmr$eval[[i]]$test_EO, lwd = 0.5)
  abline(v = which.min(res.gb$XGB.mmr$eval[[i]]$test_EO), col = "red", lwd = 0.5)
}

hist(res.gb$XGB.mmr$eval_best$EO, breaks = 10, xlab = "Iteration that Minimizes EO",
     main = "")


plot(1, type = "n", xlim = c(0, 300), ylim = c(0.3, 1.05), ylab = "EO", xlab = "Iteration",
     main = "Starting with constant")
for(i in 1:10){
  lines(res.gb$XGB.const$eval[[i]]$test_EO, lwd = 0.5)
  abline(v = which.min(res.gb$XGB.const$eval[[i]]$test_EO), col = "red", lwd = 0.5)
}

hist(res.gb$XGB.const$eval_best$EO, breaks = 10, xlab = "Iteration that Minimizes EO",
     main = "")


plot(res.gb$XGB.mmr$eval[[1]]$test_EO, xlab = "Iteration", ylab = "E/O", type = "l")
plot(res.gb$XGB.mmr$eval[[2]]$test_EO, xlab = "Iteration", ylab = "E/O", type = "l")


eo50 <- vector()
for(i in 1:100){
  eo50 <- c(eo50, res.gb$XGB.mmr$eval[[i]]$test_EO[50])
}

eo.best <- vector()
for(i in 1:100){
  eo.best <- c(eo.best, which.min(abs(res.gb$XGB.mmr$eval[[i]]$test_EO - 1)))
}
hist(eo.best, breaks = 10, xlab = "Iteration that Optimizes EO", main = "Starting with MMRpro")



setNames(data.frame(apply(best.mmr, 2, summary)), c("AUC", "BS"))

setNames(data.frame(apply(best.const, 2, summary)), c("AUC", "BS"))
