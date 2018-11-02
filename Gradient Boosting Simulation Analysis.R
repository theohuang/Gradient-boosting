## Analyzing Gradient Boosting Simulation
## Last updated: October 11, 2018

res <- vector("list", 8)
for(i in 1:100){
  load(paste(getwd(), "/Gradient Boosting/Diff10/simgbpen_", i, ".RData", sep = ""))
  for(j in 1:8){
    res[[j]] <- rbind(res[[j]], colMeans(res.gb[[j]]$perf))
  }
}

# lapply(res, function(x) lapply(x, function(y) y$perf[1, ]))

# 
# res.summ <- matrix(0, 8, 9)
# for(i in 1:8){
#   for(j in 1:821){
#     if(j == 1){
#       res.i <- colMeans(res[[j]][[i]]$perf)
#     } else{
#       res.i <- rbind(res.i, colMeans(res[[j]][[i]]$perf))
#     }
#   }
#   res.summ[i, c(1, 4, 7)] <- colMeans(res.i)
#   res.summ[i, c(2, 5, 8)] <- apply(res.i, 2, function(x) sort(x)[3])
#   res.summ[i, c(3, 6, 9)] <- apply(res.i, 2, function(x) sort(x)[97])
# }
# 


res.summ <- matrix(0, 8, 9)
for(i in 1:8){
  res.summ[i, c(1, 4, 7)] <- colMeans(res[[i]])
  res.summ[i, c(2, 5, 8)] <- apply(res[[i]], 2, function(x) sort(x)[3])
  res.summ[i, c(3, 6, 9)] <- apply(res[[i]], 2, function(x) sort(x)[97])
}

row.names(res.summ) <- names(res.gb)
colnames(res.summ) <- c("EO mean", "EO 3", "EO 97",
                        "AUC mean", "AUC 3", "AUC 97",
                        "BS mean", "BS 3", "BS 97")

save(res, res.summ, file = "Diff10_Sim_Results.RDat")