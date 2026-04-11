for(i in 0:4){
  local({
    ii <- i
    tmp <- function(X, y){
      sparse_khaos(X, y, enrichment = ii)
    }
    assign(paste0("fit", ii), tmp, envir = .GlobalEnv)
  })
}

pred <- function(obj, Xt){
  predict(obj, Xt, nugget=TRUE)
}

fit_list <- list(fit0, fit1, fit2, fit3, fit4)
prd_list <- list(pred, pred, pred, pred, pred)

library(duqling)
res <- run_sim_study(fit_list, prd_list,
                     fnames = get_paper_funcs()[1:30],
                     n_train=500, NSR=0,
                     replications=10,
                     mc_cores=5)
