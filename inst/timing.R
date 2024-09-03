library(tictoc)

cnt <- 1
while(cnt <= 100){
  pp <- sample(1:20, 1)
  dd <- sample(1:12, 1)
  qq <- min(dd-1, max(1, round(dd/2) + sample(-1:3, 1)))
  NN <- khaos:::A_size(pp,dd,qq)
  if(NN < 100000){
    cat(cnt, "-", NN, "\n")
    tic()
    tmp <- khaos:::generate_A(pp,dd,qq)
    timer <- toc(quiet=TRUE)
    tmp <- data.frame(p=pp, d=dd, q=qq, N=NN, t=timer$toc - timer$tic)
    if(cnt == 1){
      res <- tmp
    }else{
      res <- rbind(res, tmp)
    }
    cnt <- cnt + 1
  }
}

plot(log(res$N), log(res$t))
plot(res$N, res$t)

ind1 <- which(log(res$N) >= 5)
ind2 <- which(log(res$N) <= 10.3)

res_sub <- res[intersect(ind1, ind2),]
plot(log(res_sub$N), log(res_sub$t))

fit <- lm(I(log(t))~I(log(N)), data=res_sub)
fit_n <- nls(t~b1+b2*N^b3, data=res, start=list(b1=0, b2=exp(-16.6), b3=1.93))
bb <- coefficients(fit_n)

plot(res$N, res$t, ylim=c(0, 60*10))
curve(bb[1]+bb[2]*x^bb[3], add=TRUE, col="orange", lwd=2)
points(res_sub$N, res_sub$t, col="dodgerblue", pch=1)
