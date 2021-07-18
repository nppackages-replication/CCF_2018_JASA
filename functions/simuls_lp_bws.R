source("functions/simuls_functions.R")

gx.pob    = g.0.fun(c)
g.2.pob   = hessian(g.0.fun, x=c)
f0.pob    = fx(c)
s2.pob    = 1
C1.h      = C1.fun(p,v=0, K=k.fun)  
C2.h      = C2.fun(p,v=0, K=k.fun)
Vm0.pob   = C2.h*s2.pob/f0.pob 
B.h.pob   = g.2.pob/factorial(p+1);
D.h.pob   = D.fun(p, v=0, B=B.h.pob, K=k.fun);
N.h.pob   = N.fun(v=0, V=Vm0.pob);
h.mse.pob = h.fun(N.h.pob, D.h.pob, p, v=0, n)

nbws = 3
h.seq = matrix(NA,sim,nbws)
gx.us=gx.bc=se.us=se.rb=gx.lf=se.lf=gx.hh=se.hh=qz.hh=matrix(NA,sim,nbws)

# Loop
set.seed(2016)

showwhen = 1; showevery=100
for (i in 1:sim) {
  if (i==showwhen) {cat(paste("\nSimulations Completed:",i-1,"of",sim,"- Model:",model,"- x:",c,"-", Sys.time())); showwhen=showwhen+showevery}
  
  # Generate random data
  x = x.gen(n) 
  u = u.rnd(n)
  y = g.0.fun(x) + u

  lpbw      = lpbwselect(y=y, x=x, eval=c, p=p, rho=1, kernel=kernel, vce=vce, bwselect="all", interior=TRUE)
  h.mse.hat = lpbw$bws[1,2]
  h.ce.dpi  = lpbw$bws[1,6]
  h.seq[i,] = c(h.mse.pob, h.mse.hat, h.ce.dpi)

  for (j in 1:nbws) {
      lp.out = lprobust(y=y, x=x, eval=c, p=p, kernel=kernel, h=h.seq[i,j], vce=vce)
      gx.us[i,j] = lp.out$Estimate[5]
      gx.bc[i,j] = lp.out$Estimate[6]
      se.us[i,j] = lp.out$Estimate[7]
      se.rb[i,j] = lp.out$Estimate[8]
      
      fit = locfit(y~lp(x, deg = 1, h = h.seq[i,j]))
      pred = predict(fit,c(c),se.fit=TRUE)
      gx.lf[i,j] = pred$fit
      se.lf[i,j] = pred$se.fit
      
      hh = hh.fun(y=y, x=x, c=c, h=h.seq[i,j])
      gx.hh[i,j] = hh$gx.hh
      se.hh[i,j] = hh$se.hh
      qz.hh[i,j] = hh$qz.hh

  }
}

T.us = (gx.us - gx.pob) / se.us
T.bc = (gx.bc - gx.pob) / se.us
T.rb = (gx.bc - gx.pob) / se.rb
T.lf = (gx.lf - gx.pob) / se.lf
T.hh = (gx.hh - gx.pob) / se.hh

ec.us = colMeans(1*(abs(T.us)<=qz),    na.rm=T)
ec.bc = colMeans(1*(abs(T.bc)<=qz),    na.rm=T)
ec.rb = colMeans(1*(abs(T.rb)<=qz),    na.rm=T)
ec.lf = colMeans(1*(abs(T.lf)<=qz),    na.rm=T)
ec.hh = colMeans(1*(abs(T.hh)<=qz.hh), na.rm=T)

il.us = colMeans(2*qz*se.us,    na.rm=T)
il.rb = colMeans(2*qz*se.rb,    na.rm=T)
il.lf = colMeans(2*qz*se.lf,    na.rm=T)
il.hh = colMeans(2*qz.hh*se.hh, na.rm=T)

bws    = colMeans(h.seq, na.rm=T)
final.table.lp = matrix(NA,nbws,10)
colnames(final.table.lp)=c("bw", "ec.us", "ec.bc", "ec.rb", "ec.lf", "ec.hh", "il.us", "il.rb", "il.lf", "il.hh")
rownames(final.table.lp)=c("h.mse.pob", "h.mse.hat", "h.cer.dpi")

final.table.lp[,1]  = bws
final.table.lp[,2]  = ec.us
final.table.lp[,3]  = ec.bc
final.table.lp[,4]  = ec.rb
final.table.lp[,5]  = ec.lf
final.table.lp[,6]  = ec.hh
final.table.lp[,7]  = il.us
final.table.lp[,8]  = il.rb
final.table.lp[,9]  = il.lf
final.table.lp[,10] = il.hh
write.csv(final.table.lp, file=paste("simuls/jasa/output/lp_bws_m",model,"_",vce,"_n",n,"_c",k,".csv",sep=""))

parapam = cbind(h.seq, gx.us, gx.bc, gx.lf, gx.hh, se.us, se.rb, se.lf, se.hh)
colnames(parapam)=c("h.mse.pob", "h.mse.hat", "h.cer.dpi", 
                    paste("gx.us.", seq(1:nbws), sep=""),
                    paste("gx.bc.", seq(1:nbws), sep=""),
                    paste("gx.lf.", seq(1:nbws), sep=""),
                    paste("gx.hh.", seq(1:nbws), sep=""),
                    paste("se.us.", seq(1:nbws), sep=""),
                    paste("se.rb.", seq(1:nbws), sep=""),
                    paste("se.lf.", seq(1:nbws), sep=""),
                    paste("se.hh.", seq(1:nbws), sep=""))
write.csv(parapam, file=paste("output/lp_bws_all_m",model,"_",vce,"_n",n,"_c",k,".csv",sep=""))
