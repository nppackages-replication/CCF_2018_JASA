source("functions/simuls_functions.R")

gx.pob = g.0.fun(c)
g.2.pob = hessian(g.0.fun, x=c)
f0.pob  = fx(c)
s2.pob  = 1
C1.h = C1.fun(p,v=0, K=k.fun)  
C2.h = C2.fun(p,v=0, K=k.fun)
Vm0.pob = C2.h*s2.pob/f0.pob 
B.h.pob = g.2.pob/factorial(p+1);
D.h.pob = D.fun(p, v=0, B=B.h.pob, K=k.fun);
N.h.pob = N.fun(v=0, V=Vm0.pob);
h.mse.pob = h.fun(N.h.pob, D.h.pob, p, v=0, n)

h.seq = 0
gx.us1=gx.us2=gx.bc=se.us1=se.us2=se.rb=0

# Loop
set.seed(2016)

showwhen = 1; showevery=100
for (i in 1:sim) {
  if (i==showwhen) {cat(paste("\nSimulations Completed:",i-1,"of",sim,"- Model:",model,"- x:",k,"-", Sys.time())); showwhen=showwhen+showevery}
  
  # Generate random data
  x = x.gen(n) 
  u = u.rnd(n)
  y = g.0.fun(x) + u

  lpbw      = lpbwselect(y=y, x=x, eval=c, p=p, kernel=kernel, vce=vce)
  h.mse.hat = lpbw$bws[2]
  b.mse.hat = lpbw$bws[3]
  h.seq[i]  = h.mse.hat

      lp.out1 = lprobust(y=y, x=x, eval=c, p=p, kernel=kernel, h=0.5*h.mse.hat, vce=vce)
      gx.us1[i] = lp.out1$Estimate[5]
      se.us1[i] = lp.out1$Estimate[7]
      lp.out2 = lprobust(y=y, x=x, eval=c, p=p, kernel=kernel, h=0.7*h.mse.hat, vce=vce)
      gx.us2[i] = lp.out2$Estimate[5]
      se.us2[i] = lp.out2$Estimate[7]
      lp.out3 = lprobust(y=y, x=x, eval=c, p=p, kernel=kernel, h=h.mse.hat, vce=vce)
      gx.bc[i] = lp.out3$Estimate[6]
      se.rb[i] = lp.out3$Estimate[8]
}

T.us1 = (gx.us1 - gx.pob) / se.us1
T.us2 = (gx.us2 - gx.pob) / se.us2
T.rb  = (gx.bc  - gx.pob) / se.rb
ec.us1 = mean(1*(abs(T.us1)<=qz),   na.rm=T)
ec.us2 = mean(1*(abs(T.us2)<=qz),   na.rm=T)
ec.rb  = mean(1*(abs(T.rb)<=qz),    na.rm=T)
il.us1 = mean(2*qz*se.us1,  na.rm=T)
il.us2 = mean(2*qz*se.us2,  na.rm=T)
il.rb  = mean(2*qz*se.rb,   na.rm=T)

bws = mean(h.seq, na.rm=T)
final.table.lp = matrix(NA,1,7)
colnames(final.table.lp)=c("bw", "ec.us1", "ec.us2", "ec.rb", "il.us1", "il.us2", "il.rb")

final.table.lp[,1]  = bws
final.table.lp[,2]  = ec.us1
final.table.lp[,3]  = ec.us2
final.table.lp[,4]  = ec.rb
final.table.lp[,5]  = il.us1
final.table.lp[,6]  = il.us2
final.table.lp[,7]  = il.rb
write.csv(final.table.lp, file=paste("output/lp_comp_m",model,"_",vce,"_n",n,"_c",k,".csv",sep=""))
