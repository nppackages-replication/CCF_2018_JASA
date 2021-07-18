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

h.seq =gx.us=gx.bc=se.us=se.rb=matrix(NA,sim,5)
vce.list = c("hc0","hc1","hc2","hc3","nn")
# Loop
set.seed(2016)

showwhen = 1; showevery=100
for (i in 1:sim) {
  if (i==showwhen) {cat(paste("\nSimulations Completed:",i-1,"of",sim,"- Model:",model,"- x:",k,"-", Sys.time())); showwhen=showwhen+showevery}
  
  # Generate random data
  x = x.gen(n) 
  u = u.rnd(n)
  y = g.0.fun(x) + u

  for (j in 1:5) {
      lpbw      = lpbwselect(y=y, x=x, eval=c, p=p, kernel=kernel, vce=vce.list[j])
      h.mse.hat = lpbw$bws[2]
      h.seq[i,j]  = h.mse.hat
      lp.out = lprobust(y=y, x=x, eval=c, p=p, kernel=kernel, h=h.mse.hat, vce=vce.list[j])
      gx.us[i,j] = lp.out$Estimate[5]
      gx.bc[i,j] = lp.out$Estimate[6]
      se.us[i,j] = lp.out$Estimate[7]
      se.rb[i,j] = lp.out$Estimate[8]
   }
}

T.us = (gx.us - gx.pob) / se.us
T.bc = (gx.bc - gx.pob) / se.us
T.rb = (gx.bc - gx.pob) / se.rb

ec.us = colMeans(1*(abs(T.us)<=qz),    na.rm=T)
ec.bc = colMeans(1*(abs(T.bc)<=qz),    na.rm=T)
ec.rb = colMeans(1*(abs(T.rb)<=qz),    na.rm=T)

il.us = colMeans(2*qz*se.us,    na.rm=T)
il.rb = colMeans(2*qz*se.rb,    na.rm=T)

bws    = colMeans(h.seq, na.rm=T)
final.table.lp = matrix(NA,5,3)
colnames(final.table.lp)=c("bw", "ec.rb", "il.rb")
rownames(final.table.lp)=vce.list

final.table.lp[,1]  = bws
final.table.lp[,2]  = ec.rb
final.table.lp[,3]  = il.rb
write.csv(final.table.lp, file=paste("output/lp_vce_m",model,"_n",n,"_c",k,".csv",sep=""))

