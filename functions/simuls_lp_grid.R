source("functions/simuls_functions.R")

nbws = length(h.seq)
gx.us=gx.bc=se.us=se.rb=matrix(NA,sim,nbws)
gx.pob = g.0.fun(c)
s2.pob  = 1

# Loop
set.seed(2016)
showwhen = 1; showevery=100
for (i in 1:sim) {
  if (i==showwhen) {cat(paste("\nSimulations Completed:",i-1,"of",sim,"- Model:",model,"- x:",k,"-", Sys.time())); showwhen=showwhen+showevery}
  
  # Generate random data
  x = x.gen(n) 
  u = u.rnd(n)
  y = g.0.fun(x) + u

  for (j in 1:nbws) {
      lp.out = lprobust(y=y, x=x, eval=c, p=p, h=h.seq[j], vce=vce)
      gx.us[i,j] = lp.out$Estimate[5]
      gx.bc[i,j] = lp.out$Estimate[6]
      se.us[i,j] = lp.out$Estimate[7]
      se.rb[i,j] = lp.out$Estimate[8]
  }
}

T.us = (gx.us - gx.pob) / se.us
T.bc = (gx.bc - gx.pob) / se.us
T.rb = (gx.bc - gx.pob) / se.rb

ec.us = colMeans(1*(abs(T.us)<=qz),    na.rm = T)
ec.bc = colMeans(1*(abs(T.bc)<=qz),    na.rm = T)
ec.rb = colMeans(1*(abs(T.rb)<=qz),     na.rm = T)

il.us = colMeans(2*qz*se.us,    na.rm = T)
il.rb = colMeans(2*qz*se.rb,    na.rm = T)

gx.us.mean = colMeans(gx.us, na.rm = T)
gx.bc.mean = colMeans(gx.bc, na.rm = T)

final.table.lp = matrix(NA,nbws,9)
colnames(final.table.lp)=c("h", "ec.us", "ec.bc", "ec.rb", "il.us", "il.rb", "gx.us","gx.bc", "gx.pob")
final.table.lp[,1]  = h.seq
final.table.lp[,2]  = ec.us
final.table.lp[,3]  = ec.bc
final.table.lp[,4]  = ec.rb
final.table.lp[,5]  = il.us
final.table.lp[,6]  = il.rb
final.table.lp[,7] = gx.us.mean
final.table.lp[,8] = gx.bc.mean
final.table.lp[,9] = gx.pob
write.csv(final.table.lp, file=paste("output/lp_grid_m",model,"_",vce,"_n",n,"_c",k,".csv",sep=""))







