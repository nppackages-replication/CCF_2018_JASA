source("functions/simuls_functions.R")

if (model==1) {
  mu.mn    = c(0,0,0)
  Sigma.mn = c(1,1,1)
  w.mn     = c(1,0,0)
}

if (model==2) {
  mu.mn    = c(0, 1/2, 13/12)
  Sigma.mn = c(1, 2/3,  5/9)
  w.mn     = c(1/5, 1/5, 3/5)
}

if (model==3) {
  mu.mn    = c(-1, 1, 0)
  Sigma.mn = c(2/3, 2/3, 1)
  w.mn     = c(1/2, 1/2, 0)
}
if (model==4) {
  mu.mn = c(0, 3/2, 0)
  Sigma.mn = c(1, 1/3, 1)
  w.mn = c(3/4, 1/4, 0)
}

fx   = function(x) w.mn[1]*dnorm(x,mu.mn[1],Sigma.mn[1]) + w.mn[2]*dnorm(x,mu.mn[2],Sigma.mn[2]) + w.mn[3]*dnorm(x,mu.mn[3],Sigma.mn[3]) 
fx.2 = function(x) w.mn[1]*(((x-mu.mn[1])/Sigma.mn[1]^2)^2-1/Sigma.mn[1]^2)*dnorm(x,mu.mn[1],Sigma.mn[1]) +  w.mn[2]*(((x-mu.mn[2])/Sigma.mn[2]^2)^2-1/Sigma.mn[2]^2)*dnorm(x,mu.mn[2],Sigma.mn[2]) + w.mn[3]*(((x-mu.mn[3])/Sigma.mn[3]^2)^2-1/Sigma.mn[3]^2)*dnorm(x,mu.mn[3],Sigma.mn[3])

ex <- norMix(mu = mu.mn, sigma = Sigma.mn, w = w.mn)
x.gen   = function(n)   { rnorMix(n,ex) }
f0 = fx(c)

var.fun = function(K,h) { (mean((K^2)) - mean(K)^2)/(n*h^2) }
nbws = length(h.seq)
fx.us=fx.bc=se.us=se.rb=matrix(NA,sim,nbws)

# Loop
set.seed(2016)

showwhen = 1; showevery=100

for (i in 1:sim) {
  if (i==showwhen) {cat(paste("\nSimulations Completed:",i-1,"of",sim,"- Model:",model,"- x:",k,"-", Sys.time())); showwhen=showwhen+showevery}
  x = x.gen(n) 
  for (j in 1:nbws) {
    kd.out = kdrobust(x=x, eval=c, kernel=kernel, h=h.seq[j])
    fx.us[i,j] = kd.out$Estimate[5]
    fx.bc[i,j] = kd.out$Estimate[6]
    se.us[i,j] = kd.out$Estimate[7]
    se.rb[i,j] = kd.out$Estimate[8]
  }
}

T.us  = (fx.us - f0) / se.us
T.bc  = (fx.bc - f0) / se.us
T.rb  = (fx.bc - f0) / se.rb
ec.us = colMeans(1*(abs(T.us)<=q), na.rm = T)
ec.bc = colMeans(1*(abs(T.bc)<=q), na.rm = T)
ec.rb = colMeans(1*(abs(T.rb)<=q), na.rm = T)
il.us = colMeans(2*q*se.us, na.rm = T)
il.rb = colMeans(2*q*se.rb, na.rm = T)

final.table = matrix(NA,nbws,6)
colnames(final.table)=c("bw", "ec.us", "ec.bc", "ec.rb", "il.us", "il.rb")

final.table[,1] = h.seq
final.table[,2] = ec.us
final.table[,3] = ec.bc
final.table[,4] = ec.rb
final.table[,5] = il.us
final.table[,6] = il.rb
write.csv(final.table, file=paste("output/kd_grid_m",model,"_c",k,"_n",n,".csv",sep=""))
