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

#nbws = length(h.seq)
#tau.us=tau.bc=se.us=se.rb=matrix(NA,sim,nbws)
showwhen = 1; showevery=100

#var.us1=var.us2=var.bc=var.rb = matrix(0, nrow=sim, ncol=length(h.seq))
#T.us1=T.us2=T.bc=T.rb=matrix(0, nrow=sim, ncol=length(h.seq))
#ec.us1=ec.us2 = il.us1=il.us2=matrix(0, ncol=1, nrow=length(h.seq))
ec.rb = il.rb = matrix(0, ncol=length(rho.seq), nrow=length(h.seq))

showwhen = 1; showevery=10

for (i in 1:sim) {
    if (i==showwhen) {cat(paste("\nSimulations Completed:",i-1,"of",sim,"- Model:",model,"-", Sys.time())); showwhen=showwhen+showevery}
    x = x.gen(n) 
      for (j in 1:length(h.seq)) {
            for (r in 1:length(rho.seq)) {
              kd.out = kdrobust(x=x, eval=c, kernel=kernel, h=h.seq[j], rho=rho.seq[r])
              fx.bc = kd.out$Estimate[6]
              se.rb = kd.out$Estimate[8]
              T.rb = (fx.bc  - f0) / se.rb
		          ec.rb[j,r] = ec.rb[j,r] + 1*(abs(T.rb)<=q)
              il.rb[j,r] = il.rb[j,r] + 2*q*se.rb
		        }
	      }
}

aec.rb = ec.rb/sim
ail.rb = il.rb/sim

#colnames(aec.rb) = colnames(ail.rb) = paste(rho.seq)
#rownames(aec.rb) = rownames(ail.rb) = paste(h.seq)
write.csv(aec.rb, file=paste("output/kd_3d_ec_m",model,"_c",k,"_n",n,".csv",sep=""))
write.csv(ail.rb, file=paste("output/kd_3d_il_m",model,"_c",k,"_n",n,".csv",sep=""))
#x = h.seq
#y = rho.seq
#z = aec.rb
#coverage.RBC.plot <- wireframe(z, drape=TRUE, colorkey=FALSE, row.values=x, column.values=y, zlim=c(0,1), zlab=list(label="",cex=0.4), screen = list(z = -60, x = -60), xlab=list(label="Bandwidth h",cex=0.7), ylab=list(label="ratio h/b",cex=0.7), at=seq(0,2,by=0.1), main=list(label="",cex=0.8), sub=list(label="", cex=0.7), scales=list(cex=0.6, x=list(arrows=FALSE, at=c(0.5,1,1.5)), y=list(arrows=FALSE, at=c(0.5,1,1.5,2)), z=list(arrows=FALSE, at=c(0,0.5,1))))
#pdf(paste("plots/kd_3d_ec_m",model,"_eval",k,"_n",n,".pdf",sep=""))
#plot(coverage.RBC.plot)
#dev.off()

#z <- ail.rb
#il.RBC.plot <- wireframe(z, drape=TRUE, colorkey=FALSE, row.values=x, column.values=y, zlab=list(label="",cex=0.4), screen = list(z = -60, x = -60), xlab=list(label="Bandwidth h",cex=0.7), ylab=list(label="ratio h/b",cex=0.7), at=seq(0,2,by=0.1), main=list(label="",cex=0.8), sub=list(label="", cex=0.7), scales=list(cex=0.6, x=list(arrows=FALSE, at=c(0.5,1,1.5)), y=list(arrows=FALSE, at=c(0.5,1,1.5,2)), z=list(arrows=FALSE)))
#pdf(paste("plots/kd_3d_il_m",model,"_eval",k,"_n",n,".pdf",sep=""),w=5,h=5)
#plot(il.RBC.plot)
#dev.off()



