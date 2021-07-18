###############################################
##### BW GRID #################################
###############################################
c.list = c("-2","-1","0","1","2")
for (model in 1:4) {
  for (eval in 1:5) {
    x0=c.list[eval]
    output1 = read.csv(paste("output/kd_grid_m",model,"_c",eval,"_n",n,".csv",sep=""), row.names=1)
    output2 = read.csv(paste("output/kd_bws_m", model,"_c",eval,"_n",n,".csv",sep=""), row.names=1)
    
    h_seq = output1$bw; h_mse = output2$bw[1]
    ec_us = output1$ec.us;    ec_bc = output1$ec.bc;    ec_rb = output1$ec.rb
    il_us = output1$il.us;    il_rb = output1$il.rb
    d_ec = data.frame(h=h_seq, ec_us=ec_us, ec_bc=ec_bc, ec_rb=ec_rb)
    d_ec2 = melt(d_ec,id="h")
    d_il = data.frame(h=h_seq, il_us=il_us, il_rb=il_rb)
    d_il2 = melt(d_il,id="h")
    
    y.lab = ""
    x.lab = ""
    if (eval==1 | eval==4) y.lab = "Empirical Coverage"

    p <- ggplot(d_ec2, aes(x=h, y=value, group = variable, linetype = variable, color= variable))  + labs(color = "")
    p <- p + geom_line() +  coord_cartesian(xlim = c(0.1,0.7), ylim = c(0.5, 1)) + labs(x = x.lab, y = y.lab) 
    p <- p + geom_hline(aes(yintercept=0.9) ,  linetype=2) +  geom_hline(aes(yintercept=0.95) , linetype=2) +  geom_hline(aes(yintercept=1) , linetype=2)
    p <- p + geom_vline(aes(xintercept=h_mse), linetype=2, color=1) + theme(legend.justification = 'left', legend.position=c(0,0.2), legend.background = element_blank(), legend.key = element_blank())
    p <- p + scale_linetype_manual(name="", values = c("dotted","dashed", "solid"),labels=c("Undersmoothing", "Bias Corrected", "Robust Bias Corrected")) +
      scale_color_manual(name ="", values = c("blue", "red", "black"),labels=c("Undersmoothing", "Bias Corrected", "Robust Bias Corrected"))
    
    if (eval>1) {
      p <- p +  theme(legend.position="none")
    }
    if (eval==2 | eval==3 | eval==5 | eval==6) {
      p <- p + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
    }
    p
    ggsave(paste("plots/kd_grid_ec_m",model,"_c",eval,"_n",n,".pdf",sep=""), width=4, height=4)

    y.lab = ""
    x.lab = ""
    if (eval==1 | eval==4) y.lab = "Interval Length"
    p <- ggplot(d_il2, aes(x=h, y=value, group = variable, linetype = variable, color= variable))  + labs(color = "")
    p <- p + geom_line() +  coord_cartesian(xlim = c(0.1,0.7), ylim = c(0, 1)) + labs(x = x.lab, y = y.lab)
    p <- p + geom_vline(aes(xintercept=h_mse), linetype=2, color=1) + theme(legend.justification = 'left', legend.position=c(0,0.8), legend.background = element_blank(), legend.key = element_blank())
    p <- p + scale_linetype_manual(name="", values = c("dotted","solid"),labels=c("Undersmoothing", "Robust Bias Corrected")) +
      scale_color_manual(name ="", values = c("blue", "black"),labels=c("Undersmoothing", "Robust Bias Corrected"))
    
    if (eval>1) {
      p <- p +  theme(legend.position="none")
    }
    if (eval==2 | eval==3 | eval==5 | eval==6) {
      p <- p + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
    }
    
    p
    ggsave(paste("plots/kd_grid_il_m",model,"_c",eval,"_n",n,".pdf",sep=""), width=4, height=4)
  }
}

################################################
######### 3-D PLOTS ############################
################################################
h.seq = seq(0.02,0.8,0.02)
rho.seq = seq(0.2,2,0.2)
eval = c(-2,-1,0,1,2)
n=500
k=3
c = eval[k]
for (m in 1:4) {
    ec.data = read.csv(paste("output/kd_3d_ec_m",m,"_c",k,"_n",n,".csv",sep=""), row.names=1)
    il.data = read.csv(paste("output/kd_3d_il_m",m,"_c",k,"_n",n,".csv",sep=""), row.names=1)

    x <- h.seq
    y <- rho.seq
    z <- as.matrix(ec.data)
    
    colnames(z)=rownames(z)=NULL
    ec.RBC.plot <-wireframe(z,  screen = list(z = -60, x = -60), row.values=x, column.values=y, zlim=c(0.5,1),
                            zlab=list(label="",cex=0.4),
                            xlab=list(label="bandwidth h",cex=0.7),
                            ylab=list(label="ratio h/b",cex=0.7),
                            scales = list(cex=0.6, x = list(arrows=FALSE,at=c(seq(0.1,0.8,0.1))), y = list(arrows=FALSE,at=c(seq(0.4,2,0.4))), z=list(arrows=FALSE,at=c(seq(0.5,1,0.1)))),
                            drape = TRUE, colorkey = TRUE)
    pdf(paste("plots/kd_3D_ec_m",m,"_x",k,"_n",n,".pdf",sep=""),w=5,h=5)
    plot(ec.RBC.plot)
    dev.off()
    
    z <- as.matrix(il.data)
    colnames(z)=rownames(z)=NULL
   
    il.RBC.plot <-wireframe(z,  screen = list(z = -60, x = -60), row.values=x, column.values=y, zlim=c(0,2),
                            zlab=list(label="",cex=0.4),
                            xlab=list(label="bandwidth h",cex=0.7),
                            ylab=list(label="ratio h/b",cex=0.7),
                            scales = list(cex=0.6, x = list(arrows=FALSE,at=c(seq(0.1,0.8,0.1))), y = list(arrows=FALSE,at=c(seq(0.4,2,0.4))), z=list(arrows=FALSE,at=c(seq(0.4,2,0.4)))),
                            drape = TRUE, colorkey = TRUE)
    pdf(paste("plots/kd_3D_il_m",m,"_x",k,"_n",n,".pdf",sep=""),w=5,h=5)
    plot(il.RBC.plot)
    dev.off()
}

#trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1))


################################################
######### DGPS       ###########################
################################################

evalx = c(-2,-1,0,1,2)

mu.mn    = c(0,0,0)
Sigma.mn = c(1,1,1)
w.mn     = c(1,0,0)
fx   = function(x) w.mn[1]*dnorm(x,mu.mn[1],Sigma.mn[1]) + w.mn[2]*dnorm(x,mu.mn[2],Sigma.mn[2]) + w.mn[3]*dnorm(x,mu.mn[3],Sigma.mn[3]) 
x     <- evalx
y     <- fx(evalx)
data  <- data.frame(y,x)
plot2 <- ggplot(data, aes(x=evalx, y=y )) + theme_bw() + stat_function(fun = fx, colour = "red") + geom_point() + xlim(c(-3,3)) +  theme(axis.title.x=element_blank(),axis.title.y=element_blank())
plot2 
ggsave(paste("plots/kd_model",1,".pdf",sep=""), width=4, height=4)

mu.mn    = c(0, 1/2, 13/12)
Sigma.mn = c(1, 2/3,  5/9)
w.mn     = c(1/5, 1/5, 3/5)
fx   = function(x) w.mn[1]*dnorm(x,mu.mn[1],Sigma.mn[1]) + w.mn[2]*dnorm(x,mu.mn[2],Sigma.mn[2]) + w.mn[3]*dnorm(x,mu.mn[3],Sigma.mn[3]) 
x     <- evalx
y     <- fx(evalx)
data  <- data.frame(y,x)
plot2 <- ggplot(data, aes(x=evalx, y=y )) + theme_bw() + stat_function(fun = fx, colour = "red") + geom_point() + xlim(c(-3,3)) +  theme(axis.title.x=element_blank(),axis.title.y=element_blank())
plot2 
ggsave(paste("plots/kd_model",2,".pdf",sep=""), width=4, height=4)

mu.mn    = c(-1, 1, 0)
Sigma.mn = c(2/3, 2/3, 1)
w.mn     = c(1/2, 1/2, 0)
fx   = function(x) w.mn[1]*dnorm(x,mu.mn[1],Sigma.mn[1]) + w.mn[2]*dnorm(x,mu.mn[2],Sigma.mn[2]) + w.mn[3]*dnorm(x,mu.mn[3],Sigma.mn[3]) 
x     <- evalx
y     <- fx(evalx)
data  <- data.frame(y,x)
plot2 <- ggplot(data, aes(x=evalx, y=y )) + theme_bw() + stat_function(fun = fx, colour = "red") + geom_point() + xlim(c(-3,3)) +  theme(axis.title.x=element_blank(),axis.title.y=element_blank())
plot2 
ggsave(paste("plots/kd_model",3,".pdf",sep=""), width=4, height=4)

mu.mn = c(0, 3/2, 0)
Sigma.mn = c(1, 1/3, 1)
w.mn = c(3/4, 1/4, 0)
fx   = function(x) w.mn[1]*dnorm(x,mu.mn[1],Sigma.mn[1]) + w.mn[2]*dnorm(x,mu.mn[2],Sigma.mn[2]) + w.mn[3]*dnorm(x,mu.mn[3],Sigma.mn[3]) 
x     <- evalx
y     <- fx(evalx)
data  <- data.frame(y,x)
plot2 <- ggplot(data, aes(x=evalx, y=y )) + theme_bw() + stat_function(fun = fx, colour = "red") + geom_point() + xlim(c(-3,3)) +  theme(axis.title.x=element_blank(),axis.title.y=element_blank())
plot2 
ggsave(paste("plots/kd_model",4,".pdf",sep=""), width=4, height=4)