###############################################
##### BW GRID #################################
###############################################
c.list = c("-2/3","-1/3","0","1/3","2/3")
for (model in 1:6) {
  for (eval in 1:5) {
    x0=c.list[eval]
    output1 = read.csv(paste("output/lp_grid_m",model,"_",vce,"_n",n,"_c",eval,".csv",sep=""), row.names=1)
    output2 = read.csv(paste("output/lp_bws_m", model,"_",vce,"_n",n,"_c",eval,".csv",sep=""), row.names=1)
    
    h_seq = output1$h; h_mse = output2$bw[1]
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
    p <- p + geom_line() +  coord_cartesian(xlim = c(0.1,0.4), ylim = c(0.5, 1)) + labs(x = x.lab, y = y.lab) 
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
    ggsave(paste("plots/lp_grid_ec_m",model,"_c",eval,"_n",n,".pdf",sep=""), width=4, height=4)

    y.lab = ""
    x.lab = ""
    if (eval==1 | eval==4) y.lab = "Interval Length"
    p <- ggplot(d_il2, aes(x=h, y=value, group = variable, linetype = variable, color= variable))  + labs(color = "")
    p <- p + geom_line() +  coord_cartesian(xlim = c(0.1,0.4), ylim = c(0, 1)) + labs(x = x.lab, y = y.lab)
    p <- p + geom_vline(aes(xintercept=h_mse), linetype=2, color=1) + theme(legend.justification = 'left', legend.position=c(0,0.2), legend.background = element_blank(), legend.key = element_blank())
    p <- p + scale_linetype_manual(name="", values = c("dotted","solid"),labels=c("Undersmoothing", "Robust Bias Corrected")) +
      scale_color_manual(name ="", values = c("blue", "black"),labels=c("Undersmoothing", "Robust Bias Corrected"))
    
    if (eval>1) {
      p <- p +  theme(legend.position="none")
    }
    if (eval==2 | eval==3 | eval==5 | eval==6) {
      p <- p + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
    }
    p
    ggsave(paste("plots/lp_grid_il_m",model,"_c",eval,"_n",n,".pdf",sep=""), width=4, height=4)
  }
}

###############################################
##### COLOR SEGMENTS ##########################
###############################################
h.seq     = seq(0.1,0.7,0.02)
ec_test   = c(0,0.6,0.7,0.8,0.9,seq(0.9,1,length.out=length(h.seq)-5))
ec.breaks = c(-0.01,0.5,0.6,0.7,0.8,0.9,0.91,0.92,0.93,0.94,0.95,1)
ncolors <- length(ec.breaks);  
cuts <- cut(ec_test, breaks=ec.breaks)
color.list <- colorRampPalette(c('red','green'))
cols=sort(unique((color.list(ncolors)[as.numeric(cuts)])))

lb=rep(-1,length(h.seq))
ub=rep(1,length(h.seq))
df = data.frame(h=h.seq, lb  = lb, ub=ub, cuts=cuts)
p_leg <- ggplot()  +
  geom_segment(aes(x=lb, y=h, xend = ub, yend=h,  colour=cuts), data = df) +
  scale_colour_manual(name = "Coverage",  values = rev(cols))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
g2 <- function(a.gplot){
  if (!gtable::is.gtable(a.gplot))
    a.gplot <- ggplotGrob(a.gplot)
  gtable::gtable_filter(a.gplot, 'guide-box', fixed=TRUE)
}

leg1<-g_legend(p_leg)
leg2<-g2(p_leg)

pdf(paste("plots/lp_seg_leg",".pdf",sep=""), onefile=FALSE)
grid.draw(leg1) 
dev.off()

for (m in 1:6) {
  for (eval in 1:5) {
    
    model = read.csv(paste("output/lp_grid_m",m,"_",vce,"_n",n,"_c",eval,".csv",sep=""), row.names=1)
    h <- model$h
    
    il.us = model$il.us;    il.rb = model$il.rb
    ec.us = model$ec.us;    ec.rb = model$ec.rb
    gx.us = model$gx.us;    gx.bc = model$gx.bc
    bias.us = gx.us - model$gx.pob;    bias.rb = gx.bc - model$gx.pob+1
    lower.us = bias.us - il.us/2;    lower.rb = bias.rb - il.rb/2
    upper.us = bias.us + il.us/2;    upper.rb = bias.rb + il.rb/2
    
    cut_us <- cut(ec.us, breaks=ec.breaks);    cvec_us <- color.list(ncolors)[as.numeric(cut_us)]
    cut_rb <- cut(ec.rb, breaks=ec.breaks);    cvec_rb <- color.list(ncolors)[as.numeric(cut_rb)]

    seg = data.frame(id  = c(rep(0,length(lower.us)),rep(1,length(lower.rb))), 
                     lb  = c(lower.us,lower.rb),
                     ub  = c(upper.us,upper.rb), h=rep(h,2), 
                     col = c(as.character(cvec_us), as.character(cvec_rb)),
                     cut = c(as.character(cut_us),as.character(cut_rb)))
    
    coll = sort(unique(c(as.character(cvec_us), as.character(cvec_rb))))
    p <- ggplot() +  coord_cartesian( ylim = c(0.1, 0.7)) + labs(x = "", y = "Bandwidth", color="Coverage") +
      geom_segment(aes(x=lb, y=h, xend = ub, yend=h, group=id, colour=cut), data = seg) +
      geom_vline(aes(xintercept=0), linetype=2, color=1) + geom_vline(aes(xintercept=1), linetype=2, color=1) +
      scale_x_continuous(breaks=c(0,1), labels=c("US", "RBC")) +
      scale_colour_manual(name = "Coverage", values = rev(coll))

    if (eval<6) {
      p <- p +  theme(legend.position="none")
    }
    if (eval==2 | eval==3 | eval==5 | eval==6) {
      p <- p + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
    }
    
    p
    ggsave(paste("plots/lp_seg_m",m,"_c",eval,"_n",n,".pdf",sep=""))
  }
}


###############################################
##### DGP ##########################
###############################################
evalx = c(-2/3,-1/3,0,1/3,2/3)

for (model in 1:6) {
  if (model==1) g.0.fun = function(x) {sin(4*x) + 2*exp(-4*16*x^2)}
  if (model==2) g.0.fun = function(x) {2*x+2*exp(-4*16*x^2)} 
  if (model==3) g.0.fun = function(x) {0.3*exp(-4*(2*x+1)^2) + 0.7*exp(-16*(2*x-1)^2) }
  if (model==4) g.0.fun = function(x) {x+5*dnorm(10*x)}
  if (model==5) g.0.fun = function(x) {sin(3*pi*x/2)/(1+18*x^2*(sign(x)+1))}
  if (model==6) g.0.fun = function(x) {sin(  pi*x/2)/(1+ 2*x^2*(sign(x)+1))}
  
  x     <- evalx
  y     <- g.0.fun(evalx)
  data  <- data.frame(y,x)
  plot2 <- ggplot(data, aes(x=evalx, y=y )) + theme_bw() + stat_function(fun = g.0.fun, colour = "red") + geom_point() +  theme(axis.title.x=element_blank(),axis.title.y=element_blank())
  plot2 
  ggsave(paste("plots/lp_model",model,".pdf",sep=""), width=4, height=4)
}



###############################################
##### MAIN PAPER  ##########################
###############################################

c.list = c(-2/3,-1/3,0,1/3,2/3)
vce = "hc3"
model = 5
eval = 3
n=500
x0=c.list[eval]

output1 = read.csv(paste("output/lp_grid_m",model,"_",vce,"_n",n,"_c",eval,".csv",sep=""), row.names=1)
output2 = read.csv(paste("output/lp_bws_m", model,"_",vce,"_n",n,"_c",eval,".csv",sep=""), row.names=1)

h_seq = output1$h; h_mse = output2$bw[1]
ec_us = output1$ec.us;ec_bc = output1$ec.bc
ec_rb = output1$ec.rb;il_us = output1$il.us
il_rb = output1$il.rb;

d_ec = data.frame(h=h_seq, ec_us=ec_us, ec_bc=ec_bc, ec_rb=ec_rb)
d_ec2 = melt(d_ec,id="h")
d_il = data.frame(h=h_seq, il_us=il_us, il_rb=il_rb)
d_il2 = melt(d_il,id="h")

p <- ggplot(d_ec2, aes(x=h, y=value, group = variable, linetype = variable, color= variable))  + labs(color = "")
p <- p + geom_line() +  coord_cartesian(xlim = c(0.1,0.4), ylim = c(0.5, 1)) + labs(x = "Bandwidth", y = "Empirical Coverage") 
p <- p + geom_hline(aes(yintercept=0.9) ,  linetype=2) +  geom_hline(aes(yintercept=0.95) , linetype=2) +  geom_hline(aes(yintercept=1) , linetype=2)
p <- p + geom_vline(aes(xintercept=h_mse), linetype=2, color=1) + theme(legend.justification = 'left', legend.position=c(0.1,0.2), legend.background = element_blank(), legend.key = element_blank())
p <- p + scale_linetype_manual(name="", values = c("dotted","dashed", "solid"),labels=c("Undersmoothing", "Bias Corrected", "Robust Bias Corrected")) +
  scale_color_manual(name ="", values = c("blue", "red", "black"),labels=c("Undersmoothing", "Bias Corrected", "Robust Bias Corrected"))
p
ggsave(paste("plots/lp_ec_main",".pdf",sep=""), width=4, height=4)

p <- ggplot(d_il2, aes(x=h, y=value, group = variable, linetype = variable, color= variable))  + labs(color = "")
p <- p + geom_line() +  coord_cartesian(xlim = c(0.1,0.4), ylim = c(0, 1)) + labs(x = "Bandwidth", y = "Interval Length") 
p <- p + geom_vline(aes(xintercept=h_mse), linetype=2, color=1) + theme(legend.justification = 'left', legend.position=c(0.1,0.2), legend.background = element_blank(), legend.key = element_blank())
p <- p + scale_linetype_manual(name="", values = c("dotted","solid"),labels=c("Undersmoothing", "Robust Bias Corrected")) +
  scale_color_manual(name ="", values = c("blue", "black"),labels=c("Undersmoothing", "Robust Bias Corrected"))
p
ggsave(paste("plots/lp_il_main",".pdf",sep=""), width=4, height=4)


h.seq = seq(0.1,0.7,0.02)
ec_test = c(0,0.6,0.7,0.8,0.9,seq(0.9,1,length.out=length(h.seq)-5))
ec.breaks=c(-0.01,0.5,0.6,0.7,0.8,0.9,0.91,0.92,0.93,0.94,0.95,1)
ncolors <- length(ec.breaks);  
cuts <- cut(ec_test, breaks=ec.breaks)
color.list <- colorRampPalette(c('red','green'))
cols=sort(unique((color.list(ncolors)[as.numeric(cuts)])))    
    
output = read.csv(paste("output/lp_grid_m",model,"_",vce,"_n",n,"_c",eval,".csv",sep=""), row.names=1)
h <- output$h

il.us = output$il.us;    il.rb = output$il.rb
ec.us = output$ec.us;    ec.rb = output$ec.rb
gx.us = output$gx.us;    gx.bc = output$gx.bc
bias.us = gx.us - output$gx.pob;    bias.rb = gx.bc - output$gx.pob+1
lower.us = bias.us - il.us/2;    lower.rb = bias.rb - il.rb/2
upper.us = bias.us + il.us/2;    upper.rb = bias.rb + il.rb/2

cut_us <- cut(ec.us, breaks=ec.breaks);    cvec_us <- color.list(ncolors)[as.numeric(cut_us)]
cut_rb <- cut(ec.rb, breaks=ec.breaks);    cvec_rb <- color.list(ncolors)[as.numeric(cut_rb)]

seg = data.frame(id  = c(rep(0,length(lower.us)),rep(1,length(lower.rb))), 
                 lb  = c(lower.us,lower.rb),
                 ub  = c(upper.us,upper.rb), h=rep(h,2), 
                 col = c(as.character(cvec_us), as.character(cvec_rb)),
                 cut = c(as.character(cut_us),as.character(cut_rb)))

coll = sort(unique(c(as.character(cvec_us), as.character(cvec_rb))))

p <- ggplot() +  coord_cartesian( ylim = c(0.1, 0.7)) + labs(x = "", y = "Bandwidth", color="Coverage") +
  geom_segment(aes(x=lb, y=h, xend = ub, yend=h, group=id, colour=cut), data = seg) +
  geom_vline(aes(xintercept=0), linetype=2, color=1) + geom_vline(aes(xintercept=1), linetype=2, color=1) +
  scale_x_continuous(breaks=c(0,1), labels=c("US", "RBC")) +
  scale_colour_manual(name = "Coverage", values = rev(coll)) +
  theme(legend.position="none")

p
ggsave(paste("plots/lp_seg_main",".pdf",sep=""), width=4, height=4)



