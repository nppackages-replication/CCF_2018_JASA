c.list  = c(-2/3,-1/3,0,1/3,2/3)
nbws    = 3
nmodels = 6
nevals  = 5
ncols   = 10

######################################
###### EC optimal BW #################
######################################
for (m in 1:nmodels) {
model1 = read.csv(paste("output/lp_bws_m",m,"_",vce,"_n",n,"_c",1,".csv",sep=""), row.names=1)
model2 = read.csv(paste("output/lp_bws_m",m,"_",vce,"_n",n,"_c",2,".csv",sep=""), row.names=1)
model3 = read.csv(paste("output/lp_bws_m",m,"_",vce,"_n",n,"_c",3,".csv",sep=""), row.names=1)
model4 = read.csv(paste("output/lp_bws_m",m,"_",vce,"_n",n,"_c",4,".csv",sep=""), row.names=1)
model5 = read.csv(paste("output/lp_bws_m",m,"_",vce,"_n",n,"_c",5,".csv",sep=""), row.names=1)

table.reg=matrix(NA,nevals*nbws,ncols)
table.reg[,1] = formatC(c(model1$bw,model2$bw,model3$bw,model4$bw,model5$bw),  format = "f", digits = 3)
table.reg[,2] = formatC(100*c(model1$ec.us,  model2$ec.us,  model3$ec.us,  model4$ec.us,  model5$ec.us),  format = "f", digits = 1)
table.reg[,3] = formatC(100*c(model1$ec.lf,  model2$ec.lf,  model3$ec.lf,  model4$ec.lf,  model5$ec.lf),  format = "f", digits = 1)
table.reg[,4] = formatC(100*c(model1$ec.bc,  model2$ec.bc,  model3$ec.bc,  model4$ec.bc,  model5$ec.bc),  format = "f", digits = 1)
table.reg[,5] = formatC(100*c(model1$ec.hh,  model2$ec.hh,  model3$ec.hh,  model4$ec.hh,  model5$ec.hh),  format = "f", digits = 1)
table.reg[,6] = formatC(100*c(model1$ec.rb, model2$ec.rb, model3$ec.rb, model4$ec.rb, model5$ec.rb), format = "f", digits = 1)
table.reg[,7] = formatC(c(model1$il.us,  model2$il.us,  model3$il.us,  model4$il.us,  model5$il.us),      format = "f", digits = 3)
table.reg[,8] = formatC(c(model1$il.lf,  model2$il.lf,  model3$il.lf,  model4$il.lf,  model5$il.lf),      format = "f", digits = 3)
table.reg[,9] = formatC(c(model1$il.hh,  model2$il.hh,  model3$il.hh,  model4$il.hh,  model5$il.hh),      format = "f", digits = 3)
table.reg[,10]= formatC(c(model1$il.rb, model2$il.rb, model3$il.rb, model4$il.rb, model5$il.rb),     format = "f", digits = 3)

for (j in c(1,4,7,10,13)) if (table.reg[j,1]>2) table.reg[j,]=rep("-",ncols)

colnames(table.reg)= c(" ","US","Locfit","BC","HH","RBC", "US","Locfit","HH","RBC")
rownames(table.reg)= rep(c("$h_{\\tt mse}$","$\\hat{h}_{\\tt rot}$","$\\hat{h}_{\\tt dpi}$"),nevals)

table1_tex = latex(table.reg, file = paste("tables/lp_bws_m",m,"_n",n,".txt",sep=""), landscape=FALSE, n.cgroup=c(1,5,4),
                   n.rgroup=c(nbws,nbws,nbws,nbws,nbws), 
                   rgroup = c('$x=-2/3$', '$x=-1/3$', '$x=0$', '$x=1/3$', '$x=2/3$'), 
                   cgroup = c('Bandwidth','Empirical Coverage', 'Interval Length'), 
                   outer.size='scriptsize', col.just=rep('c',10), center='none', title='', table.env=FALSE)

######################################
###### BW hat performance ############
######################################

model1 = read.csv(paste("output/lp_bws_all_m",m,"_",vce,"_n",n,"_c",1,".csv",sep=""), row.names=1)
model2 = read.csv(paste("output/lp_bws_all_m",m,"_",vce,"_n",n,"_c",2,".csv",sep=""), row.names=1)
model3 = read.csv(paste("output/lp_bws_all_m",m,"_",vce,"_n",n,"_c",3,".csv",sep=""), row.names=1)
model4 = read.csv(paste("output/lp_bws_all_m",m,"_",vce,"_n",n,"_c",4,".csv",sep=""), row.names=1)
model5 = read.csv(paste("output/lp_bws_all_m",m,"_",vce,"_n",n,"_c",5,".csv",sep=""), row.names=1)

table.bw = matrix(NA,nevals*2,8)
h.mse.pob=0
h.mse.pob[1] = model1$h.mse.pob[1]
h.mse.pob[2] = model2$h.mse.pob[1]
h.mse.pob[3] = model3$h.mse.pob[1]
h.mse.pob[4] = model4$h.mse.pob[1]
h.mse.pob[5] = model5$h.mse.pob[1]
h.mse.pob=formatC(h.mse.pob, format = "f", digits = 3)
for (j in 1:5) if (h.mse.pob[j]>1) h.mse.pob[j]="-"

table.bw[1,1]   = paste(formatC(h.mse.pob[1],  format = "f", digits = 3), sep = "")
table.bw[1,2:7] = summary(model1$h.mse.hat, digits = max(3))[1:6]
table.bw[1,8]   = formatC(sd(model1$h.mse.hat),format="f",digits=3)
table.bw[2,1]   = c("-")
table.bw[2,2:7] = summary(model1$h.cer.dpi, digits = max(3))[1:6]
table.bw[2,8]   = formatC(sd(model1$h.cer.dpi),format="f",digits=3)

table.bw[3,1]   = paste(formatC(h.mse.pob[2],  format = "f", digits = 3), sep = "")
table.bw[3,2:7] = summary(model2$h.mse.hat, digits = max(3))[1:6]
table.bw[3,8]   = formatC(sd(model2$h.mse.hat),format="f",digits=3)
table.bw[4,1]   = c("-")
table.bw[4,2:7] = summary(model2$h.cer.dpi, digits = max(3))[1:6]
table.bw[4,8]   = formatC(sd(model2$h.cer.dpi),format="f",digits=3)

table.bw[5,1]   = paste(formatC(h.mse.pob[3],  format = "f", digits = 3), sep = "")
table.bw[5,2:7] = summary(model3$h.mse.hat, digits = max(3))[1:6]
table.bw[5,8]   = formatC(sd(model3$h.mse.hat),format="f",digits=3)
table.bw[6,1]   = c("-")
table.bw[6,2:7] = summary(model3$h.cer.dpi, digits = max(3))[1:6]
table.bw[6,8]   = formatC(sd(model3$h.cer.dpi),format="f",digits=3)

table.bw[7,1]   = paste(formatC(h.mse.pob[4],  format = "f", digits = 3), sep = "")
table.bw[7,2:7] = summary(model4$h.mse.hat, digits = max(3))[1:6]
table.bw[7,8]   = formatC(sd(model4$h.mse.hat),format="f",digits=3)
table.bw[8,1]   = c("-")
table.bw[8,2:7] = summary(model4$h.cer.dpi, digits = max(3))[1:6]
table.bw[8,8]   = formatC(sd(model4$h.cer.dpi),format="f",digits=3)

table.bw[9,1]   = paste(formatC(h.mse.pob[5],  format = "f", digits = 3), sep = "")
table.bw[9,2:7] = summary(model5$h.mse.hat, digits = max(3))[1:6]
table.bw[9,8]   = formatC(sd(model5$h.mse.hat),format="f",digits=3)
table.bw[10,1]   = c("-")
table.bw[10,2:7] = summary(model5$h.cer.dpi, digits = max(3))[1:6]
table.bw[10,8]   = formatC(sd(model5$h.cer.dpi),format="f",digits=3)

colnames(table.bw)=c("Pop. Par.","Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","Std. Dev.")
rownames(table.bw)=rep(c("$\\hat{h}_{\\tt rot}$","$\\hat{h}_{\\tt dpi}$"),nevals)

table1_tex = latex(table.bw, file = paste("tables/lp_bws_stats_m",m,"_n",n,".txt",sep=""), landscape=FALSE,
                   n.rgroup=c(2,2,2,2,2), rgroup = c('$x=-2/3$', '$x=-1/3$', '$x=0$', '$x=1/3$', '$x=2/3$'),  
                   outer.size='scriptsize', col.just=rep('c',8), center='none', title='', table.env=FALSE)
}
##########################################################################################################################


######################################
###### US comparison 1 (lambda) ######
######################################
n=500
vce="hc3"
c.list = c(-2/3,-1/3,0,1/3,2/3)
for (model in 1:6) {
  out=matrix(NA,5,8)
  colnames(out)= rep(c("EC","IL"),4)
  for (eval in 1:5) {
    x0=c.list[eval]
    output1 = read.csv(paste("output/lp_comp_m", model,"_",vce,"_n",n,"_c",eval,".csv",sep=""), row.names=1)
    output2 = read.csv(paste("output/lp_bws_m",  model,"_",vce,"_n",n,"_c",eval,".csv",sep=""), row.names=1)
    ec.us1 = formatC(100*output1$ec.us1,format = "f", digits = 1)
    ec.us2 = formatC(100*output1$ec.us2,format = "f", digits = 1)
    ec.rb1 = formatC(100*output1$ec.rb,format = "f", digits = 1)
    ec.rb2 = formatC(100*output2$ec.rb[3],format = "f", digits = 1)
    il.us1 = formatC(output1$il.us1,format = "f", digits = 3)
    il.us2 = formatC(output1$il.us2,format = "f", digits = 3)
    il.rb1 = formatC(output1$il.rb,format = "f", digits = 3)
    il.rb2 = formatC(output2$il.rb[3],format = "f", digits = 3)
    #h.seq  = output$bw
    out[eval,]= c(ec.us1, il.us1, ec.us2, il.us2, ec.rb1, il.rb1, ec.rb2, il.rb2)
}
assign(paste("out",model, sep=""),out)
}
out=as.matrix(rbind(out1,out2,out3,out4,out5,out6))
rownames(out)=rep(c("$x=-2/3$","$x=-1/3$","$x=0$","$x=1/3$","$x=2/3$"),6)
table1_tex = latex(out, file = paste("tables/lp_comp2","_n",n,".txt",sep=""), landscape=FALSE, 
                   n.cgroup=c(2,2,2,2), cgroup = c("US ($\\lambda=0.5$)","US ($\\lambda=0.7$)","RBC ($\\hat{h}^{\\tt rot}_{\\tt rbc}$)","RBC ($\\hat{h}^{\\tt dpi}_{\\tt rbc}$)"), 
                   n.rgroup=rep(5,6), rgroup = c('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5', 'Model 6'), 
                   outer.size='scriptsize', col.just=rep('c',8), center='none', title='', table.env=FALSE)
########################################################################################################


######################################
###### US comparison 2 (size)   ######
######################################

for (model in 1:6) {
  out=matrix(NA,5,6)
  colnames(out)= c("$h$","EC","IL","$h$","EC","IL")
  for (eval in 1:5) {
    x0=c.list[eval]
    output1 = read.csv(paste("output/lp_grid_m",model,"_",vce,"_n",n,"_c",eval,".csv",sep=""), row.names=1)
    output2 = read.csv(paste("output/lp_bws_m", model,"_",vce,"_n",n,"_c",eval,".csv",sep=""), row.names=1)
    ec.us = output1$ec.us
    ec.rb = output1$ec.rb
    il.us = output1$il.us
    il.rb = output1$il.rb
    h.seq = output1$h
    h.mse = output2$bw[1]
    n.grid=length(h.seq)
    ec.1=rep(1,n.grid); ec.9=rep(.9,n.grid); ec.95=rep(.95,n.grid)
    
    ind.h95  = which(abs(ec.us-0.949)==min(abs(ec.us-0.949)))[1]
    ind.hmse = which(abs(ec.rb-0.949)==min(abs(ec.rb-0.949)))[1]
    hmse = h.seq[ind.hmse]
    h95  = h.seq[ind.h95]
    il.sc.us = il.us[ind.h95]
    il.sc.rb = il.rb[ind.hmse]
    ec.sc.us = ec.us[ind.h95]
    ec.sc.rb = ec.rb[ind.hmse]
  
    out[eval,]= c( formatC(h95,format = "f",digits = 3),formatC(100*ec.sc.us,  format = "f", digits = 1),formatC(il.sc.us,      format = "f", digits = 3), formatC(hmse,format = "f",digits = 3),formatC(100*ec.sc.rb,  format = "f", digits = 1),formatC(il.sc.rb,  format = "f", digits = 3))
  }
  assign(paste("out",model, sep=""),out)
}
out=as.matrix(rbind(out1,out2,out3,out4,out5,out6))
rownames(out)=rep(c("$x=-2/3$","$x=-1/3$","$x=0$","$x=1/3$","$x=2/3$"),6)
table1_tex = latex(out, file = paste("tables/lp_comp","_n",n,".txt",sep=""), landscape=FALSE, 
                   n.cgroup=c(3,3), cgroup = c('US','RBC'), 
                   n.rgroup=rep(5,6), rgroup = c('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5', 'Model 6'), 
                   outer.size='scriptsize', col.just=rep('c',6), center='none', title='', table.env=FALSE)

######################################
###### VCE comparison      ###########
######################################
m=5
n=500
#vce.list = c("nn","hc0","hc1","hc2","hc3")
for (eval in 1:5) {
  x0=c.list[eval]
  out=matrix(NA,5,3)
  colnames(out)= c("$h$","EC","IL")
    output = read.csv(paste("output/lp_vce_m",m,"_n",n,"_c",eval,".csv",sep=""), row.names=1)
    out[,1] = formatC(output$bw,format = "f",digits = 3)
    out[,2] = formatC(100*output$ec.rb, format="f", digits=1)
    out[,3] = formatC(output$il.rb, format="f", digits=3)
assign(paste("out",eval, sep=""),out)
}
out=as.matrix(rbind(out1,out2,out3,out4,out5))
rownames(out)=rep(c("$HC_0$","$HC_1$","$HC_2$","$HC_3$","$NN$"),5)
table1_tex = latex(out, file = paste("tables/lp_vce","_n",n,".txt",sep=""), landscape=FALSE, 
                   #n.cgroup=c(1,3), cgroup = c('','RBC'), 
                   n.rgroup=c(5,5,5,5,5), rgroup = c("$x=-2/3$","$x=-1/3$","$x=0$","$x=1/3$","$x=2/3$"), 
                   outer.size='scriptsize', col.just=rep('c',5), center='none', title='', table.env=FALSE)


######################################
###### Table Main Paper ##############
######################################
c.list = c(-2/3,-1/3,0,1/3,2/3)
nevals=5
ncols = 11
m = 5
model1 = read.csv(paste("output/lp_bws_m",m,"_",vce,"_n",n,"_c",1,".csv",sep=""), row.names=1)
model2 = read.csv(paste("output/lp_bws_m",m,"_",vce,"_n",n,"_c",2,".csv",sep=""), row.names=1)
model3 = read.csv(paste("output/lp_bws_m",m,"_",vce,"_n",n,"_c",3,".csv",sep=""), row.names=1)
model4 = read.csv(paste("output/lp_bws_m",m,"_",vce,"_n",n,"_c",4,".csv",sep=""), row.names=1)
model5 = read.csv(paste("output/lp_bws_m",m,"_",vce,"_n",n,"_c",5,".csv",sep=""), row.names=1)
ind = 3

table.reg=matrix(NA,nevals,ncols)
table.reg[,1]  = paste(c("-2/3","-1/3","0","1/3","2/3"))
table.reg[,2]  = formatC(c(model1$bw[ind], model2$bw[ind], model3$bw[ind], model4$bw[ind], model5$bw[ind]),  format = "f", digits = 3)
table.reg[,3]  = formatC(100*c(model1$ec.us[ind], model2$ec.us[ind], model3$ec.us[ind], model4$ec.us[ind], model5$ec.us[ind]), format="f", digits=1)
table.reg[,4]  = formatC(100*c(model1$ec.lf[ind], model2$ec.lf[ind], model3$ec.lf[ind], model4$ec.lf[ind], model5$ec.lf[ind]), format="f", digits=1)
table.reg[,5]  = formatC(100*c(model1$ec.bc[ind], model2$ec.bc[ind], model3$ec.bc[ind], model4$ec.bc[ind], model5$ec.bc[ind]), format="f", digits=1)
table.reg[,6]  = formatC(100*c(model1$ec.hh[ind], model2$ec.hh[ind], model3$ec.hh[ind], model4$ec.hh[ind], model5$ec.hh[ind]), format="f", digits=1)
table.reg[,7]  = formatC(100*c(model1$ec.rb[ind], model2$ec.rb[ind], model3$ec.rb[ind], model4$ec.rb[ind], model5$ec.rb[ind]), format="f", digits=1)
table.reg[,8]  = formatC(c(model1$il.us[ind], model2$il.us[ind], model3$il.us[ind], model4$il.us[ind], model5$il.us[ind]), format="f", digits=3)
table.reg[,9]  = formatC(c(model1$il.lf[ind], model2$il.lf[ind], model3$il.lf[ind], model4$il.lf[ind], model5$il.lf[ind]), format="f", digits=3)
table.reg[,10] = formatC(c(model1$il.hh[ind], model2$il.hh[ind], model3$il.hh[ind], model4$il.hh[ind], model5$il.hh[ind]), format="f", digits=3)
table.reg[,11] = formatC(c(model1$il.rb[ind], model2$il.rb[ind], model3$il.rb[ind], model4$il.rb[ind], model5$il.rb[ind]), format="f", digits=3)
colnames(table.reg)= c("","","US","Locfit","BC","HH","RBC", "US","Locfit","HH","RBC")
#rownames(table.reg)= rep(c("$h_{\\tt mse}$","$\\hat{h}_{\\tt rot}$","$\\hat{h}_{\\tt dpi}$"),nevals)

table1_tex = latex(table.reg, file = paste("tables/lp_main_m",m,"_n",n,".txt",sep=""), landscape=FALSE, n.cgroup=c(1,1,5,4),
                     cgroup = c('Evaluation Point','Average Bandwidth', 'Empirical Coverage', 'Interval Length'), 
                     outer.size='scriptsize', col.just=rep('c',11), center='none', title='', table.env=FALSE)
  
 