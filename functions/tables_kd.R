ncols = 6
for (model in 1:4) {
  
model1 = read.csv(paste("output/kd_bws_m",model,"_c",1,"_n",n,".csv",sep=""), row.names=1)
model2 = read.csv(paste("output/kd_bws_m",model,"_c",2,"_n",n,".csv",sep=""), row.names=1)
model3 = read.csv(paste("output/kd_bws_m",model,"_c",3,"_n",n,".csv",sep=""), row.names=1)
model4 = read.csv(paste("output/kd_bws_m",model,"_c",4,"_n",n,".csv",sep=""), row.names=1)
model5 = read.csv(paste("output/kd_bws_m",model,"_c",5,"_n",n,".csv",sep=""), row.names=1)

table.den=matrix(NA,15,6)
table.den[,1]=formatC(    c(model1$bw,     model2$bw,     model3$bw,     model4$bw,     model5$bw),     format = "f", digits = 3)
table.den[,2]=formatC(100*c(model1$ec.us,  model2$ec.us,  model3$ec.us,  model4$ec.us,  model5$ec.us),  format = "f", digits = 1)
table.den[,3]=formatC(100*c(model1$ec.bc,  model2$ec.bc,  model3$ec.bc,  model4$ec.bc,  model5$ec.bc),  format = "f", digits = 1)
table.den[,4]=formatC(100*c(model1$ec.rb, model2$ec.rb, model3$ec.rb, model4$ec.rb, model5$ec.rb), format = "f", digits = 1)
table.den[,5]=formatC(    c(model1$il.us,  model2$il.us,  model3$il.us,  model4$il.us,  model5$il.us),  format = "f", digits = 3)
table.den[,6]=formatC(    c(model1$il.rb, model2$il.rb, model3$il.rb, model4$il.rb, model5$il.rb), format = "f", digits = 3)

for (j in c(1,4,7,10,13)) if (as.numeric(table.den[j,1]>10)) table.den[j,]=rep("-",ncols)
for (j in c(1,4,7,10,13)) if (table.den[j,1]==" Inf") table.den[j,]=rep("-",ncols)

colnames(table.den)= c(" ","US","BC","RBC","US","RBC")
rownames(table.den)=rep(c("$h_{\\tt mse}$","$\\hat{h}_{\\tt rot}$","$\\hat{h}_{\\tt dpi}$"),5)

table_tex = latex(table.den , file = paste("tables/kd_bws_m",model,"_n",n,".txt",sep=""), 
                  cgroup = c("Bandwidth","Empirical Coverage", "Interval Length"), n.cgroup = c(1,3,2),   
                  rgroup = c("$x=-2$", "$x=-1$", "$x=0$", "$x=1$", "$x=2$"),       n.rgroup = c(3,3,3,3,3), 
                  landscape=FALSE, center='none', col.just=rep('c',6), title='', table.env=FALSE)

### BW Table
model1 = read.csv(paste("output/kd_bws_all_m",model,"_c",1,"_n",n,".csv",sep=""), row.names=1)
model2 = read.csv(paste("output/kd_bws_all_m",model,"_c",2,"_n",n,".csv",sep=""), row.names=1)
model3 = read.csv(paste("output/kd_bws_all_m",model,"_c",3,"_n",n,".csv",sep=""), row.names=1)
model4 = read.csv(paste("output/kd_bws_all_m",model,"_c",4,"_n",n,".csv",sep=""), row.names=1)
model5 = read.csv(paste("output/kd_bws_all_m",model,"_c",5,"_n",n,".csv",sep=""), row.names=1)
table.bw = matrix(NA,10,8)

h.mse.pob=0
h.mse.pob[1] = model1$h.mse.pob[1]
h.mse.pob[2] = model2$h.mse.pob[1]
h.mse.pob[3] = model3$h.mse.pob[1]
h.mse.pob[4] = model4$h.mse.pob[1]
h.mse.pob[5] = model5$h.mse.pob[1]
for (j in 1:5) if (h.mse.pob[j]>10) h.mse.pob[j]=Inf

h.mse.pob=formatC(h.mse.pob, format = "f", digits = 3)
for (j in 1:5) if (h.mse.pob[j]==" Inf") h.mse.pob[j]="-"

table.bw[1,1]   = paste(formatC(h.mse.pob[1],  format = "f", digits = 3), sep = "")
table.bw[1,2:7] = summary(model1$h.mse.rot, digits = max(3))[1:6]
table.bw[1,8]   = format(sd(model1$h.mse.rot),format="f",digits=2,nsmall=2)
table.bw[2,1]   = c("-")
table.bw[2,2:7] = summary(model1$h.cer.dpi, digits = max(3))[1:6]
table.bw[2,8]   = format(sd(model1$h.cer.dpi),format="f",digits=2,nsmall=2)

table.bw[3,1]   = paste(formatC(h.mse.pob[2],  format = "f", digits = 3), sep = "")
table.bw[3,2:7] = summary(model2$h.mse.rot, digits = max(3))[1:6]
table.bw[3,8]   = format(sd(model2$h.mse.rot),format="f",digits=2,nsmall=2)
table.bw[4,1]   = c("-")
table.bw[4,2:7] = summary(model2$h.cer.dpi, digits = max(3))[1:6]
table.bw[4,8]   = format(sd(model2$h.cer.dpi),format="f",digits=2,nsmall=2)

table.bw[5,1]   = paste(formatC(h.mse.pob[3],  format = "f", digits = 3), sep = "")
table.bw[5,2:7] = summary(model3$h.mse.rot, digits = max(3))[1:6]
table.bw[5,8]   = format(sd(model3$h.mse.rot),format="f",digits=2,nsmall=2)
table.bw[6,1]   = c("-")
table.bw[6,2:7] = summary(model3$h.cer.dpi, digits = max(3))[1:6]
table.bw[6,8]   = format(sd(model3$h.cer.dpi),format="f",digits=2,nsmall=2)

table.bw[7,1]   = paste(formatC(h.mse.pob[4],  format = "f", digits = 3), sep = "")
table.bw[7,2:7] = summary(model4$h.mse.rot, digits = max(3))[1:6]
table.bw[7,8]   = format(sd(model4$h.mse.rot),format="f",digits=2,nsmall=2)
table.bw[8,1]   = c("-")
table.bw[8,2:7] = summary(model4$h.cer.dpi, digits = max(3))[1:6]
table.bw[8,8]   = format(sd(model4$h.cer.dpi),format="f",digits=2,nsmall=2)

table.bw[9,1]   = paste(formatC(h.mse.pob[5],  format = "f", digits = 3), sep = "")
table.bw[9,2:7] = summary(model5$h.mse.rot, digits = max(3))[1:6]
table.bw[9,8]   = format(sd(model5$h.mse.rot),format="f",digits=2,nsmall=2)
table.bw[10,1]   = c("-")
table.bw[10,2:7] = summary(model5$h.cer.dpi, digits = max(3))[1:6]
table.bw[10,8]   = format(sd(model5$h.cer.dpi),format="f",digits=2,nsmall=2)

colnames(table.bw)=c("Pop. Par.","Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","Std. Dev.")
rownames(table.bw)=rep(c("$\\hat{h}_{\\tt rot}$","$\\hat{h}_{\\tt dpi}$"),5)

table1_tex = latex(table.bw, file = paste("tables/kd_bws_stats_m",model,"_n",n,".txt",sep=""), landscape=FALSE,
                   rgroup = c("$x=-2$", "$x=-1$", "$x=0$", "$x=1$", "$x=2$"),       n.rgroup = c(2,2,2,2,2), 
                   outer.size='scriptsize', col.just=rep('c',8), center='none', title='', table.env=FALSE)
}