
table = read.table(pipe('grep "site" matrix-scan_2018-04-15.030611_of0RS4.ft'),sep="\t")

E2F4_table = table[table[,3]=="E2F4",]
E2F4_table = E2F4_table[abs(E2F4_table[,5])<1000 & abs(E2F4_table[,6])<1000,]

ETS1_table = table[table[,3]=="ETS1",]
ETS1_table = ETS1_table[abs(ETS1_table[,5])<1000 & abs(ETS1_table[,6])<1000,]

REST_table = table[table[,3]=="REST",]
REST_table = REST_table[abs(REST_table[,5])<1000 & abs(REST_table[,6])<1000,]

TCF3_table = table[table[,3]=="TCF3",]
TCF3_table = TCF3_table[abs(TCF3_table[,5])<1000 & abs(TCF3_table[,6])<1000,]

ZEB1_table = table[table[,3]=="ZEB1",]
ZEB1_table = ZEB1_table[abs(ZEB1_table[,5])<1000 & abs(ZEB1_table[,6])<1000,]


# Initialize vectors
E2F4 = rep(0,2000)
ETS1 = rep(0,2000)
REST = rep(0,2000)
TCF3 = rep(0,2000)
ZEB1 = rep(0,2000)

############# Fill up vector ##############
c_table = E2F4_table
c_vector = E2F4
#length of Binding Site
bsl = nchar(as.character(c_table[1,7]))

for( i in 1:dim(c_table)[1]) {
  line=c_table[i,]
  add = (line[,5]+1000):(line[,5]+1000+bsl)
  c_vector[add] = c_vector[add]+1
  }

E2F4 = c_vector
#

c_table = ETS1_table
c_vector = ETS1
#length of Binding Site
bsl = nchar(as.character(c_table[1,7]))

for( i in 1:dim(c_table)[1]) {
  line=c_table[i,]
  add = (line[,5]+1000):(line[,5]+1000+bsl)
  c_vector[add] = c_vector[add]+1
  }

ETS1 = c_vector
#

c_table = REST_table
c_vector = REST
#length of Binding Site
bsl = nchar(as.character(c_table[1,7]))

for( i in 1:dim(c_table)[1]) {
  line=c_table[i,]
  add = (line[,5]+1000):(line[,5]+1000+bsl)
  c_vector[add] = c_vector[add]+1
  }

REST = c_vector
#

c_table = TCF3_table
c_vector = TCF3
#length of Binding Site
bsl = nchar(as.character(c_table[1,7]))

for( i in 1:dim(c_table)[1]) {
  line=c_table[i,]
  add = (line[,5]+1000):(line[,5]+1000+bsl)
  c_vector[add] = c_vector[add]+1
  }

TCF3 = c_vector
#

c_table = ZEB1_table
c_vector = ZEB1
#length of Binding Site
bsl = nchar(as.character(c_table[1,7]))

for( i in 1:dim(c_table)[1]) {
  line=c_table[i,]
  add = (line[,5]+1000):(line[,5]+1000+bsl)
  c_vector[add] = c_vector[add]+1
  }

ZEB1 = c_vector
###########################################

NT = data.frame(E2F4=E2F4,ETS1=ETS1,REST=REST,TCF3=TCF3,ZEB1=ZEB1)
###########################################


table = read.table(pipe('grep "site" matrix-scan_2018-04-15.030853_ssr2zi.ft'),sep="\t")

E2F4_table = table[table[,3]=="E2F4",]
E2F4_table = E2F4_table[abs(E2F4_table[,5])<1000 & abs(E2F4_table[,6])<1000,]

ETS1_table = table[table[,3]=="ETS1",]
ETS1_table = ETS1_table[abs(ETS1_table[,5])<1000 & abs(ETS1_table[,6])<1000,]

REST_table = table[table[,3]=="REST",]
REST_table = REST_table[abs(REST_table[,5])<1000 & abs(REST_table[,6])<1000,]

TCF3_table = table[table[,3]=="TCF3",]
TCF3_table = TCF3_table[abs(TCF3_table[,5])<1000 & abs(TCF3_table[,6])<1000,]

ZEB1_table = table[table[,3]=="ZEB1",]
ZEB1_table = ZEB1_table[abs(ZEB1_table[,5])<1000 & abs(ZEB1_table[,6])<1000,]


# Initialize vectors
E2F4 = rep(0,2000)
ETS1 = rep(0,2000)
REST = rep(0,2000)
TCF3 = rep(0,2000)
ZEB1 = rep(0,2000)

############# Fill up vector ##############
c_table = E2F4_table
c_vector = E2F4
#length of Binding Site
bsl = nchar(as.character(c_table[1,7]))

for( i in 1:dim(c_table)[1]) {
  line=c_table[i,]
  add = (line[,5]+1000):(line[,5]+1000+bsl)
  c_vector[add] = c_vector[add]+1
  }

E2F4 = c_vector
#

c_table = ETS1_table
c_vector = ETS1
#length of Binding Site
bsl = nchar(as.character(c_table[1,7]))

for( i in 1:dim(c_table)[1]) {
  line=c_table[i,]
  add = (line[,5]+1000):(line[,5]+1000+bsl)
  c_vector[add] = c_vector[add]+1
  }

ETS1 = c_vector
#

c_table = REST_table
c_vector = REST
#length of Binding Site
bsl = nchar(as.character(c_table[1,7]))

for( i in 1:dim(c_table)[1]) {
  line=c_table[i,]
  add = (line[,5]+1000):(line[,5]+1000+bsl)
  c_vector[add] = c_vector[add]+1
  }

REST = c_vector
#

c_table = TCF3_table
c_vector = TCF3
#length of Binding Site
bsl = nchar(as.character(c_table[1,7]))

for( i in 1:dim(c_table)[1]) {
  line=c_table[i,]
  add = (line[,5]+1000):(line[,5]+1000+bsl)
  c_vector[add] = c_vector[add]+1
  }

TCF3 = c_vector
#

c_table = ZEB1_table
c_vector = ZEB1
#length of Binding Site
bsl = nchar(as.character(c_table[1,7]))

for( i in 1:dim(c_table)[1]) {
  line=c_table[i,]
  add = (line[,5]+1000):(line[,5]+1000+bsl)
  c_vector[add] = c_vector[add]+1
  }

ZEB1 = c_vector
###########################################
SH = data.frame(E2F4=E2F4,ETS1=ETS1,REST=REST,TCF3=TCF3,ZEB1=ZEB1)
###
# NT & SH
pdf("predicted_binding_sites_sameplot.pdf",width=12)
par(mfrow=c(2,3))
plot(NT[,1],type="l",col="red",xlab="ATAC-Seq peak center",ylab="# of Predicted Binding Sites",
     xaxt = "n",main=paste(colnames(NT)[1],"binding sites"),lwd=2)
lines(SH[,1],type="l",col="blue",lwd=2)
axis(1, at=c(1,500,1000,1500,2000), labels=c("-1000 bp","-500","0","500","1000 bp"))

plot(NT[,2],type="l",col="red",xlab="ATAC-Seq peak center",ylab="# of Predicted Binding Sites",
     xaxt = "n",main=paste(colnames(NT)[2],"binding sites"),lwd=2)
lines(SH[,2],type="l",col="blue",lwd=2)
axis(1, at=c(1,500,1000,1500,2000), labels=c("-1000 bp","-500","0","500","1000 bp"))

plot(NT[,3],type="l",col="red",xlab="ATAC-Seq peak center",ylab="# of Predicted Binding Sites",
     xaxt = "n",main=paste(colnames(NT)[3],"binding sites"),lwd=2)
lines(SH[,3],type="l",col="blue",lwd=2)
axis(1, at=c(1,500,1000,1500,2000), labels=c("-1000 bp","-500","0","500","1000 bp"))

plot(NT[,4],type="l",col="red",xlab="ATAC-Seq peak center",ylab="# of Predicted Binding Sites",
     xaxt = "n",main=paste(colnames(NT)[4],"binding sites"),lwd=2)
lines(SH[,4],type="l",col="blue",lwd=2)
axis(1, at=c(1,500,1000,1500,2000), labels=c("-1000 bp","-500","0","500","1000 bp"))

plot(NT[,5],type="l",col="red",xlab="ATAC-Seq peak center",ylab="# of Predicted Binding Sites",
     xaxt = "n",main=paste(colnames(NT)[5],"binding sites"),lwd=2)
lines(SH[,5],type="l",col="blue",lwd=2)
axis(1, at=c(1,500,1000,1500,2000), labels=c("-1000 bp","-500","0","500","1000 bp"))

plot.new()
legend("topleft",c("shNT","shH2AFV"),fill=c("red","blue"),bty = "n",cex=2.5)
dev.off()

