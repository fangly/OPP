sink_file<-file("all.Rout",open="wt");
sink(sink_file, type="message");
args<-commandArgs(trailingOnly = T);
#args1 should be the file name of heat map
#args2 should be name for output svg file
#args3 should be rowsidecolors
HMold<-read.table(args[1],header=TRUE,sep="\t");
x<-(dim(HMold)[2]-1);
y<-(dim(HMold)[1]);
HM<-HMold[1:y,1:x];
library(gplots);
library(RColorBrewer);
rowcols<-read.table(args[3],header=FALSE);
pdf(args[2]);
heatmap.2(as.matrix(HM),Rowv=FALSE, trace='none',col=c("white",brewer.pal(7,"PuRd")), breaks=c(0,0.01,0.1,1,5,10,25,50,100), labRow="", key=FALSE,RowSideColors=as.vector(rowcols[,1]));
dev.off();
