#!/usr/bin/env Rscript
#!/usr/bin/env Rscript
# AUTHOR = Adam M. Session
args <- commandArgs(trailingOnly = TRUE)
pdfname<-paste("./LTRFamFiles/",args[3],".",args[1],".pdf",sep="")
pdf(file=pdfname,useDingbats=F)
require(ape);require(dendextend);mybreaks<-seq(0,1,by=.01);
alignname<-paste("./",args[2],"/",args[3],".",args[1],".Aligned.fa",sep="")
ltra<-as.matrix(read.table(args[4],row.names=1,header=T))
rowMin<-function (x) 
{
code = paste("x[,", 1:(NCOL(x)), "]", sep = "", collapse = ",")
code = paste("pmin(", code, ")")
return(eval(parse(text = code)))
}


SUBFAM<-read.dna(alignname,format="fasta");SUBFAM.dist<-dist.dna(SUBFAM,model="JC",pairwise.deletion=T,as.matrix=T);diag(SUBFAM.dist)<-NA;rownames(SUBFAM.dist)<-trimws(rownames(SUBFAM.dist));colnames(SUBFAM.dist)<-trimws(rownames(SUBFAM.dist));ltra.iter<-ltra[intersect(rownames(ltra),rownames(SUBFAM.dist)),];dend.dat<-as.dendrogram(hclust(as.dist(SUBFAM.dist)));labels_colors(dend.dat) <- ltra.iter[,2][order.dendrogram(dend.dat)];for(i in names(labels_colors(dend.dat))){labels_colors(dend.dat)[i]<-ltra.iter[i,2]};dend.dat %>% set("labels_cex",.5) %>%  plot(main=paste(args[1]," Blue=A Red=B Black=Unassigned\nSee key for other colors",sep=""));
a.list<-rownames(subset(ltra,ltra[,1]==args[5]));
b.list<-rownames(subset(ltra,ltra[,1]==args[6]));
ab.list<-c(a.list,b.list);
c.list<-rownames(subset(ltra,ltra[,1]==args[7]));


AtA<-SUBFAM.dist[intersect(rownames(SUBFAM.dist),a.list),intersect(rownames(SUBFAM.dist),a.list)];AtA.filt<-AtA;
if(length(na.omit(as.numeric(AtA.filt)))>=2){hist(AtA.filt,xlim=c(0,.5),freq=F,breaks=mybreaks);AtA.median<-median(na.omit(as.numeric(AtA.filt)));AtA.filt.iter<-AtA.filt;AtA.filt.iter[is.na(AtA.filt.iter)]<-1;AtA.BH<-rowMin(AtA.filt.iter);hist(AtA.BH,xlim=c(0,.5),freq=F,breaks=mybreaks);}

AtB<-SUBFAM.dist[intersect(rownames(SUBFAM.dist),a.list),intersect(rownames(SUBFAM.dist),b.list)];AtB.filt<-AtB;
if(length(na.omit(as.numeric(AtB.filt)))>=2){hist(AtB.filt,xlim=c(0,.5),freq=F,breaks=mybreaks);AtB.median<-median(na.omit(as.numeric(AtB.filt)));AtB.filt.iter<-AtB.filt;AtB.filt.iter[is.na(AtB.filt.iter)]<-1;AtB.BH<-rowMin(AtB.filt.iter);hist(AtB.BH,xlim=c(0,.5),freq=F,breaks=mybreaks);}

BtB<-SUBFAM.dist[intersect(rownames(SUBFAM.dist),b.list),intersect(rownames(SUBFAM.dist),b.list)];BtB.filt<-BtB;
if(length(na.omit(as.numeric(BtB.filt)))>=2){hist(BtB.filt,xlim=c(0,.5),freq=F,breaks=mybreaks);BtB.median<-median(na.omit(as.numeric(BtB.filt)));BtB.filt.iter<-BtB.filt;BtB.filt.iter[is.na(BtB.filt.iter)]<-1;BtB.BH<-rowMin(BtB.filt.iter);hist(BtB.BH,xlim=c(0,.5),freq=F,breaks=mybreaks);}

ABtC<-SUBFAM.dist[intersect(rownames(SUBFAM.dist),ab.list),intersect(rownames(SUBFAM.dist),c.list)];ABtC.filt<-ABtC;
if(length(na.omit(as.numeric(ABtC.filt)))>=2){hist(ABtC.filt,xlim=c(0,.5),freq=F,breaks=mybreaks);ABtC.median<-median(na.omit(as.numeric(ABtC.filt)));ABtC.filt.iter<-ABtC.filt;ABtC.filt.iter[is.na(ABtC.filt.iter)]<-1;ABtC.BH<-rowMin(ABtC.filt.iter);hist(ABtC.BH,xlim=c(0,.5),freq=F,breaks=mybreaks);}

dev.off()
DISTANCE.out<- paste(args[3],".",args[1],".dist.tbl",sep="");
MEDIAN.out<- paste(args[3],".",args[1],".median.tbl",sep="");
BH.out<- paste(args[3],".",args[1],".BH.tbl",sep="");

AtB[AtB>=.5]<-NA;
AtBfile<-paste("./LTRFamFiles/AtB.",DISTANCE.out,sep="");
BtB[BtB>=.5]<-NA;
BtBfile<-paste("./LTRFamFiles/BtB.",DISTANCE.out,sep="");
AtA[AtA>=.5]<-NA;
AtAfile<-paste("./LTRFamFiles/AtA.",DISTANCE.out,sep="");
ABtC[ABtC>=.5]<-NA;
ABtCfile<-paste("./LTRFamFiles/ABtC.",DISTANCE.out,sep="");

write.table(x=na.omit(as.numeric(AtA)),file=AtAfile,quote=F,row=F,col=F,sep="\n");
write.table(x=na.omit(as.numeric(BtB)),file=BtBfile,quote=F,row=F,col=F,sep="\n");
write.table(x=na.omit(as.numeric(AtB)),file=AtBfile,quote=F,row=F,col=F,sep="\n");
write.table(x=na.omit(as.numeric(ABtC)),file=ABtCfile,quote=F,row=F,col=F,sep="\n");

if(length(na.omit(as.numeric(AtB.filt)))>=2){
AtB.M.file<-paste("./LTRFamFiles/AtB.",MEDIAN.out,sep="")
write.table(x=median(na.omit(as.numeric(AtB))),file=AtB.M.file,quote=F,row=F,col=F,sep="\n")

AtB.BH.file<-paste("./LTRFamFiles/AtB.",BH.out,sep="")
write.table(x=na.omit(as.numeric(AtB.BH)),file=AtB.BH.file,quote=F,row=F,col=F,sep="\n")
}

if(length(na.omit(as.numeric(BtB.filt)))>=2){
BtB.M.file<-paste("./LTRFamFiles/BtB.",MEDIAN.out,sep="")
write.table(x=median(na.omit(as.numeric(BtB))),file=BtB.M.file,quote=F,row=F,col=F,sep="\n")

BtB.BH.file<-paste("./LTRFamFiles/BtB.",BH.out,sep="")
write.table(x=na.omit(as.numeric(BtB.BH)),file=BtB.BH.file,quote=F,row=F,col=F,sep="\n",append=T)
}

if(length(na.omit(as.numeric(AtA.filt)))>=2){
AtA.M.file<-paste("./LTRFamFiles/AtA.",MEDIAN.out,sep="")
write.table(x=median(na.omit(as.numeric(AtA))),file=AtA.M.file,quote=F,row=F,col=F,sep="\n",append=T)

AtA.BH.file<-paste("./LTRFamFiles/AtA.",BH.out,sep="")
write.table(x=na.omit(as.numeric(AtA.BH)),file=AtA.BH.file,quote=F,row=F,col=F,sep="\n",append=T)
}


if(length(na.omit(as.numeric(ABtC)))>=2){
ABtC.M.file<-paste("./LTRFamFiles/ABtC.",MEDIAN.out,sep="")
write.table(x=median(na.omit(as.numeric(ABtC))),file=ABtC.M.file,quote=F,row=F,col=F,sep="\n",append=T)

ABtC.BH.file<-paste("./LTRFamFiles/ABtC.",BH.out,sep="")
write.table(x=na.omit(as.numeric(ABtC.BH)),file=ABtC.BH.file,quote=F,row=F,col=F,sep="\n",append=T)
}
