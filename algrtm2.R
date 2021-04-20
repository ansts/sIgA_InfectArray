# "Pipeline"
# filepr preprocesses .gpr files usin gprLocNorSVM to locally 
#   normalize them by SVR, the fitted bachkgrounf is recorded 
#   as B635 columns.
# cn holds the controls' sequences
# Use pepStat to get pSet and pnSet. No spots are flagged -
#   the relative difference of the duplicates does not excede 15% in v0.
#
# v0[N spots x cases] is data from pSet, 
# v[N spots x cases] - from pnSet
# vd[N spots x (cases + flag for even chip layout rows] is 
#   reordered so that the even chip rows are in the first half 
#   and odd - in the second. This separates the duplicates in 
#   two tables of the same layout, stacked on top of each other.
# v2[N1 (peptides-controls) x 2*cases] - duplicates are now separate 
#   variables and conrtrols are removed.
# vn[sames v2] - v2 after global normalization with cyclic loess 

require(pepStat)
require(limma)
require(Biobase)
require(Biostrings)
require(matrixStats)
require(mvoutlier)
require(eulerr)
require(umap)
require(reshape2)

filepr("IgA_gpr1/", sh=F)
mapFile="mapall.csv"
dirToParse="proc_IgA_gpr1"
cn=c("YPYDVPDYAG", "KEVPALTAVETGAT","GGGGG","EQKLISEEDL")
pSet=makePeptideSet(path=dirToParse, mapping.file = mapFile, use.flags = T, rm.control.list = cn, bgCorrect.method = "normexp")
pSet=summarizePeptides(pSet, summary="mean")
i=apply(pSet@assayData$exprs,1,function(x){all(!is.na(x))})
pepts=pSet@featureRange@elementMetadata@listData$peptide
v0=pSet@assayData$exprs
v0=v0[i,]
pepts=pepts[i]
v=normarr(v0, centered=F)
vn=normalizeCyclicLoess(v, iterations=1, method="affy")
#vc=ComBat(v,c(1,1,1,1,2,2,2,2))
vcnt=sweep(vn,2,apply(vn,2,median), FUN="-")
mapall=read.csv("mapall.csv",header = TRUE)
colnames(vn)=mapall[,3]

cvcnt=apply(vcnt,1,function(l){
  x=aggregate(l,by=list(as.double(dsn)), FUN=diff)
  y=aggregate(l,by=list(as.double(dsn)), FUN=mean)+2
  sqrt(sum((x/2)^2/(8*y^2)))
})
svcnt=vcnt[cvcnt<0.2,]
pepsvcnt=rownames(svcnt)

dsn<-mapall[,5]
dsn=relevel(dsn,ref="native")
calls=limfit_(dsn,svcnt+2, p=0.05)

pdf(file="Effects.pdf", width=10, height=10)
plot(as.data.frame(shinycalls), pch=16, col=cpl2(130)[rank(rowMeans(svcnt[rownames(shinycalls),]))], xlim=c(-1.5,1.5), ylim=c(-1.5,1.5))
dev.off()

calls_ind=lapply(1:6, function(i) limfit_(dsn,svcnt+2, p=0.05, coef_ind = i))
n=apply(combn(m:1,2),2,function(co){paste(dsn[co], sep="", collapse="-")})
names(calls_ind)=n[length(n):1]

calls_ind=lapply(calls_ind,function(l) l[order(l[,1], decreasing = T),])

Ags=read.csv(file="Epitopes_s.csv",quote="", stringsAsFactors = F)
colnames(Ags)=c("Seq","Org")
Ags=Ags[!(duplicated(Ags$Seq)),]
Ags=Ags[Ags$Seq!="",]
Ags=Ags[order(Ags$Seq),]
Orgs=unique(Ags$Org)
write.csv(Orgs, "Orgs.csv")
Orgs=read.csv("Orgs.csv")
Or=as.factor(Orgs[,3])
names(Or)=Orgs$x
x=Ags$Org; names(x)=Ags$Seq
pepsvcOr=as.character(Or[x[pepsvcnt]])
pepsvcOrcls=data.frame(Seq=pepsvcnt,Prev=pepsvcOr,'Fe2+'=1*(pepsvcnt %in% rownames(calls_ind$`Fe-native`)),Heme=1*(pepsvcnt %in% rownames(calls_ind$`Heme-native`)),'Low pH'=1*(pepsvcnt %in% rownames(calls_ind$`pH-native`)), stringsAsFactors = F)

plot(euler(pepsvcOrcls[pepsvcOrcls[,2]=="R",3:5]), quantities=T, main="Rare")
plot(euler(pepsvcOrcls[pepsvcOrcls[,2]=="C",3:5]), quantities=T, main="Common")
tt1=table(rowSums(pepsvcOrcls[,3:5]),pepsvcOrcls[,2])
x=chisq.test(tt1, correct = T, simulate.p.value = T)
print(x)
x$stdres

x=as.double(rownames(tt1))
y=tt1[,2]/tt1[,1]
summary(lm(y[2:4]~x[2:4]))

for (i in 2:4) {
  sink(file="Rare_Common_inCalls.txt",append = T)
  tb=table(pepsvcOrcls[,1],pepsvcOrcls[,i])
  colnames(tb)=c("Non-affected", "Affected")
  rownames(tb)=c("Common","Rare")
  print("")
  print("")
  print(colnames(pepsvcOrcls)[i])
  print(tb)
  print(chisq.test(tb))
  sink()
  pdf(file=paste(colnames(pepsvcOrcls)[i],"_incalls.pdf", collapse="",sep=""), width=5, height=5)
  plot(tb, col=1:2, main=colnames(pepsvcOrcls)[i])
  dev.off()
}


cindRC=lapply(calls_ind,function(l){
  r=rownames(l)
  si=Ags$Seq %in% r
  o=Ags$Org[si]
  s=Ags$Seq[si]
  names(o)=s
  data.frame("Effect"=l[,1],"Prevailence"=Or[o[r]])
})

r=rownames(calls)
si=Ags$Seq %in% r
o=Ags$Org[si]
s=Ags$Seq[si]
names(o)=s
cRC=cbind(calls[, c(1:6,8)], Or[o[r]])
colnames(cRC)[8]="RC"

cRCefu=sapply(c(1,2,4), function(i){
  tb=cindRC[[i]]
  boxplot(tb[,1]~tb[,2], notch=T, main=names(calls_ind)[i])
  u=wilcox.test(tb[,1]~tb[,2])
  ttb=table(tb[,1]>0,tb[,2])
  colnames(ttb)=c("Common","Rare")
  rownames(ttb)=c("Decreased","Increased")
  chi=chisq.test(ttb)
  sink(file="Rare_Common_byFactor_PosvsNeg.txt", append = T)
  print(names(cindRC)[i])
  print(ttb)
  pdf(file=paste(names(cindRC)[i],".pdf", collapse="",sep=""), width=5, height=5)
  plot(ttb, col=1:2, main=names(cindRC)[i])
  dev.off()
  print(chi)
  print(u)
  sink()
  return(list(ttb,u))
})

updown=t(sapply(cRCefu[1,],rowSums))
rownames(updown)=names(cindRC)[c(1,2,4)]
plot(as.table(updown))
sink(file="Rare_Common_byFactor_PosvsNeg.txt", append = T)
print("")
print("Distribution of increased and dicreased rectivity")
print(updown)
x=chisq.test(updown)
print(x)
print(x$stdres)
sink()

pdf(file="UpDown.pdf", width=5, height=5)
plot(as.table(updown), col=1:2, main="Increased/Decreased Reactivity by Agent")
dev.off()

chisq.test(updown[c(2,3),])

updowncls=lapply(calls_ind, function(l){
  list(Up=rownames(l)[l[,1]>0], Down=rownames(l)[l[,1]<0])
})

tud=sapply(updowncls,lengths)
tud1=tud[,c(1,2,4)]
chisq.test(t(tud1))

x=Ags$Org; names(x)=Ags$Seq
updownclsAg=lapply(updowncls,function(l1){
  lapply(l1, function(l2){
    y=cbind(l2,x[l2],as.character(Or[x[l2]]))
    rownames(y)=NULL
    y[order(y[,2]),]
  })
})

sink(file="UpDownCalls.txt")
print(updownclsAg)
sink()

topcalls=lapply(calls_ind,rownames)
names(topcalls)=names(calls_ind)

pdf(file="Venn.pdf", width=10, height=10)
plot(euler(topcalls[c(1,2,4)]), quantities=T)
dev.off()

tcnm=unique(unlist(topcalls[c(1,2,4)]))
tcnmud=t(sapply(tcnm, function(p){
  t0=rep(0,6)
  names(t0)=1:6
  x=sapply(updowncls[c(1,2,4)],function(l){
    sapply(l,function(l1){
      p %in% l1
    })
  })
  x=melt(x)
  t0[x[,3]]=1
  names(t0)=paste(x[,2],x[,1], sep="-")
  return(t0)
}))

table(apply(tcnmud,1,function(l) sum(2^((1:6)[l==1]))))
plot(euler(tcnmud[,c(1,3,5)]), main="Increased", quantities=T)
plot(euler(tcnmud[,c(2,4,6)]), main="Decreased", quantities=T)

pdf(file="UpDown_Venns.pdf", width = 16, height=10)
plot(euler(tcnmud, quantities=T))
dev.off()
                              # for increased
x=rbind(c(3,3),c(21,149-27))  #pH Fe2 overlap is almost significant
fisher.test(x)
x=rbind(c(3,1),c(21,149-25))  # heme pH overlap is significant
fisher.test(x, simulate.p.value = T)

                                # for decreased
x=rbind(c(15,34),c(26,149-75))  #pH Fe2 overlap is not significant
fisher.test(x)
x=rbind(c(20,21),c(54,149-95))  # heme pH overlap is not significant
fisher.test(x, simulate.p.value = T)
x=rbind(c(29,20),c(45,149-94))  # heme Fe overlap is trend - p=0.12
fisher.test(x, simulate.p.value = T)
