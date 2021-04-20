# A wrapper for the limma linear model based extraction of calls

limfit_<-function(fdsn,vals, p=1, coef_ind=NULL, coef_cntr=1,n_ind=4000, n_cntr=4000){
  require(limma)
  require(sva)
  require(parallel)
  m=length(unique(fdsn))
  ncl=ncol(vals)
  if (ncl!=length(fdsn)) {return("Error: ncol of design and matrix differ!")}
  clin=character(0)
  design = model.matrix(~0+fdsn)
  cn=levels(fdsn)
  colnames(design)=cn
  vv=voom(2^vals, design,span=0.3,plot=T)
  cntr=apply(combn(m:1,2),2,function(co){paste(fdsn[co], sep="", collapse="-")})
  cntr=cntr[length(cntr):1]
  cmx=makeContrasts(contrasts=cntr, levels=cn)
  fit=lmFit(vv,design)
  fit2=contrasts.fit(fit,cmx)
  fit2=eBayes(fit2, proportion = 0.1)
  #fit1=eBayes(fit, proportion = 0.1)
  tptb2=topTable(fit2,  number=n_cntr, p.value=p, adjust.method = "BH", coef = coef_ind) #,adjust.method="holm"    p.value = p, 
  #tptb1=topTable(fit1, coef=coef_ind, number=n_ind, p.value=p, adjust.method = "BH") #,adjust.method="holm"
  return(tptb2)
}  

