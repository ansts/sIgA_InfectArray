#
# Reads fin gpr file with p(ath) and writes  the fout with B(ackground) columns equal 
# to the corrected intensities of the empty chip and F columns equal to the corrected 
# intensities of the stained chip. The local normalization is based on SVR filtering
# of the background intensity determined on the bases of the lower 1/p0 quantile.
# Flags: wr(ite), sh(o)w 
# and R(ed/Green)sw(ap) - to make the green channel 
# under red column names for pepStat. The baseline is calculated 
# using SVR  
#
gprLNSVM_nosd<-function(fin, p=NULL, fout=NULL, wr=T, shw=F, Rsw=F) {
  require(parallel)
  require(reshape2)
  require(stringr)
  require(e1071)
  prt=proc.time()
  finp=paste(p,fin,sep = "")
  fcon=file(finp)
  f2l=readLines(con=fcon, n=2)
  nlns=as.double(str_extract(f2l[2], "[1-9]+"))
  fhead=readLines(con = fcon, n=nlns+2)
  f=read.delim(finp, skip=nlns+2, header=T, check.names = F, stringsAsFactors = FALSE)
  close(fcon)
  if (max(f$`Block`)>1) f=f[f$`Block`==3,]
  f=f[f$`ID`!="000000000000000",]
  dfim=data.frame(cbind(f$`Row`, f$`Column`, f$`F635 Median`))
  colnames(dfim)=c("R","C","V")
  N=nrow(dfim)
  nr=max(dfim[,1])
  nc=max(dfim[,2])
  if (shw==T){
    img0=acast(dfim, R~C, value.var = "V")
    image(img0, main=c(fin," Original"), zlim=c(0,65000), col=topo.colors(128))
  }

  r3=nr%/%10
  c3=nc%/%10
  btX=matrix(runif(nr*nc, min(dfim[,3]), max(dfim[,3])),nrow=nr)
  btX=melt(t(btX))
  x=btX[,1]
  btX[,1]=btX[,2]
  btX[,2]=x
  colnames(btX)=c("R","C","V")
  btX[1:nrow(dfim),]=dfim
  btD=dfim[dfim[,1] %in% (1:r3),]
  btD[,1]=-btD[,1]
  btU=dfim[dfim[,1] %in% ((nr-r3):nr),]
  btU[,1]=2*nr-btU[,1]+1

  btX=rbind(dfim,btD,btU)
  
  btL=btX[btX[,2] %in% (1:c3),]
  btL[,2]=-btL[,2]
  btR=btX[btX[,2] %in% ((nc-c3):nc),]
  btR[,2]=2*nc-btR[,2]+1
  
  btX=rbind(btX, btL, btR)
  nrX=max(btX[,1])
  ncX=max(btX[,2])
  NX=nrow(btX)
  
  co=FALSE 
  p0=3
  for (i in 1:NX){ 
      #print((i*100)%/%(N*2))
      rs=min(2,btX[i,1])
      cs=min(2,btX[i,2])
      p1 =min(p0,rs+cs)
      rb=max(btX[i,1]-rs,min(btX[,1]))
      rf=min(rb+2*rs+1, nr)
      cb=max(btX[i,2]-cs,min(btX[,2]))
      cf=min(cb+2*cs+1, nc)
      z0=btX[btX[,1] %in% rb:rf&btX[,2] %in% cb:cf,]
      lr=rank(z0[,3])
      z1=(z0[lr<length(lr)%/%p1,])                         # For each spot take the spots in a patch (rs.2+1)Rx(cs.2+1)C around it
      if (co==FALSE) {                                                            # and compare the spot value to the 100/p th percentile of the spots is the patch.
          btm=z1                                        # If less add this spot the btm set
          co=TRUE
      } 
      else {
      # if (dfim[i,3]<z1){                                                            
      #     btm=rbind(btm,dfim[i,])                 
      #}
          btm=rbind(btm,z1)
      }
  }
  print("bottom ready")
  frm=dfim[dfim[,1]==1|dfim[,1]==nr|dfim[,2]==1|dfim[,2]==nc,]
  btm=rbind(btm,frm)
  btm=unique(btm)

  
  lmd=svm(btm[,1:2], btm[,3], cost = 1000, gamma=3,epsilon = .001)     
  prd=predict(lmd,dfim[,1:2])
  lnew=dfim
  lnew[,3]=as.double(dfim[,3]-prd)
  print("bottom fit")

  lv=min(lnew[,3])-1                                               

  pred=prd+lv

  if (wr) {
    f$`B635`=pred
    f$`B635 Median`=pred
    f$`B635 Mean`=pred
    newp=paste("proc_",p,sep = "")                       # Write the new .grp files in 
    dir.create(newp)                                        # subfolder ...../.
    fout=paste(newp,fin,sep = "")
    fconw=file(fout, 'w')
    writeLines(fhead,con=fconw)
    write.table(f,file=fconw, row.names=F, sep = '\t')
    close(fconw)
    print(proc.time()-prt)
  }
}