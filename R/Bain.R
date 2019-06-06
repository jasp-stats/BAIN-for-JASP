
#function
Bain<-function(estimate, Sigma, grouppara = 0, jointpara = 0, n, ERr = NULL, IRr = NULL, ..., seed = 100, print = TRUE)
{

  numSP=length(estimate)
  totalRr<-list(ERr,IRr,...)
  numH<-length(totalRr)/2

  totalRr<-lapply(totalRr, function(x) if(is.null(x)){x = matrix(0,0,0)}else{x = x})

  if(length(ERr)==0&&length(IRr)==0){stop("at least one restriction matrix has to be specified")}
  if(any(lapply(totalRr,ncol)!=(numSP+1)&lapply(totalRr,ncol)!=0)){stop("number of columns in the restriction matrix is not equal to number of parameters")}

  #number of equality and/or inequality constraints in each hypothesis
  numR<-unlist(lapply(totalRr,nrow))
  j=0
  for(i in 1:(2*numH)){
    if(numR[i]==0){totalRr=totalRr[-(i-j)];j=j+1}
  }


  totalRr<-do.call(rbind,totalRr)
  J<-qr(totalRr)$rank

  ##for unit group
  if(grouppara==0){
    if(length(n)!=1){stop("n should be a number when grouppara=0")}
    if(is.list(Sigma)){stop("Sigma should be a matrix or number when grouppara=0")}
    if(nrow(as.matrix(Sigma))!=numSP||ncol(as.matrix(Sigma))!=numSP){
      stop("length of parameters and rank of covariance matrix 'Sigma' are not matched")
    }
    if(.checkcov(Sigma)==1){stop("the covariance matrix 'Sigma' you entered contains errors since it cannot exist")}

    b<-J/n
    thetacovpost<-Sigma
    thetacovprior<-thetacovpost/b
  }

  ##for multiple groups
  if(grouppara!=0){
    #if(length(n)==1){stop("n should be a vector when grouppara>0")}
    if(!is.list(Sigma)){stop("Sigma should be a list when grouppara>0")}
    if(any(unlist(lapply(Sigma,.checkcov))==1)){stop("the covariance matrix 'Sigma' you entered contains errors since it cannot exist")}

    dim_groups<-unlist(lapply(Sigma,dim))
    if(any(dim_groups!=mean(dim_groups))){
      stop("All covariance matrices in 'Sigma' should have the same dimension")
    }

    P<-length(Sigma)
    if(P!=length(n)){stop("length of sample size is not equal to the number of covariance matrices in 'Sigma'")}
    if(numSP!=grouppara*P+jointpara){stop("length of estimate is not correct for multiple groups")}
    if(any(dim_groups!=grouppara+jointpara)){stop("dim of of each covariance matrix in 'Sigma' should be grouppara+jointpara")}

    b<-rep(0,P)
    prior_cov<-list()
    for(p in 1:P){
      b[p]<-1/P*J/n[p]
      prior_cov[[p]]<-Sigma[[p]]/b[p]
    }

    inv_prior<-lapply(prior_cov,solve)
    inv_post<-lapply(Sigma,solve)

    thetacovprior<-.covmatrixfun(inv_prior,grouppara,jointpara,P)
    thetacovpost<-.covmatrixfun(inv_post,grouppara,jointpara,P)
  }

  #check about equality constraints or the comparability issue
  About=.Fortran("about",
                 as.integer(numH),
                 as.integer(numSP),
                 as.integer(numR),
                 totalRr=totalRr,
                 numARi=as.integer(rep(0,numH)),
                 error=as.integer(0)
  )

  totalRradjust=About$totalRr
  numAR=About$numARi
  error=About$error

  if(qr(totalRradjust)$rank > (qr(totalRradjust[1:sum(numR),1:numSP])$rank + sum(numAR))){
    error=2
  }

  for(i in 1:sum(numR)){
    for(j in 1:sum(numR)){
      if (all(totalRradjust[i,1:numSP]==totalRradjust[j,1:numSP])&&
          abs(totalRradjust[i,numSP+1]-totalRradjust[j,numSP+1])>0){
        error = 2
      }
      if (all(totalRradjust[i,1:numSP]==-totalRradjust[j,1:numSP])&&
          abs(totalRradjust[i,numSP+1]+totalRradjust[j,numSP+1])>0){
        error = 2
      }
    }
  }

  #Hypotheses are not comparable.
  if(error==1){stop("BaIn is not suited for the evaluation of one or more of the hypotheses,\n because the adjusted prior mean cannot be determined from R theta = r.")}
  if(error==2){stop("The informative hypotheses under evaluation are not comparable,\n and/or BaIn is not suited for the evaluation of one or more of the hypotheses,\n because the adjusted prior mean cannot be determined from R theta = r.")}

  fit=com=BF=rep(0,numH)
  fiteq=fitin=comeq=comin=rep(1,numH)
  results=rep(0,5*numH)

  for (h in 1:numH){
    ERr=IRr=constant=0
    if(numR[2*h-1]!=0){ERr<-matrix(totalRr[(sum(numR[1:(2*h-1)])-numR[2*h-1]+1):sum(numR[1:(2*h-1)]),1:(numSP+1)],numR[2*h-1],numSP+1)}
    if(numR[2*h]!=0){IRr<-matrix(totalRr[(sum(numR[1:(2*h)])-numR[2*h]+1):sum(numR[1:(2*h)]),1:(numSP+1)],numR[2*h],numSP+1);constant=IRr[,numSP+1]}

    #check redundant equality constraints.
    if(numR[2*h-1]!=0){
      if(qr(ERr)$rank<nrow(ERr)){
        stop("You have specified one or more hypotheses using equality constraints. However, for one or more of the hypotheses these equality constraints are redundant. You have to remove one or more of these constraints. Note: this will not alter your hypothesis. E.g. a=b and b=c and a=c is redundant. The same hypothesis is obtained using only a=b and b=c.")
      }
    }


    #compute the rowrank of IRr and the linear combiniation of independent constraints
    Mrank=.Fortran("mrank",
                   numR[2*h],
                   numSP,
                   rowrank=as.integer(0),
                   IRr=IRr,
                   transR=diag(0,numR[2*h],numR[2*h]),
                   constant,
                   transcon=rep(0,numR[2*h])
    )

    rowrank=Mrank$rowrank
    IRr=Mrank$IRr
    transR=Mrank$transR
    transcon=Mrank$transcon

    if(numR[2*h-1]==0){Rr=IRr}
    if(numR[2*h]==0){Rr=ERr}
    if(numR[2*h-1]!=0&&numR[2*h]!=0){Rr=rbind(ERr,IRr)}

    #parameter transformation for the estimates of theta
    thetar<-c(estimate,-1)
    betapost<-Rr[1:(numR[2*h-1]+rowrank),1:(numSP+1)]%*%thetar


    #parameter transformation for the covariance matrix of theta
    betacovpost<-Rr[1:(numR[2*h-1]+rowrank),1:numSP]%*%thetacovpost%*%t(matrix(Rr[1:(numR[2*h-1]+rowrank),1:numSP],nrow=numR[2*h-1]+rowrank,ncol=numSP))
    betacovpri<-Rr[1:(numR[2*h-1]+rowrank),1:numSP]%*%thetacovprior%*%t(matrix(Rr[1:(numR[2*h-1]+rowrank),1:numSP],nrow=numR[2*h-1]+rowrank,ncol=numSP))

    #specify prior mean
    betapri<-rep(0,numR[2*h-1]+rowrank)

    #adjust prior mean for about equality constraints
    if(numAR[h]>0){
      for(i in (numR[2*h-1]+1):(numR[2*h-1]+rowrank)){
        for(j in (sum(numR[1:(2*h)])-numR[2*h]+1):sum(numR[1:(2*h)])){
          if(all(Rr[i,1:numSP]==totalRradjust[j,1:numSP])&&Rr[i,numSP+1]!=totalRradjust[j,numSP+1])
          {betapri[i]=totalRradjust[j,numSP+1]-Rr[i,numSP+1]}
        }
      }
    }

    invbetadiagpost<-diag(solve(as.matrix(betacovpost)))
    invbetadiagpri<-diag(solve(as.matrix(betacovpri)))
    Bpost<-diag(1,numR[2*h-1]+rowrank)-solve(diag(invbetadiagpost,numR[2*h-1]+rowrank,numR[2*h-1]+rowrank))%*%solve(betacovpost)
    Bpri<-diag(1,numR[2*h-1]+rowrank)-solve(diag(invbetadiagpri,numR[2*h-1]+rowrank,numR[2*h-1]+rowrank))%*%solve(betacovpri)

    #equality constraints
    if (numR[2*h-1]>0){
      fiteq[h]<-1/sqrt((2*pi)^numR[2*h-1]*abs(det(as.matrix(betacovpost[1:numR[2*h-1],1:numR[2*h-1]]))))*exp(-1/2*(betapost[1:numR[2*h-1]]%*%solve(betacovpost[1:numR[2*h-1],1:numR[2*h-1]])%*%betapost[1:numR[2*h-1]]))
      comeq[h]<-1/sqrt((2*pi)^numR[2*h-1]*abs(det(as.matrix(betacovpri[1:numR[2*h-1],1:numR[2*h-1]]))))
    }


    #inequality constraints
    if (numR[2*h]>0){
      # function for the computation of complexity or fit for inequality constraints
      fitcom<-function(bet,invbetadiag,B,seed){
        forc=.Fortran("forc",
                      as.integer(numR[2*h-1]),
                      as.integer(numR[2*h]),
                      as.integer(rowrank),
                      bet,
                      transcon,
                      invbetadiag,
                      B,
                      transR,
                      f_or_c=as.double(0),
                      Numfc=as.integer(0),
                      as.integer(seed)
        )
        return(c(forc$f_or_c,forc$Numfc))
      }

      fitin[h]<-fitcom(betapost,invbetadiagpost,Bpost,seed)[1]
      numf<-fitcom(betapost,invbetadiagpost,Bpost,seed)[2]
      comin[h]<-fitcom(betapri,invbetadiagpri,Bpri,seed)[1]
      numc<-fitcom(betapri,invbetadiagpri,Bpri,seed)[2]
    }

    #total fit and complexity

    fit[h]<-fitin[h]*fiteq[h]
    com[h]<-comin[h]*comeq[h]

    #Bayes factor for a hypothesis vs its complement
    ifelse(numR[2*h-1]>0,BF[h]<-fit[h]/com[h],BF[h]<-(fit[h]/com[h])/((1-fit[h])/(1-com[h])))
  }


  fctable<-matrix(0,numH+1,9)
  PMPa<-c()
  PMPb<-c()

  for(h in 1:numH){
    PMPa[h]<-fit[h]/com[h]/(sum(fit/com))
    PMPb[h]<-fit[h]/com[h]/(1+sum(fit/com))
    fctable[h,]<-c(fiteq[h],fitin[h],comeq[h],comin[h],fit[h],com[h],BF[h],PMPa[h],PMPb[h])
  }
  fctable[numH+1,]<-c(rep(NA,8),1/(1+sum(fit/com)))
  fctable<-formatC(fctable, digits = 3, format = "f")

  fctable[which(fctable=="  NA",arr.ind=T)]<-"."
  fctable<-as.data.frame(fctable)
  rownames(fctable)<-c(paste("H",1:numH, sep=""),"Hu")
  colnames(fctable)<-c("f=","f>|=","c=","c>|=","f","c","BF.c","PMPa","PMPb")


  BFmatrix<-diag(1,numH)
  for(h1 in 1:numH){
    for(h2 in 1:numH){
      BFmatrix[h1,h2]<-fit[h1]/fit[h2]/(com[h1]/com[h2])
    }
  }
  rownames(BFmatrix)<-paste("H",1:numH, sep="")
  colnames(BFmatrix)<-paste("H",1:numH, sep="")


  res<-matrix(0,numH,5)
  colnames(res)<-c("fit","complexity","BF","PMPa","PMPb")
  rownames(res)<-paste("H",1:numH,sep="")

  for(h in 1:numH){
    res[h,]<-c(fit[h],com[h],BF[h],fit[h]/com[h]/sum(fit/com),fit[h]/com[h]/(1+sum(fit/com)))
  }

  if(print){
    #print result
    writeLines("Choice of b")
    writeLines(paste("J",J))
    writeLines(c("N",n),sep=" ")
    writeLines("",sep="\n")
    writeLines(c("b",formatC(b, digits = 3, format = "f")),sep=" ")
    writeLines(" ")
    cat("Estimates and covariance matrix of parameters","Estimates",sep="\n")
    writeLines(paste((formatC(estimate, digits = 3, format = "f"))),sep = " ")
    writeLines(" ")
    cat("Posterior Covariance Matrix","\n")
    write.table(data.frame(formatC(thetacovpost, digits = 3, format = "f")),col.names = FALSE,row.names = FALSE,quote = FALSE)
    cat("Prior Covariance Matrix","\n")
    write.table(data.frame(formatC(thetacovprior, digits = 3, format = "f")),col.names = FALSE,row.names = FALSE,quote = FALSE)
    writeLines(" ")
    cat("Hypothesis testing result", sep="\n")
    write.table(capture.output(fctable),col.names = FALSE,row.names = FALSE,quote = FALSE)
    writeLines(" ")
    cat("BF-matrix", sep="\n")
    write.table(capture.output(data.frame(formatC(BFmatrix, digits = 3, format = "f"))),col.names = FALSE,row.names = FALSE,quote = FALSE)
    writeLines(" ")

  }

  cl<-match.call()

  Bainres<-list(b=b,priorCov=thetacovprior, posterCov=thetacovpost,
                testResult=as.data.frame(res),
                fit_com_table=as.data.frame(fctable),
                BFmatrix=as.data.frame(BFmatrix, row.names=attributes(BFmatrix)$row.names))
  class(Bainres)<-"Bain"
  Bainres$call<-cl
  return(Bainres)
}


.checkcov<-function(Sigma){
  error=0
  if(!isTRUE(all.equal(as.matrix(Sigma),t(Sigma),tolerance=1e-10))){error=1}else{
    if(any(eigen(Sigma)$values<=0)){error=1}}
  return(error)
}


.covmatrixfun<-function(inv_cov_list,grouppara,jointpara,P){

  inv_upperleft<-lapply(inv_cov_list,function(x) x[1:grouppara,1:grouppara])
  if(jointpara>0){
    inv_upperright<-lapply(inv_cov_list,function(x) matrix(x[1:grouppara,(grouppara+1):(grouppara+jointpara)],grouppara,jointpara))
    inv_lowerleft<-lapply(inv_cov_list,function(x) matrix(x[(grouppara+1):(grouppara+jointpara),1:grouppara],jointpara,grouppara))
    inv_lowerright<-lapply(inv_cov_list,function(x) matrix(x[(grouppara+1):(grouppara+jointpara),(grouppara+1):(grouppara+jointpara)],jointpara,jointpara))
    inv_lowerright_matrix<-diag(0,jointpara)
  }

  inv_cov_total<-diag(0,P*grouppara)

  #or using cbind and rbind
  for(p in 1:P){
    inv_cov_total[((p-1)*grouppara+1):(p*grouppara),((p-1)*grouppara+1):(p*grouppara)]<-inv_upperleft[[p]]
    if(jointpara>0){
      inv_lowerright_matrix<-inv_lowerright_matrix+inv_lowerright[[p]]
    }
  }

  if(jointpara>0){
    inv_cov_total<-cbind(inv_cov_total,do.call(rbind,inv_upperright))
    inv_cov_total<-rbind(inv_cov_total,cbind(do.call(cbind,inv_lowerleft),inv_lowerright_matrix))
  }

  covmatrix<-solve(inv_cov_total)
  return(covmatrix)
}

