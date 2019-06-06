#t test
#one sample t test

Bain_ttest<-function(estimate,variance,n,nu=0,type=1){
  if(type==1&&length(estimate)==1){  # mu=nu vs mu!=nu
    Bain_res<-Bain(estimate=estimate,
                   Sigma=variance,
                   grouppara=0,
                   jointpara=1,
                   n=n,
                   ERr=matrix(c(1,nu),1,2),
                   seed=100,print=FALSE)

    res<-Bain_res$testResult["H1",c("BF","PMPb")]
    res<-cbind(res,(1-res[2]))
    names(res)<-c("BF_0u","PMP_0","PMP_u")
    cat(paste("Hypotheses H0: mu = ",nu," vs Hu: mu != ",nu,sep=""),sep="\n")
  }

  if(type==2&&length(estimate)==1){ # mu=nu vs mu>nu
    Bain_res<-Bain(estimate=estimate,
                   Sigma=variance,
                   grouppara=0,
                   jointpara=1,
                   n=n,
                   matrix(c(1,nu),1,2),
                   matrix(0,0,0),
                   matrix(0,0,0),
                   matrix(c(1,nu),1,2),
                   seed=100,print=FALSE)

    fit1<-Bain_res$testResult["H1","fit"]
    fit2<-Bain_res$testResult["H2","fit"]
    com1<-Bain_res$testResult["H1","complexity"]
    com2<-Bain_res$testResult["H2","complexity"]
    #BF1<-Bain_res$testResult["H1","BF"]
    #BF2<-Bain_res$testResult["H2","BF"]
    BF<-fit1/com1/(fit2/com2)
    PMP1<-fit1/com1/(fit1/com1+fit2/com2)
    PMP2<-1-PMP1

    res<-t(c(BF,PMP1,PMP2))
    colnames(res)<-c("BF_01","PMP_0","PMP_1")
    cat(paste("Hypotheses H0: mu = ",nu," vs H1: mu > ",nu,sep=""),sep="\n")

  }

  if(type==3&&length(estimate)==1){# mu=nu vs mu<nu
    Bain_res<-Bain(estimate=estimate,
                   Sigma=variance,
                   grouppara=0,
                   jointpara=1,
                   n=n,
                   matrix(c(1,nu),1,2),
                   matrix(0,0,0),
                   matrix(0,0,0),
                   matrix(c(-1,-nu),1,2),
                   seed=100,print=FALSE)

    fit1<-Bain_res$testResult["H1","fit"]
    fit2<-Bain_res$testResult["H2","fit"]
    com1<-Bain_res$testResult["H1","complexity"]
    com2<-Bain_res$testResult["H2","complexity"]
    #BF1<-Bain_res$testResult["H1","BF"]
    #BF2<-Bain_res$testResult["H2","BF"]
    BF<-fit1/com1/(fit2/com2)
    PMP1<-fit1/com1/(fit1/com1+fit2/com2)
    PMP2<-1-PMP1

    res<-t(c(BF,PMP1,PMP2))
    colnames(res)<-c("BF_01","PMP_0","PMP_1")
    cat(paste("Hypotheses H0: mu = ",nu," vs H1: mu < ",nu,sep=""),sep="\n")

  }

  if(type==4&&length(estimate)==1){# mu>nu vs mu<nu
    Bain_res<-Bain(estimate=estimate,
                   Sigma=variance,
                   grouppara=0,
                   jointpara=1,
                   n=n,
                   matrix(0,0,0),
                   matrix(c(1,nu),1,2),
                   matrix(0,0,0),
                   matrix(c(-1,-nu),1,2),
                   seed=100,print=FALSE)

    fit1<-Bain_res$testResult["H1","fit"]
    fit2<-Bain_res$testResult["H2","fit"]
    com1<-Bain_res$testResult["H1","complexity"]
    com2<-Bain_res$testResult["H2","complexity"]
    #BF1<-Bain_res$testResult["H1","BF"]
    #BF2<-Bain_res$testResult["H2","BF"]
    BF<-fit1/com1/(fit2/com2)
    PMP1<-fit1/com1/(fit1/com1+fit2/com2)
    PMP2<-1-PMP1

    res<-t(c(BF,PMP1,PMP2))
    colnames(res)<-c("BF_12","PMP_1","PMP_2")
    cat(paste("Hypotheses H1: mu > ",nu," vs H2: mu < ",nu,sep=""),sep="\n")
  }

  if(type==5&&length(estimate)==1){# mu=nu vs mu>nu vs mu<nu
    Bain_res<-Bain(estimate=estimate,
                   Sigma=variance,
                   grouppara=0,
                   jointpara=1,
                   n=n,
                   matrix(c(1,nu),1,2),
                   matrix(0,0,0),
                   matrix(0,0,0),
                   matrix(c(1,nu),1,2),
                   matrix(0,0,0),
                   matrix(c(-1,-nu),1,2),
                   seed=100,print=FALSE)

    fit1<-Bain_res$testResult["H1","fit"]
    fit2<-Bain_res$testResult["H2","fit"]
    fit3<-Bain_res$testResult["H3","fit"]
    com1<-Bain_res$testResult["H1","complexity"]
    com2<-Bain_res$testResult["H2","complexity"]
    com3<-Bain_res$testResult["H3","complexity"]

    BF12<-fit1/com1/(fit2/com2)
    BF13<-fit1/com1/(fit3/com3)
    BF23<-fit2/com2/(fit3/com3)
    PMP1<-fit1/com1/(fit1/com1+fit2/com2+fit3/com3)
    PMP2<-fit2/com2/(fit1/com1+fit2/com2+fit3/com3)
    PMP3<-fit3/com3/(fit1/com1+fit2/com2+fit3/com3)

    res<-t(c(BF12,BF13,BF23,PMP1,PMP2,PMP3))
    colnames(res)<-c("BF_01","BF_02","BF_12","PMP_0","PMP_1","PMP_2")
    cat(paste("Hypotheses H0: mu = ",nu," vs H1: mu > ",nu," vs H2: mu < ",nu,sep=""),sep="\n")
  }


  if(type==1&&length(estimate)==2){  # mu1=mu2 vs mu1!=mu2   or mu1-mu2=nu
    Bain_res<-Bain(estimate=estimate,
                   Sigma=list(as.matrix(variance[1]),as.matrix(variance[2])),
                   grouppara=1,
                   jointpara=0,
                   n=n,
                   ERr=matrix(c(1,-1,nu),1,3),
                   seed=100,print=FALSE)

    res<-Bain_res$testResult["H1",c("BF","PMPb")]
    res<-cbind(res,(1-res[2]))
    names(res)<-c("BF_0u","PMP_0","PMP_u")
    cat("Hypotheses H0: mu1 = mu2 vs Hu: mu1 != mu2",sep="\n")

  }

  if(type==2&&length(estimate)==2){ # mu1=mu2 vs mu1>mu2
    Bain_res<-Bain(estimate=estimate,
                   Sigma=list(as.matrix(variance[1]),as.matrix(variance[2])),
                   grouppara=1,
                   jointpara=0,
                   n=n,
                   matrix(c(1,-1,nu),1,3),
                   matrix(0,0,0),
                   matrix(0,0,0),
                   matrix(c(1,-1,nu),1,3),
                   seed=100,print=FALSE)

    fit1<-Bain_res$testResult["H1","fit"]
    fit2<-Bain_res$testResult["H2","fit"]
    com1<-Bain_res$testResult["H1","complexity"]
    com2<-Bain_res$testResult["H2","complexity"]
    #BF1<-Bain_res$testResult["H1","BF"]
    #BF2<-Bain_res$testResult["H2","BF"]
    BF<-fit1/com1/(fit2/com2)
    PMP1<-fit1/com1/(fit1/com1+fit2/com2)
    PMP2<-1-PMP1

    res<-t(c(BF,PMP1,PMP2))
    colnames(res)<-c("BF_01","PMP_0","PMP_1")
    cat("Hypotheses H0: mu1 = mu2 vs H1: mu1 > mu2",sep="\n")

  }

  if(type==3&&length(estimate)==2){# mu1=mu2 vs mu1<mu2
    Bain_res<-Bain(estimate=estimate,
                   Sigma=list(as.matrix(variance[1]),as.matrix(variance[2])),
                   grouppara=1,
                   jointpara=0,
                   n=n,
                   matrix(c(1,-1,nu),1,3),
                   matrix(0,0,0),
                   matrix(0,0,0),
                   matrix(c(-1,1,-nu),1,3),
                   seed=100,print=FALSE)

    fit1<-Bain_res$testResult["H1","fit"]
    fit2<-Bain_res$testResult["H2","fit"]
    com1<-Bain_res$testResult["H1","complexity"]
    com2<-Bain_res$testResult["H2","complexity"]
    #BF1<-Bain_res$testResult["H1","BF"]
    #BF2<-Bain_res$testResult["H2","BF"]
    BF<-fit1/com1/(fit2/com2)
    PMP1<-fit1/com1/(fit1/com1+fit2/com2)
    PMP2<-1-PMP1

    res<-t(c(BF,PMP1,PMP2))
    colnames(res)<-c("BF_01","PMP_0","PMP_1")

    cat("Hypotheses H0: mu1 = mu2 vs H1: mu1 < mu2",sep="\n")

  }

  if(type==4&&length(estimate)==2){# mu1>mu2 vs mu1<mu2
    Bain_res<-Bain(estimate=estimate,
                   Sigma=list(as.matrix(variance[1]),as.matrix(variance[2])),
                   grouppara=1,
                   jointpara=0,
                   n=n,
                   matrix(0,0,0),
                   matrix(c(1,-1,nu),1,3),
                   matrix(0,0,0),
                   matrix(c(-1,1,-nu),1,3),
                   seed=100,print=FALSE)

    fit1<-Bain_res$testResult["H1","fit"]
    fit2<-Bain_res$testResult["H2","fit"]
    com1<-Bain_res$testResult["H1","complexity"]
    com2<-Bain_res$testResult["H2","complexity"]
    #BF1<-Bain_res$testResult["H1","BF"]
    #BF2<-Bain_res$testResult["H2","BF"]
    BF<-fit1/com1/(fit2/com2)
    PMP1<-fit1/com1/(fit1/com1+fit2/com2)
    PMP2<-1-PMP1

    res<-t(c(BF,PMP1,PMP2))
    colnames(res)<-c("BF_12","PMP_1","PMP_2")

    cat("Hypotheses H1: mu1 > mu2 vs H2: mu1 < mu2",sep="\n")
  }

  if(type==5&&length(estimate)==2){# mu1=mu2 vs mu1>mu2 vs mu1<mu2
    Bain_res<-Bain(estimate=estimate,
                   Sigma=list(as.matrix(variance[1]),as.matrix(variance[2])),
                   grouppara=1,
                   jointpara=0,
                   n=n,
                   matrix(c(1,-1,nu),1,3),
                   matrix(0,0,0),
                   matrix(0,0,0),
                   matrix(c(1,-1,nu),1,3),
                   matrix(0,0,0),
                   matrix(c(-1,1,-nu),1,3),
                   seed=100,print=FALSE)

    fit1<-Bain_res$testResult["H1","fit"]
    fit2<-Bain_res$testResult["H2","fit"]
    fit3<-Bain_res$testResult["H3","fit"]
    com1<-Bain_res$testResult["H1","complexity"]
    com2<-Bain_res$testResult["H2","complexity"]
    com3<-Bain_res$testResult["H3","complexity"]
    BF12<-fit1/com1/(fit2/com2)
    BF13<-fit1/com1/(fit3/com3)
    BF23<-fit2/com2/(fit3/com3)
    PMP1<-fit1/com1/(fit1/com1+fit2/com2+fit3/com3)
    PMP2<-fit2/com2/(fit1/com1+fit2/com2+fit3/com3)
    PMP3<-fit3/com3/(fit1/com1+fit2/com2+fit3/com3)

    res<-t(c(BF12,BF13,BF23,PMP1,PMP2,PMP3))
    colnames(res)<-c("BF_01","BF_02","BF_12","PMP_0","PMP_1","PMP_2")
    cat("Hypotheses H0: mu1 = mu2 vs H1: mu1 > mu2 vs H2: mu1 < mu2" ,sep="\n")

  }


  res.x<-data.frame(formatC(as.matrix(res),digits = 3, format = "f"))
  rownames(res.x)<-rownames(res)<-""

  cl<-match.call()

  writeLines(" ")
  cat("t test result", sep="\n")
  write.table(capture.output(res.x),col.names = FALSE,row.names = FALSE,quote = FALSE)

  ttest_res<-as.data.frame(res)
  class(ttest_res)<-"Bain"
  ttest_res$call<-cl
  return(invisible(ttest_res))
}




