Bain_ttestData<-function(x,y=NULL,nu=0,type=1,paired=FALSE){
  if(is.null(y)&&!paired){
    estimate<-mean(x)
    n<-length(x)
    variance<-(sd(x))^2/n
    ttDatares<-Bain_ttest(estimate=estimate,variance=variance,n=n,nu=nu,type=type)
  }else if(is.null(y)&&paired){stop("paired test requires two groups")
  }else if(!is.null(y)&&!paired){
    estimate<-c(mean(x),mean(y))
    n1<-length(x)
    n2<-length(y)
    variance<-c((sd(x))^2/n1,(sd(y))^2/n2)
    ttDatares<-Bain_ttest(estimate=estimate,variance=variance,n=c(n1,n2),nu=nu,type=type)
  }else if(!is.null(y)&&paired){
    if(length(x)!=length(y)){stop("paired test requires equal size of two groups")}
    df<-x-y
    estimate<-mean(df)
    n<-length(df)
    variance<-(sd(df))^2/n
    ttDatares<-Bain_ttest(estimate=estimate,variance=variance,n=n,nu=nu,type=type)
  }

  cl<-match.call()
  class(ttDatares)<-"Bain"
  ttDatares$call<-cl
  return(invisible(ttDatares))
}


