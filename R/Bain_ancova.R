Bain_ancova<-function(X, dep_var, covariates, group, ERr=NULL,IRr=NULL,...){
  TRr<-list(ERr,IRr,...)

  if(!is.data.frame(X)){stop("X should be a data frame")}  #X data frame
  Xnames<-names(X)

  if(!(is.character(dep_var)&&is.character(covariates)&&is.character(group)) ){
    stop("dep_var, covariates, and group should be character")
  }
  if(length(dep_var)!=1){stop("Please specify only one dependent variable")}
  if(length(group)!=1){stop("Please specify only one group factor")}

  var_names<-c(dep_var,covariates,group)
  if(any(!var_names%in%Xnames)){stop("one of names of dep_var, covariates, and group is not in X")}

  depv<-X[,dep_var]  ##dependent variable
  groupf<-factor(X[,group]) ## group factor
  covars<-apply(as.matrix(X[,covariates]),2,function(x) x-mean(x))  ##standardized covariates
  n_covars=length(covariates)

  n<-unlist(lapply(split(depv,groupf),length)) ##sample size per group

  ##analysis
  ancovafm <-  lm(depv ~ groupf + covars -1)
  estimate<-coef(ancovafm)
  resvar<-(summary(ancovafm)$sigma)**2

  newdata<-cbind(rep(1,nrow(covars)),covars) # IMPORTANT NOW FOR EACH GROUP THE INTERCEPT IS ADDED

  split_index<-split(1:nrow(covars),groupf)
  split_newdata<-lapply(split_index, function(x) newdata[x,])
  covariance<-lapply(split_newdata, function(x) resvar * solve(t(x) %*% x))

  ##check
  if(any(unlist(lapply(TRr,ncol))!=(length(levels(groupf))+1)&unlist(lapply(TRr,ncol))!=0)){
    stop("Number of columns in ERr or IRr should equal number of groups plus 1.")
  }

  ##hypothesis
  newTRr<- lapply(TRr,.add_col,n_cov=n_covars)## add extra columns for covariates' coeffs.
  IR_chara<-paste0(",newTRr[[",1:length(newTRr),"]]",sep="",collapse = "")


  jointpara<-n_covars
  ##check
  if(jointpara!=(length(estimate)-length(n))){stop("length(n) + n_covars != length(estimate)")}

  Bain_chara<-paste0("Bain(estimate,covariance,grouppara=1,jointpara=jointpara,n=n",IR_chara,",seed=100,print=FALSE)")
  Bain_res<-eval(parse(text = Bain_chara))
  BFmatrix<-as.matrix(Bain_res$BFmatrix)


  cat("ANCOVA test result", sep="\n")
  write.table(capture.output(Bain_res$fit_com_table[,5:9]),col.names = FALSE,row.names = FALSE,quote = FALSE)

  writeLines(" ")
  cat("BF-matrix", sep="\n")
  write.table(capture.output(data.frame(formatC(BFmatrix, digits = 3, format = "f"))),col.names = FALSE,row.names = FALSE,quote = FALSE)
  writeLines(" ")

  ancovatest_res<-Bain_res$testResult

  ancova_res<-list(fit = ancovatest_res$fit, complexity = ancovatest_res$complexity,
                   BF = ancovatest_res$BF, PMPa = ancovatest_res$PMPa, PMPb = ancovatest_res$PMPb,
                   estimate_res = ancovafm)

  cl<-match.call()
  class(ancova_res)<-"Bain"
  ancova_res$call<-cl

  return(invisible(ancova_res))
}



.add_col<-function(x,n_cov){
  if(length(x)!=0) {
    if(nrow(x)!=1){
      x<-cbind(x[,-ncol(x)],matrix(0,nrow(x),n_cov),x[,ncol(x)])
    }else{
      x<-matrix(c(x[,-ncol(x)],rep(0,n_cov),x[,ncol(x)]),nrow=1)
    }
  }
  return(x)
}


