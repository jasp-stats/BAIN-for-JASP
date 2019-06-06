Bain_ancova_cm<-function(X, dep_var, covariates, group, hyp){

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

  ##possible variable names
  varnames_g<-paste("g",1:length(n),sep = "")
  varnames_gd<-paste("g",1:length(n),sep = ".")
  varnames_gn<-paste(group,1:length(n),sep = ".")
  varnames_gf<-paste(group,levels(groupf),sep = ".")

  split_temp<-unlist(strsplit(hyp,"[=<>+-]|[\\*\\/\\&\\(\\)\\{\\}\\,\\;]"))
  split_temp<-split_temp[which(grepl("[a-zA-Z]",split_temp))]
  split_temp<-unique(gsub("\\s+", "",split_temp, perl = TRUE))   ##remove blank; unique

  if(all(split_temp%in%varnames_g)){
    varnames<-varnames_g
  }else if(all(split_temp%in%varnames_gd)){
    varnames<-varnames_gd
  }else if(all(split_temp%in%varnames_gn)){
    varnames<-varnames_gn
  }else if(all(split_temp%in%varnames_gf)){
    varnames<-varnames_gf
  }else{
    stop("'hyp' contains invalid parameter(s), which does not match group factor")
  }

  TRr<-create_matrices(c(varnames,covariates),hyp)

  ##check
  if(any(unlist(lapply(TRr,ncol))!=(length(levels(groupf))+n_covars+1)&unlist(lapply(TRr,ncol))!=0)){
    stop("Number of columns in ERr or IRr should equal number of parameters plus 1.")
  }

  ##hypothesis
  IR_chara<-paste0(",TRr[[",1:length(TRr),"]]",sep="",collapse = "")


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


