Bain_anova<-function(X, dep_var=NULL, group=NULL, ERr=NULL,IRr=NULL, ...){
  TRr<-list(ERr,IRr,...)
  IR_chara<-paste0("TRr[[",1:length(TRr),"]]",sep=",",collapse = "")

  if(!is.data.frame(X)){stop("X should be a data frame")}  #X data frame

  if(is.character(dep_var)){
    if(length(dep_var)!=1){stop("Please specify only one dependent variable")}
    temp<-which(names(X)==dep_var)
    if(length(temp)!=1){stop("No variable or two or more variables has name dep_var in X")}
    depv <- X[,temp]         ###dependent variable
  }else if(is.null(dep_var)){
    depv <- X[,1]
  }else{
    stop("dep_var should be a character or NULL")
  }

  if(is.character(group)){
    temp<-which(names(X)==group)
    if(length(temp)!=1){stop("No variable or two or more variables has name group in X")}
    groupf<- factor(X[,temp])   ###dgroup factor
  }else if(is.null(group)){
    groupf <- factor(X[,2])
  }else{
    stop("group should be a character or NULL")
  }

  n<-unlist(lapply(split(depv,groupf),length))

  anovafm <-  lm(depv ~ groupf -1)
  estimate<-coef(anovafm)
  variance <- (summary(anovafm)$sigma)**2
  variance <- 1/n * variance
  covlist<- lapply(as.list(variance),matrix) ##convert variances to list of variances used in Bain

  Bain_chara<-paste0("Bain(estimate,covlist,grouppara=1,jointpara=0,n=n,",IR_chara,"seed=100,print=FALSE)")
  Bain_res<-eval(parse(text = Bain_chara))
  BFmatrix<-as.matrix(Bain_res$BFmatrix)


  cat("ANOVA test result", sep="\n")
  write.table(capture.output(Bain_res$fit_com_table[,5:9]),col.names = FALSE,row.names = FALSE,quote = FALSE)

  writeLines(" ")
  cat("BF-matrix", sep="\n")
  write.table(capture.output(data.frame(formatC(BFmatrix, digits = 3, format = "f"))),col.names = FALSE,row.names = FALSE,quote = FALSE)
  writeLines(" ")

  anovatest_res<-Bain_res$testResult

  anova_res<-list(fit = anovatest_res$fit, complexity = anovatest_res$complexity,
                  BF = anovatest_res$BF, PMPa = anovatest_res$PMPa, PMPb = anovatest_res$PMPb,
                  estimate_res = anovafm)

  cl<-match.call()
  class(anova_res)<-"Bain"
  anova_res$call<-cl
  return(invisible(anova_res))

}


