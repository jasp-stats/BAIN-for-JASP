Bain_regression_cm<-function(formula, data, hyp, standardize = FALSE){

  dependent<-model.frame(formula,data)[,1]
  predictor<-model.frame(formula,data)[,-1]
  predictor_names<-names(model.frame(formula,data))[-1]
  n<-length(dependent)

  split_temp<-unlist(strsplit(hyp,"[=<>+-]|[\\*\\/\\&\\(\\)\\{\\}\\,\\;]"))
  split_temp<-split_temp[which(grepl("[a-zA-Z_]",split_temp))]
  split_temp<-unique(gsub("\\s+", "",split_temp, perl = TRUE))   ##remove blank; unique

  if(all(split_temp%in%predictor_names)){
    varnames<-predictor_names[which(predictor_names%in%split_temp)]
  }else{
    stop("'hyp' contains invalid parameter(s), which does not match predictor names in 'formula'")
  }

  jointpara<-length(varnames)

  ###standardize
  if(!standardize){
    fit<-lm(formula,data = data)
    estimate <- coef(fit)[varnames]
    covariance <- vcov(fit)[varnames,varnames]
    reg.fit<-fit
  }else{
    sink(tempfile())
    ##Function seBeta() is from fungible package (version 1.5).
    ##It computes standardized estimates of coefficients.
    intermed <- seBeta(predictor, dependent, Nobs = n, alpha = .05, estimator = 'Normal')
    sink()
    var_num <- which(predictor_names %in% varnames)
    estimate <- intermed$CIs$estimate[var_num]
    covariance <- intermed$cov.mat[var_num,var_num]
    reg.fit<-intermed
  }

  TRr<-create_matrices(varnames,hyp)
  IR_chara<-paste0("TRr[[",1:length(TRr),"]]",sep=",",collapse = "")

  Bain_chara<-paste0("Bain(estimate,covariance,grouppara=0,jointpara=jointpara,n=n,",IR_chara,"seed=100,print=FALSE)")
  Bain_res<-eval(parse(text = Bain_chara))
  BFmatrix<-as.matrix(Bain_res$BFmatrix)

  cl<-match.call()

  cat("Regression test result", sep="\n")
  write.table(capture.output(Bain_res$fit_com_table[,5:9]),col.names = FALSE,row.names = FALSE,quote = FALSE)

  writeLines(" ")
  cat("BF-matrix", sep="\n")
  write.table(capture.output(data.frame(formatC(BFmatrix, digits = 3, format = "f"))),col.names = FALSE,row.names = FALSE,quote = FALSE)
  writeLines(" ")

  Regtest_res<-Bain_res$testResult
  #reg_res<-list(test_res = Regtest_res, estimate_res = reg.fit)
  reg_res<-list(fit = Regtest_res$fit, complexity = Regtest_res$complexity,
                BF = Regtest_res$BF, PMPa = Regtest_res$PMPa, PMPb = Regtest_res$PMPb,
                estimate_res = reg.fit)


  class(reg_res)<-"Bain"
  reg_res$call<-cl
  return(invisible(reg_res))
}

