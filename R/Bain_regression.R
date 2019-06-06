Bain_regression<-function(formula, data, ERr = NULL, IRr = NULL,...,covariates_hypo = NULL, standardize = FALSE){
  TRr<-list(ERr,IRr,...)
  IR_chara<-paste0("TRr[[",1:length(TRr),"]]",sep=",",collapse = "")

  dependent<-model.frame(formula,data)[,1]
  predictor<-model.frame(formula,data)[,-1]
  predictor_names<-names(model.frame(formula,data))[-1]
  n<-length(dependent)

  if(is.null(covariates_hypo)){ ##If null, assume all predictors are used in the hypotheses
    covariates_hypo<-names(model.frame(formula,data))[-1]
  }

  if(inherits(covariates_hypo,"formula")){  ## if formula
    covariates_hypo<-names(model.frame(covariates_hypo,data))
  }

  if(!is.character(covariates_hypo)){stop("'covariates_hypo' should be a character (vector), formula, or NULL")}
  if(!all(covariates_hypo %in% predictor_names)){
    stop("'covariates_hypo' should be a subset of predictors in 'formula'")
  }

  jointpara<-length(covariates_hypo)

  ###standardize
  if(!standardize){
    fit<-lm(formula,data = data)
    estimate <- coef(fit)[covariates_hypo]
    covariance <- vcov(fit)[covariates_hypo,covariates_hypo]
    reg.fit<-fit
  }else{
    sink(tempfile())
    ##Function seBeta() is from fungible package (version 1.5).
    ##It computes standardized estimates of coefficients.
    intermed <- seBeta(predictor, dependent, Nobs = n, alpha = .05, estimator = 'Normal')
    sink()
    covariates_hypo <- which(predictor_names %in% covariates_hypo)
    estimate <- intermed$CIs$estimate[covariates_hypo]
    covariance <- intermed$cov.mat[covariates_hypo,covariates_hypo]
    reg.fit<-intermed
  }

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

