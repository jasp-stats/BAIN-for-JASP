print.Bain <- function(x,...){
  fun_type<-as.character(x$call)[1]
  # Inputs a bain result object
  if(fun_type=="Bain"){
    print(x$fit_com_table)
  }

  # Inputs a bain t-test result object
  if(fun_type=="Bain_ttest"||fun_type=="Bain_ttestData"){
    output<-do.call(cbind,x[-length(x)])
    print(data.frame(output,row.names = ""))
  }

  # Inputs a bain anova or ancova test result object
  if(fun_type=="Bain_anova"||fun_type=="Bain_ancova"||fun_type=="Bain_regression"
     ||fun_type=="Bain_anova_cm"||fun_type=="Bain_ancova_cm"||fun_type=="Bain_regression_cm"){
    output<-data.frame(do.call(cbind,x[-c(length(x)-1,length(x))]))
    names(output)<-c("f","c","BF.c","PMPa","PMPb")
    rownames(output)<-paste("H",1:nrow(output),sep = "")
    print(output)
  }

}

