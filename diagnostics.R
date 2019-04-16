### Run Diagnostics ###
### runs test to check 
### data distribution
### normality, etc. 


run_diagnostics = function(df) {
  
  if (class(df) != "lm") {
    print("This is the Variance information:")
    print(VarCorr(df))
  }
  
  print("This is the data fit information (AIC, BIC, loglikelihood, deviance)")
  tmp <- c(AIC(df), BIC(df), logLik(df), deviance(df))
  print(tmp)
  
  par(mfrow=c(2,2))
  tmp <- cooks.distance(df)
  
  id <- tmp <= mean(tmp)*4 # find a principled reason here
  plot(tmp, xlab = "Values", ylab = "Cooks Distance", main = "Potential Outliers")
  abline(h = 4*mean(tmp, na.rm=T), lty = 2, col="red")
  text(x=1:length(tmp) + 1, y= tmp, labels=ifelse(tmp>4*mean(tmp, na.rm=T),names(tmp),""), col="blue") 
  hist((resid(df) - mean(resid(df))) / sd(resid(df)), freq = FALSE, xlab = "Residual Values", main = "Hist. of Residuals"); curve(dnorm, add = TRUE)
  
  qqnorm(residuals(df)); qqline(residuals(df))
  plot(fitted(df), residuals(df), xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs. Fitted Values");
  abline(h = 0, lty = 3, col = "gray"); abline(v = 0, lty = 3, col = "gray"); lines(smooth.spline(fitted(df), residuals(df)), col="red")
  
  plot(model.response(model.frame(df)) ~ fitted(df), xlab = "Actual", ylab = "Fitted", main = "Correlation between fitted and actual data")
  
  print(corr.test(model.response(model.frame(df)), fitted(df)))
  
  plot(df,type=c("p","smooth"), xlab = "Fitted values", ylab = "Standardized residuals")
  par(ask=TRUE)
  plot(df, sqrt(abs(resid(.)))~fitted(.), xlab = "Fitted values", ylab = "sqrt(residuals)")
}
