#' This  function is used estimate the equal valued p(Now) point and the delay and RT associated with a particular stimulus in a delayed discountig procedure
#'
#' Function to estimate the equal valued p(Now) point and the delay and RT associated with a particular stimulus in a delayed discountig procedure
#' @param QCETrialStructureList A list that specifies how the trials will be presented in the experiment.  This list specifies the selection of stimuli from stimFile.json, the ordering of stimuli, the blocking structure, etc.  If you are building a new list, then this should be NULL. If you are adding a new effect to an old list, then this should be the QCETrialStructureList that you are adding an effect to. DEFAULT = NULL
#' @param data A dataframe that contains the relevant columns.
#' @param delayColumn A string that specifies the name of the column in "data" that contains the delays.
#' @param fitColumn A string that specifies the name of the column in "data" that contains the best fit values to the RT data for each delay.
#' @param pStandardColumn A string that specifies the name of the column in "data" that contains the probability of choosing the immediate option for each delay.
#' @param unbiased An boolean that specifies whether this participant has a start point bias or not.  If yes, set to FALSE, if not set to TRUE.  DEFAULT = FALSE.
#' @param intervals An integer that specifies the number of intervals to impute. DEFAULT = 1000.
#''
#' @return the a dataframe containing the following columns: Delay; pStandard (the equal valued p(Now) point); Fit (the best guess RT for this Delay and probe)
#' @keywords delay discounting equal value p(Now) point
#' @export
#' @examples getEqualValueDelay (df.tmp, "Delay", "Fit", "pStandard", unbiased = F)

library(pracma)

# in this code pStandard is the probablity of responding now: p(now)
#return estimated delay from fit function
getEqualValueDelay <- function(data, delayColumn, fitColumn, pStandardColumn, unbiased = F, intervals = 1000) {
  data <- data[order(data[[delayColumn]]),]
  x.n <- data[[delayColumn]]
  fit <- data[[fitColumn]]
  pS <- df.tmp[[pStandardColumn]]

  approx_val.fit <- approx(x.n,fit,  n=intervals)
  approx_val.pS <- approx(x.n,pS,  approx_val.fit$x)

  indif.df.pS <- data.frame(Delay = approx_val.pS$x,pStandard = approx_val.pS$y)
  indif.df.fit <- data.frame(Delay = approx_val.fit$x,Fit = approx_val.fit$y)

  indif.df <- merge (indif.df.pS, indif.df.fit)

  if(unbiased == F) {
    df.tmp.fit <- indif.df[order(indif.df$Fit),]
    df.tmp.fit <- df.tmp.fit[complete.cases(df.tmp.fit),]
    df.tmp.fit.max <- df.tmp.fit[nrow(df.tmp.fit), ]
  } else {
    df.tmp.fit.max <- indif.df[which.min(abs(indif.df$pStandard - 0.50)),]
  }

  return(df.tmp.fit.max)
}
