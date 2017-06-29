fixNAcoefficients <- function(glmModel) {
  factor_vars <- names(glmModel$xlevels)
  touchedFactors <- FALSE
  for (factor_var in factor_vars) {
    xlevels <- grep(pattern = factor_var, x = names(coef(glmModel)), value = TRUE)
    if (!is.na(summary(is.na(glmModel$coefficients[xlevels]))["TRUE"])) {
      touchedFactors <- TRUE
  	  if (length(xlevels) > 1) {
        glmModel$coefficients["(Intercept)"] <- glmModel$coefficients["(Intercept)"] + mean(glmModel$coefficients[xlevels], na.rm=TRUE)
        glmModel$coefficients[xlevels][is.na(glmModel$coefficients[xlevels])] <- mean(glmModel$coefficients[xlevels], na.rm = TRUE)
        glmModel$coefficients[xlevels] <- scale(glmModel$coefficients[xlevels], center=TRUE, scale=FALSE)[,1]
  	  } else glmModel$coefficients[xlevels] <- 0
    }
  }
  return(glmModel)
}

fixMissingFactor <- function(factor_var, factor, glmModel) {
  # Following procedure adapted from proposal made by Sam Thomas:
  # https://stackoverflow.com/questions/31917878/handling-unknown-factor-levels-in-r-glm
  # coefficients from model
  coefficients <- coef(glmModel)[grepl(factor_var, names(coef(glmModel)))]
  # calculate coefficient for new factors
  coefficients_offset <- mean(coefficients, na.rm=TRUE)
  # insert missing coefficient
  coefficients[factor] <- coefficients_offset
  # center coefficients
  coefficients <- scale(coefficients, center=TRUE, scale=FALSE)[, 1]
  
  # replace coefficients from model
  modelcoef <- coef(glmModel)
  # add offset to intercept
  modelcoef[["(Intercept)"]] <- modelcoef[["(Intercept)"]] + coefficients_offset
  # all new coefficients
  modelcoef[names(coefficients)] <- coefficients
  glmModel$coefficients <- modelcoef
  return(glmModel)
}

predict_whileCheckingMissingFactors <- function (glmModel, newdata, type="response") {
  
  previous_na_action <- options('na.action')
  options(na.action='na.pass')
  newxreg <- model.matrix(glmModel$formula, data=newdata)
  options(na.action=previous_na_action$na.action)
  newxreg[is.na(newxreg)] <- 0
  
  predictions <- rep(NA, nrow(newdata))
  factor_vars <- names(glmModel$xlevels)
  missing_factors <- setdiff(colnames(newxreg), names(glmModel$coefficients))
  touchedFactors <- FALSE
  
  if (length(missing_factors) > 0) {
    touchedFactors <- TRUE
    factor_vars_indexes <- NULL
    for (i in 1:length(missing_factors)) {
      missing_factor <- missing_factors[i]
      tmp <- unlist(lapply(1:length(factor_vars), FUN = function(x) if (length(grep(factor_vars[x], missing_factor)) > 0) x ))
      factor_vars_indexes[[factor_vars[tmp]]] <- tmp
      glmModel <- fixMissingFactor(factor_var=factor_vars[tmp], factor=missing_factor, glmModel)
    }
    rm(tmp)
  }
  glmModel <- fixNAcoefficients(glmModel)
  glmModel$coefficients[is.nan(glmModel$coefficients)] <- 0
  # Beta0 + Beta1x...
  #glmModel$coefficients[["(Intercept)"]] + 
  pcoef <- newxreg[,] %*% glmModel$coefficients[colnames(newxreg[,])]
    #predicted response
  if (type == "response") {
    require(Brobdingnag) # using library for precision when large and small numbers go into formula
    predictions <- unlist(lapply(pcoef, function(x) {x <- as.brob(x); as.numeric(exp(x) / (1 + exp(x)))}))
  } else (print("Warning, only 'response' type prediction was implemented"))
  return(predictions)
}