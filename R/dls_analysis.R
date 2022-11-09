
#' Calculate the difference in Area-Under-the-Curve for DLS data
#'
#' This function allows you to calculate the difference in AUC between treatment and control conditions. 
#' @param data A dataframe with DLS data. One row per sample
#' @param treatments A vector of strings corresponding to a subset of the samples in 'data'.
#' @param control A string corresponding to one sample in 'data'.
#' @param metric A string specifying the DLS metric used. Must be a substring of column names in 'data', excluding the 'sample' column.
#' @export
#' @examples
#' AUC_diff(dls_data, c('treatment_1', 'treatment_2'), 'control', 'intensities')

AUC_diff <- function(data, treatments_in, control_in, metric) {
  
  dls <- data[,grep(metric,colnames(data), ignore.case = TRUE)]
  
  diff <- tibble::tibble(treatment = treatments_in, control = rep(control_in, length(treatments_in)), AUC = NaN)
  
  if (!(control_in %in% data$sample)) {
    
    stop('Control not found')
    
  }
  
  for (i in treatments_in) {
    
    if (i %in% data$sample) {
      control <- dls[data$sample == control_in,]
      treatment <- dls[data$sample == i,]
      diff$AUC[diff$treatment == i] <- sum(abs(control-treatment))
    }
    else {print(paste('Treatment ',i ,' not found'))}
    
  }
  
  return(diff)
  
}


#' Plot DLS curves
#'
#' This function allows you to plot different DLS curves on the same graph for comparison.
#' @param data A dataframe with DLS data. One row per sample
#' @param treatments A vector of strings corresponding to a subset of the samples in 'data'.
#' @param control A string corresponding to one sample in 'data'.
#' @param metric A string specifying the DLS metric used. Must be a substring of column names in 'data', excluding the 'sample' column.
#' @param size Optional. A vector of the particle sizes associated with the DLS measurements in 'data'.
#' @export
#' @examples
#' plot_DLS(dls_data, c('treatment_1', 'treatment_2'), 'control', 'intensities', size_vector)


plot_DLS <- function(data, treatments, control, metric, size){
    
    # Remove extra columns and add back columns with sample names
    dls <- data[,grep(metric,colnames(data), ignore.case = TRUE)]
    
    if (!missing(size)){
      dls <- dls[,1:length(size)]
    }

    dls$sample <- data$sample
    
    # Remove samples not in treatments or control
    dls <- dls[(dls$sample %in% treatments) | (dls$sample %in% control),]
    
    # Change to long format
    dls_long <- tidyr::gather(dls, key = "size", value = "value", grep(metric,colnames(dls), ignore.case = TRUE))
    
    # Order by sample
    dls_long <- dls_long[order(dls_long$sample),]
    
    # Add size if provided or linear vector if not
    if (missing(size)) {
      dls_long$xscale <- rep(10^(1:length(grep(metric,colnames(data), ignore.case = TRUE))/ncol(dls)), nrow(dls))
    }
    else { dls_long$xscale <- rep(size, nrow(dls))}

    # Plot one curve for each sample
    ggplot2::ggplot(dls_long, aes(x = xscale, y = value, color = sample)) + 
      geom_line() +
      labs(x = 'Size', y = metric, color = 'Treatment') +
      scale_x_continuous(trans='log10') +
      theme_bw()
  
}


#' Train linear model for titer loss over difference in AUC.
#'
#' This function allows you to train a linear model of the titer loss associated with the difference in AUC calculated from DLS data. 
#' @param AUC A dataframe with the difference in AUC for each treatment. 
#' @param titer A dataframe with the titer associated with each treatment.
#' @param control A string specifying the sample to be used as control.
#' @param intercept If TRUE, include a non-zero intercept in the linear model. Default is FALSE.
#' @export
#' @examples
#' lm_train(AUC_data, titer_data, 'control')


lm_train <- function(AUC, titer, control, intercept = FALSE) {
  
  # Calculate titer loss
  titer$titer_loss <- log10(titer$lytic + 1) - log10(titer$lytic[titer$treatment == control] + 1)
  
  # Merge AUC difference and titer loss data
  AUC_titer <- merge(AUC, titer, by = 'treatment')
  
  # Train linear model with or without intercept
  if (intercept == FALSE) {
    model <- lm(titer_loss ~ AUC - 1, data = AUC_titer)
  } else {
    model <- lm(titer_loss ~ AUC, data = AUC_titer)
  }
  
  newdata <- data.frame(AUC = seq(0, 200,len=500))
    
  #use fitted model to predict values of vs
  newdata$pred <- predict(model, newdata, interval = 'confidence')[,'fit']
  newdata$upper <- predict(model, newdata, interval = 'confidence')[,'upr']
  newdata$lower <- predict(model, newdata, interval = 'confidence')[,'lwr']
    
  plot_model <- ggplot2::ggplot()  +
    geom_point(data = AUC_titer, aes(x = AUC, y = titer_loss)) +
    geom_line(data = newdata, aes(x = AUC, y = pred)) +
    geom_ribbon(data = newdata, aes(x = AUC,ymin = lower, ymax = upper), alpha=0.3) +
    labs(x = 'AUC Difference', y = 'Difference in lytic activity') +
    theme_bw()
  
  return(list(model = model, plot = plot_model))
  
}

#' Use linear model to predict titer loss based on a difference in AUC 
#'
#' This function allows you to predict titer loss based on a trained linear model and on a difference in AUC, which can be obtained with the function 'AUC_Diff'.  
#' @param AUC A dataframe with the difference in AUC for each treatment. 
#' @param model A linear model. 
#' @export
#' @examples
#' lm_test(AUC_data, lm_model)

lm_test <- function(model, AUC) {
  
  predicted_titer <- data.frame(AUC = AUC$AUC, titer_loss = as.numeric(predict(model, AUC)))
  
  AUC_titer <- model$model
  
  newdata <- data.frame(AUC = seq(0, 200,len=500))
  
  #use fitted model to predict values of vs
  newdata$pred <- predict(model, newdata, interval = 'confidence')[,'fit']
  newdata$upper <- predict(model, newdata, interval = 'confidence')[,'upr']
  newdata$lower <- predict(model, newdata, interval = 'confidence')[,'lwr']  
  
  plot_all <- ggplot2::ggplot()  +
    geom_point(data = AUC_titer, aes(x = AUC, y = titer_loss)) +
    geom_line(data = newdata, aes(x = AUC, y = pred)) +
    geom_ribbon(data = newdata, aes(x = AUC,ymin = lower, ymax = upper), alpha=0.3) +
    geom_point(data = predicted_titer, aes(x = AUC, y = titer_loss), color = 'red') +
    labs(x = 'AUC Difference', y = 'Difference in lytic activity (log10)') +
    theme_bw()
  
  return(list(pred = predicted_titer, plot= plot_all))
  
}


