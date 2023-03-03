
#' Process DLS and titer training data, train model and predict titer loss for new DLS data
#'
#' This is a wrapper function that runs all of the analysis, replicating results from the Phage-ELF shiny app. 
#' @param data A dataframe with DLS data. One row per sample
#' @param lytic A dataframe with the titer associated with each treatment.
#' @param training_pairs A dataframe with two columns, named 'treatment' and 'control'.
#' @param testing_pairs Optional. A dataframe with two columns, named 'treatment' and 'control'.
#' @param metric A character string specifying the DLS metric to use. Default is 'Intensity'. Other metrics typically include 'Volume' and 'Number'. Column names in 'data' must contain the substrings 'intens','volume' and 'number' respectively (excluding the 'sample' column).
#' @param intercept If TRUE, include a non-zero intercept in the linear model. Default is TRUE.
#' @export
#' @examples
#' run_all(dls_data, lytic_data, training_pairs, testing pairs)

run_all <- function(data, lytic, training_pairs, testing_pairs, metric = "Intensity", intercept = TRUE) {
  
  # Create table to store AUC difference
  diff_total <- tibble::tibble(treatment = character(), control = character(), AUC = numeric())
  
  # Call AUC_diff() function to calculate AUC diff for each treatment pair
  for (i in unique(training_pairs$control)){
    
    treatments_in <- training_pairs$treatment[training_pairs$control == i]
    
    diff <- AUC_diff(data, treatments_in, i, metric)
    
    diff_total <- rbind(diff_total, diff)
    
  }
  
  # Train the model
  lm_total <- lm_train(diff_total, lytic, intercept)
  
  training_data <- lm_total$model
  training_data$treatment <- diff_total$treatment
  
  if (missing(testing_pairs)) { 
    return(list(data = training_data, model = lm_total))
  } else {
    
    # Calculate AUC diff for each test pair
    diff_test <- tibble::tibble(treatment = character(), control = character(), AUC = numeric())
    
    for (i in unique(testing_pairs$control)){ # Debug from here
      
      treatments_in <- testing_pairs$treatment[testing_pairs$control == i]
      
      diff <- AUC_diff(data, treatments_in, i, metric)
      
      diff_test <- rbind(diff_test, diff)
      
    }
    
    pred <- lm_test(lm_total, diff_test)
    
    return(list(data = training_data, model = lm_total, predictions = pred))
    
  }

}

#' Calculate the difference in Area-Under-the-Curve for DLS data
#'
#' This function allows you to calculate the difference in AUC between treatment and control conditions. 
#' @param data A dataframe with DLS data. One row per sample
#' @param treatments A vector of strings corresponding to a subset of the samples in 'data'.
#' @param control A string corresponding to one sample in 'data'.
#' @param metric A character string specifying the DLS metric to use. Default is 'Intensity'. Other metrics typically include 'Volume' and 'Number'. Column names in 'data' must contain the substrings 'intens','volume' and 'number' respectively (excluding the 'sample' column).
#' @export
#' @examples
#' AUC_diff(dls_data, c('treatment_1', 'treatment_2'), 'control', 'Intensity')

AUC_diff <- function(data, treatments_in, control_in, metric = "Intensity") {
  
  names(data)[names(data) == 'Sample.Name'] <- 'sample'
  
  if (metric == 'Intensity') {
    match.metric <- 'intens'
  } else {match.metric <- metric}
 
  dls <- data[,grep(match.metric,colnames(data), ignore.case = TRUE)]
  
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
#' @param metric A character string specifying the DLS metric to use. Default is 'Intensity'. Other metrics typically include 'Volume' and 'Number'. Column names in 'data' must contain the substrings 'intens','volume' and 'number' respectively (excluding the 'sample' column).
#' @param size Optional. A vector of the particle sizes associated with the DLS measurements in 'data'.
#' @export
#' @examples
#' plot_DLS(dls_data, c('treatment_1', 'treatment_2'), 'control', 'Intensity', size_vector)

plot_DLS <- function(data, treatments, control, metric = 'Intensity', size){
  
    names(data)[names(data) == 'Sample.Name'] <- 'sample'
    
    if (metric == 'Intensity') {
      match.metric <- "intens"
    } else {match.metric <- metric}
    
    # Remove extra columns and add back columns with sample names
    dls <- data[,grep(match.metric,colnames(data), ignore.case = TRUE)]
    
    if (!missing(size)){
      dls <- dls[,1:length(size)]
    }

    dls$sample <- data$sample
    
    # Remove samples not in treatments or control
    dls <- dls[(dls$sample %in% treatments) | (dls$sample %in% control),]
    
    # Change to long format
    dls_long <- tidyr::gather(dls, key = "size", value = "value", grep(match.metric,colnames(dls), ignore.case = TRUE))
    
    # Order by sample
    dls_long <- dls_long[order(dls_long$sample),]
    
    # Add size if provided or linear vector if not
    if (missing(size)) {
      dls_long$xscale <- rep(10^(1:length(grep(match.metric,colnames(data), ignore.case = TRUE))/ncol(dls)), nrow(dls))
    }
    else { dls_long$xscale <- rep(size, nrow(dls))}

    # Plot one curve for each sample
    ggplot2::ggplot(dls_long, ggplot2::aes(x = xscale, y = value, color = sample)) + 
      ggplot2::geom_line() +
      ggplot2::labs(x = bquote("D"[h]~"(nm)"), y = paste(metric, '(%)'), color = 'Treatment') +
      ggplot2::scale_x_continuous(trans='log10') +
      ggplot2::theme_classic() + 
      ggplot2::theme(axis.text = ggplot2::element_text(size = 13), axis.title = ggplot2::element_text(size = 14), 
                     legend.text = ggplot2::element_text(size = 14), legend.title = ggplot2::element_blank()) +
      ggplot2::annotation_logticks(sides = 'b')
  
}

#' Train linear model for titer loss over difference in AUC.
#'
#' This function allows you to train a linear model of the titer loss associated with the difference in AUC calculated from DLS data. 
#' @param AUC A dataframe with the difference in AUC for each treatment. 
#' @param titer A dataframe with the titer associated with each treatment.
#' @param intercept If TRUE, include a non-zero intercept in the linear model. Default is TRUE.
#' @export
#' @examples
#' lm_train(AUC_data, titer_data)

lm_train <- function(AUC, titer, intercept = TRUE) {
  
  AUC$titer_loss <- NA
  
  for (n in 1:nrow(AUC)){
    
    control <- titer$lytic[titer$treatment == AUC$control[n]]
    treatment <- titer$lytic[titer$treatment == AUC$treatment[n]]
    AUC$titer_loss[n] <- log10(treatment+1) - log10(control+1)
    
  }
  
  # Train linear model with or without intercept
  if (intercept == FALSE) {
    model <- lm(titer_loss ~ AUC - 1, data = AUC)
  } else {
    model <- lm(titer_loss ~ AUC, data = AUC)
  }
  
  return(model)
  
}

#' Plot linear model for titer loss over difference in AUC.
#'
#' This function allows you to plot the training data and the trained linear model of the titer loss 
#' associated with the difference in AUC calculated from DLS data. 
#' @param model A trained linear model. Output of lm_train.
#' @export
#' @examples
#' plot_model(AUC_data, titer_data, 'control')

plot_model <- function(model) {

  AUC_titer <- model$model
  
  newdata <- data.frame(AUC = seq(0, 200,len=500))
  
  #use fitted model to predict values of vs
  newdata$pred <- predict(model, newdata, interval = 'confidence')[,'fit']
  newdata$upper <- predict(model, newdata, interval = 'confidence')[,'upr']
  newdata$lower <- predict(model, newdata, interval = 'confidence')[,'lwr']
  
  ggplot2::ggplot()  +
    ggplot2::geom_point(data = AUC_titer, ggplot2::aes(x = AUC, y = titer_loss)) +
    ggplot2::geom_line(data = newdata, ggplot2::aes(x = AUC, y = pred)) +
    ggplot2::geom_ribbon(data = newdata, ggplot2::aes(x = AUC,ymin = lower, ymax = upper), alpha=0.3) +
    ggplot2::labs(x = 'AUC Difference (%)', y = 'Difference in lytic activity (log10)') +
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text = ggplot2::element_text(size = 13), axis.title = ggplot2::element_text(size = 14))
  
  
}

#' Use linear model to predict titer loss based on a difference in AUC 
#'
#' This function allows you to predict titer loss based on a trained linear model and on a difference in AUC, which can be obtained with the function 'AUC_Diff'.  
#' @param AUC A dataframe with the difference in AUC for each treatment. Output of AUC_diff().
#' @param model A trained linear model. Output of lm_train().
#' @export
#' @examples
#' lm_test(lm_model, AUC_data)

lm_test <- function(model, AUC) {
  
  predicted_titer <- data.frame(treatment = AUC$treatment, AUC = AUC$AUC, titer_loss = as.numeric(predict(model, AUC)))
  
  return(predicted_titer)
  
}

#' Plot titer loss predictions
#'
#' This function plots titer loss predictions with the training dat and model  
#' @param AUC A dataframe with the difference in AUC for each treatment. Output of AUC_diff().
#' @param model A trained linear model. Output of lm_train().
#' @export
#' @examples
#' plot_test_data(lm_model, AUC_data)

plot_test_data <- function(model, AUC) {
  
  predicted_titer <- data.frame(AUC = AUC$AUC, titer_loss = as.numeric(predict(model, AUC)))
  
  AUC_titer <- model$model
  
  newdata <- data.frame(AUC = seq(0, 200,len=500))
  
  #use fitted model to predict values of vs
  newdata$pred <- predict(model, newdata, interval = 'confidence')[,'fit']
  newdata$upper <- predict(model, newdata, interval = 'confidence')[,'upr']
  newdata$lower <- predict(model, newdata, interval = 'confidence')[,'lwr']  
  
  ggplot2::ggplot()  +
    ggplot2::geom_point(data = AUC_titer, ggplot2::aes(x = AUC, y = titer_loss, color = 'Training data')) +
    ggplot2::geom_line(data = newdata, ggplot2::aes(x = AUC, y = pred)) +
    ggplot2::geom_ribbon(data = newdata, ggplot2::aes(x = AUC,ymin = lower, ymax = upper), alpha=0.3) +
    ggplot2::geom_point(data = predicted_titer, ggplot2::aes(x = AUC, y = titer_loss, color = 'Predictions')) +
    ggplot2::labs(x = 'AUC Difference (%)', y = 'Difference in lytic activity (log10)') +
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text = ggplot2::element_text(size = 13), axis.title = ggplot2::element_text(size = 14), 
                   legend.text = ggplot2::element_text(size = 14), legend.title = ggplot2::element_blank())
  
}


