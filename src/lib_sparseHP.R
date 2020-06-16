library(ggplot2)
library(quantreg)
library(pracma)
library(CVXR)
library(gurobi)
library(coronavirus)
library(moments)
library(dplyr)
library(tidyr)
library(tidyverse)

########################################################
#
#
# Collection of Filters
#
#
########################################################
# R function for L2 trend filtering
# library CVXR is needed 
# Inputs: 
#         y: n*1 vector of time series
#         lambda: penalization parameter
# Outputs:       
#         betahat: n*1 vector of trend estimates
# [Note. There is a closed-form solution for HP filter. However, the CVXR format is used here]

l2tf = function(y=y,lambda=lambda){
  
  beta = Variable(length(y))
  objective = Minimize(p_norm(y - beta)^2 +
                         lambda * p_norm(diff(x = beta, differences = 2))^2)
  prb = Problem(objective)
  betahat = solve(prb)$getValue(beta)   
  fidelity = sum((y - betahat)^2)
  
  return(list(betahat=betahat, fidelity=fidelity))
  
  # Clear space
  rm(beta, objective, prb, betahat)
}



# R function for L1 trend filtering
# library CVXR is needed
# Inputs: 
#         y: n*1 vector of time series
#         lambda: penalization parameter
# Outputs:       
#         betahat: n*1 vector of trend estimates

l1tf = function(y=y,lambda=lambda){
  
  beta = Variable(length(y))
  objective = Minimize(p_norm(y - beta)^2 +
                         lambda * p_norm(diff(x = beta, differences = 2), 1))
  prb = Problem(objective)
  betahat = solve(prb)$getValue(beta)   
  fidelity = sum((y - betahat)^2)
  
  return(list(betahat=betahat, fidelity=fidelity))
  
  
  # Clear space
  rm(beta, objective, prb, betahat)
}


# R function for square-root L1 trend filtering
# library CVXR is needed
# Inputs: 
#         y: n*1 vector of time series
#         lambda: penalization parameter
# Outputs:       
#         betahat: n*1 vector of trend estimates

sqrt_l1tf = function(y=y,lambda=lambda){
  
  beta = Variable(length(y))
  objective = Minimize(p_norm(y - beta) +
                         lambda * p_norm(diff(x = beta, differences = 2), 1))
  prb = Problem(objective)
  betahat = solve(prb)$getValue(beta)   
  fidelity = sum((y - betahat)^2)
  
  return(list(betahat=betahat, fidelity=fidelity))
  
  
  # Clear space
  rm(beta, objective, prb, betahat)
}


# R function to solve L0-constrained trend filtering
# library gurobi is needed
# Inputs: 
# 
#        
#  
#        
# Outputs:       
#        

l0tfc <- function(y=y,l0constraint=l0constraint,M=M,l2penalty=l2penalty,...){
  
  n <- length(y)
  D_matrix <- matrix(0,nrow=(n-2),ncol=n) 
  
  for (i in 1:(n-2))
  {
    D_matrix[i,i:(i+2)]=c(1,-2,1)        
  }                    
  
  A1 <- cbind(D_matrix,-M*diag(n-2))
  A2 <- cbind(-D_matrix,-M*diag(n-2))
  A3 <- matrix(c(rep(0,n),rep(1,(n-2))),nrow=1)
  A_matrix <- rbind(A1,A2,A3)
  A_c <- ncol(A_matrix)
  A_r <- nrow(A_matrix)
  
  b_vector <- c(rep(0,(A_r-1)),l0constraint)
  
  c_vector <- c(-2*y,rep(0,(n-2)))
  
  Q_matrix <- matrix(0,nrow=A_c,ncol=A_c)
  Q_matrix[1:n,1:n] <- diag(n) + l2penalty*(t(D_matrix)%*%D_matrix)
  
  x_ub <- c(rep(max(y),n),rep(1,(n-2)))
  x_lb <- c(rep(min(y),n),rep(0,(n-2)))
  
  x_type <- rep('C',n)
  z_type <- rep('B',(n-2))
  
  model <- list()
  
  model$A          <- A_matrix
  model$rhs        <-  b_vector  
  #  model$sense      <- '<'
  model$sense      <- c(rep('<',A_r-1),"<=")
  model$obj        <- c_vector
  model$Q          <- Q_matrix 
  model$modelsense <- 'min'
  model$lb         <- x_lb
  model$ub         <- x_ub
  model$vtype      <- c(x_type,z_type)
  
  params <- list(OutputFlag=0)
  
  result <- gurobi(model, params)
  
  return(result)
  
  # Clear space
  rm(model, result, params)
}


########################################################
#
#
# Collection of Tuning Parameter Selections
#
#
########################################################




########################################################
### Select the tuning parameter for HP filter (l2tf) ###
########################################################
l2tf_tuning <-function(y, fidelity, ld_set, country.name, file_name){
  n_ld_set = length(ld_set)
  result = data.frame(ld_set, obj_val = rep(NA, n_ld_set), abs_dist = rep(NA, n_ld_set))
  result_all <- foreach(i = 1:n_ld_set, .combine = rbind ) %dopar% l2tf(y=y,lambda=ld_set[i])
  result[,2] = unlist(result_all[,2])
  result[,3] = abs(result[,2] - fidelity)
  min_index = which.min(result[,3])
  print(result[min_index,])
  
  data = data.frame(lambda = result[,1], abs_dist = result[,3])
  p = ggplot(data, aes(x=lambda, y=abs_dist)) +
    geom_line() +
    ggtitle(paste0(country.name, " : Tuning Parameter (HP)")) +
    ylab("Absolute Distance of Fidelities") +
    xlab(bquote(lambda)) +
    theme(axis.title = element_text(size = 15)) +
    geom_vline(xintercept = result[min_index,1], linetype="dashed", color="red")
  print(p)
  gg_file_name = paste("../results/tuning_",file_name,".pdf",sep="")
  ggsave(gg_file_name, width=4, height=4, units="in")

  return(list(opt_ld=result[min_index,1] , result=result))
}

########################################################
### Select the tuning parameter for L1 filter (l1tf) ###
########################################################
l1tf_tuning <-function(y, fidelity, ld_set, country.name, file_name){

  n_ld_set = length(ld_set)
  result = data.frame(ld_set, obj_val = rep(NA, n_ld_set), abs_dist = rep(NA, n_ld_set))
  result_all <- foreach(i = 1:n_ld_set, .combine = rbind ) %dopar% l1tf(y=y,lambda=ld_set[i])
  result[,2] = unlist(result_all[,2])
  result[,3] = abs(result[,2] - fidelity)
  min_index = which.min(result[,3])
  print(result[min_index,])
  
  data = data.frame(lambda = result[,1], abs_dist = result[,3])
  p = ggplot(data, aes(x=lambda, y=abs_dist)) +
    geom_line() +
    ggtitle(paste0(country.name, ": Tuning Parameter (L1)")) +
    ylab("Absolute Distance of Fidelities") +
    xlab(bquote(lambda)) +
    theme(axis.title = element_text(size = 15)) +
    geom_vline(xintercept = result[min_index,1], linetype="dashed", color="red")
  print(p)
  gg_file_name = paste("../results/tuning_",file_name,".pdf",sep="")
  ggsave(gg_file_name, width=4, height=4, units="in")

  return(list(opt_ld=result[min_index,1] , result=result))
}



##################################################################
### Select the tuning parameter for SQRT_L1 filter (sqrt_l1tf) ###
##################################################################
sqrt_l1tf_tuning <-function(y, fidelity, ld_set, country.name, file_name){
  
  n_ld_set = length(ld_set)
  result = data.frame(ld_set, obj_val = rep(NA, n_ld_set), abs_dist = rep(NA, n_ld_set))
  result_all <- foreach(i = 1:n_ld_set, .combine = rbind ) %dopar% sqrt_l1tf(y=y,lambda=ld_set[i])
  result[,2] = unlist(result_all[,2])
  result[,3] = abs(result[,2] - fidelity)
  min_index = which.min(result[,3])
  print(result[min_index,])
  
  data = data.frame(lambda = result[,1], abs_dist = result[,3])
  p = ggplot(data, aes(x=lambda, y=abs_dist)) +
    geom_line() +
    ggtitle(paste0(country.name, ": Tuning Parameter (SQRT_L1)")) +
    ylab("Absolute Distance of Fidelities") +
    xlab(bquote(lambda)) +
    theme(axis.title = element_text(size = 15)) +
    geom_vline(xintercept = result[min_index,1], linetype="dashed", color="red")
  print(p)
  gg_file_name = paste("../results/tuning_",file_name,".pdf",sep="")
  ggsave(gg_file_name, width=4, height=4, units="in")
  
  return(list(opt_ld=result[min_index,1] , result=result))
}


# R functions for Leave-one-out Cross-validations for L0 constraint and L2 penalty
#
#
# l0tfc_cv        : Calculate all the paths of the obj values for the set of tuning parameters (k,ld)
# l0tfc_single    : Return the objective function values using the output of l0tfc_leave
# l0tfc_leave     : Estimat betaHat for a single pair of tuning parameters (k,ld) when one observatino is dropped (leave)
#

l0tfc_cv <- function(y, l0constraint_set, M, l2penalty_set){
  n <- length(y)
  n_l0 <- length(l0constraint_set)
  n_l2 <- length(l2penalty_set)
  n_grid <- n_l0*n_l2
  
  result_set <- expand.grid(l2penalty_set, l0constraint_set)
  result_set <- cbind(result_set, rep(NA, n_grid))
  colnames(result_set) <- c('l2', 'l0', 'obj_val')
  
  result_set[,3] <- foreach(i = 1:n_grid, .combine = c ) %dopar% l0tfc_cv_single(y=y, l0constraint=result_set[i,2], M=M, l2penalty=result_set[i,1])
  #result_set[,3] <- foreach(i = 1:n_grid, .combine = c ) %do% l0tfc_cv_single(y=y, l0constraint=result_set[i,2], M=M, l2penalty=result_set[i,1])
  
  return(result_set)
}

l0tfc_cv_single <- function(y, l0constraint, M, l2penalty){
  
  n <- length(y)
  betaHat_s = rep(NA, n)                    # Leave-one-out estimates for betaHat

  for (leave.s in c(1:n) ){
    betaHat_s[leave.s] = l0tfc_leave(y=y, l0constraint=l0constraint,M=M,l2penalty=l2penalty, leave=leave.s)   
  }
  
  obj_val = t(betaHat_s - y) %*% (betaHat_s - y)
  
  return(obj_val)     
}   


l0tfc_leave <- function(y, l0constraint, M=M, l2penalty, leave){
  
  n <- length(y)
  D_matrix <- matrix(0,nrow=(n-2),ncol=n) 
  
  for (i in 1:(n-2))
  {
    D_matrix[i,i:(i+2)]=c(1,-2,1)        
  }                    
  
  A1 <- cbind(D_matrix,-M*diag(n-2))
  A2 <- cbind(-D_matrix,-M*diag(n-2))
  A3 <- matrix(c(rep(0,n),rep(1,(n-2))),nrow=1)
  A_matrix <- rbind(A1,A2,A3)
  A_c <- ncol(A_matrix)
  A_r <- nrow(A_matrix)
  
  b_vector <- c(rep(0,(A_r-1)),l0constraint)
  
  c_vector <- c(-2*y,rep(0,(n-2)))
  c_vector_leave <- c_vector
  c_vector_leave[leave] <- 0 
  
  Q_matrix <- matrix(0,nrow=A_c,ncol=A_c)
  Q_matrix_leave = Q_matrix
  Q_matrix_leave[1:n,1:n] = diag(n)
  Q_matrix_leave[leave,leave] = 0
  Q_matrix_leave[1:n,1:n] = Q_matrix_leave[1:n,1:n] + l2penalty*(t(D_matrix)%*%D_matrix)
  
  x_ub <- c(rep(max(y),n),rep(1,(n-2)))
  x_lb <- c(rep(min(y),n),rep(0,(n-2)))
  
  x_type <- rep('C',n)
  z_type <- rep('B',(n-2))
  
  model <- list()
  
  model$A          <- A_matrix
  model$rhs        <-  b_vector  
  model$sense      <- c(rep('<',A_r-1),'<=')
  model$obj        <- c_vector_leave
  model$Q          <- Q_matrix_leave 
  model$modelsense <- 'min'
  model$lb         <- x_lb
  model$ub         <- x_ub
  model$vtype      <- c(x_type,z_type)
  
  params <- list(OutputFlag=0, Threads=6, TimeLimit=300)
  
  result <- gurobi(model, params)
  
  return(result$x[leave])
  #return(result)
  
  # Clear space
  rm(model, result, params)
}

########################################################################
#
#
# COVID-19 Data extraction
#
#
########################################################################
data_extract <- function(country.name, final.date, save.csv=FALSE){
  update_dataset()
  
  data("coronavirus")
  output_file_name = paste0("covid_",country.name,"_",max(coronavirus$date),".csv")
  #########################################################
  # Change the country.name to match those in the dataset
  #########################################################
  if (country.name =="South Korea") {
    country.name = "Korea, South"  
  } else if (country.name == "UK") {
    country.name = "United Kingdom"
  }
  
  data_out <- coronavirus %>%
    filter(country==country.name) %>%
    filter(date <= final.date) %>%
    select(date, type, cases) %>%
    group_by(date, type) %>%
    summarize(total_cases = sum(cases)) %>%
    pivot_wider(names_from = type,
                values_from = total_cases)
  
  
  if (save.csv){
    write.csv(data_out, file=output_file_name)
  } 
  
  return(data_out)  
  
}

######################################################################
####  R function to generate variables                            ####
####  Inputs: population size                                     #### 
####          daily confirmed cases                               #### 
####          daily recovered cases                               #### 
####          daily deaths                                        ####  
####          date  ["Date" class]                                ####  
####          minimum number of cumulative cases                  ####
####  Output: date and key variables [Variable for analysis]      ####  
######################################################################
gen_vars = function(date=date,
                    confirmed=confirmed,
                    death=death,
                    recovered=recovered,
                    population=population,
                    min_cum_case=min_cum_case){
  
  # setting up and normalizing by population size  
  Np = population 
  cum_case = cumsum(confirmed)/Np
  cum_death = cumsum(death)/Np
  cum_recov = cumsum(recovered)/Np
  cum_susc = 1 - cum_case
  cum_infec = cum_case - cum_recov - cum_death
  
  # sample size    
  n = length(cum_case)  
  seq_tmp = 1:n
  
  # first period          
  
  date_tmp = date[cumsum(confirmed) > min_cum_case]     
  first_period = seq_tmp[date==date_tmp[1]]
  
  # date          
  day =  date[(first_period+1):n]  
  
  # contructing lags          
  cum_infec_lag = cum_infec[first_period:(n-1)]
  cum_case_lag = cum_case[first_period:(n-1)]
  cum_susc_lag = cum_susc[first_period:(n-1)]   
  cum_recov_lag = cum_recov[first_period:(n-1)]   
  
  # contructing changes          
  change_recov = cum_recov[(first_period+1):n] - cum_recov_lag
  change_infec = cum_infec[(first_period+1):n] - cum_infec_lag
  change_case = cum_case[(first_period+1):n] - cum_case_lag
  
  # contructing growth rates        
  growth_case = change_case/cum_case_lag
  growth_infec = change_infec/cum_infec_lag
  
  # contructing main variables 
  
  DeltaC = change_case
  I_lag = cum_infec_lag
  S_lag = cum_susc_lag
  Y = DeltaC/(I_lag*S_lag)
  
  vars = data.frame(
    date = day,
    DeltaC = DeltaC,
    I_lag = I_lag,      
    S_lag = S_lag,
    Y = Y
  )      
  
  return(vars)          
  
}  


########################################################################
#
#
# Functions to draw graphs
#
#
########################################################################
draw_fit_R0_SHP_L1 <- function(data, SHP_kink_dates, L1_kink_dates, country.name, file_name){
  p = ggplot(data = data) +
    geom_line(mapping = aes(x = date, y = R), color = 'grey50') +
    labs(x = "Date", y = bquote(R[0])) +
    geom_line(mapping = aes(x = date, y = filter_Rhat), color = 'blue', size = 0.7) + 
    geom_line(mapping = aes(x = date, y = L1_filter_Rhat), linetype="twodash", color = 'red', size = 0.7) +
    geom_vline(xintercept = SHP_kink_dates, col="blue") +
    geom_vline(xintercept = L1_kink_dates, linetype = "twodash", col="red") +
    ggtitle(bquote(.(country.name) ~ ": " ~ R[0] ~ "Fit")) +
    theme_minimal()
  
  if (lockdown_date[1] != 0){
    p = p + geom_vline(xintercept = date[lockdown_date], col="orange") 
  }
  
  gg_file_name = paste("../results/R0_",file_name,".pdf",sep="")
  ggsave(gg_file_name, width=6, height=4, units="in")
}

draw_resid_R0_SHP_L1 <- function(data, SHP_kink_dates, L1_kink_dates, country.name, file_name){
  p = ggplot(data = data) +
    geom_line(mapping = aes(x = date, y = 0), color = 'grey50') +
    labs(x = "Date", y = bquote(R[0])) +
    geom_line(mapping = aes(x = date, y = filter_Rresid), color = 'blue', size = 0.7) + 
    geom_line(mapping = aes(x = date, y = L1_filter_Rresid), linetype="twodash", color = 'red', size = 0.7) +
    geom_vline(xintercept = SHP_kink_dates, col="blue") +
    geom_vline(xintercept = L1_kink_dates, linetype = "twodash", col="red") +
    ggtitle(bquote(.(country.name) ~ ": Residuals from "~ R[0] ~"Fit")) +
    theme_minimal()
  
  if (data$lockdown_date[1] != 0){
    p = p + geom_vline(xintercept = date[lockdown_date], col="orange") 
  }
  
  gg_file_name = paste("../results/R0_resid_",file_name,".pdf",sep="")
  ggsave(gg_file_name, width=6, height=4, units="in")
}

draw_fit_beta_SHP_L1 <- function(data, SHP_kink_dates, L1_kink_dates, country.name, file_name){
  p = ggplot(data = data) +
    geom_line(mapping = aes(x = date, y = ly), color = 'grey50') +
    labs(x = "Date", y = bquote(log~beta[t])) +
    geom_line(mapping = aes(x = date, y = filter_yhat), color = 'blue', size = 0.7) + 
    geom_line(mapping = aes(x = date, y = L1_filter_yhat), linetype = "twodash", color = 'red', size = 0.7) +
    geom_vline(xintercept = SHP_kink_dates, col="blue") +
    geom_vline(xintercept = L1_kink_dates, linetype = "twodash", col="red") +
    ggtitle(bquote(.(country.name) ~ ": " ~ log ~ beta[t] ~ "Fit")) +
    theme_minimal()
  
  if (data$lockdown_date[1] != 0){
    p = p + geom_vline(xintercept = date[lockdown_date], col="orange") 
  }
  
  gg_file_name = paste("../results/log_beta_",file_name,".pdf",sep="")
  ggsave(gg_file_name, width=6, height=4, units="in")
}

draw_resid_beta_SHP_L1 <- function(data, SHP_kink_dates, L1_kink_dates, country.name, file_name){
  p = ggplot(data = data) +
    geom_line(mapping = aes(x = date, y = 0), color = 'grey50') +
    labs(x = "Date", y = bquote(log~beta[t])) +
    geom_line(mapping = aes(x = date, y = filter_resid), color = 'blue', size = 0.7) + 
    geom_line(mapping = aes(x = date, y = L1_filter_resid), linetype = "twodash", color = 'red', size = 0.7) +
    geom_vline(xintercept = SHP_kink_dates, col="blue") +
    geom_vline(xintercept = L1_kink_dates, linetype = "twodash", col="red") +
    ggtitle(bquote(.(country.name) ~ ": Residuals from "~ log ~ beta[t] ~ "Fit")) +
    theme_minimal()
  
  if (data$lockdown_date[1] != 0){
    p = p + geom_vline(xintercept = date[lockdown_date], col="orange") 
  }
  
  gg_file_name = paste("../results/log_beta_resid_",file_name,".pdf",sep="")
  ggsave(gg_file_name, width=6, height=4, units="in")
}


#########################################################################################
#
# Single Model Draws
#
#########################################################################################


draw_fit_R0 <- function(data, kink_dates=NULL, country.name, file_name){
  p = ggplot(data = data) +
    geom_line(mapping = aes(x = date, y = R), color = 'grey50') +
    labs(x = "Date",  y = bquote(R[0])) +
    geom_line(mapping = aes(x = date, y = lRhat), color = 'red', size = 0.7) +
    geom_line(mapping = aes(x = date, y = filter_Rhat), color = 'blue', size = 0.7) + 
    geom_vline(xintercept = date[lockdown_date], col="orange") +
    ggtitle(bquote(.(country.name) ~ ": " ~ R[0] ~ "Fit")) +
    theme_minimal()
  if ( !is.null(kink_dates)) {
    p = p + geom_vline(xintercept = kink_dates, col="blue", size = 0.35) 
  }
  gg_file_name = paste("../results/R0_",file_name,".pdf",sep="")
  ggsave(gg_file_name, width=6, height=4, units="in")
}

draw_resid_R0 <- function(data, country.name, file_name){
  ggplot(data = data) +
    geom_line(mapping = aes(x = date, y = 0), color = 'grey50') +
    labs(x = "Date", y = bquote(R[0])) +
    geom_line(mapping = aes(x = date, y = lRresid), color = 'red', size = 0.7) +
    geom_line(mapping = aes(x = date, y = filter_Rresid), color = 'blue', size = 0.7) + 
    geom_vline(xintercept = date[lockdown_date], col="orange") +
    ggtitle(bquote(.(country.name) ~ ": Residuals from "~ R[0] ~"Fit")) +
    theme_minimal()
  
  gg_file_name = paste("../results/R0_resid_",file_name,".pdf",sep="")
  ggsave(gg_file_name, width=6, height=4, units="in")
}

draw_fit_beta <- function(data, kink_dates=NULL, country.name, file_name){
  p = ggplot(data = data) +
    geom_line(mapping = aes(x = date, y = ly), color = 'grey50') +
    labs(x = "Date", y = bquote(log~beta[t])) +
    geom_line(mapping = aes(x = date, y = lyhat), color = 'red', size = 0.7) +
    geom_line(mapping = aes(x = date, y = filter_yhat), color = 'blue', size = 0.7) + 
    geom_vline(xintercept = date[lockdown_date], col="orange") +
    ggtitle(bquote(.(country.name) ~ ": " ~ log ~ beta[t] ~ "Fit")) +
    theme_minimal()
  if ( !is.null(kink_dates)) {
    p = p + geom_vline(xintercept = kink_dates, col="blue", size = 0.35) 
  }
  
  gg_file_name = paste("../results/log_beta_",file_name,".pdf",sep="")
  ggsave(gg_file_name, width=6, height=4, units="in")
}

draw_resid_beta <- function(data, country.name, file_name){
  ggplot(data = data) +
    geom_line(mapping = aes(x = date, y = 0), color = 'grey50') +
    labs(x = "Date", y = bquote(log~beta[t])) +
    geom_line(mapping = aes(x = date, y = lresid), color = 'red', size = 0.7) +
    geom_line(mapping = aes(x = date, y = filter_resid), color = 'blue', size = 0.7) + 
    geom_vline(xintercept = date[lockdown_date], col="orange") +
    ggtitle(bquote(.(country.name) ~ ": Residuals from "~ log ~ beta[t] ~ "Fit")) +
    theme_minimal()
  
  gg_file_name = paste("../results/log_beta_resid_",file_name,".pdf",sep="")
  ggsave(gg_file_name, width=6, height=4, units="in")
}

draw_legend <- function(data, method_name1, method_name2, file_name){

  df2 <- data %>%
    select(date, lresid, filter_resid) %>%
    gather(key = "variable", value = "value", -date)
  
  ggplot(df2, aes(x = date, y = value)) + 
    geom_line(aes(color = variable)) + 
    theme(legend.position="bottom") +
    #scale_colour_discrete(labels = c(method_name1, method_name2), values = c("blue", "red")) +
    scale_color_manual(values = c("blue", "red"), labels = c(method_name1, method_name2)) +
    theme(legend.title=element_blank()) 
  
  gg_file_name = paste("../results/legend_",method_name1,"_",method_name2,".pdf",sep="")
  ggsave(gg_file_name, width=6, height=4, units="in")
}

draw_legend2 <- function(data, method_name1, method_name2, file_name){
  
  
  df2 <- data %>%
    select(date, lresid, filter_resid) %>%
    gather(key = "variable", value = "value", -date)
  
  ggplot(df2, aes(x = date, y = value)) + 
    geom_line(aes(color = variable, linetype= variable)) + 
    theme(legend.position="bottom") +
    #scale_colour_discrete(labels = c(method_name1, method_name2), values = c("blue", "red")) +
    scale_colour_manual(values = c("blue", "red"), labels = c(method_name1, method_name2)) +
    theme(legend.title=element_blank()) 
  
  gg_file_name = paste("../results/legend_",method_name1,"_",method_name2,".pdf",sep="")
  ggsave(gg_file_name, width=6, height=4, units="in")
}




