rm(list=ls())   # initialization

library(ggplot2)
library(quantreg)
library(pracma)
library(CVXR)
library(gurobi)
library(moments)
library(viridis)
library(doMC)
library(parallel)

n_cores = detectCores()
registerDoMC(cores=n_cores)
#registerDoMC(cores=5)

source('lib_sparseHP.R')

######################################################
### UK specific parameters                         ###
######################################################

country.name = "UK"     # The country name should match with the one in the 'coronavirus' R package
final.date = "2020-06-08"   # Fianl date of the data: "YYYY-MM-DD". Default is today
if (is.null(final.date)) final.date=Sys.Date()

Np = 328200000          # total population in US
min_cum_case = 100      # The first date of the analysis begins when the number of cumulative cases reaches 100
lockdown = "2020-03-24" # March 24th: Lockdown measures begin
ghat = 1/18             # gamma

bwl = 3                 # backward window length for moving average

output_file_name = paste0("../results/",country.name,"_",final.date,".out")
sink(output_file_name)


######################################################
### load data                                      ###
######################################################

#df = data_extract(country.name, final.date)
df = read.csv(paste0("../data/covid_",country.name,"_2020-06-08.csv"))

######################################################
### Constructing key variables                     ###
######################################################

BetaData = gen_vars(date=as.Date(df$date),
                    confirmed=df$confirmed,
                    death=df$death,
                    recovered=df$recovered,
                    population=Np,
                    min_cum_case=min_cum_case)

date = BetaData$date
DeltaC = BetaData$DeltaC
I_lag = BetaData$I_lag
S_lag = BetaData$S_lag
Y = BetaData$Y

# N.B. The first observation of Y=beta(t) is at "first_period+1" by variable construction.

Y = movavg(Y, bwl, type="s")

# simple moving average of Y_t using back window length = bwl (e.g. Y_7 = (Y_7+Y_6+Y_5)/3 if bwl=3)

ly = log(Y)

skew = skewness(ly)
kurt = kurtosis(ly)
M = max(abs(diff(x = ly, differences = 2)))

######################################################
### construct lockdown variables                   ###
######################################################

n = length(ly)   # sample size
date_seq = 1:n

lockdown_date = date_seq[date==lockdown]     # define lockdown_date to be consistent with date_seq
lockdown_indicator = (date_seq > lockdown_date)
lockdown_X = (date_seq-lockdown_date)*lockdown_indicator


# Piecewise time trend: pre-lockdown: constant; post-lockdown: linear time trend

   lrg    = lm(ly~lockdown_X)
 lresid   = lrg$residuals
 lyhat    = lrg$fitted.values
 lRhat    = exp(lyhat)/ghat
 R        = Y/ghat
 lRresid  = R - lRhat

######################################################
### Selection of Method for Filtering              ###
######################################################

# options below will be turned off if zero is selected and on if one is selected
 selection_L2L0c = 1
 selection_L2 = 0
 selection_L1 = 1
 selection_SQRT_L1= 0

 draw_L2L0c_L1_togher = 1

# switch for cross-validation

  cross_validation = 0
  tuning_selection = 1

############################################### ERASE LATER ###################################################
#  fidelity = 1.78015                # fidelity from L2/L0 results
############################################### ERASE LATER ###################################################

#####################################################
### Sparse HP                                     ###
### L2/L0-constrained Trend Filtering             ###
#####################################################

if  (selection_L2L0c == 1){
  cat("------------------------------------- \n")
  cat("Sparse HP \n")
  cat("------------------------------------- \n")

  file_name = paste0("SparseHP_",country.name,"_",final.date)

  # Set the tuning parameters
    if (cross_validation == 0){
      l0constraint = 4
      l2penalty = 32
    }
    else{
      l0constraint_set = c(2:4)
      l2penalty_set = c(2^(seq(0,5,1)))

      n_l0 = length(l0constraint_set)
      n_l1 = length(l2penalty_set)
      cv_file_name = paste0(file_name,'_CV')
      #cv_file_name = paste0(file_name,'_CV','_',l0constraint_set)
      tic("cv")
      result <- l0tfc_cv(y=ly, l0constraint_set, M, l2penalty_set)
      comp_time <- toc()

      min_l2_l0 = result[which.min(result$obj_val),c(1:2)]

      l2penalty = min_l2_l0[1,1]
      l0constraint = min_l2_l0[1,2]

      cv_save=list(result=result, comp_time=comp_time)
      save.image(paste0('../results/',cv_file_name,'.RData'))

      # draw the cross-validation result graph
      result[,2] = as.factor(result[,2])
      f = ggplot(data=result, aes(x=l2, y=obj_val, group=l0, color=l0) ) +
        geom_line() +
        scale_color_viridis(discrete = TRUE) +
        ggtitle(paste0(country.name, ": Leave-one-out Cross-validation")) +
        ylab("Objective Function") +
        xlab(bquote(lambda)) +
        theme(axis.title = element_text(size = 15)) +
        guides(color=guide_legend(title=bquote(kappa))) +
        theme(legend.title = element_text(size=14)) +
        #scale_x_continuous(trans = "log2", breaks=c(0.125,0.25,0.5,1,2,4,8,16,32)) +
        scale_x_continuous(trans = "log2") +
        geom_vline(xintercept = min_l2_l0$l2, linetype="dashed", color="red")

      print(f)
      ggsave(paste0('../results/',cv_file_name,"_Graph.pdf"))

    }

    # Estimate L2/L0-constrained filtering
    l0_results = l0tfc(y=ly,l0constraint=l0constraint,M=M,l2penalty=l2penalty)
        l0_obj = l0_results$objval
        l0_est = l0_results$x

       betaHat = l0_est[1:length(ly)]
          zHat = l0_est[(length(ly)+1):length(l0_est)]
         SHP_kink_dates = date[zHat==1] + 1          # Double difference makes one lag. Should add 1 more day
         cat('\n\n Kink dates: \n')
         cat('---------------------------------------------------------------- \n')
         print(SHP_kink_dates)
         cat('---------------------------------------------------------------- \n')
         cat('\n\n \n')
         SHP_kink_dates=c(SHP_kink_dates)

         filter_yhat   = betaHat
         filter_resid  = ly - betaHat
         filter_Rhat   = exp(betaHat)/ghat
         filter_Rresid = R - filter_Rhat
         fidelity      = sum(filter_resid^2)
         cat('---------------------------------------------------------------- \n')
         cat('\n\n fidelity : ',fidelity,'\n' )
         cat('---------------------------------------------------------------- \n')
         cat('\n\n \n')



        data = data.frame(
           date,
           lockdown_date,
           R,            # R_t original scale
           ly,           # log(Y_t) = beta_t
           # linear fit reulsts
           lyhat,        # fitted value of linear reg, log scale
           lRhat,        # fitted value of linear reg, R_t scale
           lresid,       # residual of linear fit, log scale
           lRresid,      # residual of linear fit, R_t scale
           # filtering results
           filter_yhat,  # fitted value, log scale
           filter_Rhat,  # fitted value, R_t scale
           filter_resid, # residuals, log scale
           filter_Rresid # residuals, R_t scale
          #kinks         # kink dates
         )

        if (draw_L2L0c_L1_togher==0) {
         draw_fit_R0(data, country.name, file_name)
         draw_resid_R0(data, country.name, file_name)
         draw_fit_beta(data, country.name, file_name)
         draw_resid_beta(data, country.name, file_name)
    }
}

  #####################################################
  ### L1                                            ###
  #####################################################

  if  (selection_L1 == 1){
    cat("------------------------------------- \n")
    cat("L1 \n")
    cat("------------------------------------- \n")

    file_name = paste0("L1_",country.name,"_",final.date)

    # Set the tuning parameters
    if (tuning_selection == 0){
      ld = 0.8                 # User can set the tuning parameter lambda
    }
    else{
      ld_set = seq(0,10,0.1)
      l1_ld_out <- l1tf_tuning(y=ly, fidelity, ld_set, country.name, file_name)
      ld = l1_ld_out$opt_ld
    }

    l1_results = l1tf(y=ly,lambda=ld)
    betaHat = l1_results$betahat
    L1_kink_dates =  date[abs(diff(x=betaHat, differences=2))>1e-6] + 1   # Double difference makes one lag. Should add 1 more day
    cat('\n\n Kink dates: \n')
    cat('---------------------------------------------------------------- \n')
    print(L1_kink_dates)
    cat('---------------------------------------------------------------- \n')
    cat('\n\n \n')


    L1_filter_yhat   = betaHat
    L1_filter_resid  = ly - betaHat
    L1_filter_Rhat   = exp(betaHat)/ghat
    L1_filter_Rresid = R - filter_Rhat

    if (draw_L2L0c_L1_togher == 1) {
      data = cbind(data,
                   L1_filter_yhat,
                   L1_filter_resid,
                   L1_filter_Rhat,
                   L1_filter_Rresid
                   )
      #data$L1_kink_dates = as.Date(L1_kink_dates)
      draw_fit_R0_SHP_L1(data, SHP_kink_dates, L1_kink_dates, country.name, paste0('SHP_',file_name))
      draw_resid_R0_SHP_L1(data, SHP_kink_dates, L1_kink_dates,  country.name, paste0('SHP_',file_name))
      draw_fit_beta_SHP_L1(data, SHP_kink_dates, L1_kink_dates, country.name, paste0('SHP_',file_name))
      draw_resid_beta_SHP_L1(data, SHP_kink_dates, L1_kink_dates, country.name, paste0('SHP_',file_name))
    } else {
    data = data.frame(
      date,
      lockdown_date,
      R,            # R_t original scale
      ly,           # log(Y_t) = beta_t
      # SHP-filter fit reulsts
      lyhat,        # fitted value of linear reg, log scale
      lRhat,        # fitted value of linear reg, R_t scale
      lresid,       # residual of linear fit, log scale
      lRresid,      # residual of linear fit, R_t scale
      # L1-filtering results
      filter_yhat = L1_filter_yhat,  # fitted value, log scale
      filter_Rhat = L1_filter_resid,  # fitted value, R_t scale
      filter_resid = L1_filter_Rhat, # residuals, log scale
      filter_Rresid =L1_filter_Rresid# residuals, R_t scale
    )

    draw_fit_R0(data, country.name, file_name)
    draw_resid_R0(data, country.name, file_name)
    draw_fit_beta(data, country.name, file_name)
    draw_resid_beta(data, country.name, file_name)
  }

}

#####################################################
### HP                                            ###
#####################################################

if  (selection_L2 == 1){
  cat("------------------------------------- \n")
  cat("HP \n")
  cat("------------------------------------- \n")

  file_name = paste0("HP_",country.name,"_",final.date)

  # Set the tuning parameters
  if (cross_validation == 0){
    ld = 100
  }
  else{
    ld_set = seq(120,200,5)
    l2_ld_out <- l2tf_tuning(y=ly, fidelity, ld_set, country.name, file_name)
    ld = l2_ld_out$opt_ld
  }
  l2_results = l2tf(y=ly,lambda=ld)
  betaHat = l2_results$betahat

  filter_yhat   = betaHat
  filter_resid  = ly - betaHat
  filter_Rhat   = exp(betaHat)/ghat
  filter_Rresid = R - filter_Rhat

  data = data.frame(
    date,
    lockdown_date,
    R,            # R_t original scale
    ly,           # log(Y_t) = beta_t
    # linear fit reulsts
    lyhat,        # fitted value of linear reg, log scale
    lRhat,        # fitted value of linear reg, R_t scale
    lresid,       # residual of linear fit, log scale
    lRresid,      # residual of linear fit, R_t scale
    # filtering results
    filter_yhat,  # fitted value, log scale
    filter_Rhat,  # fitted value, R_t scale
    filter_resid, # residuals, log scale
    filter_Rresid # residuals, R_t scale
  )

  draw_fit_R0(data, country.name, file_name)
  draw_resid_R0(data, country.name, file_name)
  draw_fit_beta(data, country.name, file_name)
  draw_resid_beta(data, country.name, file_name)
}


#####################################################
### SQRT L1                                            ###
#####################################################

if  (selection_SQRT_L1 == 1){
  cat("------------------------------------- \n")
  cat("SQRT L1 \n")
  cat("------------------------------------- \n")

  file_name = paste0("SQRT_L1_",country.name,"_",final.date)

  # Set the tuning parameters
  if (cross_validation == 0){
    ld = 1.9                  # User can set the tuning parameter lambda
  }
  else{
    ld_set = seq(0,10,0.1)
    sqrt_l1_ld_out <- sqrt_l1tf_tuning(y=ly, fidelity, ld_set, country.name, file_name)
    ld = sqrt_l1_ld_out$opt_ld
  }


  sqrt_l1_results = sqrt_l1tf(y=ly,lambda=ld)
  betaHat = sqrt_l1_results$betahat
  kink_dates =  date[abs(diff(x=betaHat, differences=2))>1e-6]
  cat('\n\n Kink dates: \n')
  cat('---------------------------------------------------------------- \n')
  print(kink_dates)
  cat('---------------------------------------------------------------- \n')
  cat('\n\n \n')
  filter_yhat   = betaHat
  filter_resid  = ly - betaHat
  filter_Rhat   = exp(betaHat)/ghat
  filter_Rresid = R - filter_Rhat

  data = data.frame(
    date,
    lockdown_date,
    R,            # R_t original scale
    ly,           # log(Y_t) = beta_t
    # linear fit reulsts
    lyhat,        # fitted value of linear reg, log scale
    lRhat,        # fitted value of linear reg, R_t scale
    lresid,       # residual of linear fit, log scale
    lRresid,      # residual of linear fit, R_t scale
    # filtering results
    filter_yhat,  # fitted value, log scale
    filter_Rhat,  # fitted value, R_t scale
    filter_resid, # residuals, log scale
    filter_Rresid # residuals, R_t scale
  )

  draw_fit_R0(data, country.name, file_name)
  draw_resid_R0(data, country.name, file_name)
  draw_fit_beta(data, country.name, file_name)
  draw_resid_beta(data, country.name, file_name)
}

  cat("------------------------------------- \n")
  cat("Calculate the Growth Rate of R_0(t) \n")
  cat("------------------------------------- \n")

  g.logbeta = 100*diff(filter_yhat, differences=1) / filter_yhat[1:(n-1)]
  g.R = 100*diff(filter_Rhat, differences=1) / filter_Rhat[1:(n-1)]
  print(factor(round(g.R,2)))
  print(factor(round(diff(g.R,1),2)))
sink()