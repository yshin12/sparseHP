######################################################
### Construct Variables and First Look             ###
### Parametric mdoel fitting                       ###
### Updated on 6 June 2020                         ###
######################################################

rm(list=ls())   # initialization

library(ggplot2)
library(quantreg)
library(pracma)
library(CVXR)
library(gurobi)
library(coronavirus)
library(moments)


source('lib_sparseHP.R')               # call library  that contains required functions


######################################################
### US specific parameters                         ###  
######################################################

country.name = "US"         # The country name should match with the one in the 'coronavirus' R package
final.date = "2020-06-08"   # Fianl date of the data: "YYYY-MM-DD". Default is today
if (is.null(final.date)) final.date=Sys.Date()

Np = 328200000          # total population in US
min_cum_case = 100      # The first date of the analysis begins when the number of cumulative cases reaches 100
lockdown = "2020-03-30" # March 30th: 30 statewide stay-at-home orders
ghat = 1/18             # gamma   

bwl = 3                 # backward window length for moving average


file_name = paste0(country.name,"_",final.date)  

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
                  
######################################################
### Plots of key variables                         ###
######################################################     

data = data.frame(
        date = date,
         ts1 = DeltaC,
         ts2 = I_lag,
         ts3 = S_lag,
         ts4 = Y
     )

# plot of DeltaC
 
 ggplot(data, aes(x=date, y=ts1)) +
   geom_line(color = 'maroon') +
   xlab("Date") + 
   ylab(bquote(Delta~ C[t])) +
   geom_vline(xintercept = date[date==lockdown], col="orange") +
   scale_y_continuous(labels = scales::scientific) +
   ggtitle("Daily positives") +
   theme_minimal()
 
 gg_file_name = paste("../results/DeltaC",file_name,".pdf",sep="")
 
 ggsave(gg_file_name, width=6, height=4, units="in")
 
# plot of I_lag
 
 ggplot(data, aes(x=date, y=ts2)) +
   geom_line(color = 'maroon') +
   xlab("Date") + 
   ylab(bquote(I[t-1])) +
   geom_vline(xintercept = date[date==lockdown], col="orange") +
   ggtitle("Lagged cumulative infectives") +
   theme_minimal()
 
 gg_file_name = paste("../results/I_lag",file_name,".pdf",sep="")
 
 ggsave(gg_file_name, width=6, height=4, units="in")
 
# plot of S_lag
 
 ggplot(data, aes(x=date, y=ts3)) +
   geom_line(color = 'maroon') +
   xlab("Date") + 
   ylab(bquote(S[t-1])) +
   geom_vline(xintercept = date[date==lockdown], col="orange") +
   ggtitle("Lagged cumulative susceptibles") +
   theme_minimal()
 
 gg_file_name = paste("../results/S_lag",file_name,".pdf",sep="")
 
 ggsave(gg_file_name, width=6, height=4, units="in")
 
 
# plot of Y
 
 ggplot(data, aes(x=date, y=ts4)) +
   geom_line(color = 'maroon') +
   xlab("Date") + 
   ylab(bquote(beta[t])) +
   geom_vline(xintercept = date[date==lockdown], col="orange") +
   ggtitle(bquote("Time-varying " ~ beta[t])) +
   theme_minimal()
 
 gg_file_name = paste("../results/beta",file_name,".pdf",sep="")
 
 ggsave(gg_file_name, width=6, height=4, units="in")

 # plot of log(Y)
 
 ggplot(data, aes(x=date, y=log(ts4))) +
   geom_line(color = 'maroon') +
   xlab("Date") + 
   ylab(bquote(log ~ beta[t] )) +
   geom_vline(xintercept = date[date==lockdown], col="orange") +
   ggtitle(bquote("Time-varying " ~ log ~ beta[t])) +
   theme_minimal()
 
 gg_file_name = paste("../results/logbeta",file_name,".pdf",sep="")
 
 ggsave(gg_file_name, width=6, height=4, units="in")
 
 
 
######################################################
### construct lockdown variables                   ###
######################################################          
 
 n = length(ly)   # sample size 
 date_seq = 1:n
 
      lockdown_date = date_seq[date==lockdown]     # define lockdown_date to be consistent with date_seq
 lockdown_indicator = (date_seq > lockdown_date)
         lockdown_X = (date_seq-lockdown_date)*lockdown_indicator
 
 
 # Piecewise time trend: pre-lockdown: constant; post-lockdown: linear time trend
 
 lrg = lm(ly~lockdown_X)
 resid = lrg$residuals
 lyhat = lrg$fitted.values
   
######################################################
### plot of log(beta_t) & parametric model         ###
######################################################     
 
 data = data.frame(
   date = date,
   value = ly,
   yhat = lyhat,
   resid = resid,
   y_level = Y/ghat,
   yhat_level = exp(lyhat)/ghat,
   resid_level = Y/ghat - exp(lyhat)/ghat
 )

 ggplot(data, aes(x=date, y=value)) +
   geom_line(color="gray50") + 
   geom_line(mapping = aes(x = date, y = yhat), color = 'red', size = 0.7) +
   xlab("Date") + 
   ylab(bquote(log ~ beta[t])) +
   geom_vline(xintercept = date[date==lockdown], col="orange") +   
   ggtitle("Parametric fitting of " ~ log ~ beta[t]) +
   theme_minimal() 
 
 
 gg_file_name = paste("../results/logbeta_para",file_name,".pdf",sep="")
 
 ggsave(gg_file_name, width=6, height=4, units="in")
 
# Plot of R_0(t) = beta(t)/gamma
 
 ggplot(data, aes(x=date, y=y_level)) +
   geom_line(color="gray50") + 
   geom_line(mapping = aes(x = date, y = yhat_level), color = 'red', size = 0.7) +
   xlab("Date") + 
   ylab(bquote(R[0])) +
   geom_vline(xintercept = date[date==lockdown], col="orange") +   
   ggtitle(bquote("Parametric fitting of " ~ R[0](t) ~ "=" ~ beta[t] / gamma)) +
   theme_minimal() 
 
 
 gg_file_name = paste("../results/R0_para",file_name,".pdf",sep="")
 
 ggsave(gg_file_name, width=6, height=4, units="in")
 
######################################################
### plot of residuals from the parametric model    ###
######################################################       
 
 ggplot(data, aes(x=date, y=resid)) +
   geom_line(color = 'red', size = 0.7)  +
   geom_area(fill="red", alpha=0.3) +
   xlab("Date") + 
   ylab(bquote(log ~ beta[t])) +
   ggtitle(bquote("Residuals from parametric fitting of " ~ log ~ beta[t])) +
   theme_minimal()
 
 gg_file_name = paste("../results/logbeta_resid",file_name,".pdf",sep="")
 
 ggsave(gg_file_name, width=6, height=4, units="in")
 
 ggplot(data, aes(x=date, y=resid_level)) +
   geom_line(color="red", size = 0.7) + 
   geom_area(fill="red", alpha=0.3) +
   xlab("Date") + 
   ylab(bquote(R[0])) +
   ggtitle(bquote("Residuals from parametric fitting of " ~ R[0](t) ~ "=" ~ beta[t] / gamma)) +
   theme_minimal()

 gg_file_name = paste("../results/R0_resid",file_name,".pdf",sep="")
 
 ggsave(gg_file_name, width=6, height=4, units="in")
 
 
 