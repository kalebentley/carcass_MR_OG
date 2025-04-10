#---------------------------------------------------------------------------------------------------------- -
#                                                   ----
#---------------------------------------------------------------------------------------------------------- -

# Load functions
sapply(FUN = source, paste(getwd(), "functions", list.files("functions"), sep="/"))
logit<-function(x){log(x/(1-x))}

# Install/Load packages
package_list<-c("here", "gplots", "tidyverse", "dataRetrieval", "RColorBrewer", "MuMIn", "openxlsx" 
                , "RMark" , "R2jags", "tidybayes", "MCMCvis", "patchwork", "ggthemes", "glue", "coda") 
install_or_load_pack(package_list)

# Specify model files of interest
models<-c(
    "sst_model.r"
  # , "stt_model.r"
  # , "tst_model.r"
   ,"ttt_model.r"
  # , "timeseries_ttt_model.r",
  # , "timeseries_ttt_model_days.r",
  # , "timeseries_ttt_model_mvn_v2.r",#indexing requires same number of periods for all strata; Wishart Prior
  # , "timeseries_ttt_model_mvn_v3.r",#indexing requires same number of periods for all strata; only loops through  periods with non-0's for R and T; Wishart Prior..SOMETHING WRONG WITH RESULTS
  # , "timeseries_ttt_model_mvn_v4.r",#indexing requires same number of periods for all strata only loops through  periods with non-0's for R and T; Non-Wishart Prior (prior on hyperpriors of Sigma)..SOMETHING WRONG WITH RESULTS
  # , "timeseries_ttt_model_mvn_v5_kale.r", #indexing requires same number of periods for all strata but only loops through  periods with non-0's for R and T; Wishart Prior; mvn prior for first state of log-delta (same intercept for each strata!) ..SOMETHING WRONG WITH RESULTS
  # , "timeseries_ttt_model_v2index.r",#indexing requires same number of periods for all strata; only loops through non-0's for R and R Wishart Prior..SOMETHING WRONG WITH RESULTS
  # , "timeseries_ttt_model_mvn_v3_Nsuper.r",#indexing requires same number of periods for all strata but only loops through non-0's for R and T; Wishart Prior, N super is hierarchical prior...SOMETHING WRONG WITH RESULTS
  # ,  "timeseries_ttt_model_mvn_v2_Nsuper.r",#indexing requires same number of periods for all strata and no 0's for R and T; Wishart Prior; Nsuper is hierarchical prior
  # , "timeseries_ttt_model_mvn_v2_Nsuper_days.r"#indexing requires same number of periods for all strata; Wishart Prior
)

# Load M-R dataset via saved .rds file
location<-c("NFLewis")
species<-c("Chinook")
JS_stats_ALL<-readRDS(glue("data/Lewis_2017_JS_stats_threeGroups_noPooling.rds"))

# Collapse 3rd dimension of JS_stats into a single grouping using "collapse_third_dimension" function
(collapsed_df <- collapse_third_dimension(JS_stats_ALL))
(year<-year(min(collapsed_df$dates)))

# Define "good" JS vectors of data based on rule-sets
good_s_T_mat_1 = apply(as.matrix(collapsed_df$zi + collapsed_df$mi),2,function(x) c(which(x>0),rep(NA,length(which(x==0)))))
good_s_T_mat_2 = apply(apply(good_s_T_mat_1,2,function(x) ifelse(x%in%c(2:(nrow(collapsed_df)-1)),x,NA)),2,function(x) sort(x,na.last=TRUE))
good_s_T_final = apply(good_s_T_mat_2,2,function(x) length(x[!is.na(x)]))  

good_s_R_mat_1<-apply(as.matrix(collapsed_df$Ri),2,function(x) c(which(x>0),rep(NA,length(which(x==0)))))
good_s_R_mat_2=apply(apply(good_s_R_mat_1,2,function(x) ifelse(x%in%c(1:(nrow(collapsed_df)-1)),x,NA)),2,function(x) sort(x,na.last=TRUE))
good_s_R_final<-apply(good_s_R_mat_2,2,function(x) length(x[!is.na(x)]))

# Create data object for JAGS model
MR_data<-
  list(
    a=0.5
  , b=0.5
  , alpha = 0.5
  , LL = sum(collapsed_df$ui)
  , UL = sum(collapsed_df$ui) * 20
  , s = nrow(collapsed_df)
  #, num_strata = dim(n)[2]
  , uTot = sum(collapsed_df$ui)
  , n = collapsed_df$ni
  , R = collapsed_df$Ri
  , u = collapsed_df$ui
  , r = collapsed_df$ri
   ,m = collapsed_df$mi
  , T = collapsed_df$zi + collapsed_df$mi
  , days=c(NA,as.numeric(diff(collapsed_df$dates))) #number of days between last period and current period...for first period, value is NA
  , good_s_T_mat = good_s_T_mat_2
  , good_s_R_mat = good_s_R_mat_2
  , good_s_T = good_s_T_final
  , good_s_R = good_s_R_final
  , logit_phi_1_mu_prior = logit(exp(log(0.5)/mean(c(NA,as.numeric(diff(collapsed_df$dates))),na.rm=TRUE)))
  )

if(is.null(dim(MR_data$R)[2])==TRUE){
  MR_data[lapply(MR_data,class)=="matrix"]<-lapply(MR_data[lapply(MR_data,class)=="matrix"],as.vector)
  MR_data$alpha<-rep(1,MR_data$s)
}

# Create inits for uTot
uTot<-MR_data$uTot
inits_func<-function(){
  inits<-list(
    "Ntot" = uTot * runif(1,1,5)
  )
}

# Run JS models of interest
for (i in 1:length(models)){
  start.time<-Sys.time()
  print(start.time)
  out<- 
    jags.parallel(
      data = MR_data
      , model.file = glue("models/{models[i]}")
      , n.chains= 4
      , n.thin= 40
      , n.burnin= 30000
      , n.iter = 120000 
      , parameters.to.save=c("p",
                           "phi",
                           "logit_phi",
                           "pent",
                           "log_delta",
                           "delta",
                           "Nsuper",
                           "Bstar",
                           "rho_phi",
                           "rho_p",
                           "rho_pent",
                           "sigma_phi_process",
                           "sigma_p_process",
                           "sigma_pent_process"
      )
      #, inits = inits_func
      # , n.chains = n_chains
      # , n.iter = n_iter
      # , n.burnin = n_burn
      # , n.thin = n_thin
      # , parameters.to.save=params_to_save
    )
  end.time<-Sys.time()
  model.duration<-end.time-start.time
  print(model.duration)
  print(out,dig=3)
  print(i)
  write.csv(x = out$BUGSoutput$summary, file = glue("output/JS-summary-results_{location}_{species}_{year}_{str_remove(models[i], '\\\\.r$')}.csv"))
}


#plots
library(RColorBrewer)
pdf(paste(Sys.Date(), models[i],"p_phi_pent.pdf",sep=""))
cols<-brewer.pal(length(win.data$uTot),"Blues")
par(mfrow=c(1,1))
x<-matrix(rep(1:win.data$s, length(win.data$uTot)),ncol=length(win.data$uTot))
#par(mfrow=c(2,1),mar=c(4,4,4,4))
matplot(out$BUGSoutput$median$p,type="l",col=cols[1:length(win.data$uTot)],lty=1,lwd=2,ylim=c(0,1),xlim=c(1,win.data$s), ylab="Capture")
legend("topright",legend=names(win.data$uTot),lty=rep(1,length(win.data$uTot)),lwd=2,col=cols[1:length(win.data$uTot)])
#matplot(y=win.data$m/win.data$T,x=x,type="l",col=cols[1:length(win.data$uTot)],lty=2,lwd=2,ylim=c(0,1),xlim=c(1,win.data$s),ylab ="RAW (m/T)")
#legend("topright",legend=names(win.data$uTot),lty=rep(2,length(win.data$uTot)),lwd=2,col=cols[1:length(win.data$uTot)])
#par(mfrow=c(2,1),mar=c(4,4,4,4))
matplot(out$BUGSoutput$median$phi,type="l",col=cols[1:length(win.data$uTot)],lty=1,lwd=2,ylim=c(0,1),xlim=c(1,c(win.data$s-1)),ylab="Survival")
legend("topright",legend=names(win.data$uTot),lty=rep(1,length(win.data$uTot)),lwd=2,col=cols[1:length(win.data$uTot)])
#matplot(y=win.data$r/win.data$R,x=x,type="l",col=cols[1:length(win.data$uTot)],lty=2,lwd=2,ylim=c(0,1),xlim=c(1,c(win.data$s-1)),ylab ="RAW (r/R)")
#legend("topright",legend=names(win.data$uTot),lty=rep(2,length(win.data$uTot)),lwd=2,col=cols[1:length(win.data$uTot)])
#par(mfrow=c(1,1))
matplot(out$BUGSoutput$median$pent,type="l",col=cols[1:length(win.data$uTot)],lty=1,lwd=2,ylim=c(0,1),ylab="Pent")
legend("topright",legend=names(win.data$uTot),lty=1,lwd=2,col=cols[1:length(win.data$uTot)])
dev.off()