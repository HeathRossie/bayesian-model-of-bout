rm(list = ls())
# import library
library(ggplot2)
library(dplyr)
library(rstan)
library(gridExtra)

### defining functions ------------------------------------------------------------------------------------#
#-----response generation according to Brackney et al. (2011)-----#
bout.sim <- function(VR.value, b.bout.VR, b.bout.VI,  PRP, ref.per, q.VR = 5, q.VI = 3,  w.bout.VR = 0.15, w.bout.VI = 0.25){
  
  bins <- function(data, bin){
    (floor(data / bin * 10^0) /10^0)
  }
  
  dat.para <- data.frame(VR.value = VR.value,
                         b.bout.VR = b.bout.VR,
                         b.bout.VI = b.bout.VI,
                         q.VR = q.VR,
                         q.VI = q.VI,
                         w.bout.VR = w.bout.VR,
                         w.bout.VI = w.bout.VI,
                         ref.per = ref.per)
  
  # simulation
  dat.sim <- NULL
  for(i in 1:5){
    IRT <- NULL
    n <- 1
    trial.num <- 0
    trial <- NULL
    
    while(trial.num < 81){
      # PRP (PRPs are not assumed in default)
      if(rbinom(n = 1, size = 1, prob = 1 / dat.para$VR.value) == 1){
        IRT[n] <- sample(PRP, 1, replace = T)
        trial.num <- trial.num + 1
        trial[n] <- trial.num
        n <- n + 1
      }
      
      if(rbinom(n = 1, size = 1, prob = 1 / dat.para$q.VR) == 1){
        # sampling between-bout IRT
        pro.IRT <- rexp(n = 1, rate = 1 / dat.para$b.bout.VR)
        
        # resampling if the sampled IRT exceeded refractory periods
        while(pro.IRT <= dat.para$ref.per){
          pro.IRT <- rexp(n = 1, rate = 1 / dat.para$b.bout.VR)
        }
        IRT[n] <- pro.IRT
        trial[n] <- trial.num
        n <- n + 1
      } else {
        # sampling within-bout IRT
        pro.IRT <- rexp(n = 1, rate = 1 / dat.para$w.bout.VR)
        
        # resampling if the sampled IRT exceeded refractory periods
        while(pro.IRT <= dat.para$ref.per){
          pro.IRT <- rexp(n = 1, rate = 1 / dat.para$w.bout.VR)
        }
        
        IRT[n] <- pro.IRT
        trial[n] <- trial.num
        n <- n + 1
      }
    }
    dat.sim[[i]] <- data.frame(session = i,
                               IRT = IRT,
                               trial = trial)
  }
  
  VR.IRT <- do.call(rbind, dat.sim)
  
  # binning process
  list.VR <- lapply(split(VR.IRT, list(VR.IRT$trial, VR.IRT$session), drop = T), function(data){
    data$timer <- cumsum(data$IRT)
    data$bin <- bins(data$timer, 1)
    content <- data.frame(bin = unique(data$bin),
                          resp = as.vector(table(data$bin)))
    container <- data.frame(bin = 0:max(data$bin),
                            resp = 0)
    container$resp[content$bin + 1] <- as.integer(content$resp)
    data.frame(session =unique(data$session),
               container,
               trial = unique(data$trial))
  })
  
  VR.bin<- data.frame(do.call(rbind, list.VR))
  
  # calculating VI value
  VI.value <- 60 / ((60 / mean(VR.IRT$IRT)) / dat.para$VR.value)
  
  #VI simulation
  dat.sim <- NULL
  for(i in 1:5){
    IRT <- NULL
    n <- 1
    trial.num <- 0
    trial <- NULL
    
    while(trial.num < 81){
      
      # PRP (PRPs are not assumed in default)
      IRI <- rexp(n = 1, rate = 1 / VI.value)
      interval <- 0
      IRT[n] <- sample(PRP, 1, replace = T)
      interval <- interval + IRT[n]
      trial[n] <- trial.num
      n <- n + 1
      
      while(interval <= IRI){
        if(rbinom(n = 1, size = 1, prob = 1 / dat.para$q.VI) == 1){
          
          # sampling between-bout IRT
          pro.IRT <- rexp(n = 1, rate = 1 / dat.para$b.bout.VI)
          
          # resampling if the sampled IRT exceeded refractory periods
          while(pro.IRT <= dat.para$ref.per){
            pro.IRT <- rexp(n = 1, rate = 1 / dat.para$b.bout.VI)
          }
          IRT[n] <- pro.IRT
          interval <- interval + pro.IRT
          trial[n] <- trial.num
          n <- n + 1
        } else {
          # sampling within-bout IRT
          pro.IRT <- rexp(n = 1, rate = 1 / dat.para$w.bout.VI)
          
          # resampling if the sampled IRT exceeded refractory periods
          while(pro.IRT <= dat.para$ref.per){
            pro.IRT <- rexp(n = 1, rate = 1 / dat.para$w.bout.VI)
          }
          IRT[n] <- pro.IRT
          interval <- interval + pro.IRT
          trial[n] <- trial.num
          n <- n + 1
        }
      }
      trial.num <- trial.num + 1
    }
    dat.sim[[i]] <- data.frame(session = i,
                               IRT = IRT,
                               trial = trial)
  }
  
  VI.IRT <- do.call(rbind, dat.sim)
  
  # binning process
  list.VI <- lapply(split(VI.IRT, list(VI.IRT$trial, VI.IRT$session), drop = T), function(data){
    data$timer <- cumsum(data$IRT)
    data$bin <- bins(data$timer, 1)
    content <- data.frame(bin = unique(data$bin),
                          resp = as.vector(table(data$bin)))
    container <- data.frame(bin = 0:max(data$bin),
                            resp = 0)
    container$resp[content$bin + 1] <- as.integer(content$resp)
    data.frame(session = unique(data$session),
               container,
               trial = unique(data$trial))
  })
  
  VI.bin <- data.frame(do.call(rbind, list.VI))
  
  dat.sim <- list(VR.IRT = VR.IRT,
                  VR.bin = VR.bin[VR.bin$trial != 0, ],
                  VI.IRT = VI.IRT,
                  VI.bin = VI.bin[VI.bin$trial != 0, ])
  
  return(dat.sim)
}

#-----modified plot function for visualzing results-----#
# function to process data for log-survivor plots
log.sur <- function(IRT.bin, bin, cond="VI"){
  
  bins <- function(data, bin){
    (floor(data / bin * 10^0) /10^0)
  }
  
  IRT.bin <- ((bins(IRT.bin, bin) + 1) / 10)* 2
  
  
  content <- data.frame(ser.num = sort(unique(IRT.bin))* 5,
                        IRT.bin = sort(unique(IRT.bin)),
                        num = as.vector(table(IRT.bin)))
  
  container <- data.frame(IRT.bin = seq(0.2, max(IRT.bin), 0.2),
                          num = 0)
  
  container$num[content$ser.num] <- content$num
  container$num <- container$num / sum(container$num)
  
  sur.plot <- data.frame(IRT.bin = container$IRT.bin,
                         sur.fun = 1 - cumsum(container$num),
                         cond = cond)
  
  return(sur.plot)
}

# plot function for log-survivor plots
plot_ls <- function(sur.plot, ymin=-5, ymax=0, xmin=0, xmax=100){
  pl <- ggplot(data = sur.plot, aes(x = IRT.bin, y = log10(sur.fun), linetype=cond)) +
    geom_line() + 
    xlab(NULL) + ylab(NULL) + TM + ylim(values=c(ymin,ymax)) +  xlim(values=c(xmin,xmax)) +
    theme(axis.text = element_text(colour="black"), 
          axis.ticks.length =  unit(.2, "cm"),
          legend.position = "none"
    )
  
  geom_line() + return(pl)
}

plot_mr <- function(respVR, binVR, respVI, binVI){
  datVR <- as.vector(tapply(respVR, binVR, mean))[1:used_bin]
  datVI <- as.vector(tapply(respVI, binVI, mean))[1:used_bin]
  
  ggplot() +
    geom_line(aes(x=1:used_bin, y=datVR)) + 
    geom_line(aes(x=1:used_bin, y=datVI), lty=2) +
    xlab(NULL) + ylab(NULL) + TM + ylim(values=c(0,4)) +
    scale_x_continuous(breaks = c(0,25,50,75)) + 
    theme(axis.text = element_text(colour="black"), 
          axis.ticks.length =  unit(.2, "cm"))
}

### simulation part ------------------------------------------------------------------------------------#

### generating data and fitting the model
bout_simulator <- function(b.bout){
  
  # setting a seed
  set.seed(1234)
  d_sim <- bout.sim(VR.value = 30, b.bout.VR = b.bout, b.bout.VI = b.bout, PRP = 0, ref.per = 0.08,
                    q.VR = 14.29, q.VI = 5.57,  w.bout.VR = 0.25, w.bout.VI = 0.35)
  
  d_VR <- d_sim$VR.bin
  d_VI <- d_sim$VI.bin
  d_VR$cond <- "VR"
  d_VI$cond <- "VI"
  d_input <- rbind(d_VR, d_VI)
  d_input$bin <- d_input$bin + 1
  used_bin <- 75
  d_input <- d_input[d_input$bin <= used_bin,]
  
  d.list <- list(N_VI=nrow(d_input[d_input$cond=="VI",]),
                 N_VR=nrow(d_input[d_input$cond=="VR",]),
                 resp_VI=d_input[d_input$cond=="VI",]$resp,
                 resp_VR=d_input[d_input$cond=="VR",]$resp,
                 time_VI=d_input[d_input$cond=="VI",]$bin,
                 time_VR=d_input[d_input$cond=="VR",]$bin,
                 N_time_VI=used_bin,
                 N_time_VR=used_bin)
  
  # fitting the model with STAN
  md <- stan_model('model.stan')
  fit <- vb(md, d.list, seed=1234, algorithm="meanfield")
  res <- rstan::extract(fit)
  
  return(list(d_sim = d_sim, input = d_input, res = res))
}

# Note that PRP = 0 means that post-reinforcement pause is not assumed in the simulation
# though PRP was sampled emprical distribution of the real data in the simulation of the article, 


### perform simulation
# set working directory to the location where the model code was saved
# setwd("~/Desktop") 
b.bout <- c(1, 2, 4, 8)
sim1 <- bout_simulator(b.bout[1])
sim2 <- bout_simulator(b.bout[2])
sim3 <- bout_simulator(b.bout[3])
sim4 <- bout_simulator(b.bout[4])

### Visualizing the results ------------------------------------------------------------------------------------#
### set the theme of plots
TM <- theme_classic() + theme(axis.text = element_text(family="Times New Roman",colour="black", size=rel(3)))

### depicting log-survivor plots
p1 <- log.sur(sim4$d_sim$VR.IRT$IRT, bin=0.08, cond="VR") %>% 
  rbind(., log.sur(sim4$d_sim$VI.IRT$IRT, 0.08, cond="VI")) %>% 
  mutate(., cond=factor(as.character(cond), levels=c("VR", "VI"))) %>% 
  plot_ls(., xmin=0, xmax=30, ymin=-3)
p2 <- log.sur(sim3$d_sim$VR.IRT$IRT, 0.08, cond="VR") %>% 
  rbind(., log.sur(sim3$d_sim$VI.IRT$IRT, 0.08, cond="VI")) %>% 
  mutate(., cond=factor(as.character(cond), levels=c("VR", "VI"))) %>% 
  plot_ls(., xmin=0, xmax=30, ymin=-3)
p3 <- log.sur(sim2$d_sim$VR.IRT$IRT, 0.08, cond="VR") %>% 
  rbind(., log.sur(sim2$d_sim$VI.IRT$IRT, 0.08, cond="VI")) %>% 
  mutate(., cond=factor(as.character(cond), levels=c("VR", "VI"))) %>% 
  plot_ls(., xmin=0, xmax=30, ymin=-3)
p4 <- log.sur(sim1$d_sim$VR.IRT$IRT, 0.08, cond="VR") %>% 
  rbind(., log.sur(sim1$d_sim$VI.IRT$IRT, 0.08, cond="VI")) %>% 
  mutate(., cond=factor(as.character(cond), levels=c("VR", "VI"))) %>% 
  plot_ls(., xmin=0, xmax=30, ymin=-3)
grid.arrange(p1, p2, p3, p4,ncol=1)

### the response rates as a function of elasped times in IRIs 
# equivalent figures of Figure 5ab in the article
# solid lines and dashed lines denote simulation data and means of posterior predictive distribution
# Note that bitonic function would not appear if the PRP=0 because of lack of post-reinforcement pauses
X <- 1:75
EXtract <- function(d, cond) d[d$cond==cond,]
Posmean <- function(d, cond){
  if(cond == "VR") pm <- apply(d$pred_VR2, 2, mean)
  if(cond == "VI") pm <- apply(d$pred_VI2, 2, mean)
  return(pm)
} 

# VR
grid.arrange(
  ggplot() + 
    geom_line(aes(X, y=tapply(sim4$d_sim$VR.bin$resp, sim4$d_sim$VR.bin$bin, mean)[1:75])) +
    geom_line(aes(X, y=tapply(Posmean(sim4$res, "VR"), EXtract(sim4$input, "VR")$bin, mean)), linetype=2) +
    xlab(NULL) + ylab(NULL) + TM + 
    scale_y_continuous(limits=c(0,4), breaks = c(0,1,2,3,4)) + 
    scale_x_continuous(breaks = c(0,25,50,75)) + 
    theme(axis.text = element_text(colour="black"), 
          axis.ticks.length =  unit(.2, "cm"))
  ,
  ggplot() + 
    geom_line(aes(X, y=tapply(sim3$d_sim$VR.bin$resp, sim3$d_sim$VR.bin$bin, mean)[1:75])) +
    geom_line(aes(X, y=tapply(Posmean(sim3$res, "VR"), EXtract(sim3$input, "VR")$bin, mean)), linetype=2) +
    xlab(NULL) + ylab(NULL) + TM + 
    scale_y_continuous(limits=c(0,4), breaks = c(0,1,2,3,4)) + 
    scale_x_continuous(breaks = c(0,25,50,75)) + 
    theme(axis.text = element_text(colour="black"), 
          axis.ticks.length =  unit(.2, "cm"))
  ,
  ggplot() + 
    geom_line(aes(X, y=tapply(sim2$d_sim$VR.bin$resp, sim2$d_sim$VR.bin$bin, mean)[1:75])) +
    geom_line(aes(X, y=tapply(Posmean(sim2$res, "VR"), EXtract(sim2$input, "VR")$bin, mean)), linetype=2) +
    xlab(NULL) + ylab(NULL) + TM + 
    scale_y_continuous(limits=c(0,4), breaks = c(0,1,2,3,4)) + 
    scale_x_continuous(breaks = c(0,25,50,75)) + 
    theme(axis.text = element_text(colour="black"), 
          axis.ticks.length =  unit(.2, "cm"))
  ,
  ggplot() + 
    geom_line(aes(1:75, y=tapply(sim1$d_sim$VR.bin$resp, sim1$d_sim$VR.bin$bin, mean)[1:75])) +
    geom_line(aes(X, y=tapply(Posmean(sim1$res, "VR"), EXtract(sim1$input, "VR")$bin, mean)), linetype=2) +
    xlab(NULL) + ylab(NULL) + TM + 
    scale_y_continuous(limits=c(0,4), breaks = c(0,1,2,3,4)) + 
    scale_x_continuous(breaks = c(0,25,50,75)) + 
    theme(axis.text = element_text(colour="black"), 
          axis.ticks.length =  unit(.2, "cm"))
  , nrow=4)


# VI
grid.arrange(
  ggplot() + 
    geom_line(aes(X, y=tapply(sim4$d_sim$VI.bin$resp, sim4$d_sim$VI.bin$bin, mean)[1:75])) +
    geom_line(aes(X, y=tapply(Posmean(sim4$res, "VI"), EXtract(sim4$input, "VI")$bin, mean)), linetype=2) +
    xlab(NULL) + ylab(NULL) + TM + 
    scale_y_continuous(limits=c(0,4), breaks = c(0,1,2,3,4)) + 
    scale_x_continuous(breaks = c(0,25,50,75)) + 
    theme(axis.text = element_text(colour="black"), 
          axis.ticks.length =  unit(.2, "cm"))
  ,
  ggplot() + 
    geom_line(aes(X, y=tapply(sim3$d_sim$VI.bin$resp, sim3$d_sim$VI.bin$bin, mean)[1:75])) +
    geom_line(aes(X, y=tapply(Posmean(sim3$res, "VI"), EXtract(sim3$input, "VI")$bin, mean)), linetype=2) +
    xlab(NULL) + ylab(NULL) + TM + 
    scale_y_continuous(limits=c(0,4), breaks = c(0,1,2,3,4)) + 
    scale_x_continuous(breaks = c(0,25,50,75)) + 
    theme(axis.text = element_text(colour="black"), 
          axis.ticks.length =  unit(.2, "cm"))
  ,
  ggplot() + 
    geom_line(aes(X, y=tapply(sim2$d_sim$VI.bin$resp, sim2$d_sim$VI.bin$bin, mean)[1:75])) +
    geom_line(aes(X, y=tapply(Posmean(sim2$res, "VI"), EXtract(sim2$input, "VI")$bin, mean)), linetype=2) +
    xlab(NULL) + ylab(NULL) + TM +
    scale_y_continuous(limits=c(0,4), breaks = c(0,1,2,3,4)) + 
    scale_x_continuous(breaks = c(0,25,50,75)) + 
    theme(axis.text = element_text(colour="black"), 
          axis.ticks.length =  unit(.2, "cm"))
  ,
  ggplot() + 
    geom_line(aes(X, y=tapply(sim1$d_sim$VI.bin$resp, sim1$d_sim$VI.bin$bin, mean)[1:75])) +
    geom_line(aes(X, y=tapply(Posmean(sim1$res, "VI"), EXtract(sim1$input, "VI")$bin, mean)), linetype=2) +
    xlab(NULL) + ylab(NULL) + TM + 
    scale_y_continuous(limits=c(0,4), breaks = c(0,1,2,3,4)) +
    scale_x_continuous(breaks = c(0,25,50,75)) + 
    theme(axis.text = element_text(colour="black"), 
          axis.ticks.length =  unit(.2, "cm"))
  , nrow=4)

### depictng mean response rates from simulation data and posterior predictive distribution
#  equivalent figures of Figure 5c in the article
get_rdat <- function(d){
  mrr <- function(r){
    rate <- sum(r)/length(r)
    rate <- rate 
    return(rate)
  }
  rdat <- data.frame(resp = c(tapply(d$input$resp,d$input$cond, mrr),
                              mean(d$res$pred_VI2),
                              mean(d$res$pred_VR2)
  ),
  minute = c(tapply(d$input$resp,d$input$cond, mean),
             mean(d$res$pred_VR2),
             mean(d$res$pred_VR2)
  ),
  cond = factor(rep(c("VI","VR"),2), levels=c("VR","VI")),
  data = factor(rep(c("real","pred"),each=2),levels=c("real","pred")))
  rownames(rdat) <- NULL
  return(rdat)
}

plt <- function(d){
  pl <- ggplot(d) + 
    geom_bar(aes(x=interaction(data,cond), y=resp, fill=data), stat = "identity", colour="black")  +
    scale_fill_manual(values=c("gray", "black")) +
    scale_x_discrete(label=c("VR", "VR", "VI", "VI")) +
    scale_y_continuous(expand = c(0,0), limits=c(0,4), breaks=c(0,1,2,3,4))  + 
    xlab(NULL) + ylab(NULL) + TM + 
    theme(
          axis.text.x = element_text(colour="black", size=rel(.7)), 
          axis.text.y = element_text(colour="black"), 
          axis.ticks.length =  unit(.2, "cm"),
          axis.ticks.x = element_blank(),
          plot.margin= unit(c(2, 1, 1, 1), "lines"))
  return(pl)
}

grid.arrange(plt(get_rdat(sim4)),
             plt(get_rdat(sim3)),
             plt(get_rdat(sim2)),
             plt(get_rdat(sim1)),
             nrow=4)

### means of parameter estimates 
# equivalent figures to Figure 6ab
# lambda
rpdat <- data.frame(
  cond = factor(c("yorked-VI","VR30"), levels=c("VR30","yorked-VI")),
  resp1 = c(mean(sim4$res$mean_lambda_VI),mean(sim4$res$mean_lambda_VR)),
  resp2 = c(mean(sim3$res$mean_lambda_VI),mean(sim3$res$mean_lambda_VR)),
  resp3 = c(mean(sim2$res$mean_lambda_VI),mean(sim2$res$mean_lambda_VR)),
  resp4 = c(mean(sim1$res$mean_lambda_VI),mean(sim1$res$mean_lambda_VR)),
  lower1 = c(quantile(sim4$res$mean_lambda_VI, 0.025), quantile(sim4$res$mean_lambda_VR, 0.025)),
  lower2 = c(quantile(sim3$res$mean_lambda_VI, 0.025), quantile(sim3$res$mean_lambda_VR, 0.025)),
  lower3 = c(quantile(sim2$res$mean_lambda_VI, 0.025), quantile(sim2$res$mean_lambda_VR, 0.025)),
  lower4 = c(quantile(sim1$res$mean_lambda_VI, 0.025), quantile(sim1$res$mean_lambda_VR, 0.025)),
  upper1 = c(quantile(sim4$res$mean_lambda_VI, 0.975), quantile(sim4$res$mean_lambda_VR, 0.975)),
  upper2 = c(quantile(sim3$res$mean_lambda_VI, 0.975), quantile(sim3$res$mean_lambda_VR, 0.975)),
  upper3 = c(quantile(sim2$res$mean_lambda_VI, 0.975), quantile(sim2$res$mean_lambda_VR, 0.975)),
  upper4 = c(quantile(sim1$res$mean_lambda_VI, 0.975), quantile(sim1$res$mean_lambda_VR, 0.975))
)

p1 <- ggplot(rpdat) + 
  geom_bar(aes(x=cond, y=resp1, fill=cond), stat = "identity", colour="black") + 
  geom_errorbar(aes(x=cond, ymin=lower1, ymax=upper1), width=.2, colour="black") +
  scale_fill_manual(values=c("gray", "gray15")) + 
  scale_y_continuous(expand = c(0,0), limits=c(0,4), breaks=c(0,1,2,3,4))  +
  xlab(NULL) + ylab(NULL) + TM +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="white"),
        axis.ticks.length =  unit(.2, "cm"),
        axis.ticks.x = element_blank(),
        plot.margin= unit(c(1, 1, -1, 1), "lines"))

p2 <- ggplot(rpdat) + 
  geom_bar(aes(x=cond, y=resp2, fill=cond), stat = "identity", colour="black") + 
  geom_errorbar(aes(x=cond, ymin=lower2, ymax=upper2), width=.2, colour="black") +
  scale_fill_manual(values=c("gray", "gray15")) + 
  scale_y_continuous(expand = c(0,0), limits=c(0,4), breaks=c(0,1,2,3,4))  +
  xlab(NULL) + ylab(NULL) + TM +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="white"),
        axis.ticks.length =  unit(.2, "cm"),
        axis.ticks.x = element_blank(),
        plot.margin= unit(c(1, 1, -1, 1), "lines"))

p3 <- ggplot(rpdat) + 
  geom_bar(aes(x=cond, y=resp3, fill=cond), stat = "identity", colour="black") + 
  geom_errorbar(aes(x=cond, ymin=lower3, ymax=upper3), width=.2, colour="black") +
  scale_fill_manual(values=c("gray", "gray15")) + 
  scale_y_continuous(expand = c(0,0), limits=c(0,4), breaks=c(0,1,2,3,4))  +
  xlab(NULL) + ylab(NULL) + TM +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="white"),
        axis.ticks.length =  unit(.2, "cm"),
        axis.ticks.x = element_blank(),
        plot.margin= unit(c(1, 1, -1, 1), "lines"))

p4 <- ggplot(rpdat) + 
  geom_bar(aes(x=cond, y=resp4, fill=cond), stat = "identity", colour="black") + 
  geom_errorbar(aes(x=cond, ymin=lower4, ymax=upper4), width=.2, colour="black") +
  scale_fill_manual(values=c("gray", "gray15")) + 
  scale_y_continuous(expand = c(0,0), limits=c(0,4), breaks=c(0,1,2,3,4))  +
  xlab(NULL) + ylab(NULL) + TM +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="white"),
        axis.ticks.length =  unit(.2, "cm"),
        axis.ticks.x = element_blank(),
        plot.margin= unit(c(1, 1, -1, 1), "lines"))

grid.arrange(p1,p2,p3,p4, ncol=1) 

# q
rpdat <- data.frame(
  cond = factor(c("yorked-VI","VR30"), levels=c("VR30","yorked-VI")),
  resp1 = c(mean(sim4$res$mean_q_VI),mean(sim4$res$mean_q_VR)),
  resp2 = c(mean(sim3$res$mean_q_VI),mean(sim3$res$mean_q_VR)),
  resp3 = c(mean(sim2$res$mean_q_VI),mean(sim2$res$mean_q_VR)),
  resp4 = c(mean(sim1$res$mean_q_VI),mean(sim1$res$mean_q_VR)),
  lower1 = c(quantile(sim4$res$mean_q_VI, 0.025), quantile(sim4$res$mean_q_VR, 0.025)),
  lower2 = c(quantile(sim3$res$mean_q_VI, 0.025), quantile(sim3$res$mean_q_VR, 0.025)),
  lower3 = c(quantile(sim2$res$mean_q_VI, 0.025), quantile(sim2$res$mean_q_VR, 0.025)),
  lower4 = c(quantile(sim1$res$mean_q_VI, 0.025), quantile(sim1$res$mean_q_VR, 0.025)),
  upper1 = c(quantile(sim4$res$mean_q_VI, 0.975), quantile(sim4$res$mean_q_VR, 0.975)),
  upper2 = c(quantile(sim3$res$mean_q_VI, 0.975), quantile(sim3$res$mean_q_VR, 0.975)),
  upper3 = c(quantile(sim2$res$mean_q_VI, 0.975), quantile(sim2$res$mean_q_VR, 0.975)),
  upper4 = c(quantile(sim1$res$mean_q_VI, 0.975), quantile(sim1$res$mean_q_VR, 0.975))
)

p1 <- ggplot(rpdat) + 
  geom_bar(aes(x=cond, y=resp1, fill=cond), stat = "identity", colour="black") + 
  geom_errorbar(aes(x=cond, ymin=lower1, ymax=upper1), width=.2, colour="black") +
  scale_fill_manual(values=c("gray", "gray15")) + 
  scale_y_continuous(expand = c(0,0), limits=c(0,1), breaks=c(0,.25,.5,.75,1))  +
  xlab(NULL) + ylab(NULL) + TM +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="white"),
        axis.ticks.length =  unit(.2, "cm"),
        axis.ticks.x = element_blank(),
        plot.margin= unit(c(1, 1, -1, 1), "lines"))

p2 <- ggplot(rpdat) + 
  geom_bar(aes(x=cond, y=resp2, fill=cond), stat = "identity", colour="black") + 
  geom_errorbar(aes(x=cond, ymin=lower2, ymax=upper2), width=.2, colour="black") +
  scale_fill_manual(values=c("gray", "gray15")) + 
  scale_y_continuous(expand = c(0,0), limits=c(0,1), breaks=c(0,.25,.5,.75,1))  +
  xlab(NULL) + ylab(NULL) + TM +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="white"),
        axis.ticks.length =  unit(.2, "cm"),
        axis.ticks.x = element_blank(),
        plot.margin= unit(c(1, 1, -1, 1), "lines"))

p3 <- ggplot(rpdat) + 
  geom_bar(aes(x=cond, y=resp3, fill=cond), stat = "identity", colour="black") + 
  geom_errorbar(aes(x=cond, ymin=lower3, ymax=upper3), width=.2, colour="black") +
  scale_fill_manual(values=c("gray", "gray15")) + 
  scale_y_continuous(expand = c(0,0), limits=c(0,1), breaks=c(0,.25,.5,.75,1))  +
  xlab(NULL) + ylab(NULL) + TM +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="white"),
        axis.ticks.length =  unit(.2, "cm"),
        axis.ticks.x = element_blank(),
        plot.margin= unit(c(1, 1, -1, 1), "lines"))

p4 <- ggplot(rpdat) + 
  geom_bar(aes(x=cond, y=resp4, fill=cond), stat = "identity", colour="black") + 
  geom_errorbar(aes(x=cond, ymin=lower4, ymax=upper4), width=.2, colour="black") +
  scale_fill_manual(values=c("gray", "gray15")) + 
  scale_y_continuous(expand = c(0,0), limits=c(0,1), breaks=c(0,.25,.5,.75,1))  +
  xlab(NULL) + ylab(NULL) + TM +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="white"),
        axis.ticks.length =  unit(.2, "cm"),
        axis.ticks.x = element_blank(),
        plot.margin= unit(c(1, 1, -1, 1), "lines"))

grid.arrange(p1,p2,p3,p4, ncol=1) # これがbです
