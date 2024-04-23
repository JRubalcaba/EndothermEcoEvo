######## R Code for the article ######################################################################
## Rubalcaba, J.G. (in review) "The evolution of endothermy as a result of life-history optimization"
## jg.rubalcaba@gmail.com
###################################################################################################### 

######## Temperature Functions ######## 

# temperature dependence of energy assimilation
phi_Tb <- function(Tb, T_A =	10332, T_AL =	50000, T_AH = 90000, T_L = 283.15, T_H = 313.15, T_ref =	293.15){
  phi <- exp(T_A/T_ref - T_A/(Tb+273)) / (1 + exp(T_AL/(Tb+273) - T_AL/T_L) + exp(T_AH/T_H - T_AH/(Tb+273)))
  return(phi)
}
# energy costs of endothermy
theta_Tb <- function(Tb, Ta, r=1){ 
  q <- exp(r*(Tb-Ta))/(1+exp(r*(Tb-Ta))) * (Tb-Ta)
  return(q)
}

######## Analytical solution energy budget : fitness + growth ######## 

balance_function <- function(t0, tmat, tfin, Tb, Ta, mu, lambda, kappa, b, v, a, beta, c){
  phi = phi_Tb(Tb)
  mu_phi = mu * phi
  lambda_phi = lambda * phi
  kappa_theta = kappa * theta_Tb(Tb,Ta)
  
  linf <- (mu_phi - kappa_theta) / lambda_phi
  linf[linf<0]<-0
  k <- lambda_phi*(1-b)
  
  t <- t0:tfin
  T_t <- (1-v) * (t - t0) - v/a * (log(1 + exp(-a*(t-tmat))) - log(1 + exp(-a*(t0-tmat))))
  
  l <- linf * (1 - exp(-k*T_t))
  
  p_t <- v / (1 + exp(-a*(t-tmat)))
  integral <- p_t * l^2 * (linf - l)
  integral_t0 <- integral[1]
  integral_tfin <- integral[tfin]
  
  f <- 3 * beta * k/c * (integral_tfin-integral_t0)
  
  return(cbind(l,f,t))
}

Tb <- seq(0, 50,length.out=100) # body temperature range

mu = 0.05
lambda = 0.001
tmat=365
tfin=365*5

v = 0.9
a = 1
t0 = 0

beta = 1
c = 1
b = 2/3

w0=1
Ta=37
Eassim_Emaint <- mu * w0^0.667 * phi_Tb(Tb) - lambda * w0 * phi_Tb(Tb)

Ta=37
kappa=0.012
Ethermo_37 <- kappa * w0^0.667 * theta_Tb(Tb, Ta, r=4)
Ta=20
kappa=0.011
Ethermo_20_001 <- kappa * w0^0.667 * theta_Tb(Tb, Ta, r=4)
kappa=0.012
Ethermo_20_005 <- kappa * w0^0.667 * theta_Tb(Tb, Ta, r=4)
kappa=0.013
Ethermo_20_010 <- kappa * w0^0.667 * theta_Tb(Tb, Ta, r=4)
kappa=0.014
Ethermo_20_015 <- kappa * w0^0.667 * theta_Tb(Tb, Ta, r=4)
kappa=0.015
Ethermo_20_020 <- kappa * w0^0.667 * theta_Tb(Tb, Ta, r=4)

out_37 <- out_20_0.001 <- out_20_0.005 <- out_20_0.010 <- out_20_0.015 <- out_20_0.020 <- as.data.frame(array(NA, dim=c(100, 3)))
for(i in 1:100){
  Ta=37
  kappa=0.012
  out <- balance_function(t0, tmat, tfin, Tb[i], Ta, mu, lambda, kappa, b, v, a, beta, c)
  out_37[i,] <- out[tfin,]
  
  Ta=20
  kappa=0.011
  out <- balance_function(t0, tmat, tfin, Tb[i], Ta, mu, lambda, kappa, b, v, a, beta, c)
  out_20_0.001[i,] <- out[tfin,]
  kappa=0.012
  out <- balance_function(t0, tmat, tfin, Tb[i], Ta, mu, lambda, kappa, b, v, a, beta, c)
  out_20_0.005[i,] <- out[tfin,]
  kappa=0.013
  out <- balance_function(t0, tmat, tfin, Tb[i], Ta, mu, lambda, kappa, b, v, a, beta, c)
  out_20_0.010[i,] <- out[tfin,]
  kappa=0.014
  out <- balance_function(t0, tmat, tfin, Tb[i], Ta, mu, lambda, kappa, b, v, a, beta, c)
  out_20_0.015[i,] <- out[tfin,]
  kappa=0.015
  out <- balance_function(t0, tmat, tfin, Tb[i], Ta, mu, lambda, kappa, b, v, a, beta, c)
  out_20_0.020[i,] <- out[tfin,]
}

require(ggplot2)
require(RColorBrewer)
df_fig1A <- data.frame(Eassim_Emaint, Ethermo_37, Tb)
p1A <- ggplot() + theme_classic() + ylab("Energy requirements") + xlab("Body temperature") + ylim(-0.001, 0.4) +
  theme(axis.title = element_text(size=14), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12, colour="black")) +
  geom_line(df_fig1A, mapping=aes(y=Eassim_Emaint, x=Tb), size=0.8, col="black") +
  geom_line(df_fig1A, mapping=aes(y=Ethermo_37, x=Tb), size=0.8, col="#D6604D")
#
df_fig1B <- data.frame(out_37, Tb)
p1B <- ggplot() + theme_classic() + ylab("Scope for reproduction") + xlab("Body temperature") + 
  theme(axis.title = element_text(size=14), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12, colour="black")) +
  geom_line(df_fig1B, mapping=aes(y=V2, x=Tb), size=0.8, col="black")
#
df_fig1C <- data.frame(Eassim_Emaint=rep(Eassim_Emaint,5),
                       Ethermo=c(Ethermo_20_001,Ethermo_20_005,Ethermo_20_010,Ethermo_20_015,Ethermo_20_020),
                       k_value=c(sort(rep(c("1","2","3","4","5"),100))),
                       Tb=rep(Tb,5))
p1C <- ggplot() + theme_classic() + ylab("Energy requirements") + xlab("Body temperature") + ylim(-0.001, 0.4) +
  theme(legend.position = "none",
        axis.title = element_text(size=14), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12, colour="black")) +
  geom_line(df_fig1C, mapping=aes(y=Eassim_Emaint, x=Tb), size=0.8, col="black") +
  geom_line(df_fig1C, mapping=aes(y=Ethermo, x=Tb, col=k_value), size=0.8) +
  scale_colour_manual(values=rev(brewer.pal(6, "RdBu")))
#
df_fig1D <- data.frame(Wfin=c(out_20_0.001$V2,out_20_0.005$V2,out_20_0.010$V2,out_20_0.015$V2,out_20_0.020$V2),
                       k_value=c(sort(rep(c("1","2","3","4","5"),100))),
                       Tb=rep(Tb,5))
p1D <- ggplot() + theme_classic() + ylab("Scope for reproduction") + xlab("Body temperature") + 
  theme(legend.position = "none",
        axis.title = element_text(size=14), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12, colour="black")) +
  geom_line(df_fig1D, mapping=aes(y=(Wfin), x=Tb, col=k_value), size=0.8) +
  scale_colour_manual(values=rev(brewer.pal(6, "RdBu")))

######## Fig. 2 ----

require(ggpubr)
ggarrange(p1A, p1B, p1C, p1D, ncol=2, nrow=2) # 800x600


######## Evolutionary dynamics ######## 

EvolDyn <-function(params, time){
  # calculations t=0
  W <- numeric(time)
  N <- numeric(time)
  Tb_i <- numeric(time)
  f_mean <- numeric(time)
  l_mean <- numeric(time)
  vA_Tb_N <- numeric(time)
  vP_Tb_N <- numeric(time)
  Tenv <- numeric(time)
  S_coef <- numeric(time)
  
  # genetics
  Tb0 = params$Tb
  vP_Tb_N[1] = params$vP_Tb
  h2 = params$h2
  vmut = params$vmut
  
  # population parameters
  Wmax = params$Wmax
  N0 = params$N0
  K0 = params$K0
  NC = params$NC
  
  vA_Tb_N[1] <- vP_Tb_N[1] * h2 # initial additive genetic variance
  vE_Tb <- vP_Tb_N[1] - vA_Tb_N[1] # initial environmental variance
  
  N[1] <- N0 # Initial population size
  Tb_i[1] <- Tb0
  
  Tb <- rnorm(N[1], Tb_i[1], sqrt(vP_Tb_N[1])) # Population values for trait, with N[1] individuals...
  
  # Ta <- seq(params$Ta_max, params$Ta_min, length.out=time)
  Ta <- params$Ta_min + (params$Ta_max-params$Ta_min) * exp(0.1*(time/4 - 1:time)) / (1 + exp(0.1*(time/4 - 1:time)))
  
  # Ta <- params$Ta_min + (params$Ta_max-params$Ta_min)/2 + (params$Ta_max-params$Ta_min)/2 * sin(2+0.01 * 1:time)
  
  # Run model
  for(i in 2:time){
    
    if(N[i-1] < NC){
      break
    }
    
    # Covariance Tb vs Fitness
    Tb_seq <- seq(min(Tb,na.rm=T), max(Tb,na.rm=T), length.out=20)
    out_tb <- array(NA, dim=c(20, 3))
    for(j in 1:20){
      out <- balance_function(t0=0, tmat=params$tmat, tfin=params$tfin, Tb=Tb_seq[j], 
                              Ta=Ta[i], mu=params$mu, lambda=params$lambda, kappa=params$kappa, 
                              b=params$b, v=params$v, a=1, beta=1, c=1)
      out_tb[j,] <- out[params$tfin,]
    }
    
    f_mean[i] <- mean(out_tb[,2], na.rm=T)
    l_mean[i] <- mean(out_tb[,1], na.rm=T)
    
    #### Genetic variance components
    h2 <- params$h2 # vA_Tb_N[i-1] / vP_Tb_N[i-1]
    vA_Tb_N[i] <- vA_Tb_N[i-1] * (1-(1/(2*N[i-1]))) #loss of va due to drift https://lukejharmon.github.io/pcm/chapter3_bmintro/ 
    vE_Tb <- vP_Tb_N[i-1] - vA_Tb_N[i]
    vM_Tb <- vmut * vE_Tb #vA_Tb_N[i]# mutational variance ~1e-3 1e-2 VE, Johnson and Barton 2005 Philos trans
    
    vP_Tb_N[i] <- vA_Tb_N[i] + vE_Tb + vM_Tb  
    
    ### Selection coefficient
    S <- cov(out_tb[,2], Tb_seq, use="complete") / f_mean[i]
    
    ### Breeder equation
    R <- h2 * S
    
    # New population numbers
    newTb <- mean(Tb,na.rm=T) + R
    Tb <- rnorm(N[i-1], newTb, sqrt(vP_Tb_N[i]))
    Tb_i[i] <- newTb
    
    #population growth
    K <- round(rnorm(1, K0, K0*0.000010))
    
    W <- f_mean[i] * params$r
    N[i-1] <- ifelse(N[i-1] > K, K, N[i-1])
    N[i] <- round(N[i-1] + ((log(W)*N[i-1]) * (1 - (N[i-1]/K))))
    
    print(i)
  }
  
  out <- data.frame(Tb_i, vA_Tb_N, vP_Tb_N, N, l_mean, f_mean, Ta, S)
  return(out[-1,])
} 

kappa <- 0.02 
params <- list(Ta_max=37, Ta_min=25, mu=0.05, lambda=0.001, kappa=kappa, v=0.9, b=2/3, tmat=365/2, tfin=365*1,
               Tb=37, vP_Tb=2, h2=0.5,                                    
               r=40, N0=100, K0=1e3, NC=10, vmut=1e-3) 
out <- EvolDyn(params, time=300)

Eassim <- params$mu*out$l_mean^params$b*phi_Tb(out$Tb_i)
Emainten <- params$lambda*out$l_mean*phi_Tb(out$Tb_i)
Ethermo <- params$kappa*out$l_mean^params$b*theta_Tb(out$Tb_i,out$Ta,r=10)
Ta <- params$Ta_min + (params$Ta_max-params$Ta_min) * exp(0.1*(299/4 - 1:299)) / (1 + exp(0.1*(299/4 - 1:299)))

Ethermo_simul <- Tb_simul <-l_simul <- array(NA, dim=c(299,10))
for(i in 1:10){
  out <- EvolDyn(params, time=300)
  Tb_simul[,i] <- out$Tb_i
  l_simul[,i] <- out$l_mean
  Ethermo_simul[,i] <- params$kappa*out$l_mean^params$b*theta_Tb(out$Tb_i,out$Ta,r=10)
  points(out$Tb_i)
}

plot(Tb_simul[,1], ylim=c(10,40))
for(i in 2:10) points(Tb_simul[,i])
plot(Ethermo_simul[,1], ylim=c(0,1))
for(i in 2:10) points(Ethermo_simul[,i])
plot(l_simul[,1])
for(i in 2:10) points(l_simul[,i])
plot(Tb_simul[,1] ~ l_simul[,1])
for(i in 2:10) points(Tb_simul[,i] ~ l_simul[,i])

df_fig3 <- data.frame(Tb=as.numeric(Tb_simul), Ethermo=as.numeric(Ethermo_simul), l=as.numeric(l_simul),
                      simul=sort(rep(1:10,299)),
                      Ta=rep(Ta,10), time=rep(1:299,10))

p3A <- ggplot() + theme_classic() + ylab("Temperature") + xlab("time") + ylim(15,40) +
  theme(axis.title = element_text(size=14), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12, colour="black")) +
  geom_line(df_fig3, mapping=aes(y=Tb, x=time, group=simul), size=1, col="black") +
  geom_line(df_fig3, mapping=aes(y=Ta, x=time), size=1, col="#4393C3") 

p3B <- ggplot() + theme_classic() + ylab("Thermoregulatory costs") + xlab("time") + ylim(-0.001, 1) +
  theme(axis.title = element_text(size=14), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12, colour="black")) +
  geom_line(df_fig3, mapping=aes(y=Ethermo, x=time, group=simul), size=1, col="black") 

######## Fig. 3 ----

ggarrange(p3A, p3B, ncol=1, nrow=2) 



