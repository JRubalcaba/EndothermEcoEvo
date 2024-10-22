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
# energy costs of thermoregulation
theta_Tb <- function(Tb, Ta, hmin=1, hmax=4, B=12, gamma=1){
  if(length(Tb)==1 & length(Ta)==1){
    Thmin <- B/hmin+Ta
    Thmax <- B/hmax+Ta
    
    h <- B/(Tb-Ta)
    h[Tb>=Thmin] <- hmin
    h[Tb<Thmax] <- hmax
    
    Qtot <- Qt <- h * (Tb-Ta) - B
    Qt[Tb>Thmax] <- Qtot[Tb>Thmax] 
    Qt[Tb<Thmin] <- -gamma*Qtot[Tb<Thmin] 
  }
  if(length(Ta)>1 & length(Tb)==1){
    Thmin <- Tb-B/hmin
    Thmax <- Tb-B/hmax

    h <- B/(Tb-Ta)
    h[Ta<=Thmin] <- hmin
    h[Ta>Thmax] <- hmax
    
    Qtot <- Qt <- h * (Tb-Ta) - B
    Qt[Ta<Thmin] <- Qtot[Ta<Thmin] 
    Qt[Ta>Thmax] <- -gamma*Qtot[Ta>Thmax] 
  }
  if(length(Tb)>1 & length(Ta)==1){
    Thmin <- B/hmin+Ta
    Thmax <- B/hmax+Ta
    
    h <- B/(Tb-Ta)
    h[Tb>=Thmin] <- hmin
    h[Tb<Thmax] <- hmax

    Qtot <- Qt <- h * (Tb-Ta) - B
    Qt[Tb>Thmax] <- Qtot[Tb>Thmax] 
    Qt[Tb<Thmin] <- -gamma*Qtot[Tb<Thmin] 
  }
  return(Qt)
}

## Details on thermoregulation costs model
Tb=37   # Body temperature
Ta=10:40# Air temperature
hmin=1  # Min and Max heat transfer coefficients: determine the range within 
hmax=4  #       which the animal can thermoregulate by adjusting the thermal conductivity 
        #       without incurring energy costs
B=12    # Basal metabolic heat production   
gamma=1 # Cooling rate-energy conversion coefficient

Qmin <- hmin*(Tb-Ta) # Min and max heat dissipation rates given 
Qmax <- hmax*(Tb-Ta) # min and max heat transfer coefficients

require(ggplot2)
ggplot() + theme_classic() + ylab("Heat dissipation") + xlab("Air temperature") +
  ylim(-5,40) +
  geom_line(mapping=aes(y=Qmin, x=Ta)) + geom_line(mapping=aes(y=Qmax, x=Ta)) +
  geom_hline(yintercept = B, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = Tb, linetype="dashed")

Thmin <- Tb-B/hmin # air temperatures at which hmin / hmax are reached
Thmax <- Tb-B/hmax
 
h <- B/(Tb-Ta) # heat transfer coefficient required to maintain heat production at basal
h[Ta<=Thmin] <- hmin # constrain the required h within hmin and hmax
h[Ta>Thmax] <- hmax 

p1 <- ggplot() + theme_classic() + ylab("Heat transfer coefficient") + xlab("Air temperature") +
  ylim(0,5) + xlim(10,40) +
  geom_line(mapping=aes(y=h, x=Ta), size=1.2) +
  geom_hline(yintercept = hmin, linetype="dashed") +
  geom_hline(yintercept = hmax, linetype="dashed")

Qtot <- Qt <- h * (Tb-Ta) - B 
Qt[Ta<Thmin] <- Qtot[Ta<Thmin] # Energy cost of heating 
Qt[Ta>Thmax] <- -gamma*Qtot[Ta>Thmax] # Energy cost of cooling 

p2 <- ggplot() + theme_classic() + ylab("Energy costs of thermoregulation") + xlab("Air temperature") +
  ylim(-5,40) +
  geom_line(mapping=aes(y=Qmin, x=Ta), lty=3) + 
  geom_line(mapping=aes(y=Qmax, x=Ta), lty=3) +
  geom_line(mapping=aes(y=Qt+B, x=Ta), size=1.2) +
  geom_hline(yintercept = B, linetype="dashed") +
  geom_vline(xintercept = Tb, linetype="dashed")

require(ggpubr)
ggarrange(p1,p2)

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
Ta1=25
Ta2=10

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
Eassim_Emaint <- mu * w0^0.667 * phi_Tb(Tb) - lambda * w0 * phi_Tb(Tb)

kappa=0.013
Ethermo_37 <- kappa * w0^0.667 * theta_Tb(Tb, Ta1)

kappa=0.013
Ethermo_20_001 <- kappa * w0^0.667 * theta_Tb(Tb, Ta2)
kappa=0.014
Ethermo_20_005 <- kappa * w0^0.667 * theta_Tb(Tb, Ta2)
kappa=0.015
Ethermo_20_010 <- kappa * w0^0.667 * theta_Tb(Tb, Ta2)
kappa=0.016
Ethermo_20_015 <- kappa * w0^0.667 * theta_Tb(Tb, Ta2)
kappa=0.017
Ethermo_20_020 <- kappa * w0^0.667 * theta_Tb(Tb, Ta2)

out_37 <- out_20_0.001 <- out_20_0.005 <- out_20_0.010 <- out_20_0.015 <- out_20_0.020 <- as.data.frame(array(NA, dim=c(100, 3)))
for(i in 1:100){
  kappa=0.013
  out <- balance_function(t0, tmat, tfin, Tb[i], Ta1, mu, lambda, kappa, b, v, a, beta, c)
  out_37[i,] <- out[tfin,]
  
  kappa=0.013
  out <- balance_function(t0, tmat, tfin, Tb[i], Ta2, mu, lambda, kappa, b, v, a, beta, c)
  out_20_0.001[i,] <- out[tfin,]
  kappa=0.014
  out <- balance_function(t0, tmat, tfin, Tb[i], Ta2, mu, lambda, kappa, b, v, a, beta, c)
  out_20_0.005[i,] <- out[tfin,]
  kappa=0.015
  out <- balance_function(t0, tmat, tfin, Tb[i], Ta2, mu, lambda, kappa, b, v, a, beta, c)
  out_20_0.010[i,] <- out[tfin,]
  kappa=0.016
  out <- balance_function(t0, tmat, tfin, Tb[i], Ta2, mu, lambda, kappa, b, v, a, beta, c)
  out_20_0.015[i,] <- out[tfin,]
  kappa=0.017
  out <- balance_function(t0, tmat, tfin, Tb[i], Ta2, mu, lambda, kappa, b, v, a, beta, c)
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
  geom_line(df_fig1A, mapping=aes(y=Eassim_Emaint, x=Tb), linewidth=0.8, col="black") +
  geom_line(df_fig1A, mapping=aes(y=Ethermo_37, x=Tb), linewidth=0.8, col="#D6604D")
#
df_fig1B <- data.frame(out_37, Tb)
p1B <- ggplot() + theme_classic() + ylab("Scope for reproduction") + xlab("Body temperature") + 
  theme(axis.title = element_text(size=14), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=12, colour="black")) +
  geom_line(df_fig1B, mapping=aes(y=V2, x=Tb), linewidth=0.8, col="black")
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
  
  # Ta_mean <- seq(params$Ta_max, params$Ta_min, length.out=time)
  # Ta_mean <- params$Ta_min + (params$Ta_max-params$Ta_min)/2 + (params$Ta_max-params$Ta_min)/2 * sin(2+0.01 * 1:time)
  Ta_mean <- params$Ta_min + (params$Ta_max-params$Ta_min) * exp(0.05*(time/3 - 1:time)) / (1 + exp(0.05*(time/3 - 1:time)))
  Ta <- Ta_mean + rnorm(time, 0, params$Ta_SD)
  
  # plot(Ta)
  
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
    
    if(f_mean[i] == 0){
      break
    }
    
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
  out$Tb_i[which(out$Tb_i == 0)] <- NA
  return(out[-1,])
} 

kappa <- 0.005 # 0.005, 0.03; branching: 0.0183 (SD=1, vP_Tb=4)
params <- list(Ta_max=25, Ta_min=10, Ta_SD=2,
               mu=0.05, lambda=0.001, kappa=kappa, v=0.9, b=2/3, tmat=365/2, tfin=365*1,
               Tb=37, vP_Tb=4, h2=0.5,                                    
               r=50, N0=500, K0=1e3, NC=100, vmut=1e-3) 
out <- EvolDyn(params, time=500)

Eassim <- params$mu*out$l_mean^params$b*phi_Tb(out$Tb_i)
Emainten <- params$lambda*out$l_mean*phi_Tb(out$Tb_i)
theta <- numeric(499)
for(i in 1:499) theta[i]<-theta_Tb(out$Tb_i[i],out$Ta[i])
Ethermo <- params$kappa*out$l_mean^params$b*theta
Ta <- params$Ta_min + (params$Ta_max-params$Ta_min) * exp(0.05*(499/3 - 1:499)) / (1 + exp(0.05*(499/3 - 1:499)))

plot(out$Tb_i, pch=20, col="blue", ylim=c(20,40))

reps=20
Ethermo_simul <- Tb_simul <- Ta_simul <- l_simul <- array(NA, dim=c(499,reps))
for(i in 1:reps){
  out <- EvolDyn(params, time=500)
  Tb_simul[,i] <- out$Tb_i
  l_simul[,i] <- out$l_mean
  Ta_simul[,i] <- out$Ta
  theta <- numeric(499)
  for(j in 1:499) theta[j]<-theta_Tb(out$Tb_i[j],out$Ta[j])
  Ethermo_simul[,i] <- params$kappa*out$l_mean^params$b*theta
  points(out$Tb_i)
}

plot(Tb_simul[,1], ylim=c(10,40))
for(i in 2:reps) points(Tb_simul[,i])
plot(Tb_simul[which(Ta_simul[,1]<15),1] ~ Ta_simul[which(Ta_simul[,1]<15),1], ylim=c(10,40))
for(i in 2:reps) points(Tb_simul[which(Ta_simul[,i]<15),i] ~ Ta_simul[which(Ta_simul[,i]<15),i])
plot(Ethermo_simul[,1], ylim=c(0,1))
for(i in 2:reps) points(Ethermo_simul[,i])
plot(l_simul[,1])
for(i in 2:reps) points(l_simul[,i])
plot(Tb_simul[,1] ~ l_simul[,1])
for(i in 2:reps) points(Tb_simul[,i] ~ l_simul[,i])

df_fig3 <- data.frame(Tb=as.numeric(Tb_simul), Ethermo=as.numeric(Ethermo_simul), l=as.numeric(l_simul),
                      simul=sort(rep(1:reps,499)),
                      Ta=rep(Ta,10), time=rep(1:499,reps))

p3A <- ggplot() + theme_classic() + ylab("Temperature") + xlab("time") + ylim(20,40) +
  theme(axis.title = element_text(size=14), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12, colour="black")) +
  geom_line(df_fig3, mapping=aes(y=Tb, x=time, group=simul), size=1, col="black") 
  # geom_line(df_fig3, mapping=aes(y=Ta, x=time), size=1, col="#4393C3") 

p3B <- ggplot() + theme_classic() + ylab("Thermoregulatory costs") + xlab("time") + ylim(-0.001, 0.5) +
  theme(axis.title = element_text(size=14), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12, colour="black")) +
  geom_line(df_fig3, mapping=aes(y=Ethermo, x=time, group=simul), size=1, col="black") 

######## Fig. 3 ----

ggarrange(p3A, p3B, ncol=1, nrow=2) #600x400

##### Temperature fluctuations, heterotermy vs homeothermy

# FIRST set gamma = 0.2 in the theta function

Ta_warm <- 25
Ta_cold <- 10
SD <- 1
kappa_homeo <- 0.005
kappa_hetero <- 0.03

params_warm_homeo <- list(Ta_max=Ta_warm, Ta_min=Ta_warm, Ta_SD=SD,
                          mu=0.05, lambda=0.001, kappa=kappa_homeo, v=0.9, b=2/3, tmat=365/2, tfin=365*1,
                          Tb=37, vP_Tb=4, h2=0.5,                                    
                          r=50, N0=500, K0=1e3, NC=100, vmut=1e-3) 
out_warm_homeo <- EvolDyn(params_warm_homeo, time=2000)

params_warm_hetero <- list(Ta_max=Ta_warm, Ta_min=Ta_warm, Ta_SD=SD,
                          mu=0.05, lambda=0.001, kappa=kappa_hetero, v=0.9, b=2/3, tmat=365/2, tfin=365*1,
                          Tb=37, vP_Tb=4, h2=0.5,                                    
                          r=50, N0=500, K0=1e3, NC=100, vmut=1e-3) 
out_warm_hetero <- EvolDyn(params_warm_hetero, time=2000)

params_cold_homeo <- list(Ta_max=Ta_cold, Ta_min=Ta_cold, Ta_SD=SD,
                          mu=0.05, lambda=0.001, kappa=kappa_homeo, v=0.9, b=2/3, tmat=365/2, tfin=365*1,
                          Tb=37, vP_Tb=4, h2=0.5,                                    
                          r=50, N0=500, K0=1e3, NC=100, vmut=1e-3) 
out_cold_homeo <- EvolDyn(params_cold_homeo, time=2000)

params_cold_hetero <- list(Ta_max=Ta_cold, Ta_min=Ta_cold, Ta_SD=SD,
                           mu=0.05, lambda=0.001, kappa=kappa_hetero, v=0.9, b=2/3, tmat=365/2, tfin=365*1,
                           Tb=37, vP_Tb=4, h2=0.5,                                    
                           r=50, N0=500, K0=1e3, NC=100, vmut=1e-3) 
out_cold_hetero <- EvolDyn(params_cold_hetero, time=2000)

require(ggplot2)
ggplot() + theme_classic() + ylab("Body temperature (ºC)") + xlab("Ambient temperature (ºC)") +
  geom_point(mapping=aes(y=out_warm_homeo$Tb_i, x=out_warm_homeo$Ta), col="red") +
  geom_smooth(mapping=aes(y=out_warm_homeo$Tb_i, x=out_warm_homeo$Ta), method="lm", col="darkred") +
  geom_point(mapping=aes(y=out_cold_homeo$Tb_i, x=out_cold_homeo$Ta), col="blue") +
  geom_smooth(mapping=aes(y=out_cold_homeo$Tb_i, x=out_cold_homeo$Ta), method="lm", col="darkblue") 

ggplot() + theme_classic() + ylab("Body temperature (ºC)") + xlab("Ambient temperature (ºC)") +
  geom_point(mapping=aes(y=out_warm_hetero$Tb_i, x=out_warm_hetero$Ta), col="red") +
  geom_smooth(mapping=aes(y=out_warm_hetero$Tb_i, x=out_warm_hetero$Ta), method="lm", col="darkred") +
  geom_point(mapping=aes(y=out_cold_hetero$Tb_i, x=out_cold_hetero$Ta), col="blue") +
  geom_smooth(mapping=aes(y=out_cold_hetero$Tb_i, x=out_cold_hetero$Ta), method="lm", col="darkblue") 

summary(lm(out_warm_homeo$Tb_i ~ out_warm_homeo$Ta))
summary(lm(out_warm_hetero$Tb_i ~ out_warm_hetero$Ta))

summary(lm(out_cold_homeo$Tb_i ~ out_cold_homeo$Ta))
summary(lm(out_cold_hetero$Tb_i ~ out_cold_hetero$Ta))

### Vector field Tb vs Ta 

n_vectors = 30 # number of vectors
size = 0.05 # vector size

Tb_seq <- seq(20, 42, length.out=n_vectors)
Ta_seq <- seq(7, 30, length.out=n_vectors)

kappa <- 0.03 # 0.0215 (branching) 
params <- list(Ta_max=25, Ta_min=10, Ta_SD=1,
               mu=0.05, lambda=0.001, kappa=kappa, v=0.9, b=2/3, tmat=365/2, tfin=365*1,
               Tb=37, vP_Tb=4, h2=0.5,
               r=50, N0=500, K0=1e3, NC=100, vmut=1e-3)

dW_dTa <- dW_dTb <- array(NA, dim=c(n_vectors,n_vectors))
for(i in 1:n_vectors){
  Ta <- Ta_seq[i]
  for(j in 1:n_vectors){
    Tb <- Tb_seq[j]  
    Tb_diff <- seq(Tb-size, Tb+size, length.out=20)
    out_tb <- array(NA, dim=c(20, 3))
    for(k in 1:20){
      out <- balance_function(t0=0, tmat=params$tmat, tfin=params$tfin, Tb=Tb_diff[k], 
                              Ta=Ta, mu=params$mu, lambda=params$lambda, kappa=params$kappa, 
                              b=params$b, v=params$v, a=1, beta=1, c=1)
      out_tb[k,] <- out[params$tfin,]
    }
    dW_dTb[i,j] <- cov(out_tb[,2],Tb_diff)#/mean(out_tb[,2])
    
    # Ta_diff <- seq(Ta-0.1, Ta+0.1, length.out=20)
    # out_ta <- array(NA, dim=c(20, 3))
    # for(k in 1:20){
    #   out <- balance_function(t0=0, tmat=params$tmat, tfin=params$tfin, Tb=Tb, 
    #                           Ta=Ta_diff[k], mu=params$mu, lambda=params$lambda, kappa=params$kappa, 
    #                           b=params$b, v=params$v, a=1, beta=1, c=1)
    #   out_ta[k,] <- out[tfin,]
    # }
    # dW_dTa[i,j] <- cov(out_ta[,2],Tb_diff)/mean(out_ta[,2])
  }
}

df_vector <- data.frame(Ta=rep(Ta_seq,n_vectors), Tb=sort(rep(Tb_seq,n_vectors)), 
                        dW_dTa=as.numeric(dW_dTa), dW_dTb=as.numeric(dW_dTb),
                        module=sqrt(as.numeric(dW_dTa)^2+as.numeric(dW_dTb)^2))
df_vector$sign <- df_vector$dW_dTb
df_vector$sign[df_vector$sign<0] <- -1
df_vector$sign[df_vector$sign>0] <- 1

require(grid)
require(ggplot2)
# First load "functions multiple scales.R"

ggplot(df_vector, aes(y = Tb, x = Ta)) + theme_bw() + # 500x400
  theme(panel.grid = element_blank(), legend.position="none",
        axis.text = element_text(colour="black")) +
  ylab("Body temperature") + xlab("Ambient temperature") +
  xlim(range(Ta_seq)) + ylim(range(Tb_seq)) +
  geom_segment(df_vector[df_vector$sign>0,], mapping=aes(yend = Tb+dW_dTb*15, xend = Ta, 
                                                         colour=(dW_dTb)),
               arrow = arrow(length = unit(0.3, "cm"), type="closed") ) +
  scale_colour_continuous(low="pink", high="darkred") +
  new_scale_color() +
  geom_segment(df_vector[df_vector$sign<0,], mapping=aes(yend = Tb+dW_dTb*5, xend = Ta, 
                                                         colour=-(dW_dTb)),
               arrow = arrow(length = unit(0.3, "cm"), type="closed")) +
  scale_colour_continuous(low="lightblue", high="darkblue") 

###### Body size dynamics

kappa <- 0.03 # branching: 0.0189

Ta <- seq(25, 10, length.out=200)
Ethermo_simul <- Tb_simul <- l_simul <- numeric(length(Ta))
for(i in 1:length(Ta)){
  params <- list(Ta_max=Ta[i], Ta_min=Ta[i], Ta_SD = 0,
                 mu=0.05, lambda=0.001, kappa=kappa, v=0.9, b=2/3, tmat=365/2, tfin=365*1,
                 Tb=37,  vP_Tb=4, h2=0.5,                                    
                 r=50, N0=500, K0=1e3, NC=100, vmut=1e-3) 

  out <- EvolDyn(params, time=50)
  Tb_simul[i] <- out$Tb_i[nrow(out)]
  l_simul[i] <- out$l_mean[nrow(out)]
  Ethermo_simul[i] <- params$kappa*out$l_mean[nrow(out)]^params$b*theta_Tb(out$Tb_i[nrow(out)],Ta[i])
}

ggplot() + theme_classic() + ylab("Body tempeature") + xlab("Body size") + ylim(20,40) + xlim(0,15) +
  theme(axis.title = element_text(size=14), 
        axis.text = element_text(size=12, colour="black")) +
  geom_point(mapping=aes(y=Tb_simul, x=l_simul, colour=Ta), size=3) +
  scale_color_distiller(palette="RdBu") # 500x400

ggplot() + theme_classic() + ylab("Body tempeature") + xlab("Ambient temperature") + ylim(25,40) + xlim(18,37) +
  theme(axis.title = element_text(size=14), 
        axis.text = element_text(size=12, colour="black")) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(mapping=aes(y=Tb_simul, x=Ta, size=round(l_simul)^2)) 
