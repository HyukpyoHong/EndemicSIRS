ODE_SIRS_prop_vacc_v3 <- function(Time, State, Pars){
  # This function simulate the time course of the mathematical model (see Fig. 1C and supplementary methods) for given parameters:
  # Ntot: total number of the population
  # beta.param: infectious rate
  # gamma.param: recovery rate
  # omega_RL: waining rate of the infection-preventing immunity (R -> S_L)
  # omega_LH: waining rate of the severity-reducing immunity (S_L -> S_H)
  # v: daily vaccination rate among susceptible populations.
  # h_S: probability of experiencing severe disases when an individual from S_H get infected.
  # l_S: probability of experiencing severe disases when an individual from S_L get infected. 
  # Because of the severity-reducing immunity of S_L population, l_S < h_S.
  
  with(as.list(c(State, Pars)), {
    dS_H = -v * S_H - beta.param * S_H * (I_M + I_S) / Ntot + omega_LH * S_L
    dS_L = -v * S_L - beta.param * S_L * (I_M + I_S) / Ntot - omega_LH * S_L + omega_RL * R
    dI_M = beta.param / Ntot * (I_M + I_S) * ((1-h_S)*S_H + (1-l_S)*S_L) - gamma.param * I_M
    dI_S = beta.param / Ntot * (I_M + I_S) * (h_S*S_H + l_S*S_L) - gamma.param * I_S
    dR = v * (S_H + S_L) + gamma.param * (I_M + I_S) - omega_RL * R
    return(list(c(dS_H, dS_L, dI_M, dI_S, dR)))
  })
}

SS_SIRS_prop_vacc_v3 <- function(Ntot, beta.param, gamma.param, omega_RL, omega_LH, v, h_S, l_S){
  # This function calculate the steady-state values of the mathematical model (see Fig. 1C and supplementary methods) for given parameters:
  # Ntot: total number of the population
  # beta.param: infectious rate
  # gamma.param: recovery rate
  # omega_RL: waining rate of the infection-preventing immunity (R -> S_L)
  # omega_LH: waining rate of the severity-reducing immunity (S_L -> S_H)
  # v: daily vaccination rate among susceptible populations.
  # h_S: probability of experiencing severe disases when an individual from S_H get infected.
  # l_S: probability of experiencing severe disases when an individual from S_L get infected. 
  # Because of the severity-reducing immunity of S_L population, l_S < h_S.
  
  if(gamma.param/beta.param < omega_RL/(omega_RL+v)){
    Stot_inf = gamma.param/beta.param * Ntot
    Itot_inf = 1/(beta.param*(omega_RL + gamma.param)) * (omega_RL*beta.param - omega_RL*gamma.param - v * gamma.param) * Ntot
    R_inf = Ntot - Stot_inf - Itot_inf
    S_L_inf = omega_RL * R_inf/(v + omega_LH + beta.param * Itot_inf / Ntot)
    S_H_inf = Stot_inf - S_L_inf
    I_S_inf = beta.param * Itot_inf / (Ntot *gamma.param) *((h_S*S_H_inf + l_S*S_L_inf))
    I_M_inf = beta.param * Itot_inf / (Ntot *gamma.param) *(((1-h_S)*S_H_inf + (1-l_S)*S_L_inf))
  }else{
    Stot_inf = gamma.param/beta.param * Ntot
    Itot_inf = 0
    R_inf = Ntot - Stot_inf - Itot_inf
    S_L_inf = omega_RL * R_inf/(v + omega_LH)
    S_H_inf = Stot_inf - S_L_inf
    I_S_inf = 0
    I_M_inf = 0
  }
  return(list(S_H_inf = S_H_inf, S_L_inf=S_L_inf, I_S_inf=I_S_inf, I_M_inf=I_M_inf, R_inf=R_inf))
}


# install.packages("deSolve") # Before the first execution of the code, the package for numerical ODE solver, 'deSolve', has to be installed once.
# library("deSolve") # When R session is newly open, deSolve package should be called. 

####### ===================== FIGURE production ======================== ######

#### FIG 2 (R0 and original parameter) and S2 (R0 and altered gamma, w1, and w2) ===========

S_H_init = 5*10^7
S_L_init = 0
I_M_init = 4750
I_S_init = 250
R_init = 0
yinit = c(S_H = S_H_init, S_L = S_L_init, I_M = I_M_init, I_S = I_S_init, R = R_init)

# for FIG 2
omega_RLtmp = 1/365 * 1
omega_LHtmp = 1/365 * 1/3 
gam = 1/10

betalist = seq(from = gam, by = 0.01, to = 16*gam)
vacc.rate.list0 = c(0, 0.0001, 0.001, 0.01)


# for FIG S2
# omega_RLtmp = 1/365 * 1/1.5
# omega_LHtmp = 1/365 * 1/6 
# gam = 1/4

par(mfrow = c(3, 4))
par(mar = c(2,2,2,2))
for(ii in 1:3){
  for(imm_efficacy in c(0.8,0.9,0.95,1)){
    
    ISlist = matrix(NA, nrow = length(betalist), ncol = length(vacc.rate.list0))
    h_S = 0.05
    l_S = h_S * (1-imm_efficacy)
    for(bbid in 1:length(betalist)){
      for(vacc.id in 1:length(vacc.rate.list0)){
        vacc = vacc.rate.list0[vacc.id]
        bb = betalist[bbid]
        out1 = SS_SIRS_prop_vacc_v3(sum(yinit), beta.param = bb, gamma.param = gam, omega_RL = omega_RLtmp, omega_LH = omega_LHtmp, v = vacc,
                                    h_S = h_S, l_S = l_S)
        if(ii == 1){
          # 1. Daily infection cases
          ylim.max = 0.25
          ISlist[bbid, vacc.id] = gam * (out1$I_M_inf + out1$I_S_inf) / sum(yinit) * 100
          
        }else if(ii == 2){
          # 2. ratio for progression of severe disease, i.e., I_S / (I_M + I_S)
          ylim.max = 4
          if((out1$I_M_inf + out1$I_S_inf) == 0){
            ISlist[bbid, vacc.id] = NA
          }else{
            ISlist[bbid, vacc.id] = out1$I_S_inf/(out1$I_M_inf + out1$I_S_inf) * 100
          }
        }else if(ii == 3){
          # 1. Daily severe cases
          ylim.max = 30*10^(-4)
          ISlist[bbid, vacc.id] = gam * out1$I_S_inf / sum(yinit) * 100
        }else{
          if (out1$S2_inf == 0){
            ISlist[bbid, vacc.id] = NA
          }else{
            ISlist[bbid, vacc.id] = out1$S_H_inf / out1$S_L_inf * 100
          }
        }
      }
    }
    
    plot(betalist / gam, ISlist[,1], type = "l", col = 1, 
         ylim = c(0, ylim.max), ylab = "", 
         # xlim = c(0, max(betalist / gam)), lwd = 2, cex.lab = 2)
         xlim = c(1, 8), lwd = 1, cex.lab = 2)#, log = "x" # , xaxt="n": remove xtick
    # axis(side=1, at=c(1,2,4,8,16))
    collist = matrix(c(0,0,0,
                       204,121,167,
                       0,114,175,
                       82,167,221,
                       0,158,115,
                       213,94,0)/255, ncol = 3, byrow=TRUE)
    for(colid in 2:length(vacc.rate.list0)){
      lines(betalist / gam, ISlist[,colid], col = rgb(collist[colid,1],collist[colid,2],collist[colid,3]), lwd = 1)
    }
  }
}

#### FIG S1 (Rv and original parameter) ===========

# for FIG S1
omega_RLtmp = 1/365 * 1
omega_LHtmp = 1/365 * 1/3 
gam = 1/10

par(mfrow = c(3, 4), mar = c(2,2,2,2), oma = c(2, 2, 2, 0))

for(ii in 1:3){
  for(imm_efficacy in c(0.8,0.9,0.95,1)){
    betalist = seq(from = gam, by = 0.01, to = 60*gam)
    ISlist = matrix(NA, nrow = length(betalist), ncol = length(vacc.rate.list0))
    h_S = 0.05
    l_S = h_S * (1-imm_efficacy)
    for(bbid in 1:length(betalist)){
      for(vacc.id in 1:length(vacc.rate.list0)){
        vacc = vacc.rate.list0[vacc.id]
        bb = betalist[bbid]
        out1 = SS_SIRS_prop_vacc_v3(sum(yinit), beta.param = bb, gamma.param = gam, omega_RL = omega_RLtmp, omega_LH = omega_LHtmp, v = vacc,
                                    h_S = h_S, l_S = l_S)
        if(ii == 1){
          # 1. Daily infection cases
          ylim.max = 0.25
          ISlist[bbid, vacc.id] = gam * (out1$I_M_inf + out1$I_S_inf) / sum(yinit) * 100
          
        }else if(ii == 2){
          # 2. ratio for progression of severe disease, i.e., I_S / (I_M + I_S)
          ylim.max = 4
          if((out1$I_M_inf + out1$I_S_inf) == 0){
            ISlist[bbid, vacc.id] = NA
          }else{
            ISlist[bbid, vacc.id] = out1$I_S_inf/(out1$I_M_inf + out1$I_S_inf) * 100
          }
        }else if(ii == 3){
          # 1. Daily severe cases
          ylim.max = 30*10^(-4)
          ISlist[bbid, vacc.id] = gam * out1$I_S_inf / sum(yinit) * 100
        }else{
          if (out1$S2_inf == 0){
            ISlist[bbid, vacc.id] = NA
          }else{
            ISlist[bbid, vacc.id] = out1$S_H_inf / out1$S_L_inf * 100
          }
        }
      }
    }
    
    plot(betalist / gam * (omega_RLtmp/(vacc.rate.list0[1] + omega_RLtmp)), ISlist[,1], type = "l", col = 1, 
         ylim = c(0, ylim.max), ylab = "", xlim = c(1, 8), lwd = 1, cex.lab = 2)
    collist = matrix(c(0,0,0,
                       204,121,167,
                       0,114,175,
                       82,167,221,
                       0,158,115,
                       213,94,0)/255, ncol = 3, byrow=TRUE)
    for(colid in 2:length(vacc.rate.list0)){
      lines(betalist / gam * (omega_RLtmp/(vacc.rate.list0[colid] + omega_RLtmp)), ISlist[,colid], col = rgb(collist[colid,1],collist[colid,2],collist[colid,3]), lwd = 1)
    }
  }
}
# mtext(text = "Basic reproduction number (R0)", cex = 1.5, side = 1, outer = TRUE)
# mtext(text = "Daily severe cases (%)", cex = 1, side = 2, adj = 0, outer = TRUE)
# mtext(text = "Severity rate (%)", cex = 1, side = 2, adj = 0.5, outer = TRUE)
# mtext(text = "Daily infection cases (%)", cex = 1, side = 2, adj = 1, outer = TRUE)

#### Figure for time course -- FIGURE 3 in the paper ==

# Time course of daily infection cases -- FIGURE 3A, 3D, S3 in the paper ==
par(mfrow = c(1,3))
gam = 1/10
vacc = 0.001
omega_RLparam = 1/365
omega_LHparam = 1/365 * 1/3
h_S = 0.05
l_S = 0.0025
times = seq(from = 0, to = 5000, by = 0.5)
betalist = c(seq(from = 0.16, by = 0.02, to = 0.3))
xlim_list = matrix(c(0, 2000, 
                     0, 2000, 
                     0, 2000), ncol = 2, byrow = TRUE)
ylim_list = matrix(c(0, 2.5, 
                     0, 1, 
                     0, 1), ncol = 2, byrow = TRUE)

for(imm_init_id in 1:3){
  imm_init = c(0.1, 0.5, 0.8)[imm_init_id]
  yinit_trans_v2 = c(S_H = (1-imm_init)*5*10^7, S_L = 0, I_M = 4750, I_S = 250, R = imm_init*5*10^7)

  Imat = matrix(NA, nrow = length(times), ncol = length(betalist))
  
  for(bid in 1:length(betalist)){
    bb = betalist[bid]
    pars_trans_v2 = c(beta.param = bb, gamma.param = gam, omega_RL = omega_RLparam, 
                      omega_LH = omega_LHparam, v = vacc, h_S, l_S, Ntot = sum(yinit_trans_v2))
    out_trans = as.data.frame(ode(yinit_trans_v2, times, ODE_SIRS_prop_vacc_v3, pars_trans_v2))
    Imat[, bid] = gam *(out_trans$I_M + out_trans$I_S)/sum(yinit_trans_v2) * 100
  }
  
  for(bid in 1:length(betalist)){
    if(bid == 1){
      plot(times,  Imat[, bid], type = "l", col =  rgb(bid/length(betalist), 0, 0,maxColorValue = 1), 
           xlim = xlim_list[imm_init_id,], ylim = ylim_list[imm_init_id,], 
           xlab = "Time (day)", ylab = paste("Daily infection cases (Initial R: ",imm_init*100,"%)", sep = ""), cex.lab = 1.5)
    }else{
      lines(times, Imat[, bid], type = "l", col =  rgb(bid/length(betalist), 0, 0,maxColorValue = 1))
    }
  }
}

# Time course of severity rate -- FIGURE 3B,3E, S3 in the paper ==
par(mfrow = c(1,3))
gam = 1/10
vacc = 0.001
omega_RLparam = 1/365
omega_LHparam = 1/365 * 1/3
h_S = 0.05
l_S = 0.0025
times = seq(from = 0, to = 5000, by = 0.5)
betalist = c(seq(from = 0.16, by = 0.02, to = 0.3))
xlim_list = matrix(c(0, 2000, 
                     0, 2000, 
                     0, 2000), ncol = 2, byrow = TRUE)
ylim_list = matrix(c(0, 5, 
                     0, 5, 
                     0, 5), ncol = 2, byrow = TRUE)

for(imm_init_id in 1:3){
  imm_init = c(0.1, 0.5, 0.8)[imm_init_id]
  yinit_trans_v2 = c(S_H = (1-imm_init)*5*10^7, S_L = 0, I_M = 4750, I_S = 250, R = imm_init*5*10^7)
  
  Imat = matrix(NA, nrow = length(times), ncol = length(betalist))
  
  for(bid in 1:length(betalist)){
    bb = betalist[bid]
    pars_trans_v2 = c(beta.param = bb, gamma.param = gam, omega_RL = omega_RLparam, 
                      omega_LH = omega_LHparam, v = vacc, h_S, l_S, Ntot = sum(yinit_trans_v2))
    out_trans = as.data.frame(ode(yinit_trans_v2, times, ODE_SIRS_prop_vacc_v3, pars_trans_v2))
    Imat[, bid] = out_trans$I_S/(out_trans$I_M + out_trans$I_S) * 100
  }
  
  for(bid in 1:length(betalist)){
    if(bid == 1){
      plot(times,  Imat[, bid], type = "l", col =  rgb(bid/length(betalist), 0, 0,maxColorValue = 1), 
           xlim = xlim_list[imm_init_id,], ylim = ylim_list[imm_init_id,], 
           xlab = "Time (day)", ylab = paste("Severity rate (Initial R: ",imm_init*100,"%)", sep = ""), cex.lab = 1.5)
    }else{
      lines(times, Imat[, bid], type = "l", col =  rgb(bid/length(betalist), 0, 0,maxColorValue = 1))
    }
  }
}



#### Figure for sever infectious -- FIGURE 3C in the paper

# Time course of daily sever cases -- FIGURE 3C,3F,3H, S3 in the paper ==
par(mfrow = c(1,3))
gam = 1/10
vacc = 0.001
omega_RLparam = 1/365
omega_LHparam = 1/365 * 1/3
h_S = 0.05
l_S = 0.0025
times = seq(from = 0, to = 5000, by = 0.5)
betalist = c(seq(from = 0.16, by = 0.02, to = 0.3))
xlim_list = matrix(c(0, 2000, 
                     0, 2000, 
                     0, 2000), ncol = 2, byrow = TRUE)
ylim_list = matrix(c(0, 1.1, 
                     0, 0.25, 
                     0, 0.12), ncol = 2, byrow = TRUE)

for(imm_init_id in 1:3){
  imm_init = c(0.1, 0.5, 0.8)[imm_init_id]
  yinit_trans_v2 = c(S_H = (1-imm_init)*5*10^7, S_L = 0, I_M = 4750, I_S = 250, R = imm_init*5*10^7)
  
  Imat = matrix(NA, nrow = length(times), ncol = length(betalist))
  
  for(bid in 1:length(betalist)){
    bb = betalist[bid]
    pars_trans_v2 = c(beta.param = bb, gamma.param = gam, omega_RL = omega_RLparam, 
                      omega_LH = omega_LHparam, v = vacc, h_S, l_S, Ntot = sum(yinit_trans_v2))
    out_trans = as.data.frame(ode(yinit_trans_v2, times, ODE_SIRS_prop_vacc_v3, pars_trans_v2))
    Imat[, bid] = out_trans$I_S/sum(yinit_trans_v2) * 100
  }
  
  for(bid in 1:length(betalist)){
    if(bid == 1){
      plot(times,  Imat[, bid], type = "l", col =  rgb(bid/length(betalist), 0, 0,maxColorValue = 1), 
           xlim = xlim_list[imm_init_id,], ylim = ylim_list[imm_init_id,], 
           xlab = "Time (day)", ylab = paste("Daily severe cases (Initial R: ",imm_init*100,"%)", sep = ""), cex.lab = 1.5)
    }else{
      lines(times, Imat[, bid], type = "l", col =  rgb(bid/length(betalist), 0, 0,maxColorValue = 1))
    }
  }
}



# Time to reach steay state -- FIGURE 3G in the paper ==
par(mfrow = c(1,1))
gam = 1/10
vacc = 0.001
omega_RLparam = 1/365
omega_LHparam = 1/365 * 1/3
h_S = 0.05
l_S = 0.0025
times = seq(from = 0, to = 5000, by = 0.5)
betalist = c(seq(from = 0.16, by = 0.02, to = 0.3))
xlim_list = matrix(c(0, 2000, 
                     0, 2000, 
                     0, 2000), ncol = 2, byrow = TRUE)
ylim_list = matrix(c(0, 1.1, 
                     0, 0.25, 
                     0, 0.12), ncol = 2, byrow = TRUE)

ss_thres = 0.3

imm_init = 0.8
yinit_trans_v2 = c(S_H = (1-imm_init)*5*10^7, S_L = 0, I_M = 4750, I_S = 250, R = imm_init*5*10^7)

Imat = matrix(NA, nrow = length(times), ncol = length(betalist))
I_res_mat = matrix(NA, nrow = length(times), ncol = length(betalist))
TTS_list = rep(NA, length(betalist))

for(bid in 1:length(betalist)){
  bb = betalist[bid]
  pars_trans_v2 = c(beta.param = bb, gamma.param = gam, omega_RL = omega_RLparam, 
                    omega_LH = omega_LHparam, v = vacc, h_S, l_S, Ntot = sum(yinit_trans_v2))
  
  out_trans = as.data.frame(ode(yinit_trans_v2, times, ODE_SIRS_prop_vacc_v3, pars_trans_v2))
  Imat[, bid] = out_trans$I_S
  out1 = SS_SIRS_prop_vacc_v3(sum(yinit_trans_v2), beta.param = bb, gamma.param = gam, omega_RL = omega_RLtmp, omega_LH = omega_LHtmp, v = vacc,
                              h_S = h_S, l_S = l_S)
  IS_SS = out1$I_S_inf
  I_res_mat[,bid] = abs(Imat[,bid] - IS_SS) / IS_SS
  for(bid in 1:length(betalist)){
    TTS_list[bid] = max(which(I_res_mat[,bid] > ss_thres))
  }
}

plot(betalist*10, TTS_list/365, xlab = "R0", ylab = "Time to steady state (year)", 
     ylim = c(1,12), cex.lab = 1.3, cex =2, pch = 16, col = rgb((1:8)/8, 0, 0,maxColorValue = 1))



