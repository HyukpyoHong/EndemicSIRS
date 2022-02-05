#### Functions
1. ODE_SIRS_prop_vacc_v3(Time, State, Pars)
 - This function simulate the time course of the modified SIRS model (see Fig. 1C and supplementary methods) for given parameters.
 - Detail description is included in the funciton file as comments.

2. SS_SIRS_prop_vacc_v3(Ntot, beta.param, gamma.param, omega_RL, omega_LH, v, h_S, l_S)
 - This function calculate the steady-state values of the modified SIRS model (see Fig. 1C and supplementary methods) for given parameters.
 - Detail description is included in the funciton file as comments.

Based on the above two functions, the rest part of the code generate the figures in the main text and supplementary materials. One can generates the figures by seqeuntially running the relevant section of the codes.