# This code implements agent-based simulation of growing and dividing cells and stochastic gene expression.
# a lineage of cells (like mother machine) is run and inside the cells and a Pol2 based model of non-bursty expression of a single gene is modeled
# Some figures are produced. 
# The results are used in Figure 3 of Sun, Bowman et al 2019 (bioRxiv, 754788), code by Vahid Shahrezaei 2019

rm(list=ls())


pom1_posterior <- read.csv("~/Box Sync/Fish project figures/Final push/pom1/prefit_pom1_posterior.csv", header=TRUE)
wt_posterior <- read.csv("~/Box Sync/Fish project figures/Final push/pom1/prefit_wt_posterior.csv", header=TRUE)


# Cell-Nuclear Noisy linear Map parameters for pom1 deleted cells from ABC inference on the cell and nuclear size data
# NLM parameters for cell size
a = 0.892279 # a of noisy linear map, a=1 is adder, a=0 is sizer
b = 5.595312 # b of noisy linear map, average birth size is b/(2-a)
noise_1 = 1.191693 # noise parameter in the added size
noise_2 = 0.07747589 # noise parameter of division
noise_3 = 0.1045408 # noise in growth rate
# NLM parameters for nucleaus size
a_N = 0.4675869 # a of noisy linear map, a=1 is adder, a=0 is sizer
b_N = 3.951028 # b of noisy linear map, average birth size is b/(2-a)
noise_1_N = 0.3910837 # noise parameter in the added size
noise_2_N = 0.03177691 # noise parameter of division
corr_LN = 0.5 #correlation between the birth cell size and nuclear size estimated from the time-lapse imaging. 

N_CellCycle = 5000 # Number of cell cycles to simulate
d_p = log(2)/2.5 # Pombe Growth rate (h^-1); which is equal to protein dilution rate

# Transcription parameters
d_m = log(2)/0.6 # mRNA decay rate
v_m = log(2)/0.02 # transcription rate
 

# Pol II model parameters
v_Pol2 = 4000*d_p/(b/(2-a)) # rate of production (translation) of Pol2 proportional to cell size picked to achive average of 4000 Pol2 based on litrature estimates
t_Pol2 = 2/(b_N/(2-a_N)) # rate of transport of Pol2 proportional to nuclear size pickd to be fast 
K_A_D = 10*(b_N/(2-a_N)) # K_A (association rate) of Pol2 binding times total DNA

L_birth = array(dim = N_CellCycle) # the array of initial sizes
L_division = array(dim = N_CellCycle) # the array of division sizes
gamma = array(dim = N_CellCycle) # the array of growth rates
T_CellCycle = array(dim = N_CellCycle) # the array of cell cycle times

L_birth = array(dim = N_CellCycle) # the array of initial sizes
L_division = array(dim = N_CellCycle) # the array of division sizes
L_birth_N = array(dim = N_CellCycle) # the array of initial nuclear sizes
L_division_N = array(dim = N_CellCycle) # the array of division nuclear sizes



# Setting up the cell cycle and nuclear size parameters for the whole lineage of cells
Noise_1 = rnorm(N_CellCycle, mean = 0, sd = noise_1)
gamma = d_p*rtnorm(N_CellCycle, mean = 1, a = 0, b= Inf, sd = noise_3)
Noise_1_N = rnorm(N_CellCycle, mean = 0, sd = noise_1_N)
MVN = rbvnorm(N_CellCycle, mean1 = 0.5, mean2 = 0.5, sd1 = noise_2, sd2 = noise_2_N, cor = corr_LN)
MVN[MVN<0.2] = 0.2
MVN[MVN>0.8] = 0.8
Noise_2 = MVN[,1]
Noise_2_N = MVN[,2]

L_div_0 = 2*b/(2-a)
L_div_0_N = 2*b_N/(2-a_N)

for (i in 1:N_CellCycle){
  L_birth[i] = L_div_0*Noise_2[i]
  L_division[i] = a*L_birth[i] + b + Noise_1[i]
  if (L_division[i]<L_birth[i]) 
    { L_division[i] = L_birth[i] + L_div_0/2}
  L_birth_N[i] = L_div_0_N*Noise_2_N[i]
  L_division_N[i] = a_N*L_birth_N[i] + b_N + Noise_1_N[i]
  L_div_0 = L_division[i]
  L_div_0_N = L_division_N[i]
}
T_CellCycle = 1/gamma*log(L_division/L_birth)  


## Main loop

N_t = 100 # number of time-steps per CellCycle
mRNA = 1.5*v_m/d_m # initial condition
Pol2_N =  v_Pol2/d_p * b/(2-a) # initial condition
Pol2_N_b = K_A_D/(K_A_D + b_N/(2-a_N))*Pol2_N # initial condition
Pol2 = 0 # initial condition

t_tot = array(dim = N_CellCycle*N_t)
M_tot = array(dim = N_CellCycle*N_t)
L_tot = array(dim = N_CellCycle*N_t)
L_N_tot = array(dim = N_CellCycle*N_t)
Pol2_N_b_tot = array(dim = N_CellCycle*N_t)



for(i in 1:N_CellCycle){
    t = 0
    L = L_birth[i]
    L_N = L_birth_N[i]
    j = 0
    mRNA = rbinom(1, mRNA, Noise_2[i]) #binomial partioning of mRNAs
    Pol2 = rbinom(1, round(Pol2), Noise_2[i]) #partitioning of cytosolic Pol2
    Pol2_N = rbinom(1, round(Pol2_N - Pol2_N_b), Noise_2_N[i]) #partitioning of free nuclear Pol2 is nuclear voloume dependent
    Pol2_N_b = rbinom(1, round(Pol2_N_b), 0.5) #partitioning of bound Pol2 equally divided as on DNA
    Pol2_N = Pol2_N_b + Pol2_N  
    while (t < T_CellCycle[i]){
      a_d_m = mRNA* d_m # propensity of mRNA decay
      a_v_m = v_m * Pol2_N_b/(v_Pol2*(b/(2-a)))*d_p # propensity of transcription in Pol2 model is propotional to the normalised levle of Pol2_N_b 
      a_tot = a_d_m + a_v_m
      tau = rexp(1, a_tot)
      if (t+tau <(j+1)*T_CellCycle[i]/N_t) {
        t = t + tau
        r = runif(1, 0, a_tot)
        if (r < a_d_m) {
          mRNA = mRNA - 1}
        else {
          mRNA = mRNA + 1}
      }
      else {
        j = j+1
        t = j*T_CellCycle[i]/N_t 
        L = L * exp(gamma[i]*T_CellCycle[i]/N_t)
        L_N = L_N * exp(gamma[i]*T_CellCycle[i]/N_t)
        Pol2_0 = Pol2
        Pol2 = Pol2 + (v_Pol2*L- Pol2*t_Pol2*L_N)*T_CellCycle[i]/N_t  # Descretised Pol2 ODE: is produced by a rate proptional to cell size and transported in by a rate proportional to nuclear size
        Pol2_N = Pol2_N + Pol2_0*t_Pol2*L_N*T_CellCycle[i]/N_t # Descretised Pol2_N ODE: transported from cytosol to nucleus with a rate proportional to nuclear size
        Pol2_N_b = K_A_D/(K_A_D + L_N)*Pol2_N  #Assuming fast binding to DNA compared to cell cycle using the quasi equilibrium approximation
        t_tot[(i-1)*N_t + j]= sum(T_CellCycle[1:i-1]) +  j* T_CellCycle[i]/N_t
        M_tot[(i-1)*N_t + j]=mRNA
        L_tot[(i-1)*N_t + j]=L
        L_N_tot[(i-1)*N_t + j]=L_N
        Pol2_N_b_tot [(i-1)*N_t + j]= Pol2_N_b
      }
    }
}

# sample N_CellCycle observations from the whole lineage ignoring the first 10 cell cycles.   
index = sample(10*N_t:N_CellCycle*N_t, N_CellCycle, replace = FALSE)

par(mfrow = c(2, 2))
plot(L_tot[index], M_tot[index],xlab =  "Length", ylab = "mRNA")
plot(L_tot[index], Pol2_N_b_tot[index], xlab =  "Length", ylab = "DNA bound PolII")
plot(L_tot[index], M_tot[index]/L_tot[index], xlab =  "Length", ylab = "mRNA concentration", main = "pom1 rpb1")
plot(L_N_tot[index], M_tot[index]/L_N_tot[index], xlab =  "Nuclear Area", ylab = "mRNA/Nuclear area")



cat(cor(L_tot[index], M_tot[index]/L_tot[index]))

cat(cor(L_N_tot[index], M_tot[index]/L_N_tot[index]))

