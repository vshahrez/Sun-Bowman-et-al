# This code implements agent-based simulation of growing and dividing cells and stochastic gene expression.
# a lineage of cells (like mother machine) is run and inside the cells a TASEP model of non-bursty expression of a single gene is modeled
# some statistics of model output is recorded for differenent choices of model parameters and three different models of initation, elongation or termination scaling with cell size
# The results are used in Figure 3 of Sun, Bowman et al 2019 (bioRxiv, 754788), code by Vahid Shahrezaei 2019

rm(list=ls())


# Noisy linear map parameters for WT cells from ABC inference on the data
a = 0.6191072 # a = 1 would be adder, pombe is between an adder and sizer
b = 8.493817 # b of noisy linear map in micrometers, average birth size is b/(2-a)
noise_1 = 0.9984188 # noise parameter in the added size
noise_2 = 0.04342351 # noise parameter of division
noise_3 = 0.1463003 # noise in growth rate

N_CellCycle = 300 # Number of cell cycles to simulate
d_p = log(2)/2.5 # Pombe Growth rate (h^-1); which is equal to protein dilution rate

# Model 1 Initiation scaling, 2 Elongation scaling, 3 Termination scaling 
model = 1

L_birth = array(dim = N_CellCycle) # the array of initial sizes
L_division = array(dim = N_CellCycle) # the array of division sizes
gamma = array(dim = N_CellCycle) # the array of growth rates
T_CellCycle = array(dim = N_CellCycle) # the array of cell cycle times

Output = array(dim = c(200, 17, 3)) # the main output: contains 17 features for 200 parameters in the 3 models of initiation, elongation or termination scaling 
colnames(Output) = c('L', 'M', 'N', 'cor', 'slope', 'slope_e', 'In', 'El', 'Tr', 'cor_threshold', 'P', 'cor_P', 'Gene_1_s', 'Gene_1_l', 'Gene_2_s', 'Gene_2_l', 'cor_freq')


# simulation parameters and variables
param = c(0, 0, 0)

N_t = 100 # number of time-steps per CellCycle
L_G = 40 # number of TASEP sites
GENE = array(0, dim = L_G)
t_tot = array(dim = N_CellCycle*N_t+1)
M_tot = array(dim = N_CellCycle*N_t+1)
L_tot = array(dim = N_CellCycle*N_t+1)
Nascent_tot = array(dim = N_CellCycle*N_t+1)
Pol_tot = array(dim = N_CellCycle*N_t+1)
GENE_tot = array(dim = c(N_CellCycle*N_t+1, L_G))


# Main loop 3 models, 200 parameters
for (model in (1:3)){
q = 0
while (q < 201){

# prior over the parameter sets  
param[1] = runif(1, min = 0.001, max = 0.1)
param[2] = runif(1, min = 0.001, max = 0.1)
param[3] = runif(1, min = 0.001, max = 0.1)
# Transcroption parameters
L_G = 20 # number of TASEP sites
d_m = 1/0.6 # mRNA decay rate
In = 1/param[1] # Initiation rate
El = 1/param[2]*(L_G -1) # Elongation rate (for per cite multipy by (L_G - 1)
Tr = 1/param[3] # Termination rate


# Setting up the cell cycle parameters for the whole lineage of cells
Noise_1 = rnorm(N_CellCycle, mean = 0, sd = noise_1)
Noise_2 = rtnorm(N_CellCycle, mean = 0.5, a = 0, b = 1, sd = noise_2)
gamma = d_p*rtnorm(N_CellCycle, mean = 1, a = 0, b= Inf, sd = noise_3)
L_div_0 = 2*b/(2-a)
for (i in 1:N_CellCycle){
  L_birth[i] = L_div_0*Noise_2[i]
  L_division[i] = a*L_birth[i] + b + Noise_1[i]
  L_div_0 = L_division[i]
}
T_CellCycle = 1/gamma*log(L_division/L_birth)  


## Main cell cycle loop - first run for 5 cell cycles to check if the parameters chosen satisfy the conditions set 

mRNA = 0

for(i in 1:5){
    t = 0
    L = L_birth[i]
    j = 0
    mRNA = rbinom(1, mRNA, Noise_2[i]) #binomial partiioning of mRNAs
    
    # stchastic gene expression during a cell cycle
    while (t < T_CellCycle[i]){
      a_d_m = mRNA* d_m # propensity of mRNA decay
      a_In = (1-GENE[1])*In # propensity of initiation
      a_Tr = GENE[L_G]*Tr # propensity of termination
      Elon = GENE[1:L_G-1]*(1-GENE[2:L_G]) 
      a_El = sum(Elon)*El # propensity of elongation
      if (model == 1) a_In = a_In *L/(1.5*b/(2-a))
      if (model == 2) a_El = a_El *L/(1.5*b/(2-a))
      if (model ==3) a_tr = a_Tr *L/(1.5*b/(2-a))
      a_tot = a_d_m + a_In + a_Tr + a_El
      tau = rexp(1, a_tot)
      if (t+tau <(j+1)*T_CellCycle[i]/N_t) {
        t = t + tau
        r = runif(1, 0, a_tot)
        if (r < a_d_m) {
          mRNA = mRNA - 1}
        else if (r < a_d_m + a_El) {
          index_GENE = which(Elon == 1)
          GENE[index_GENE[ceiling((r - a_d_m)/a_El*sum(Elon))]]= 0
          GENE[index_GENE[ceiling((r - a_d_m)/a_El*sum(Elon))]+1]= 1
        }
        else if (r < a_d_m + a_El+ a_In) {
          GENE[1] = 1
        }
        else {
          GENE[L_G] = 0
          mRNA = mRNA + 1
        }
      }
      else{
        j = j+1
        t = j*T_CellCycle[i]/N_t 
        L = L * exp(gamma[i]*T_CellCycle[i]/N_t)
        t_tot[(i-1)*N_t + j]= sum(T_CellCycle[1:i-1])+  j* T_CellCycle[i]/N_t
        M_tot[(i-1)*N_t + j]=mRNA
        L_tot[(i-1)*N_t + j]=L
        Nascent_tot[(i-1)*N_t + j]=sum(GENE*seq(1:L_G))/L_G
        Pol_tot[(i-1)*N_t + j]=sum(GENE)
        GENE_tot[(i-1)*N_t + j,]= GENE
      }
    }
  
    GENE = sample(0:1, L_G, replace = TRUE)*GENE # nascent RNAs remain in daughter with prob = 0.5
}

# check if the nascent intensity is about 5 and mean mRNA numbers about 25 (these set as values comparable to what is observed for differnet genes)
# then continue with simulating the rest of the cell cycles for this parameter choice else pick a new random parameter set
if ((abs(mean(Nascent_tot[N_t:(5*N_t)])-5)<6) & ((abs(mean(M_tot[N_t:(5*N_t)]) - 25) < 5)))
{ 
for(i in 6:N_CellCycle){
  t = 0
  L = L_birth[i]
  j = 0
  mRNA = rbinom(1, mRNA, Noise_2[i]) #binomial partiioning of mRNAs
  
  # stchastic gene expression during a cell cycle
  while (t < T_CellCycle[i]){
    a_d_m = mRNA* d_m # propensity of mRNA decay
    a_In = (1-GENE[1])*In # propensity of initiation
    a_Tr = GENE[L_G]*Tr # propensity of termination
    Elon = GENE[1:L_G-1]*(1-GENE[2:L_G]) 
    a_El = sum(Elon)*El # propensity of elongation
    if (model == 1) a_In = a_In *L/(1.5*b/(2-a))
    if (model == 2) a_El = a_El *L/(1.5*b/(2-a))
    if (model ==3) a_tr = a_Tr *L/(1.5*b/(2-a))
    a_tot = a_d_m + a_In + a_Tr + a_El
    tau = rexp(1, a_tot)
    if (t+tau <(j+1)*T_CellCycle[i]/N_t) {
      t = t + tau
      r = runif(1, 0, a_tot)
      if (r < a_d_m) {
        mRNA = mRNA - 1}
      else if (r < a_d_m + a_El) {
        index_GENE = which(Elon == 1)
        GENE[index_GENE[ceiling((r - a_d_m)/a_El*sum(Elon))]]= 0
        GENE[index_GENE[ceiling((r - a_d_m)/a_El*sum(Elon))]+1]= 1
      }
      else if (r < a_d_m + a_El+ a_In) {
        GENE[1] = 1
      }
      else {
        GENE[L_G] = 0
        mRNA = mRNA + 1
      }
    }
    else{
      j = j+1
      t = j*T_CellCycle[i]/N_t 
      L = L * exp(gamma[i]*T_CellCycle[i]/N_t)
      t_tot[(i-1)*N_t + j]= sum(T_CellCycle[1:i-1])+  j* T_CellCycle[i]/N_t
      M_tot[(i-1)*N_t + j]=mRNA
      L_tot[(i-1)*N_t + j]=L
      Nascent_tot[(i-1)*N_t + j]=sum(GENE*seq(1:L_G))/L_G
      Pol_tot[(i-1)*N_t + j]=sum(GENE)
      GENE_tot[(i-1)*N_t + j,]= GENE
    }
  }
  
  GENE = sample(0:1, L_G, replace = TRUE)*GENE # nascent RNAs remain in daughter with prob = 0.5
}

# sample N_CellCycle*5 observations from the whole lineage ignoring the first 10 cell cycles.   
index = sample(10*N_t:N_CellCycle*N_t, N_CellCycle*5, replace = FALSE)


mean_N = mean(Nascent_tot[index])
mean_M = mean(M_tot[index])
 
cor_n = cor(L_tot[index], Nascent_tot[index])
mean_L = mean(L_tot[index])
mean_Pol = mean(Pol_tot[index])

index_1 = index[L_tot[index] < mean_L]
index_2 = index[L_tot[index] > mean_L]
a1 = lm(M_tot[index_1]~L_tot[index_1])
a2 = lm(M_tot[index_2]~L_tot[index_2])

# Wrting the summery results for each model/parameter choice
Output[q, 1, model] = mean_L # mean cell length
Output[q, 2, model] = mean_M # mean mRNA number
Output[q, 3, model] = mean_N # mean Nascent intensity
Output[q, 4, model] = cor_n # correlation between cell length and nascent intensity
Output[q, 5, model] = (summary(a1)$coef[2,1]-summary(a2)$coef[2,1])/mean_M # normalised differnce between the regression slope of the first half and second half of the data, if it is close to zere the relationship between mRNAs and cell length is linear
Output[q, 6, model] = (summary(a1)$coef[2,1]-summary(a1)$coef[2,2]-summary(a2)$coef[2,1]-summary(a2)$coef[2,2])/mean_M #like above but taking into account the confidence intervals of the slopes
Output[q, 7, model] = param[1] # initiation rate
Output[q, 8, model] = param[2] # elongation rate
Output[q, 9, model] = param[3] # termination rate
Nascent_tot[Nascent_tot<2.5] = 0 
cor_n = cor(L_tot[index], Nascent_tot[index]) 
Output[q, 10, model] = cor_n # correlation between cell length and nascent intensity for spots brighter than the threshold of 2.5 as used in data analysis
Output[q, 11, model] = mean_Pol # mean number of Pol2 on the gene
Output[q, 12, model] = cor(L_tot[index], Pol_tot[index]) # correlation between cell length and number of Pol2 on the gene
GENE_1 = colMeans(GENE_tot[index_1,])
GENE_2 = colMeans(GENE_tot[index_2,])
Output[q, 13, model] = sum(GENE_1[1:(L_G/2)]) # mean occupancy of the 5' end of the gene in small cells
Output[q, 14, model] = sum(GENE_1[(L_G/2 + 1): L_G]) # mean occupancy of the 3' end of the gene in small cells
Output[q, 15, model] = sum(GENE_2[1:(L_G/2)]) # mean occupancy of the 5' end of gene in large cells
Output[q, 16, model] = sum(GENE_2[(L_G/2 + 1): L_G]) # mean occupancy of the 3' end of gene in large cells

Nascent_Freq = array(10)
Nascent_Freq[1] = mean(Nascent_tot[which((L_tot < 5.5))]>2.5)
Nascent_Freq[2] = mean(Nascent_tot[which((L_tot > 5.5) & (L_tot < 6.5))]>2.5)
Nascent_Freq[3] = mean(Nascent_tot[which((L_tot > 6.5) & (L_tot < 7.5))]>2.5)
Nascent_Freq[4] = mean(Nascent_tot[which((L_tot > 7.5) & (L_tot < 8.5))]>2.5)
Nascent_Freq[5] = mean(Nascent_tot[which((L_tot > 8.5) & (L_tot < 9.5))]>2.5)
Nascent_Freq[6] = mean(Nascent_tot[which((L_tot > 9.5) & (L_tot < 10.5))]>2.5)
Nascent_Freq[7] = mean(Nascent_tot[which((L_tot > 10.5) & (L_tot < 11.5))]>2.5)
Nascent_Freq[8] = mean(Nascent_tot[which((L_tot > 11.5) & (L_tot < 12.5))]>2.5)
Nascent_Freq[9] = mean(Nascent_tot[which((L_tot > 12.5) & (L_tot < 13.5))]>2.5)
Nascent_Freq[10] = mean(Nascent_tot[which(L_tot > 13.5)]>2.5)
Output[q, 17, model] = cor(Nascent_Freq, 5:14) # correlation between cell length and Nascent frequency using the binning of cell size 
q = q+1
cat(q)
}
}
}

save(Output, file = "~/PATH/Tasep_output.RData")


