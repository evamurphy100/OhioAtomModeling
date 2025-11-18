### Factor analysis on Ohio county and zcta level data in 2020
library(dplyr)
library(readxl)
library(zctaCrosswalk)
library(writexl)
library(tigris)

# Load Atom data
atom_data = read_excel("Atom_data.xlsx")
# Load the crosswalk data
crosswalk = read_excel("Data_redo/Ohio_Crosswalk.xlsx")
crosswalk = crosswalk[,c(1,5)]

########## Make map from atom to county ##########

#### make a vector of the number of atoms in each county
#### and then a longer vector that is the indices for each atom for each county
n.cnty = length(unique.fips)
atom.per.cnty = aggregate(atom_data$Atom, by = list(atom_data$numeric.fips), FUN = length)
atom.per.cnty = atom.per.cnty$x
atom.per.cnty = c(0, atom.per.cnty)
cty.atom = rep(NA,sum(atom.per.cnty))
for(i in 1:n.cnty){
  #save the atom index corresponding to cnty = i
  cty.atom[(sum(atom.per.cnty[1:(i)])+1):(sum(atom.per.cnty[1:(i)])+atom.per.cnty[i+1])] = unique(atom_data$Atom[which(atom_data$numeric.fips==i)])
}


### map for what atoms go to what counties
acmap <- atom_data$numeric.fips

########## Make map from atom to zcta ##########

#### make a vector of the number of atoms in each zcta
#### and then a longer vector that is the indices for each atom for each zcta
n.zcta = length(unique.ZCTA)
atom.per.zcta = aggregate(atom_data$Atom, by = list(atom_data$numeric.ZCTA), FUN = length)
atom.per.zcta = atom.per.zcta$x
atom.per.zcta = c(0, atom.per.zcta)
zcta.atom = rep(NA,sum(atom.per.zcta))
for(i in 1:n.zcta){
  #save the atom index corresponding to zcta = i
  zcta.atom[(sum(atom.per.zcta[1:(i)])+1):(sum(atom.per.zcta[1:(i)])+atom.per.zcta[i+1])] = unique(atom_data$Atom[which(atom_data$numeric.ZCTA==i)])
}

### map for what atoms go to what ZCTA
azmap <- atom_data$numeric.ZCTA


###### Make weights ###########

### make the atom weights that are constant over time
## assume the weights are the proportion of population within 
## each atom relative to county
## each row a county and column an atom
n.cnty = length(unique.fips)
n.atom = nrow(atom_data)
wts.AC <- matrix(0,n.cnty,n.atom)
for(i in 1:n.cnty){
  Ind = which(atom_data$numeric.fips == i)
  wts.AC[i,Ind] = atom_data$Atom_pop[Ind]/atom_data$County_pop[Ind]
}

### assume the weights are the proportion of population within 
### each atom relative to county
### each row a county and column an atom
n.zcta = length(unique.ZCTA)
wts.AZ <- matrix(0, n.zcta, n.atom)
for(i in 1:n.zcta){
  Ind = which(atom_data$numeric.ZCTA == i)
  wts.AZ[i,Ind] = atom_data$Atom_pop[Ind]/atom_data$ZCTA_pop[Ind]
}

############ ZCTA level data ############ 

## Load naloxone data
naloxone.data = read_excel("Naloxone.xlsx")
naloxone.data = na.omit(naloxone.data)

## Combine naloxone with the zcta data
naloxone.data = full_join(naloxone.data, crosswalk, by = "Zip")
naloxone.data = naloxone.data %>%
  group_by(ZCTA) %>%
  summarize(Total = sum(Nalox_count, na.rm = F))
## Order the naloxone data to correspond to the ordering of atom.data
naloxone.data = full_join(naloxone.data, unique(atom_data[,c(2,9)]), by = c("ZCTA" = "zcta"))
naloxone.data = naloxone.data[order(naloxone.data$numeric.ZCTA), ]

## Save the zcta counts in a new variable
y.N = naloxone.data$Total[1:n.zcta]


############ County level data ############ 

## Load the death counts 
death.data = read.csv("all_death_data_07-23.csv")
death.data = death.data %>%
  filter(DeathYear == 2020)
death.data = death.data[,c(2,3)]
death.data = full_join(death.data, cnty_pop, by = "CountyName")
death.data = death.data[,-c(1,4)]
death.data = death.data[c("county", "Deaths")]
death.data = death.data %>%
  mutate(county = as.numeric(county))
## Order the death data to correspond to the ordering of atom.data
death.data = full_join(death.data, unique(atom_data[,c(3,8)]), by = "county")
death.data = death.data[order(death.data$numeric.fips),]

y.D = as.numeric(death.data$Deaths[1:n.cnty])

## Load treatment counts
treatment.data = read_excel("Kline Opioid Dec2022.xlsx")
treatment.data = treatment.data[-90,c(2,4,5)]
colnames(treatment.data) = c("CountyName", "u20", "o20")
treatment.data <- treatment.data %>%
  mutate(u20 = ifelse(u20 == "*", NA, as.numeric(u20)),
         o20 = ifelse(o20 == "*", 0, as.numeric(o20)))
treatment.data = full_join(treatment.data, cnty_pop, by = "CountyName")
treatment.data = treatment.data[,2:4]
treatment.data = treatment.data[c("county", "u20", "o20")]
treatment.data = treatment.data %>%
  mutate(county = as.numeric(county))
## Order the treatment data to correspond to the ordering of atom.data
treatment.data = full_join(treatment.data, unique(atom_data[,c(3,8)]), by = "county")
treatment.data = treatment.data[order(treatment.data$numeric.fips),]
treatment.data = treatment.data[1:n.cnty, ]

y.T = treatment.data$u20 + treatment.data$o20

## Adress censored data
censC = 1*(is.na(treatment.data$u20))*(1-is.na(treatment.data$o20))
censA = 1*(1-is.na(treatment.data$u20))*(is.na(treatment.data$o20))
censBoth = 1*(is.na(treatment.data$u20))*(is.na(treatment.data$o20))
censNone = 1*(1-is.na(treatment.data$u20))*(1-is.na(treatment.data$o20))
cenlb = rep(NA,length(y.T))
cenub = rep(NA,length(y.T))
cenlb[which(censC==1)]=treatment.data$o20[which(censC==1)]
cenub[which(censC==1)]=treatment.data$o20[which(censC==1)]+9
cenlb[which(censA==1)]=treatment.data$u20[which(censA==1)]
cenub[which(censA==1)]=treatment.data$u20[which(censA==1)]+9
cenlb[which(censBoth==1)]=1
cenub[which(censBoth==1)]=18
cenlb[which(is.na(cenlb))]=y.T[which(is.na(cenlb))]
cenub[which(is.na(cenub))]=y.T[which(is.na(cenub))]

cen=censA+censC+censBoth
bd=cbind(cenlb,cenub)

### plug in values for censored right now
y.T[which(censC==1)] = treatment.data$o20[which(censC==1)]+5
y.T[which(censA==1)] = treatment.data$u20[which(censA==1)]+5
y.T[which(censBoth==1)] = 10
#y.T = y.T[1:n.cnty]

## Load HCV data
HCV.data = read_excel("HCV_2020.xlsx")
HCV.data = full_join(HCV.data, cnty_pop, by = "CountyName")
HCV.data = HCV.data[,3:4]
HCV.data = HCV.data[c("county", "Total")]
HCV.data = HCV.data %>%
  mutate(county = as.numeric(county))
## Order the HCV data to correspond to the ordering of atom.data
HCV.data = full_join(HCV.data, unique(atom_data[,c(3,8)]), by = "county")
HCV.data = HCV.data[order(HCV.data$numeric.fips),]

y.C = HCV.data$Total[1:88]

## Load HIV data
HIV.data = read_excel("HIV_2020.xlsx")
HIV.data = full_join(HIV.data, cnty_pop, by = "CountyName")
HIV.data = HIV.data[,c(2,3)]
HIV.data = HIV.data[c("county", "Total")]
HIV.data = HIV.data %>%
  mutate(county = as.numeric(county))
## Order the HIV data to correspond to the ordering of atom.data
HIV.data = full_join(HIV.data, unique(atom_data[,c(3,8)]), by = "county")
HIV.data = HIV.data[order(HIV.data$numeric.fips),]

y.I = HIV.data$Total[1:n.cnty]


## Make atom adjacency matrix
#library(corrplot)
#library(sf)
#library(dplyr)
#library(tidyverse)

adj.zcta.mat = as.matrix(read.csv("adj_zcta_mat.csv", row.names = 1, check.names = FALSE))

county_adj_matrix = read.csv("OhioAdjacency.csv", header = F)  # County adjacency matrix
zcta_county_mapping = read_excel("zipcode_county_matrix.xlsx")  # ZCTA-to-county mapping

# Initialize the Adjacency Matrix
num_atoms <- as.numeric(nrow(atom_data))
atom_adj_matrix <- matrix(0, nrow = num_atoms, ncol = num_atoms)

# Loop through each pair of atoms
for (i in 1:num_atoms) {
  for (j in 1:num_atoms) {
    if (i != j) {
      county_i <- atom_data$county[i]
      county_j <- atom_data$county[j]
      zcta_i <- atom_data$zcta[i]
      zcta_j <- atom_data$zcta[j]
      
      # Check if atoms share the same county and if their ZCTAs are neighbors
      if (county_i == county_j && adj.zcta.mat[as.character(zcta_i), as.character(zcta_j)] == 1) {
        atom_adj_matrix[i, j] <- 1
      }
      
      # Check if atoms are from different counties but share the same ZCTA
      if (county_i != county_j && zcta_i == zcta_j) {
        atom_adj_matrix[i, j] <- 1
      }
    }
  }
}

# Handle Special Cases: Neighboring counties and neighboring ZCTAs 
# that spill over into a neighboring county
for (i in 1:num_atoms) {
  for (j in 1:num_atoms) {
    if (i != j && atom_adj_matrix[i, j] == 0) {
      county_i <- atom_data$numeric.fips[i]
      county_j <- atom_data$numeric.fips[j]
      zcta_i <- atom_data$zcta[i]
      zcta_j <- atom_data$zcta[j]
      
      # Check if counties are neighbors and ZCTAs are neighbors
      if (county_adj_matrix[county_i, county_j] == 1 &&
          adj.zcta.mat[as.character(zcta_i), as.character(zcta_j)] == 1) {
        
        # Ensure there is no other atom with zcta_i in county_j and zcta_j in county_i
        if (!any(atom_data$zcta == zcta_i & atom_data$numeric.fips == county_j) &&
            !any(atom_data$zcta == zcta_j & atom_data$numeric.fips == county_i)) {
          atom_adj_matrix[i, j] <- 1
        }
      }
    }
  }
}


# Create a data frame with unique fips and corresponding county populations
Tot.Pop.county = atom_data[!duplicated(atom_data$numeric.fips), c("numeric.fips", "County_pop")]
Tot.Pop.county = Tot.Pop.county[order(Tot.Pop.county$numeric.fips), ]
# Extract the County_pop in the order of numeric.fips
Tot.Pop.county = Tot.Pop.county$County_pop
# Expected count standardized to the state level
ET = rep(0,n.cnty)
ED = rep(0,n.cnty)
EC = rep(0,n.cnty)
EI = rep(0,n.cnty)
for(i in 1:n.cnty){
    ET[i] =  sum(y.T)/sum(Tot.Pop.county)*Tot.Pop.county[i]
    ED[i] =  sum(y.D)/sum(Tot.Pop.county)*Tot.Pop.county[i]
    EC[i] =  sum(y.C)/sum(Tot.Pop.county)*Tot.Pop.county[i]
    EI[i] =  sum(y.I)/sum(Tot.Pop.county)*Tot.Pop.county[i]
}

# Create a data frame with unique zcta and corresponding zcta populations
Tot.Pop.zcta = atom_data[!duplicated(atom_data$numeric.ZCTA), c("numeric.ZCTA", "ZCTA_pop")]
Tot.Pop.zcta = Tot.Pop.zcta[order(Tot.Pop.zcta$numeric.ZCTA), ]
# Extract the ZCTA_pop in the order of numeric.ZCTA
Tot.Pop.zcta = Tot.Pop.zcta$ZCTA_pop
# Expected count standardized to the state level
EN =rep(0,n.zcta)
for(i in 1:n.zcta){
  EN[i] =  sum(y.N)/sum(Tot.Pop.zcta)*Tot.Pop.zcta[i]
}


############################################################
###########         Set up for NIMBLE          ############
############################################################
library(nimble)
library(coda)
atom_adj_matrix = as.data.frame(atom_adj_matrix)
num = colSums(atom_adj_matrix)
## ICAR variance
adj.ICAR = NULL
for(j in 1:num_atoms){
  adj.ICAR = c(adj.ICAR,which(atom_adj_matrix[j,]==1)) ##takes only the neighbors
}
adj.ICAR = as.vector(adj.ICAR)
num = as.vector(num)
weights.ICAR = 1+0*adj.ICAR

load.mat = as.vector(c(1,0,0,0,0))

model_code=nimbleCode({
  #county data
  for (i in 1:n.cnty) {
    y.D[i] ~ dpois(ED[i]*lambdaD[i])
    cen[i] ~ dinterval(y.T[i],bd[i,1:2])
    y.T[i] ~ dpois(ET[i]*lambdaT[i])
    y.C[i] ~ dpois(EC[i]*lambdaC[i])
    y.I[i] ~ dpois(EI[i]*lambdaI[i])
    log(lambdaT[i]) <- muT + inprod(wts.AC[i,1:num_atoms], load.mat[1]* U[1:num_atoms]) + VT[i]
    log(lambdaD[i]) <- muD + inprod(wts.AC[i,1:num_atoms], load.mat[2]* U[1:num_atoms]) + VD[i]
    log(lambdaC[i]) <- muC + inprod(wts.AC[i,1:num_atoms], load.mat[3]* U[1:num_atoms]) + VC[i]
    log(lambdaI[i]) <- muI + inprod(wts.AC[i,1:num_atoms], load.mat[4]* U[1:num_atoms]) + VI[i]
    VD[i] ~ dnorm(0,tau.VD) 
    VT[i] ~ dnorm(0,tau.VT) 
    VC[i] ~ dnorm(0,tau.VC) 
    VI[i] ~ dnorm(0,tau.VI) 
  }
  #ZCTA data
  for (i in 1:n.zcta) {
    y.N[i] ~ dpois(EN[i]*lambdaN[i])
    log(lambdaN[i]) <- muN + inprod(wts.AZ[i,1:num_atoms], load.mat[5]* U[1:num_atoms]) + VN[i]
    VN[i] ~ dnorm(0,tau.VN) 
  }
  
  # ICAR prior for the spatial factors
  U[1:num_atoms] ~ dcar_normal(adj[], weights[], num[], tau = tau.U, zero_mean=1)
  tau.U ~ dgamma(0.5,0.5)
  
  #Priors
  load.mat[2] ~ dnorm(0,sd=10)
  load.mat[3] ~ dnorm(0,sd=10)
  load.mat[4] ~ dnorm(0,sd=10)
  load.mat[5] ~ dnorm(0,sd=10)
  
  muT ~ dflat()
  muD ~ dflat()
  muC ~ dflat()
  muI ~ dflat()
  muN ~ dflat()
  
  tau.VD ~ dgamma(0.5,0.5)
  tau.VT ~ dgamma(0.5,0.5)
  tau.VC ~ dgamma(0.5,0.5)
  tau.VI ~ dgamma(0.5,0.5)
  tau.VN ~ dgamma(0.5,0.5)
  
})

############################################################
#########              Call NIMBLE                 ###########
############################################################

mod_constants=list(bd=bd,num_atoms=num_atoms,num=num,adj=adj.ICAR,weights=weights.ICAR,n.cnty=n.cnty,
                   n.zcta=n.zcta,wts.AC=wts.AC,wts.AZ=wts.AZ)

mod_data=list(cen=cen,y.D=y.D,y.T=y.T,y.C=y.C,y.I=y.I,y.N=y.N,
              ET=ET,ED=ED,EC=EC,EI=EI,EN=EN)

mod_inits=list(U=rep(0,num_atoms),load.mat=load.mat,muT=0,muD=0,muC=0,muI=0,muN=0,
               VD=rep(0,n.cnty),VT=rep(0,n.cnty),VC=rep(0,n.cnty),VI=rep(0,n.cnty),VN=rep(0,n.zcta),
               y.T=floor(.5*(bd[,1]+bd[,2])))

# Build the model.
nim_model <- nimbleModel(model_code,mod_constants,mod_data,mod_inits)
compiled_model <- compileNimble(nim_model,resetFunctions = TRUE)

# Set up samplers.
mcmc_conf <- configureMCMC(nim_model,monitors=c("muT","muD","muC","muI","muN","load.mat",
                                                "U","VD","VT","VC","VI", "VN", "lambdaT","lambdaD",
                                                "lambdaC","lambdaI", "lambdaN"),useConjugacy = TRUE)

mod_mcmc<-buildMCMC(mcmc_conf)
compiled_mcmc<-compileNimble(mod_mcmc, project = nim_model,resetFunctions = TRUE)


# Run the model 
MCS=50000
st<-Sys.time()
samples=runMCMC(compiled_mcmc,inits=mod_inits,
                nchains = 1, nburnin=floor(MCS/2),niter = MCS,samplesAsCodaMCMC = TRUE,thin=50,
                summary = FALSE, WAIC = FALSE,progressBar=TRUE) 
Sys.time()-st

save(samples,file="Samples/ModOutput50000.Rda")
