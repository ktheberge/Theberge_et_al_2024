#### H. GAMMARUS ####
#Updated code used for the 2024 publication

rm(list=ls(all=TRUE))
library(tidyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(Hmisc)

setwd('C:/Users/kaalt/Desktop/Lobster Manuscript/Theberge_et_al_review3/Lobster_Plots_r3')

#Could adjust v-notching compliance here
vcomp <- seq(0,1, by = 0.1) # 0.9 (Mazur 2019) ~0.7 (Murphy 2018)
vcompliance <- vcomp[11]

#Adjust this value according to protection of interest:  =0 if no female protection. =1 if toss berried fems back. =4 if v notch fem.
Vmaxlast <- 1  

#### Population simulation parameters ####
Npop <- 100        #initial pop size
Tmax <- 200        #max years
Tfishing <- 100    #when we start fishing
Ngroup <- 4        #Number of separate population groups #Unnotched, males, all notched, all females
Amax <- c(35, 35,0,0)  #40-5. max cohort age across all groups (Sundelof et al 2015)
Vmax <- c(0,1,4,Vmaxlast)  #use this only to calculate total egg production

#Values of interest for every level of fishing
mu_f = seq(0, 0.9, by = 0.1)   #Varying levels of fishing pressure
SPR <- array(dim=c(length(Vmax),length(mu_f)), data=0)
RelN <- array(dim=c(length(Vmax),length(mu_f)), data=0)
Yield <- array(dim=c(length(Vmax),length(mu_f)), data=0)
Sexratio <- array(dim=c(length(Vmax),length(mu_f)), data=0)   
fert <- seq(0, 1, by = 0.1)   #proportion of eggs fertilized
eggsfish <- array(dim=c(length(Vmax),length(mu_f)), data=NA)    #number of eggs

#### Life history parameters ####
initialsize <- c(56, 56,0,0) #mm size at "recruitment" to our population model ((53 mm from 2020 americanus assessment: age 5))
matsize <- c(92.5,92.5,0,0) #mm length at 50% maturation.  (Coleman et al. 2023)
a_mat <- -40.59 #14.96 #parameter for maturity ogive       (Coleman et al. 2023)
b_mat <- 0.43 #-0.16   #parameter for maturity ogive       (Coleman et al. 2023)
prop <- c(0.5, 0.5,0,0)     #initial proportions of each sex at recruitment (birth) 

#TL to CL conversion coefficients (Sordalen et al 2022)
  #female:  CL = 0.3842 * TL – 8.45  #alternative values CL = 0.382 x TL - 12.569 (Pavicic et al. 2021)
  #male: TL = CL = 0.3987 * TL – 10.42  #alternative values CL = 0.405 x TL - 17.167 (Pavicic et al. 2021)
inter <- c(-8.45,-10.42)  #c(-6.84,-10.91)
slop <- c(0.3842,0.3987)  #c(0.38,0.40)

#regulations
fsmin <- slop[1] * 250 + inter[1] #(was 86.9, now 87.6)
fsmax <- slop[1] * 320 + inter[1] #(was 114.76, now 114.494)
msmin <- slop[2] * 250 + inter[2] #(was 88.6, now 89.255)
msmax <- slop[2] * 320 + inter[2] #(was 117.09, now 117.164)

#Weight - Length (Sordalen et al 2022) #need to do for CL not TL
#a_w <- c( -12.00057, -1.30364)
#b_w <- c(3.28428, 0.25540)

#natural mortality by sex and size
Mfem <- 250*slop[1]+inter[1]  #250 mm TL is cut off size for size based M.
Mmal <- 250*slop[2]+inter[2]
sizecutoff <- c(Mfem,Mmal)  #in mm CL
#Fernandez-Chacon et al. 2021 - Biological Conservation Estimate of lobster survival - from three reserves at T2 (after protection, so can estimate M w/o fishing).
Flodevigen <- c(0.732751975, 0.846274514, 0.573599506, 0.657874245) #big female, small female, big male, small male
Bolaerne <- c(0.680054367, 0.668539401, 0.473879818, 0.643383665)
Hvaler <- c(0.737909049, 0.840175394, 0.541784944, 0.748808103)
MPAs <- rbind(Flodevigen, Bolaerne, Hvaler)
survival <- colMeans(MPAs) # survival of females > 250 mm total length, females < 250 mm total length, males > 250 mm total length, males < 250 males mm total length
natmort <- -log(survival)  # Annual survival = exp(-M)   #ln(Annual Survival) = - M

natmortsm <- c(natmort[2],natmort[4]) #sm female, male
natmortlg <- c(natmort[1],natmort[3]) #lg female, male

#assuming recruitment follows what was used in Bannister and Addison 1986 since there is no better info available    ####
alpha <- 100000
beta <- 1
kappa <- 0.851  

#### Set up vectors to store age-specific life history characters and fishing probability ####
# cohort "age", time, group #

W <- array(dim=c(Amax[1],Ngroup), data=NA)        #individual weight-at-age
L <- array(dim=c(Amax[1],Ngroup), data=NA)        #length-at-age #switched from NA to 0
mu <- array(dim=c(Amax[1],Ngroup), data=0)        #mortality-at-age (this is constant in initial case)
select <- array(dim=c(Amax[1],Ngroup), data=NA)   #probability of being caught at each age, given size-at-age
Fishing <- array(dim=c(Amax[1],Ngroup), data=0)   #Fishing mortality at age
M <- array(dim=c(Amax[1],Ngroup), data=NA)        #natural mortality by age, sex #Fernandez-Chacon 2020

N <- array(dim=c(Amax[1],Tmax,Ngroup), data=0)    #population numbers in each year class over time
N_agef <- array(dim=c(Amax[1],length(mu_f),length(Vmax)), data=0)  #subset of population for age structure plots
N_agem <- array(dim=c(Amax[1],length(mu_f),length(Vmax)), data=0)
E <- array(dim=c(1,Tmax), data=0)                 #number of eggs at each time
P <- array(dim=c(1,Tmax), data=0)                 #number of larvae at each time
Catch <- array(dim=c(Amax[1],Tmax,Ngroup), data=0) #Catch number or biomass, given selectivity
pmat <- array(dim=c(Amax[1],Ngroup), data=0)       #probability of maturation each year, note it is a matrix so could change with population density if we want to consider fisheries-induced changes in maturation over time

pmolt <- array(dim=c(Amax[1],Ngroup), data=0)      #probability of molting at a given length each year #switched from NA to 0
mi <- array(dim=c(Amax[1],Ngroup), data=0)         #molt increment at length
spnfreq <- 0.5                                      #female annual spawning freq by size (Sordalen)
eggers <- array(dim=c(Amax[1],Tmax), data=0)        #all females that are eligible to be v notched 

#### MOLTING FREQUENCY BY LENGTH ####
#from Sordalen et al 2022
#females
moltclf <- seq(from=56,to=143, by=1)
moltpf <- c(0.954437536,0.951874652,0.949055923,0.945927518,0.942784535,
            0.939530258,0.936100065,0.932807186,0.92935697,0.92538549,
            0.921234435,0.916775608,0.911843396,0.906923762,0.902283444,
            0.896988184,0.891412574,0.886096953,0.879958055,0.873675353,
            0.866951426,0.859425173,0.852112224,0.84520645,0.83731222,
            0.829186494,0.821435598,0.812478453,0.802718743,0.793064894,
            0.784099387,0.774080974,0.763527978,0.753681431,0.742565172,
            0.731202672,0.720503912,0.707339922,0.69384065,0.682299597,
            0.669446448,0.656278599,0.644048897,0.630305036,0.616653674,
            0.602488109,0.587108185,0.5726753,0.559467477,0.546174765,
            0.531626731,0.516898653,0.503461685,0.488632416,0.472730276,
            0.45791236,0.443497331,0.42889552,0.416457773,0.401268319,
            0.38625381,0.374399442,0.359800654, 0.345440052,0.334156281,
            0.320325453,0.306783159,0.296190702,0.283264728,0.270663335,
            0.260848703,0.251267347,0.239636408,0.22835544,0.219612972,
            0.209039508,0.198820958,0.190929407,0.181416957,0.172253565,
            0.165199068,0.156721073,0.148577818,0.142326205,0.134833107,
            0.127654331,0.122156777,0.118162142)
moltfem<- cbind(moltclf,moltpf)

#males
moltclm <- seq(from=56,to=142, by=1)
moltpm <- c(0.998031793,0.997853824,0.997674756,0.997458729,0.997230019,
            0.996999087,0.996748869,0.996449909,0.996127833,0.995805314,
            0.995418397,0.995004883,0.994589388,0.994088931,0.99355836,
            0.993023449,0.992376624,0.991629659,0.990855772,0.99003498,
            0.989121857,0.988222255,0.987163733,0.985995874,0.984841857,
            0.983593863,0.982140663,0.98051032,0.978913532,0.977047655,
            0.974972418,0.972934678,0.970545759,0.96791238,0.965320566,
            0.962273078,0.958944566,0.955661957,0.95179224,0.947606036,
            0.943470757,0.939028606,0.933856402,0.928179295,0.922634705,
            0.916169568,0.909149658,0.902291393,0.894288662,0.885698819,
            0.877310038,0.867522502,0.857145936,0.847024491,0.835228116,
            0.822886978,0.810874752,0.796904772,0.782497381,0.768517259,
            0.753918104,0.737259109,0.719949946,0.703470663,0.68475285,
            0.665626889,0.649639633,0.626155551,0.604382422,0.586642373,
            0.564143359,0.541395128,0.52305956,0.500058221,0.477056652,
            0.458720535,0.440495522,0.417954552,0.395728968,0.378257643,
            0.356877946,0.336017375,0.319787075,0.30012358,0.281123685,
            0.266480619,0.255806675)
moltmal<- cbind(moltclm,moltpm)

#### MOLT INCREMENT ####   
#from Sordalen et al 2022
mi_int <- c(13.31382,  9.52874)
mi_slo <- c(-0.05861, 0.02223)


###########################
#### Life History LOOP ####
###########################

for(g in 1:2) { 
  #for both sexes, where g =1 is unnotched female and g = 2 is male
  #N[cohort "age" (yr), time (yr), group (sex)]
  
  N[1:Amax[1],1,g] = Npop*prop[g]     #initial population size    
  L[1,g] = initialsize[g]      #initial individual size when entering the population model 
  
  #size-specific annual growth: molt probability and molt increment
  for (a in 1:(Amax[1]-1)) {   #for each cohort
    if (a <= Amax[1]) {          
      
      #GROWTH
      #proportion molting given premolt size and sex  #Sordalen manuscript
      #females
      if(g==1){
        for (len in 1:88){
          if((moltfem[len,1]) == ceiling(L[a,1])){
            pmolt[a,1] <- moltfem[len,2]
          }}}
      #males
      if(g==2){
        for(len in 1:87){
          if((moltmal[len,1] ) == ceiling(L[a,2])){
            #print(L[a, g ])
            #print(moltmal[len, 1])
            pmolt[a,2] <- moltmal[len,2] 
          }}}
      
      #Molt increment & final growth
      mi[a,g] <- L[a,g]*mi_slo[g]+mi_int[g]
      
      L[a+1,g] <- L[a,g] + (pmolt[a,g]*mi[a,g]) #Fogarty 1995 Factor book ch 6
      #print(pmolt[a,g])
      #print(mi[a, g])
      
      #WEIGHT
      #W[a,g]=  a_w[g]*L[a,g]^b_w[g]  
      
      #MATURATION   
      pmat[a,g]=1/(1+exp(-(a_mat+b_mat*L[a,g])))   #(Walker 2005, Coleman et al. 2023)

      #FISHERY SELECTIVITY (Sordalen 2019)
      if(g==1){
        if (L[a, g] > fsmin & L[a, g] < fsmax) {        #Females  
          select[a,g] = 1
        } else select[a,g] = 0
      }
      if(g==2){
        if (L[a, g] > msmin & L[a, g] < msmax) {         #Males 
          select[a,g] = 1   
        }else select[a,g] = 0
      }
      
      #NATURAL MORTALITY    (Fernandez-Chacon et al. 2020)
      if (L[a,g] <= sizecutoff[g]){
        M[a,g] <- natmortsm[g]
      }
      if (L[a,g]> sizecutoff[g]){
        M[a,g] <- natmortlg[g]
      }
    } #a if loop
  } #a for loop
} #g

#define female egg production as a function of her length 
eggs <- (468.17*L[, 1])-33004 #(Agnalt 2008)
eggs[1:2] <- 0 #otherwise this is a negative number  
#note these coefficients are for mass-at-age, can be adjusted based on what kind of data are available. original: #estimate of mass-specific egg production (in grams)

########################################
#### Population Dynamics Simulation ####
########################################

###Simulate population dynamics, start from an arbitrary population size and let the population reach a stable age distribution, then start fishing. The population will reach a new, fished, steady state (stable age dist).
## *Note: For the base case, recruitment is based only on Female abundance, and assuming 50:50 offspring sex ratio. We can change this to look at fertilization limitation or sex ratio variation
for(v in (1:length(Vmax))){
  notched <- array(dim=c(Amax[1],Tmax,4), data=0)  #this year's notched females by compliance; there are 4 stacked matrices to show the notched lobs slowing losing their notch
  lost <- array(dim=c(Amax[1],Tmax), data=NA)         #the notched females that lose their notch and return to the unnotched popn.
  
for (fish in 1:length(mu_f)) {     #loop over all levels of fishing mortality from (mu_f)
  for(t in 1:(Tmax-1)) {
    
    #SPAWNING
    N[1,t,4] <- N[1,t,1]
    E[t]=sum(N[, t, 4]*pmat[, 1]*eggs*spnfreq[1] )         #assuming spawning depends ONLY on mature females
    P[t]= E[t]*fert[10]                                    #assumes all eggs are fertilized (for now) #ALTER VALUE FOR SPERM/FERTILIZATION LIMITATION HERE
    
    for(g in 1:Ngroup) {                                    #here g is sex, where g = 1 is unnotched female, and g = 2 is male  
      
      N[1,t+1,1]=(alpha*P[t])/(1+(P[t]/kappa)^beta)*prop[1] # Bannister and Addison 1986
      N[1,t+1,2]=(alpha*P[t])/(1+(P[t]/kappa)^beta)*prop[2] # for unnotched females and then males
      
      #FISHING MORTALITY, given selectivity function

      for (age in 1:(Amax[1]-1)) {
        
        if (Tfishing < t) {
          
          Fishing[age,g] = select[age,g]*mu_f[fish]
          Fishing[Amax[g], g] = select[Amax[g], g]*mu_f[fish] 
          
          # BERRIED FEMALE PROTECTIONS #
          # Vmax=0: none; Vmax=1: one year; Vmax=4 is how long v-notch lasts.
          
          eggers[age,t] <- N[age,t,1]*pmat[age,1]*spnfreq[1]    #all berried females
          
          if (Vmax[v] == 0){ 
            N[age,t,3] <- 0
            N[age,t,4] <- N[age,t,1]
          }
          if (Vmax[v] == 1){      #assuming only female that are caught were thrown back and not fished
            notched[age,t,1] <- eggers[age,t]*select[age,1]*mu_f[fish]  
            N[age,t,3] <- notched[age,t,1]    #nat mort is included in survival eqs
          }
          if (Vmax[v] == 4){      #assuming only females that are caught (function of selectivity & fishing pressure) can be notched
            notched[age,t,1] <- eggers[age,t]*vcompliance*select[age,1]*mu_f[fish]  #from 0-1 by 0.1; so 1:11 . this year's notched females with compliance of 0.90
            N[age,t,3] <- notched[age,t,1]+notched[age,t,2]+notched[age,t,3]+notched[age,t,4] #nat mort is included in survival eqs
          }
        } else {
          Fishing[age,g] = 0
          N[age,t,3] <- 0
          N[age,t,4] <- N[age,t,1]
        } 
        
        #CATCH AND SURVIVAL
        
        Catch[,t+1,g] <- N[,t,g]*(1-exp(-Fishing[,g]))  #Note this is catch numbers, not biomass, assumes natural mortality occurs later in the year than fishing
        #Catch of all the ages, the number of individuals in every age class * fishing mortality of every age class (Fishing=Fishing mortality)
        #Survival to next year
        N[age+1,t+1,2] <- N[age,t,2]*exp(-M[age,2]-Fishing[age,2]) #surviving MALES in each group enter the next age class in the following year, all fish get to spawn before mortality 		 
        
        if(Vmax[v] == 0){    
          N[age+1,t+1,1] <-N[age,t,1] * exp(-M[age,1]-Fishing[age,1]) 
          N[age+1,t+1,4] <- N[age+1,t+1,1] 
        }
        if(Vmax[v] == 1){    #Protect eggers without a notch, so add them back into total female popn without Fishing mortality.
          N[age+1,t+1,1] <- ((N[age,t,1] - (notched[age,t,1])) * exp(-M[age,1]-Fishing[age,1])) + (N[age,t,3]* exp(-M[age,1]))  
          N[age+1,t+1,4] <- N[age+1,t+1,1]  #Calc these to get next year's egg production
        }
        if(Vmax[v] == 4){       #Move lobs through notched matrix until v-notch wears out, and they join unnotched popn
          notched[age+1,t+1,2] <- notched[age,t,1]*exp(-M[age,1]) 
          notched[age+1,t+1,3] <- notched[age,t,2]*exp(-M[age,1])  
          notched[age+1,t+1,4] <- notched[age,t,3]*exp(-M[age,1])  
          lost[age,t] <- notched[age,t,4]        #when leave notched matrix; will return to popn in survival equation.
          #unnotched + ones that lost v-notch can now be fished
          N[age+1,t+1,1] <- ((N[age,t,1] - (notched[age,t,1])) * exp(-M[age,1]-Fishing[age,1])) + (lost[age,t]*exp(-M[age,1]-Fishing[age,1]))  
          N[age+1,t+1,3] <- (N[age,t,3]-lost[age,t])*exp(-M[age,1])
          N[age+1,t+1,4] <- N[age+1,t+1,1] + N[age+1,t+1,3] #Calc these to get next year's egg production
        }
      } #age
    } #g
  } #t
  
  #REFERENCE POINTS AND POPULATION CALCULATIONS OVER FISHING PRESSURES
  RelN[v, fish] <- sum(N[-c(1:4) ,Tfishing+50, c(2,4)])/sum(N[-c(1:4) ,Tfishing-10,c(2,4)])  
  SPR[v,fish] <- E[Tfishing+50]/E[Tfishing-10]                   #cohort lifetime production
  Yield[v,fish] <- sum(N[-35,Tfishing+50,c(2,4)]*(1-exp(-Fishing[-35,-c(3:4)])))    #All fish, all sexes #All fish, all sexes
  Sexratio[v,fish] <- sum(N[, Tfishing+50, 4]*pmat[,1])/(sum(N[, Tfishing+50, 2]*pmat[,2])+sum(N[, Tfishing+50, 4]*pmat[,1]))  #percentage of MATURE females in MATURE popn
  eggsfish[v,fish] <- E[Tfishing+50]  #calculate egg production for each protection scenario here
  N_agef[,fish,v] <- N[,Tfishing+50,4]  #for age structure plots, females
  N_agem[,fish,v] <- N[,Tfishing+50,2]  #for age structure plots, males
  } #end fishing mortality loop
} #end Vmax loop

## POPULATION OVER TIME PLOTS ##


#pal <- c("#000000","#969696")
#tripal <- c("#000000","#969696","#BDBDBD")

colors1 <- c("#CC6677","#88CCEE","#DDCC77","#117733","#332288","#AA4499",
                      "#44AA99","#999933","#882255","#661100","#6699CC","#888888")

colors2 <- c("#F7CAC9","#92A8D1","#F7786B","#034F84")

colors3 <- c("#FF9146","#6ECDDC","#414D55","#FFB9AA","#96A0AF","#FF645A","#646E7D","#007D91","#6ECDDC")

pal <- c(colors2[3],colors2[4])    #f,m
tripal <- colors3[1:3] #protections

plot(colSums(N[,1:Tmax-1,4]), type="l", ylim=c(0,2e+5),col=pal[1],las=1,lwd=2,xlab="Time (years)",ylab="")
lines(colSums(N[,,2]), type="l", lty=2,las=1,lwd=2,col=pal[2])
title(ylab="Relative Population Size", line=4, cex.lab=1.2)
legend("topright",                    
       legend = c("Females", "Males"),
       col = pal,lwd=2,
       lty=c(1,2))

plot(colSums(N[,,3]), type="l", ylim=c(0,1e+5),col=pal[1],las=1,lwd=2,xlab="Time (years)",ylab="")
plot(colSums(N[,,1]), type="l", ylim=c(0,2e+5),col=pal[1],las=1,lwd=2,xlab="Time (years)",ylab="")


###############
#### PLOTS ####
###############

## GROWTH ##

ageplot <- c(5:39) #add ages to plots since recruit to the model at age ~5
Lbind <- cbind(ageplot,L)
Wbind <- cbind(ageplot,W)
mibind <- cbind(ageplot,mi)
pmoltbind <- cbind(ageplot,pmolt)
pmatbind <- cbind(ageplot,pmat)

#Individual growth L
pdf("gam_growth.pdf",   width=5, height=5,  )
matplot(Lbind[,2:3],x=Lbind[,1], xlim = c(0,40), ylim = c(0,180),type="l",lwd=2, las=1, xlab="Age (years)", ylab="Carapace Length (mm)", col=pal)
minor.tick(nx = 2, ny = 5,  
           tick.ratio = 0.5) 
legend("bottomright",                    
       legend = c("Females", "Males"),
       col = pal,lwd=2,
       lty=c(1,2))
dev.off()

#Molt increment
mibind <- mibind[-35,]
pdf("gam_mi.pdf",   width=5, height=5,  )
matplot(mibind[1:Amax[1]-1,2:3],x=mibind[,1], xlim = c(0,40), ylim = c(0,16), xlab="Age (years)",ylab="Molt Increment (mm)",type="l",lwd=2, las=1, col=pal)
legend("bottomleft",                    
       legend = c("Females", "Males"),
       col = pal,lwd=2,
       lty=c(1,2))
dev.off()

#proportion molting
pdf("gam_propmolt.pdf",   width=5, height=5,  )
matplot(x=pmoltbind[,1],y=pmoltbind[1:Amax[1],2],xlab="Age (years)",ylab="Proportion Molting", type="l",lwd=2, las=1, col=pal, ylim=c(0:1),xlim=c(0,30))
lines(x=pmoltbind[1:11,1],y=pmoltbind[1:11,3],col=pal[2],lty=2,lwd=2)
legend("topright",                    
       legend = c("Females", "Males"),
       col = pal,lwd=2,
       lty=c(1,2))
dev.off()

#Weight-Length relationship
# matplot(Lbind[,2:3],W[,2:3], xlab="Length",ylab="Weight", type="l",lwd=2, las=1, col=pal)
# legend("bottomright",                    
#        legend = c("Females", "Males"),
#        col = pal,lwd=2,
#        lty=c(1,2))

#Probability mature
matplot(pmatbind[,2:3],x=pmatbind[,1], type="l", lwd=2, las=1, ylab ="Probability Mature", xlab = "Age (years)", col=pal, lty=1, xlim = c(0,30))
abline(h=0.5,lty=3)

matplot(pmatbind[,2:3],x=L[,1:2], type="l", lwd=2, las=1, ylab ="Probability Mature", xlab = "Carapace Length", col=pal, lty=1, xlim = c(0,130))
abline(h=0.5,lty=3)
abline(h=96,lty=3)
abline(v=92.5, lty=3)

## REPRODUCTION & SEX RATIO ##

##Function for values over levels of fishing pressure

ratiofig <- function(ylbl,      #eg. y label: Sex Ratio (Females:Males) or Ratio of Unprotected to Protected Females
                     value,    #eg. name of table Sexratio, Vratio, etc.
                     name,     #name of image file "name.pdf"
                     leg     #legend location
) {
  first <- value
  first <- as.data.frame(first[-4,])
  row.names(first)<-c("No Female Protection", "1-year Protection","4-year Protection")
  colnames(first)<-c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")
  first$Scenario <- c("No Female Protection", "1-year Protection","4-year Protection")
  #return(first)
  ratio_melt <- melt(first, id.vars = 'Scenario')
  pdf(name,   width=6, height=4,  )
  rf <-ggplot(ratio_melt, aes(x=variable,y=value,color=Scenario,group=Scenario))+
    theme_bw()+
    geom_line(linetype='dashed')+
    geom_point(stat= 'identity',size=2)+
    scale_color_manual(breaks=c("No Female Protection", "1-year Protection","4-year Protection"),values=tripal)+
    labs(x='Fishing Pressure',y=ylbl)+
    ylim(0,1.01)+
    ggtitle("H. gammarus")+
    theme(plot.title = element_text(color="black", size=14, face="italic"))+
    theme(legend.position = c(0.8, leg),legend.background = element_rect(fill = "white"))
  return(rf)
}
ratiofig(value=SPR,ylbl = "SPR",name="gam_SPR.pdf",leg=0.8)
dev.off()
ratiofig(value=Sexratio, ylbl="Proportion Female",name="gam_SexRatio.pdf",leg=0.8)
dev.off()

##### Egg production plot ####

#by fishing pressure
eggsfish1 <- eggsfish[-4,]
row.names(eggsfish1)<-c("No Female Protection", "1-year Protection","4-year Protection")
colnames(eggsfish1)<-c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")
alleggs <- array(dim=c(3,length(mu_f)), data=NA)
row.names(alleggs)<-c("No Female Protection", "1-year Protection","4-year Protection")
colnames(alleggs)<-c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")
for(f in 1:length(mu_f)){
  alleggs[,f] <- eggsfish1[,f]/(eggsfish1[1,1])
} 
alleggs <- as.data.frame(alleggs)
alleggs$Scenario <- c("No Female Protection", "1-year Protection","4-year Protection")
dategg <- melt(alleggs,id.vars="Scenario")

pdf("gam_eggprod.pdf",   width=6, height=4,  )
ggplot(dategg,aes(x=variable,y=value,color=Scenario,group=Scenario))+
  geom_line(linetype='dashed')+
  geom_point(stat = "identity",size=2)+
  xlab("Fishing Pressure")+
  ylab("Egg Production Proportional to Unfished Production")+
  ylim(0,1)+
  scale_color_manual("Scenario",breaks=c("No Female Protection", "1-year Protection","4-year Protection"),values=tripal)+
  theme_bw()+
  ggtitle("H. gammarus")+
  theme(plot.title = element_text(color="black", size=14, face="italic"))+
  theme(legend.position = c(0.8, 0.8),legend.background = element_rect(fill = "white"))
dev.off()

#### Age structure ####

#Set up function to make plots
agefig <- function(ylbl,      #eg. y label: fished or unfished abundance 
                   fish,     #eg. what level of fishing pressure use 1-10 for 0-0.9
                   v,       #eg. vmax =0,1,4 but use 1,2,3 to select
                   name,     #name of image file "name.pdf"
                   zoom1,
                   zoom2
) {
  age_structure <- cbind(N_agef[,fish,v], N_agem[,fish,v])
  u<- as.data.frame(age_structure)
  names(u) <- c("females", "males")
  u$age<-c(5:(Amax[1]+4))
  age<-c(5:(Amax[1]+4))
  dat1 <- melt(u,id.vars="age")
  
  pdf(name,   width=5, height=5,  )
  ff <- ggplot(dat1) + 
    aes(age) + 
    coord_flip()+  
    theme_bw()+
    xlab("Relative Age (years)")+
    ylab(ylbl)+
    xlim(zoom1,zoom2)+
    geom_bar(data =  dat1[dat1[["variable"]]=="females",], aes(y = -value, fill="females"), stat="identity", fill=pal[1])+ 
    geom_bar(data =  dat1[dat1[["variable"]]=="males",], aes(y = value, fill="males"), stat="identity", fill=pal[2])+
    scale_fill_manual("Sex",breaks=c("Females", "Males"),values=pal)+
    theme_classic()+
    theme(legend.position = c(0.8, 0.8),legend.background = element_rect(fill = "white"),axis.text.x = element_blank(),axis.ticks = element_blank())
  ff
}

# No fishing, F = 0.0
a00 <- agefig(ylbl="Unfished Abundance", fish=1,v=1,name="gam_age_nofish.pdf",zoom1 = 4,zoom2 = 40)    #without fishing
a00
dev.off()

#Full at F = 0.9
a09 <- agefig(ylbl="Fished Abundance at F = 0.9", fish=10,v=1,name="gam_age_v0_f9.pdf",zoom1 = 4,zoom2 = 40)    #with fishing, vmax=0
a09
dev.off()
a19 <- agefig(ylbl="Fished Abundance at F = 0.9", fish=10,v=2,name="gam_age_v1_f9.pdf",zoom1 = 4,zoom2 = 40)    #with fishing, vmax=1
a19
dev.off()
a49 <- agefig(ylbl="Fished Abundance at F = 0.9", fish=10,v=3,name="gam_age_v4_f9.pdf",zoom1 = 4,zoom2 = 40)    #with fishing, vmax=4
a49
dev.off()

#Zoom in on slot limit at F = 0.9
# a09z <- agefig(ylbl="Fished Abundance at F = 0.9 & no female protection", fish=10,v=1,name="gam_age_v0_f9_zoom.pdf",zoom1 = 7,zoom2 = 25)    #with fishing, vmax=0
# a09z
# dev.off()
# a19z <- agefig(ylbl="Fished Abundance at F = 0.9 & 1-year female protection", fish=10,v=2,name="gam_age_v1_f9_zoom.pdf",zoom1 = 7,zoom2 = 25)    #with fishing, vmax=1
# dev.off()
# a49z <- agefig(ylbl="Fished Abundance at F = 0.9 & 4-year female protection", fish=10,v=3,name="gam_age_v4_f9_zoom.pdf",zoom1 = 7,zoom2 = 25)    #with fishing, vmax=4
# dev.off()


#Full at F = 0.3
a03 <- agefig(ylbl="Fished Abundance at F = 0.3", fish=4,v=1,name="gam_age_v0_f3.pdf",zoom1 = 4,zoom2 = 40)    #with fishing, vmax=0
a03
dev.off()
a13 <- agefig(ylbl="Fished Abundance at F = 0.3", fish=4,v=2,name="gam_age_v1_f3.pdf",zoom1 = 4,zoom2 = 40)    #with fishing, vmax=1
a13
dev.off()
a43 <- agefig(ylbl="Fished Abundance at F = 0.3", fish=4,v=3,name="gam_age_v4_f3.pdf",zoom1 = 4,zoom2 = 40)    #with fishing, vmax=4
a43
dev.off()

#Zoom in on slot limit at F = 0.3
# a03z <- agefig(ylbl="Fished Abundance at F = 0.3 & no female protection", fish=4,v=1,name="gam_age_v0_f3_zoom.pdf",zoom1 = 7,zoom2 = 25)    #with fishing, vmax=0
# dev.off()
# a13z <- agefig(ylbl="Fished Abundance at F = 0.3 & 1-year female protection", fish=4,v=2,name="gam_age_v1_f3_zoom.pdf",zoom1 = 7,zoom2 = 25)    #with fishing, vmax=1
# dev.off()
# a43z <- agefig(ylbl="Fished Abundance at F = 0.3 & 4-year female protection", fish=4,v=3,name="gam_age_v4_f3_zoom.pdf",zoom1 = 7,zoom2 = 25)    #with fishing, vmax=4
# dev.off()


# all together
pdf("gam_age_all_f3andf9.pdf",   width=10, height=15,  )
age3and9 <- ggarrange(a00,a00,a03,a09,a13,a19,a43,a49, labels=c("No Fishing","No Fishing","No female protection","No female protection", "1-year female protection", "1-year female protection","4-year female protection","4-year female protection"), 
                      ncol = 2, nrow = 4, common.legend = TRUE, legend = "bottom")
age3and9
annotate_figure(age3and9, top = text_grob("H. gammarus", face = "italic", color = "black", size = 22))
dev.off()

#### Size Structure Plots ####

##Mature individuals only
#Set up function to make plots
sizefig <- function(
  fish,  #eg. 2-10 for 0.1-0.9
  v,      #eg. 1-3 for 0,1,4  
  name     #name of image file "name.pdf"
  
) {
  fishedsteady_f=N_agef[ ,fish,v]
  fishedsteady_m=N_agem[ ,fish,v]
  
  binsize=10
  bins= seq(56, 190, by = binsize)
  Lindex<-matrix(nrow=Amax[1], ncol= 2)
  for (g in 1:2){
    for (i in 1:(length(bins)-1)) {
      for (a in 1:Amax[1]) {
        if(L[a,g] >= bins[i] & L[a,g] < bins[i+1]) {
          Lindex[a,g] <- i 
        }}}}
  sizeind=(Lindex*binsize)+55
  
  fishedstructure_f<-tapply(fishedsteady_f*pmat[, 1], sizeind[,1], sum)
  fishedstructure_m<-tapply(fishedsteady_m*pmat[, 2], sizeind[,2], sum)
   
  fished_b <- as.data.frame(cbind(fishedstructure_f,fishedstructure_m))
  fished_b$Size_Class <- as.numeric(row.names(fished_b)) 
  colnames(fished_b) <- c('Females','Males','Size_Class')
  f<-melt(fished_b,id.vars="Size_Class")
  f$Size_Class <- as.factor(f$Size_Class)
  binnames=c("56-65","66-75","76-85","86-95","96-105","106-115", "116-125",
              "126-135","136-145","146-155","155-165","166-175","176+")
  
  pdf(name,   width=6, height=4,  )
  ggplot(data=f, aes(fill=variable,x=Size_Class,y=value))+
    geom_bar(stat='identity',position = 'dodge')+
    xlab("Size (mm CL)")+
    ylab("Abundance of Mature Individuals")+
    scale_fill_manual(name="Sex",breaks=c("Females","Males"),values=pal)+
    scale_x_discrete(labels=binnames)+
    ylim(0,14500)+ 
    geom_vline(xintercept = 3.5, color="red") +
    geom_vline(xintercept = 7.5, color="red") +
    theme_classic()+
    theme(legend.position = c(0.85, 0.85),legend.background = element_rect(fill = "white"),
          axis.text.x=element_text(angle=45,hjust=1,vjust=1, size=20),
          axis.text.y=element_text(size=18),
          axis.title = element_text(size=20))
}

#F = 0.0
s00 <- sizefig(fish = 1,v=1,name = "gam_size_nofish.pdf") #unfished
s00
dev.off()

#F = 0.1
sizefig(fish = 2,v=1,name="gam_size_v0_f1.pdf") #vmax=0
dev.off()
sizefig(fish = 2, v=2,name="gam_size_v1_f1.pdf") #vmax=1
dev.off()
sizefig(fish = 2, v=3,name="gam_size_v4_f1.pdf") #vmax=4
dev.off()

#F = 0.3
s03 <- sizefig(fish = 4,v=1,name="gam_size_v0_f3.pdf") #vmax=0
s03
dev.off()
s13 <- sizefig(fish = 4, v=2,name="gam_size_v1_f3.pdf") #vmax=1
s13
dev.off()
s43 <- sizefig(fish = 4, v=3,name="gam_size_v4_f3.pdf") #vmax=4
s43
dev.off()

#F = 0.5
sizefig(fish = 6,v=1,name="gam_size_v0_f5.pdf") #vmax=0
dev.off()
sizefig(fish = 6, v=2,name="gam_size_v1_f5.pdf") #vmax=1
dev.off()
sizefig(fish = 6, v=3,name="gam_size_v4_f5.pdf") #vmax=4
dev.off()

#F = 0.7
sizefig(fish = 8,v=1,name="gam_size_v0_f7.pdf") #vmax=0
dev.off()
sizefig(fish = 8, v=2,name="gam_size_v1_f7.pdf") #vmax=1
dev.off()
sizefig(fish = 8, v=3,name="gam_size_v4_f7.pdf") #vmax=4
dev.off()

#F = 0.9
s09 <- sizefig(fish = 10,v=1,name="gam_size_v0_f9.pdf") #vmax=0
s09
dev.off()
s19 <- sizefig(fish = 10, v=2,name="gam_size_v1_f9.pdf") #vmax=1
s19
dev.off()
s49 <- sizefig(fish = 10, v=3,name="gam_size_v4_f9.pdf") #vmax=4
s49
dev.off()

#all in one
pdf("gam_size_f3_and_f9.pdf",  width = 15, height = 15,  )
f3and9 <- ggarrange(s03+theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()),
                    s09+theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
                              axis.text.y = element_blank(), axis.title.y = element_blank()),
                    s13+theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                    s19+theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
                              axis.text.y = element_blank(), axis.title.y = element_blank()),
                    s43+theme(axis.title.y = element_blank()),
                    s49+theme(axis.text.y = element_blank(), axis.title.y = element_blank()), 
                    labels=c("    No female protection","    No female protection", 
                             " 1-year female protection", "1-year female protection",
                             " 4-year female protection","4-year female protection"),
                    ncol=2, nrow=3, hjust=-1.1, align = "v",
                    common.legend = TRUE, legend = "bottom", font.label = list(size=18))
f3and9
annotate_figure(f3and9, top = text_grob("F = 0.3                                                                  F = 0.9",
                                        face = "italic", color = "black", size = 22))
dev.off()


###Size Structure Plot With immature individuals
#Set up function to make plots
sizefig2 <- function(
  fish,  #eg. 2-10 for 0.1-0.9
  v,      #eg. 1-3 for 0,1,4
  name     #name of image file "name.pdf"
) {
  fishedsteady_f=N_agef[ ,fish,v]
  fishedsteady_m=N_agem[ ,fish,v]
  
  binsize=10
  bins= seq(56, 190, by = binsize)
  Lindex<-matrix(nrow=Amax[1], ncol= 2)
  for (g in 1:2){
  for (i in 1:(length(bins)-1)) {
    for (a in 1:Amax[1]) {
      if(L[a,g] >= bins[i] & L[a,g] < bins[i+1]) {
        Lindex[a,g] <- i 
      }}}}
  sizeind=(Lindex*binsize)+55
  
  fishedstructure_f<-tapply(fishedsteady_f*pmat[, 1], sizeind[,1], sum)
  fishedstructure_m<-tapply(fishedsteady_m*pmat[, 2], sizeind[,2], sum)
  
  fished_b <- as.data.frame(cbind(fishedstructure_f,fishedstructure_m))
  fished_b$Size_Class <- as.numeric(row.names(fished_b)) 
  colnames(fished_b) <- c('Mature_Females','Mature_Males','Size_Class')
  f<-melt(fished_b,id.vars="Size_Class")
  f$Size_Class <- as.factor(f$Size_Class)
  
  fishedall_f<-tapply(fishedsteady_f, sizeind[,1], sum)
  fishedall_m<-tapply(fishedsteady_m, sizeind[,2], sum)
  fished_ball <- as.data.frame(cbind(fishedall_f,fishedall_m))
  fished_ball$Size_Class <- as.numeric(row.names(fished_ball)) 
  colnames(fished_ball) <- c('Juvenile_Females','Juvenile_Males','Size_Class')
  fall<-melt(fished_ball,id.vars="Size_Class")
  fall$Size_Class <- as.factor(f$Size_Class)
  binnames=c("56-65","66-75","76-85","86-95","96-105","106-115", "116-125",
             "126-135","136-145","146-155","155-165","166-175","176+")

  pdf(name,   width=6, height=4,  )
  ggplot(NULL, aes())+
    geom_bar(data=fall, aes(fill=variable,x=Size_Class,y=value),stat='identity',position = 'dodge')+
    geom_bar(data=f, aes(fill=variable,x=Size_Class,y=value),stat='identity',position = 'dodge')+
    xlab("Size (mm CL)")+
    ylab("Total Abundance")+
    scale_x_discrete(labels=binnames)+
    ylim(0,8e+04)+
    geom_vline(xintercept = 3.5, color="red") +
    geom_vline(xintercept = 7.5, color="red") +
    theme_classic()+
    #ggtitle("H. gammarus")+
    theme(plot.title = element_text(color="black", size=18, face="italic"))+
    scale_fill_manual(name="Group",labels=c("Juvenile Females","Juvenile Males","Mature Females","Mature Males"),values = colors2)+
    theme(legend.position = c(0.86, 0.77),legend.background = element_rect(fill = "white"), 
          legend.title=element_text(size=20,face="bold"), legend.text=element_text(size=18),
          #axis.text.x=element_text(angle=45,hjust=1,vjust=1, size=18),
          axis.text.x = element_blank(),axis.title.x = element_blank(),
          axis.text.y=element_text(size=18), 
          axis.title = element_text(size=20))
  
}

#Unfished
sim00 <- sizefig2(fish = 1,v=1,name = "gam_size_im_nofish.pdf") 
sim00
dev.off()


#stack with mature only plot
pdf("gam_size_nofish_im_mat.pdf",  width = 15, height = 14,  ) 
nofish_im_mat <- ggarrange(sim00,s00,
                           ncol=1, nrow=2, common.legend = TRUE, legend = "bottom",
                           labels=c("C","D"), hjust=-8, align = "hv", font.label = list(size=20))
nofish_im_mat
annotate_figure(nofish_im_mat, top = text_grob("H. gammarus",
                                               face = "italic", color = "black", size = 26))
dev.off()

###Average Size for Mature Females and Males across all scenarios in one plot
withLf <-array(dim=c(length(mu_f),3), data=NA)
withLm <-array(dim=c(length(mu_f),3), data=NA)

for (vi in 1:3){
  for (fi in 1:10){
    withLf[fi,vi] <- sum((N_agef[,fi,vi]*pmat[, 1])*L[,1])/(sum(N_agef[,fi,vi]*pmat[, 1]))
    withLm[fi,vi] <- sum((N_agem[,fi,vi]*pmat[, 1])*L[,1])/(sum(N_agem[,fi,vi]*pmat[, 1]))
  }
}

colnames(withLf)<-c("No Female Protection", "1-year Protection","4-year Protection")
row.names(withLf)<-c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")
colnames(withLm)<-c("No Female Protection", "1-year Protection","4-year Protection")
row.names(withLm)<-c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")
withLf<-as.data.frame(withLf)
withLm<-as.data.frame(withLm)
withLf$Sex<-"Females"
withLm$Sex<-"Males"

withboth <- rbind(withLf,withLm)
wb<-melt(withboth,id.vars="Sex")
wb$Fishing <- rep(c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),6)

sizefig3<-function(v,name){
  pdf(name,   width=6, height=4,  )
  ggplot(wb[v,],aes(color=Sex,group=Sex,x=Fishing,y=value))+
    geom_line(linetype = "dashed")+
    geom_point(size=2)+
    xlab("Fishing Pressure")+
    ylab("Mean Length (mm CL) of Mature Individuals")+
    scale_color_manual(name="Sex",breaks=c("Females","Males"),values=pal)+
    theme_bw()+
    ggtitle("H. gammarus")+
    theme(plot.title = element_text(color="black", size=14, face="italic"))+
    theme(legend.position = c(0.85, 0.85),legend.background = element_rect(fill = "white"))
}

sizefig3(v=1:20,name = "gam_size_0avg_disparity.pdf") #no female protections
dev.off()
sizefig3(v=21:40,name = "gam_size_1avg_disparity.pdf") #vmax=1
dev.off()
sizefig3(v=41:60,name = "gam_size_4avg_disparity.pdf") #vmax=4
dev.off()


#Average difference in size females-males

dif <- withLf[,1:3]-withLm[,1:3]
dif <- as.data.frame(dif)
dif$Fishing <- c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")
datdif <- melt(dif,id.vars="Fishing")

pdf("gam_size_diff.pdf",   width=6, height=4,  )
ggplot(datdif,aes(x=Fishing,y=value,color=variable,group=variable))+
  geom_point(stat = "identity",size=0)+
  geom_point(size=2)+
  geom_line(linetype='dashed')+
  geom_hline(yintercept=0,color="red")+
  xlab("Fishing Pressure")+
  ylab("Mean Length Difference (mm CL)")+
  scale_color_manual("Scenario",breaks=c("No Female Protection", "1-year Protection","4-year Protection"),values=tripal)+
  theme_bw()+
  ggtitle("H. gammarus")+
  theme(plot.title = element_text(color="black", size=14, face="italic"))+
  theme(legend.position = c(0.17, 0.3),legend.background = element_rect(fill = "white"))
dev.off()

per <- withLm[,1:3]/withLf[,1:3]
per <- as.data.frame(per)
per$Fishing <- c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")
datper <- melt(per,id.vars="Fishing")

#END

