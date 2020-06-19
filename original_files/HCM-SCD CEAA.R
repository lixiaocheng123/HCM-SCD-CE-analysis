# Bayesian Models for Cost-Effectiveness Analysis
# Loads packages

 S <- 4 # Number of health states
        # s = 1 Healthy health state
        # s = 2 Stroke HCM-Related Health state
        # s = 3 SCD(Sudden Cardiac Death) Health state
        # s = 4 DAC (Death All Causes) Health state
        # EXISTING METHOD OF SCD-RISK-PREDICTION (EMSRP):Current Method of SCD Risk Prediction
        # HCM - SCD RISK PREDICTION MODEL(HSRPM)
 
 J <-10 # Number of years of follow up
 
 # Now load the observered data on transitions among the states for the two treatments
 Model.file="HCMmodel.txt"              # Specifies file with the model defining observed data
 inits1=source("hcminits1.txt")$value   # Loads the initial values for the first chain
 inits2=source("hcminits2.txt")$value   # Loads the initial values of the second chain
 inits=list(inits1,inits2)              # Combines into a list files with inital values
 dataBugs=source("HCMdata.txt")$value   # Loads observed data
 
 # Now run the MCMC
 library(R2OpenBUGS)
 params=c("lambda.0","lambda.1")        # Defines parameters to save
 n.iter <- 10000
 n.burnin<- 5000
 n.thin <- floor((n.iter- n.burnin)/500)
 n.chain=2
 debug=FALSE
 mm1 <- bugs(data=dataBugs,inits=inits,parameters.to.save = params,
             model.file = Model.file, n.chains=n.chain,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC = TRUE)
 print(mm1,digits = 3)
 attach.bugs(mm1)
# Analyze MCMC Convergence with CODA
# Check if all entries in Rhat component of Bugs output are less than 1.1
# all(mm1$summary[,"Rhat"]< 1.1)

# Convert bugs output for coda and create an mcmc object
 hcm<- bugs(data=dataBugs,inits=inits,parameters.to.save = params,
            model.file = Model.file,codaPkg=TRUE, n.chains=n.chain,n.iter=n.iter,n.burnin=n.burnin,n.thin=n.thin,DIC = TRUE)
 hcm.coda<-read.bugs(hcm)
 
 # Now we run the Markov model from R
 start <- c(1000,0,0,0)                 # Analysis for virtual cohort of 1000, individuals
                                        # NB All patients enter the model from the first state "Healthy"
 
 #Determine the Markov transitions
 m.0 <- m.1 <- array(NA,c(n.sims,S,(J+1)))
 for(s in 1:S){
        m.0[,s,1] <- start[s]
        m.1[,s,1] <- start[s]
 }
 
 #NB
 #     BUGS only outputs matrices for lambda.0 and lambda.1 with simulations for the "random" part
 #     ie only the first 2 rows, as the last two are deterministically defined as c(0,0,1,1)
 #     because once a patient is in SCD,and DAC, they can't move away. So there is the need to 
 #     reconstruct a full matrix with S rows and S columns for each MCMC simulations. This is done by 
 #     defining  new arrays lamda0 and lamda1 and then stacking up the simulated values for the first (S-2)
 #     rows saved in lambda.0[i,,] and lambda.1[i,,] for MCMC simulations i with a row vector
 #     containing (S-2) 0s and then two 1's, ie c(0,0, 1,1)

 lamda0=lamda1=array(NA, c(n.sims,S,S))
 for (i in 1:n.sims) {
         lamda0[i,,]=rbind(rbind(lambda.0[i,,],c(0,0,1,1)),c(0,0,1,1))
         lamda1[i,,]=rbind(rbind(lambda.1[i,,],c(0,0,1,1)),c(0,0,1,1))
         for (j in 2: (J+1)) {
                 for (s in 1:S) {
                         # Now use lamda0,and lamda1, for the matrix multiplication
                         m.0[i,s,j] <- sum(m.0[i,,j-1]*lamda0[i,,s])
                         m.1[i,s,j] <- sum(m.1[i,,j-1]*lamda1[i,,s])
                 }
                 
         }
         
 }
 # Now we draw barplot of the number of people in each state at each time point during follow up
 
 par(mfrow=c(1,2))
 barplot(apply(m.0,c(2,3),sum),names.arg=seq(0,10),space=.2,xlab="Virtual follow up",
         ylab="Proportion of patients in each state",main="EMSRP",ylim = c(0,200000000))
 
 barplot(apply(m.1,c(2,3),sum),names.arg=seq(0,10),space=.2,xlab="Virtual follow up",
         ylab="Proportion of patients in each state",main="HSRPM",ylim = c(0,200000000)) 
 # Plot trace and density for all mcmc chains
plot(hcm.coda, trace = TRUE, density = TRUE, smooth = FALSE)
 
 # Run the economic analysis
 #Now we define the benefits
 utility.score.0<-utility.score.1<-array(NA,c(n.sims,2,J)) # Defines measures of accumulated untility
                                                           # under each treatment
 dec.rate<-0.35                    # Defines utility decreament rate to apply when a non-fatal HCM event occurs at j >0
 utility.score.0[,1,]<-rep(0.637,J) # Utility for occupying the state "Healthy" under treatment t=0
 utility.score.1[,1,]<-rep(0.637,J) # Utility for occupyingthe state "Healty" under treatment t= 1
 utility.score.0[,2,1]<-dec.rate*utility.score.0[1,1,1] # Utility for occupying state "Stroke-HCM Related" under treatment t=0
 utility.score.1[,2,1]<-dec.rate*utility.score.1[1,1,1] # Utility for occupying state "Stroke-HCM Related" under treatment t=1
  for (i in 1:n.sims) {
         for (j in 2:J) {
                 utility.score.0[i,2,j]<-utility.score.0[i,2,j-1]*(1- dec.rate)
                 utility.score.1[i,2,j]<-utility.score.1[i,2,j-1]*(1- dec.rate)
         }
         
 }
 # We now compute QALY's accumulated under each treatment for each year of follow up
 Qal.0<-Qal.1<-matrix(NA,n.sims,J)
 for (i in 1:n.sims) {
         for (j in 1:J) {
                 Qal.0[i,j]<-(m.0[i,1,j]%*%utility.score.0[i,1,j] 
                                + m.0[i,2,j]%*%utility.score.0[i,2,j])/m.0[1,1,1]
                 Qal.1[i,j]<-(m.1[i,1,j]%*%utility.score.1[i,1,j] 
                                + m.1[i,2,j]%*%utility.score.1[i,2,j])/m.1[1,1,1]
                 
         }
         
 }
# Now sum values across all time points, and create matrix effectiveness
 eff<-array(NA,c(n.sims,2,J))
 eff[,1,]<-apply(Qal.0, 1,sum)
 eff[,2,]<-apply(Qal.1, 1,sum)
 
 # We define the annual cost for each non-fatal health state under each treatment
 unit.cost.0 <-c(4792,22880)
 unit.cost.1 <-c(4812,22880) 
 
 #Create a holding cost variable to track yearly (j>0)accumulated cost under each treatment
 cost.0<-cost.1<-matrix(NA,n.sims,J)
 for (i in 1: n.sims) {
         for (j in 2:(J+1)) {
                 cost.0[i,j-1]<-(m.0[i,S,j]+ m.0[i,(S-1),j])*(unit.cost.0%*%m.0[i,1:(S-2),j])/sum(m.0[i,1:(S-2),j])
                 + unit.cost.0%*%m.0[i,1:(S-2),j]
                 cost.1[i,j-1]<-(m.0[i,S,j]+ m.0[i,(S-1),j])*(unit.cost.1%*%m.0[i,1:(S-2),j])/sum(m.0[i,1:(S-2),j])
                 + unit.cost.1%*%m.0[i,1:(S-2),j]
                 
         }
         
 }
 
 
# We now derive a general formulation to apply discount to cost and benefits
 rate.b <- 0.035	# discount rate for benefits (3.5%)
 rate.c <- 0.035	# discount rate for costs (3.5%)
 # Defines the discount factors
 disc.b <- numeric(); disc.c <- numeric()
 disc.b[1] <- 1; disc.c[1] <- 1
 for (j in 2:J) {
         disc.b[j] <- (1+rate.b)^(j-1)
         disc.c[j] <- (1+rate.c)^(j-1)
 }
 disc.cost.0 <- disc.eff.0 <- disc.cost.1 <- disc.eff.1 <- matrix(NA,n.sims,J)
 for (j in 1:J) {
         disc.cost.0[,j] <- cost.0[,j]/disc.c[j]
         disc.cost.1[,j] <- cost.1[,j]/disc.c[j]
         disc.eff.0[,j] <- eff[,1,j]/disc.b[j]
         disc.eff.1[,j] <- eff[,2,j]/disc.b[j]
 }
 
 # Now sum the values across all time points and create matrix of costs
 c <- matrix(NA,n.sims,2)
 c[,1] <- apply(disc.cost.0,1,sum)
 c[,2] <- apply(disc.cost.1,1,sum)
 
 # Sum all discounted valees of effectiveness and create a matrix of discounted effectiveness
 e <- matrix(NA,n.sims,2)
 e[,1] <- apply(disc.eff.0,1,sum)
 e[,2] <- apply(disc.eff.1,1,sum)
 
 
 # Cost-effectiveness analysis
 library(BCEA)
 ints <- c("EMSRP","HSRPM")
 m <- bcea(e,c,ref=2,interventions=ints,Kmax=25000)
 
 contour2(m,25000)
 plot(m)
 summary(m)
 