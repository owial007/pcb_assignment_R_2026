## Below is the code for the assignment
## Each section code is given under individual bracket
## The ipynb notebook for the code is available at:-
#Reading the files
#abundance file
print("Abundance File:")
abund<-t(read.csv("abundance.csv",row.names=1,header =1))
View(head(abund))

#metadata file
print("Metadata File:")
metadata<-read.csv("metadata.csv",row.names=1,header =1)
View(head(metadata))

#taxonomy file
print("Taxonomy File:")
tax<-read.csv("taxonomy.csv",row.names=1,header =1)
View(head(tax))

#Part A: Species Prevalence and Abundance
#A1) The 5 most prevelant species are:
#Checking for null alues in the abundance matrix
cat(sum(is.na(abund)))
#extracting the prevalant columns that have abundance more than 0
prev_num<-colSums(abund>0)
#Calculating the relative prevelance of each prevalent column
prev_perc<-colSums(abund>0)/nrow(abund) * 100

#The top 5 species based on count
top_5_prev<-names(prev_num[order(-prev_num,names(prev_num))][1:5])
cat("Top 5 Prevalent species based on Count:\n",top_5_prev)

#the top 5 species based on percentage
top_5_prev<-names(prev_perc[order(-prev_perc,names(prev_perc))][1:5])
cat("\nTop 5 Prevalent species based on Percentage:\n",top_5_prev)
#They both are same

#Printing the genus and speices of the top 5 most relevent features
cat("\n\n")
View(tax[top_5_prev,'species'])


#A2) The relative Abundances of these 5 species

#calculating the total sum for each sample
sample_sum<-rowSums(abund)

#Calculating the relative abundance matrix 
rel_abundance<-sweep(abund,1,sample_sum,"/")*100

#Extracting the reletive abunda of the top5 species
top_5_rel_abundance <-rel_abundance[,top_5_prev]
#View(head(top_5_rel_abundance))

#Calculating the mean and median for each of the top 5 prevalent species
top_5_rel_abd_sum<-data.frame(
  species=colnames(top_5_rel_abundance),
  mean_rel_abund=colMeans(top_5_rel_abundance),
  median_rel_abund=apply(top_5_rel_abundance,2,median)
)
rownames(top_5_rel_abd_sum)<-NULL
View(t(t(top_5_rel_abd_sum)))

#Part B: Descriptive statistics and group comparison
#Checking the null values in the metadata
cat(sum(is.na(metadata)))

#matching the rownames of metadata and relative abundance
metadata<-metadata[match(rownames(top_5_rel_abundance), rownames(metadata)), ]

#B1: Extracting the mean,median and IQR for each species in CRC
crc<-metadata$disease=='CRC'
healthy<-metadata$disease=='healthy'
#extracting CRFC and healthy rows
top_5_crc<-top_5_rel_abundance[crc,]
top_5_healthy<-top_5_rel_abundance[healthy,]
#defining the metric cal function
sum_grp<-function(data){
  data.frame(
    species=colnames(data),
    mean = colMeans(data),
    median = apply(data,2,median),
    IQR = apply(data,2,IQR),
    std = apply(data,2,sd),
    row.names=NULL
  )
}

#Extracting values for CRC
crc_sum<-sum_grp(top_5_crc)
View(t(t(crc_sum)))

hel_sum<-sum_grp(top_5_healthy)
View(t(t(hel_sum)))

#B2: Getting the difference between the median for the two groups
diff_pp<-crc_sum$median-hel_sum$median
diff_abs<-abs(crc_sum$median-hel_sum$median)
diff_category<-ifelse(diff_pp>0,"CRC",ifelse(diff_pp<0,"Healthy","Equal"))

b2_df<-data.frame(species=crc_sum$species,
                  difference=diff_pp,
                  Absoulte=diff_abs,
                  Category=diff_category)
print("For Median Difference")
View(t(t(b2_df)))

#B3: Getting the difference between the standard dev for the two groups
diff_pp<-crc_sum$std-hel_sum$std
diff_abs<-abs(crc_sum$std-hel_sum$std)
diff_category<-ifelse(diff_pp>0,"CRC",ifelse(diff_pp<0,"Healthy","Equal"))

b2_df<-data.frame(species=crc_sum$species,
                  difference=diff_pp,
                  Absoulte=diff_abs,
                  Category=diff_category)
print("For Std. Dev Difference")
View(t(t(b2_df)))


#Part C: Probability and independence
#C1: Getting the probablities
#probablity of belonging to CRC group
prob_crc<- sum(metadata$disease=='CRC')/nrow(metadata)
cat("Probability of belonging to CRC group:\n",prob_crc)
#probability of being a female
prob_fem <- sum(metadata$gender=='female')/nrow(metadata)
cat("\nProbability of being a female:\n",prob_fem)
#probablity of being a female and belonging to CRC group
prob_fem_crc<-sum(metadata$gender=='female' & metadata$disease=='CRC')/nrow(metadata)
cat("\nProbability of being a female and belonging to CRC group:\n",prob_fem_crc)

#C2: Calculating the marginal probablities and performing chi2 test
prob_fem_crc_marg <- prob_crc * prob_fem
comp_df<- data.frame(
  values=c("Joint Prob. P(CRC and Female)","Marginal Prob. P(CRC)*P(Female)","Difference"),
  prob= round(c(prob_fem_crc,prob_fem_crc_marg,abs(prob_fem_crc-prob_fem_crc_marg)),4)
)
View(comp_df)

#performing the chi2 test
obs <- table(disease=metadata$disease,gender=metadata$gender)
chi <- chisq.test(obs)
chi

#C3: Comparing the simulation based probability
#Theoretical probabilities for exactly 4 CRCs adn atmot 4CRCs
sample=10
prob_exact_4 <- dbinom(4,size=sample, prob = prob_crc)
prob_atmost_4 <- pbinom(4,size=sample,prob=prob_crc)

#Simulating over 1000 groups of 10
set.seed(123)
prob_sim <- rbinom(1000,size=sample,prob=prob_crc)
sim_exact_4 <- mean(prob_sim==4)
sim_atmost_4 <- mean(prob_sim<=4)

c3_df<- data.frame(
  probability = c("P(Exact 4)","P(atmost 4)"),
  theoretical_prob = round(c(prob_exact_4,prob_atmost_4),4),
  simulated_prob = round(c(sim_exact_4,sim_atmost_4),4)
)
View(c3_df)

#Part D: Normal Distribution and Simulation
#D1:Meand and Std of observed age and simulated age

#observed age mean
mean_age<-mean(metadata$age)
sd_age <- mean(metadata$age)

#simulated age mean
set.seed(123)
sim_age <- rnorm(1000,mean=mean_age,sd=sd_age)

sim_age_mean <- mean(sim_age)
sim_age_sd <- sd(sim_age)

#comparison between the given age and simulated age
d1_df<- data.frame(
  age_source = c("Given","Simulated"),
  mean = round(c(mean_age,sim_age_mean),2),
  sd = round(c(sd_age,sim_age_sd),2)
)
cat("D1:\n\n")
View(d1_df)

#D2: Density and probabilty for age at 50 years

#density
den_50 <- dnorm(50,mean = mean_age, sd = sd_age)

#probability
prob_atmost_50 <- pnorm(50,mean = mean_age, sd = sd_age)

#comparison between the density and probability
d2_df<- data.frame(
  prob_source = c("Density at age 50","Prob. of age <=50"),
  val = round(c(den_50,prob_atmost_50),2)
)
cat("D2:\n\n")
View(d2_df)

#D3: Getting the theoretical prob for ages 40> and <=60.

#probability of observed proportions
prob_40_60 <- sum(metadata$age>=40 & metadata$age<=60)/nrow(metadata$age)

#theoretical probability
theor_prob_40_60 <- pnorm(60,mean_age, sd_age) - pnorm(40,mean_age,sd_age)

#comparison
d3_df<- data.frame(
  prob_source = c("Obs. Prob","Theoretical Prob."),
  val = round(c(prob_40_60,theor_prob_40_60),2)
)
cat("D3:\n\n")
View(d3_df)

#D4: Checking whether the provided data follows normality

s_w <- shapiro.test(metadata$age)

s_w_test <- unname(s_w$statistic)
s_w_pval <- s_w$p.value

cat("D4:\n\n")
cat("Shapiro-Wilk T-stat:",s_w_test)
cat("\nShapiro-Wilk P-val:",s_w_pval)