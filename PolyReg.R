################################################################################
############### Data Analysis ##################################################
library(moments)
library(BSDA)
library(leaps)
## Data sets
alloy<-read.csv(file="alloy2323093.csv", header=TRUE)
head(alloy)

### Question 1.
### a.
boxplot(alloy, xlab="Alloy", ylab="Megapascals", main="Vickers hardness tests score",
        col="skyblue2")

hist(alloy$V1, main="Vickers hardness tests score", xlab="Alloy", col="skyblue2", 
     breaks=5 )

print(skewness(alloy))
### b.

### c.
mean(alloy$V1)
sd(alloy$V1)
n=length(alloy$V1)
n

## Let X be the distribution of alloy hardness
## Z=(x-mean)/sd
## when x < 850
z850=(850-mean(alloy$V1))/sd(alloy$V1)
z880=(880-mean(alloy$V1))/sd(alloy$V1)
z840=(840-mean(alloy$V1))/sd(alloy$V1)
z860=(860-mean(alloy$V1))/sd(alloy$V1)

# i) P(X<850)
mosaic::xpnorm(z850)
1-0.06876616

# ii) P(X>880)
mosaic::xpnorm(z880)

# iii) P(840<=X<=860)
mosaic::xpnorm(z840)
mosaic::xpnorm(z860)
0.7656787-(1-0.8057)

# d. CI of the mean
ul<-mean(alloy$V1)+1.96*sd(alloy$V1)/sqrt(n)
ll<-mean(alloy$V1)-1.96*sd(alloy$V1)/sqrt(n)
c(ll, ul)

## Hypothesis test group mean
## Two sided z-test
##  H0 mean=860
##  H1 mean≠860
##t.test(alloy$V1, mu = 860, alternative = "two.sided")
z.test(alloy$V1, sigma.x=12.6031, alternative="two.sided", mu=860,
       conf.level=0.95)

### Question 2 ################################################################

cola<-read.csv("cola2323093.csv", header = TRUE)
head(cola)

summary(cola)
c(sd(cola$temperature), sd(cola$volume))
c(IQR(cola$temperature), IQR(cola$volume))

# a.
cola.lm<-lm(volume~temperature, data=cola)
par(mfrow=c(1,1))
plot(cola$temperature, cola$volume, xlab="Temperature (in degrees celcius)", 
     col="blue2",
     ylab="Volume (in 1000 litres)", main="Temperature vs cola volume produced")
abline(cola.lm$coef)
summary(cola.lm)

# b.
#### Assumptions on ϵi mean=0 and var=sigma^2 fir i=1,2,...n 
#### Errors are independent of each other.

par(mfrow=c(2,2))
plot(cola.lm, col="skyblue3")
cor(cola$temperature, cola$volume)

# c. Suggest a quadratic model curved residuals.

# d. Quadratic model

par(mfrow=c(1,1))
plot(cola$temperature, cola$volume, xlab="Temperature (in degrees celcius)", 
     col="blue2",
     ylab="Volume (in 1000 litres)", main="Temperature vs cola volume produced") 
cola2.lm <- lm(volume ~ temperature + I(temperature^2), data=cola)

#use model to get predicted values
pred <- predict(cola2.lm)
ix <- sort(cola$temperature, index.return=T)$ix
#add polynomial curve to plot
lines(cola$temperature[ix], pred[ix], lwd=2)

# checking model assumptions
par(mfrow=c(2,2))
plot(cola2.lm, col="skyblue3")
summary(cola2.lm)
c(AIC(cola.lm), AIC(cola2.lm))
c(BIC(cola.lm), BIC(cola2.lm))

## Model Y=171.16206 + 17.56727x - 0.34179x^2
y<-171.16206 + 17.56727*20 - 0.34179*(20^2)
y
summary(cola2.lm)

#### Question 3. ###############################################################

magnia<-read.csv("magnia2323093.csv", header = TRUE)
colnames(magnia)<-c("TPSA","SAacc","H_050","MLOGP","RDCHI","GATS1p","nN",
                    "C_040","LC50")
head(magnia)
str(magnia)
### Data summary
summary(magnia)
par(mfrow=c(3,3))
### Histogram of the distribution of variables
hist(magnia$LC50, main="Distribution of LC50", xlab="LC50", col="skyblue") ######
hist(magnia$TPSA, main="Distribution of TPSA", xlab="TPSA", col="skyblue")
hist(magnia$SAacc, main="Distribution of SAacc", xlab="SAacc", col="skyblue")
hist(magnia$H_050, main="Distribution of H_050", xlab="H_050", col="skyblue")
hist(magnia$MLOGP, main="Distribution of MLOGP", xlab="MLOGP", col="skyblue")
hist(magnia$RDCHI, main="Distribution of RDCHI", xlab="RDCHI", col="skyblue")
hist(magnia$GATS1p, main="Distribution of GATS1p", xlab="GATS1p", col="skyblue")
hist(magnia$nN, main="Distribution of nN", xlab="nN", col="skyblue")
hist(magnia$C_040, main="Distribution of C_040", xlab="C_040", col="skyblue")

##### Boxplot
boxplot(magnia$LC50, main="Distribution of LC50", xlab="LC50", col="skyblue") ###
boxplot(magnia$TPSA, main="Distribution of TPSA", xlab="TPSA", col="skyblue")
boxplot(magnia$SAacc, main="Distribution of SAacc", xlab="SAacc", col="skyblue")
boxplot(magnia$H_050, main="Distribution of H_050", xlab="H_050", col="skyblue")
boxplot(magnia$MLOGP, main="Distribution of MLOGP", xlab="MLOGP", col="skyblue")
boxplot(magnia$RDCHI, main="Distribution of RDCHI", xlab="RDCHI", col="skyblue")
boxplot(magnia$GATS1p, main="Distribution of GATS1p", xlab="GATS1p", col="skyblue")
boxplot(magnia$nN, main="Distribution of nN", xlab="nN", col="skyblue")
boxplot(magnia$C_040, main="Distribution of C_040", xlab="C_040", col="skyblue")

par(mfrow=c(1,1))
####### Plot matrix to show multiple relationships
pairs(~ TPSA + SAacc + H_050 + MLOGP + RDCHI + GATS1p + nN + C_040 + LC50,
      data=magnia,
      main="Response variable (LC50)  versus explanatory variables", col="skyblue3")

magnia.lm<-lm(LC50~ .,data=magnia)
summary(magnia.lm)

### Perform step-wise backwards regression######################################
slm1<-step(magnia.lm)

## final model backward selection AIC
##                    LC50 = TPSA + SAacc + MLOGP + RDCHI + GATS1p + nN
summary(slm1)

#####################backward Model selection with BIC ############################
n=length(magnia$TPSA)
slm1.1<-step(magnia.lm, k=log(n))
summary(slm1.1)
#########             LC50 = TPSA + SAacc + MLOGP + RDCHI + GATS1p + nN

### Comparing AIC and BIC for both full model and final model 
AIC(magnia.lm)
BIC(magnia.lm)
AIC(slm1.1)
BIC(slm1.1)

########### Forward selection #################################################
# This is the intercept only model
magnia.lm2<-lm(LC50~1,data=magnia)
summary(magnia.lm2)
slm2<-step(magnia.lm2, scope=formula(magnia.lm), direction='forward')
summary(slm2)
AIC(slm2)
BIC(slm2)
## Final model forward selection
## LC50 = MLOGP + TPSA + SAacc + nN + GATS1p + RDCHI


### Model check
par(mfrow=c(2,2))
plot(slm2, col="skyblue3")

#### Dictation and removal of outliers/influential obs
cooksD <- cooks.distance(slm2)
influential <- cooksD[(cooksD > (3 * mean(cooksD, na.rm = TRUE)))]
influential

names_of_influential <- names(influential)
outliers <- magnia[names_of_influential,]
magnia_without_outliers <- magnia %>% anti_join(outliers)
magnia.reduced <- lm(LC50 ~ MLOGP + TPSA + SAacc + nN + GATS1p + RDCHI, 
             data = magnia_without_outliers)
summary(magnia.reduced)
plot(magnia.reduced, col="skyblue2")
AIC(magnia.reduced)
BIC(magnia.reduced)






