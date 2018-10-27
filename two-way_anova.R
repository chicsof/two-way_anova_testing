##########################################################
########### TWO WAY ANOVA IN R WITH INTERACTION #########
#########################################################

#built in data set
head(warpbreaks)
summary(warpbreaks)

#standart model
model1 <- aov(breaks ~ wool + tension, data = warpbreaks)
#we can see that tension is significant to the breaks
summary(model1)

#add an interaction manually, product of wool and tension
model2 <- aov(breaks ~ wool + tension + wool:tension, data = warpbreaks)
#we see that the interaction of wool with tension (combination) is fairly significant 
summary(model2)

#check for all interactions, should return the same model
model3 <- aov(breaks ~ wool * tension, data = warpbreaks)
summary(model3)

############################## POISON EXAMPLE MANUALLY, from scrach :'] #######################################################

###fix the path for you###
#this is a balanced test since for every poisson the same amount of tests/treatments were used
servivolTimeDS <- read.csv("~/Projects/statistical_analysis/two-way_anova_testing/poison_balanced.tsv", sep="")

#we can see the boxplot to get an idea of the varience
boxplot(survivolTime~ treatment* poison ,data=servivolTimeDS)

summary(servivolTimeDS)

#H0: type of poison has no affect
#H0: type o treatment has no affect
#H0: combination of poison and type has no affect

############################ inspect the mean servivol time for each combination #######################################
library(plyr)
#per poison
meanPerPoison <-ddply(servivolTimeDS, .(poison), summarize, mean=mean(survivolTime))
meanPerPoison
#per treatment
meanPerTreatment <-ddply(servivolTimeDS, .(treatment), summarize, mean=mean(survivolTime))
meanPerTreatment
#for each possion each treatment
meanPoisonTreat <-ddply(servivolTimeDS, .(treatment, poison), summarize, mean=mean(survivolTime))
meanPoisonTreat
#total
meanTimeForAll <- mean(servivolTimeDS$survivolTime)
meanTimeForAll

##########################Sum square of first factor(poison)#############################################################

#this is given by calculating the squared difference of the grand mean to the mean for each 
#poison and then suming the result

#for robustness we make a function that selects the mean for a given poison and does the calculations
sumOfSquaresForPoisonF <- function(poisonGiven){
  #we multiply by 4 because we have 4 treatments per poison
  4*(((subset(meanPerPoison,poison==poisonGiven,select="mean"))[1,]-meanTimeForAll)^2)
}
#then we apply this function for poion 1,2,3 and sum the result
sumOfSquaresForPoison <- sum(mapply(sumOfSquaresForPoisonF, c(1,2,3)))


##########################Sum square of second factor(treatment)######################################################
#as above but for treatments this time

sumOfSquaresForTreatmentsF <- function(treatmentGiven){
  #we multiply by 3 because we have 3 poisons per treatment
  3*(((subset(meanPerTreatment,treatment==treatmentGiven,select="mean"))[1,]-meanTimeForAll)^2)
}


sumOfSquaresForTreatments <- sum(mapply(sumOfSquaresForTreatmentsF, c("A", "B", "C", "D")))
sumOfSquaresForTreatments 
########################## Sum square within error ######################################################################

#This is the square sum for each servivaltime in our dataset minus the avarage for that poison and treatment

sumOfSquaresWithErrorF <- function(treatmentGiven, poisonGiven){
  subTrPoi <- subset(servivolTimeDS, treatment==treatmentGiven & poison== poisonGiven)
  meanOfTrPoi <- subset(meanPoisonTreat, treatment==treatmentGiven& poison== poisonGiven,select="mean")
  sumTrPoi <- sum((subTrPoi$survivolTime - meanOfTrPoi)^2)
  return(sumTrPoi)
}
#to get all combinations A1, A2,A3,B1....we need to use outer product
x <- factor(c("A", "B", "c", "D"))
y <- c(1,2,3)
product <- expand.grid(x, y)

#apply the function to all possible combinations and sum them up
sumOfSquaresWithError <- sum(mapply( sumOfSquaresWithErrorF,treatmentGiven= as.character(product[1,"Var1"]),poisonGiven=product[1,"Var2"]))
sumOfSquaresWithError
########################## Sum of Square Total ################################
sumSquareTotal <- sum((servivolTimeDS$survivolTime - meanTimeForAll)^2) 
sumSquareTotal

######################### sum square of both factors ############################################
#sum of both factors is given by  
sumOfSquareBothFactors <- sumSquareTotal - sumOfSquaresWithError -sumOfSquaresForTreatments - sumOfSquaresForPoison

######################### calculating the degrees of freedom for each sum of squares ############

#for first factor (poison)
dfFirstFactor <- 3-1
#for second factor (treatmeant)
dfSecondFactor <- 4-1
#for within error we add up n-1 of each treatment for each poison so:
#(4-1) a treat ment for a poison *3 one treatment for each poison *4 each treatment for each poison
dfWithinError <- (4-1)*3*4
#sum of both squares, we multiply df of first and second
dfSumOfBoth <- dfFirstFactor * dfSecondFactor
#total degree of freedoms, this is the sum of all of them
dfTotal <- dfFirstFactor + dfSecondFactor + dfWithinError + dfSumOfBoth
dfTotal

######################### calculating the mean square of sum of square within error ############

#we will need this to calculate the f-scores which will allow us to draw our conclusions for each H0
#this is the sum of squares within error devided by its degrees of freedom so:
meanSquareOfSumWithinError <- sumOfSquaresWithError/dfWithinError
meanSquareOfSumWithinError

######################### H0: poison does not affect the survivol time ############################

#we need to calculate the F-score for this whcich is meansquareof1rstfactor/meansquareWithinError
#so we need the mean square of 1rst factor:
meanSquareOfFirstFactor <- sumOfSquaresForPoison/dfFirstFactor 
FscoreForPoison <- meanSquareOfFirstFactor/meanSquareOfSumWithinError
# F(dfFirstFactor ,dfWithinError) = FscoreForPoison p<0.5 for a 95%confidence interval,
#we can find the critical value for df of numerator dfFirstFactor and dfWithinError demonimator from the F disctribution
cvForPoison <- qf(.95, df1=dfFirstFactor, df2=dfWithinError) 
cvForPoison
#FscoreForPoison falls in the rejection aeria and so we can reject he H0 and accept that the poison does affect the
#servival time

####### you can repeat this for the rest of the H0's ##############################################
### or we can be sensible and use..... r #####

#however you need to ensure your factors are actually of type factor!!!
servivolTimeDS$poison <- as.factor(servivolTimeDS$poison)
modelForServivolTime <- lm(survivolTime ~ treatment * poison, data = servivolTimeDS)
#type does not matter since our test is balanced
install.packages("car")
library(car)
Anova(modelForServivolTime, type="III")

#from the above we can reject that poison and treatment has no effect and accept that the interaction has no effect
