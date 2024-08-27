#LIBRARIES
library(dplyr)
library(ggplot2)

#data paths
data_path = "./data/fish_size/"
save_figure_path = "./data/fish_size/figures/" 
save_table_path = "./data/fish_size/tables/" 

#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#

fish_size = read.csv(paste0(data_path, "biometric_data.csv"))
str(fish_size)

#fish size comparisons

## Muddyfoot ##
muddy.lm.size.roach = lm(Total_length ~ Treatment, data = subset(tagsize, Species == 'Roach' & Lake == 'Muddyfoot'))
summary(muddy.lm.size.roach)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   18.9933     0.4469    42.5   <2e-16 ***
#   TreatmentMix  -0.8362     0.6432    -1.3    0.205   

muddy.lm.size.perch = lm(Total_length ~ Treatment, data = subset(tagsize, Species == 'Perch' & Lake == 'Muddyfoot'))
summary(muddy.lm.size.perch)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   20.4867     0.3716  55.132   <2e-16 ***
#   TreatmentMix  -1.0467     0.5255  -1.992   0.0562 .


## BT ##
BT.lm.size.roach = lm(Total_length ~ Treatment, data = subset(tagsize, Species == 'Roach' & Lake == 'BT'))
summary(BT.lm.size.roach)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   17.2400     0.3563  48.389   <2e-16 ***
#   TreatmentMix  -0.4933     0.5039  -0.979    0.336   

BT.lm.size.perch = lm(Total_length ~ Treatment, data = subset(tagsize, Species == 'Perch' & Lake == 'BT'))
summary(BT.lm.size.perch)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   20.1600     0.8036   25.09   <2e-16 ***
#   TreatmentMix  -0.6133     1.1364   -0.54    0.594  

## Cow Paradise ##
Cow.lm.size.roach = lm(Total_length ~ Treatment, data = subset(tagsize, Species == 'Roach' & Lake == 'Cow Paradise'))
summary(Cow.lm.size.roach)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   15.8600     0.3516  45.107   <2e-16 ***
#   TreatmentMix   0.4300     0.4972   0.865    0.399  

Cow.lm.size.perch = lm(Total_length ~ Treatment, data = subset(tagsize, Species == 'Perch' & Lake == 'Cow Paradise'))
summary(Cow.lm.size.perch)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   19.4150     0.5027   38.62   <2e-16 ***
#   TreatmentMix  -0.3700     0.7109   -0.52    0.606    
