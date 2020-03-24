data_full = read.csv(file = '../data/database.csv')
ind  = which(data_full$OriginalID == 'K_Adamberg_Lactococcus_lactis_growth')
name = data_full$OriginalID[ind]
std_mu_max = data_full$SpecificTraitValue[ind]
temp = data_full$AmbientTemp[ind]
minimal_data_example = data.frame('id' = name, 'temp' = temp, 'std_mu_max' = std_mu_max)
write.csv(minimal_data_example, file = '../data/minimal_example_data.csv')
