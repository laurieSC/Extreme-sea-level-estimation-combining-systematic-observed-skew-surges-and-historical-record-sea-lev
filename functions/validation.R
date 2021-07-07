# --------------------------------------------------------------------------------------------------------------------------------------

# Author : Laurie Saint Criq (lauriesaintcriq@gmail.com)

# This file compute the Monte Carlo simulations with a selected method for a given site and save the results in the folder "./output/method".

# --------------------------------------------------------------------------------------------------------------------------------------

# Choose the site of observations and the method to implement
site = "" # "Brest", "Dunkerque", "La Rochelle" or "Saint Nazaire"
method = ""   # "1", "2", "3", "4", "5" or "6"

# loading of the function validation
source("./functions/validationFunctions.R")

# data loading
load(paste("./data/",site,"_data.RData",sep=""))

# --------------------------------------------------------------------------------------------------------------------------------------

validation(method, site)
