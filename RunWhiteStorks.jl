# Script to run the model for the white storks

using DataFrames

include("migrationConnectivityModel.jl")

##  Initialisation of the system  ##
##  The distance matrices and distribution of energy supply could be imported from real-world systems  ##

distanceMatrix_west = Array(readtable("distanceMatrix_west.csv", header=false))
distanceMatrix_east = Array(readtable("distanceMatrix_east.csv", header=false))

energySupply_winter = Array(readtable("energySupply_winter.csv", header=false))      # energy supply is derived from NDVI at the corresponding season
energySupply_summer = Array(readtable("energySupply_summer.csv", header=false))
energySupply_summer = reshape(energySupply_summer,1, 1081)

#energySupply_winter = zeros(length(distanceMatrix_west[:,1])) + 1000      # It is assumed here that the energy supply is homogeneously distributed across the species range
#energySupply_summer = zeros(length(distanceMatrix_west[1,:])) + 1000
#energySupply_summer = reshape(energySupply_summer,1, 1081)


##  Model parameters  ##

alph = 0.01    # Cost of movement
bet = 5        # Growth rate
lambd = 1      # Strength of the conformity bias


result = migraConnectivityModel(distanceMatrix_west, distanceMatrix_east, energySupply_winter, energySupply_summer, alph, bet, lambd)

writecsv("modelOutput_whiteStorks_west.csv", result[1])
writecsv("modelOutput_whiteStorks_east.csv", result[2])
