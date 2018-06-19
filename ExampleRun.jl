# Script to run an example simulation

include("migrationConnectivityModel.jl")

##  Initialisation of the system  ##
##  The distance matrices and distribution of energy supply could be imported from real-world systems  ##

Nw = 2          # Number of wintering sites
Ns = 2          # Number of breeding sites
energySupply_winter = zeros(Nw) + 100
energySupply_summer = zeros(Ns) + 100

distanceMatrix_west = zeros(Nw, Ns)
distanceMatrix_west[:,1] = [1, 2]
distanceMatrix_west[:,2] = [2, 3]

distanceMatrix_east = zeros(Nw, Ns)
distanceMatrix_east[:,1] = [3, 2]
distanceMatrix_east[:,2] = [2, 1]


##  Model parameters  ##

alph = 0.01    # Cost of movement
bet = 5        # Growth rate
lambd = 1      # Strength of the conformity bias


result = migraConnectivityModel(distanceMatrix_west, distanceMatrix_east, energySupply_winter, energySupply_summer, alph, bet, lambd)

writecsv("modelOutput_pathway1.csv", result[1])
writecsv("modelOutput_pathway2.csv", result[2])
