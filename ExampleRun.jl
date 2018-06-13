# Script to run an example simulation

include("migrationConnectivityModel.jl")

##  Initialisation of the system  ##
##  The distance matrices and distribution of energy supply could be imported from real-world systems  ##

Nw = 4          # Number of wintering sites
Ns = 3          # Number of breeding sites
energySupply_winter = zeros(Nw) + 100
energySupply_summer = zeros(Ns) + 100
energySupply_summer = reshape(energySupply_summer,1, Ns)

energySupply_winter[4] = 130
energySupply_summer[2] = 120
energySupply_summer[3] = 200

distanceMatrix_west = zeros(Nw, Ns)
distanceMatrix_west[:,1] = [1, 2, 3, 3.5]
distanceMatrix_west[:,2] = [2, 3, 4, 4.5]
distanceMatrix_west[:,3] = [3, 4, 5, 5.5]

distanceMatrix_east = zeros(Nw, Ns)
distanceMatrix_east[:,1] = [5, 4, 3, 3.5]
distanceMatrix_east[:,2] = [4, 3, 2, 2.5]
distanceMatrix_east[:,3] = [3, 2, 1, 1.5]


##  Model parameters  ##

alph = 0.7    # Cost of movement
bet = 10        # Growth rate
lambd = 1      # Strength of the conformity bias


result = migraConnectivityModel(distanceMatrix_west, distanceMatrix_east, energySupply_winter, energySupply_summer, alph, bet, lambd)

result

# Average distance travelled by individuals of the species
averageDist = (sum(result[1] .* distanceMatrix_west) + sum(result[2] .* distanceMatrix_east)) / (sum(result[1]) + sum(result[2]))

# Average distance travelled by individuals around the migratory divide
averageDist = (sum(result[1] .* distanceMatrix_west, 1)[2] + sum(result[2] .* distanceMatrix_east, 1)[2]) / (sum(result[1],1)[2] + sum(result[2],1)[2])


pathwaysWest = result[1]
pathwaysWest[pathwaysWest.<0.1] = 0
pathwaysEast = result[2]
pathwaysEast[pathwaysEast.<0.1] = 0
pathways = pathwaysWest + pathwaysEast
prevalenceWest = sum(pathwaysWest, 1) ./ sum(pathways, 1)

sum(pathwaysWest, 1)
sum(pathwaysEast, 1)


writecsv("modelOutput_pathway1.csv", result[1])
writecsv("modelOutput_pathway2.csv", result[2])








Nw = 2          # Number of wintering sites
Ns = 1          # Number of breeding sites
energySupply_winter = zeros(Nw) + 100
energySupply_summer = zeros(Ns) + 100
energySupply_summer = reshape(energySupply_summer,1, Ns)

energySupply_winter[1] = 150

distanceMatrix_west = zeros(Nw, Ns)
distanceMatrix_west[:,1] = [1, 2]
distanceMatrix_east = zeros(Nw, Ns)
distanceMatrix_east[:,1] = [2, 1]

alph = 0.1    # Cost of movement
bet = 5        # Growth rate
lambd = 2      # Strength of the conformity bias
result = migraConnectivityModel(distanceMatrix_west, distanceMatrix_east, energySupply_winter, energySupply_summer, alph, bet, lambd)

result

# Average distance travelled by individuals of the species
averageDist = (sum(result[1] .* distanceMatrix_west) + sum(result[2] .* distanceMatrix_east)) / (sum(result[1]) + sum(result[2]))

pathwaysWest = result[1]
pathwaysWest[pathwaysWest.<0.1] = 0
pathwaysEast = result[2]
pathwaysEast[pathwaysEast.<0.1] = 0
pathways = pathwaysWest + pathwaysEast
prevalenceWest = sum(pathwaysWest, 1) ./ sum(pathways, 1)

sum(pathwaysWest, 1)
sum(pathwaysEast, 1)
