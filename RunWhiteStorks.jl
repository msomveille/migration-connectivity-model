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

alph = 0.01   # Cost of movement
bet = 5       # Growth rate
lambd = 1     # Strength of the conformity bias

result = migraConnectivityModel(distanceMatrix_west, distanceMatrix_east, energySupply_winter, energySupply_summer, alph, bet, lambd)



##  Summary statistics  ##

# Average distance travelled by individuals of the species
averageDist = (sum(result[1] .* distanceMatrix_west) + sum(result[2] .* distanceMatrix_east)) / (sum(result[1]) + sum(result[2]))

# Average distance travelled by individuals around the migratory divide
switchingPops = (longitude_breedingHexagons .< 16) .* (longitude_breedingHexagons .> 12) .* (latitude_breedingHexagons .> 41)
averageDist_switchingPops = (sum(sum(result[1] .* distanceMatrix_west, 1)[switchingPops]) + sum(sum(result[2] .* distanceMatrix_east, 1)[switchingPops])) / (sum(sum(result[1],1)[switchingPops]) + sum(sum(result[2],1)[switchingPops]))
averageDist_switchingPops

# Number of western and eastern migrants around the migratory divide
westernMigrants = sum(sum(result[1],1)[switchingPops])
easternMigrants = sum(sum(result[2],1)[switchingPops])




# Location of the migratory divide (i.e. longitude of the eastern breeding site dominated by western migrants)
pathwaysWest = result[1]
pathwaysWest[pathwaysWest.<0.1] = 0
pathwaysEast = result[2]
pathwaysEast[pathwaysEast.<0.1] = 0
pathways = pathwaysWest + pathwaysEast
prevalenceWest = sum(pathwaysWest, 1) ./ sum(pathways, 1)
migraDivide = maximum(longitude_breedingHexagons[(latitude_breedingHexagons .> 41) .* (prevalenceWest .> 0.66)])

writecsv("modelOutput_whiteStorks_west.csv", result[1])
writecsv("modelOutput_whiteStorks_east.csv", result[2])




##   Run the model for increasing values of lambda  ##

coordinates_breedingHexagons = Array(readtable("coordinates_breedingHexagons.csv", header=false))
longitude_breedingHexagons = coordinates_breedingHexagons[:,1]
latitude_breedingHexagons = coordinates_breedingHexagons[:,2]
longitude_breedingHexagons = reshape(longitude_breedingHexagons, 1, 1081)
latitude_breedingHexagons = reshape(latitude_breedingHexagons, 1, 1081)

migraDivide = zeros(20)
for i in 1:20
    result = migraConnectivityModel(distanceMatrix_west, distanceMatrix_east, energySupply_winter, energySupply_summer, alph, bet, lambd)

    pathwaysWest = result[1]
    pathwaysWest[pathwaysWest.<0.1] = 0
    pathwaysEast = result[2]
    pathwaysEast[pathwaysEast.<0.1] = 0
    pathways = pathwaysWest + pathwaysEast
    prevalenceWest = sum(pathwaysWest, 1) ./ sum(pathways, 1)

    migraDivide[i] = maximum(longitude_breedingHexagons[(latitude_breedingHexagons .> 41) .* (prevalenceWest .> 0.66)])

    lambd = lambd + 0.1
end

writecsv("migraDivide_output.csv", migraDivide)
