
function migraConnectivityModel(distanceMatrix_1, distanceMatrix_2, supply_winter, supply_summer, alpha, beta, lambda)

    Nw = length(distanceMatrix_1[:,1])  # Number of wintering sites
    Ns = length(distanceMatrix_1[1,:])   # Number of breeding sites

    # Initialise the pathways
    pathways1 = zeros(Nw, Ns) + 0.01
    pathways2 = zeros(Nw, Ns) + 0.01
    pathways = pathways1 + pathways2

    # Conformity function
    function conformity(x, lambda)
	   L = ((x ./ (1 .- x)).^lambda) ./ (1 .+ ((x ./ (1 .- x)).^lambda))
	   L[x .== 1] = 1
	   return(L)
    end

    # Compute the migration cost for each pathway
    migrationCost_1 = exp(-alpha.*distanceMatrix_1)
    migrationCost_2 = exp(-alpha.*distanceMatrix_2)

    # Run the model
    diff = 1000
    t = 1
    while sum(diff) > 1

        # Compute the competition component
        compet = exp(-(sum(pathways, 1) ./ supply_summer)) .* exp(-(sum(pathways, 2) ./ supply_winter))

        # Compute the conformity component for the first pathway
        conform_1 = sum(pathways1,1) ./ sum(pathways, 1)
        conform_1[sum(pathways, 1) .== 0] = 0
        conform_1 = conformity(conform_1, lambda)

        # Compute the conformity component for the second pathway
        conform_2 = 1 .- conform_1

        # Compute the population dynamics for each pathway
        pathways1b = beta .* ((pathways1 ./ sum(pathways1,1)) .* sum(pathways,1) .* conform_1) .* migrationCost_1 .* compet
        pathways1b[isnan(pathways1b)] = 0
        pathways2b = beta .* ((pathways2 ./ sum(pathways2,1)) .* sum(pathways,1) .* conform_2) .* migrationCost_2 .* compet
        pathways2b[isnan(pathways2b)] = 0

        diff = sum((pathways1b .- pathways1).^2 .+ (pathways2b - pathways2).^2)
        # Measure of difference between the state of the system at time t and t+1, which is used to stop the run when this difference is small enough
        #print(diff)
        pathways1 = pathways1b
        pathways2 = pathways2b
        pathways = pathways1 .+ pathways2
        t=t+1
        print(t-1)
    end
    print(t-1)
    return((pathways1, pathways2))

end
