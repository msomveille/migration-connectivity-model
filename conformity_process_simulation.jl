# Script to simulate conformity

using StatsBase

res = zeros(100)

for k in 1:length(res)

# Local population with individuals' ID and behaviour (0: first behaviour, 1: second behaviour)
popSize = 1000
pop_id = collect(1:1:popSize)
pop_behaviours = zeros(length(pop_id))
pop_behaviours[sample(pop_id, Int(0.9*popSize), replace=false)] = 1

# Sample groups of migrating individuals in the local population
groupSize = 5
pop = collect(1:1:popSize)
groups = zeros(length(pop))
i = 1
while length(pop) > 0
    s = zeros(length(pop))
    s[sample(1:length(s), groupSize, replace=false)] = 1
    pop2 = zeros(length(pop_id))
    pop2[pop[s .== 1]] = 1
    groups = hcat(groups, pop2)
    pop = pop[s .== 0]
    i=i+1
end
groups = groups[:,2:i]

for i in 1:length(groups[1,:])
    if sum(pop_behaviours[groups[:,i] .== 1]) > (groupSize/2)
        pop_behaviours[groups[:,i] .== 1] .= 1
    end
    if sum(pop_behaviours[groups[:,i] .== 1]) < (groupSize/2)
        pop_behaviours[groups[:,i] .== 1] .= 0
    end
    if sum(pop_behaviours[groups[:,i] .== 1]) == (groupSize/2)
        pop_behaviours[groups[:,i] .== 1] .= sample([1,0],1)
    end
end

res[k] = sum(pop_behaviours)

end

mean(res)





res = vector()

pop.behav = c(rep("W", 20), rep("E", 80))

for(i in 1:length(spl)){
	if(length(which(pop.behav[spl[[i]]] == "W")) > 5){
		pop.behav[spl[[i]]] <- rep("W", 10)
	}
	if(length(which(pop.behav[spl[[i]]] == "W")) < 5){
		pop.behav[spl[[i]]] <- rep("E", 10)
	}
	if(length(which(pop.behav[spl[[i]]] == "W")) == 5){
		pop.behav[spl[[i]]] <- rep(sample(c("W","E"),1), 10)
	}
}

res = c(res, length(which(pop.behav == "W")))
