


pop = 1:100
spl <- list()
i=1
while(length(pop) > 0){
	s = sample(1:length(pop), 10)
	spl[[i]] <- pop[s]
	pop = pop[-s]
	i=i+1
}

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


x = seq(0,1,0.1)
y = c(0, 0.004, 0.018, 0.096, 0.268, 0.499, 0.738, 0.905, 0.982, 0.999, 1)
y2 = c(0, 0.009, 0.058, 0.164, 0.319, 0.499, 0.682, 0.837, 0.944, 0.992, 1)



lambda=2
yy = ((x/(1-x))^lambda) / (1+((x/(1-x))^lambda))

plot(x, x, type='l')
points(x, y, type='l', col="red")
points(x, y2, type='l', col="orange")

points(x, yy, type='l', col="blue")


