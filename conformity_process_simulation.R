
prevals <- seq(0,1,0.1)
y <- vector()
for(k in 1:length(prevals)){

popSize = 1000
groupSize = 10
pop = 1:popSize
spl <- list()
i=1
while(length(pop) > 0){
	s = sample(1:length(pop), groupSize)
	spl[[i]] <- pop[s]
	pop = pop[-s]
	i=i+1
}


prevalW = prevals[k]
pop.behav = c(rep("W", popSize*prevalW), rep("E", (popSize-(popSize*prevalW))))

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


y[k] = length(which(pop.behav == "W")) / popSize

}



lambda=3
yy = ((prevals/(1-prevals))^lambda) / (1+((prevals/(1-prevals))^lambda))
yy[length(yy)] <- 1

plot(prevals, prevals, type='l')
points(prevals, y, type='l', col="red")
points(prevals, yy, type='l', col="blue")

