#### Version 2. Chris Gaynor. 19th October 2016
# Script for selecting a diverse set of MAGIC population founders

#### Read marker data ----
#Final geno data is a data frame with markers coded as 0 and 1
geno = read.table("Triticeae_DArT_markers_CerealsDB.csv",sep=",",header=T)
#Remove any blank columns
tmp = apply(geno,2,function(x) !all(is.na(x)))
geno = geno[,tmp]
rm(tmp)
#Set line names as row names
row.names(geno) = geno[,1]
geno = geno[,-1]

#Create marker matrix with uniquely numbered alleles
#Allows for more efficient computation of total alleles
#"by" is the maximum number of alleles per loci
tmp = seq(from=0,by=2,length.out=ncol(geno)) 
uniqueGeno = sweep(geno,2,tmp,"+")
rm(tmp)


#### Calculate best lines ----
library(GA)
nLines = 16 #Change this value for the number of desired individuals

#Functions for the genetic algorithm, used by ga()
population = function(object,...){
  #Sets initial solutions
  #Solutions are random samples of nLines
  population = matrix(0,nrow=object@popSize,ncol=object@nBits)
  for(i in 1:nrow(population)){
    population[i,sample.int(object@nBits,nLines)]=1
  }
  return(population)
}
crossover = function(object, parents, ...){
  #Performs crossovers between solutions
  #Crossovers are constrained to maintain solutions with nLines
  fitness = object@fitness[parents]
  parents = object@population[parents, ,drop=FALSE]
  children = matrix(as.double(NA), nrow=2, ncol=ncol(parents))
  fitnessChildren = rep(NA, 2)
  diffSite = which(abs(parents[1,]-parents[2,])==1)
  n=length(diffSite)/2
  nCrossOver=sample(0:n,size=1)
  if(nCrossOver==0){ #No crossover
    children = parents
    fitnessChildren = fitness
  }else if(nCrossOver==n){ #All crossover
    children[1:2, ] = parents[2:1, ]
    fitnessChildren[1:2] = fitness[2:1]
  }else{
    children[1,] = parents[1,]
    children[2,] = parents[2,]
    site1 = diffSite[which(parents[1,diffSite]==1)]
    site0 = diffSite[which(parents[1,diffSite]==0)]
    children[1,site1[1:n]] = 0
    children[1,site0[1:n]] = 1
    children[2,site1[1:n]] = 1
    children[2,site0[1:n]] = 0
  }
  out = list(children = children, fitness = fitnessChildren)
  return(out)
}
mutation = function(object,parent,...){ 
  #Mutates solutions
  #Mutations are constrained to nLines
  mutate = as.vector(object@population[parent, ])
  i = which(mutate==0)
  j = which(mutate==1)
  n = 1 
  i = sample(i, size=n)
  j = sample(j, size=n)
  mutate[c(i,j)] = abs(mutate[c(i,j)]-1)
  return(mutate)
}
objAllele = function(x,geno){
  #Objective function for number of alleles
  #Requires all alleles uniquely coded
  x = as.logical(x)
  tmp = unlist(geno[x,])
  value = length(na.omit(unique(tmp)))
  return(value)
}
objDiv = function(x,geno){
  #Objective function for Nei's genetic diversity
  #Requires biallelic markers coded as 0 and 1
  x = as.logical(x)
  p = colMeans(geno[x,],na.rm=TRUE)
  value = sum(2*p*(1-p),na.rm=TRUE)/ncol(geno)
  return(value)
}

#Run genetic algorithm for number of alleles
genAlgAllele = ga(type="binary",fitness=objAllele,geno=uniqueGeno,nBits=nrow(geno),
                  popSize=nrow(geno)*3, #Number of possibilities tested each iteration
                  run=20, #Stop when this number of iteration fails to achieve a better solution
                  maxiter=500, #Stop when reaching this number of iterations
                  population=population,
                  crossover=crossover,
                  mutation=mutation,
                  maxFitness=objAllele(rep(1,nrow(uniqueGeno)),uniqueGeno))
#Get names of selected lines and print summary information
tmp = as.logical(genAlgAllele@solution[1,])
takeAllele = row.names(geno)[tmp]
cat("\nSelected lines:",takeAllele,"\n")
cat("Total number of alleles:",objAllele(tmp,geno=uniqueGeno),"\n")
cat("Diversity:",objDiv(tmp,geno=geno),"\n\n")
rm(tmp)

#Run genetic algorithm for diversity
genAlgDiv = ga(type="binary",fitness=objDiv,geno=geno,nBits=nrow(geno),
               popSize=nrow(geno)*3, #Number of possibilities tested each iteration
               run=20, #Stop when this number of iteration fails to achieve a better solution
               maxiter=500, #Stop when reaching this number of iterations
               population=population,
               crossover=crossover,
               mutation=mutation)
#Get names of selected lines and print summary information
tmp = as.logical(genAlgDiv@solution[1,])
takeDiv = row.names(geno)[tmp]
cat("\nSelected lines:",takeDiv,"\n")
cat("Total number of alleles:",objAllele(tmp,geno=uniqueGeno),"\n")
cat("Diversity:",objDiv(tmp,geno=geno),"\n\n")
rm(tmp)
