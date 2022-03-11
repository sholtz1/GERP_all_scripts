R scripts to get format fitness, allele grequency, and GERP data and to combine the three for analysis. 

###Beetle Fitness#####
Starting file "F2_census.csv"

run Rscript Format_beetle_fitness
This script calcultes fitness from groth rate data. 
Final file "Beetle_fitness_filtered.delim" 
Landsape is the number of the replicate from the Experiment
Location is wether the sampled pool of beetle came from the core or the edge of the landscape
Trement S= structured C= controll or shuffeled
Fitness is the average growth rate of beetles for that sampled location and landscape

Landscape	Location	Treatment	Fitness
21		core		S		6.679167
21		edge		S		8.500000
22		core		S		4.933333
22		edge		S		5.016667


### Allele frequencies ###
starting file "Allele_frequencies.delim" 
This comes from the end of the Pooledseq_to_allelefrq pipeline that is run on the cluster
run Rscript Format_allele_frqs