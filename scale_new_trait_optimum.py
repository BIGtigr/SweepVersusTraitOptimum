# -*- coding: utf-8 -*-
"""
@author: neda
"""

'''
For the simulations minimum fitness, maximum fitness and standard deviation are 0.5, 4.5, and 0.3.
Using the below script you can set the new trait optimum for populations with different number of loci
and loci with different effect sizes
'''

#given the starting frequency and effect size computes the phenotype of one locus
def compute_trait_phenotype(minor_al_freq,eff_size):
    major_al_freq = 1-minor_al_freq
    d = 0 #d=0 specifies co-dominance
    pheno = eff_size*(minor_al_freq-major_al_freq) + 2*d*minor_al_freq*major_al_freq
    return pheno

##########################################################
#adjust the new triat optimum for different number of loci
##########################################################

#compute the contribution of one locus to the phenotype with allele frequency of 0.05 and effect size of 0.04 
pheno_al = compute_trait_phenotype(0.05, 0.04)

#compute the mean phenotype of populations with different number of loci
loci_num = [100,50,20,10]
initial_pheno = []
for i in loci_num:
    initial_pheno.append(pheno_al*i)

#compute the new optimum given 1.1 shift for all scenarios  
new_pheno = [item+1.1 for item in initial_pheno] 

#################################################################
#adjust the new triat optimum for loci of different effect sizes
#################################################################

#compute the mean phenotype of populations with 100 loci with allele frequency of 0.05 and different effect sizes
eff_size = [0.04,0.08,0.2,0.4]
initial_pheno2 = []
for i in eff_size:
    initial_pheno2.append(compute_trait_phenotype(0.05, i)*100)

#compute the new optimum given 1.1 shift for all scenarios  
new_pheno2 = [item+1.1 for item in initial_pheno2] 



