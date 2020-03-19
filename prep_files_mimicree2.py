# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 17:51:34 2020

@author: Neda
"""
import numpy as np

'''
#######################################################
#1.make haplotype file 
#this file will NOT be used for simulations directly
#but will be used for estimation of allele frequencies
#and selection of selected loci 
#######################################################

#The phased haplotype data is available from https://academic.oup.com/gbe/article/11/4/1345/5454723 (deposited in Dryad)
#you can convert the above data to a tab-delimited file containing 193 columns: chromosome, position, reference, major/minor alleles and 189 phased haplotype
#then (as done below) convert the the haplotype data to mimicree format (including only chromosomes 2 and 3) and make population size of 450 and 9000 diploid individuals
'''

#A. population of 450 diploids 
InputFile = open('all_biallelicSNPs_q50.hapmat','r')
OutputFile = open('FlLines_FreeBayes_biallelicSNPs_q50_450Ne_No4X.mimhap','w')

base = ['A','T','G','C']
linenum = 0
for line in InputFile:
    linenum += 1
    line = line.rstrip()
    if linenum == 1:
        header = line
    else:
        cols = line.split('\t')
        if cols[0] != 'X':
            biallele = list(set([item for item in cols[3:]]))
            if len(biallele) > 1:
                if biallele[0] in base and biallele[1] in base:
                    count1,count2 = cols[3:].count(biallele[0]),cols[3:].count(biallele[1])
                    if count1 > count2: 
                        ancestral,derived = biallele[0],biallele[1]
                        OutputFile.write('\t'.join(cols[0:3])+'\t'+ ancestral+'/'+derived+'\t'+' '.join([item+item for item in cols[3:]*2]+[item+item for item in cols[3:75]])+'\n') 
                    if count1 < count2: 
                        ancestral,derived = biallele[1],biallele[0] 
                        OutputFile.write('\t'.join(cols[0:3])+'\t'+ ancestral+'/'+derived+'\t'+' '.join([item+item for item in cols[3:]*2]+[item+item for item in cols[3:75]])+'\n')
InputFile.close()
OutputFile.flush()

#B. population of 9000 diploids
InputFile = open('all_biallelicSNPs_q50.hapmat','r')
OutputFile = open('FlLines_FreeBayes_biallelicSNPs_q50_9000Ne_No4X.mimhap','w')

base = ['A','T','G','C']
linenum = 0
for line in InputFile:
    linenum += 1
    line = line.rstrip()
    if linenum == 1:
        header = line
    else:
        cols = line.split('\t')
        if cols[0] != 'X':
            biallele = list(set([item for item in cols[3:]]))
            if len(biallele) > 1:
                if biallele[0] in base and biallele[1] in base:
                    count1,count2 = cols[3:].count(biallele[0]),cols[3:].count(biallele[1])
                    if count1 > count2: 
                        ancestral,derived = biallele[0],biallele[1]
                        OutputFile.write('\t'.join(cols[0:3])+'\t'+ ancestral+'/'+derived+'\t'+' '.join([item+item for item in cols[3:]*47]+[item+item for item in cols[3:120]])+'\n') 
                    if count1 < count2: 
                        ancestral,derived = biallele[1],biallele[0] 
                        OutputFile.write('\t'.join(cols[0:3])+'\t'+ ancestral+'/'+derived+'\t'+' '.join([item+item for item in cols[3:]*47]+[item+item for item in cols[3:120]])+'\n')
InputFile.close()
OutputFile.flush()

'''
########################################################
#2. select beneficial alleles with frequency of 0.5
# and prepare a file containing selected sites for both 
#sweep (sel.txt) and quantitative trait (effect_size.txt) 
#simulations and prepare haplotype file (.mimhap file)
########################################################

The below script will generate required files for population of 450
for a population of 9000 just specify Ne = 9000 and input and output files

'''

#A. population of 450 diploids

#for polymorphic positions (only for chromosome 2 and 3) save the frequency of major and minor alleles from haplotype data (step 1, above)
#the information for each chromosome arm is stored in a separate list. Each SNP is stored as [position,  1st allele (major), 2nd allele (minor), frequency of 1st allele, frequency of 2nd allele]
#first store info about frequency of major and minor alleles for all SNPs
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/FlLines_FreeBayes_biallelicSNPs_q50_450Ne_No4X.mimhap','r')

Ne = 450
haplo_AF_2l,haplo_AF_2r,haplo_AF_3l,haplo_AF_3r = [],[],[],[]
chrom = {'2L':haplo_AF_2l, '2R':haplo_AF_2r, '3L':haplo_AF_3l, '3R':haplo_AF_3r}

for line in InputFile:
    line = line.rstrip()
    cols = line.split('\t')
    Allele1,Allele2 = cols[3].split('/')
    if cols[0] in chrom:
        chrom[cols[0]].append([int(cols[1]),Allele1, Allele2,round(cols[4].count(Allele1)/(Ne*2),4), round(cols[4].count(Allele2)/(Ne*2),4)]) 
InputFile.close()  

########################################################################
#then randomly pick 100 loci, 25 in each chr arm with p0 of 0.05+/-0.005
########################################################################

'''
Here I chose 25 loci in each chromosome arm: a total of 100 loci
for simulations of 10, 20 and 50 loci you just need to change the number of
selected loci in each chromosome arm
'''

pos_2l = np.random.choice([item[0] for item in haplo_AF_2l if 0.045 <= item[4] <= 0.055],25,replace=False)
pos_2r = np.random.choice([item[0] for item in haplo_AF_2r if 0.045 <= item[4] <= 0.055],25,replace=False)
pos_3l = np.random.choice([item[0] for item in haplo_AF_3l if 0.045 <= item[4] <= 0.055],25,replace=False)
pos_3r = np.random.choice([item[0] for item in haplo_AF_3r if 0.045 <= item[4] <= 0.055],25,replace=False)

#to include only selected sites in the haplotypes file
chrom = {'2L':pos_2l, '2R':pos_2r, '3L':pos_3l, '3R':pos_3r}

#select the selected SNPs from the haplotype file (from step 1) and store it as a new file (this file will be used for simulations)
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/mimicree2_sims/sims_proposal/FlLines_FreeBayes_biallelicSNPs_q50_450Ne_No4X.mimhap','r')
OutputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/mimicree2_sims/sims_proposal/Ne450_100loci_0.05p0_0.08s/FlLines_FreeBayes_biallelicSNPs_q50_450Ne_No4X_100loci.mimhap','w')
for line in InputFile:
    line = line.rstrip()
    cols = line.split('\t')
    if cols[0] in chrom:
        if int(cols[1]) in chrom[cols[0]]:
            OutputFile.write(line+'\n')
InputFile.close()
OutputFile.flush()

########################################
#store the info about selected sites 
#######################################

selected_sites_2l = [['2L',item[0],item[1]+'/'+item[2],item[4]] for item in haplo_AF_2l if item[0] in pos_2l]
selected_sites_2r = [['2R',item[0],item[1]+'/'+item[2],item[4]] for item in haplo_AF_2r if item[0] in pos_2r]
selected_sites_3l = [['3L',item[0],item[1]+'/'+item[2],item[4]] for item in haplo_AF_3l if item[0] in pos_3l]
selected_sites_3r = [['3R',item[0],item[1]+'/'+item[2],item[4]] for item in haplo_AF_3r if item[0] in pos_3r]

##########################################################
#write the selected sites in a file for sweep simulations
##########################################################

'''
Here I chose s=0.08 but for different selection coefficients (0.02, 0.05 and 0.1)
you just need to change s
'''
#save the selected sites in sel.txt file for sweep simulations
OutputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/mimicree2_sims/sims_proposal/Ne450_100loci_0.5p0/sel.txt','w')

s = 0.08 #selection coefficient
OutputFile.write('[s]'+'\n')
for item in selected_sites_2l+selected_sites_2r+selected_sites_3l+selected_sites_3r:
    OutputFile.write('\t'.join([str(thing) for thing in item[0:3]])+'\t'+str(s)+'\t'+'0.5'+'\n') #0.5 is for co-dominance
OutputFile.flush()


########################################################################
#write the selected sites in a file for quantitative trait simulations
########################################################################

'''
Here I chose effect size of 0.04 but for different effect size (0.08, 0.2 and 0.4)
you just need to change eff_size
'''
#save the selected sites in effect_size.txt file for quantitative trait simulations
OutputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne200_sel_0.5-4.5--2.5-0.3/effect_size.txt','w')

eff_size = 0.04
for item in selected_sites_2l+selected_sites_2r+selected_sites_3l+selected_sites_3r:
    OutputFile.write('\t'.join([item[0],str(item[1]),item[2].split('/')[1]+'/'+item[2].split('/')[0],str(eff_size),'0'])+'\n') 
OutputFile.flush()

          
