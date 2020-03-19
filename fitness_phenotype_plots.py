# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 11:08:18 2020

@author: Neda
"""

import numpy as np
import matplotlib.pyplot as plt

#get the fitness of sweep simulations and categorize them based on timepoint and replicate
#then compute mean and standard deviation from the median fitness of all replicates         
def get_fitness_SimBatch_w(InputFile_f,sim_nums):
    gen_index = {0:0,10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
    fitness = [[[] for _ in xrange(sim_nums)] for _ in xrange(len(gen_index))]
    for line in InputFile_f:
        line = line.rstrip()
        cols = line.split('\t')
        if int(cols[0]) <= sim_nums:
            gen,rep=gen_index[int(cols[1])],int(cols[0])
            fitness[gen][rep-1].append(np.log10(float(cols[5])))
    #compute mean and standard deviation of fitness from median of replicates for each timepoint
    mean_of_median_fit, std_of_mean_fit = [], []
    for tmp in fitness:
        mean_of_median_fit.append(np.mean([np.median(r) for r in tmp]))
        std_of_mean_fit.append(np.std([np.median(r) for r in tmp]))
    #compute the area around the mean for plotting std    
    std_fit_down = [tmp-(std_of_mean_fit[ind]/2) for ind,tmp in enumerate(mean_of_median_fit)]
    std_fit_up = [tmp+(std_of_mean_fit[ind]/2) for ind,tmp in enumerate(mean_of_median_fit)]
    return (mean_of_median_fit,std_fit_up,std_fit_down)
    InputFile_f.close()
    
#get the phenotype of quantitaive trait simulations and categorize them based on timepoint and replicate
#then compute mean and standard deviation from the median phenotype of all replicates 
#note that for the sake of plotting the phenotypes are scaled so they are not negative         
def get_fitness_SimBatch_qff(InputFile_f,sim_nums,scale_fac):
    gen_index = {0:0,10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
    phenotype = [[[] for _ in xrange(sim_nums)] for _ in xrange(len(gen_index))]
    for line in InputFile_f:
        line = line.rstrip()
        cols = line.split('\t')
        if int(cols[0]) <= sim_nums:
            gen,rep=gen_index[int(cols[1])],int(cols[0])
            phenotype[gen][rep-1].append(float(cols[4])-scale_fac)
    #compute mean and standard deviation of phenotype from median of replicates for each timepoint
    mean_of_median_fit, std_of_mean_fit = [], []
    for tmp in phenotype:
        mean_of_median_fit.append(np.mean([np.median(r) for r in tmp]))
        std_of_mean_fit.append(np.std([np.median(r) for r in tmp]))
    #compute the area around the mean for plotting std    
    std_fit_down = [tmp-(std_of_mean_fit[ind]/2) for ind,tmp in enumerate(mean_of_median_fit)]
    std_fit_up = [tmp+(std_of_mean_fit[ind]/2) for ind,tmp in enumerate(mean_of_median_fit)]
    return (mean_of_median_fit,std_fit_up,std_fit_down)
    InputFile_f.close()

'''
###########
#Figure2
###########
'''

###############################
#get fitness data for sweep
#fitness is log 10 transformed
##############################

#w 450
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/recom_rate_noX/1000sims/fitness_result.txt','r')  
mean_fit_w_450_100, std_fit_w_450_100_up, std_fit_w_450_100_down = get_fitness_SimBatch_w(InputFile_f,500)

#w 9k
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.08s/recom_rate_noX/500sims/fitness_result.txt','r')  
mean_fit_w_9k_100, std_fit_w_9k_100_up, std_fit_w_9k_100_down = get_fitness_SimBatch_w(InputFile_f,500)

##############################################################
#get phenotype for QTL
#phenotype is rescaled by the initial phenotype of population
##############################################################

#QTL 450
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/recom_rate_noX/1000sims/GenoPheno.gpf','r')  
loci_num, eff_size, p0 = 100, 0.04, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_450_100, std_fit_qtl_450_100_up, std_fit_qtl_450_100_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

#QTL 9k
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--2.5-0.3/recom_rate_noX/500sims/GenoPheno.gpf','r')  
loci_num, eff_size, p0 = 100, 0.04, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_9k_100, std_fit_qtl_9k_100_up, std_fit_qtl_9k_100_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

#########
#plot
#########

mean_450 = [mean_fit_w_450_100,mean_fit_qtl_450_100]
mean_9k = [mean_fit_w_9k_100,mean_fit_qtl_9k_100]

std_450_up = [std_fit_w_450_100_up,std_fit_qtl_450_100_up]
std_450_down = [std_fit_w_450_100_down,std_fit_qtl_450_100_down]

std_9k_up = [std_fit_w_9k_100_up,std_fit_qtl_9k_100_up]
std_9k_down = [std_fit_w_9k_100_down,std_fit_qtl_9k_100_down]
           
gen_index = {0:0,10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
color_450 = ['turquoise','#e699b2']
color_9k = ['teal','#cc3366']
color_label = [['turquoise','teal'],['#e699b2','#cc3366']]
label_450 = ['450 sweep', '450 trait optimum']
label_9k = ['9000 sweep', '9000 trait optimum']

fig , ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True,figsize=(14,5),dpi=100, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.3, wspace=0.1)
for rep in range(2):
    ax=plt.subplot(1,2,rep+1)
    plt.plot(np.arange(0,(len(gen_index))), mean_450[rep], color = color_450[rep],linewidth = 3, label = label_450[rep])
    plt.fill_between(np.arange(0,(len(gen_index))), std_450_up[rep], std_450_down[rep], color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(np.arange(0,(len(gen_index))), mean_9k[rep], color = color_9k[rep],linewidth = 3, linestyle = '--', label = label_9k[rep])
    plt.fill_between(np.arange(0,(len(gen_index))), std_9k_up[rep], std_9k_down[rep], color = color_9k[rep], alpha = 0.5, edgecolor = None, linewidth=0.0)
    ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
    ax.tick_params(axis='y', which='both' ,right='off', left='on')
    plt.xlim(-0.2,(len(gen_index))-0.8) 
    plt.ylim(0,)
    plt.xticks(np.arange(0,len(gen_index),1), [i for i in sorted(gen_index.keys())],fontsize=12, va='top')
    plt.yticks(fontsize=12, va='top')
    plt.xlabel('Generation', size = 14,labelpad=15)
    if rep == 0:
        ax.set_ylabel('Fitness', size = 14,labelpad=15)
        ax.legend(loc = 3,ncol=1, borderaxespad=0., bbox_to_anchor=(0.733, 1))
    if rep == 1:
        ax.set_ylabel('Phenotype', size = 14)
        ax.legend(loc=3,ncol=1, borderaxespad=0.,bbox_to_anchor=(0.632, 1))
plt.tight_layout()
plt.savefig('sweep_qtl_meanPhenoFitness_sweepLog10_qtlPhenoRescaled_defaultParameters.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_qtl_meanPhenoFitness_sweepLog10_qtlPhenoRescaled_defaultParameters.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''
##########
#Figure S5
##########
 '''
 
##########################
#get log 10 fitness for w 
##########################

#w 450, loci 10, 20, 50, 100
sim = 'Ne450_10loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_450_10, std_fit_w_450_10_up, std_fit_w_450_10_down = get_fitness_SimBatch_w(InputFile_f,500)

sim =  'Ne450_20loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_450_20, std_fit_w_450_20_up, std_fit_w_450_20_down = get_fitness_SimBatch_w(InputFile_f,500)

sim = 'Ne450_50loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_450_50, std_fit_w_450_50_up, std_fit_w_450_50_down = get_fitness_SimBatch_w(InputFile_f,500)

sim =  'Ne450_100loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/1000sims/fitness_result.txt' %sim,'r')  
mean_fit_w_450_100, std_fit_w_450_100_up, std_fit_w_450_100_down = get_fitness_SimBatch_w(InputFile_f,500)

#w 9000, loci 10, 20, 50, 100
sim = 'Ne9000_10loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_9k_10, std_fit_w_9k_10_up, std_fit_w_9k_10_down = get_fitness_SimBatch_w(InputFile_f,500)

sim =  'Ne9000_20loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_9k_20, std_fit_w_9k_20_up, std_fit_w_9k_20_down = get_fitness_SimBatch_w(InputFile_f,500)

sim = 'Ne9000_50loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_9k_50, std_fit_w_9k_50_up, std_fit_w_9k_50_down = get_fitness_SimBatch_w(InputFile_f,500)

sim =  'Ne9000_100loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_9k_100, std_fit_w_9k_100_up, std_fit_w_9k_100_down = get_fitness_SimBatch_w(InputFile_f,500)

#########
#plot
#########
    
mean_450 = [mean_fit_w_450_10,mean_fit_w_450_20,mean_fit_w_450_50, mean_fit_w_450_100]
mean_9k = [mean_fit_w_9k_10,mean_fit_w_9k_20,mean_fit_w_9k_50, mean_fit_w_9k_100]

std_450_up = [std_fit_w_450_10_up,std_fit_w_450_20_up,std_fit_w_450_50_up, std_fit_w_450_100_up]
std_450_down = [std_fit_w_450_10_down,std_fit_w_450_20_down,std_fit_w_450_50_down, std_fit_w_450_100_down]

std_9k_up = [std_fit_w_9k_10_up,std_fit_w_9k_20_up,std_fit_w_9k_50_up, std_fit_w_9k_100_up]
std_9k_down = [std_fit_w_9k_10_down,std_fit_w_9k_20_down,std_fit_w_9k_50_down, std_fit_w_9k_100_down]

gen_index = {0:0,10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
color_label = ['turquoise','teal']
loci_num = [10,20,50,100]

fig , ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True,figsize=(14,10),dpi=100, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.2, wspace=0.1)
for rep in range(4):
    ax=plt.subplot(2,2,rep+1)
    plt.plot(np.arange(0,(len(gen_index))), mean_450[rep], color = 'turquoise',linewidth = 3, label = '450')
    plt.fill_between(np.arange(0,(len(gen_index))), std_450_up[rep], std_450_down[rep], color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(np.arange(0,(len(gen_index))), mean_9k[rep], color = 'teal',linewidth = 3, label = '9000', linestyle = '--')
    plt.fill_between(np.arange(0,(len(gen_index))), std_9k_up[rep], std_9k_down[rep], color = 'teal', alpha = 0.3, edgecolor = None, linewidth=0.0)
    ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
    ax.tick_params(axis='y', which='both' ,right='off', left='on')
    plt.xlim(-0.2,(len(gen_index))-0.8) 
    plt.ylim(0,)
    plt.yticks(fontsize=12, va='top')    
    plt.xticks(np.arange(0,len(gen_index),1), [i for i in sorted(gen_index.keys())],fontsize=12, va='top')
    ax.set_title("%i loci" %(loci_num[rep]),fontsize=12)
    if rep == 1:
        ax.legend(loc = 3,ncol=2, borderaxespad=0., bbox_to_anchor=(0.65, 1))
    if rep == 0 or rep == 2:
        ax.set_ylabel('Fitness', size = 14,labelpad=15)
    if rep == 2 or rep == 3:
        plt.xlabel('Generation', size = 14,labelpad=15)
    if rep == 0 or rep == 1:plt.ylim(-0.02, 0.7)
    if rep == 2 or rep == 3:plt.ylim(-0.02, 3)
plt.tight_layout()
plt.savefig('sweep_10to100loci_0.08s_fitnesslog10.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_10to100loci_0.08s_fitnesslog10s.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''
#########
#Figure 5
#########
'''

################################
#get rescaled phenotype for QTL 
################################

#qff 450, 10, 20, 50, 100 loci
sim = 'sel_0.5-4.5-0.74-0.3_eff0.04_10loci_movingTO'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 10, 0.04, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_450_10, std_fit_qtl_450_10_up, std_fit_qtl_450_10_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

sim = 'sel_0.5-4.5-0.38-0.3_eff0.04_20loci_movingTO'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 20, 0.04, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_450_20, std_fit_qtl_450_20_up, std_fit_qtl_450_20_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

sim = 'sel_0.5-4.5--0.7-0.3_eff0.04_50loci_movingTO'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 50, 0.04, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_450_50, std_fit_qtl_450_50_up, std_fit_qtl_450_50_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

sim = 'sel_0.5-4.5--2.5-0.3'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/recom_rate_noX/1000sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 100, 0.04, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_450_100, std_fit_qtl_450_100_up, std_fit_qtl_450_100_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

#qff 9000, 10, 20, 50, 100 loci
sim = 'Ne9000_sel_0.5-4.5-0.74-0.3_eff0.04_10loci_movingTO'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 10, 0.04, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_9k_10, std_fit_qtl_9k_10_up, std_fit_qtl_9k_10_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

sim = 'Ne9000_sel_0.5-4.5-0.38-0.3_eff0.04_20loci_movingTO'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 20, 0.04, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_9k_20, std_fit_qtl_9k_20_up, std_fit_qtl_9k_20_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

sim = 'Ne9000_sel_0.5-4.5--0.7-0.3_eff0.04_50loci_movingTO'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 50, 0.04, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_9k_50, std_fit_qtl_9k_50_up, std_fit_qtl_9k_50_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

sim = 'Ne9000_sel_0.5-4.5--2.5-0.3'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/recom_rate_noX/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 100, 0.04, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_9k_100, std_fit_qtl_9k_100_up, std_fit_qtl_9k_100_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

#########
#plot
#########
    
mean_450 = [mean_fit_qtl_450_10,mean_fit_qtl_450_20,mean_fit_qtl_450_50, mean_fit_qtl_450_100]
mean_9k = [mean_fit_qtl_9k_10,mean_fit_qtl_9k_20,mean_fit_qtl_9k_50, mean_fit_qtl_9k_100]

std_450_up = [std_fit_qtl_450_10_up,std_fit_qtl_450_20_up,std_fit_qtl_450_50_up, std_fit_qtl_450_100_up]
std_450_down = [std_fit_qtl_450_10_down,std_fit_qtl_450_20_down,std_fit_qtl_450_50_down, std_fit_qtl_450_100_down]

std_9k_up = [std_fit_qtl_9k_10_up,std_fit_qtl_9k_20_up,std_fit_qtl_9k_50_up, std_fit_qtl_9k_100_up]
std_9k_down = [std_fit_qtl_9k_10_down,std_fit_qtl_9k_20_down,std_fit_qtl_9k_50_down, std_fit_qtl_9k_100_down]

gen_index = {0:0,10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
color_label = ['#e699b2','#cc3366']
loci_num = [10,20,50,100]

TO = [item-item+1.1 for item in [0.74,0.38,-0.7,-2.5]]

fig , ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True,figsize=(14,10),dpi=100, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.2, wspace=0.1)
for rep in range(4):
    ax=plt.subplot(2,2,rep+1)
    plt.hlines(y = TO[rep],xmin = -0.2,xmax=len(gen_index), colors = 'gold', linestyles = 'solid',linewidth = 1)
    plt.plot(np.arange(0,(len(gen_index))), mean_450[rep], color = '#e699b2',linewidth = 3, label = '450')
    plt.fill_between(np.arange(0,(len(gen_index))), std_450_up[rep], std_450_down[rep], color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(np.arange(0,(len(gen_index))), mean_9k[rep], color = '#cc3366',linewidth = 3, linestyle = '--', label = '9000')
    plt.fill_between(np.arange(0,(len(gen_index))), std_9k_up[rep], std_9k_down[rep], color = '#cc3366', alpha = 0.3, edgecolor = None, linewidth=0.0)
    ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
    ax.tick_params(axis='y', which='both' ,right='off', left='on')
    plt.xlim(-0.2,(len(gen_index))-0.8) 
    plt.ylim(-0.05)
    plt.yticks(fontsize=12, va='top')    
    plt.xticks(np.arange(0,len(gen_index),1), [i for i in sorted(gen_index.keys())],fontsize=12, va='top')
    ax.set_title("%i loci" %(loci_num[rep]),fontsize=12)
    if rep == 1:
        ax.legend(bbox_to_anchor=(0.65, 1), loc=3,ncol=2, borderaxespad=0.)    
    if rep == 0 or rep == 2:
        ax.set_ylabel('Phenotype', size = 14,labelpad=15)
    if rep == 2 or rep == 3:
        plt.xlabel('Generation', size = 14,labelpad=15)
plt.tight_layout()
plt.savefig('QTL_10to100loci_0.04EffSize_rescaledPheno.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('QTL_10to100loci_0.04EffSize_rescaledPheno.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''
#########
#Figure 7
#########
 '''

################################
#get log10 fitness for w 
################################

#w 450, 0.02,0.05,0.08,0.1
sim = 'Ne450_100loci_0.05p0_0.02s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_450_02, std_fit_w_450_02_up, std_fit_w_450_02_down = get_fitness_SimBatch_w(InputFile_f,500)

sim =  'Ne450_100loci_0.05p0_0.05s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_450_05, std_fit_w_450_05_up, std_fit_w_450_05_down = get_fitness_SimBatch_w(InputFile_f,500)

sim = 'Ne450_100loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/1000sims/fitness_result.txt' %sim,'r')  
mean_fit_w_450_08, std_fit_w_450_08_up, std_fit_w_450_08_down = get_fitness_SimBatch_w(InputFile_f,500)

sim =  'Ne450_100loci_0.05p0_0.1s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_450_1, std_fit_w_450_1_up, std_fit_w_450_1_down = get_fitness_SimBatch_w(InputFile_f,500)

#w 9000, 0.02,0.05,0.08,0.1
sim = 'Ne9000_100loci_0.05p0_0.02s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_9k_02, std_fit_w_9k_02_up, std_fit_w_9k_02_down = get_fitness_SimBatch_w(InputFile_f,500)

sim =  'Ne9000_100loci_0.05p0_0.05s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_9k_05, std_fit_w_9k_05_up, std_fit_w_9k_05_down = get_fitness_SimBatch_w(InputFile_f,500)

sim = 'Ne9000_100loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_9k_08, std_fit_w_9k_08_up, std_fit_w_9k_08_down = get_fitness_SimBatch_w(InputFile_f,500)

sim =  'Ne9000_100loci_0.05p0_0.1s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
mean_fit_w_9k_1, std_fit_w_9k_1_up, std_fit_w_9k_1_down = get_fitness_SimBatch_w(InputFile_f,500)

##########
#plot
#########
    
mean_450 = [mean_fit_w_450_02,mean_fit_w_450_05,mean_fit_w_450_08, mean_fit_w_450_1]
mean_9k = [mean_fit_w_9k_02,mean_fit_w_9k_05,mean_fit_w_9k_08, mean_fit_w_9k_1]

std_450_up = [std_fit_w_450_02_up,std_fit_w_450_05_up,std_fit_w_450_08_up, std_fit_w_450_1_up]
std_450_down = [std_fit_w_450_02_down,std_fit_w_450_05_down,std_fit_w_450_08_down, std_fit_w_450_1_down]

std_9k_up = [std_fit_w_9k_02_up,std_fit_w_9k_05_up,std_fit_w_9k_08_up, std_fit_w_9k_1_up]
std_9k_down = [std_fit_w_9k_02_down,std_fit_w_9k_05_down,std_fit_w_9k_08_down, std_fit_w_9k_1_down]

gen_index = {0:0,10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
color_label = ['turquoise','teal']
s = [0.02,0.05,0.08,0.1]

fig , ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True,figsize=(14,10),dpi=100, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.2, wspace=0.1)
for rep in range(4):
    ax=plt.subplot(2,2,rep+1)
    plt.plot(np.arange(0,(len(gen_index))), mean_450[rep], color = 'turquoise',linewidth = 3, label = '450')
    plt.fill_between(np.arange(0,(len(gen_index))), std_450_up[rep], std_450_down[rep], color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(np.arange(0,(len(gen_index))), mean_9k[rep], color = 'teal',linewidth = 3, linestyle = '--', label = '9000')
    plt.fill_between(np.arange(0,(len(gen_index))), std_9k_up[rep], std_9k_down[rep], color = 'teal', alpha = 0.3, edgecolor = None, linewidth=0.0)
    ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
    ax.tick_params(axis='y', which='both' ,right='off', left='on')
    plt.xlim(-0.2,(len(gen_index))-0.8) 
    plt.ylim(0,)
    plt.yticks(fontsize=12, va='top')    
    plt.xticks(np.arange(0,len(gen_index),1), [i for i in sorted(gen_index.keys())],fontsize=12, va='top')
    ax.set_title("s = %.2f" %(s[rep]),fontsize=12)
    if rep == 1:
        ax.legend(bbox_to_anchor=(0.65, 1), loc=3,ncol=2, borderaxespad=0.)    
    if rep == 0 or rep == 2:
        ax.set_ylabel('Fitness', size = 14,labelpad=15)
    if rep == 2 or rep == 3:
        plt.xlabel('Generation', size = 14,labelpad=15)
    if rep == 2 or rep == 3:plt.ylim(-0.02, 4)
plt.tight_layout()
plt.savefig('sweep_100loci_0.02to0.1s_log10fitness.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_100loci_0.02to0.1s_log10fitness.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''
############
#Figure S11
############
 '''

################################
#get rescaled phenotype for QTL
################################

#qff 450, 100 loci, different effect size: 0.04, 0.08, 0.2, 0.4
sim = 'sel_0.5-4.5--2.5-0.3'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/recom_rate_noX/1000sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 100, 0.04, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_450_04, std_fit_qtl_450_04_up, std_fit_qtl_450_04_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

sim = 'sel_0.5-4.5--6.1-0.3_eff0.08_100loci_movingTO'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 100, 0.08, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_450_08, std_fit_qtl_450_08_up, std_fit_qtl_450_08_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

sim = 'sel_0.5-4.5--16.9-0.3_eff0.2_100loci_movingTO'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 100, 0.2, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_450_2, std_fit_qtl_450_2_up, std_fit_qtl_450_2_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

sim = 'sel_0.5-4.5--34.9-0.3_eff0.4_100loci_movingTO'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 100, 0.4, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_450_4, std_fit_qtl_450_4_up, std_fit_qtl_450_4_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

#qff 9000, 100 loci, different effect size: 0.04, 0.08, 0.2, 0.4
sim = 'Ne9000_sel_0.5-4.5--2.5-0.3'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/recom_rate_noX/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 100, 0.04, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_9k_04, std_fit_qtl_9k_04_up, std_fit_qtl_9k_04_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

sim = 'Ne9000_sel_0.5-4.5--6.1-0.3_eff0.08_100loci_movingTO'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 100, 0.08, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_9k_08, std_fit_qtl_9k_08_up, std_fit_qtl_9k_08_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

sim = 'Ne9000_sel_0.5-4.5--16.9-0.3_eff0.2_100loci_movingTO'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 100, 0.2, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_9k_2, std_fit_qtl_9k_2_up, std_fit_qtl_9k_2_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

sim = 'Ne9000_sel_0.5-4.5--34.9-0.3_eff0.4_100loci_movingTO'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/GenoPheno.gpf' %sim,'r')  
loci_num, eff_size, p0 = 100, 0.4, 0.05
pop_opt = loci_num*(p0-(1-p0))*eff_size
mean_fit_qtl_9k_4, std_fit_qtl_9k_4_up, std_fit_qtl_9k_4_down = get_fitness_SimBatch_qff(InputFile_f,500,pop_opt)

##########
#plot
##
#######
    
mean_450 = [mean_fit_qtl_450_04,mean_fit_qtl_450_08,mean_fit_qtl_450_2, mean_fit_qtl_450_4]
mean_9k = [mean_fit_qtl_9k_04,mean_fit_qtl_9k_08,mean_fit_qtl_9k_2, mean_fit_qtl_9k_4]

std_450_up = [std_fit_qtl_450_04_up,std_fit_qtl_450_08_up,std_fit_qtl_450_2_up, std_fit_qtl_450_4_up]
std_450_down = [std_fit_qtl_450_04_down,std_fit_qtl_450_08_down,std_fit_qtl_450_2_down, std_fit_qtl_450_4_down]

std_9k_up = [std_fit_qtl_9k_04_up,std_fit_qtl_9k_08_up,std_fit_qtl_9k_2_up, std_fit_qtl_9k_4_up]
std_9k_down = [std_fit_qtl_9k_04_down,std_fit_qtl_9k_08_down,std_fit_qtl_9k_2_down, std_fit_qtl_9k_4_down]

gen_index = {0:0,10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
color_label = ['#e699b2','#cc3366']
s = [0.04,0.08,0.2,0.4]
TO = [item-item+1.1 for item in [-2.5, -6.1,-16.9,-34.9]]

fig , ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True,figsize=(14,10),dpi=100, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.2, wspace=0.1)
for rep in range(4):
    ax=plt.subplot(2,2,rep+1)
    plt.hlines(y = TO[rep],xmin = -0.2,xmax=(len(gen_index)*2)+len(gen_index), colors = 'gold', linestyles = 'solid',linewidth = 1)
    plt.plot(np.arange(0,(len(gen_index))), mean_450[rep], color = '#e699b2',linewidth = 3, label = '450')
    plt.fill_between(np.arange(0,(len(gen_index))), std_450_up[rep], std_450_down[rep], color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(np.arange(0,(len(gen_index))), mean_9k[rep], color = '#cc3366',linewidth = 3, linestyle = '--', label = '9000')
    plt.fill_between(np.arange(0,(len(gen_index))), std_9k_up[rep], std_9k_down[rep], color = '#cc3366', alpha = 0.3, edgecolor = None, linewidth=0.0)
    ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
    ax.tick_params(axis='y', which='both' ,right='off', left='on')
    plt.xlim(-0.2,(len(gen_index))-0.8) 
    plt.ylim(-0.05,)
    plt.yticks(fontsize=12, va='top')    
    plt.xticks(np.arange(0,len(gen_index),1), [i for i in sorted(gen_index.keys())],fontsize=12, va='top')
    ax.set_title("effect size = %.2f" %(s[rep]),fontsize=12)
    if rep == 1:
        ax.legend(bbox_to_anchor=(0.653, 1), loc=3,ncol=2, borderaxespad=0.)    
    if rep == 0 or rep == 2:
        ax.set_ylabel('Phenotype', size = 14,labelpad=15)
    if rep == 2 or rep == 3:
        plt.xlabel('Generation', size = 14,labelpad=15)
plt.tight_layout()
plt.savefig('qtl_100loci_0.04to0.4effectsize_phenotypeScaled.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('qtl_100loci_0.04to0.4effectsize_phenotypeScaled.pdf', dpi=300,format='pdf', bbox_inches = 'tight')
