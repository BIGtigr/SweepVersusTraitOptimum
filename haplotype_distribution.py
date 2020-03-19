# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 11:22:08 2020

@author: Neda
"""
import gzip
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#get info about the selected sites as a single haplotype for sweep simulations
def get_selected_haplo_w(infile):
    sel_haplo =''
    linenum = 0
    for line in infile:
        linenum += 1
        line = line.rstrip()
        if linenum == 1:
            header = line
        else:
            cols = line.split('\t')        
            sel_haplo+=(cols[2].split('/')[1])
    return sel_haplo
    infile.close()  

#get info about the selected sites as a single haplotype for quantitative trait simulations
def get_selected_haplo_qtl(infile):
    sel_haplo_qff =''
    for line in infile:
        line = line.rstrip()
        cols = line.split('\t')        
        sel_haplo_qff+=(cols[2].split('/')[0])
    return sel_haplo_qff
    infile.close()  

#extract haplotype information from the haplotypes stored during simulations
def get_haplo_info(path, num_sims, num_haplos, gen):
    haplo_all = []
    for filename in glob.glob(path):
        if gen in filename:
            with gzip.open(filename, 'rb') as f:
                print 'processing',filename
                haplo = ['' for _ in range(2*num_haplos)]
                linenum = 0
                for line in f:
                    linenum += 1
                    line=line.rstrip()
                    cols=line.split('\t')
                    if linenum == 1 : header = line
                    if linenum > 1:
                        all_haplo_tmp=cols[4].replace(' ','') 
                        for ind,h in enumerate(all_haplo_tmp):
                            haplo[ind]+=h
                haplo_all.append(haplo)
            f.close() 
    return haplo_all

#count the number of selected loci in haplotypes 
def compare_score_simBatch(ref_haplo, com_haplo, num_loci):
    score_all = []
    for s in com_haplo:
        tmp_scored_haplo = []
        for h in s:
            score = 0
            for i in range(num_loci):
                if h[i] == ref_haplo[i]:
                    score+=1
            tmp_scored_haplo.append(score)
        score_all.append(tmp_scored_haplo)
    return score_all

#this function distribute the values of each replicate into bins and gets the mean and standard deviation of all replicates for each bin 
# bins is a list or array of bin boundaries
def mean_binned_HaploScore(score_base,bins):
    #bin the haplotype score of all individuals for all replicates
    hist_score = []
    for ind,rep in enumerate(score_base):
        hist_score.append([rep.count(i)/float(len(rep)) for i in bins])
    #for each bin transpose the data so that for each bin a total of replicates item will be present 
    trans_hist_score = np.transpose(hist_score)   
    #then for each timepoint get the mean for each bin
    mean_hist_score=[np.mean(b) for b in trans_hist_score]
    #compute standard deviation of the mean for each bin
    std_hist_score=[np.std(b) for b in trans_hist_score]
    #compute the area around the mean for plotting std
    std_hist_score_up, std_hist_score_down = [], []
    for ind,b in enumerate(mean_hist_score):
        std_hist_score_down.append(b-(std_hist_score[ind]/2))
        std_hist_score_up.append(b+(std_hist_score[ind]/2))
    return (mean_hist_score,std_hist_score_up,std_hist_score_down)   
    
'''
##########
#Figure 4
##########
'''
#sweep 450
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/sel.txt','r')   
selected_haplo_w_450 = get_selected_haplo_w(InputFile)

#sweep 9000
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.08s/sel.txt','r')   
selected_haplo_w_9k = get_selected_haplo_w(InputFile)

#get the haplotypes for F0, F40, F80, F140 for Ne450 
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/recom_rate_noX/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=450

gen='g0'
haplo_base_w_450 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_450 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_w_450 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_450 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_450 = compare_score_simBatch(selected_haplo_w_450, haplo_base_w_450, 100)                
score_evol40_w_450 = compare_score_simBatch(selected_haplo_w_450, haplo_evol40_w_450, 100)                
score_evol80_w_450 = compare_score_simBatch(selected_haplo_w_450, haplo_evol80_w_450, 100)                
score_evol140_w_450 = compare_score_simBatch(selected_haplo_w_450, haplo_evol140_w_450, 100) 

mean_score_base_w_450, std_score_base_w_450_up, std_score_base_w_450_down = mean_binned_HaploScore(score_base_w_450,bins=np.arange(0,101,1))
mean_score_evol40_w_450, std_score_evol40_w_450_up, std_score_evol40_w_450_down = mean_binned_HaploScore(score_evol40_w_450,bins=np.arange(0,101,1))
mean_score_evol80_w_450, std_score_evol80_w_450_up, std_score_evol80_w_450_down = mean_binned_HaploScore(score_evol80_w_450,bins=np.arange(0,101,1))
mean_score_evol140_w_450, std_score_evol140_w_450_up, std_score_evol140_w_450_down = mean_binned_HaploScore(score_evol140_w_450,bins=np.arange(0,101,1))

#sweep Ne9000
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.08s/recom_rate_noX/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=9000

gen='g0'
haplo_base_w_9k = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_9k = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_w_9k = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_9k = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_9k = compare_score_simBatch(selected_haplo_w_9k, haplo_base_w_9k, 100)                
score_evol40_w_9k = compare_score_simBatch(selected_haplo_w_9k, haplo_evol40_w_9k, 100)                
score_evol80_w_9k = compare_score_simBatch(selected_haplo_w_9k, haplo_evol80_w_9k, 100)                
score_evol140_w_9k = compare_score_simBatch(selected_haplo_w_9k, haplo_evol140_w_9k, 100) 

mean_score_base_w_9k, std_score_base_w_9k_up, std_score_base_w_9k_down = mean_binned_HaploScore(score_base_w_9k,bins=np.arange(0,101,1))
mean_score_evol40_w_9k, std_score_evol40_w_9k_up, std_score_evol40_w_9k_down = mean_binned_HaploScore(score_evol40_w_9k,bins=np.arange(0,101,1))
mean_score_evol80_w_9k, std_score_evol80_w_9k_up, std_score_evol80_w_9k_down = mean_binned_HaploScore(score_evol80_w_9k,bins=np.arange(0,101,1))
mean_score_evol140_w_9k, std_score_evol140_w_9k_up, std_score_evol140_w_9k_down = mean_binned_HaploScore(score_evol140_w_9k,bins=np.arange(0,101,1))

#get selected sites as a single haplotype
#QTL 450
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/effect_size.txt','r')   
selected_haplo_qtl_450 = get_selected_haplo_qtl(InputFile)

#QTL 9000
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--2.5-0.3/effect_size.txt','r')   
selected_haplo_qtl_9k = get_selected_haplo_qtl(InputFile)

#get the haplotypes for F0, F40, F140 for Ne450
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/recom_rate_noX/2ndRun/evolved_haplos/*"
num_sims=50
num_haplos=450

gen='g0'
haplo_base_qtl_450 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_qtl_450 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_qtl_450 = get_haplo_info(path, num_sims, num_haplos, gen)

#now count the selected loci in haplos 
score_base_qtl_450 = compare_score_simBatch(selected_haplo_qtl_450, haplo_base_qtl_450, 100)                
score_evol40_qtl_450 = compare_score_simBatch(selected_haplo_qtl_450, haplo_evol40_qtl_450, 100)                
score_evol140_qtl_450 = compare_score_simBatch(selected_haplo_qtl_450, haplo_evol140_qtl_450, 100)                

mean_score_base_qtl_450, std_score_base_qtl_450_up, std_score_base_qtl_450_down = mean_binned_HaploScore(score_base_qtl_450,bins=np.arange(0,101,1))
mean_score_evol40_qtl_450, std_score_evol40_qtl_450_up, std_score_evol40_qtl_450_down = mean_binned_HaploScore(score_evol40_qtl_450,bins=np.arange(0,101,1))
mean_score_evol140_qtl_450, std_score_evol140_qtl_450_up, std_score_evol140_qtl_450_down = mean_binned_HaploScore(score_evol140_qtl_450,bins=np.arange(0,101,1))

#get the haplotypes for F0, F40, F140 for Ne9000
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--2.5-0.3/recom_rate_noX/2ndRun/evolved_haplos/*"
num_sims=50
num_haplos=9000

gen='g0'
haplo_base_qtl_9k = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_qtl_9k = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_qtl_9k = get_haplo_info(path, num_sims, num_haplos, gen)

#now count the selected loci in haplos 
score_base_qtl_9k = compare_score_simBatch(selected_haplo_qtl_9k, haplo_base_qtl_9k, 100)                
score_evol40_qtl_9k = compare_score_simBatch(selected_haplo_qtl_9k, haplo_evol40_qtl_9k, 100)                
score_evol140_qtl_9k = compare_score_simBatch(selected_haplo_qtl_9k, haplo_evol140_qtl_9k, 100)                
    
mean_score_base_qtl_9k, std_score_base_qtl_9k_up, std_score_base_qtl_9k_down = mean_binned_HaploScore(score_base_qtl_9k,bins=np.arange(0,101,1))
mean_score_evol40_qtl_9k, std_score_evol40_qtl_9k_up, std_score_evol40_qtl_9k_down = mean_binned_HaploScore(score_evol40_qtl_9k,bins=np.arange(0,101,1))
mean_score_evol140_qtl_9k, std_score_evol140_qtl_9k_up, std_score_evol140_qtl_9k_down = mean_binned_HaploScore(score_evol140_qtl_9k,bins=np.arange(0,101,1))
    
########
#plot
########

base = [mean_score_base_w_450,mean_score_base_w_9k,mean_score_base_qtl_450,mean_score_base_qtl_9k]
F40 = [mean_score_evol40_w_450,mean_score_evol40_w_9k,mean_score_evol40_qtl_450,mean_score_evol40_qtl_9k]
F80 = [mean_score_evol80_w_450,mean_score_evol80_w_9k]
F140 = [mean_score_evol140_w_450,mean_score_evol140_w_9k,mean_score_evol140_qtl_450,mean_score_evol140_qtl_9k]

base_up = [std_score_base_w_450_up,std_score_base_w_9k_up,std_score_base_qtl_450_up,std_score_base_qtl_9k_up]
F40_up = [std_score_evol40_w_450_up,std_score_evol40_w_9k_up,std_score_evol40_qtl_450_up,std_score_evol40_qtl_9k_up]
F80_up = [std_score_evol80_w_450_up,std_score_evol80_w_9k_up]
F140_up = [std_score_evol140_w_450_up,std_score_evol140_w_9k_up,std_score_evol140_qtl_450_up,std_score_evol140_qtl_9k_up]

base_down = [std_score_base_w_450_down,std_score_base_w_9k_down,std_score_base_qtl_450_down,std_score_base_qtl_9k_down]
F40_down = [std_score_evol40_w_450_down,std_score_evol40_w_9k_down,std_score_evol40_qtl_450_down,std_score_evol40_qtl_9k_down]
F80_down = [std_score_evol80_w_450_down,std_score_evol80_w_9k_down]
F140_down = [std_score_evol140_w_450_down,std_score_evol140_w_9k_down,std_score_evol140_qtl_450_down,std_score_evol140_qtl_9k_down]

model = ['sweep','sweep','trait optimum','trait optimum']
pop_size = [450,9000,450,9000]
color_label = ['#7fcdbb','#41b6c4','#1d91c0','#253494']
x = np.arange(0,101,1)
 
fig , ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True,figsize=(18,10),dpi=100, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.2, wspace=0.2)
for rep in range(4):
    ax=plt.subplot(2,2,rep+1)
    plt.plot(x,base[rep],color = color_label[0], linewidth = 3)#, alpha = 0.6
    plt.fill_between(x, base_up[rep], base_down[rep], color = color_label[0], alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(x,F40[rep],color = color_label[1], linewidth = 3)
    plt.fill_between(x, F40_up[rep], F40_down[rep], color = color_label[1], alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(x,F140[rep],color = color_label[3], linewidth = 3)
    plt.fill_between(x, F140_up[rep], F140_down[rep], color = color_label[3], alpha = 0.3, edgecolor = None, linewidth=0.0)
    ax.tick_params(axis='x', which='both' ,bottom='off', top='off')
    ax.tick_params(axis='y', which='both' ,right='off', left='on') 
    ax.set_title("%s, N = %i" %(model[rep],pop_size[rep]),fontsize=12)
    plt.xticks(fontsize=12, va='top')
    plt.yticks(fontsize=12, va='top') 
    plt.ylim(0,0.25)
    if rep == 0 or rep == 1:
        plt.plot(x,F80[rep],color = color_label[2], linewidth = 3)
        plt.fill_between(x, F80_up[rep], F80_down[rep], color = color_label[2], alpha = 0.3, edgecolor = None, linewidth=0.0)
    if rep == 1:
        texts = ['Gen 0','Gen 40','Gen 80','Gen 140']
        patches = [mpatches.Patch(color=color_label[i], label="{:s}".format(texts[i]) ) for i in range(len(texts))]
        ax.legend(handles=patches,loc='upper left',bbox_to_anchor=(0.77, 1),frameon=False)    
    if rep == 0 or rep == 2:
        ax.set_ylabel('Proportion', size = 12, labelpad=15)
    if rep == 2 or rep == 3:
        plt.xlabel('No. of selected loci per haplotype', size = 12, labelpad=15)        
plt.savefig('sweep_qtl_selectedLociPerHaplo_proportion_defaultParameters.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_qtl_selectedLociPerHaplo_proportion_defaultParameters.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''
##########
#Figure 9
##########   
'''

###########
#sweep 450
###########

#sweep 450, 0.02 s
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.02s/sel.txt','r')   
selected_haplo_w_450_02 = get_selected_haplo_w(InputFile)

#get the haplotypes for F0, F20, F40, F80, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.02s/recom_rate_noX/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=450

gen='g0'
haplo_base_w_450_02 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_w_450_02 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_450_02 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_w_450_02 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_450_02 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_450_02 = compare_score_simBatch(selected_haplo_w_450_02, haplo_base_w_450_02, 100)                
score_evol20_w_450_02 = compare_score_simBatch(selected_haplo_w_450_02, haplo_evol20_w_450_02, 100)                
score_evol40_w_450_02 = compare_score_simBatch(selected_haplo_w_450_02, haplo_evol40_w_450_02, 100)                
score_evol80_w_450_02 = compare_score_simBatch(selected_haplo_w_450_02, haplo_evol80_w_450_02, 100) 
score_evol140_w_450_02 = compare_score_simBatch(selected_haplo_w_450_02, haplo_evol140_w_450_02, 100) 

mean_score_base_w_450_02, std_score_base_w_450_02_up, std_score_base_w_450_02_down = mean_binned_HaploScore(score_base_w_450_02,bins=np.arange(0,101,1))
mean_score_evol20_w_450_02, std_score_evol20_w_450_02_up, std_score_evol20_w_450_02_down = mean_binned_HaploScore(score_evol20_w_450_02,bins=np.arange(0,101,1))
mean_score_evol40_w_450_02, std_score_evol40_w_450_02_up, std_score_evol40_w_450_02_down = mean_binned_HaploScore(score_evol40_w_450_02,bins=np.arange(0,101,1))
mean_score_evol80_w_450_02, std_score_evol80_w_450_02_up, std_score_evol80_w_450_02_down = mean_binned_HaploScore(score_evol80_w_450_02,bins=np.arange(0,101,1))
mean_score_evol140_w_450_02, std_score_evol140_w_450_02_up, std_score_evol140_w_450_02_down = mean_binned_HaploScore(score_evol140_w_450_02,bins=np.arange(0,101,1))

#sweep 450, 0.05 s
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.05s/sel.txt','r')   
selected_haplo_w_450_05 = get_selected_haplo_w(InputFile)

#get the haplotypes for F0, F20, F40, F80, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.05s/recom_rate_noX/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=450

gen='g0'
haplo_base_w_450_05 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_w_450_05 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_450_05 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_w_450_05 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_450_05 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_450_05 = compare_score_simBatch(selected_haplo_w_450_05, haplo_base_w_450_05, 100)                
score_evol20_w_450_05 = compare_score_simBatch(selected_haplo_w_450_05, haplo_evol20_w_450_05, 100)                
score_evol40_w_450_05 = compare_score_simBatch(selected_haplo_w_450_05, haplo_evol40_w_450_05, 100)                
score_evol80_w_450_05 = compare_score_simBatch(selected_haplo_w_450_05, haplo_evol80_w_450_05, 100) 
score_evol140_w_450_05 = compare_score_simBatch(selected_haplo_w_450_05, haplo_evol140_w_450_05, 100) 

mean_score_base_w_450_05, std_score_base_w_450_05_up, std_score_base_w_450_05_down = mean_binned_HaploScore(score_base_w_450_05,bins=np.arange(0,101,1))
mean_score_evol20_w_450_05, std_score_evol20_w_450_05_up, std_score_evol20_w_450_05_down = mean_binned_HaploScore(score_evol20_w_450_05,bins=np.arange(0,101,1))
mean_score_evol40_w_450_05, std_score_evol40_w_450_05_up, std_score_evol40_w_450_05_down = mean_binned_HaploScore(score_evol40_w_450_05,bins=np.arange(0,101,1))
mean_score_evol80_w_450_05, std_score_evol80_w_450_05_up, std_score_evol80_w_450_05_down = mean_binned_HaploScore(score_evol80_w_450_05,bins=np.arange(0,101,1))
mean_score_evol140_w_450_05, std_score_evol140_w_450_05_up, std_score_evol140_w_450_05_down = mean_binned_HaploScore(score_evol140_w_450_05,bins=np.arange(0,101,1))

#sweep 450, 0.08 s
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/sel.txt','r')   
selected_haplo_w_450_08 = get_selected_haplo_w(InputFile)

#get the haplotypes for F0, F20, F40, F80, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/recom_rate_noX/4thRun/evolved_haplos/*"

num_sims=50
num_haplos=450

gen='g0'
haplo_base_w_450_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_w_450_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_450_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_w_450_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_450_08 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_450_08 = compare_score_simBatch(selected_haplo_w_450_08, haplo_base_w_450_08, 100)                
score_evol20_w_450_08 = compare_score_simBatch(selected_haplo_w_450_08, haplo_evol20_w_450_08, 100)                
score_evol40_w_450_08 = compare_score_simBatch(selected_haplo_w_450_08, haplo_evol40_w_450_08, 100)                
score_evol80_w_450_08 = compare_score_simBatch(selected_haplo_w_450_08, haplo_evol80_w_450_08, 100) 
score_evol140_w_450_08 = compare_score_simBatch(selected_haplo_w_450_08, haplo_evol140_w_450_08, 100) 

mean_score_base_w_450_08, std_score_base_w_450_08_up, std_score_base_w_450_08_down = mean_binned_HaploScore(score_base_w_450_08,bins=np.arange(0,101,1))
mean_score_evol20_w_450_08, std_score_evol20_w_450_08_up, std_score_evol20_w_450_08_down = mean_binned_HaploScore(score_evol20_w_450_08,bins=np.arange(0,101,1))
mean_score_evol40_w_450_08, std_score_evol40_w_450_08_up, std_score_evol40_w_450_08_down = mean_binned_HaploScore(score_evol40_w_450_08,bins=np.arange(0,101,1))
mean_score_evol80_w_450_08, std_score_evol80_w_450_08_up, std_score_evol80_w_450_08_down = mean_binned_HaploScore(score_evol80_w_450_08,bins=np.arange(0,101,1))
mean_score_evol140_w_450_08, std_score_evol140_w_450_08_up, std_score_evol140_w_450_08_down = mean_binned_HaploScore(score_evol140_w_450_08,bins=np.arange(0,101,1))

#sweep 450, 0.1 s
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.1s/sel.txt','r')   
selected_haplo_w_450_1 = get_selected_haplo_w(InputFile)

#get the haplotypes for F0, F20, F40, F80, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.1s/recom_rate_noX/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=450

gen='g0'
haplo_base_w_450_1 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_w_450_1 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_450_1 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_w_450_1 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_450_1 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_450_1 = compare_score_simBatch(selected_haplo_w_450_1, haplo_base_w_450_1, 100)                
score_evol20_w_450_1 = compare_score_simBatch(selected_haplo_w_450_1, haplo_evol20_w_450_1, 100)                
score_evol40_w_450_1 = compare_score_simBatch(selected_haplo_w_450_1, haplo_evol40_w_450_1, 100)                
score_evol80_w_450_1 = compare_score_simBatch(selected_haplo_w_450_1, haplo_evol80_w_450_1, 100) 
score_evol140_w_450_1 = compare_score_simBatch(selected_haplo_w_450_1, haplo_evol140_w_450_1, 100) 

mean_score_base_w_450_1, std_score_base_w_450_1_up, std_score_base_w_450_1_down = mean_binned_HaploScore(score_base_w_450_1,bins=np.arange(0,101,1))
mean_score_evol20_w_450_1, std_score_evol20_w_450_1_up, std_score_evol20_w_450_1_down = mean_binned_HaploScore(score_evol20_w_450_1,bins=np.arange(0,101,1))
mean_score_evol40_w_450_1, std_score_evol40_w_450_1_up, std_score_evol40_w_450_1_down = mean_binned_HaploScore(score_evol40_w_450_1,bins=np.arange(0,101,1))
mean_score_evol80_w_450_1, std_score_evol80_w_450_1_up, std_score_evol80_w_450_1_down = mean_binned_HaploScore(score_evol80_w_450_1,bins=np.arange(0,101,1))
mean_score_evol140_w_450_1, std_score_evol140_w_450_1_up, std_score_evol140_w_450_1_down = mean_binned_HaploScore(score_evol140_w_450_1,bins=np.arange(0,101,1))

#########
#sweep 9k
#########

#sweep 9k, 0.02 s
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.02s/sel.txt','r')   
selected_haplo_w_9k_02 = get_selected_haplo_w(InputFile)

#get the haplotypes for F0, F20, F40, F80, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.02s/recom_rate_noX/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=9000

gen='g0'
haplo_base_w_9k_02 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_w_9k_02 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_9k_02 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_w_9k_02 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_9k_02 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_9k_02 = compare_score_simBatch(selected_haplo_w_9k_02, haplo_base_w_9k_02, 100)                
score_evol20_w_9k_02 = compare_score_simBatch(selected_haplo_w_9k_02, haplo_evol20_w_9k_02, 100)                
score_evol40_w_9k_02 = compare_score_simBatch(selected_haplo_w_9k_02, haplo_evol40_w_9k_02, 100)                
score_evol80_w_9k_02 = compare_score_simBatch(selected_haplo_w_9k_02, haplo_evol80_w_9k_02, 100) 
score_evol140_w_9k_02 = compare_score_simBatch(selected_haplo_w_9k_02, haplo_evol140_w_9k_02, 100) 

mean_score_base_w_9k_02, std_score_base_w_9k_02_up, std_score_base_w_9k_02_down = mean_binned_HaploScore(score_base_w_9k_02,bins=np.arange(0,101,1))
mean_score_evol20_w_9k_02, std_score_evol20_w_9k_02_up, std_score_evol20_w_9k_02_down = mean_binned_HaploScore(score_evol20_w_9k_02,bins=np.arange(0,101,1))
mean_score_evol40_w_9k_02, std_score_evol40_w_9k_02_up, std_score_evol40_w_9k_02_down = mean_binned_HaploScore(score_evol40_w_9k_02,bins=np.arange(0,101,1))
mean_score_evol80_w_9k_02, std_score_evol80_w_9k_02_up, std_score_evol80_w_9k_02_down = mean_binned_HaploScore(score_evol80_w_9k_02,bins=np.arange(0,101,1))
mean_score_evol140_w_9k_02, std_score_evol140_w_9k_02_up, std_score_evol140_w_9k_02_down = mean_binned_HaploScore(score_evol140_w_9k_02,bins=np.arange(0,101,1))

#sweep 9k, 0.05 s
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.05s/sel.txt','r')   
selected_haplo_w_9k_05 = get_selected_haplo_w(InputFile)

#get the haplotypes for F0, F20, F40, F80, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.05s/recom_rate_noX/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=9000

gen='g0'
haplo_base_w_9k_05 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_w_9k_05 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_9k_05 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_w_9k_05 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_9k_05 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_9k_05 = compare_score_simBatch(selected_haplo_w_9k_05, haplo_base_w_9k_05, 100)                
score_evol20_w_9k_05 = compare_score_simBatch(selected_haplo_w_9k_05, haplo_evol20_w_9k_05, 100)                
score_evol40_w_9k_05 = compare_score_simBatch(selected_haplo_w_9k_05, haplo_evol40_w_9k_05, 100)                
score_evol80_w_9k_05 = compare_score_simBatch(selected_haplo_w_9k_05, haplo_evol80_w_9k_05, 100) 
score_evol140_w_9k_05 = compare_score_simBatch(selected_haplo_w_9k_05, haplo_evol140_w_9k_05, 100) 

mean_score_base_w_9k_05, std_score_base_w_9k_05_up, std_score_base_w_9k_05_down = mean_binned_HaploScore(score_base_w_9k_05,bins=np.arange(0,101,1))
mean_score_evol20_w_9k_05, std_score_evol20_w_9k_05_up, std_score_evol20_w_9k_05_down = mean_binned_HaploScore(score_evol20_w_9k_05,bins=np.arange(0,101,1))
mean_score_evol40_w_9k_05, std_score_evol40_w_9k_05_up, std_score_evol40_w_9k_05_down = mean_binned_HaploScore(score_evol40_w_9k_05,bins=np.arange(0,101,1))
mean_score_evol80_w_9k_05, std_score_evol80_w_9k_05_up, std_score_evol80_w_9k_05_down = mean_binned_HaploScore(score_evol80_w_9k_05,bins=np.arange(0,101,1))
mean_score_evol140_w_9k_05, std_score_evol140_w_9k_05_up, std_score_evol140_w_9k_05_down = mean_binned_HaploScore(score_evol140_w_9k_05,bins=np.arange(0,101,1))

#sweep 9k, 0.08 s
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.08s/sel.txt','r')   
selected_haplo_w_9k_08 = get_selected_haplo_w(InputFile)

#get the haplotypes for F0, F20, F40, F80, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.08s/recom_rate_noX/4thRun/evolved_haplos/*"

num_sims=50
num_haplos=9000

gen='g0'
haplo_base_w_9k_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_w_9k_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_9k_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_w_9k_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_9k_08 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_9k_08 = compare_score_simBatch(selected_haplo_w_9k_08, haplo_base_w_9k_08, 100)                
score_evol20_w_9k_08 = compare_score_simBatch(selected_haplo_w_9k_08, haplo_evol20_w_9k_08, 100)                
score_evol40_w_9k_08 = compare_score_simBatch(selected_haplo_w_9k_08, haplo_evol40_w_9k_08, 100)                
score_evol80_w_9k_08 = compare_score_simBatch(selected_haplo_w_9k_08, haplo_evol80_w_9k_08, 100) 
score_evol140_w_9k_08 = compare_score_simBatch(selected_haplo_w_9k_08, haplo_evol140_w_9k_08, 100) 

mean_score_base_w_9k_08, std_score_base_w_9k_08_up, std_score_base_w_9k_08_down = mean_binned_HaploScore(score_base_w_9k_08,bins=np.arange(0,101,1))
mean_score_evol20_w_9k_08, std_score_evol20_w_9k_08_up, std_score_evol20_w_9k_08_down = mean_binned_HaploScore(score_evol20_w_9k_08,bins=np.arange(0,101,1))
mean_score_evol40_w_9k_08, std_score_evol40_w_9k_08_up, std_score_evol40_w_9k_08_down = mean_binned_HaploScore(score_evol40_w_9k_08,bins=np.arange(0,101,1))
mean_score_evol80_w_9k_08, std_score_evol80_w_9k_08_up, std_score_evol80_w_9k_08_down = mean_binned_HaploScore(score_evol80_w_9k_08,bins=np.arange(0,101,1))
mean_score_evol140_w_9k_08, std_score_evol140_w_9k_08_up, std_score_evol140_w_9k_08_down = mean_binned_HaploScore(score_evol140_w_9k_08,bins=np.arange(0,101,1))

#sweep 9k, 0.1 s
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.1s/sel.txt','r')   
selected_haplo_w_9k_1 = get_selected_haplo_w(InputFile)

#get the haplotypes for F0, F20, F40, F80, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.1s/recom_rate_noX/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=9000

gen='g0'
haplo_base_w_9k_1 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_w_9k_1 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_9k_1 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_w_9k_1 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_9k_1 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_9k_1 = compare_score_simBatch(selected_haplo_w_9k_1, haplo_base_w_9k_1, 100)                
score_evol20_w_9k_1 = compare_score_simBatch(selected_haplo_w_9k_1, haplo_evol20_w_9k_1, 100)                
score_evol40_w_9k_1 = compare_score_simBatch(selected_haplo_w_9k_1, haplo_evol40_w_9k_1, 100)                
score_evol80_w_9k_1 = compare_score_simBatch(selected_haplo_w_9k_1, haplo_evol80_w_9k_1, 100) 
score_evol140_w_9k_1 = compare_score_simBatch(selected_haplo_w_9k_1, haplo_evol140_w_9k_1, 100) 

mean_score_base_w_9k_1, std_score_base_w_9k_1_up, std_score_base_w_9k_1_down = mean_binned_HaploScore(score_base_w_9k_1,bins=np.arange(0,101,1))
mean_score_evol20_w_9k_1, std_score_evol20_w_9k_1_up, std_score_evol20_w_9k_1_down = mean_binned_HaploScore(score_evol20_w_9k_1,bins=np.arange(0,101,1))
mean_score_evol40_w_9k_1, std_score_evol40_w_9k_1_up, std_score_evol40_w_9k_1_down = mean_binned_HaploScore(score_evol40_w_9k_1,bins=np.arange(0,101,1))
mean_score_evol80_w_9k_1, std_score_evol80_w_9k_1_up, std_score_evol80_w_9k_1_down = mean_binned_HaploScore(score_evol80_w_9k_1,bins=np.arange(0,101,1))
mean_score_evol140_w_9k_1, std_score_evol140_w_9k_1_up, std_score_evol140_w_9k_1_down = mean_binned_HaploScore(score_evol140_w_9k_1,bins=np.arange(0,101,1))

###########
#QTL 450
###########

#QTL 450, 0.04 effect size
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/effect_size.txt','r')   
selected_haplo_qtl_450_04 = get_selected_haplo_qtl(InputFile)

#get the haplotypes for F0, F40, F70, F110, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/recom_rate_noX/4thRun/evolved_haplos/*"

num_sims=50
num_haplos=450

gen='g0'
haplo_base_qtl_450_04 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_qtl_450_04 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_qtl_450_04 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_qtl_450_04 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_qtl_450_04 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_qtl_450_04 = compare_score_simBatch(selected_haplo_qtl_450_04, haplo_base_qtl_450_04, 100)                
score_evol20_qtl_450_04 = compare_score_simBatch(selected_haplo_qtl_450_04, haplo_evol20_qtl_450_04, 100)                
score_evol40_qtl_450_04 = compare_score_simBatch(selected_haplo_qtl_450_04, haplo_evol40_qtl_450_04, 100)                
score_evol80_qtl_450_04 = compare_score_simBatch(selected_haplo_qtl_450_04, haplo_evol80_qtl_450_04, 100) 
score_evol140_qtl_450_04 = compare_score_simBatch(selected_haplo_qtl_450_04, haplo_evol140_qtl_450_04, 100) 

mean_score_base_qtl_450_04, std_score_base_qtl_450_04_up, std_score_base_qtl_450_04_down = mean_binned_HaploScore(score_base_qtl_450_04,bins=np.arange(0,101,1))
mean_score_evol20_qtl_450_04, std_score_evol20_qtl_450_04_up, std_score_evol20_qtl_450_04_down = mean_binned_HaploScore(score_evol20_qtl_450_04,bins=np.arange(0,101,1))
mean_score_evol40_qtl_450_04, std_score_evol40_qtl_450_04_up, std_score_evol40_qtl_450_04_down = mean_binned_HaploScore(score_evol40_qtl_450_04,bins=np.arange(0,101,1))
mean_score_evol80_qtl_450_04, std_score_evol80_qtl_450_04_up, std_score_evol80_qtl_450_04_down = mean_binned_HaploScore(score_evol80_qtl_450_04,bins=np.arange(0,101,1))
mean_score_evol140_qtl_450_04, std_score_evol140_qtl_450_04_up, std_score_evol140_qtl_450_04_down = mean_binned_HaploScore(score_evol140_qtl_450_04,bins=np.arange(0,101,1))

#QTL 450, 0.08 effect size
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--6.1-0.3_eff0.08_100loci_movingTO/effect_size.txt','r')   
selected_haplo_qtl_450_08 = get_selected_haplo_qtl(InputFile)

#get the haplotypes for F0, F40, F70, F110, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--6.1-0.3_eff0.08_100loci_movingTO/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=450

gen='g0'
haplo_base_qtl_450_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_qtl_450_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_qtl_450_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_qtl_450_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_qtl_450_08 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_qtl_450_08 = compare_score_simBatch(selected_haplo_qtl_450_08, haplo_base_qtl_450_08, 100)                
score_evol20_qtl_450_08 = compare_score_simBatch(selected_haplo_qtl_450_08, haplo_evol20_qtl_450_08, 100)                
score_evol40_qtl_450_08 = compare_score_simBatch(selected_haplo_qtl_450_08, haplo_evol40_qtl_450_08, 100)                
score_evol80_qtl_450_08 = compare_score_simBatch(selected_haplo_qtl_450_08, haplo_evol80_qtl_450_08, 100) 
score_evol140_qtl_450_08 = compare_score_simBatch(selected_haplo_qtl_450_08, haplo_evol140_qtl_450_08, 100) 

mean_score_base_qtl_450_08, std_score_base_qtl_450_08_up, std_score_base_qtl_450_08_down = mean_binned_HaploScore(score_base_qtl_450_08,bins=np.arange(0,101,1))
mean_score_evol20_qtl_450_08, std_score_evol20_qtl_450_08_up, std_score_evol20_qtl_450_08_down = mean_binned_HaploScore(score_evol20_qtl_450_08,bins=np.arange(0,101,1))
mean_score_evol40_qtl_450_08, std_score_evol40_qtl_450_08_up, std_score_evol40_qtl_450_08_down = mean_binned_HaploScore(score_evol40_qtl_450_08,bins=np.arange(0,101,1))
mean_score_evol80_qtl_450_08, std_score_evol80_qtl_450_08_up, std_score_evol80_qtl_450_08_down = mean_binned_HaploScore(score_evol80_qtl_450_08,bins=np.arange(0,101,1))
mean_score_evol140_qtl_450_08, std_score_evol140_qtl_450_08_up, std_score_evol140_qtl_450_08_down = mean_binned_HaploScore(score_evol140_qtl_450_08,bins=np.arange(0,101,1))

#QTL 450, 0.2 effect size
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--16.9-0.3_eff0.2_100loci_movingTO/effect_size.txt','r')   
selected_haplo_qtl_450_2 = get_selected_haplo_qtl(InputFile)

#get the haplotypes for F0, F40, F70, F110, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--16.9-0.3_eff0.2_100loci_movingTO/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=450

gen='g0'
haplo_base_qtl_450_2 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_qtl_450_2 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_qtl_450_2 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_qtl_450_2 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_qtl_450_2 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_qtl_450_2 = compare_score_simBatch(selected_haplo_qtl_450_2, haplo_base_qtl_450_2, 100)                
score_evol20_qtl_450_2 = compare_score_simBatch(selected_haplo_qtl_450_2, haplo_evol20_qtl_450_2, 100)                
score_evol40_qtl_450_2 = compare_score_simBatch(selected_haplo_qtl_450_2, haplo_evol40_qtl_450_2, 100)                
score_evol80_qtl_450_2 = compare_score_simBatch(selected_haplo_qtl_450_2, haplo_evol80_qtl_450_2, 100) 
score_evol140_qtl_450_2 = compare_score_simBatch(selected_haplo_qtl_450_2, haplo_evol140_qtl_450_2, 100) 

mean_score_base_qtl_450_2, std_score_base_qtl_450_2_up, std_score_base_qtl_450_2_down = mean_binned_HaploScore(score_base_qtl_450_2,bins=np.arange(0,101,1))
mean_score_evol20_qtl_450_2, std_score_evol20_qtl_450_2_up, std_score_evol20_qtl_450_2_down = mean_binned_HaploScore(score_evol20_qtl_450_2,bins=np.arange(0,101,1))
mean_score_evol40_qtl_450_2, std_score_evol40_qtl_450_2_up, std_score_evol40_qtl_450_2_down = mean_binned_HaploScore(score_evol40_qtl_450_2,bins=np.arange(0,101,1))
mean_score_evol80_qtl_450_2, std_score_evol80_qtl_450_2_up, std_score_evol80_qtl_450_2_down = mean_binned_HaploScore(score_evol80_qtl_450_2,bins=np.arange(0,101,1))
mean_score_evol140_qtl_450_2, std_score_evol140_qtl_450_2_up, std_score_evol140_qtl_450_2_down = mean_binned_HaploScore(score_evol140_qtl_450_2,bins=np.arange(0,101,1))

#QTL 450, 0.4 effect size
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--34.9-0.3_eff0.4_100loci_movingTO/effect_size.txt','r')   
selected_haplo_qtl_450_4 = get_selected_haplo_qtl(InputFile)

#get the haplotypes for F0, F40, F70, F110, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--34.9-0.3_eff0.4_100loci_movingTO/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=450

gen='g0'
haplo_base_qtl_450_4 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_qtl_450_4 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_qtl_450_4 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_qtl_450_4 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_qtl_450_4 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_qtl_450_4 = compare_score_simBatch(selected_haplo_qtl_450_4, haplo_base_qtl_450_4, 100)                
score_evol20_qtl_450_4 = compare_score_simBatch(selected_haplo_qtl_450_4, haplo_evol20_qtl_450_4, 100)                
score_evol40_qtl_450_4 = compare_score_simBatch(selected_haplo_qtl_450_4, haplo_evol40_qtl_450_4, 100)                
score_evol80_qtl_450_4 = compare_score_simBatch(selected_haplo_qtl_450_4, haplo_evol80_qtl_450_4, 100) 
score_evol140_qtl_450_4 = compare_score_simBatch(selected_haplo_qtl_450_4, haplo_evol140_qtl_450_4, 100) 

mean_score_base_qtl_450_4, std_score_base_qtl_450_4_up, std_score_base_qtl_450_4_down = mean_binned_HaploScore(score_base_qtl_450_4,bins=np.arange(0,101,1))
mean_score_evol20_qtl_450_4, std_score_evol20_qtl_450_4_up, std_score_evol20_qtl_450_4_down = mean_binned_HaploScore(score_evol20_qtl_450_4,bins=np.arange(0,101,1))
mean_score_evol40_qtl_450_4, std_score_evol40_qtl_450_4_up, std_score_evol40_qtl_450_4_down = mean_binned_HaploScore(score_evol40_qtl_450_4,bins=np.arange(0,101,1))
mean_score_evol80_qtl_450_4, std_score_evol80_qtl_450_4_up, std_score_evol80_qtl_450_4_down = mean_binned_HaploScore(score_evol80_qtl_450_4,bins=np.arange(0,101,1))
mean_score_evol140_qtl_450_4, std_score_evol140_qtl_450_4_up, std_score_evol140_qtl_450_4_down = mean_binned_HaploScore(score_evol140_qtl_450_4,bins=np.arange(0,101,1))

#########
#QTL 9000
#########

#QTL 9k, 0.04 effect size
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--2.5-0.3/effect_size.txt','r')   
selected_haplo_qtl_9k_04 = get_selected_haplo_qtl(InputFile)

#get the haplotypes for F0, F40, F70, F110, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--2.5-0.3/recom_rate_noX/4thRun/evolved_haplos/*"

num_sims=50
num_haplos=9000

gen='g0'
haplo_base_qtl_9k_04 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_qtl_9k_04 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_qtl_9k_04 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_qtl_9k_04 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_qtl_9k_04 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_qtl_9k_04 = compare_score_simBatch(selected_haplo_qtl_9k_04, haplo_base_qtl_9k_04, 100)                
score_evol20_qtl_9k_04 = compare_score_simBatch(selected_haplo_qtl_9k_04, haplo_evol20_qtl_9k_04, 100)                
score_evol40_qtl_9k_04 = compare_score_simBatch(selected_haplo_qtl_9k_04, haplo_evol40_qtl_9k_04, 100)                
score_evol80_qtl_9k_04 = compare_score_simBatch(selected_haplo_qtl_9k_04, haplo_evol80_qtl_9k_04, 100) 
score_evol140_qtl_9k_04 = compare_score_simBatch(selected_haplo_qtl_9k_04, haplo_evol140_qtl_9k_04, 100) 

mean_score_base_qtl_9k_04, std_score_base_qtl_9k_04_up, std_score_base_qtl_9k_04_down = mean_binned_HaploScore(score_base_qtl_9k_04,bins=np.arange(0,101,1))
mean_score_evol20_qtl_9k_04, std_score_evol20_qtl_9k_04_up, std_score_evol20_qtl_9k_04_down = mean_binned_HaploScore(score_evol20_qtl_9k_04,bins=np.arange(0,101,1))
mean_score_evol40_qtl_9k_04, std_score_evol40_qtl_9k_04_up, std_score_evol40_qtl_9k_04_down = mean_binned_HaploScore(score_evol40_qtl_9k_04,bins=np.arange(0,101,1))
mean_score_evol80_qtl_9k_04, std_score_evol80_qtl_9k_04_up, std_score_evol80_qtl_9k_04_down = mean_binned_HaploScore(score_evol80_qtl_9k_04,bins=np.arange(0,101,1))
mean_score_evol140_qtl_9k_04, std_score_evol140_qtl_9k_04_up, std_score_evol140_qtl_9k_04_down = mean_binned_HaploScore(score_evol140_qtl_9k_04,bins=np.arange(0,101,1))

#QTL 9k, 0.08 effect size
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--6.1-0.3_eff0.08_100loci_movingTO/effect_size.txt','r')   
selected_haplo_qtl_9k_08 = get_selected_haplo_qtl(InputFile)

#get the haplotypes for F0, F40, F70, F110, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--6.1-0.3_eff0.08_100loci_movingTO/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=9000

gen='g0'
haplo_base_qtl_9k_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_qtl_9k_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_qtl_9k_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_qtl_9k_08 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_qtl_9k_08 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_qtl_9k_08 = compare_score_simBatch(selected_haplo_qtl_9k_08, haplo_base_qtl_9k_08, 100)                
score_evol20_qtl_9k_08 = compare_score_simBatch(selected_haplo_qtl_9k_08, haplo_evol20_qtl_9k_08, 100)                
score_evol40_qtl_9k_08 = compare_score_simBatch(selected_haplo_qtl_9k_08, haplo_evol40_qtl_9k_08, 100)                
score_evol80_qtl_9k_08 = compare_score_simBatch(selected_haplo_qtl_9k_08, haplo_evol80_qtl_9k_08, 100) 
score_evol140_qtl_9k_08 = compare_score_simBatch(selected_haplo_qtl_9k_08, haplo_evol140_qtl_9k_08, 100) 

mean_score_base_qtl_9k_08, std_score_base_qtl_9k_08_up, std_score_base_qtl_9k_08_down = mean_binned_HaploScore(score_base_qtl_9k_08,bins=np.arange(0,101,1))
mean_score_evol20_qtl_9k_08, std_score_evol20_qtl_9k_08_up, std_score_evol20_qtl_9k_08_down = mean_binned_HaploScore(score_evol20_qtl_9k_08,bins=np.arange(0,101,1))
mean_score_evol40_qtl_9k_08, std_score_evol40_qtl_9k_08_up, std_score_evol40_qtl_9k_08_down = mean_binned_HaploScore(score_evol40_qtl_9k_08,bins=np.arange(0,101,1))
mean_score_evol80_qtl_9k_08, std_score_evol80_qtl_9k_08_up, std_score_evol80_qtl_9k_08_down = mean_binned_HaploScore(score_evol80_qtl_9k_08,bins=np.arange(0,101,1))
mean_score_evol140_qtl_9k_08, std_score_evol140_qtl_9k_08_up, std_score_evol140_qtl_9k_08_down = mean_binned_HaploScore(score_evol140_qtl_9k_08,bins=np.arange(0,101,1))

#QTL 9k, 0.2 effect size
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--16.9-0.3_eff0.2_100loci_movingTO/effect_size.txt','r')   
selected_haplo_qtl_9k_2 = get_selected_haplo_qtl(InputFile)

#get the haplotypes for F0, F40, F70, F110, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--16.9-0.3_eff0.2_100loci_movingTO/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=9000

gen='g0'
haplo_base_qtl_9k_2 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_qtl_9k_2 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_qtl_9k_2 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_qtl_9k_2 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_qtl_9k_2 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_qtl_9k_2 = compare_score_simBatch(selected_haplo_qtl_9k_2, haplo_base_qtl_9k_2, 100)                
score_evol20_qtl_9k_2 = compare_score_simBatch(selected_haplo_qtl_9k_2, haplo_evol20_qtl_9k_2, 100)                
score_evol40_qtl_9k_2 = compare_score_simBatch(selected_haplo_qtl_9k_2, haplo_evol40_qtl_9k_2, 100)                
score_evol80_qtl_9k_2 = compare_score_simBatch(selected_haplo_qtl_9k_2, haplo_evol80_qtl_9k_2, 100) 
score_evol140_qtl_9k_2 = compare_score_simBatch(selected_haplo_qtl_9k_2, haplo_evol140_qtl_9k_2, 100) 

mean_score_base_qtl_9k_2, std_score_base_qtl_9k_2_up, std_score_base_qtl_9k_2_down = mean_binned_HaploScore(score_base_qtl_9k_2,bins=np.arange(0,101,1))
mean_score_evol20_qtl_9k_2, std_score_evol20_qtl_9k_2_up, std_score_evol20_qtl_9k_2_down = mean_binned_HaploScore(score_evol20_qtl_9k_2,bins=np.arange(0,101,1))
mean_score_evol40_qtl_9k_2, std_score_evol40_qtl_9k_2_up, std_score_evol40_qtl_9k_2_down = mean_binned_HaploScore(score_evol40_qtl_9k_2,bins=np.arange(0,101,1))
mean_score_evol80_qtl_9k_2, std_score_evol80_qtl_9k_2_up, std_score_evol80_qtl_9k_2_down = mean_binned_HaploScore(score_evol80_qtl_9k_2,bins=np.arange(0,101,1))
mean_score_evol140_qtl_9k_2, std_score_evol140_qtl_9k_2_up, std_score_evol140_qtl_9k_2_down = mean_binned_HaploScore(score_evol140_qtl_9k_2,bins=np.arange(0,101,1))

#QTL 9k, 0.4 effect size
#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--34.9-0.3_eff0.4_100loci_movingTO/effect_size.txt','r')   
selected_haplo_qtl_9k_4 = get_selected_haplo_qtl(InputFile)

#get the haplotypes for F0, F40, F70, F110, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--34.9-0.3_eff0.4_100loci_movingTO/2ndRun/evolved_haplos/*"

num_sims=50
num_haplos=9000

gen='g0'
haplo_base_qtl_9k_4 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g20'
haplo_evol20_qtl_9k_4 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_qtl_9k_4 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g80'
haplo_evol80_qtl_9k_4 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_qtl_9k_4 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_qtl_9k_4 = compare_score_simBatch(selected_haplo_qtl_9k_4, haplo_base_qtl_9k_4, 100)                
score_evol20_qtl_9k_4 = compare_score_simBatch(selected_haplo_qtl_9k_4, haplo_evol20_qtl_9k_4, 100)                
score_evol40_qtl_9k_4 = compare_score_simBatch(selected_haplo_qtl_9k_4, haplo_evol40_qtl_9k_4, 100)                
score_evol80_qtl_9k_4 = compare_score_simBatch(selected_haplo_qtl_9k_4, haplo_evol80_qtl_9k_4, 100) 
score_evol140_qtl_9k_4 = compare_score_simBatch(selected_haplo_qtl_9k_4, haplo_evol140_qtl_9k_4, 100) 

mean_score_base_qtl_9k_4, std_score_base_qtl_9k_4_up, std_score_base_qtl_9k_4_down = mean_binned_HaploScore(score_base_qtl_9k_4,bins=np.arange(0,101,1))
mean_score_evol20_qtl_9k_4, std_score_evol20_qtl_9k_4_up, std_score_evol20_qtl_9k_4_down = mean_binned_HaploScore(score_evol20_qtl_9k_4,bins=np.arange(0,101,1))
mean_score_evol40_qtl_9k_4, std_score_evol40_qtl_9k_4_up, std_score_evol40_qtl_9k_4_down = mean_binned_HaploScore(score_evol40_qtl_9k_4,bins=np.arange(0,101,1))
mean_score_evol80_qtl_9k_4, std_score_evol80_qtl_9k_4_up, std_score_evol80_qtl_9k_4_down = mean_binned_HaploScore(score_evol80_qtl_9k_4,bins=np.arange(0,101,1))
mean_score_evol140_qtl_9k_4, std_score_evol140_qtl_9k_4_up, std_score_evol140_qtl_9k_4_down = mean_binned_HaploScore(score_evol140_qtl_9k_4,bins=np.arange(0,101,1))

#######
#plot 
#######

base = [mean_score_base_w_450_02,mean_score_base_w_9k_02,mean_score_base_qtl_450_04,mean_score_base_qtl_9k_04,
        mean_score_base_w_450_05,mean_score_base_w_9k_05,mean_score_base_qtl_450_08,mean_score_base_qtl_9k_08,
        mean_score_base_w_450_08,mean_score_base_w_9k_08,mean_score_base_qtl_450_2,mean_score_base_qtl_9k_2,
        mean_score_base_w_450_1,mean_score_base_w_9k_1,mean_score_base_qtl_450_4,mean_score_base_qtl_9k_4]
F20 = [mean_score_evol20_w_450_02,mean_score_evol20_w_9k_02,mean_score_evol20_qtl_450_04,mean_score_evol20_qtl_9k_04,
        mean_score_evol20_w_450_05,mean_score_evol20_w_9k_05,mean_score_evol20_qtl_450_08,mean_score_evol20_qtl_9k_08,
        mean_score_evol20_w_450_08,mean_score_evol20_w_9k_08,mean_score_evol20_qtl_450_2,mean_score_evol20_qtl_9k_2,
        mean_score_evol20_w_450_1,mean_score_evol20_w_9k_1,mean_score_evol20_qtl_450_4,mean_score_evol20_qtl_9k_4]
F40 = [mean_score_evol40_w_450_02,mean_score_evol40_w_9k_02,mean_score_evol40_qtl_450_04,mean_score_evol40_qtl_9k_04,
        mean_score_evol40_w_450_05,mean_score_evol40_w_9k_05,mean_score_evol40_qtl_450_08,mean_score_evol40_qtl_9k_08,
        mean_score_evol40_w_450_08,mean_score_evol40_w_9k_08,mean_score_evol40_qtl_450_2,mean_score_evol40_qtl_9k_2,
        mean_score_evol40_w_450_1,mean_score_evol40_w_9k_1,mean_score_evol40_qtl_450_4,mean_score_evol40_qtl_9k_4]
F80 = [mean_score_evol80_w_450_02,mean_score_evol80_w_9k_02,mean_score_evol80_qtl_450_04,mean_score_evol80_qtl_9k_04,
        mean_score_evol80_w_450_05,mean_score_evol80_w_9k_05,mean_score_evol80_qtl_450_08,mean_score_evol80_qtl_9k_08,
        mean_score_evol80_w_450_08,mean_score_evol80_w_9k_08,mean_score_evol80_qtl_450_2,mean_score_evol80_qtl_9k_2,
        mean_score_evol80_w_450_1,mean_score_evol80_w_9k_1,mean_score_evol80_qtl_450_4,mean_score_evol80_qtl_9k_4]
F140 = [mean_score_evol140_w_450_02,mean_score_evol140_w_9k_02,mean_score_evol140_qtl_450_04,mean_score_evol140_qtl_9k_04,
        mean_score_evol140_w_450_05,mean_score_evol140_w_9k_05,mean_score_evol140_qtl_450_08,mean_score_evol140_qtl_9k_08,
        mean_score_evol140_w_450_08,mean_score_evol140_w_9k_08,mean_score_evol140_qtl_450_2,mean_score_evol140_qtl_9k_2,
        mean_score_evol140_w_450_1,mean_score_evol140_w_9k_1,mean_score_evol140_qtl_450_4,mean_score_evol140_qtl_9k_4]

base_up = [std_score_base_w_450_02_up,std_score_base_w_9k_02_up,std_score_base_qtl_450_04_up,std_score_base_qtl_9k_04_up,
        std_score_base_w_450_05_up,std_score_base_w_9k_05_up,std_score_base_qtl_450_08_up,std_score_base_qtl_9k_08_up,
        std_score_base_w_450_08_up,std_score_base_w_9k_08_up,std_score_base_qtl_450_2_up,std_score_base_qtl_9k_2_up,
        std_score_base_w_450_1_up,std_score_base_w_9k_1_up,std_score_base_qtl_450_4_up,std_score_base_qtl_9k_4_up]
F20_up = [std_score_evol20_w_450_02_up,std_score_evol20_w_9k_02_up,std_score_evol20_qtl_450_04_up,std_score_evol20_qtl_9k_04_up,
        std_score_evol20_w_450_05_up,std_score_evol20_w_9k_05_up,std_score_evol20_qtl_450_08_up,std_score_evol20_qtl_9k_08_up,
        std_score_evol20_w_450_08_up,std_score_evol20_w_9k_08_up,std_score_evol20_qtl_450_2_up,std_score_evol20_qtl_9k_2_up,
        std_score_evol20_w_450_1_up,std_score_evol20_w_9k_1_up,std_score_evol20_qtl_450_4_up,std_score_evol20_qtl_9k_4_up]
F40_up = [std_score_evol40_w_450_02_up,std_score_evol40_w_9k_02_up,std_score_evol40_qtl_450_04_up,std_score_evol40_qtl_9k_04_up,
        std_score_evol40_w_450_05_up,std_score_evol40_w_9k_05_up,std_score_evol40_qtl_450_08_up,std_score_evol40_qtl_9k_08_up,
        std_score_evol40_w_450_08_up,std_score_evol40_w_9k_08_up,std_score_evol40_qtl_450_2_up,std_score_evol40_qtl_9k_2_up,
        std_score_evol40_w_450_1_up,std_score_evol40_w_9k_1_up,std_score_evol40_qtl_450_4_up,std_score_evol40_qtl_9k_4_up]
F80_up = [std_score_evol80_w_450_02_up,std_score_evol80_w_9k_02_up,std_score_evol80_qtl_450_04_up,std_score_evol80_qtl_9k_04_up,
        std_score_evol80_w_450_05_up,std_score_evol80_w_9k_05_up,std_score_evol80_qtl_450_08_up,std_score_evol80_qtl_9k_08_up,
        std_score_evol80_w_450_08_up,std_score_evol80_w_9k_08_up,std_score_evol80_qtl_450_2_up,std_score_evol80_qtl_9k_2_up,
        std_score_evol80_w_450_1_up,std_score_evol80_w_9k_1_up,std_score_evol80_qtl_450_4_up,std_score_evol80_qtl_9k_4_up]
F140_up = [std_score_evol140_w_450_02_up,std_score_evol140_w_9k_02_up,std_score_evol140_qtl_450_04_up,std_score_evol140_qtl_9k_04_up,
        std_score_evol140_w_450_05_up,std_score_evol140_w_9k_05_up,std_score_evol140_qtl_450_08_up,std_score_evol140_qtl_9k_08_up,
        std_score_evol140_w_450_08_up,std_score_evol140_w_9k_08_up,std_score_evol140_qtl_450_2_up,std_score_evol140_qtl_9k_2_up,
        std_score_evol140_w_450_1_up,std_score_evol140_w_9k_1_up,std_score_evol140_qtl_450_4_up,std_score_evol140_qtl_9k_4_up]

base_down = [std_score_base_w_450_02_down,std_score_base_w_9k_02_down,std_score_base_qtl_450_04_down,std_score_base_qtl_9k_04_down,
        std_score_base_w_450_05_down,std_score_base_w_9k_05_down,std_score_base_qtl_450_08_down,std_score_base_qtl_9k_08_down,
        std_score_base_w_450_08_down,std_score_base_w_9k_08_down,std_score_base_qtl_450_2_down,std_score_base_qtl_9k_2_down,
        std_score_base_w_450_1_down,std_score_base_w_9k_1_down,std_score_base_qtl_450_4_down,std_score_base_qtl_9k_4_down]
F20_down = [std_score_evol20_w_450_02_down,std_score_evol20_w_9k_02_down,std_score_evol20_qtl_450_04_down,std_score_evol20_qtl_9k_04_down,
        std_score_evol20_w_450_05_down,std_score_evol20_w_9k_05_down,std_score_evol20_qtl_450_08_down,std_score_evol20_qtl_9k_08_down,
        std_score_evol20_w_450_08_down,std_score_evol20_w_9k_08_down,std_score_evol20_qtl_450_2_down,std_score_evol20_qtl_9k_2_down,
        std_score_evol20_w_450_1_down,std_score_evol20_w_9k_1_down,std_score_evol20_qtl_450_4_down,std_score_evol20_qtl_9k_4_down]
F40_down = [std_score_evol40_w_450_02_down,std_score_evol40_w_9k_02_down,std_score_evol40_qtl_450_04_down,std_score_evol40_qtl_9k_04_down,
        std_score_evol40_w_450_05_down,std_score_evol40_w_9k_05_down,std_score_evol40_qtl_450_08_down,std_score_evol40_qtl_9k_08_down,
        std_score_evol40_w_450_08_down,std_score_evol40_w_9k_08_down,std_score_evol40_qtl_450_2_down,std_score_evol40_qtl_9k_2_down,
        std_score_evol40_w_450_1_down,std_score_evol40_w_9k_1_down,std_score_evol40_qtl_450_4_down,std_score_evol40_qtl_9k_4_down]
F80_down = [std_score_evol80_w_450_02_down,std_score_evol80_w_9k_02_down,std_score_evol80_qtl_450_04_down,std_score_evol80_qtl_9k_04_down,
        std_score_evol80_w_450_05_down,std_score_evol80_w_9k_05_down,std_score_evol80_qtl_450_08_down,std_score_evol80_qtl_9k_08_down,
        std_score_evol80_w_450_08_down,std_score_evol80_w_9k_08_down,std_score_evol80_qtl_450_2_down,std_score_evol80_qtl_9k_2_down,
        std_score_evol80_w_450_1_down,std_score_evol80_w_9k_1_down,std_score_evol80_qtl_450_4_down,std_score_evol80_qtl_9k_4_down]
F140_down = [std_score_evol140_w_450_02_down,std_score_evol140_w_9k_02_down,std_score_evol140_qtl_450_04_down,std_score_evol140_qtl_9k_04_down,
        std_score_evol140_w_450_05_down,std_score_evol140_w_9k_05_down,std_score_evol140_qtl_450_08_down,std_score_evol140_qtl_9k_08_down,
        std_score_evol140_w_450_08_down,std_score_evol140_w_9k_08_down,std_score_evol140_qtl_450_2_down,std_score_evol140_qtl_9k_2_down,
        std_score_evol140_w_450_1_down,std_score_evol140_w_9k_1_down,std_score_evol140_qtl_450_4_down,std_score_evol140_qtl_9k_4_down]
       
model = ['sweep','sweep','trait optimum','trait optimum']
pop_size = [450,9000,450,9000]
color_label = ['#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494']
x = np.arange(0,101,1)

fig , ax = plt.subplots(nrows=4, ncols=4, sharex=True, sharey=True,figsize=(26,26),dpi=100, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.2, wspace=0.2)
for rep in range(16):
    ax=plt.subplot(4,4,rep+1)
    plt.plot(x,base[rep],color = color_label[0], linewidth = 3)
    plt.fill_between(x, base_up[rep], base_down[rep], color = color_label[0], alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(x,F20[rep],color = color_label[1], linewidth = 3)
    plt.fill_between(x, F20_up[rep], F20_down[rep], color = color_label[1], alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(x,F40[rep],color = color_label[2], linewidth = 3)
    plt.fill_between(x, F40_up[rep], F40_down[rep], color = color_label[2], alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(x,F80[rep],color = color_label[3], linewidth = 3)
    plt.fill_between(x, F80_up[rep], F80_down[rep], color = color_label[3], alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(x,F140[rep],color = color_label[4], linewidth = 3)
    plt.fill_between(x, F140_up[rep], F140_down[rep], color = color_label[4], alpha = 0.3, edgecolor = None, linewidth=0.0)
    ax.tick_params(axis='x', which='both' ,bottom='off', top='off')
    ax.tick_params(axis='y', which='both' ,right='off', left='on') 
    ylim = [0.25]*16
    plt.ylim(0,ylim[rep])
    xlim = [40]*4+[100,100,40,40]*3
    plt.xlim(0,xlim[rep])
    plt.xticks(fontsize=12, va='top')
    plt.yticks(fontsize=12, va='top') 
    if rep == 0:
        ax.set_title("%s, N = %i" %(model[rep],pop_size[rep]),fontsize=24)
        plt.ylabel('s=0.02/eff_size=0.04',fontsize=24, labelpad = 20)#rotation=0
    if rep == 1:ax.set_title("%s, N = %i" %(model[rep],pop_size[rep]),fontsize=24)
    if rep == 2:ax.set_title("%s, N = %i" %(model[rep],pop_size[rep]),fontsize=24)
    if rep == 3:ax.set_title("%s, N = %i" %(model[rep],pop_size[rep]),fontsize=24)
    if rep == 4:plt.ylabel('s=0.05/eff_size=0.08',fontsize=24, labelpad = 20)
    if rep == 8:plt.ylabel('s=0.08/eff_size=0.2',fontsize=24, labelpad = 20)
    if rep == 12:plt.ylabel('s=0.1/eff_size=0.4',fontsize=24, labelpad = 20)
    if rep == 3:
        texts = ['Gen 0','Gen 20','Gen 40','Gen 80','Gen 140']
        patches = [mpatches.Patch(color=color_label[i], label="{:s}".format(texts[i]) ) for i in range(len(texts))]
        ax.legend(handles=patches,loc='upper left',bbox_to_anchor=(0.67, 1),frameon=False)    
fig.text(0.05, 0.5, 'Proportion', va='center', rotation='vertical',size = 24) #ylabel 
fig.text(0.5, 0.07, 'No. of selected loci per haplotype', ha='center',size = 24) #xlabel
plt.savefig('sweep_qtl_diffEffSize_selectedLociPerHaplo_proportion.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_qtl_diffEffSize_selectedLociPerHaplo_proportiony.pdf', dpi=300,format='pdf', bbox_inches = 'tight')
       
  