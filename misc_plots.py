# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 10:29:19 2020

@author: Neda
"""
import gzip
import glob
import operator
import matplotlib.patches as mpatches
import numpy as np
import matplotlib.pyplot as plt

#This function requires the sel.txt (for sweep simulation) and the output file from w function of mimicree2
#it will compute and store the allele frequency of seleted allele
def get_AF_w(infile_s, infile_af): 
    chr_pos_s =[]
    linenum = 0
    for line in infile_s:
        linenum += 1
        line = line.rstrip()
        if linenum == 1:
            header = line
        else:
            cols = line.split('\t')        
            chr_pos_s.append(cols[:])
    infile_s.close()  

    base_index = {'A':0, 'T':1, 'C':2, 'G':3}
    AF_sim = []
    for line in infile_af:
        line = line.rstrip()
        cols = line.split('\t')
        for thing in chr_pos_s:
            if cols[0] == thing[0]:
                if cols[1] == thing[1]:
                    unselect, select = base_index[thing[2].split('/')[0]],base_index[thing[2].split('/')[1]]
                    snp_AF = []
                    for ind,rep in enumerate(cols[3:]):
                        AF = int(rep.split(':')[select]) / float(int(rep.split(':')[select])+int(rep.split(':')[unselect]))
                        snp_AF.append(AF)
                    AF_sim.append(snp_AF)
    infile_af.close()
    return AF_sim

#This function requires the effect_size.txt (for quantitative and the output file from qff function of mimicree2
#it will compute and store the allele frequency of seleted allele
def get_AF_qff(infile_s, infile_af): 
    chr_pos_s =[]
    for line in infile_s:
        line = line.rstrip()
        cols = line.split('\t')        
        chr_pos_s.append(cols[:])
    infile_s.close()  

    base_index = {'A':0, 'T':1, 'C':2, 'G':3}
    AF_sim = []
    for line in infile_af:
        line = line.rstrip()
        cols = line.split('\t')
        for thing in chr_pos_s:
            if cols[0] == thing[0]:
                if cols[1] == thing[1]:
                    unselect, select = base_index[thing[2].split('/')[1]],base_index[thing[2].split('/')[0]]
                    snp_AF = []
                    for ind,rep in enumerate(cols[3:]):
                        AF = int(rep.split(':')[select]) / float((int(rep.split(':')[select])+int(rep.split(':')[unselect])))
                        snp_AF.append(AF)
                    AF_sim.append(snp_AF)
    infile_af.close()
    return AF_sim

def process_AF_sim(AF_sim,loci_num,replicates):
    #classify the frequency of different loci and replicates based on generations/timepoints
    #replicates = 50
    timepoints = 15   
    timepoint_index = []
    for i in range(0,timepoints):
        timepoint_index.append(np.arange(i,(replicates*timepoints),timepoints))
    
    AF_sim_classified = [list() for _ in xrange(timepoints)]
    for item in AF_sim: #each item is one locus and for each locus we have replicate*timepoint entries
        for ind,timepoint in enumerate(timepoint_index):
            for t in timepoint:
                AF_sim_classified[ind].append(item[t])
    
    #to compute the AFC for every timepoint compared to F0
    #first divide 5000 data points for each timepoint (100loci * 50 sims) into separate sims
    #sim_num = 50
    #loci_num = 100   
    sim_index = []
    for i in np.arange(0,loci_num*replicates,loci_num):
        sim_index.append([i,i+loci_num])
    
    #get AFC for every simulation between specific timepoint compared to F0
    AFC_sims_alltimepoints=[]
    for item in AF_sim_classified[1:]:
        alltimepoints = [map(operator.sub, item[s[0]:s[1]], AF_sim_classified[0][s[0]:s[1]]) for s in sim_index]
        AFC_sims_alltimepoints.append(alltimepoints)
    
    #separate AFC for every timepoint
    separate_afc = []
    for tp in AFC_sims_alltimepoints:
        tmp_af = []
        for sim in tp:
            for af in sim:
                tmp_af.append(af)
        separate_afc.append(tmp_af)
    
    #get 5% upper tail for each timepoint and set them as threshold
    AFC_thresh = []
    for tp in separate_afc:
        AFC_thresh.append(round(np.percentile(tp,95),3)) 
    
    return(AF_sim_classified,AFC_sims_alltimepoints,AFC_thresh)

#store the allele freqeuncies in specific simulations that were sun until generation 5000
def process_AF_sim_F5000(AF_sim,loci_num):
    #classify the frequency of different loci and replicates based on generations/timepoints
    replicates = 1
    timepoints = 501   
    timepoint_index = []
    for i in range(0,timepoints):
        timepoint_index.append(np.arange(i,(replicates*timepoints),timepoints))
    
    AF_sim_classified = [list() for _ in xrange(timepoints)]
    for item in AF_sim: #each item is one locus and for each locus we have replicate*timepoint entries
        for ind,timepoint in enumerate(timepoint_index):
            for t in timepoint:
                AF_sim_classified[ind].append(item[t])
    
    #to compute the AFC for every timepoint compared to F0
    #first divide 5000 data points for each timepoint (100loci * 50 sims) into separate sims
    sim_num = 1
    #loci_num = 100   
    sim_index = []
    for i in np.arange(0,loci_num*sim_num,loci_num):
        sim_index.append([i,i+loci_num])
    
    #get AFC for every simulation between specific timepoint compared to F0
    AFC_sims_alltimepoints=[]
    for item in AF_sim_classified[1:]:
        alltimepoints = [map(operator.sub, item[s[0]:s[1]], AF_sim_classified[0][s[0]:s[1]]) for s in sim_index]
        AFC_sims_alltimepoints.append(alltimepoints)
    
    #separate AFC for every timepoint
    separate_afc = []
    for tp in AFC_sims_alltimepoints:
        tmp_af = []
        for sim in tp:
            for af in sim:
                tmp_af.append(af)
        separate_afc.append(tmp_af)
    
    #get 5% upper tail for each timepoint and set them as threshold
    AFC_thresh = []
    for tp in separate_afc:
        AFC_thresh.append(round(np.percentile(tp,95),3)) 
    
    return(AF_sim_classified,AFC_sims_alltimepoints,AFC_thresh)

#get the fitness of sweep simulations and categorize them based on timepoint and replicate
#then compute variance from the median fitness of all replicates         
def get_fitnessVar_SimBatch_w(InputFile_f,sim_nums):
    gen_index = {0:0,10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
    fitness = [[[] for _ in xrange(sim_nums)] for _ in xrange(len(gen_index))]
    for line in InputFile_f:
        line = line.rstrip()
        cols = line.split('\t')
        if int(cols[0]) <= sim_nums:
            gen,rep=gen_index[int(cols[1])],int(cols[0])
            fitness[gen][rep-1].append(np.log10(float(cols[5])))
    #compute variance of fitness from median of replicates for each timepoint
    mean_of_varaince_fit = []
    for tmp in fitness:
        mean_of_varaince_fit.append(np.mean([np.var(r) for r in tmp]))
    return mean_of_varaince_fit
    InputFile_f.close()

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

#extract haplotype information from the stored haplotypes from simulations
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

#count the number of the selected loci in each haplotype
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

def set_color_flier(col):
    flierprops = dict( marker = '.',color = col, markerfacecolor = col, markeredgecolor = col, markersize =3)   
    return flierprops

def set_color_box(col):    
    colors=[ col] * 100
    fullset_color_box = colors
    fullset_color_whisker = [ col] * 80
    for patch, col in zip(box['boxes'], fullset_color_box):
        patch.set_facecolor(col)
        patch.set_edgecolor(col)
    for patch, col in zip(box['whiskers'], fullset_color_whisker):
        patch.set_color(col)
        patch.set_linestyle('-')

'''
###########
#Figure S1
###########
'''

#450
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.0s/sel.txt','r')   
InputFile_af= gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.0s/recom_rate_noX/500sims/AF_results.sync','rb')

AF_sim_neu_450 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_neu_450, AFC_sim_alltmpnt_neu_450, AFC_thresh_neu_450 = process_AF_sim(AF_sim_neu_450,100,500)

#9000
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.0s/sel.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.0s/recom_rate_noX/500sims/AF_results.sync','rb')

AF_sim_neu_9k = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_neu_9k, AFC_sim_alltmpnt_neu_9k, AFC_thresh_neu_9k = process_AF_sim(AF_sim_neu_9k,100,500)

######
#plot 
######

data_450 = [AF_sim_class_neu_450]
data_9k = [AF_sim_class_neu_9k]

gen_index = {0:0,10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
color = ['darkkhaki','#585528']
color_label = ['darkkhaki','#585528']

fig , ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True,figsize=(9,5),dpi=100, facecolor='w', edgecolor='k')
n = 0
for rep in range(1):
    ax=plt.subplot(1,1,rep+1)
    set_color_flier(color[n])
    box=ax.boxplot(data_450[rep],positions = np.arange(1,(len(gen_index)*2)+len(gen_index),3), vert = True,showcaps = False, showfliers = True, patch_artist=True, flierprops=set_color_flier(color[n]))
    set_color_box(color[n])
    n+=1
    for element in ['medians']:plt.setp(box[element], color='black')    
    set_color_flier(color[n])
    plt.plot(np.arange(4,(len(gen_index)*2)+len(gen_index),3),AFC_thresh_neu_450,'*', color = 'red', markeredgecolor = 'red')
    box=ax.boxplot(data_9k[rep],positions = np.arange(2,(len(gen_index)*2)+len(gen_index),3), vert = True,showcaps = False, showfliers = True, patch_artist=True, flierprops=set_color_flier(color[n]))
    set_color_box(color[n])
    n+=1
    for element in ['medians']:plt.setp(box[element], color='black')
    plt.plot(np.arange(5,(len(gen_index)*2)+len(gen_index),3),AFC_thresh_neu_9k,'*',color = 'red', markeredgecolor = 'red')
    ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
    ax.tick_params(axis='y', which='both' ,right='off', left='on')
    plt.xlim(0,(len(gen_index)*2)+len(gen_index))   
    plt.ylim(-0.02, 1.02)
    plt.xticks(np.arange(1,(len(gen_index)*2)+len(gen_index),3), [i for i in sorted(gen_index.keys())],fontsize=12, va='top')
    plt.yticks(fontsize=12, va='top')
    plt.xlabel('Generation', size = 14,labelpad=15)
    plt.ylabel('Allele frequency', size = 14,labelpad=15)
    texts = ['450','9000']
    patches = [mpatches.Patch(color=color_label[i], label="{:s}".format(texts[i]) ) for i in range(len(texts))]
    ax.legend(handles=patches,bbox_to_anchor=(0.69, 1), loc=3,ncol=2, borderaxespad=0.)    
plt.savefig('AF_drift_defaultParameters_withMarker.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('AF_drift_defaultParameters_withMarker.pdf', dpi=300,format='pdf', bbox_inches = 'tight')


'''
##########
#Figure S3
##########
'''

#QTL 450
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/effect_size.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/recom_rate_noX/gen5000/AF.sync.gz','rb')  

AF_sim_450_qtl = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_450_qtl, AFC_sim_alltmpnt_450_qtl, AFC_thresh_450_qtl = process_AF_sim_F5000(AF_sim_450_qtl,100)

#using this indexing we'll get all 100 loci in the first simulation
AF_sim_class_450_qtl_sub = [list() for _ in range(100)]
for tp in AF_sim_class_450_qtl:
    for ind,loci in enumerate(tp):
        AF_sim_class_450_qtl_sub[ind].append(loci)

#plot the results only until generation 2500 because all alleles are fixed or lost at this point      
gen_index = {}
n=0
for item in np.arange(0,2510,10):
    gen_index[item] = n
    n+=1
    
#log-transform the x axis
fig , ax = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True,figsize=(18,4.5),dpi=100, facecolor='w', edgecolor='k')
for ind,item in enumerate(AF_sim_class_450_qtl_sub[:20]):
    plt.plot(np.log10(range(1,len(item[:250])+1)),item[:250], '-', color = '#e699b2', markerfacecolor = '#e699b2', markeredgecolor = '#e699b2', alpha = 1)
ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
ax.tick_params(axis='y', which='both' ,right='off', left='on')
plt.xlim(-0.02,np.log10(len(gen_index)+10))   
sel_xtick_label = [index for index,item in enumerate(np.arange(0,2510,100)) if item in [0,100,200,300,400,500,1000,2000]]
sel_xtick = [x for index,x in enumerate(np.log10(np.arange(1,len(gen_index)+1,10))) if index in sel_xtick_label]
plt.xticks(sel_xtick, [0,100,200,300,400,500,1000,2000],fontsize=12, va='top')
plt.yticks(fontsize=12, va='top')
plt.ylim(-0.02, 1.02)
plt.ylabel('Allele frequency', size = 14,labelpad=15)
plt.xlabel('Generation', size = 14,labelpad=15)
plt.tight_layout()
plt.savefig('AF_1simSample_Ne450_5kgen_2500genShown_20loci_log10.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('AF_1simSample_Ne450_5kgen_2500genShown_20loci_log10.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''
###########
#Figure S4
###########
'''

#####################
#get simulation data
#####################

#sweep 450
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/sel.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/recom_rate_noX/1000sims/AF_results.sync','rb')  

AF_sim_450_w = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_450_w, AFC_sim_alltmpnt_450_w, AFC_thresh_450_w = process_AF_sim(AF_sim_450_w,100,500)

#sweep 9k
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.08s/sel.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.08s/recom_rate_noX/500sims/AF_results.sync','rb') 

AF_sim_9k_w = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_9k_w, AFC_sim_alltmpnt_9k_w, AFC_thresh_9k_w = process_AF_sim(AF_sim_9k_w,100,500)

#qtl 450
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/effect_size.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/recom_rate_noX/1000sims/AF.sync.gz','rb')   

AF_sim_450_qtl = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_450_qtl, AFC_sim_alltmpnt_450_qtl, AFC_thresh_450_qtl = process_AF_sim(AF_sim_450_qtl,100,500)

#qtl 9k
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--2.5-0.3/effect_size.txt','r')  
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--2.5-0.3/recom_rate_noX/500sims/AF.sync.gz','rb')  

AF_sim_9k_qtl = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_9k_qtl, AFC_sim_alltmpnt_9k_qtl, AFC_thresh_9k_qtl = process_AF_sim(AF_sim_9k_qtl,100,500)

######
#plot 
######

#using this indexing we'll get all 100 loci in the first simulation
AF_sim_class_450_sub = [[tp[i] for tp in AF_sim_class_450_w] for i in np.arange(0,50000,500)]
AF_sim_class_9k_sub = [[tp[i] for tp in AF_sim_class_9k_w] for i in np.arange(0,50000,500)]
AF_sim_class_450_qtl_sub = [[tp[i] for tp in AF_sim_class_450_qtl] for i in np.arange(0,50000,500)]
AF_sim_class_9k_qtl_sub = [[tp[i] for tp in AF_sim_class_9k_qtl] for i in np.arange(0,50000,500)]

data_toplot = [AF_sim_class_450_sub,AF_sim_class_9k_sub,AF_sim_class_450_qtl_sub,AF_sim_class_9k_qtl_sub]

gen_index = {0:0,10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
color = ['turquoise','teal','#e699b2','#cc3366']
model = ['sweep','sweep','trait optimum','trait optimum']
pop_size = [450,9000,450,9000]

fig , ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True,figsize=(18,10),dpi=100, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.3, wspace=0.1)
for rep in range(4):
    ax=plt.subplot(2,2,rep+1)
    for ind,item in enumerate(data_toplot[rep]):
        plt.plot(range(1,len(item)+1),item, '-', color = color[rep], markerfacecolor = color[rep], markeredgecolor = color[rep], alpha = 0.5)
    ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
    ax.tick_params(axis='y', which='both' ,right='off', left='on')
    ax.set_title("%s, N = %i" %(model[rep],pop_size[rep]),fontsize=12)
    plt.xlim(0.8,len(gen_index)+0.1)   
    plt.xticks(np.arange(1,len(gen_index)+1,1), [i for i in sorted(gen_index.keys())],fontsize=12, va='top')
    plt.yticks(fontsize=12, va='top')
    plt.ylim(-0.02, 1.02)
    if rep == 0 or rep == 2:ax.set_ylabel('Allele frequency', size = 14,labelpad=15)
    if rep == 2 or rep == 3:ax.set_xlabel('Generation', size = 14,labelpad=15)
plt.tight_layout()
plt.savefig('sweep_qtl_AF_1simSample_defaultParameters_outliers.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_qtl_AF_1simSample_defaultParameters_outliers.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''
###########
#Figure S8
##########
'''

#################################################################
#get fitness for sweep simulations with different number of loci
#################################################################

sim='Ne450_10loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
var_fit_450_10 = get_fitnessVar_SimBatch_w(InputFile_f,500)

sim='Ne450_20loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
var_fit_450_20 = get_fitnessVar_SimBatch_w(InputFile_f,500)

sim='Ne450_50loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/fitness_result.txt' %sim,'r')  
var_fit_450_50 = get_fitnessVar_SimBatch_w(InputFile_f,500)

sim='Ne450_100loci_0.05p0_0.08s'
InputFile_f = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/1000sims/fitness_result.txt' %sim,'r')  
var_fit_450_100 = get_fitnessVar_SimBatch_w(InputFile_f,500)

############################################
#get number of selected sites per haplotype 
#for each replicate and then compute varaince
#in each replicate
############################################

#get selected sites as a single haplotype
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_10loci_0.05p0_0.08s/sel.txt','r')   
selected_haplo_w_450_10 = get_selected_haplo_w(InputFile)

#get the haplotypes for F0, F40, F70, F110, F140
path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_10loci_0.05p0_0.08s/recom_rate_noX/2ndrun/evolved_haplos/*"

num_sims=50
num_haplos=450

gen='g0'
haplo_base_w_450_10 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_450_10 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g70'
haplo_evol70_w_450_10 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g110'
haplo_evol110_w_450_10 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_450_10 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_450_10 = compare_score_simBatch(selected_haplo_w_450_10, haplo_base_w_450_10, 10)                
score_evol40_w_450_10 = compare_score_simBatch(selected_haplo_w_450_10, haplo_evol40_w_450_10, 10)                
score_evol70_w_450_10 = compare_score_simBatch(selected_haplo_w_450_10, haplo_evol70_w_450_10, 10)                
score_evol110_w_450_10 = compare_score_simBatch(selected_haplo_w_450_10, haplo_evol110_w_450_10, 10) 
score_evol140_w_450_10 = compare_score_simBatch(selected_haplo_w_450_10, haplo_evol140_w_450_10, 10) 

var_score_base_w_450_10 = np.mean([np.var(rep) for rep in score_base_w_450_10])
var_score_evol40_w_450_10 = np.mean([np.var(rep) for rep in score_evol40_w_450_10])
var_score_evol70_w_450_10 = np.mean([np.var(rep) for rep in score_evol70_w_450_10])
var_score_evol110_w_450_10 = np.mean([np.var(rep) for rep in score_evol110_w_450_10])
var_score_evol140_w_450_10 = np.mean([np.var(rep) for rep in score_evol140_w_450_10])

#sweep 450 20 loci
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_20loci_0.05p0_0.08s/sel.txt','r')   
selected_haplo_w_450_20 = get_selected_haplo_w(InputFile)

path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_20loci_0.05p0_0.08s/recom_rate_noX/2ndRun/evolved_haplos/*"

gen='g0'
haplo_base_w_450_20 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_450_20 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g70'
haplo_evol70_w_450_20 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g110'
haplo_evol110_w_450_20 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_450_20 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_450_20 = compare_score_simBatch(selected_haplo_w_450_20, haplo_base_w_450_20, 20)                
score_evol40_w_450_20 = compare_score_simBatch(selected_haplo_w_450_20, haplo_evol40_w_450_20, 20)                
score_evol70_w_450_20 = compare_score_simBatch(selected_haplo_w_450_20, haplo_evol70_w_450_20, 20)                
score_evol110_w_450_20 = compare_score_simBatch(selected_haplo_w_450_20, haplo_evol110_w_450_20, 20) 
score_evol140_w_450_20 = compare_score_simBatch(selected_haplo_w_450_20, haplo_evol140_w_450_20, 20) 

var_score_base_w_450_20 = np.mean([np.var(rep) for rep in score_base_w_450_20])
var_score_evol40_w_450_20 = np.mean([np.var(rep) for rep in score_evol40_w_450_20])
var_score_evol70_w_450_20 = np.mean([np.var(rep) for rep in score_evol70_w_450_20])
var_score_evol110_w_450_20 = np.mean([np.var(rep) for rep in score_evol110_w_450_20])
var_score_evol140_w_450_20 = np.mean([np.var(rep) for rep in score_evol140_w_450_20])

#sweep 450 50 loci
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_50loci_0.05p0_0.08s/sel.txt','r')   
selected_haplo_w_450_50 = get_selected_haplo_w(InputFile)

path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_50loci_0.05p0_0.08s/recom_rate_noX/2ndRun/evolved_haplos/*"

gen='g0'
haplo_base_w_450_50 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_450_50 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g70'
haplo_evol70_w_450_50 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g110'
haplo_evol110_w_450_50 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_450_50 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_450_50 = compare_score_simBatch(selected_haplo_w_450_50, haplo_base_w_450_50, 50)                
score_evol40_w_450_50 = compare_score_simBatch(selected_haplo_w_450_50, haplo_evol40_w_450_50, 50)                
score_evol70_w_450_50 = compare_score_simBatch(selected_haplo_w_450_50, haplo_evol70_w_450_50, 50)                
score_evol110_w_450_50 = compare_score_simBatch(selected_haplo_w_450_50, haplo_evol110_w_450_50, 50) 
score_evol140_w_450_50 = compare_score_simBatch(selected_haplo_w_450_50, haplo_evol140_w_450_50, 50) 

var_score_base_w_450_50 = np.mean([np.var(rep) for rep in score_base_w_450_50])
var_score_evol40_w_450_50 = np.mean([np.var(rep) for rep in score_evol40_w_450_50])
var_score_evol70_w_450_50 = np.mean([np.var(rep) for rep in score_evol70_w_450_50])
var_score_evol110_w_450_50 = np.mean([np.var(rep) for rep in score_evol110_w_450_50])
var_score_evol140_w_450_50 = np.mean([np.var(rep) for rep in score_evol140_w_450_50])

#sweep 450 100 loci
InputFile = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/sel.txt','r')   
selected_haplo_w_450_100 = get_selected_haplo_w(InputFile)

path = "/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/recom_rate_noX/3rdRun/evolved_haplos/*"

gen='g0'
haplo_base_w_450_100 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g40'
haplo_evol40_w_450_100 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g70'
haplo_evol70_w_450_100 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g110'
haplo_evol110_w_450_100 = get_haplo_info(path, num_sims, num_haplos, gen)

gen='g140'
haplo_evol140_w_450_100 = get_haplo_info(path, num_sims, num_haplos, gen)

score_base_w_450_100 = compare_score_simBatch(selected_haplo_w_450_100, haplo_base_w_450_100, 100)                
score_evol40_w_450_100 = compare_score_simBatch(selected_haplo_w_450_100, haplo_evol40_w_450_100, 100)                
score_evol70_w_450_100 = compare_score_simBatch(selected_haplo_w_450_100, haplo_evol70_w_450_100, 100)                
score_evol110_w_450_100 = compare_score_simBatch(selected_haplo_w_450_100, haplo_evol110_w_450_100, 100) 
score_evol140_w_450_100 = compare_score_simBatch(selected_haplo_w_450_100, haplo_evol140_w_450_100, 100) 

var_score_base_w_450_100 = np.mean([np.var(rep) for rep in score_base_w_450_100])
var_score_evol40_w_450_100 = np.mean([np.var(rep) for rep in score_evol40_w_450_100])
var_score_evol70_w_450_100 = np.mean([np.var(rep) for rep in score_evol70_w_450_100])
var_score_evol110_w_450_100 = np.mean([np.var(rep) for rep in score_evol110_w_450_100])
var_score_evol140_w_450_100 = np.mean([np.var(rep) for rep in score_evol140_w_450_100])

######################
#plot the variance of
#fitness and no. of 
#selected loci per haplotype
######################

haplo_w_450_10 = [var_score_base_w_450_10,var_score_evol40_w_450_10,var_score_evol70_w_450_10,var_score_evol110_w_450_10,var_score_evol140_w_450_10]
haplo_w_450_20 = [var_score_base_w_450_20,var_score_evol40_w_450_20,var_score_evol70_w_450_20,var_score_evol110_w_450_20,var_score_evol140_w_450_20]
haplo_w_450_50 = [var_score_base_w_450_50,var_score_evol40_w_450_50,var_score_evol70_w_450_50,var_score_evol110_w_450_50,var_score_evol140_w_450_50]
haplo_w_450_100 = [var_score_base_w_450_100,var_score_evol40_w_450_100,var_score_evol70_w_450_100,var_score_evol110_w_450_100,var_score_evol140_w_450_100]

linestyle = [':','-.','--','-']*14
labels = ['450, 10 loci', '450, 20 loci','450, 50 loci','450, 100 loci']

fig , ax = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True,figsize=(8,10),dpi=100, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.4, wspace=0.1)

#plot variance in number of selected alleles on haplotypes
ax=plt.subplot(2,1,1)
plt.plot([0,4,7,11,14],haplo_w_450_10, linestyle = linestyle[0] , color = 'turquoise', linewidth = 3, label = labels[0])
plt.plot([0,4,7,11,14],haplo_w_450_20, linestyle = linestyle[1] , color = 'turquoise', linewidth = 3, label = labels[1])
plt.plot([0,4,7,11,14],haplo_w_450_50, linestyle = linestyle[2] , color = 'turquoise', linewidth = 3, label = labels[2])
plt.plot([0,4,7,11,14],haplo_w_450_100, linestyle = linestyle[3] , color = 'turquoise', linewidth = 3, label = labels[3])
ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
ax.tick_params(axis='y', which='both' ,right='off', left='on')
plt.xlim(-0.2,14.2)   
plt.xticks([0,4,7,11,14], ['0','40','70','110','140'],fontsize=12, va='top')
plt.yticks(fontsize=12, va='top') 
ax.legend(bbox_to_anchor=(0.477, 1), loc=3,ncol=2, borderaxespad=0.)
ax.set_ylabel('Variance of no. of selected loci per haplotype', size = 14,labelpad=20)

#plot fitness variance
ax=plt.subplot(2,1,2)
plt.plot(np.arange(0,15,1), var_fit_450_10, linestyle = linestyle[0] , color = 'turquoise', linewidth = 3, label = labels[0])
plt.plot(np.arange(0,15,1),var_fit_450_20, linestyle = linestyle[1] , color = 'turquoise', linewidth = 3, label = labels[1])
plt.plot(np.arange(0,15,1),var_fit_450_50, linestyle = linestyle[2] , color = 'turquoise', linewidth = 3, label = labels[2])
plt.plot(np.arange(0,15,1),var_fit_450_100, linestyle = linestyle[3] , color = 'turquoise', linewidth = 3, label = labels[3])
plt.xlim(-0.2,14.2)   
plt.xticks([0,4,7,11,14], ['0','40','70','110','140'],fontsize=12, va='top')
ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
ax.tick_params(axis='y', which='both' ,right='off', left='on')
plt.yticks(fontsize=12, va='top') 
ax.set_ylabel('Variance of fitness', size = 14,labelpad=20)
plt.xlabel('Generation', size = 14,labelpad=15)
plt.tight_layout()
plt.savefig('sweep_Ne450_diffNoLoci_varFitness_varNoLociPerHaplo.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_Ne450_diffNoLoci_varFitness_varNoLociPerHaplo.pdf', dpi=300,format='pdf', bbox_inches = 'tight')
