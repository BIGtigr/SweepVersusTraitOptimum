# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 10:44:19 2020

@author: Neda
"""
import gzip
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

#compute Jaccrd index for a specific timepoint
def compute_Jaccard_index(AF_sim_tmp,sim_reps, tmp,AFC_thresh_neu):
    def get_index(tmp,sim_reps):
        tmp_ind = {10:1, 20:2, 30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10, 110:11, 120:12, 130:13, 140:14}   
        timepoint = 15
        F0_rep_index = list(np.arange(0,timepoint*sim_reps,timepoint))
        Flater_rep_index = list(np.arange(tmp_ind[tmp],timepoint*sim_reps,timepoint))
        get_index = sorted(F0_rep_index+Flater_rep_index)
        return get_index

    #determine the index of AF based on timepoint
    index_toget = get_index(tmp,sim_reps)
    #get the AF based on index
    AF_sim = [[locus[ind] for ind in index_toget] for locus in AF_sim_tmp]

    #convert eary 10 simulations to 1 experiment of 10 replicates 
    num_reps = 10
    rep_index = {}
    n = 0
    for item in np.arange(0,sim_reps*2+1-num_reps*2,num_reps*2):
        rep_index[n] = [item, item+num_reps*2]
        n += 1
    
    sims_AFC = [list() for _ in xrange(int(sim_reps/num_reps))]
    for i in range(len(AF_sim)):
        for rep in rep_index:
            sims_AFC[rep].append(map(operator.sub,[AF_sim[i][ind] for ind in np.arange(rep_index[rep][0]+1,rep_index[rep][1],2)],[AF_sim[i][ind] for ind in np.arange(rep_index[rep][0],rep_index[rep][1],2)]))

    #determine which loci has >neutral in each replicate
    tmp_ind = {10:1, 20:2, 30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10, 110:11, 120:12, 130:13, 140:14}   
    AFC = AFC_thresh_neu[tmp_ind[tmp]-1] 
    ALLsims_risingreps = []
    for s in range(int(sim_reps/num_reps)): 
        onesim_risingreps = [list() for _ in xrange(10)]
        for ind1,i in enumerate(sims_AFC[s]):
            for ind2,item in enumerate(i):
                if item >= AFC: #selected rising
                    onesim_risingreps[ind2].append(ind1)  
        ALLsims_risingreps.append(onesim_risingreps)

    #get Jaccard index
    jaccard_index_ALLsims = []
    for s in range(len(ALLsims_risingreps)):
        jaccard_index_sim = [list() for _ in range(10)]
        for i in range(10):
            for j in range(10):
                total = len(set(ALLsims_risingreps[s][i]+ALLsims_risingreps[s][j]))
                intersect = len(set(ALLsims_risingreps[s][i]).intersection(ALLsims_risingreps[s][j]))
                if total == 0:jaccard_index_sim[i].append(0)
                if total != 0:jaccard_index_sim[i].append(round(intersect/float(total),2))
        jaccard_index_ALLsims.append(jaccard_index_sim)
    
    #remove the repetitive numbers
    #exclude 1s
    reduced_jaccard_index = []
    for ind,item in enumerate(jaccard_index_ALLsims):
        reduced_jaccard_index.append([rep[ind2+1:] for ind2,rep in enumerate(item)])
    
    #get mean and standard deviation for jaccard index for every simulation
    jaccard_mean = np.mean([np.mean([j for rep in sim for j in rep]) for sim in reduced_jaccard_index])
    jaccard_std = np.std([np.mean([j for rep in sim for j in rep]) for sim in reduced_jaccard_index])
    std_jac_down = jaccard_mean-(jaccard_std/2)
    std_jac_up = jaccard_mean+(jaccard_std/2)
    #jaccard_median = np.median([np.median([j for rep in sim for j in rep]) for sim in reduced_jaccard_index])
    return (jaccard_mean,std_jac_up,std_jac_down) 

#this function runs jaccard index for multiple timepoints 
#and computes the mean and standard variation    
def concat_Jaccard_index(AF_sim_tmp,sim_reps, timepoint,AFC_thresh_neu):
    mean_jac, std_jac_up, std_jac_down = [],[],[]
    for tmpnt in timepoint:
        tmp_mean, tmp_std_up, tmp_std_down = compute_Jaccard_index(AF_sim_tmp,sim_reps, tmpnt,AFC_thresh_neu)
        mean_jac.append(tmp_mean)
        std_jac_up.append(tmp_std_up)
        std_jac_down.append(tmp_std_down)
    return (mean_jac, std_jac_up, std_jac_down)

'''
###########
#Figure 3
##########
'''

#sweep 450
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/sel.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/recom_rate_noX/1000sims/AF_results.sync','rb')  

AF_sim_450_w = get_AF_w(InputFile_s,InputFile_af)

#sweep 9k
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.08s/sel.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.08s/recom_rate_noX/500sims/AF_results.sync','rb') 

AF_sim_9k_w = get_AF_w(InputFile_s,InputFile_af)

#QTL 450
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/effect_size.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/recom_rate_noX/1000sims/AF.sync.gz','rb')   

AF_sim_450_qtl = get_AF_qff(InputFile_s,InputFile_af)

#QTL 9k
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--2.5-0.3/effect_size.txt','r')  
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--2.5-0.3/recom_rate_noX/500sims/AF.sync.gz','rb')  

AF_sim_9k_qtl = get_AF_qff(InputFile_s,InputFile_af)

###########
#drift sims
###########

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

########################
#compute Jccard index
########################

mean_Jac_450_neu,std_Jac_450_neu_up,std_Jac_450_neu_down = concat_Jaccard_index(AF_sim_neu_450,500, np.arange(10,150,10),[0]*14)
mean_Jac_9k_neu,std_Jac_9k_neu_up,std_Jac_9k_neu_down = concat_Jaccard_index(AF_sim_neu_9k,500, np.arange(10,150,10),[0]*14)

mean_Jac_450_w,std_Jac_450_w_up,std_Jac_450_w_down = concat_Jaccard_index(AF_sim_450_w,500, np.arange(10,150,10),AFC_thresh_neu_450)
mean_Jac_9k_w,std_Jac_9k_w_up,std_Jac_9k_w_down = concat_Jaccard_index(AF_sim_9k_w,500, np.arange(10,150,10),AFC_thresh_neu_9k)

mean_Jac_450_qtl,std_Jac_450_qtl_up,std_Jac_450_qtl_down = concat_Jaccard_index(AF_sim_450_qtl,500, np.arange(10,150,10),AFC_thresh_neu_450)
mean_Jac_9k_qtl,std_Jac_9k_qtl_up,std_Jac_9k_qtl_down = concat_Jaccard_index(AF_sim_9k_qtl,500, np.arange(10,150,10),AFC_thresh_neu_9k)

######
#plot
######

gen_index = {10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
color_label = ['turquoise','teal','#e699b2','#cc3366']
x = np.arange(1,len(np.arange(10,150,10))+1,1)

fig , ax = plt.subplots(sharex=True, sharey=True,figsize=(8,5),dpi=100, facecolor='w', edgecolor='k')
plt.plot(x,mean_Jac_450_w, color = 'turquoise',linewidth = 3)
plt.fill_between(x, std_Jac_450_w_up, std_Jac_450_w_down, color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
plt.plot(x,mean_Jac_9k_w, color = 'teal',linewidth = 3)
plt.fill_between(x, std_Jac_9k_w_up, std_Jac_9k_w_down, color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
plt.plot(x,mean_Jac_450_qtl, color = '#e699b2',linewidth = 3)
plt.fill_between(x, std_Jac_450_qtl_up, std_Jac_450_qtl_down, color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
plt.plot(x,mean_Jac_9k_qtl, color = '#cc3366',linewidth = 3)
plt.fill_between(x, std_Jac_9k_qtl_up, std_Jac_9k_qtl_down, color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
plt.plot(x,mean_Jac_450_neu, color = 'lightgrey',linewidth = 3)
plt.fill_between(x, std_Jac_450_neu_up, std_Jac_450_neu_down, color = 'lightgrey', alpha = 0.3, edgecolor = None, linewidth=0.0)
plt.plot(x,mean_Jac_9k_neu, color = 'grey',linewidth = 3, linestyle = '--')
plt.fill_between(x, std_Jac_9k_neu_up, std_Jac_9k_neu_down, color = 'lightgrey', alpha = 0.3, edgecolor = None, linewidth=0.0)
ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
ax.tick_params(axis='y', which='both' ,right='off', left='on')
plt.xlim(0,len(x)+1)   
plt.ylim(0,1.1)
plt.xticks(x, [i for i in sorted(gen_index.keys())],fontsize=12, va='top')
plt.yticks(fontsize=12, va='top')
plt.xlabel('Generation', size = 14,labelpad=15)
plt.ylabel('Jaccard similarity index', size = 14,labelpad=15)
texts = ['450 sweep','9000 sweep','450 trait optimum','9000 trait optimum']
patches = [mpatches.Patch(color=color_label[i], label="{:s}".format(texts[i]) ) for i in range(len(texts))]
ax.legend(handles=patches,bbox_to_anchor=(0.421, 1), loc=3,ncol=2, borderaxespad=0.)    
plt.tight_layout()
plt.savefig('sweep_qtl_JaccardIndex_defaultParameters.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_qtl_JaccardIndex_defaultParameters.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''
###########
#Figure S9A
###########
'''

####################################
#get drift and compute Jaccard index
####################################

#450
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.0s/sel.txt','r')   
InputFile_af= gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.0s/recom_rate_noX/500sims/AF_results.sync','rb')

AF_sim_neu_450 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_neu_450, AFC_sim_alltmpnt_neu_450, AFC_thresh_neu_450 = process_AF_sim(AF_sim_neu_450,100,500)
mean_Jac_450_neu,std_Jac_450_neu_up,std_Jac_450_neu_down = concat_Jaccard_index(AF_sim_neu_450,500, np.arange(10,150,10),[0]*14)

#9000
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.0s/sel.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.0s/recom_rate_noX/500sims/AF_results.sync','rb')

AF_sim_neu_9k = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_neu_9k, AFC_sim_alltmpnt_neu_9k, AFC_thresh_neu_9k = process_AF_sim(AF_sim_neu_9k,100,500)
mean_Jac_9k_neu,std_Jac_9k_neu_up,std_Jac_9k_neu_down = concat_Jaccard_index(AF_sim_neu_9k,500, np.arange(10,150,10),[0]*14)

####################################
#get AF and compute Jaccard index
####################################

#qff 450, 10, 20, 50, 100 loci
sim = 'sel_0.5-4.5-0.74-0.3_eff0.04_10loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_450_qtl_10,std_Jac_450_qtl_10_up,std_Jac_450_qtl_10_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

sim = 'sel_0.5-4.5-0.38-0.3_eff0.04_20loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_450_qtl_20,std_Jac_450_qtl_20_up,std_Jac_450_qtl_20_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

sim = 'sel_0.5-4.5--0.7-0.3_eff0.04_50loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_450_qtl_50,std_Jac_450_qtl_50_up,std_Jac_450_qtl_50_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

sim = 'sel_0.5-4.5--2.5-0.3'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/recom_rate_noX/1000sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_450_qtl_100,std_Jac_450_qtl_100_up,std_Jac_450_qtl_100_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

#qff 9000, 10, 20, 50, 100 loci
sim = 'Ne9000_sel_0.5-4.5-0.74-0.3_eff0.04_10loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_9k_qtl_10,std_Jac_9k_qtl_10_up,std_Jac_9k_qtl_10_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

sim = 'Ne9000_sel_0.5-4.5-0.38-0.3_eff0.04_20loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_9k_qtl_20,std_Jac_9k_qtl_20_up,std_Jac_9k_qtl_20_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

sim = 'Ne9000_sel_0.5-4.5--0.7-0.3_eff0.04_50loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_9k_qtl_50,std_Jac_9k_qtl_50_up,std_Jac_9k_qtl_50_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

sim = 'Ne9000_sel_0.5-4.5--2.5-0.3'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/recom_rate_noX/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_9k_qtl_100,std_Jac_9k_qtl_100_up,std_Jac_9k_qtl_100_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

#sweep 450, loci 10, 20, 50, 100
sim = 'Ne450_10loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_450_w_10,std_Jac_450_w_10_up,std_Jac_450_w_10_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

sim =  'Ne450_20loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_450_w_20,std_Jac_450_w_20_up,std_Jac_450_w_20_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

sim = 'Ne450_50loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_450_w_50,std_Jac_450_w_50_up,std_Jac_450_w_50_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

sim =  'Ne450_100loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/1000sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_450_w_100,std_Jac_450_w_100_up,std_Jac_450_w_100_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

#sweep 9000, loci 10, 20, 50, 100
sim = 'Ne9000_10loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_9k_w_10,std_Jac_9k_w_10_up,std_Jac_9k_w_10_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

sim =  'Ne9000_20loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_9k_w_20,std_Jac_9k_w_20_up,std_Jac_9k_w_20_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

sim = 'Ne9000_50loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_9k_w_50,std_Jac_9k_w_50_up,std_Jac_9k_w_50_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

sim =  'Ne9000_100loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_9k_w_100,std_Jac_9k_w_100_up,std_Jac_9k_w_100_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

#######
#plot
#######

mean_Jac = [[mean_Jac_450_w_10,mean_Jac_450_w_20, mean_Jac_450_w_50, mean_Jac_450_w_100,mean_Jac_9k_w_10,mean_Jac_9k_w_20, mean_Jac_9k_w_50, mean_Jac_9k_w_100],
       [mean_Jac_450_qtl_10,mean_Jac_450_qtl_20, mean_Jac_450_qtl_50, mean_Jac_450_qtl_100,mean_Jac_9k_qtl_10,mean_Jac_9k_qtl_20, mean_Jac_9k_qtl_50, mean_Jac_9k_qtl_100]]

std_Jac_up = [[std_Jac_450_w_10_up,std_Jac_450_w_20_up, std_Jac_450_w_50_up, std_Jac_450_w_100_up,std_Jac_9k_w_10_up,std_Jac_9k_w_20_up, std_Jac_9k_w_50_up, std_Jac_9k_w_100_up],
       [std_Jac_450_qtl_10_up,std_Jac_450_qtl_20_up, std_Jac_450_qtl_50_up, std_Jac_450_qtl_100_up,std_Jac_9k_qtl_10_up,std_Jac_9k_qtl_20_up, std_Jac_9k_qtl_50_up, std_Jac_9k_qtl_100_up]]

std_Jac_down = [[std_Jac_450_w_10_down,std_Jac_450_w_20_down, std_Jac_450_w_50_down, std_Jac_450_w_100_down,std_Jac_9k_w_10_down,std_Jac_9k_w_20_down, std_Jac_9k_w_50_down, std_Jac_9k_w_100_down],
       [std_Jac_450_qtl_10_down,std_Jac_450_qtl_20_down, std_Jac_450_qtl_50_down, std_Jac_450_qtl_100_down,std_Jac_9k_qtl_10_down,std_Jac_9k_qtl_20_down, std_Jac_9k_qtl_50_down, std_Jac_9k_qtl_100_down]]

gen_index = {10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
color_label = ['turquoise']*4+['teal']*4+['#e699b2']*4+['#cc3366']*4
color_shade = ['grey']*4+['teal']*4+['grey']*4+['#cc3366']*4
alpha = [0.3]*4+[0.1]*4+[0.3]*4+[0.1]*4

linestyle = [':','-.','--','-']*4
labels = ['450, 10 loci', '450, 20 loci','450, 50 loci','450, 100 loci','9000, 10 loci', '9000, 20 loci','9000, 50 loci','9000, 100 loci']*2
x = np.arange(1,len(np.arange(10,150,10))+1,1)

fig , ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True,figsize=(18,5),dpi=100, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.3, wspace=0.1)
n = 0
for rep in range(2):
    m = 0
    ax=plt.subplot(1,2,rep+1)
    plt.plot(x,mean_Jac_450_neu, color = 'gold',linewidth = 3)
    plt.fill_between(x, std_Jac_450_neu_up, std_Jac_450_neu_down, color = 'gold', alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(x,mean_Jac_9k_neu, color = 'gold',linewidth = 3, linestyle = '--')
    plt.fill_between(x, std_Jac_9k_neu_up, std_Jac_9k_neu_down, color = 'gold', alpha = 0.3, edgecolor = None, linewidth=0.0)
    for item in mean_Jac[rep]:    
        plt.plot(x,item, color = color_label[n],linewidth = 3,linestyle =linestyle[n], label = labels[n])
        plt.fill_between(x, std_Jac_up[rep][m], std_Jac_down[rep][m], color = color_shade[n], alpha = alpha[n], edgecolor = None, linewidth=0.0)
        n+=1
        m+=1
    ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
    ax.tick_params(axis='y', which='both' ,right='off', left='on')
    plt.xlim(0,len(x)+1)   
    plt.ylim(-0.1,1.1)
    plt.xticks(x, [i for i in sorted(gen_index.keys())],fontsize=12, va='top')
    plt.yticks(fontsize=12, va='top')
    plt.xlabel('Generation', size = 14,labelpad=15)
    if rep == 0: plt.ylabel('Jaccard similarity index', size = 14,labelpad=15)
    ax.legend(bbox_to_anchor=(0.113, 1), loc=3,ncol=4, borderaxespad=0.)    
plt.tight_layout()
plt.savefig('sweep_qtl_DiffLociNum_JaccardIndex.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_qtl_DiffLociNum_JaccardIndex.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''
###########
#figure S9B
###########
'''

####################################
#get drift and compute Jaccard index
####################################

#450
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.0s/sel.txt','r')   
InputFile_af= gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.0s/recom_rate_noX/500sims/AF_results.sync','rb')

AF_sim_neu_450 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_neu_450, AFC_sim_alltmpnt_neu_450, AFC_thresh_neu_450 = process_AF_sim(AF_sim_neu_450,100,500)
mean_Jac_450_neu,std_Jac_450_neu_up,std_Jac_450_neu_down = concat_Jaccard_index(AF_sim_neu_450,500, np.arange(10,150,10),[0]*14)

#9000
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.0s/sel.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.0s/recom_rate_noX/500sims/AF_results.sync','rb')

AF_sim_neu_9k = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_neu_9k, AFC_sim_alltmpnt_neu_9k, AFC_thresh_neu_9k = process_AF_sim(AF_sim_neu_9k,100,500)
mean_Jac_9k_neu,std_Jac_9k_neu_up,std_Jac_9k_neu_down = concat_Jaccard_index(AF_sim_neu_9k,500, np.arange(10,150,10),[0]*14)

####################################
#get AF and compute Jaccard index
####################################

#sweep 450 0.02,0.05,0.08,0.1 s
sim = 'Ne450_100loci_0.05p0_0.02s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_450_w_02,std_Jac_450_w_02_up,std_Jac_450_w_02_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

sim =  'Ne450_100loci_0.05p0_0.05s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_450_w_05,std_Jac_450_w_05_up,std_Jac_450_w_05_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

sim = 'Ne450_100loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/1000sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_450_w_08,std_Jac_450_w_08_up,std_Jac_450_w_08_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

sim =  'Ne450_100loci_0.05p0_0.1s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_450_w_1,std_Jac_450_w_1_up,std_Jac_450_w_1_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

#sweep 9000, 0.02,0.05,0.08,0.1 s
sim = 'Ne9000_100loci_0.05p0_0.02s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_9k_w_02,std_Jac_9k_w_02_up,std_Jac_9k_w_02_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

sim =  'Ne9000_100loci_0.05p0_0.05s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_9k_w_05,std_Jac_9k_w_05_up,std_Jac_9k_w_05_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

sim = 'Ne9000_100loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_9k_w_08,std_Jac_9k_w_08_up,std_Jac_9k_w_08_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

sim =  'Ne9000_100loci_0.05p0_0.1s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim = get_AF_w(InputFile_s,InputFile_af)
mean_Jac_9k_w_1,std_Jac_9k_w_1_up,std_Jac_9k_w_1_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

#qff 450, 100 loci, different effect size: 0.04, 0.08, 0.2, 0.4 
sim = 'sel_0.5-4.5--2.5-0.3'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/recom_rate_noX/1000sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_450_qtl_04,std_Jac_450_qtl_04_up,std_Jac_450_qtl_04_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

sim = 'sel_0.5-4.5--6.1-0.3_eff0.08_100loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_450_qtl_08,std_Jac_450_qtl_08_up,std_Jac_450_qtl_08_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

sim = 'sel_0.5-4.5--16.9-0.3_eff0.2_100loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_450_qtl_2,std_Jac_450_qtl_2_up,std_Jac_450_qtl_2_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

sim = 'sel_0.5-4.5--34.9-0.3_eff0.4_100loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_450_qtl_4,std_Jac_450_qtl_4_up,std_Jac_450_qtl_4_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_450)

#qff 9000, 100 loci, different effect size: 0.04, 0.08, 0.2, 0.4 
sim = 'Ne9000_sel_0.5-4.5--2.5-0.3'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/recom_rate_noX/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_9k_qtl_04,std_Jac_9k_qtl_04_up,std_Jac_9k_qtl_04_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

sim = 'Ne9000_sel_0.5-4.5--6.1-0.3_eff0.08_100loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_9k_qtl_08,std_Jac_9k_qtl_08_up,std_Jac_9k_qtl_08_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

sim = 'Ne9000_sel_0.5-4.5--16.9-0.3_eff0.2_100loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_9k_qtl_2,std_Jac_9k_qtl_2_up,std_Jac_9k_qtl_2_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

sim = 'Ne9000_sel_0.5-4.5--34.9-0.3_eff0.4_100loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim = get_AF_qff(InputFile_s,InputFile_af)
mean_Jac_9k_qtl_4,std_Jac_9k_qtl_4_up,std_Jac_9k_qtl_4_down = concat_Jaccard_index(AF_sim,500, np.arange(10,150,10),AFC_thresh_neu_9k)

#######
#plot
#######

mean_Jac = [[mean_Jac_450_w_02,mean_Jac_450_w_05, mean_Jac_450_w_08, mean_Jac_450_w_1,mean_Jac_9k_w_02,mean_Jac_9k_w_05, mean_Jac_9k_w_08, mean_Jac_9k_w_1],
       [mean_Jac_450_qtl_04,mean_Jac_450_qtl_08, mean_Jac_450_qtl_2, mean_Jac_450_qtl_4,mean_Jac_9k_qtl_04,mean_Jac_9k_qtl_08, mean_Jac_9k_qtl_2, mean_Jac_9k_qtl_4]]

std_Jac_up = [[std_Jac_450_w_02_up,std_Jac_450_w_05_up, std_Jac_450_w_08_up, std_Jac_450_w_1_up,std_Jac_9k_w_02_up,std_Jac_9k_w_05_up, std_Jac_9k_w_08_up, std_Jac_9k_w_1_up],
       [std_Jac_450_qtl_04_up,std_Jac_450_qtl_08_up, std_Jac_450_qtl_2_up, std_Jac_450_qtl_4_up,std_Jac_9k_qtl_04_up,std_Jac_9k_qtl_08_up, std_Jac_9k_qtl_2_up, std_Jac_9k_qtl_4_up]]

std_Jac_down = [[std_Jac_450_w_02_down,std_Jac_450_w_05_down, std_Jac_450_w_08_down, std_Jac_450_w_1_down,std_Jac_9k_w_02_down,std_Jac_9k_w_05_down, std_Jac_9k_w_08_down, std_Jac_9k_w_1_down],
       [std_Jac_450_qtl_04_down,std_Jac_450_qtl_08_down, std_Jac_450_qtl_2_down, std_Jac_450_qtl_4_down,std_Jac_9k_qtl_04_down,std_Jac_9k_qtl_08_down, std_Jac_9k_qtl_2_down, std_Jac_9k_qtl_4_down]]

gen_index = {10:1,20:2,30:3,40:4,50:5,60:6,70:7,80:8,90:9,100:10,110:11,120:12,130:13,140:14}
color_label = ['turquoise']*4+['teal']*4+['#e699b2']*4+['#cc3366']*4
linestyle = [':','-.','-','--']*2+['-',':','-.','--']*2
labels = ['450, s=0.02', '450, s=0.05','450, s=0.08','450, s=0.1','9000, s=0.02', '9000, s=0.05','9000, s=0.08','9000, s=0.1',
          '450, eff_size=0.04', '450, eff_size=0.08','450, eff_size=0.2','450, eff_size=0.4','9000, eff_size=0.04', '9000, eff_size=0.08','9000, eff_size=0.2','9000, eff_size=0.4']
color_shade = ['grey']*4+['teal']*4+['grey']*4+['#cc3366']*4
alpha = [0.3]*4+[0.1]*4+[0.3]*4+[0.1]*4
x = np.arange(1,len(np.arange(10,150,10))+1,1)

fig , ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True,figsize=(18,5),dpi=100, facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace=0.3, wspace=0.1)
n,m = 0,0
for rep in range(2):
    m=0
    ax=plt.subplot(1,2,rep+1)
    plt.plot(x,mean_Jac_450_neu, color = 'gold',linewidth = 3)
    plt.fill_between(x, std_Jac_450_neu_up, std_Jac_450_neu_down, color = 'gold', alpha = 0.3, edgecolor = None, linewidth=0.0)
    plt.plot(x,mean_Jac_9k_neu, color = 'gold',linewidth = 3, linestyle = '--')
    plt.fill_between(x, std_Jac_9k_neu_up, std_Jac_9k_neu_down, color = 'gold', alpha = 0.3, edgecolor = None, linewidth=0.0)
    for item in mean_Jac[rep]:    
        plt.plot(x,item, color = color_label[n],linewidth = 3,linestyle =linestyle[n], label = labels[n])
        plt.fill_between(x, std_Jac_up[rep][m], std_Jac_down[rep][m], color = color_shade[n], alpha = alpha[n], edgecolor = None, linewidth=0.0)
        n+=1
        m+=1
    ax.tick_params(axis='x', which='both' ,bottom='on', top='off')
    ax.tick_params(axis='y', which='both' ,right='off', left='on')
    plt.xlim(0,len(x)+1)  
    plt.ylim(-0.1,1.1)
    plt.xticks(x, [i for i in sorted(gen_index.keys())],fontsize=12, va='top')
    plt.yticks(fontsize=12, va='top')
    plt.xlabel('Generation', size = 14,labelpad=15)
    if rep == 0: 
        plt.ylabel('Jaccard similarity index', size = 14,labelpad=15)
        ax.legend(bbox_to_anchor=(0.337, 1), loc=3,ncol=3, borderaxespad=0.)
    if rep == 1: ax.legend(bbox_to_anchor=(0.153, 1), loc=3,ncol=3, borderaxespad=0.)
plt.tight_layout()
plt.savefig('sweep_qtl_DiffEffSize_JaccardIndex.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_qtl_DiffEffSize_JaccardIndex.pdf', dpi=300,format='pdf', bbox_inches = 'tight')
