# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 15:24:03 2020

@author: Neda
"""
import gzip
import operator
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import pchip

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

#this funcgtion will convert the output of functions get_AF_w or get_AF_qff (each item is one locus and for each locus we have replicate*timepoint entries)
# to AF_sim_rep where each item (timepoint) has separated replicates based on replicate
def process_AF_SimBatch(AF_sim,loci_num,replicates):
    #classify the frequency of different loci and replicates based on timepoints
    timepoints = 15   #I know it's hard coded, you can change it if you need. 
    timepoint_index = []
    for i in range(0,timepoints):
        timepoint_index.append(np.arange(i,(replicates*timepoints),timepoints))
    
    AF_sim_classified = [list() for _ in xrange(timepoints)]
    for item in AF_sim: #each item in AF_sim is one locus and for each locus we have replicate*timepoint entries
        for ind,timepoint in enumerate(timepoint_index):
            for t in timepoint:
                AF_sim_classified[ind].append(item[t])
    
    #to classify the AF for every timepoint based on separate sim
    #first divide 5000 data points for each timepoint (100loci * 500 sims) into separate sims  
    sim_index = []
    for i in np.arange(0,loci_num*replicates,loci_num):
        sim_index.append([i,i+loci_num])
    
    #for each item (timepoint) separate replicates based on replicate
    AF_sims_reps=[]
    for item in AF_sim_classified: #each item in AF_sim_classified is one timepoint with loci_num*replicates entries
        oneRep = [item[s[0]:s[1]] for s in sim_index]
        AF_sims_reps.append(oneRep) 
    
    return(AF_sim_classified,AF_sims_reps)

def process_AF_sim(AF_sim,loci_num,replicates):
    #classify the frequency of different loci and replicates based on generations/timepoints
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

#this function will keep the loci with AFC more than neutrality. 
#the input data is allele frequency NOT allele frequency change 
#this will also count rising and swept loci
def AF_driftCorrected_RisingSweptLoci(AF_sims_reps,AFC_thresh_neu):
    AF_drifted,rising_all,swept_all = [AF_sims_reps[0]],[],[]
    for ind,timepoint in enumerate(AF_sims_reps[1:]):
        drifted_tmp,rising_tmp, swept_tmp = [],[],[]
        for ind2,rep in enumerate(timepoint):
            drifted_rep = [l for ind3,l in enumerate(rep) if l >= AF_sims_reps[0][ind2][ind3]+AFC_thresh_neu[ind]]  
            rising,swept = len(drifted_rep), len([i for i in drifted_rep if i >= 0.9])
            rising_tmp.append(rising)
            swept_tmp.append(swept)
            drifted_tmp.append(drifted_rep) 
        AF_drifted.append(drifted_tmp) 
        rising_all.append(rising_tmp)    
        swept_all.append(swept_tmp)
    return AF_drifted,rising_all,swept_all

#to divide a list of items into bins
def hist_prop(freqs, bins):#freqs - list of floats; bins is a list or array of bin boundaries
    histolist = []
    if len(freqs) == 0:histolist=[0]*(len(bins)-1)
    if len(freqs) > 0:
        for i in range(len(bins[1:])):
            histolist.append(len([j for j in freqs if bins[i] < j <= bins [i+1]])/float(len(freqs)))
    return histolist

#this function distribute the values of each replicate into bins and gets the mean and standard deviation of all replicates for each bin 
# bins is a list or array of bin boundaries
def mean_binned_AF(AF_sim_rep_drifted,bins):
    #get the number of loci in each bin for all timepoints and replicates  
    hist_AF = []
    for ind,tmp in enumerate(AF_sim_rep_drifted[1:]):
        hist_tmp = []    
        for ind2,rep in enumerate(tmp):
            hist_tmp.append(hist_prop(rep,bins))
        hist_AF.append(hist_tmp)
    
    #for each timepoint transpose the data so that for each bin a total of replicates item will be present 
    trans_hist_AF = [np.transpose(tmp) for tmp in hist_AF]    
    #then for each timepoint get the mean for each bin
    mean_hist_AF=[[np.mean(tmp[i]) for i in range(len(bins)-1)] for tmp in trans_hist_AF]
    #for each timepoint compute standard deviation of the mean for each bin
    std_hist_AF=[[np.std(tmp[i]) for i in range(len(bins)-1)] for tmp in trans_hist_AF]
    #compute the area around the mean for plotting std
    std_hist_AF_up, std_hist_AF_down = [], []
    for ind,tmp in enumerate(mean_hist_AF):
        std_hist_AF_down.append([i-(std_hist_AF[ind][ind2]/2) for ind2,i in enumerate(tmp)])
        std_hist_AF_up.append([i+(std_hist_AF[ind][ind2]/2) for ind2,i in enumerate(tmp)])
    return (mean_hist_AF,std_hist_AF,std_hist_AF_up,std_hist_AF_down)    

# this function masks the axes
def mask_XYaxes():
    if i in np.arange(6,51,4) or i in np.arange(7,52,4):
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='off', top='off', labelbottom = False)
        ax1.tick_params(axis='y', which='both' ,right='off', left='off', labelleft = False)
    if i in np.arange(5,50,4):
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='off', top='off', labelbottom = False)
        ax1.tick_params(axis='y', which='both' ,right='off', left='on')
    if i in np.arange(8,53,4):
        ax1.spines['top'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='off', top='off', labelbottom = False)
        ax1.tick_params(axis='y', which='both' ,right='off', left='off', labelleft = False)
    if i == 1:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(True)
        ax1.tick_params(axis='x', which='both' ,bottom='off', top='off', labelbottom = False)
        ax1.tick_params(axis='y', which='both' ,right='off', left='on')
    if i == 2 or i == 3:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='off', top='off', labelbottom = False)
        ax1.tick_params(axis='y', which='both' ,right='off', left='off', labelleft = False)
    if i == 4:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='off', top='off', labelbottom = False)
        ax1.tick_params(axis='y', which='both' ,right='off', left='off', labelleft = False)
    if i == 53:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(True)
        ax1.spines['left'].set_visible(True)
        ax1.tick_params(axis='x', which='both' ,bottom='on', top='off')
        ax1.tick_params(axis='y', which='both' ,right='off', left='on')
    if i == 54 or i == 55:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(True)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='on', top='off')
        ax1.tick_params(axis='y', which='both' ,right='off', left='off', labelleft = False)
    if i == 56:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(True)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='on', top='off')
        ax1.tick_params(axis='y', which='both' ,right='off', left='off', labelleft = False)

# this function masks the axes        
def mask_XYaxes_sub():
    if i in [1,5,9]:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(True)
        ax1.tick_params(axis='x', which='both' ,bottom='off', top='off', labelbottom = False)
        ax1.tick_params(axis='y', which='both' ,right='off', left='on')
    if i in [14,15,16]:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(True)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='on', top='off')
        ax1.tick_params(axis='y', which='both' ,right='off', left='off', labelleft = False)
    if i == 13:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(True)
        ax1.spines['left'].set_visible(True)
        ax1.tick_params(axis='x', which='both' ,bottom='on', top='off')
        ax1.tick_params(axis='y', which='both' ,right='off', left='on')
    if i in [2,3,4,6,7,8,10,11,12]:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='off', top='off', labelbottom = False)
        ax1.tick_params(axis='y', which='both' ,right='off', left='off', labelleft = False)

'''        
###################
#A. figure 1 and S2
###################
'''

#####################
#get simulation data
#####################

#sweep 450
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/sel.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne450_100loci_0.05p0_0.08s/recom_rate_noX/1000sims/AF_results.sync','rb')  

AF_sim_450_w = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_450_w, AF_sim_rep_450_w = process_AF_SimBatch(AF_sim_450_w,100,500)

#sweep 9k
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.08s/sel.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/Ne9000_100loci_0.05p0_0.08s/recom_rate_noX/500sims/AF_results.sync','rb') 

AF_sim_9k_w = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_9k_w, AF_sim_rep_9k_w = process_AF_SimBatch(AF_sim_9k_w,100,500)

#qtl 450
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/effect_size.txt','r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/sel_0.5-4.5--2.5-0.3/recom_rate_noX/1000sims/AF.sync.gz','rb')   

AF_sim_450_qtl = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_450_qtl, AF_sim_rep_450_qtl = process_AF_SimBatch(AF_sim_450_qtl,100,500)

#qtl 9k
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--2.5-0.3/effect_size.txt','r')  
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/Ne9000_sel_0.5-4.5--2.5-0.3/recom_rate_noX/500sims/AF.sync.gz','rb')  

AF_sim_9k_qtl = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_9k_qtl, AF_sim_rep_9k_qtl = process_AF_SimBatch(AF_sim_9k_qtl,100,500)

################
#get drift data
################

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
    
################################
#correct for drift
#get # of rising and sweep loci
################################

#sweep
AF_sim_rep_450_w_drifted, w_450_NumRisingLoci, w_450_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_w,AFC_thresh_neu_450)
AF_sim_rep_9k_w_drifted, w_9k_NumRisingLoci, w_9k_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_w,AFC_thresh_neu_9k)

#qtl
AF_sim_rep_450_qtl_drifted, qtl_450_NumRisingLoci, qtl_450_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_qtl,AFC_thresh_neu_450)
AF_sim_rep_9k_qtl_drifted, qtl_9k_NumRisingLoci, qtl_9k_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_qtl,AFC_thresh_neu_9k)

#compute median AF
median_AF_sim_rep_450_w_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_w_drifted[1:]]   
median_AF_sim_rep_9k_w_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_w_drifted[1:]]   

median_AF_sim_rep_450_qtl_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_qtl_drifted[1:]]   
median_AF_sim_rep_9k_qtl_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_qtl_drifted[1:]]   

######################################################
#get binned mean and standard deviation of replicates
######################################################

mean_binned_450_w,std_binned_450_w,std_binned_450_w_up,std_binned_450_w_down = mean_binned_AF(AF_sim_rep_450_w_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_w,std_binned_9k_w,std_binned_9k_w_up,std_binned_9k_w_down =mean_binned_AF(AF_sim_rep_9k_w_drifted,bins=np.arange(0,1.05,0.05))

mean_binned_450_qtl,std_binned_450_qtl,std_binned_450_qtl_up,std_binned_450_qtl_down =mean_binned_AF(AF_sim_rep_450_qtl_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_qtl,std_binned_9k_qtl,std_binned_9k_qtl_up,std_binned_9k_qtl_down =mean_binned_AF(AF_sim_rep_9k_qtl_drifted,bins=np.arange(0,1.05,0.05))

################
#plot Figure 1
###############

#plot selected generations F10, F50, F100, F140
ind_toplot = [0,4,9,13]
#check how to set max value of Y axis
for i in ind_toplot:
    print max([max(std_binned_450_w_up[i]),max(std_binned_9k_w_up[i]),max(std_binned_450_qtl_up[i]),max(std_binned_9k_qtl_up[i])])

yMax=[1,0.5,0.5,0.5] 
timepoints = len(ind_toplot)
gen_label =  ['Gen %i'%tm for tm in np.arange(10,150,10)]
nbins = 20
bin_edge =np.array([float(i)/nbins for i in range(nbins+1)])
x = 0.5*(bin_edge[1:]+bin_edge[:-1])
xx = np.linspace(0, 1, 300)

fig, ax = plt.subplots(nrows=4, ncols=2, figsize=(10,10), facecolor='w', edgecolor='k')  
fig.subplots_adjust(hspace=0.15, wspace=0.1)
#first plot sweep
for ind,i in enumerate(np.arange(1,timepoints*2,2)):
    ax1=plt.subplot(4,2,i)
    pch_up = pchip(x,std_binned_450_w_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_450_w_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, label = '450 sweep')
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_drifted[ind_toplot[ind]]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_9k_w_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle = '--', label = '9000 sweep')
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_drifted[ind_toplot[ind]]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    if round(np.mean(w_450_NumSweep[ind_toplot[ind]])) == 0 and round(np.mean(w_9k_NumSweep[ind_toplot[ind]])) == 0:
        ax1.annotate('%i' %int(round(np.mean(w_450_NumRisingLoci[ind_toplot[ind]]))), xy=(0.7, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(w_9k_NumRisingLoci[ind_toplot[ind]]))), xy=(0.7, 0.45),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if round(np.mean(w_450_NumSweep[ind_toplot[ind]])) != 0 or round(np.mean(w_9k_NumSweep[ind_toplot[ind]])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(w_450_NumSweep[ind_toplot[ind]])))), xy=(0.7, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(w_9k_NumSweep[ind_toplot[ind]])))), xy=(0.7, 0.45),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
    if i < 7:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(True)
        ax1.tick_params(axis='x', which='both' ,bottom='off', top='off', labelbottom = False)
        ax1.tick_params(axis='y', which='both' ,right='off', left='on')
    if i == 1 :
        ax1.legend(bbox_to_anchor=(0.62, 1), loc=3,ncol=1, borderaxespad=0.) 
    if i == 7:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='on', top='off')
        ax1.tick_params(axis='y', which='both' ,right='off', left='on')
#then plot qtl
for ind,i in enumerate(np.arange(2,timepoints*2+1,2)):
    ax1=plt.subplot(4,2,i)
    pch_up = pchip(x,std_binned_450_qtl_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_450_qtl_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, label = '450 trait optimum')
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_drifted[ind_toplot[ind]]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_9k_qtl_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle = '--', label = '9000 trait optimum')
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_drifted[ind_toplot[ind]]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    ax1.annotate(gen_label[ind_toplot[ind]], xy=(0.85, 0.8), xytext=(0.85, 0.8),textcoords='axes fraction',xycoords='axes fraction', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_NumSweep[ind_toplot[ind]])) == 0 and round(np.mean(qtl_9k_NumSweep[ind_toplot[ind]])) == 0:
        ax1.annotate('%i' %int(round(np.mean(qtl_450_NumRisingLoci[ind_toplot[ind]]))), xy=(0.8, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(qtl_9k_NumRisingLoci[ind_toplot[ind]]))), xy=(0.8, 0.45),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_NumSweep[ind_toplot[ind]])) != 0 or round(np.mean(qtl_9k_NumSweep[ind_toplot[ind]])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(qtl_450_NumSweep[ind_toplot[ind]])))), xy=(0.8, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(qtl_9k_NumSweep[ind_toplot[ind]])))), xy=(0.8, 0.45),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
    plt.yticks([])
    if i < 7:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='off', top='off', labelbottom = False)
        ax1.tick_params(axis='y', which='both' ,right='off', left='on')
    if i == 2:
        color_label = ['#e699b2','#cc3366']
        ax1.legend(bbox_to_anchor=(0.475, 1), loc=3,ncol=1, borderaxespad=0.)   
    if i == 8:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='on', top='off')
        ax1.tick_params(axis='y', which='both' ,right='off', left='off')
fig.text(0.05, 0.5, 'Proportion', va='center', rotation='vertical',size = 18) #ylabel Count Probability
fig.text(0.5, 0.05, 'Allele frequency', ha='center',size = 18) #xlabel
plt.savefig('sweep_qtl_AF_RisingLociNum_smoothCurve_defaultParameters_sub.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_qtl_AF_RisingLociNum_smoothCurve_defaultParameters_sub.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

################
#plot figure S2 
###############

#check how to set max value of Y axis
for i in range(14):
    print max([max(std_binned_450_w_up[i]),max(std_binned_9k_w_up[i]),max(std_binned_450_qtl_up[i]),max(std_binned_9k_qtl_up[i])])

yMax=[1,0.8,0.6]+[0.5]*11 
timepoints = 14
gen_label = ['Gen %i'%tm for tm in np.arange(10,150,10)]
nbins = 20
bin_edge =np.array([float(i)/nbins for i in range(nbins+1)])
x = 0.5*(bin_edge[1:]+bin_edge[:-1])
xx = np.linspace(0, 1, 300)

fig, ax = plt.subplots(nrows=14, ncols=2, figsize=(10,25), facecolor='w', edgecolor='k')  
fig.subplots_adjust(hspace=0.15, wspace=0.1)
#first plot sweep
for ind,i in enumerate(np.arange(1,timepoints*2,2)):
    ax1=plt.subplot(14,2,i)
    pch_up = pchip(x,std_binned_450_w_up[ind])
    pch_down = pchip(x,std_binned_450_w_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w[ind])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, label = '450 sweep')
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_drifted[ind]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_up[ind])
    pch_down = pchip(x,std_binned_9k_w_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w[ind])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle = '--', label = '9000 sweep')
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_drifted[ind]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    if round(np.mean(w_450_NumSweep[ind])) == 0 and round(np.mean(w_9k_NumSweep[ind])) == 0:
        ax1.annotate('%i' %int(round(np.mean(w_450_NumRisingLoci[ind]))), xy=(0.7, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(w_9k_NumRisingLoci[ind]))), xy=(0.7, 0.45),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if round(np.mean(w_450_NumSweep[ind])) != 0 or round(np.mean(w_9k_NumSweep[ind])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_NumRisingLoci[ind]))),int(round(np.mean(w_450_NumSweep[ind])))), xy=(0.7, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_NumRisingLoci[ind]))),int(round(np.mean(w_9k_NumSweep[ind])))), xy=(0.7, 0.45),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
    if i < 27:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(True)
        ax1.tick_params(axis='x', which='both' ,bottom='off', top='off', labelbottom = False)
        ax1.tick_params(axis='y', which='both' ,right='off', left='on')
    if i == 1 :
        ax1.legend(bbox_to_anchor=(0.64, 1), loc=3,ncol=1, borderaxespad=0.) 
    if i == 27:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='on', top='off')
        ax1.tick_params(axis='y', which='both' ,right='off', left='on')
#then plot qtl
for ind,i in enumerate(np.arange(2,timepoints*2+1,2)):
    ax1=plt.subplot(14,2,i)
    pch_up = pchip(x,std_binned_450_qtl_up[ind])
    pch_down = pchip(x,std_binned_450_qtl_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl[ind])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, label = '450 trait optimum')
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_drifted[ind]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_up[ind])
    pch_down = pchip(x,std_binned_9k_qtl_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl[ind])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle = '--', label = '9000 trait optimum')
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_drifted[ind]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    ax1.annotate(gen_label[ind], xy=(0.85, 0.8), xytext=(0.85, 0.8),textcoords='axes fraction',xycoords='axes fraction', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_NumSweep[ind])) == 0 and round(np.mean(qtl_9k_NumSweep[ind])) == 0:
        ax1.annotate('%i' %int(round(np.mean(qtl_450_NumRisingLoci[ind]))), xy=(0.8, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(qtl_9k_NumRisingLoci[ind]))), xy=(0.8, 0.45),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_NumSweep[ind])) != 0 or round(np.mean(qtl_9k_NumSweep[ind])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_NumRisingLoci[ind]))),int(round(np.mean(qtl_450_NumSweep[ind])))), xy=(0.8, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_NumRisingLoci[ind]))),int(round(np.mean(qtl_9k_NumSweep[ind])))), xy=(0.8, 0.45),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
    plt.yticks([])
    if i < 27:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='off', top='off', labelbottom = False)
        ax1.tick_params(axis='y', which='both' ,right='off', left='on')
    if i == 2:
        color_label = ['#e699b2','#cc3366']
        ax1.legend(bbox_to_anchor=(0.475, 1), loc=3,ncol=1, borderaxespad=0.)   
    if i == 28:
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_visible(False)
        ax1.tick_params(axis='x', which='both' ,bottom='on', top='off')
        ax1.tick_params(axis='y', which='both' ,right='off', left='off')
fig.text(0.05, 0.5, 'Proportion', va='center', rotation='vertical',size = 18) #ylabel Count Probability
fig.text(0.5, 0.1, 'Allele frequency', ha='center',size = 18) #xlabel
plt.savefig('sweep_qtl_AF_RisingLociNum_smoothCurve_defaultParameters.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_qtl_AF_RisingLociNum_smoothCurve_defaultParameters.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''                
##########
#Figure S6
##########
'''

################
#get drift data
###############

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

##################
#get simulation AF 
##################

#sweep 450, loci 10, 20, 50, 100
sim = 'Ne450_10loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_450_w_10 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_450_w_10, AF_sim_rep_450_w_10 = process_AF_SimBatch(AF_sim_450_w_10,10,500)

sim =  'Ne450_20loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_450_w_20 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_450_w_20, AF_sim_rep_450_w_20 = process_AF_SimBatch(AF_sim_450_w_20,20,500)

sim = 'Ne450_50loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_450_w_50 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_450_w_50, AF_sim_rep_450_w_50 = process_AF_SimBatch(AF_sim_450_w_50,50,500)

sim =  'Ne450_100loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/1000sims/AF_results.sync' %sim,'rb')  
AF_sim_450_w_100 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_450_w_100, AF_sim_rep_450_w_100 = process_AF_SimBatch(AF_sim_450_w_100,100,500)

#sweep 9000, loci 10, 20, 50, 100
sim = 'Ne9000_10loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_9k_w_10 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_9k_w_10, AF_sim_rep_9k_w_10 = process_AF_SimBatch(AF_sim_9k_w_10,10,500)

sim =  'Ne9000_20loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_9k_w_20 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_9k_w_20, AF_sim_rep_9k_w_20 = process_AF_SimBatch(AF_sim_9k_w_20,20,500)

sim = 'Ne9000_50loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_9k_w_50 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_9k_w_50, AF_sim_rep_9k_w_50 = process_AF_SimBatch(AF_sim_9k_w_50,50,500)

sim =  'Ne9000_100loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_9k_w_100 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_9k_w_100, AF_sim_rep_9k_w_100 = process_AF_SimBatch(AF_sim_9k_w_100,100,500)

################################
#correct for drift
#get # of rising and sweep loci
################################

#sweep 450, loci 10, 20, 50, 100
AF_sim_rep_450_w_10_drifted, w_450_10_NumRisingLoci, w_450_10_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_w_10,AFC_thresh_neu_450)
AF_sim_rep_450_w_20_drifted, w_450_20_NumRisingLoci, w_450_20_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_w_20,AFC_thresh_neu_450)
AF_sim_rep_450_w_50_drifted, w_450_50_NumRisingLoci, w_450_50_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_w_50,AFC_thresh_neu_450)
AF_sim_rep_450_w_100_drifted, w_450_100_NumRisingLoci, w_450_100_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_w_100,AFC_thresh_neu_450)

#sweep 9k, loci 10, 20, 50, 100
AF_sim_rep_9k_w_10_drifted, w_9k_10_NumRisingLoci, w_9k_10_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_w_10,AFC_thresh_neu_9k)
AF_sim_rep_9k_w_20_drifted, w_9k_20_NumRisingLoci, w_9k_20_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_w_20,AFC_thresh_neu_9k)
AF_sim_rep_9k_w_50_drifted, w_9k_50_NumRisingLoci, w_9k_50_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_w_50,AFC_thresh_neu_9k)
AF_sim_rep_9k_w_100_drifted, w_9k_100_NumRisingLoci, w_9k_100_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_w_100,AFC_thresh_neu_9k)

#compute median AF
#sweep 450, loci 10, 20, 50, 100
median_AF_sim_rep_450_w_10_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_w_10_drifted[1:]]   
median_AF_sim_rep_450_w_20_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_w_20_drifted[1:]]   
median_AF_sim_rep_450_w_50_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_w_50_drifted[1:]]   
median_AF_sim_rep_450_w_100_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_w_100_drifted[1:]]   

#sweep 9k, loci 10, 20, 50, 100
median_AF_sim_rep_9k_w_10_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_w_10_drifted[1:]]   
median_AF_sim_rep_9k_w_20_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_w_20_drifted[1:]]   
median_AF_sim_rep_9k_w_50_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_w_50_drifted[1:]]   
median_AF_sim_rep_9k_w_100_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_w_100_drifted[1:]]   

######################################################
#get binned mean and standard deviation of replicates
######################################################

#sweep 450, loci 10, 20, 50, 100
mean_binned_450_w_10,std_binned_450_w_10,std_binned_450_w_10_up,std_binned_450_w_10_down = mean_binned_AF(AF_sim_rep_450_w_10_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_450_w_20,std_binned_450_w_20,std_binned_450_w_20_up,std_binned_450_w_20_down = mean_binned_AF(AF_sim_rep_450_w_20_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_450_w_50,std_binned_450_w_50,std_binned_450_w_50_up,std_binned_450_w_50_down = mean_binned_AF(AF_sim_rep_450_w_50_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_450_w_100,std_binned_450_w_100,std_binned_450_w_100_up,std_binned_450_w_100_down = mean_binned_AF(AF_sim_rep_450_w_100_drifted,bins=np.arange(0,1.05,0.05))

#sweep 9k, loci 10, 20, 50, 100
mean_binned_9k_w_10,std_binned_9k_w_10,std_binned_9k_w_10_up,std_binned_9k_w_10_down = mean_binned_AF(AF_sim_rep_9k_w_10_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_w_20,std_binned_9k_w_20,std_binned_9k_w_20_up,std_binned_9k_w_20_down = mean_binned_AF(AF_sim_rep_9k_w_20_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_w_50,std_binned_9k_w_50,std_binned_9k_w_50_up,std_binned_9k_w_50_down = mean_binned_AF(AF_sim_rep_9k_w_50_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_w_100,std_binned_9k_w_100,std_binned_9k_w_100_up,std_binned_9k_w_100_down = mean_binned_AF(AF_sim_rep_9k_w_100_drifted,bins=np.arange(0,1.05,0.05))

######
#plot 
######

#check how to set max value of Y axis
for i in range(14):
    print max([max(std_binned_450_w_10_up[i]),max(std_binned_450_w_20_up[i]),max(std_binned_450_w_50_up[i]),max(std_binned_450_w_50_up[i]),
               max(std_binned_9k_w_10_up[i]),max(std_binned_9k_w_20_up[i]),max(std_binned_9k_w_50_up[i]),max(std_binned_9k_w_100_up[i])])

yMax = [1.1,0.8,0.8,0.8,0.6,0.5,0.5,0.5,0.5,0.5,0.6,0.6,0.8,1]
sim_label = ['450, 10loci','450, 20loci','450, 50loci','450, 100loci','9000, 10loci','9000, 20loci','9000, 50loci','9000, 100loci']
linestyle = [':','-.','--','-']*14
timepoints = 14
gen_label = ['Gen %i'%tm for tm in np.arange(10,150,10)]
nbins = 20
bin_edge =np.array([float(i)/nbins for i in range(nbins+1)])
x = 0.5*(bin_edge[1:]+bin_edge[:-1])
xx = np.linspace(0, 1, 300)
          
fig, ax = plt.subplots(nrows=14, ncols=4, figsize=(20,25), facecolor='w', edgecolor='k')  
fig.subplots_adjust(hspace=0.15, wspace=0.1)
for ind,i in enumerate(np.arange(1,timepoints*4,4)): #10 loci
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_w_10_up[ind])
    pch_down = pchip(x,std_binned_450_w_10_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w_10[ind])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[0])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_10_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_10_drifted[ind]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_10_up[ind])
    pch_down = pchip(x,std_binned_9k_w_10_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w_10[ind])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[4])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_10_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_10_drifted[ind]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    if i in [41,45,49,53]:
        if round(np.mean(w_450_10_NumSweep[ind])) == 0 and round(np.mean(w_9k_10_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(w_450_10_NumRisingLoci[ind]))), xy=(0.55, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(w_9k_10_NumRisingLoci[ind]))), xy=(0.55, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
        if round(np.mean(w_450_10_NumSweep[ind])) != 0 or round(np.mean(w_9k_10_NumSweep[ind])) != 0:
            ax1.annotate('%i/%i' %(int(round(np.mean(w_450_10_NumRisingLoci[ind]))),int(round(np.mean(w_450_10_NumSweep[ind])))), xy=(0.55, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_10_NumRisingLoci[ind]))),int(round(np.mean(w_9k_10_NumSweep[ind])))), xy=(0.55, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    else:
        if round(np.mean(w_450_10_NumSweep[ind])) == 0 and round(np.mean(w_9k_10_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(w_450_10_NumRisingLoci[ind]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(w_9k_10_NumRisingLoci[ind]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
        if round(np.mean(w_450_10_NumSweep[ind])) != 0 or round(np.mean(w_9k_10_NumSweep[ind])) != 0:
            ax1.annotate('%i/%i' %(int(round(np.mean(w_450_10_NumRisingLoci[ind]))),int(round(np.mean(w_450_10_NumSweep[ind])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_10_NumRisingLoci[ind]))),int(round(np.mean(w_9k_10_NumSweep[ind])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if i == 1:ax1.legend(bbox_to_anchor=(0.57, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
for ind,i in enumerate(np.arange(2,timepoints*4,4)): #20 loci
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_w_20_up[ind])
    pch_down = pchip(x,std_binned_450_w_20_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w_20[ind])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[1])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_20_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_20_drifted[ind]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_20_up[ind])
    pch_down = pchip(x,std_binned_9k_w_20_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w_20[ind])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[5])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_20_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_20_drifted[ind]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    if i in [42,46,50,54]:
        if round(np.mean(w_450_20_NumSweep[ind])) == 0 and round(np.mean(w_9k_20_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(w_450_20_NumRisingLoci[ind]))), xy=(0.55, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(w_9k_20_NumRisingLoci[ind]))), xy=(0.55, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
        if round(np.mean(w_450_20_NumSweep[ind])) != 0 or round(np.mean(w_9k_20_NumSweep[ind])) != 0:
            ax1.annotate('%i/%i' %(int(round(np.mean(w_450_20_NumRisingLoci[ind]))),int(round(np.mean(w_450_20_NumSweep[ind])))), xy=(0.55, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_20_NumRisingLoci[ind]))),int(round(np.mean(w_9k_20_NumSweep[ind])))), xy=(0.55, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    else:
        if round(np.mean(w_450_20_NumSweep[ind])) == 0 and round(np.mean(w_9k_20_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(w_450_20_NumRisingLoci[ind]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(w_9k_20_NumRisingLoci[ind]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
        if round(np.mean(w_450_20_NumSweep[ind])) != 0 or round(np.mean(w_9k_20_NumSweep[ind])) != 0:
            ax1.annotate('%i/%i' %(int(round(np.mean(w_450_20_NumRisingLoci[ind]))),int(round(np.mean(w_450_20_NumSweep[ind])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_20_NumRisingLoci[ind]))),int(round(np.mean(w_9k_20_NumSweep[ind])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if i == 2:ax1.legend(bbox_to_anchor=(0.57, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
for ind,i in enumerate(np.arange(3,timepoints*4,4)): #50 loci
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_w_50_up[ind])
    pch_down = pchip(x,std_binned_450_w_50_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w_50[ind])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[2])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_50_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_50_drifted[ind]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_50_up[ind])
    pch_down = pchip(x,std_binned_9k_w_50_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w_50[ind])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[6])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_50_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_50_drifted[ind]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    if i in [43,47,51,55]:
        if round(np.mean(w_450_50_NumSweep[ind])) == 0 and round(np.mean(w_9k_50_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(w_450_50_NumRisingLoci[ind]))), xy=(0.6, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(w_9k_50_NumRisingLoci[ind]))), xy=(0.6, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
        if round(np.mean(w_450_50_NumSweep[ind])) != 0 or round(np.mean(w_9k_50_NumSweep[ind])) != 0:
            ax1.annotate('%i/%i' %(int(round(np.mean(w_450_50_NumRisingLoci[ind]))),int(round(np.mean(w_450_50_NumSweep[ind])))), xy=(0.6, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_50_NumRisingLoci[ind]))),int(round(np.mean(w_9k_50_NumSweep[ind])))), xy=(0.6, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    else:
        if round(np.mean(w_450_50_NumSweep[ind])) == 0 and round(np.mean(w_9k_50_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(w_450_50_NumRisingLoci[ind]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(w_9k_50_NumRisingLoci[ind]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
        if round(np.mean(w_450_50_NumSweep[ind])) != 0 or round(np.mean(w_9k_50_NumSweep[ind])) != 0:
            ax1.annotate('%i/%i' %(int(round(np.mean(w_450_50_NumRisingLoci[ind]))),int(round(np.mean(w_450_50_NumSweep[ind])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_50_NumRisingLoci[ind]))),int(round(np.mean(w_9k_50_NumSweep[ind])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if i == 3:ax1.legend(bbox_to_anchor=(0.57, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
for ind,i in enumerate(np.arange(4,timepoints*4+1,4)): #100 loci
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_w_100_up[ind])
    pch_down = pchip(x,std_binned_450_w_100_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w_100[ind])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[3])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_100_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_100_drifted[ind]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_100_up[ind])
    pch_down = pchip(x,std_binned_9k_w_100_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w_100[ind])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[7])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_100_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_100_drifted[ind]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    ax1.annotate(gen_label[ind], xy=(0.85, 0.8), xytext=(0.85, 0.8),textcoords='axes fraction',xycoords='axes fraction', fontweight = 'bold', size = 14)
    if i in [4,8,12,16,20]:
        if round(np.mean(w_450_100_NumSweep[ind])) == 0 and round(np.mean(w_9k_100_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(w_450_100_NumRisingLoci[ind]))), xy=(0.75, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(w_9k_100_NumRisingLoci[ind]))), xy=(0.75, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
        if round(np.mean(w_450_100_NumSweep[ind])) != 0 or round(np.mean(w_9k_100_NumSweep[ind])) != 0:
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_100_NumRisingLoci[ind]))),int(round(np.mean(w_450_100_NumSweep[ind])))), xy=(0.75, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_100_NumRisingLoci[ind]))),int(round(np.mean(w_9k_100_NumSweep[ind])))), xy=(0.75, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    else:
        if round(np.mean(w_450_100_NumSweep[ind])) == 0 and round(np.mean(w_9k_100_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(w_450_100_NumRisingLoci[ind]))), xy=(0.65, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(w_9k_100_NumRisingLoci[ind]))), xy=(0.65, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
        if round(np.mean(w_450_100_NumSweep[ind])) != 0 or round(np.mean(w_9k_100_NumSweep[ind])) != 0:
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_100_NumRisingLoci[ind]))),int(round(np.mean(w_450_100_NumSweep[ind])))), xy=(0.65, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_100_NumRisingLoci[ind]))),int(round(np.mean(w_9k_100_NumSweep[ind])))), xy=(0.65, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if i == 4:ax1.legend(bbox_to_anchor=(0.55, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
fig.text(0.08, 0.5, 'Proportion', va='center', rotation='vertical',size = 24) #ylabel Count Probability
fig.text(0.5, 0.1, 'Allele frequency', ha='center',size = 24) #xlabel
plt.savefig('sweep_diffNumLoci_AF_RisingLociNum_smoothCurve.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_diffNumLoci_AF_RisingLociNum_smoothCurve.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''
#################
#Figure 6 and S7
#################
'''

################
#get drift data
###############

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

##################
#get simulation AF 
##################

#qff 450, 10, 20, 50, 100 loci
sim = 'sel_0.5-4.5-0.74-0.3_eff0.04_10loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_450_qtl_10 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_450_qtl_10, AF_sim_rep_450_qtl_10 = process_AF_SimBatch(AF_sim_450_qtl_10,10,500)

sim = 'sel_0.5-4.5-0.38-0.3_eff0.04_20loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_450_qtl_20 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_450_qtl_20, AF_sim_rep_450_qtl_20 = process_AF_SimBatch(AF_sim_450_qtl_20,20,500)

sim = 'sel_0.5-4.5--0.7-0.3_eff0.04_50loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_450_qtl_50 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_450_qtl_50, AF_sim_rep_450_qtl_50 = process_AF_SimBatch(AF_sim_450_qtl_50,50,500)

sim = 'sel_0.5-4.5--2.5-0.3'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/recom_rate_noX/1000sims/AF.sync.gz' %sim,'rb')  
AF_sim_450_qtl_100 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_450_qtl_100, AF_sim_rep_450_qtl_100 = process_AF_SimBatch(AF_sim_450_qtl_100,100,500)

#qff 9000, 10, 20, 50, 100 loci, effect size: 0.04 moving TO
sim = 'Ne9000_sel_0.5-4.5-0.74-0.3_eff0.04_10loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_9k_qtl_10 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_9k_qtl_10, AF_sim_rep_9k_qtl_10 = process_AF_SimBatch(AF_sim_9k_qtl_10,10,500)

sim = 'Ne9000_sel_0.5-4.5-0.38-0.3_eff0.04_20loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_9k_qtl_20 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_9k_qtl_20, AF_sim_rep_9k_qtl_20 = process_AF_SimBatch(AF_sim_9k_qtl_20,20,500)

sim = 'Ne9000_sel_0.5-4.5--0.7-0.3_eff0.04_50loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_9k_qtl_50 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_9k_qtl_50, AF_sim_rep_9k_qtl_50 = process_AF_SimBatch(AF_sim_9k_qtl_50,50,500)

sim = 'Ne9000_sel_0.5-4.5--2.5-0.3'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/recom_rate_noX/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_9k_qtl_100 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_9k_qtl_100, AF_sim_rep_9k_qtl_100 = process_AF_SimBatch(AF_sim_9k_qtl_100,100,500)

################################
#correct for drift
#get # of rising and sweep loci
################################

#QTL 450, loci 10, 20, 50, 100
AF_sim_rep_450_qtl_10_drifted, qtl_450_10_NumRisingLoci, qtl_450_10_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_qtl_10,AFC_thresh_neu_450)
AF_sim_rep_450_qtl_20_drifted, qtl_450_20_NumRisingLoci, qtl_450_20_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_qtl_20,AFC_thresh_neu_450)
AF_sim_rep_450_qtl_50_drifted, qtl_450_50_NumRisingLoci, qtl_450_50_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_qtl_50,AFC_thresh_neu_450)
AF_sim_rep_450_qtl_100_drifted, qtl_450_100_NumRisingLoci, qtl_450_100_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_qtl_100,AFC_thresh_neu_450)

#QTL 9k, loci 10, 20, 50, 100
AF_sim_rep_9k_qtl_10_drifted, qtl_9k_10_NumRisingLoci, qtl_9k_10_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_qtl_10,AFC_thresh_neu_9k)
AF_sim_rep_9k_qtl_20_drifted, qtl_9k_20_NumRisingLoci, qtl_9k_20_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_qtl_20,AFC_thresh_neu_9k)
AF_sim_rep_9k_qtl_50_drifted, qtl_9k_50_NumRisingLoci, qtl_9k_50_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_qtl_50,AFC_thresh_neu_9k)
AF_sim_rep_9k_qtl_100_drifted, qtl_9k_100_NumRisingLoci, qtl_9k_100_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_qtl_100,AFC_thresh_neu_9k)

#compute median AF
#QTL 450, loci 10, 20, 50, 100
median_AF_sim_rep_450_qtl_10_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_qtl_10_drifted[1:]]   
median_AF_sim_rep_450_qtl_20_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_qtl_20_drifted[1:]]   
median_AF_sim_rep_450_qtl_50_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_qtl_50_drifted[1:]]   
median_AF_sim_rep_450_qtl_100_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_qtl_100_drifted[1:]]   

#QTL 9k, loci 10, 20, 50, 100
median_AF_sim_rep_9k_qtl_10_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_qtl_10_drifted[1:]]   
median_AF_sim_rep_9k_qtl_20_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_qtl_20_drifted[1:]]   
median_AF_sim_rep_9k_qtl_50_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_qtl_50_drifted[1:]]   
median_AF_sim_rep_9k_qtl_100_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_qtl_100_drifted[1:]]   

######################################################
#get binned mean and standard deviation of replicates
######################################################

#QTL 450, loci 10, 20, 50, 100
mean_binned_450_qtl_10,std_binned_450_qtl_10,std_binned_450_qtl_10_up,std_binned_450_qtl_10_down = mean_binned_AF(AF_sim_rep_450_qtl_10_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_450_qtl_20,std_binned_450_qtl_20,std_binned_450_qtl_20_up,std_binned_450_qtl_20_down = mean_binned_AF(AF_sim_rep_450_qtl_20_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_450_qtl_50,std_binned_450_qtl_50,std_binned_450_qtl_50_up,std_binned_450_qtl_50_down = mean_binned_AF(AF_sim_rep_450_qtl_50_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_450_qtl_100,std_binned_450_qtl_100,std_binned_450_qtl_100_up,std_binned_450_qtl_100_down = mean_binned_AF(AF_sim_rep_450_qtl_100_drifted,bins=np.arange(0,1.05,0.05))

#QTL 9k, loci 10, 20, 50, 100
mean_binned_9k_qtl_10,std_binned_9k_qtl_10,std_binned_9k_qtl_10_up,std_binned_9k_qtl_10_down = mean_binned_AF(AF_sim_rep_9k_qtl_10_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_qtl_20,std_binned_9k_qtl_20,std_binned_9k_qtl_20_up,std_binned_9k_qtl_20_down = mean_binned_AF(AF_sim_rep_9k_qtl_20_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_qtl_50,std_binned_9k_qtl_50,std_binned_9k_qtl_50_up,std_binned_9k_qtl_50_down = mean_binned_AF(AF_sim_rep_9k_qtl_50_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_qtl_100,std_binned_9k_qtl_100,std_binned_9k_qtl_100_up,std_binned_9k_qtl_100_down = mean_binned_AF(AF_sim_rep_9k_qtl_100_drifted,bins=np.arange(0,1.05,0.05))

###############
#plot Figure 6
###############

#check how to set max value of Y axis
for i in range(14):
    print max([max(std_binned_450_qtl_10_up[i]),max(std_binned_450_qtl_20_up[i]),max(std_binned_450_qtl_50_up[i]),max(std_binned_450_qtl_50_up[i]),
               max(std_binned_9k_qtl_10_up[i]),max(std_binned_9k_qtl_20_up[i]),max(std_binned_9k_qtl_50_up[i]),max(std_binned_9k_qtl_100_up[i])])

yMax = [1.1,1.1,1.1,1,1,0.8,0.8,0.8,0.8,0.8,0.5,0.5,0.5,0.5]
timepoints = 14
gen_label = ['Gen %i'%tm for tm in np.arange(10,150,10)]
sim_label = ['450, 10loci','450, 20loci','450, 50loci','450, 100loci','9000, 10loci','9000, 20loci','9000, 50loci','9000, 100loci']
linestyle = [':','-.','--','-']*14
nbins = 20
bin_edge =np.array([float(i)/nbins for i in range(nbins+1)])
x = 0.5*(bin_edge[1:]+bin_edge[:-1])
xx = np.linspace(0, 1, 300)
          
fig, ax = plt.subplots(nrows=14, ncols=4, figsize=(20,25), facecolor='w', edgecolor='k')  
fig.subplots_adjust(hspace=0.15, wspace=0.1)
for ind,i in enumerate(np.arange(1,timepoints*4,4)): #10 loci
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_qtl_10_up[ind])
    pch_down = pchip(x,std_binned_450_qtl_10_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl_10[ind])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[0])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_10_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_10_drifted[ind]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_10_up[ind])
    pch_down = pchip(x,std_binned_9k_qtl_10_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl_10[ind])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[4])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_10_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_10_drifted[ind]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    if round(np.mean(qtl_450_10_NumSweep[ind])) == 0 and round(np.mean(qtl_9k_10_NumSweep[ind])) == 0:
        ax1.annotate('%i' %int(round(np.mean(qtl_450_10_NumRisingLoci[ind]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(qtl_9k_10_NumRisingLoci[ind]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_10_NumSweep[ind])) != 0 or round(np.mean(qtl_9k_10_NumSweep[ind])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_10_NumRisingLoci[ind]))),int(round(np.mean(qtl_450_10_NumSweep[ind])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_10_NumRisingLoci[ind]))),int(round(np.mean(qtl_9k_10_NumSweep[ind])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if i == 1:ax1.legend(bbox_to_anchor=(0.57, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
for ind,i in enumerate(np.arange(2,timepoints*4,4)): #20 loci
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_qtl_20_up[ind])
    pch_down = pchip(x,std_binned_450_qtl_20_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl_20[ind])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[1])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_20_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_20_drifted[ind]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_20_up[ind])
    pch_down = pchip(x,std_binned_9k_qtl_20_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl_20[ind])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[5])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_20_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_20_drifted[ind]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    if i in [42,46,50,54]:
        if round(np.mean(qtl_450_20_NumSweep[ind])) == 0 and round(np.mean(qtl_9k_20_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(qtl_450_20_NumRisingLoci[ind]))), xy=(0.55, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(qtl_9k_20_NumRisingLoci[ind]))), xy=(0.55, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
        if round(np.mean(qtl_450_20_NumSweep[ind])) != 0 or round(np.mean(qtl_9k_20_NumSweep[ind])) != 0:
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_20_NumRisingLoci[ind]))),int(round(np.mean(qtl_450_20_NumSweep[ind])))), xy=(0.55, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_20_NumRisingLoci[ind]))),int(round(np.mean(qtl_9k_20_NumSweep[ind])))), xy=(0.55, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    else:
        if round(np.mean(qtl_450_20_NumSweep[ind])) == 0 and round(np.mean(qtl_9k_20_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(qtl_450_20_NumRisingLoci[ind]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(qtl_9k_20_NumRisingLoci[ind]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
        if round(np.mean(qtl_450_20_NumSweep[ind])) != 0 or round(np.mean(qtl_9k_20_NumSweep[ind])) != 0:
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_20_NumRisingLoci[ind]))),int(round(np.mean(qtl_450_20_NumSweep[ind])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_20_NumRisingLoci[ind]))),int(round(np.mean(qtl_9k_20_NumSweep[ind])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if i == 2:ax1.legend(bbox_to_anchor=(0.57, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
for ind,i in enumerate(np.arange(3,timepoints*4,4)): #50 loci
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_qtl_50_up[ind])
    pch_down = pchip(x,std_binned_450_qtl_50_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl_50[ind])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[2])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_50_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_50_drifted[ind]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_50_up[ind])
    pch_down = pchip(x,std_binned_9k_qtl_50_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl_50[ind])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[6])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_50_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_50_drifted[ind]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    if round(np.mean(qtl_450_50_NumSweep[ind])) == 0 and round(np.mean(qtl_9k_50_NumSweep[ind])) == 0:
        ax1.annotate('%i' %int(round(np.mean(qtl_450_50_NumRisingLoci[ind]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(qtl_9k_50_NumRisingLoci[ind]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_50_NumSweep[ind])) != 0 or round(np.mean(qtl_9k_50_NumSweep[ind])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_50_NumRisingLoci[ind]))),int(round(np.mean(qtl_450_50_NumSweep[ind])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_50_NumRisingLoci[ind]))),int(round(np.mean(qtl_9k_50_NumSweep[ind])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if i == 3:ax1.legend(bbox_to_anchor=(0.57, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
for ind,i in enumerate(np.arange(4,timepoints*4+1,4)): #100 loci
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_qtl_100_up[ind])
    pch_down = pchip(x,std_binned_450_qtl_100_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl_100[ind])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[3])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_100_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_100_drifted[ind]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_100_up[ind])
    pch_down = pchip(x,std_binned_9k_qtl_100_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl_100[ind])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[7])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_100_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_100_drifted[ind]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    ax1.annotate(gen_label[ind], xy=(0.85, 0.8), xytext=(0.85, 0.8),textcoords='axes fraction',xycoords='axes fraction', fontweight = 'bold', size = 14)
    if i in [4,8,12,16,20]:
        if round(np.mean(qtl_450_100_NumSweep[ind])) == 0 and round(np.mean(qtl_9k_100_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(qtl_450_100_NumRisingLoci[ind]))), xy=(0.75, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(qtl_9k_100_NumRisingLoci[ind]))), xy=(0.75, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
        if round(np.mean(qtl_450_100_NumSweep[ind])) != 0 or round(np.mean(qtl_9k_100_NumSweep[ind])) != 0:
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_100_NumRisingLoci[ind]))),int(round(np.mean(qtl_450_100_NumSweep[ind])))), xy=(0.75, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_100_NumRisingLoci[ind]))),int(round(np.mean(qtl_9k_100_NumSweep[ind])))), xy=(0.75, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    else:
        if round(np.mean(qtl_450_100_NumSweep[ind])) == 0 and round(np.mean(qtl_9k_100_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(qtl_450_100_NumRisingLoci[ind]))), xy=(0.65, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(qtl_9k_100_NumRisingLoci[ind]))), xy=(0.65, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
        if round(np.mean(qtl_450_100_NumSweep[ind])) != 0 or round(np.mean(qtl_9k_100_NumSweep[ind])) != 0:
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_100_NumRisingLoci[ind]))),int(round(np.mean(qtl_450_100_NumSweep[ind])))), xy=(0.65, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_100_NumRisingLoci[ind]))),int(round(np.mean(qtl_9k_100_NumSweep[ind])))), xy=(0.65, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if i == 4:ax1.legend(bbox_to_anchor=(0.55, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
fig.text(0.08, 0.5, 'Proportion', va='center', rotation='vertical',size = 24) #ylabel Count Probability
fig.text(0.5, 0.1, 'Allele frequency', ha='center',size = 24) #xlabel
plt.savefig('QTL_diffNumLoci_AF_RisingLociNum_smoothCurve.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('QTL_diffNumLoci_AF_RisingLociNum_smoothCurve.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

###############
#plot Figure S7
###############

#plot selected generations F10, F50, F100, F140
ind_toplot = [0,4,9,13]

#check how to set max value of Y axis
for i in ind_toplot:
    print max([max(std_binned_450_qtl_10_up[i]),max(std_binned_450_qtl_20_up[i]),max(std_binned_450_qtl_50_up[i]),max(std_binned_450_qtl_50_up[i]),
               max(std_binned_9k_qtl_10_up[i]),max(std_binned_9k_qtl_20_up[i]),max(std_binned_9k_qtl_50_up[i]),max(std_binned_9k_qtl_100_up[i])])

yMax = [1.1,1,0.8,0.5]
timepoints = len(ind_toplot)
gen_label = ['Gen %i'%tm for tm in np.arange(10,150,10)]
sim_label = ['450, 10loci','450, 20loci','450, 50loci','450, 100loci','9000, 10loci','9000, 20loci','9000, 50loci','9000, 100loci']
linestyle = [':','-.','--','-']*14
nbins = 20
bin_edge =np.array([float(i)/nbins for i in range(nbins+1)])
x = 0.5*(bin_edge[1:]+bin_edge[:-1])
xx = np.linspace(0, 1, 300)
          
fig, ax = plt.subplots(nrows=4, ncols=4, figsize=(17,10), facecolor='w', edgecolor='k')  
fig.subplots_adjust(hspace=0.15, wspace=0.1)
for ind,i in enumerate(np.arange(1,timepoints*4,4)): #10 loci
    ax1=plt.subplot(4,4,i)
    pch_up = pchip(x,std_binned_450_qtl_10_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_450_qtl_10_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl_10[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[0])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_10_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_10_drifted[ind_toplot[ind]]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_10_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_9k_qtl_10_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl_10[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[4])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_10_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_10_drifted[ind_toplot[ind]]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes_sub()
    if round(np.mean(qtl_450_10_NumSweep[ind_toplot[ind]])) == 0 and round(np.mean(qtl_9k_10_NumSweep[ind_toplot[ind]])) == 0:
        ax1.annotate('%i' %int(round(np.mean(qtl_450_10_NumRisingLoci[ind_toplot[ind]]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(qtl_9k_10_NumRisingLoci[ind_toplot[ind]]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_10_NumSweep[ind_toplot[ind]])) != 0 or round(np.mean(qtl_9k_10_NumSweep[ind_toplot[ind]])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_10_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(qtl_450_10_NumSweep[ind_toplot[ind]])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_10_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(qtl_9k_10_NumSweep[ind_toplot[ind]])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if i == 1:ax1.legend(bbox_to_anchor=(0.5, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
for ind,i in enumerate(np.arange(2,timepoints*4,4)): #20 loci
    ax1=plt.subplot(4,4,i)
    pch_up = pchip(x,std_binned_450_qtl_20_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_450_qtl_20_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl_20[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[1])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_20_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_20_drifted[ind_toplot[ind]]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_20_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_9k_qtl_20_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl_20[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[5])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_20_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_20_drifted[ind_toplot[ind]]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes_sub()
    if i in [10,14]:
        if round(np.mean(qtl_450_20_NumSweep[ind_toplot[ind]])) == 0 and round(np.mean(qtl_9k_20_NumSweep[ind_toplot[ind]])) == 0:
            ax1.annotate('%i' %int(round(np.mean(qtl_450_20_NumRisingLoci[ind_toplot[ind]]))), xy=(0.55, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(qtl_9k_20_NumRisingLoci[ind_toplot[ind]]))), xy=(0.55, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
        if round(np.mean(qtl_450_20_NumSweep[ind_toplot[ind]])) != 0 or round(np.mean(qtl_9k_20_NumSweep[ind_toplot[ind]])) != 0:
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_20_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(qtl_450_20_NumSweep[ind_toplot[ind]])))), xy=(0.55, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_20_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(qtl_9k_20_NumSweep[ind_toplot[ind]])))), xy=(0.55, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    else:
        if round(np.mean(qtl_450_20_NumSweep[ind_toplot[ind]])) == 0 and round(np.mean(qtl_9k_20_NumSweep[ind_toplot[ind]])) == 0:
            ax1.annotate('%i' %int(round(np.mean(qtl_450_20_NumRisingLoci[ind_toplot[ind]]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(qtl_9k_20_NumRisingLoci[ind_toplot[ind]]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
        if round(np.mean(qtl_450_20_NumSweep[ind_toplot[ind]])) != 0 or round(np.mean(qtl_9k_20_NumSweep[ind_toplot[ind]])) != 0:
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_20_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(qtl_450_20_NumSweep[ind_toplot[ind]])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_20_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(qtl_9k_20_NumSweep[ind_toplot[ind]])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if i == 2:ax1.legend(bbox_to_anchor=(0.5, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
for ind,i in enumerate(np.arange(3,timepoints*4,4)): #50 loci
    ax1=plt.subplot(4,4,i)
    pch_up = pchip(x,std_binned_450_qtl_50_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_450_qtl_50_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl_50[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[2])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_50_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_50_drifted[ind_toplot[ind]]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_50_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_9k_qtl_50_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl_50[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[6])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_50_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_50_drifted[ind_toplot[ind]]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes_sub()
    if round(np.mean(qtl_450_50_NumSweep[ind_toplot[ind]])) == 0 and round(np.mean(qtl_9k_50_NumSweep[ind_toplot[ind]])) == 0:
        ax1.annotate('%i' %int(round(np.mean(qtl_450_50_NumRisingLoci[ind_toplot[ind]]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(qtl_9k_50_NumRisingLoci[ind_toplot[ind]]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_50_NumSweep[ind_toplot[ind]])) != 0 or round(np.mean(qtl_9k_50_NumSweep[ind_toplot[ind]])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_50_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(qtl_450_50_NumSweep[ind_toplot[ind]])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_50_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(qtl_9k_50_NumSweep[ind_toplot[ind]])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if i == 3:ax1.legend(bbox_to_anchor=(0.5, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
for ind,i in enumerate(np.arange(4,timepoints*4+1,4)): #100 loci
    ax1=plt.subplot(4,4,i)
    pch_up = pchip(x,std_binned_450_qtl_100_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_450_qtl_100_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl_100[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[3])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_100_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_100_drifted[ind_toplot[ind]]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_100_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_9k_qtl_100_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl_100[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[7])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_100_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_100_drifted[ind_toplot[ind]]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes_sub()
    ax1.annotate(gen_label[ind_toplot[ind]], xy=(0.85, 0.8), xytext=(0.85, 0.8),textcoords='axes fraction',xycoords='axes fraction', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_100_NumSweep[ind_toplot[ind]])) == 0 and round(np.mean(qtl_9k_100_NumSweep[ind_toplot[ind]])) == 0:
        ax1.annotate('%i' %int(round(np.mean(qtl_450_100_NumRisingLoci[ind_toplot[ind]]))), xy=(0.65, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(qtl_9k_100_NumRisingLoci[ind_toplot[ind]]))), xy=(0.65, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_100_NumSweep[ind_toplot[ind]])) != 0 or round(np.mean(qtl_9k_100_NumSweep[ind_toplot[ind]])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_100_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(qtl_450_100_NumSweep[ind_toplot[ind]])))), xy=(0.65, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_100_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(qtl_9k_100_NumSweep[ind_toplot[ind]])))), xy=(0.65, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if i == 4:ax1.legend(bbox_to_anchor=(0.55, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
fig.text(0.08, 0.5, 'Proportion', va='center', rotation='vertical',size = 24) #ylabel Count Probability
fig.text(0.5, 0.05, 'Allele frequency', ha='center',size = 24) #xlabel
plt.savefig('QTL_diffNumLoci_AF_RisingLociNum_smoothCurve_sub.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('QTL_diffNumLoci_AF_RisingLociNum_smoothCurve_sub.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''
##################
#Figure 8 and S10
#################
'''
################
#get drift data
###############

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

##################
#get simulation AF 
##################

#sweep 450 0.02,0.05,0.08,0.1 s
sim = 'Ne450_100loci_0.05p0_0.02s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_450_w_02 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_450_w_02, AF_sim_rep_450_w_02 = process_AF_SimBatch(AF_sim_450_w_02,100,500)

sim =  'Ne450_100loci_0.05p0_0.05s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_450_w_05 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_450_w_05, AF_sim_rep_450_w_05 = process_AF_SimBatch(AF_sim_450_w_05,100,500)

sim = 'Ne450_100loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/1000sims/AF_results.sync' %sim,'rb')  
AF_sim_450_w_08 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_450_w_08, AF_sim_rep_450_w_08 = process_AF_SimBatch(AF_sim_450_w_08,100,500)

sim =  'Ne450_100loci_0.05p0_0.1s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_450_w_1 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_450_w_1, AF_sim_rep_450_w_1 = process_AF_SimBatch(AF_sim_450_w_1,100,500)

#sweep 9000, 0.02,0.05,0.08,0.1 s
sim = 'Ne9000_100loci_0.05p0_0.02s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_9k_w_02 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_9k_w_02, AF_sim_rep_9k_w_02 = process_AF_SimBatch(AF_sim_9k_w_02,100,500)

sim =  'Ne9000_100loci_0.05p0_0.05s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_9k_w_05 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_9k_w_05, AF_sim_rep_9k_w_05 = process_AF_SimBatch(AF_sim_9k_w_05,100,500)

sim = 'Ne9000_100loci_0.05p0_0.08s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_9k_w_08 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_9k_w_08, AF_sim_rep_9k_w_08 = process_AF_SimBatch(AF_sim_9k_w_08,100,500)

sim =  'Ne9000_100loci_0.05p0_0.1s'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/sel.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/%s/recom_rate_noX/500sims/AF_results.sync' %sim,'rb')  
AF_sim_9k_w_1 = get_AF_w(InputFile_s,InputFile_af)
AF_sim_class_9k_w_1, AF_sim_rep_9k_w_1 = process_AF_SimBatch(AF_sim_9k_w_1,100,500)

################################
#correct for drift
#get # of rising and sweep loci
################################

#sweep 450, 0.02,0.05,0.08,0.1 s
AF_sim_rep_450_w_02_drifted, w_450_02_NumRisingLoci, w_450_02_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_w_02,AFC_thresh_neu_450)
AF_sim_rep_450_w_05_drifted, w_450_05_NumRisingLoci, w_450_05_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_w_05,AFC_thresh_neu_450)
AF_sim_rep_450_w_08_drifted, w_450_08_NumRisingLoci, w_450_08_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_w_08,AFC_thresh_neu_450)
AF_sim_rep_450_w_1_drifted, w_450_1_NumRisingLoci, w_450_1_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_w_1,AFC_thresh_neu_450)

#sweep 9k, 0.02,0.05,0.08,0.1 s
AF_sim_rep_9k_w_02_drifted, w_9k_02_NumRisingLoci, w_9k_02_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_w_02,AFC_thresh_neu_9k)
AF_sim_rep_9k_w_05_drifted, w_9k_05_NumRisingLoci, w_9k_05_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_w_05,AFC_thresh_neu_9k)
AF_sim_rep_9k_w_08_drifted, w_9k_08_NumRisingLoci, w_9k_08_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_w_08,AFC_thresh_neu_9k)
AF_sim_rep_9k_w_1_drifted, w_9k_1_NumRisingLoci, w_9k_1_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_w_1,AFC_thresh_neu_9k)

#compute median AF
#sweep 450, 0.02,0.05,0.08,0.1 s
median_AF_sim_rep_450_w_02_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_w_02_drifted[1:]]   
median_AF_sim_rep_450_w_05_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_w_05_drifted[1:]]   
median_AF_sim_rep_450_w_08_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_w_08_drifted[1:]]   
median_AF_sim_rep_450_w_1_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_w_1_drifted[1:]]   

#sweep 9k, 0.02,0.05,0.08,0.1 s
median_AF_sim_rep_9k_w_02_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_w_02_drifted[1:]]   
median_AF_sim_rep_9k_w_05_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_w_05_drifted[1:]]   
median_AF_sim_rep_9k_w_08_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_w_08_drifted[1:]]   
median_AF_sim_rep_9k_w_1_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_w_1_drifted[1:]]   

######################################################
#get binned mean and standard deviation of replicates
######################################################

#sweep 450, 0.02,0.05,0.08,0.1 s
mean_binned_450_w_02,std_binned_450_w_02,std_binned_450_w_02_up,std_binned_450_w_02_down = mean_binned_AF(AF_sim_rep_450_w_02_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_450_w_05,std_binned_450_w_05,std_binned_450_w_05_up,std_binned_450_w_05_down = mean_binned_AF(AF_sim_rep_450_w_05_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_450_w_08,std_binned_450_w_08,std_binned_450_w_08_up,std_binned_450_w_08_down = mean_binned_AF(AF_sim_rep_450_w_08_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_450_w_1,std_binned_450_w_1,std_binned_450_w_1_up,std_binned_450_w_1_down = mean_binned_AF(AF_sim_rep_450_w_1_drifted,bins=np.arange(0,1.05,0.05))

#sweep 9k, 0.02,0.05,0.08,0.1 s
mean_binned_9k_w_02,std_binned_9k_w_02,std_binned_9k_w_02_up,std_binned_9k_w_02_down = mean_binned_AF(AF_sim_rep_9k_w_02_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_w_05,std_binned_9k_w_05,std_binned_9k_w_05_up,std_binned_9k_w_05_down = mean_binned_AF(AF_sim_rep_9k_w_05_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_w_08,std_binned_9k_w_08,std_binned_9k_w_08_up,std_binned_9k_w_08_down = mean_binned_AF(AF_sim_rep_9k_w_08_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_w_1,std_binned_9k_w_1,std_binned_9k_w_1_up,std_binned_9k_w_1_down = mean_binned_AF(AF_sim_rep_9k_w_1_drifted,bins=np.arange(0,1.05,0.05))

#################
#plot Figure S10
#################

#check how to set max value of Y axis
for i in range(14):
    print max([max(std_binned_450_w_02_up[i]),max(std_binned_450_w_05_up[i]),max(std_binned_450_w_08_up[i]),max(std_binned_450_w_08_up[i]),
               max(std_binned_9k_w_02_up[i]),max(std_binned_9k_w_05_up[i]),max(std_binned_9k_w_08_up[i]),max(std_binned_9k_w_1_up[i])])

yMax = [1.1,1.1,1.1,1,1,0.8,0.6,0.8,0.8,0.8,0.6,0.6,0.6,0.8]
timepoints = 14
gen_label = ['Gen %i'%tm for tm in np.arange(10,150,10)]
sim_label = ['450, s=0.02','450, s=0.05','450, s=0.08','450, s=0.1','9000, s=0.02','9000, s=0.05','9000, s=0.08','9000, s=0.1']
linestyle = [':','-.','-','--']*14  
nbins = 20
bin_edge =np.array([float(i)/nbins for i in range(nbins+1)])
x = 0.5*(bin_edge[1:]+bin_edge[:-1])
xx = np.linspace(0, 1, 300)
          
fig, ax = plt.subplots(nrows=14, ncols=4, figsize=(20,25), facecolor='w', edgecolor='k')  
fig.subplots_adjust(hspace=0.15, wspace=0.1)
for ind,i in enumerate(np.arange(1,timepoints*4,4)): #s 0.02
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_w_02_up[ind])
    pch_down = pchip(x,std_binned_450_w_02_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w_02[ind])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[0])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_02_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_02_drifted[ind]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_02_up[ind])
    pch_down = pchip(x,std_binned_9k_w_02_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w_02[ind])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[4])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_02_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_02_drifted[ind]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    if round(np.mean(w_450_02_NumSweep[ind])) == 0 and round(np.mean(w_9k_02_NumSweep[ind])) == 0:
        ax1.annotate('%i' %int(round(np.mean(w_450_02_NumRisingLoci[ind]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(w_9k_02_NumRisingLoci[ind]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if round(np.mean(w_450_02_NumSweep[ind])) != 0 or round(np.mean(w_9k_02_NumSweep[ind])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_02_NumRisingLoci[ind]))),int(round(np.mean(w_450_02_NumSweep[ind])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_02_NumRisingLoci[ind]))),int(round(np.mean(w_9k_02_NumSweep[ind])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if i == 1:ax1.legend(bbox_to_anchor=(0.57, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
for ind,i in enumerate(np.arange(2,timepoints*4,4)): #s 0.05
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_w_05_up[ind])
    pch_down = pchip(x,std_binned_450_w_05_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w_05[ind])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[1])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_05_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_05_drifted[ind]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_05_up[ind])
    pch_down = pchip(x,std_binned_9k_w_05_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w_05[ind])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[5])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_05_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_05_drifted[ind]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    if round(np.mean(w_450_05_NumSweep[ind])) == 0 and round(np.mean(w_9k_05_NumSweep[ind])) == 0:
        ax1.annotate('%i' %int(round(np.mean(w_450_05_NumRisingLoci[ind]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(w_9k_05_NumRisingLoci[ind]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if round(np.mean(w_450_05_NumSweep[ind])) != 0 or round(np.mean(w_9k_05_NumSweep[ind])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_05_NumRisingLoci[ind]))),int(round(np.mean(w_450_05_NumSweep[ind])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_05_NumRisingLoci[ind]))),int(round(np.mean(w_9k_05_NumSweep[ind])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if i == 2:ax1.legend(bbox_to_anchor=(0.57, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
for ind,i in enumerate(np.arange(3,timepoints*4,4)): #s 0.08
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_w_08_up[ind])
    pch_down = pchip(x,std_binned_450_w_08_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w_08[ind])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[2])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_08_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_08_drifted[ind]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_08_up[ind])
    pch_down = pchip(x,std_binned_9k_w_08_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w_08[ind])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[6])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_08_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_08_drifted[ind]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    if i in [43,47,51,55]:
        if round(np.mean(w_450_08_NumSweep[ind])) == 0 and round(np.mean(w_9k_08_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(w_450_08_NumRisingLoci[ind]))), xy=(0.6, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(w_9k_08_NumRisingLoci[ind]))), xy=(0.6, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
        if round(np.mean(w_450_08_NumSweep[ind])) != 0 or round(np.mean(w_9k_08_NumSweep[ind])) != 0:
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_08_NumRisingLoci[ind]))),int(round(np.mean(w_450_08_NumSweep[ind])))), xy=(0.6, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_08_NumRisingLoci[ind]))),int(round(np.mean(w_9k_08_NumSweep[ind])))), xy=(0.6, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    else:
        if round(np.mean(w_450_08_NumSweep[ind])) == 0 and round(np.mean(w_9k_08_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(w_450_08_NumRisingLoci[ind]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(w_9k_08_NumRisingLoci[ind]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
        if round(np.mean(w_450_08_NumSweep[ind])) != 0 or round(np.mean(w_9k_08_NumSweep[ind])) != 0:
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_08_NumRisingLoci[ind]))),int(round(np.mean(w_450_08_NumSweep[ind])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_08_NumRisingLoci[ind]))),int(round(np.mean(w_9k_08_NumSweep[ind])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if i == 3:ax1.legend(bbox_to_anchor=(0.57, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
for ind,i in enumerate(np.arange(4,timepoints*4+1,4)): #s 0.1
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_w_1_up[ind])
    pch_down = pchip(x,std_binned_450_w_1_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w_1[ind])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[3])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_1_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_1_drifted[ind]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_1_up[ind])
    pch_down = pchip(x,std_binned_9k_w_1_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w_1[ind])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[7])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_1_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_1_drifted[ind]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    ax1.annotate(gen_label[ind], xy=(0.7, 0.8), xytext=(0.7, 0.8),textcoords='axes fraction',xycoords='axes fraction', fontweight = 'bold', size = 14)
    if i in [4,8,12,16,20]:
        if round(np.mean(w_450_1_NumSweep[ind])) == 0 and round(np.mean(w_9k_1_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(w_450_1_NumRisingLoci[ind]))), xy=(0.6, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(w_9k_1_NumRisingLoci[ind]))), xy=(0.6, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
        if round(np.mean(w_450_1_NumSweep[ind])) != 0 or round(np.mean(w_9k_1_NumSweep[ind])) != 0:
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_1_NumRisingLoci[ind]))),int(round(np.mean(w_450_1_NumSweep[ind])))), xy=(0.6, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_1_NumRisingLoci[ind]))),int(round(np.mean(w_9k_1_NumSweep[ind])))), xy=(0.6, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    else:
        if round(np.mean(w_450_1_NumSweep[ind])) == 0 and round(np.mean(w_9k_1_NumSweep[ind])) == 0:
            ax1.annotate('%i' %int(round(np.mean(w_450_1_NumRisingLoci[ind]))), xy=(0.6, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i' %int(round(np.mean(w_9k_1_NumRisingLoci[ind]))), xy=(0.6, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
        if round(np.mean(w_450_1_NumSweep[ind])) != 0 or round(np.mean(w_9k_1_NumSweep[ind])) != 0:
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_1_NumRisingLoci[ind]))),int(round(np.mean(w_450_1_NumSweep[ind])))), xy=(0.6, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
            ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_1_NumRisingLoci[ind]))),int(round(np.mean(w_9k_1_NumSweep[ind])))), xy=(0.6, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if i == 4:ax1.legend(bbox_to_anchor=(0.55, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
    if yMax[ind] == 0.5: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.1))
fig.text(0.08, 0.5, 'Proportion', va='center', rotation='vertical',size = 24) #ylabel Count Probability
fig.text(0.5, 0.1, 'Allele frequency', ha='center',size = 24) #xlabel
plt.savefig('sweep_diffS_AF_RisingLociNum_smoothCurve.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_diffS_AF_RisingLociNum_smoothCurve.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

###############
#plot Figure 8
###############

#plot selected generations F10, F50, F100, F140
ind_toplot = [0,4,9,13]

#check how to set max value of Y axis
for i in ind_toplot:
    print max([max(std_binned_450_w_02_up[i]),max(std_binned_450_w_05_up[i]),max(std_binned_450_w_08_up[i]),max(std_binned_450_w_08_up[i]),
               max(std_binned_9k_w_02_up[i]),max(std_binned_9k_w_05_up[i]),max(std_binned_9k_w_08_up[i]),max(std_binned_9k_w_1_up[i])])

yMax = yMax = [1.1,1,0.8,0.8]
timepoints = len(ind_toplot)
gen_label = ['Gen %i'%tm for tm in np.arange(10,150,10)]
sim_label = ['450, s=0.02','450, s=0.05','450, s=0.08','450, s=0.1','9000, s=0.02','9000, s=0.05','9000, s=0.08','9000, s=0.1']
linestyle = [':','-.','-','--']*14  
nbins = 20
bin_edge =np.array([float(i)/nbins for i in range(nbins+1)])
x = 0.5*(bin_edge[1:]+bin_edge[:-1])
xx = np.linspace(0, 1, 300)
          
fig, ax = plt.subplots(nrows=4, ncols=4, figsize=(17,10), facecolor='w', edgecolor='k')  
fig.subplots_adjust(hspace=0.15, wspace=0.1)
for ind,i in enumerate(np.arange(1,timepoints*4,4)): #s 0.02
    ax1=plt.subplot(4,4,i)
    pch_up = pchip(x,std_binned_450_w_02_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_450_w_02_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w_02[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[0])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_02_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_02_drifted[ind_toplot[ind]]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_02_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_9k_w_02_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w_02[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[4])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_02_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_02_drifted[ind_toplot[ind]]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes_sub()
    if round(np.mean(w_450_02_NumSweep[ind_toplot[ind]])) == 0 and round(np.mean(w_9k_02_NumSweep[ind_toplot[ind]])) == 0:
        ax1.annotate('%i' %int(round(np.mean(w_450_02_NumRisingLoci[ind_toplot[ind]]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(w_9k_02_NumRisingLoci[ind_toplot[ind]]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if round(np.mean(w_450_02_NumSweep[ind_toplot[ind]])) != 0 or round(np.mean(w_9k_02_NumSweep[ind_toplot[ind]])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_02_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(w_450_02_NumSweep[ind_toplot[ind]])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_02_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(w_9k_02_NumSweep[ind_toplot[ind]])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if i == 1:ax1.legend(bbox_to_anchor=(0.5, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
for ind,i in enumerate(np.arange(2,timepoints*4,4)): #s 0.05
    ax1=plt.subplot(4,4,i)
    pch_up = pchip(x,std_binned_450_w_05_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_450_w_05_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w_05[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[1])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_05_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_05_drifted[ind_toplot[ind]]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_05_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_9k_w_05_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w_05[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[5])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_05_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_05_drifted[ind_toplot[ind]]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes_sub()
    if round(np.mean(w_450_05_NumSweep[ind_toplot[ind]])) == 0 and round(np.mean(w_9k_05_NumSweep[ind_toplot[ind]])) == 0:
        ax1.annotate('%i' %int(round(np.mean(w_450_05_NumRisingLoci[ind_toplot[ind]]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(w_9k_05_NumRisingLoci[ind_toplot[ind]]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if round(np.mean(w_450_05_NumSweep[ind_toplot[ind]])) != 0 or round(np.mean(w_9k_05_NumSweep[ind_toplot[ind]])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_05_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(w_450_05_NumSweep[ind_toplot[ind]])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_05_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(w_9k_05_NumSweep[ind_toplot[ind]])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if i == 2:ax1.legend(bbox_to_anchor=(0.5, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
for ind,i in enumerate(np.arange(3,timepoints*4,4)): #s 0.08
    ax1=plt.subplot(4,4,i)
    pch_up = pchip(x,std_binned_450_w_08_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_450_w_08_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w_08[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[2])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_08_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_08_drifted[ind_toplot[ind]]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_08_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_9k_w_08_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w_08[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[6])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_08_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_08_drifted[ind_toplot[ind]]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes_sub()
    if round(np.mean(w_450_08_NumSweep[ind_toplot[ind]])) == 0 and round(np.mean(w_9k_08_NumSweep[ind_toplot[ind]])) == 0:
        ax1.annotate('%i' %int(round(np.mean(w_450_08_NumRisingLoci[ind_toplot[ind]]))), xy=(0.7, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(w_9k_08_NumRisingLoci[ind_toplot[ind]]))), xy=(0.7, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if round(np.mean(w_450_08_NumSweep[ind_toplot[ind]])) != 0 or round(np.mean(w_9k_08_NumSweep[ind_toplot[ind]])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_08_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(w_450_08_NumSweep[ind_toplot[ind]])))), xy=(0.7, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_08_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(w_9k_08_NumSweep[ind_toplot[ind]])))), xy=(0.7, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if i == 3:ax1.legend(bbox_to_anchor=(0.5, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
for ind,i in enumerate(np.arange(4,timepoints*4+1,4)): #s 0.1
    ax1=plt.subplot(4,4,i)
    pch_up = pchip(x,std_binned_450_w_1_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_450_w_1_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_w_1[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = 'turquoise', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[3])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_w_1_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_w_1_drifted[ind_toplot[ind]]),min_diff_y,'*', color = 'turquoise',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_w_1_up[ind_toplot[ind]])
    pch_down = pchip(x,std_binned_9k_w_1_down[ind_toplot[ind]])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'teal', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_w_1[ind_toplot[ind]])
    plt.plot(xx, pch(xx), color = 'teal', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[7])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_w_1_drifted[ind_toplot[ind]])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_w_1_drifted[ind_toplot[ind]]),min_diff_y,'*', color = 'teal',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes_sub()
    ax1.annotate(gen_label[ind_toplot[ind]], xy=(0.6, 0.8), xytext=(0.6, 0.8),textcoords='axes fraction',xycoords='axes fraction', fontweight = 'bold', size = 14)
    if round(np.mean(w_450_1_NumSweep[ind_toplot[ind]])) == 0 and round(np.mean(w_9k_1_NumSweep[ind_toplot[ind]])) == 0:
        ax1.annotate('%i' %int(round(np.mean(w_450_1_NumRisingLoci[ind_toplot[ind]]))), xy=(0.6, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(w_9k_1_NumRisingLoci[ind_toplot[ind]]))), xy=(0.6, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if round(np.mean(w_450_1_NumSweep[ind_toplot[ind]])) != 0 or round(np.mean(w_9k_1_NumSweep[ind_toplot[ind]])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_450_1_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(w_450_1_NumSweep[ind_toplot[ind]])))), xy=(0.6, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = 'turquoise', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(w_9k_1_NumRisingLoci[ind_toplot[ind]]))),int(round(np.mean(w_9k_1_NumSweep[ind_toplot[ind]])))), xy=(0.6, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = 'teal', fontweight = 'bold', size = 14)
    if i == 4:ax1.legend(bbox_to_anchor=(0.55, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
fig.text(0.08, 0.5, 'Proportion', va='center', rotation='vertical',size = 24) #ylabel Count Probability
fig.text(0.5, 0.05, 'Allele frequency', ha='center',size = 24) #xlabel
plt.savefig('sweep_diffS_AF_RisingLociNum_smoothCurve_sub.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('sweep_diffS_AF_RisingLociNum_smoothCurve_sub.pdf', dpi=300,format='pdf', bbox_inches = 'tight')

'''
###########
#Figure S12
###########
'''

################
#get drift data
###############

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

##################
#get simulation AF 
##################

#qff 450 different effect size: 0.04, 0.08, 0.2, 0.4
sim = 'sel_0.5-4.5--2.5-0.3'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/recom_rate_noX/1000sims/AF.sync.gz' %sim,'rb')  
AF_sim_450_qtl_04 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_450_qtl_04, AF_sim_rep_450_qtl_04 = process_AF_SimBatch(AF_sim_450_qtl_04,100,500)

sim = 'sel_0.5-4.5--6.1-0.3_eff0.08_100loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_450_qtl_08 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_450_qtl_08, AF_sim_rep_450_qtl_08 = process_AF_SimBatch(AF_sim_450_qtl_08,100,500)

sim = 'sel_0.5-4.5--16.9-0.3_eff0.2_100loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_450_qtl_2 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_450_qtl_2, AF_sim_rep_450_qtl_2 = process_AF_SimBatch(AF_sim_450_qtl_2,100,500)

sim = 'sel_0.5-4.5--34.9-0.3_eff0.4_100loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_450_qtl_4 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_450_qtl_4, AF_sim_rep_450_qtl_4 = process_AF_SimBatch(AF_sim_450_qtl_4,100,500)

#qff 9000, different effect size: 0.04, 0.08, 0.2, 0.4 
sim = 'Ne9000_sel_0.5-4.5--2.5-0.3'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/recom_rate_noX/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_9k_qtl_04 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_9k_qtl_04, AF_sim_rep_9k_qtl_04 = process_AF_SimBatch(AF_sim_9k_qtl_04,100,500)

sim = 'Ne9000_sel_0.5-4.5--6.1-0.3_eff0.08_100loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_9k_qtl_08 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_9k_qtl_08, AF_sim_rep_9k_qtl_08 = process_AF_SimBatch(AF_sim_9k_qtl_08,100,500)

sim = 'Ne9000_sel_0.5-4.5--16.9-0.3_eff0.2_100loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_9k_qtl_2 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_9k_qtl_2, AF_sim_rep_9k_qtl_2 = process_AF_SimBatch(AF_sim_9k_qtl_2,100,500)

sim = 'Ne9000_sel_0.5-4.5--34.9-0.3_eff0.4_100loci_movingTO'
InputFile_s = open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/effect_size.txt' %sim,'r')   
InputFile_af = gzip.open('/Volumes/Temp2/Dsim_base_hotcage/heterogeneity_paper/new_analysis/CS_proposal_sims/mimicree2_sims/qff/%s/500sims/AF.sync.gz' %sim,'rb')  
AF_sim_9k_qtl_4 = get_AF_qff(InputFile_s,InputFile_af)
AF_sim_class_9k_qtl_4, AF_sim_rep_9k_qtl_4 = process_AF_SimBatch(AF_sim_9k_qtl_4,100,500)

################################
#correct for drift
#get # of rising and sweep loci
################################

#QTL 450, different effect size: 0.04, 0.08, 0.2, 0.4 
AF_sim_rep_450_qtl_04_drifted, qtl_450_04_NumRisingLoci, qtl_450_04_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_qtl_04,AFC_thresh_neu_450)
AF_sim_rep_450_qtl_08_drifted, qtl_450_08_NumRisingLoci, qtl_450_08_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_qtl_08,AFC_thresh_neu_450)
AF_sim_rep_450_qtl_2_drifted, qtl_450_2_NumRisingLoci, qtl_450_2_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_qtl_2,AFC_thresh_neu_450)
AF_sim_rep_450_qtl_4_drifted, qtl_450_4_NumRisingLoci, qtl_450_4_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_450_qtl_4,AFC_thresh_neu_450)

#QTL 9k, different effect size: 0.04, 0.08, 0.2, 0.4 
AF_sim_rep_9k_qtl_04_drifted, qtl_9k_04_NumRisingLoci, qtl_9k_04_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_qtl_04,AFC_thresh_neu_9k)
AF_sim_rep_9k_qtl_08_drifted, qtl_9k_08_NumRisingLoci, qtl_9k_08_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_qtl_08,AFC_thresh_neu_9k)
AF_sim_rep_9k_qtl_2_drifted, qtl_9k_2_NumRisingLoci, qtl_9k_2_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_qtl_2,AFC_thresh_neu_9k)
AF_sim_rep_9k_qtl_4_drifted, qtl_9k_4_NumRisingLoci, qtl_9k_4_NumSweep = AF_driftCorrected_RisingSweptLoci(AF_sim_rep_9k_qtl_4,AFC_thresh_neu_9k)

#compute median AF
#QTL 450, different effect size: 0.04, 0.08, 0.2, 0.4 
median_AF_sim_rep_450_qtl_04_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_qtl_04_drifted[1:]]   
median_AF_sim_rep_450_qtl_08_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_qtl_08_drifted[1:]]   
median_AF_sim_rep_450_qtl_2_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_qtl_2_drifted[1:]]   
median_AF_sim_rep_450_qtl_4_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_450_qtl_4_drifted[1:]]   

#QTL 9k, different effect size: 0.04, 0.08, 0.2, 0.4 
median_AF_sim_rep_9k_qtl_04_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_qtl_04_drifted[1:]]   
median_AF_sim_rep_9k_qtl_08_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_qtl_08_drifted[1:]]   
median_AF_sim_rep_9k_qtl_2_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_qtl_2_drifted[1:]]   
median_AF_sim_rep_9k_qtl_4_drifted = [[np.nanmedian(rep) for rep in tmp] for tmp in AF_sim_rep_9k_qtl_4_drifted[1:]]   

######################################################
#get binned mean and standard deviation of replicates
######################################################

#QTL 450, different effect size: 0.04, 0.08, 0.2, 0.4 
mean_binned_450_qtl_04,std_binned_450_qtl_04,std_binned_450_qtl_04_up,std_binned_450_qtl_04_down = mean_binned_AF(AF_sim_rep_450_qtl_04_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_450_qtl_08,std_binned_450_qtl_08,std_binned_450_qtl_08_up,std_binned_450_qtl_08_down = mean_binned_AF(AF_sim_rep_450_qtl_08_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_450_qtl_2,std_binned_450_qtl_2,std_binned_450_qtl_2_up,std_binned_450_qtl_2_down = mean_binned_AF(AF_sim_rep_450_qtl_2_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_450_qtl_4,std_binned_450_qtl_4,std_binned_450_qtl_4_up,std_binned_450_qtl_4_down = mean_binned_AF(AF_sim_rep_450_qtl_4_drifted,bins=np.arange(0,1.05,0.05))

#QTL 9k, different effect size: 0.04, 0.08, 0.2, 0.4 
mean_binned_9k_qtl_04,std_binned_9k_qtl_04,std_binned_9k_qtl_04_up,std_binned_9k_qtl_04_down = mean_binned_AF(AF_sim_rep_9k_qtl_04_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_qtl_08,std_binned_9k_qtl_08,std_binned_9k_qtl_08_up,std_binned_9k_qtl_08_down = mean_binned_AF(AF_sim_rep_9k_qtl_08_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_qtl_2,std_binned_9k_qtl_2,std_binned_9k_qtl_2_up,std_binned_9k_qtl_2_down = mean_binned_AF(AF_sim_rep_9k_qtl_2_drifted,bins=np.arange(0,1.05,0.05))
mean_binned_9k_qtl_4,std_binned_9k_qtl_4,std_binned_9k_qtl_4_up,std_binned_9k_qtl_4_down = mean_binned_AF(AF_sim_rep_9k_qtl_4_drifted,bins=np.arange(0,1.05,0.05))

#######
#plot
########

#check how to set max value of Y axis
for i in range(14):
    print max([max(std_binned_450_qtl_04_up[i]),max(std_binned_450_qtl_08_up[i]),max(std_binned_450_qtl_2_up[i]),max(std_binned_450_qtl_2_up[i]),
               max(std_binned_9k_qtl_04_up[i]),max(std_binned_9k_qtl_08_up[i]),max(std_binned_9k_qtl_2_up[i]),max(std_binned_9k_qtl_4_up[i])])

yMax = [1.1,1.1,1,1,1,1,1,0.8,0.8,0.8,0.8,0.6,0.6,0.6]
timepoints = 14
gen_label = ['Gen %i'%tm for tm in np.arange(10,150,10)]
sim_label = ['450, eff_size=0.04','450, eff_size=0.08','450, eff_size=0.2','450, eff_size=0.4','9000, eff_size=0.04','9000, eff_size=0.08','9000, eff_size=0.2','9000, eff_size=0.4']
linestyle = ['-',':','-.','--']*14  
nbins = 20
bin_edge =np.array([float(i)/nbins for i in range(nbins+1)])
x = 0.5*(bin_edge[1:]+bin_edge[:-1])
xx = np.linspace(0, 1, 300)
          
fig, ax = plt.subplots(nrows=14, ncols=4, figsize=(20,25), facecolor='w', edgecolor='k')  
fig.subplots_adjust(hspace=0.15, wspace=0.1)
for ind,i in enumerate(np.arange(1,timepoints*4,4)): #0.04 effect size
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_qtl_04_up[ind])
    pch_down = pchip(x,std_binned_450_qtl_04_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl_04[ind])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[0])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_04_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_04_drifted[ind]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_04_up[ind])
    pch_down = pchip(x,std_binned_9k_qtl_04_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl_04[ind])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[4])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_04_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_04_drifted[ind]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    if round(np.mean(qtl_450_04_NumSweep[ind])) == 0 and round(np.mean(qtl_9k_04_NumSweep[ind])) == 0:
        ax1.annotate('%i' %int(round(np.mean(qtl_450_04_NumRisingLoci[ind]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(qtl_9k_04_NumRisingLoci[ind]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_04_NumSweep[ind])) != 0 or round(np.mean(qtl_9k_04_NumSweep[ind])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_04_NumRisingLoci[ind]))),int(round(np.mean(qtl_450_04_NumSweep[ind])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_04_NumRisingLoci[ind]))),int(round(np.mean(qtl_9k_04_NumSweep[ind])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if i == 1:ax1.legend(bbox_to_anchor=(0.4, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
for ind,i in enumerate(np.arange(2,timepoints*4,4)): #0.08 effect size
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_qtl_08_up[ind])
    pch_down = pchip(x,std_binned_450_qtl_08_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl_08[ind])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[1])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_08_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_08_drifted[ind]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_08_up[ind])
    pch_down = pchip(x,std_binned_9k_qtl_08_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl_08[ind])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[5])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_08_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_08_drifted[ind]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    if round(np.mean(qtl_450_08_NumSweep[ind])) == 0 and round(np.mean(qtl_9k_08_NumSweep[ind])) == 0:
        ax1.annotate('%i' %int(round(np.mean(qtl_450_08_NumRisingLoci[ind]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(qtl_9k_08_NumRisingLoci[ind]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_08_NumSweep[ind])) != 0 or round(np.mean(qtl_9k_08_NumSweep[ind])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_08_NumRisingLoci[ind]))),int(round(np.mean(qtl_450_08_NumSweep[ind])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_08_NumRisingLoci[ind]))),int(round(np.mean(qtl_9k_08_NumSweep[ind])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if i == 2:ax1.legend(bbox_to_anchor=(0.4, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
for ind,i in enumerate(np.arange(3,timepoints*4,4)):#0.2 effect size
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_qtl_2_up[ind])
    pch_down = pchip(x,std_binned_450_qtl_2_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl_2[ind])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[2])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_2_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_2_drifted[ind]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_2_up[ind])
    pch_down = pchip(x,std_binned_9k_qtl_2_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl_2[ind])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[6])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_2_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_2_drifted[ind]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    if round(np.mean(qtl_450_2_NumSweep[ind])) == 0 and round(np.mean(qtl_9k_2_NumSweep[ind])) == 0:
        ax1.annotate('%i' %int(round(np.mean(qtl_450_2_NumRisingLoci[ind]))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(qtl_9k_2_NumRisingLoci[ind]))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_2_NumSweep[ind])) != 0 or round(np.mean(qtl_9k_2_NumSweep[ind])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_2_NumRisingLoci[ind]))),int(round(np.mean(qtl_450_2_NumSweep[ind])))), xy=(0.8, 0.7),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_2_NumRisingLoci[ind]))),int(round(np.mean(qtl_9k_2_NumSweep[ind])))), xy=(0.8, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if i == 3:ax1.legend(bbox_to_anchor=(0.4, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
for ind,i in enumerate(np.arange(4,timepoints*4+1,4)): #0.4 effect size
    ax1=plt.subplot(14,4,i)
    pch_up = pchip(x,std_binned_450_qtl_4_up[ind])
    pch_down = pchip(x,std_binned_450_qtl_4_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = 'grey', alpha = 0.3, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_450_qtl_4[ind])
    plt.plot(xx, pch(xx), color = '#e699b2', linewidth = 3, linestyle =linestyle[i-1], label = sim_label[3])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_450_qtl_4_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_450_qtl_4_drifted[ind]),min_diff_y,'*', color = '#e699b2',markersize = 15, markeredgecolor = 'black')
    pch_up = pchip(x,std_binned_9k_qtl_4_up[ind])
    pch_down = pchip(x,std_binned_9k_qtl_4_down[ind])
    plt.fill_between(xx, pch_up(xx), pch_down(xx), color = '#cc3366', alpha = 0.1, edgecolor = None, linewidth=0.0)
    pch = pchip(x,mean_binned_9k_qtl_4[ind])
    plt.plot(xx, pch(xx), color = '#cc3366', linewidth = 3, linestyle =linestyle[i-1],label = sim_label[7])
    diff = [abs(item-np.nanmean(median_AF_sim_rep_9k_qtl_4_drifted[ind])) for item in xx]
    min_diff_y = pch(xx)[diff.index(min(diff))]
    plt.plot(np.nanmean(median_AF_sim_rep_9k_qtl_4_drifted[ind]),min_diff_y,'*', color = '#cc3366',markersize = 15, markeredgecolor = 'black')
    mask_XYaxes()
    ax1.annotate(gen_label[ind], xy=(0.8, 0.8), xytext=(0.8, 0.8),textcoords='axes fraction',xycoords='axes fraction', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_4_NumSweep[ind])) == 0 and round(np.mean(qtl_9k_4_NumSweep[ind])) == 0:
        ax1.annotate('%i' %int(round(np.mean(qtl_450_4_NumRisingLoci[ind]))), xy=(0.65, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i' %int(round(np.mean(qtl_9k_4_NumRisingLoci[ind]))), xy=(0.65, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if round(np.mean(qtl_450_4_NumSweep[ind])) != 0 or round(np.mean(qtl_9k_4_NumSweep[ind])) != 0:
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_450_4_NumRisingLoci[ind]))),int(round(np.mean(qtl_450_4_NumSweep[ind])))), xy=(0.65, 0.65),textcoords='axes fraction',xycoords='axes fraction', color = '#e699b2', fontweight = 'bold', size = 14)
        ax1.annotate('%i$\it{(%i)}$' %(int(round(np.mean(qtl_9k_4_NumRisingLoci[ind]))),int(round(np.mean(qtl_9k_4_NumSweep[ind])))), xy=(0.65, 0.5),textcoords='axes fraction',xycoords='axes fraction', color = '#cc3366', fontweight = 'bold', size = 14)
    if i == 4:ax1.legend(bbox_to_anchor=(0.4, 1), loc=3,ncol=1, borderaxespad=0.)  
    plt.xlim(0,1)
    plt.ylim(-0.02,yMax[ind])
    if yMax[ind] in [1.1,1,0.8,0.6]: plt.yticks(np.arange(0, yMax[ind]+0.1, 0.2))
fig.text(0.08, 0.5, 'Proportion', va='center', rotation='vertical',size = 24) #ylabel Count Probability
fig.text(0.5, 0.1, 'Allele frequency', ha='center',size = 24) #xlabel
plt.savefig('QTL_diffEffSize_AF_RisingLociNum_smoothCurve.png', dpi=300,format='png', bbox_inches = 'tight')
plt.savefig('QTL_diffEffSize_AF_RisingLociNum_smoothCurve.pdf', dpi=300,format='pdf', bbox_inches = 'tight')
