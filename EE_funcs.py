# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 13:09:34 2020

@author: aw
"""
import numpy as np #for poisson distribution
import sys 
import math
import random # for random number generators


#does a Benjamini hochberg correction on a list of P-values
#sorts the list of p-values, then computes the maximal rank for
#which the p-value is lower than the BH threshold. Returns the p-value 
#corresponding to this rank, all p-values lower than  or equal to that value
# will be considered significant
def benj_hoch(plist, fdr):
    plist_sorted=sorted(plist)
    n=len(plist_sorted)
    #the maximal index of the p-value that is below the BH threshold
    maxindex=-1
    for rank in range(0, n):
        BH_thresh=fdr*((rank+1)/n) #rank must run from 1 to n
        if plist_sorted[rank]<BH_thresh:
            maxindex=rank
    #the maximal p value that is considered significant
    if maxindex> -1:
        return plist_sorted[maxindex]
    else: #this is the case where no p value is below the BH threshold, in this case all p-values are 
        #nonsignificant   
        return 0


#the Hamming distance of two strings
def hamdist(str1, str2):
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

#determines the smallest hamming distance between a string and an array of strings
def mindist(str, strarr):
    mind=len(str)
    for x in strarr:
        currdist=hamdist(x, str)
        if currdist < mind:
            mind=currdist        
    return mind

# Kimura's formula for the fixation of a mutation under selection and drift,
# adapted for haploid populations
def pfix(N, s):
    if(s != 0):
        try: 
            x=(1-math.exp(-2*s))/(1-math.exp(-2*N*s))
        #can happen if s<0 and N large, in which case the denominator becomes too large to handle and pfix should be zero
        except OverflowError: 
            x=0
        return x
    else:
        return 1/N


#calculate the fixation probabilities for an entire network according
#to Kimura's fixation probability
def pfix_gnet_Kimura(gnet, N):
    for g in gnet.vs:
        sumprob=0
        
        for e in g.out_edges():
            tmp=pfix(N, e["fitdiff"])
            e["pfix"]=tmp
            sumprob+=tmp
            
        #when fitness differences are negative and with large 
        #population sizes, then all fixation probabilities may be zero
        #this means we are at a peak, and we simply set all
        #fixation probabilities to zero, without normalization
        if sumprob==0:
            for e in g.out_edges():
                e["pfix"]=0
        else:
            #normalize fixation probabilities to one
            for e in g.out_edges():
                e["pfix"]/=sumprob
                
                

            
# a random walk of L steps, using Kimura fixation probabilities, starting
# from a genotype g, assumes that fixation probabilities are already associated 
# with each edge of the graph
# returns an array with three entries:
# 1. the fitness at the end of the random walk, 
# 2. a dictionary with the number of edges of a specific type traversed, 
# 3. the actual length of the walk, which is important if a peak has been reached
# and the walk needed to be aborted early 
def ranwalk_Kimura(gnet, g, L):   
    
    #the random walker we start with
    ranwalker=g
    currfit=ranwalker['fit'] #the current fitness of the random walker
    
    #edges traversed
    e_trav={}    
    e_trav['EE_edges']=0 #the number of evolvability enhancing mutations along this random walk
    e_trav['EN_edges']=0 #the number of evolvability reducing mutations along this random walk
    e_trav['ED_edges']=0 #the number of evolvability reducing mutations along this random walk
    #and a finer grouping of edges traversed
    e_trav['ben_EE_edges']=0
    e_trav['ben_EN_edges']=0 
    e_trav['ben_ED_edges']=0
    e_trav['neut_EE_edges']=0 
    e_trav['neut_EN_edges']=0
    e_trav['neut_ED_edges']=0             
    e_trav['del_EE_edges']=0
    e_trav['del_EN_edges']=0
    e_trav['del_ED_edges']=0   
    
    #the actual number of steps will be L+1, so that we admit the case of a 
    #random walk with zero steps, i.e., one that starts at a peak
    for step in range(0, L+1):
        outedges=ranwalker.out_edges()
        pfixarr=[]
        for e in ranwalker.out_edges():
            pfixarr.append(e['pfix'])
        
        
        #a potential problem here is that all fixation probabilites may be 
        #zero, in which case we have reached a peak and the random walk should
        #abort or it will get hung
        if np.sum(pfixarr)==0:
            return [currfit, e_trav, step]
        
        n_neigh=len(pfixarr)
        ran_neigh=random.randrange(n_neigh)
        while random.random()>=pfixarr[ran_neigh]:
            ran_neigh=random.randrange(n_neigh)
            
        #to determine whether the mutation we chose was 
        #evovability-enhancing
        if outedges[ran_neigh]['evolv'][-2:]=='EE':
            e_trav['EE_edges']+=1
        if outedges[ran_neigh]['evolv'][-2:]=='EN':
            e_trav['EN_edges']+=1
        if outedges[ran_neigh]['evolv'][-2:]=='ED':
            e_trav['ED_edges']+=1
        
        if outedges[ran_neigh]['evolv']=='ben_EE':
            e_trav['ben_EE_edges']+=1
        if outedges[ran_neigh]['evolv']=='ben_EN':
            e_trav['ben_EN_edges']+=1
        if outedges[ran_neigh]['evolv']=='ben_ED':
            e_trav['ben_ED_edges']+=1
        if outedges[ran_neigh]['evolv']=='neut_EE':
            e_trav['neut_EE_edges']+=1
        if outedges[ran_neigh]['evolv']=='neut_EN':
            e_trav['neut_EN_edges']+=1
        if outedges[ran_neigh]['evolv']=='neut_ED':
            e_trav['neut_ED_edges']+=1
        if outedges[ran_neigh]['evolv']=='del_EE':
            e_trav['del_EE_edges']+=1
        if outedges[ran_neigh]['evolv']=='del_EN':
            e_trav['del_EN_edges']+=1
        if outedges[ran_neigh]['evolv']=='del_ED':
            e_trav['del_ED_edges']+=1
        
        
        #now change the state of the random walker to the target of
        #the chosen edge, note that the vs function returns an edge sequence (in this case with one member),
        #not an edge, hence the subscript[0]
        ranwalker=gnet.vs(outedges[ran_neigh].target)[0]
        currfit=ranwalker['fit']
        
    return [currfit, e_trav, step]





# like ranwalk_Kimura, but returns the entire walk, for plotting
# returns
# 1. the fitness of the current random walker(excluding the starting genotype passed into the function) 
# 2. a list of strings containing 'EE', 'EN' ,or 'ED' depending on whether the 
#    step leading to the current random walker was evolvability enhancing, neutral, or reducing
def ranwalk_Kimura_verbose(gnet, g, L):   
    fitvalarr=[]
    evolvarr=[]    
    
    
    #the random walker we start with
    ranwalker=g
    currfit=ranwalker['fit'] #the current fitness of the random walker
    #edges traversed
    e_trav={}    
    e_trav['EE_edges']=0 #the number of evolvability enhancing mutations along this random walk
    e_trav['EN_edges']=0 #the number of evolvability reducing mutations along this random walk
    e_trav['ED_edges']=0 #the number of evolvability reducing mutations along this random walk
    #and a finer grouping of edges traversed
    e_trav['ben_EE_edges']=0
    e_trav['ben_EN_edges']=0 
    e_trav['ben_ED_edges']=0
    e_trav['neut_EE_edges']=0 
    e_trav['neut_EN_edges']=0
    e_trav['neut_ED_edges']=0             
    e_trav['del_EE_edges']=0
    e_trav['del_EN_edges']=0
    e_trav['del_ED_edges']=0   
    
    #the actual number of steps will be L+1, so that we admit the case of a 
    #random walk with zero steps, i.e., one that starts at a peak
    for step in range(0, L+1):
        outedges=ranwalker.out_edges()
        
        pfixarr=[]
        for e in ranwalker.out_edges():
            pfixarr.append(e['pfix'])
        
       
        
        #a potential problem here is that all fixation probabilites may be 
        #zero, in which case we have reached a peak and the random walk should
        #abort or it will get hung
        if np.sum(pfixarr)==0:
            return [fitvalarr, evolvarr]
        
        n_neigh=len(pfixarr)
        ran_neigh=random.randrange(n_neigh)
        while random.random()>=pfixarr[ran_neigh]:
            ran_neigh=random.randrange(n_neigh)
            
        #to determine whether the mutation we chose was 
        #evovability-enhancing
        if outedges[ran_neigh]['evolv'][-2:]=='EE':
            e_trav['EE_edges']+=1
            evolvarr.append('EE')
        if outedges[ran_neigh]['evolv'][-2:]=='EN':
            e_trav['EN_edges']+=1
            evolvarr.append('EN')
        if outedges[ran_neigh]['evolv'][-2:]=='ED':
            e_trav['ED_edges']+=1
            evolvarr.append('ED')    
        
        
        if outedges[ran_neigh]['evolv']=='ben_EE':
            e_trav['ben_EE_edges']+=1
        if outedges[ran_neigh]['evolv']=='ben_EN':
            e_trav['ben_EN_edges']+=1
        if outedges[ran_neigh]['evolv']=='ben_ED':
            e_trav['ben_ED_edges']+=1
        if outedges[ran_neigh]['evolv']=='neut_EE':
            e_trav['neut_EE_edges']+=1
        if outedges[ran_neigh]['evolv']=='neut_EN':
            e_trav['neut_EN_edges']+=1
        if outedges[ran_neigh]['evolv']=='neut_ED':
            e_trav['neut_ED_edges']+=1
        if outedges[ran_neigh]['evolv']=='del_EE':
            e_trav['del_EE_edges']+=1
        if outedges[ran_neigh]['evolv']=='del_EN':
            e_trav['del_EN_edges']+=1
        if outedges[ran_neigh]['evolv']=='del_ED':
            e_trav['del_ED_edges']+=1
        
        
        #now change the state of the random walker to the target of
        #the chosen edge, note that the vs function returns an edge sequence (in this case with one member),
        #not an edge, hence the subscript[0]
        ranwalker=gnet.vs(outedges[ran_neigh].target)[0]
        currfit=ranwalker['fit']
        fitvalarr.append(currfit)
        
    return [fitvalarr, evolvarr]
        
      



# determines the mean and sdev of the mean fitness of all 1-mutant neighbors of a sequence
# minus those sequences whose residue at position pos differs from that residue in the  
# wild-type, this can be important for multiallelic variants where one wants to have 
# a specific residue at a specific position
    
# seqfit is a 1-d dictionary (key, value)=(seq, fitness)
# vararr is a 1-d dictionary (key, value)=(seq, fitness variance)  

#note that the s.dev of the mean of n independent rvs xi with variance var i
#calculates as follows var (sum xi) = sum (var xi), thus s.dev (sum xi)=sqrt (sum var xi)
#to get the standard deviation of the mean, divide sqrt (sum var xi) by n

#returns an array [mean fitness, sdev of mean fitness, n(number of neighbors)]
def mfit_sdevm_1mut_excludepos(wt, seqfit, varfit, pos):
    
    fitarr=[]
    vararr=[]
    for seq in seqfit.keys():
        if(hamdist(seq, wt)==1 and seq[pos] == wt[pos]):
            #print(seq, wt, seqfit[seq])
            fitarr.append(seqfit[seq])
            vararr.append(varfit[seq])                
    if len(fitarr)==0:
        return [math.nan, math.nan, 0]
    else:
        #now calculate the average 
        sdev_of_mfit=0
        for var in vararr:
            sdev_of_mfit += var
        sdev_of_mfit = (1/len(vararr))*np.sqrt(sdev_of_mfit)
        #print("elements in fitarr :", len(fitarr))
        return [ np.mean(fitarr), sdev_of_mfit, len(fitarr) ]

# determines the mean and sdev of the fitness of all 1-mutant neighbors of a sequence
# minus those sequences whose residue at position pos differs from that residue in the  
# wild-type, this can be important for multiallelic variants where one wants to have 
# a specific residue at a specific position

#passes a data structure neigh that contains
#all sequences in the neighborhood of the passed sequence for computational
#efficiency
def mfit_1mut_excludepos_v3(wt, neigh, seqfit, pos):   
    fitarr=[]
    for seq in neigh[wt]:
        if seq[pos] == wt[pos]:
            fitarr.append(seqfit[seq])
    return [ np.mean(fitarr), np.std(fitarr), len(fitarr) ]


    
    

 