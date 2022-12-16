# -*- coding: utf-8 -*-

"""
Analyzes the yeast tRNA fitness landscape of Domingo et al. (Nature 2018).
Requires as an input file the file RNA_landscape_data.txt, which is identical to
Table S1 from this paper, but modified to 
eliminate the table caption, and to preserve only the following columns: 
seq, num_vars. pos_vars, ntd_vars, ntd_wt, fitness, SE.

Asks for all pairs of wild-type and mutants
with two-sided one sample t tests whether the difference 
between the mean fitness of the 
1-neighbors of the mutant and the mean fitness of the 1-neighbors of the wild type
differs significantly from zero.

Also computes for every wt-mut pair, the fraction of the
neighbors of the wt that are beneficial (have greater fitness than the wild-type),
and the fraction of neighbors of the mutant that have greater fitness than
the mutant, as well as the average benefits of these beneficial mutations. 

The output file RNA_EEmutations_out.txt contains the following columns
1. wt sequence
2. mutant (mt) sequence
3. position at which they differ
4. genotye of wt at that position
5. genotype of mt at that position
6. fitness of wt
7. fitness of mt
8. difference del_f=(f_mt-f_wt)
9. mean fitness of neighbors of wt: f_n_wt
10.mean fitness of neighbors of mt: f_n_wt
11.difference in mean fitness of the two: del_f_n=f_n_mt-f_n_wt
12.fdr value
13.p value threshold (at that fdr) for the t-test of the null hypothesis that del_f_n>del_f
14. test result: +1 (-1) if del_f_n significantly greater (smaller) than del_f, 0 otherwise
15. fraction of neighbors of wt that are beneficial
16. fraction of neighbors of mt that are beneficial
17. mean fitness of beneficial neighbors of wt
18. mean fitness of beneficial neighbors of mt
19.p value threshold (at that fdr) for the t-test of the null hypothesis that del_f_n>0
20. test result: +1 (-1) if del_f_n significantly greater (smaller) than 0, 0 otherwise




"""

import EE_funcs

import pandas as pd
from pathlib import Path


import numpy as np
from scipy import stats

import igraph as gr

#to be able to exit the program with sys.exit("error message")
import sys

#########################################################
### extract RNA data 
#########################################################

fitdatfile="RNA_landscape_data.txt"
pathstr="./"
filepath = pathstr + fitdatfile
infile = Path(filepath)

df =  pd.read_csv(infile, sep='\t')
fitdatall=df.T.to_dict('list')

#two dictionaries indexed by sequence that will hold the fitness and its standard error 
fit={}
varfit={}
for n in fitdatall.keys():
    seq=fitdatall[n][0]
    varsites=seq[0]+seq[1]+seq[5]+seq[26]+seq[42]+seq[45]+seq[65]+seq[68]+seq[69]+seq[70]
    #fitness is defined in the log domain
    logfit=fitdatall[n][6]
    fit[varsites]= logfit
    
    #transform standard error into a variance
    #data is based on six measurements, as reported by Domingo et al., Nature 2018  
    tmpserrfit=fitdatall[n][7]
    tmpsdevfit=tmpserrfit*np.sqrt(6) 
    tmpvarfit=tmpsdevfit**2
    varfit[varsites] = tmpvarfit
    

# a sorted list of all variable sites
allseq=sorted(fit.keys())

#mapping the index of variable sites from 0 to 9 onto actual tRNA positions
posvarseq_to_tRNApos={}
posvarseq_to_tRNApos[0]=1
posvarseq_to_tRNApos[1]=2
posvarseq_to_tRNApos[2]=6
posvarseq_to_tRNApos[3]=27
posvarseq_to_tRNApos[4]=43
posvarseq_to_tRNApos[5]=46
posvarseq_to_tRNApos[6]=66
posvarseq_to_tRNApos[7]=69
posvarseq_to_tRNApos[8]=70
posvarseq_to_tRNApos[9]=71


#####################################################################
#now write a header line for the main output file
######################################################################

#the false discovery rate(s) (fdr) to be used for the Benjamini-Hochberg
#multiple testing correction
fdrarr=[0.01]


outfilename="RNA_EEmutations_out.txt"
print("will write data output to", outfilename)
outfile = open(outfilename, 'w')
print("wt\tmut\tposdiff\tg_wt\tg_mut1\t", end ='', file = outfile)
print("f_wt\tf_mut\tdelfit\t", end ='', file = outfile)
print("meanf_wt_neigh\tmeanf_mut_neigh\tdelmeanfit\t", end ='', file = outfile)

for i in range(0,len(fdrarr)):
    print('fdr_'+str(i), end ='\t', file = outfile)

    print('pthresh_del_neigh_fit_del_fit_'+str(i), end ='\t', file = outfile)
    print('del_neigh_fit_del_fit_'+str(i), end ='\t', file = outfile)
    #fraction of beneficial neighbors of both wt and mutant
    print('frac_neigh_wt_ben'+str(i), end ='\t', file = outfile)
    print('frac_neigh_mut_ben'+str(i), end ='\t', file = outfile)
    #mean benefit of beneficial neighbors of both wt and mutant
    print('mean_ben_neigh_wt'+str(i), end ='\t', file = outfile)
    print('mean_ben_neigh_mut'+str(i), end ='\t', file = outfile)
    
    print('pthresh_del_neigh_fit_0_'+str(i), end ='\t', file = outfile)
    print('del_neigh_fit_0_'+str(i), end ='\t', file = outfile)

print('', end ='\n', file = outfile)

outfile.close()


 
########################################################################################
#### core computation block to prepare input for t-tests and to compute test statistics
######################################################################################## 

print("now beginning statistical testing")

#dicts that hold p values of test of the null hypotheses that the difference between the mean fitness of the 
#1-neighbors of the mutant and the mean fitness of the 1-neighbors of the wild type
#equals (i) the difference of fitness between mutant and wild type, or (ii) 0
p_del_neigh_fit_del_fit={}
p_del_neigh_fit_0={} 
 
#position at which wild type and mutant differ
diffpos={}

#mean fitness values, sdev, and n for 1 neighborhoods of sequences
mfit_1neigh={}
sdev_1neigh={}
n_1neigh={} #number of sequences in one neighborhood


#compute mean and standard deviation of the neighbors of all sequences 
#by looping over all wild type sequences, we also cycle over all
#possible mutants
ctr=0
for wt in allseq:
    mfit_1neigh[wt]= {}
    sdev_1neigh[wt]= {}
    n_1neigh[wt]={}
    #there are 10 variable positions in the domingo data set
    for pos in range(0, 10):
        #compute mean and standard deviation of fitness of neighbors of wt
        #under the exclusion of variants at the position pos
        fitstat_wt_posconst=EE_funcs.mfit_sdevm_1mut_excludepos(wt, fit, varfit, pos)
        mfit_1neigh[wt][pos]= fitstat_wt_posconst[0]
        sdev_1neigh[wt][pos]= fitstat_wt_posconst[1]
        n_1neigh[wt][pos]= fitstat_wt_posconst[2]
    ctr +=1
        
#compute input for t-test
#cycle over all sequences in the data set and their 1-mutant neighbors
#notice that it is necessary to cycle over all sequence pairs, and not just i<j
ctr=0
for wt in allseq:
    
        
        p_del_neigh_fit_del_fit[wt]={}
        p_del_neigh_fit_0[wt]={} 
                
        diffpos[wt]={}        
        
        for mut1 in allseq:
            if EE_funcs.hamdist(wt,mut1)==1:
               
                #determine at which residue the mutant differs from the wild-type
                #needed below to hold that residue constant
                for pos in range(len(mut1)):
                    if wt[pos] != mut1[pos]:
                        diffpos[wt][mut1]=pos
                        break
                  
                #now use a one-sample t-test to find out whether the difference between the mean fitness of the 
                #1-neighbors of the mutant and the mean fitness of the 1-neighbors of the wild type
                #differs significantly from the single number fit[mut1]-fit[wt]. 
                tmpmean=mfit_1neigh[mut1][pos]-mfit_1neigh[wt][pos]
                #the standard deviation of the difference of two random variables
                tmpsdev=np.sqrt((sdev_1neigh[mut1][pos])**2 + (sdev_1neigh[mut1][pos])**2 )
                #sample size is based on the actual number of neighbors, 
                #since the random variate is a difference of two random variates that may have different sample sizes, 
                #just use the smaller sample size to render the test conservative
                tmpn=np.min([n_1neigh[wt][pos], n_1neigh[mut1][pos]])
                tmpstat=stats.ttest_ind_from_stats(mean1=tmpmean, std1=tmpsdev, nobs1=tmpn, mean2=fit[mut1]-fit[wt], std2=0, nobs2=2, equal_var=False)
                p_del_neigh_fit_del_fit[wt][mut1]=tmpstat.pvalue 
                
                #now use the same one-sample t-test to find out whether the difference between the mean fitness of the 
                #1-neighbors of the mutant and the mean fitness of the 1-neighbors of the wild type
                #differs significantly from zero. 
                tmpstat=stats.ttest_ind_from_stats(mean1=tmpmean, std1=tmpsdev, nobs1=tmpn, mean2=0, std2=0, nobs2=2, equal_var=False)
                p_del_neigh_fit_0[wt][mut1]=tmpstat.pvalue 
                
                ctr=ctr+1

#once all p-values are calculated, perform a BH FDR correction
#for this, first write all p values into an array
#do this separately for the different kinds of tests above
p_del_neigh_fit_del_fit_allp=[]
p_del_neigh_fit_0_allp=[]

for wt in p_del_neigh_fit_del_fit.keys():
    for mut1 in p_del_neigh_fit_del_fit[wt].keys():
        #order of p values does not matter, since we will be sorting later
        p_del_neigh_fit_del_fit_allp.append(p_del_neigh_fit_del_fit[wt][mut1])
        p_del_neigh_fit_0_allp.append(p_del_neigh_fit_0[wt][mut1])
        

print("computing p-value thresholds for Benjamini-Hochberg FDR correction")
#a dictionary keyed by fdr value that holds the p value threshold for that value
pthresh_del_neigh_fit_del_fit={}
pthresh_del_neigh_fit_0={}
#note that it is in principle possible to use an array of multiple fdr values in 
#the program allows multiple fdr values to be used 
for fdr in fdrarr:
    pthresh_del_neigh_fit_del_fit[fdr] = EE_funcs.benj_hoch(p_del_neigh_fit_del_fit_allp, fdr)
    pthresh_del_neigh_fit_0[fdr] = EE_funcs.benj_hoch(p_del_neigh_fit_0_allp, fdr)
    
#######################################################
#analysis of t-test results and output
#######################################################
print("writing results to output file")
outfile = open(outfilename, 'a')
tmpctr=0
for wt in allseq:
    for mut1 in allseq: 
        if EE_funcs.hamdist(wt,mut1)==1:
            tmpctr+=1            
            delfit=fit[mut1]-fit[wt]
            
            #the position at which wild type and mutant differ
            pos=diffpos[wt][mut1]
            del_mfit_1neigh = mfit_1neigh[mut1][pos]-mfit_1neigh[wt][pos]
                
            print(wt, end ='\t', file = outfile)
            print(mut1, end ='\t', file = outfile)
            #the actual position that differs between wt and mut1
            print(posvarseq_to_tRNApos[pos], end ='\t', file = outfile)
            print(wt[pos], end ='\t', file = outfile)
            print(mut1[pos], end ='\t', file = outfile)
              
            print('{0:.4f}\t'.format(fit[wt]), end ='', file = outfile)
            print('{0:.4f}\t'.format(fit[mut1]), end ='', file = outfile)
            print('{0:.4f}\t'.format(delfit), end ='', file = outfile)
    
            print('{0:.4f}\t'.format(mfit_1neigh[wt][pos]), end ='', file = outfile)
            print('{0:.4f}\t'.format(mfit_1neigh[mut1][pos]), end ='', file = outfile)
            print('{0:.4f}\t'.format(del_mfit_1neigh), end ='', file = outfile)
     
        
            #now print the current fdr value for which the statistical analysis is 
            #to be conducted
            for fdr in pthresh_del_neigh_fit_del_fit.keys():
                print('{0:.4f}\t'.format(fdr), end ='', file = outfile)  
                
                
                #first output the p value threshold for this fdr
                print('{0:.4f}\t'.format(pthresh_del_neigh_fit_del_fit[fdr]), end ='', file = outfile)

                #the following is relevant to identify beneficial or neutral EE mutations (fit[mut1]-fit[wt]>=0):
                #CASE 1.1:  difference in mean fitness of the neighbors of the wild-type 
                #and the mutant is significantly different from delfit                
                if p_del_neigh_fit_del_fit[wt][mut1] <= pthresh_del_neigh_fit_del_fit[fdr]: 
                    #if mean fitness difference of neighbors exceeds delfit
                    #write +1, if delfit>0, this means that the mutation is evolvability enhancing
                    if del_mfit_1neigh>delfit:
                        print('1\t', end ='', file = outfile)
                    #in the opposite case, write -1, if delfit>0, the mutation is evolvability reducing
                    elif del_mfit_1neigh<delfit:
                        print('-1\t', end ='', file = outfile) 
                    #sanity check: no fitness difference despite a significant p-value: impossible
                    else:
                        print('error_aw in main')
                        sys.exit()
                #CASE 1.2:  difference in mean fitness of the neighbors of the wild-type 
                #and the mutant is NOT significantly different from delfit  
                if p_del_neigh_fit_del_fit[wt][mut1] > pthresh_del_neigh_fit_del_fit[fdr]:                    
                    print('0\t', end ='', file = outfile)
                      
                #note the following block is independent of fdr, and will result in the same
                #output of the fraction of beneficial neighbors and their mean benefit for each fdr
                #creates redundant output but done to keep data on each fdr level in one block of 
                #the output file
            
                #now cycle over all 1-neighbors of the wt except mut1 and ask whether its 
                #fitness is greater than the wild type
                neigh_ctr=0
                #the fraction of beneficial mutations
                frac_ben_wt_neighs=0
                #the mean benefit of those mutations that are beneficial
                mean_ben_wt_neighs=0
                for tmpseq in allseq:
                    #the third and last condition requires that the candidate neighbor tmpseq
                    #does not differ from the wildtype at the same residue as mut1, 
                    #this is important for loci with more than two alleles, and it is included here
                    #for consistency to how neighborhoods are determined in the first pass
                    #through the data above
                    if EE_funcs.hamdist(tmpseq, wt)==1 and tmpseq != mut1 and tmpseq[pos]==wt[pos]:
                        neigh_ctr+=1
                        if fit[tmpseq]>fit[wt]: 
                            frac_ben_wt_neighs+=1
                            mean_ben_wt_neighs=fit[tmpseq]-fit[wt]
                #fraction of neighbors of the wild type that are beneficial
                frac_ben_wt_neighs /= neigh_ctr 
                #mean fitness benefit of the beneficial mutants                           
                mean_ben_wt_neighs /= neigh_ctr   
                
                #do the same for the mutants (examine all neighbors except the wt)
                neigh_ctr=0
                frac_ben_mut1_neighs=0
                mean_ben_mut1_neighs=0
                for tmpseq in allseq:
                    if EE_funcs.hamdist(tmpseq, mut1)==1 and tmpseq != wt and tmpseq[pos]==mut1[pos]:
                        neigh_ctr+=1
                        if fit[tmpseq]>fit[mut1]: 
                            frac_ben_mut1_neighs += 1
                            mean_ben_mut1_neighs=fit[tmpseq]-fit[mut1]
                frac_ben_mut1_neighs /= neigh_ctr
                mean_ben_mut1_neighs /= neigh_ctr   
                            
                
                print('{0:.4f}\t'.format(frac_ben_wt_neighs), end ='', file = outfile) 
                print('{0:.4f}\t'.format(frac_ben_mut1_neighs), end ='', file = outfile)
                #the mean benefit of beneficial mutations is only meaningful when 
                #there are beneficial mutations
                if frac_ben_wt_neighs>0:
                    print('{0:.4f}\t'.format(mean_ben_wt_neighs), end ='', file = outfile) 
                else:
                    print('NaN\t', end ='', file = outfile)
                if frac_ben_mut1_neighs>0:
                    print('{0:.4f}\t'.format(mean_ben_mut1_neighs), end ='', file = outfile) 
                else:
                    print('NaN\t', end ='', file = outfile)
        
                #the following is relevant to identify deleterious EE mutations (fit[mut1]-fit[wt]<0):
                #CASE 2.1:  difference in mean fitness of the neighbors of the wild-type 
                #and the mutant is significantly different from zero              
                print('{0:.4f}\t'.format(pthresh_del_neigh_fit_del_fit[fdr]), end ='', file = outfile)
                if p_del_neigh_fit_0[wt][mut1] <= pthresh_del_neigh_fit_0[fdr]: 
                    #if mean fitness difference of neighbors exceeds zero
                    #write +1, if delfit<0, this means that the mutation is evolvability enhancing
                    if del_mfit_1neigh>0:
                        print('1\t', end ='', file = outfile) 
                    #if mean fitness difference of neighbors <0
                    #write -1, if delfit<0, this means that the mutation is evolvability reducing
                    elif del_mfit_1neigh<0:
                        print('-1\t', end ='', file = outfile) 
                    #sanity check: no fitness difference despite a significant p-value: impossible
                    else:
                        print('error_aw in main')
                        sys.exit()
                #CASE 2.2:  difference in mean fitness of the neighbors of the wild-type 
                #and the mutant is NOT significantly different from zero 
                if p_del_neigh_fit_0[wt][mut1] > pthresh_del_neigh_fit_0[fdr]: 
                        print('0\t', end ='', file = outfile)
                
                
                        
            print('\n', end ='', file = outfile)

outfile.close()
    
            

