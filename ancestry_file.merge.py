##################################################################################
## python2.7
##################################################################################
## Description: calibrate and merge ancestry components in acnestry files
##                sort ancestry components and estimate supporting ratio
##                #calibrate and merge P files for contml function in phylip
##################################################################################
## Usage:        python2.7 ancestry_file.merge.py -h
##################################################################################
## Output:        prefix.logfile
##                prefix.consensus.filelist
##                prefix.conflic.filelist
##                prefix.merge.ancestry
##                prefix.supporting_ratio.txt
##################################################################################

import argparse, sys
import pandas
import numpy as np
from scipy.stats.stats import pearsonr
from scipy import stats
import glob, os
import time
import itertools
from decimal import Decimal
import random

##################################################################################
#######################           ancestry files           #######################
##################################################################################
def find_files(root_path, K):
    filelist = []
    for root,subdirs,files in os.walk(root_path):
        if len(files) > 0:
            for file in files:
                if len(file.split('.')) >= 2:
                    if (file.split('.')[-1] == 'ancestry') & (file.split('.')[-2] == str(K)): 
                        filelist.append(os.path.join(root,file))
                    else:
                        pass
                else:
                    pass
        else:
            pass
    return filelist

##################################################################################
#######################                match               #######################
##################################################################################
def compare_files(ref_data, target_ancestry_file_name, K, out):
    match = {}  ## pairs of ancestry components
    ancref = ref_data.copy()
    ancref.index = ancref[0]
    anctar = pandas.read_csv(target_ancestry_file_name,sep='\s+',header=None)
    anctar.index = anctar[0]
    overlapsample = list(set(ancref.index) & set(anctar.index))

    for col1 in list(ancref.columns)[2:]:
        largest_cor = 0.0
        second_cor = 0.0
        for col2 in list(anctar.columns)[2:]:
            cor = pearsonr(ancref.loc[overlapsample,col1], anctar.loc[overlapsample,col2])[0]
            if cor > largest_cor:
                largest_cor, second_cor = cor, largest_cor
                match[col1] = col2
            elif (cor <= largest_cor) & (cor > second_cor):
                second_cor = cor
            else:
                pass
        if (abs(largest_cor) - abs(second_cor)) < 0.5:
            del match[col1]
        else:
            pass

    if len(set(match.keys())) > len(set(match.values())):
        multimatch = [(x,y) for x,y in list(itertools.combinations(match.keys(),2)) if match[x]==match[y]]
        multimatch = list(set(itertools.chain.from_iterable(multimatch)))
        match = {x:y for x,y in match.items() if not x in multimatch}
    else:
        pass

    if len(match.keys()) == K:
        anctar = anctar[[0,1] + [match[x] for x in list(ancref.columns)[2:]]]
        anctar.columns = list(ancref.columns)
        with open(out+'.consensus.filelist','a') as fw:
            fw.write(target_ancestry_file_name+'\n')
        return anctar
    else:
        with open(out+'.conflict.filelist','a') as fw:
            fw.write(target_ancestry_file_name+'\n')
        return match

##################################################################################
#######################                 merge              #######################
##################################################################################
def merge_files(ref_data, target_ancestry_file_list, K, out):
    n = len(target_ancestry_file_list)
    result = map(compare_files, [ref_data]*n, target_ancestry_file_list, [K]*n, [out]*n)
    ## estimate supporting ratio for each component
    matching = [x for x in result if isinstance(x,dict)]
    counting = list(itertools.chain(*[x.keys() for x in matching])) + range(2,K+2)*(len(target_ancestry_file_list)-len(matching))
    counting = pandas.DataFrame(pandas.value_counts(counting),columns=['count'])
    counting['ratio'] = counting['count']*1.0 / len(target_ancestry_file_list)
    ## merge consensus ancestry file
    ancestry = [y for y in result if isinstance(y,pandas.DataFrame)]
    merging = pandas.concat(ancestry,ignore_index=True)
    merging = merging.groupby([0,1]).mean().reset_index()
    return counting, merging

##################################################################################
#######################                 sort               #######################
##################################################################################
def output_ancestry(path, K, ref_data, order, out):
    tar_anc_filelist = find_files(path, K)
    counts, merge = merge_files(ref_data, tar_anc_filelist, K, out)
    ## sort populations
    if order != 'NO':
        poplist = pandas.read_csv(order,header=None)
        poplist = list(poplist[0])
    else:
        poplist = list(ref_data[1].unique())
    poplist = [x for x in poplist if x in list(merge[1])]
    merge[1] = merge[1].astype('category')
    merge[1].cat.reorder_categories(poplist, inplace=True)
    merge.sort_values(1,ascending=True,inplace=True)
    ## sort individuals within each population
    remerge = pandas.DataFrame()
    for pop in poplist:
        temp = merge[merge[1]==pop].copy()
        mean = list(temp.iloc[:,2:].mean())
        sort = list(np.argsort(mean)+2)
        sort.reverse()
        temp.sort_values(by=sort,ascending=True,inplace=True)
        remerge = pandas.concat([remerge,temp],ignore_index=True)
    merge = remerge.copy()
    ## sort components
    popmean = merge.groupby([1]).mean()
    resort = pandas.DataFrame(columns=['component','pop','index'])
    for comp in list(popmean.columns):
        resort.loc[comp,'index'] = list(popmean[comp]).index(popmean[comp].max())
    resort['component'] = resort.index
    resort['pop'] = resort['index'].apply(lambda x: list(popmean.index)[x])
    resort.sort_values(by=['index'],ascending=True,inplace=True)
    resortlist = list(resort['component'])
    merge = merge[[0,1]+resortlist]
    merge.columns = [0,1] + range(2,K+2)
    merge.to_csv(out+'.merge.ancestry',sep='\t',header=None,index=None)
    ## update supporting ratio
    resort['count'] = resort['component'].apply(lambda x: counts.loc[x,'count'])
    resort['ratio'] = resort['component'].apply(lambda x: counts.loc[x,'ratio'])
    resort['component'] = ['comp'+str(i) for i in range(1,K+1)]
    resort.drop('index',axis=1,inplace=True)
    resort.columns = ['component','represent_pop','support_counts','support_ratio']
    resort.to_csv(out+'.supporting_ratio.txt',sep='\t',index=None)
    ## finish
    return merge

##################################################################################
#######################      calibrate & merge P file      #######################
##################################################################################
'''
def calibrate_merge_P_file(target_ancestry_file_list, merged_ancestry_data, K, out):
    with open(out+'.contml.phy','w') as fw:
        pass
    ancref = merged_ancestry_data.copy()
    ancref.index = ancref[0]
    for anctarfile in target_ancestry_file_list:
        match = {}  ## pairs of ancestry components
        anctar = pandas.read_csv(anctarfile,sep='\s+',header=None)
        anctar.index = anctar[0]
        overlapsample = list(set(ancref.index) & set(anctar.index))

        for col1 in list(ancref.columns)[2:]:
            largest_cor = 0.0
            for col2 in list(anctar.columns)[2:]:
                cor = pearsonr(ancref.loc[overlapsample,col1], anctar.loc[overlapsample,col2])[0]
                if cor > largest_cor:
                    largest_cor = cor
                    match[col1] = col2
                else:
                    pass
        resort = [match[x]-2 for x in list(ancref.columns)[2:]]

        pdata = pandas.read_csv(os.path.splitext(anctarfile)[0]+'.P',sep='\s+',header=None)
        pdata = pdata[resort]
        pdata.columns = ['comp'+str(i) + ' '*(6-len(str(i))) for i in range(1,K+1)]
        #pdata.rename(columns=lambda x: x+' '*(10-len(str(x))),inplace=True)
        pdata.round(decimals=6).to_csv(os.path.splitext(anctarfile)[0]+'.calibrate.P',sep=' ',header=None,index=None)

        samplenum = pdata.shape[1]; markernum = pdata.shape[0]
        pdata = pdata.applymap(lambda x: str(Decimal(str(x)).quantize(Decimal('0.000000')))).T
        pdata.insert(0,'component',pdata.index)

        with open(out+'.contml.phy','a') as fw:
            fw.write('    '+str(samplenum)+'    '+str(markernum)+'\n')
            fw.write(' '.join(['2']*markernum)+'\n')
        #pdata = pdata.sample(frac=1)
        pdata.apply(lambda x: ' '.join(x),axis=1).to_csv(out+'.contml.phy',header=None,index=None,mode='a')
        with open(out+'.contml.phy','a') as fw:
            fw.write('\n')
'''

##################################################################################
#######################                 main               #######################
##################################################################################
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", type=str, required = True, \
                        help="/path/for/ancestry/files/from/AncestryPainter.pl, ancestry file names: prefix.K.ancestry")
    parser.add_argument("--ref", type=str, required = True, \
                        help="/reference/ancestry.file")
    parser.add_argument("--k", type=int, required=True, \
                        help="Number of ancestry component.")
    parser.add_argument("--order", type=str, required=False, default='NO', \
                        help="Population order list.")
    #parser.add_argument("--p", type=str, required=False, choices=['T','F'], default='F', \
    #                    help="To merge P files for phylip (contml), with each P file in the same floder of the corresponding ancestry files. better to randomize the input order with option J in phylip, and it is advisable to do jumbling as least 10 times.")
    parser.add_argument("--out", type=str, required=False, default='out', \
                        help="prefix of output")
    args = parser.parse_args()

    with open(args.out+'.logfile','w') as fw:
        fw.write(time.strftime("%Y-%m-%d %X",time.localtime())+'\n')
        fw.write(os.getcwd()+'\n\n')
        fw.write('python '+sys.argv[0]+'\n')
        fw.write('       --path '+args.path+'\n')
        fw.write('       --ref '+args.ref+'\n')
        fw.write('       --k '+str(args.k)+'\n')
        fw.write('       --order '+args.order+'\n')
        #fw.write('       --p '+args.p+'\n')
        fw.write('       --out '+args.out+'\n')

    ref_data = pandas.read_csv(args.ref,sep='\s+',header=None)
    with open(args.out+'.consensus.filelist','w') as fw:
        pass
    with open(args.out+'.conflict.filelist','w') as fw:
        pass
    ancmerge = output_ancestry(args.path, args.k, ref_data, args.order, args.out)

    '''
    if args.p == 'T':
        anclist = pandas.read_csv(args.out+'.consensus.filelist',header=None)
        anclist = list(anclist[0])
        calibrate_merge_P_file(anclist, ancmerge, args.k, args.out)
    else:
        pass
    '''

    print('Done.\nHave a Nice Day!')

if __name__ == '__main__':
    main()
