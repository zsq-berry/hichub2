import hicstraw
import pandas as pd
import numpy as np
from functools import reduce
from pybedtools import BedTool
import os
import pyarrow.feather as feather

import leidenalg
from itertools import permutations 
import igraph as ig
from pyranges import PyRanges
from scipy.stats import wilcoxon
import random
random.seed(10)
import multiprocessing

from scipy.stats import nbinom
import statsmodels.api as sm
from scipy.stats import nbinom
from scipy.stats import norm
from pymc import ZeroInflatedNegativeBinomial
import scipy.stats as st
from scipy import special
import pyBigWig
from itertools import product
import argparse
from optparse import OptionParser
import warnings
warnings.filterwarnings("ignore")

def get_chrom_names(hic_file):
    hic_info = hicstraw.HiCFile(hic_file)
    chrom_names = hic_info.getChromosomes()
    chrom_names = pd.DataFrame([t.name for t in chrom_names])[0]
    chrom_names = chrom_names[~chrom_names.str.contains('All|ALL|Y|M|_')]
    return chrom_names
    
def implement_hicstraw(hic, data_type, norm, resolution, chrom):
    result = hicstraw.straw(data_type,norm, hic, chrom, chrom, 'BP', resolution)
    if len(result)>0:
        bin1 = pd.DataFrame([t.binX for t in result])
        bin2 = pd.DataFrame([t.binY for t in result])
        counts = pd.DataFrame([t.counts for t in result])
        df_temp = pd.DataFrame(data={'bin1':bin1[0],'bin2':bin2[0], 'score':counts[0]})
        return df_temp
    else: return []

def get_hic_by_chromosome(chrom):
    hic_mtx = []
    for i in range(len(hic)):
        mtx = implement_hicstraw(hic[i], data_type, norm, resolution, chrom)
        mtx.columns = ['bin1','bin2',name[i]]
        if 'chr' in chrom: seqname = chrom
        if 'chr' not in chrom: seqname = 'chr'+chrom
        mtx['#chr'] = seqname
        mtx = mtx.reindex(columns=['#chr','bin1','bin2',name[i]])
        hic_mtx.append(mtx)
    hic_merge = reduce(lambda  left,right: 
                       pd.merge(left,right,on=['#chr','bin1','bin2'],how='outer'), hic_mtx).fillna(0)
    hic_merge['distance'] = (hic_merge['bin2'] - hic_merge['bin1'])/10**4
#     hic_merge = distance_correction(hic_merge)
    hic_merge = hic_merge[(hic_merge['distance']<=200)]
    print(seqname + ' finished')
    return hic_merge

def get_normed_hic(chrom):
    data_type = ['observed','oe']
    mtx1 = implement_hicstraw(hic[i], data_type[0], norm, resolution, chrom)
    mtx2 = implement_hicstraw(hic[i], data_type[1], norm, resolution, chrom)
    mtx1.columns = ['bin1','bin2',name[i]+'-'+data_type[0]]
    mtx2.columns = ['bin1','bin2',name[i]+'-'+data_type[1]]
    mtx1[name[i]+'-'+data_type[0]] = np.round(mtx1[name[i]+'-'+data_type[0]],4)
    mtx2[name[i]+'-'+data_type[1]] = np.round(mtx2[name[i]+'-'+data_type[1]],4)
    if 'chr' in chrom: seqname = chrom
    if 'chr' not in chrom: seqname = 'chr'+chrom
    mtx1['#chr'] = seqname
    mtx2['#chr'] = seqname
    mtx = mtx1.merge(mtx2,on=['#chr','bin1','bin2'])
    mtx = mtx.reindex(columns=['#chr','bin1','bin2',name[i]+'-'+data_type[0],name[i]+'-'+data_type[1]])
    mtx['distance'] = (mtx['bin2'] - mtx['bin1'])/10**4
    mtx = mtx[(mtx['distance']<=200)]
#     mtx.to_csv(path+'hic_data_'+ name[i] +'.txt',mode='a',sep='\t', index=None)
    print(seqname + ' finished')
    return mtx



def Stitch_Region(region,resolution,gap_size):
    former = region[1][range(len(region[1])-1)].values
    latter = region[1][range(1,len(region[1]))].values
    gap = np.where((latter - former)/resolution>gap_size)[0]
    region_stitch = pd.DataFrame({'chr':region[0][0],'start':region[1][np.append(0,(gap+1))].values,
                                  'end':region[1][np.append(gap+1,len(region))-1].values})#+resolution
    return region_stitch

def PageRank_Filter(hic_mtx,col,ratio):
    graph = ig.Graph.DataFrame(hic_mtx,directed=None,use_vids=False)

    graph.vs['pagerank'] = graph.pagerank(weights=graph.es[col])
    global_median = np.percentile(graph.vs['pagerank'], ratio)

    idx = np.where(graph.vs['pagerank']>global_median)[0]
    bin0 = pd.DataFrame(graph.vs['name'])[0][idx]
    
    hic_mtx = hic_mtx[(hic_mtx['reg1'].isin(bin0)) & (hic_mtx['reg2'].isin(bin0))]
    return hic_mtx

def Permutation(name):
    perm = pd.DataFrame(permutations(range(len(name)),2))
    diff = pd.DataFrame()
    for i in range(len(perm)):
        temp = pd.DataFrame(hic_score[name[perm.iloc[i,0]]]-
                            hic_score[name[perm.iloc[i,1]]])
        temp.columns = [name[perm.iloc[i,0]]+'-'+name[perm.iloc[i,1]]]
        diff = pd.concat([diff,temp],axis=1)
    return(diff)
    
    
def Community_Detection(hic_all,weight,res):
    hic_chr = hic_all[['reg1','reg2',weight]]
#     hic_chr = PageRank_Filter(hic_chr,weight,50)
    graph = ig.Graph.DataFrame(hic_chr,directed=None,use_vids=False)
    community = leidenalg.find_partition(graph, leidenalg.CPMVertexPartition,resolution_parameter = res)
    subgraphs = community.subgraphs()
    return(subgraphs)

def community_multilevel(hic_all,weight):
    hic_chr = hic_all[['reg1','reg2',weight]]
    hic_chr = PageRank_Filter(hic_chr,weight,30)
    graph = ig.Graph.DataFrame(hic_chr,directed=None,use_vids=False)
    community = graph.community_multilevel(return_levels=True,weights=graph.es[weight])
    subgraphs = community[0].subgraphs()
    return(subgraphs)

def Island_Capture(subgraphs):
    hub_all = pd.DataFrame()
    for j in range(len(subgraphs)):
        subgraph = subgraphs[j]
        size = len(subgraph.vs)
        edge = len(subgraph.es)
        if size >= 3:
            region = pd.DataFrame([i.split(':', 1) for i in subgraph.vs['name']])
            region[1] = region[1].astype(int)
            hub = Stitch_Region(region,resolution,gap_size)
            hub = hub[hub['end']-hub['start']>=1*10**4]
            if len(hub)>1:
                perm = pd.DataFrame(permutations(list(range(len(hub))),2))
                hub0 = pd.concat([hub.iloc[perm[0].values,:].reset_index(drop=True),
                           hub.iloc[perm[1].values,:].reset_index(drop=True)],axis=1)
                hub0 = pd.concat([pd.concat([hub,hub],axis=1),hub0])
            if len(hub)<=1:hub0 = pd.concat([hub,hub],axis=1)
            hub0.columns = ['chr1','start1','end1','chr2','start2','end2']
            hub0 = hub0[hub0['start1']<=hub0['start2']]
            hub_all = pd.concat([hub_all,hub0])
        else: break
    return(hub_all.reset_index(drop=True))

def pagerank_median(hic_all,weight):
    hic_chr = hic_all[['reg1','reg2',weight]]
    graph = ig.Graph.DataFrame(hic_chr,directed=None,use_vids=False)
    graph.vs['pagerank'] = graph.pagerank(weights=graph.es[weight])
    global_median = np.percentile(graph.vs['pagerank'], 50)
    return global_median

def run_community_detection(i):
    subgraphs = Community_Detection(hic_select,name[i],res_community)
#     subgraphs = community_multilevel(hic_select,pair[i],name[int(i/(len(name)-1))])
    print(name[i] + ' finished with ' + str(len(subgraphs)) + ' subgraphs')
    return subgraphs

def run_island_capture(i):
    cisland = Island_Capture(subgraphs[i])
    cisland.to_csv(path + 'cisland/' +pair[i]+'_pair_hub.bed',index=None,sep='\t')
    print(pair[i] + ' finished with ' + str(len(cisland)) + ' hichubs')

def set_resolution(num_pixels):
    thres = 5*10**7
    res_lim = 0.9
    if num_pixels >= thres:
        return res_lim
    if num_pixels < thres:
        return res_lim - 0.1*(thres - num_pixels)/10**7



def significant_score(hic_distance):
    score = hic_distance['score']
    chrom = hic_distance['#chr'].unique()[0]
    distance = hic_distance['distance'].unique()[0]
    genome_size = pd.read_csv('/home/szhu/annotation/hg38.genome',header=None,sep='\t')
    size = int(genome_size[1][genome_size[0]==chrom]/10**4)
    score0 = np.concatenate((score,np.repeat(0,int(size)- distance - len(score))))
    psi = len(score)/len(score0)
    score = score0[score0<np.percentile(score,95)]

    length = max(score0)+1
    x = np.arange(0, length, 1)
    mu=score.mean()
    sigma=score.std()
    n=mu**2/(sigma**2-mu)
    p=mu/sigma**2
    y = ZeroInfNegBinom(n, mu, psi, x)

    s = []; s0 = 0
    significance = pd.DataFrame({'score':x})
    significance['sig_score'] = 0
    for j in range(0,int(length)):
        s0 = s0 + y[j]
        significance.loc[j,'sig_score'] = -np.log10(1-s0)
    if (len(significance['sig_score'].unique())>2):
        significance.loc[significance['sig_score'].isna(),'sig_score'] = significance['sig_score'].max()
        significance.loc[significance['sig_score']==np.inf,'sig_score'] = np.sort(significance['sig_score'].unique())[-2]
    if (len(significance['sig_score'].unique())<=2):
        significance.loc[significance['sig_score']==np.inf,'sig_score'] = 0
    hic_sig = hic_distance.merge(significance)
    return(hic_sig)

def extend_hichub(cisland_chr):
    hic_pixel = pd.DataFrame()
    for i in range(len(cisland_chr)):
        left = (cisland_chr.loc[i,'end1'] - cisland_chr.loc[i,'start1'])/resolution
        right = (cisland_chr.loc[i,'end2'] - cisland_chr.loc[i,'start2'])/resolution
        left_bin = cisland_chr.loc[i,'start1'] + list(range(0,(int(left))*resolution,resolution))
        right_bin = cisland_chr.loc[i,'start2'] + list(range(0,(int(right))*resolution,resolution))
        temp = pd.DataFrame(product(left_bin,right_bin))
        temp['idx'] = i
        hic_pixel = pd.concat([hic_pixel,temp])
    hic_pixel['#chr'] = cisland_chr.iloc[0,0]
    hic_pixel.columns = ['bin1','bin2','idx','#chr']
    return(hic_pixel)

def ZeroInfNegBinom(a, m, psi, x):
    pmf = special.binom(x + a - 1, x) * (a / (m + a))**a * (m / (m + a))**x
    pmf[0] = (1 - psi) + pmf[0]
    pmf[1:] =  psi * pmf[1:]
    pmf /= pmf.sum()
    return pmf

def significant_hic(hic_distance):
    score = hic_distance['score']
    chrom = hic_distance['#chr'].unique()[0]
    distance = hic_distance['distance'].unique()[0]
    genome_size = pd.read_csv('/home/szhu/annotation/hg38.genome',header=None,sep='\t')
    size = int(genome_size[1][genome_size[0]==chrom]/10**4)
    score0 = np.concatenate((score,np.repeat(0,int(size)- distance - len(score))))
    psi = len(score)/len(score0)
#     score = score0[score0<np.percentile(score,99.5)]

    length = max(score)
    j=max(score)
    x = np.arange(0, length, 1)
    mu=score.mean()
    sigma=score.std()
    n=mu**2/(sigma**2-mu)
    p=mu/sigma**2
#         y=nbinom.pmf(x, n, p)
    y = ZeroInfNegBinom(n, mu, psi, x)

    pval = 0.95;s = 0
    for j in range(0,len(x)):
        s = s + y[j]
        if(s >= pval):break
    hic_cut = hic_distance[hic_distance['score']>=j]
    return hic_cut

def find_hichub(chrom):
    hic_score = implement_hicstraw(hic,'observed', Norm, resolution, chrom)

    if len(hic_score)>0:
        hic_score['distance'] = (hic_score['bin2'] - hic_score['bin1'])/10**4
        hic_score['#chr'] = chrom
        hic_score = hic_score[(hic_score['distance']<=300)]
        hic_sig = hic_score.groupby('distance').apply(significant_score).droplevel('distance').reset_index()
        hic_select = hic_sig[hic_sig['sig_score']> -np.log10(0.01)]
        hic_select = hic_select.reset_index()

        hic_select['reg1'] = hic_select['#chr']+':'+hic_select['bin1'].astype(str)
        hic_select['reg2'] = hic_select['#chr']+':'+hic_select['bin2'].astype(str)

        subgraphs = Community_Detection(hic_select,'score',res_community)
        print(chrom + ' finished with ' + str(len(subgraphs)) + ' subgraphs')
        cisland_chr = Island_Capture(subgraphs)
        
        hic_pixel = extend_hichub(cisland_chr)
        hic_pixel = hic_pixel.merge(hic_sig,how='left').fillna(0)
        cisland_chr['sig_score'] = hic_pixel.groupby('idx')['sig_score'].mean()
        
        print(chrom + ' finished with ' + str(len(cisland_chr)) + ' hichubs')
    return cisland_chr

desc="Single state hichub"

parser = OptionParser(description=desc)

parser.add_option("-p", "--path", action="store", type="string",
                dest="path", help="Input/output path", metavar="<file>")

parser.add_option("-n", "--name", action="store", type="string",
                dest="name", help="Hic file name", metavar="<file>")

parser.add_option("-r", "--resolution", action="store", type="int",
                dest="res", help="Resolution for the output hic file.", metavar="<int>")


(opt, args) = parser.parse_args()
path =  opt.path
name = opt.name
resolution = opt.res

print("Signle state hichub for: %s" % name)
print("Resolution: %s" % resolution)

data_type = 'oe'
Norm = 'NONE'
gap_size = 2
res_community = 0.45
hic = path + '/' + name +'.hic'
chrom_names = get_chrom_names(hic)

import click
@click.group()
@click.command()

def single():

    pool = multiprocessing.Pool()
    pool = multiprocessing.Pool(processes=12)
    outputs = pool.map(find_hichub, chrom_names)
    pool.terminate()
    pool.join()
    hichub = pd.concat(outputs)
    hichub.to_csv(path + '/' + name +'_hichub.bed',index=None,sep='\t')

    if __name__ == "__main__":
        single()
