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
import warnings
warnings.filterwarnings("ignore")
import multiprocessing
import time

from scipy.stats import nbinom
import statsmodels.api as sm
import numpy as np
from scipy.stats import nbinom
from scipy.stats import norm
from pymc import ZeroInflatedNegativeBinomial
import scipy.stats as st
from scipy import special
import pyBigWig
from itertools import product
from itertools import combinations

import qnorm
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from pybedtools import BedTool
import igraph
import argparse
from optparse import OptionParser

random.seed(10)

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

def Stitch_Region(region,resolution,gap_size):
    former = region[1][range(len(region[1])-1)].values
    latter = region[1][range(1,len(region[1]))].values
    gap = np.where((latter - former)/resolution>gap_size)[0]
    region_stitch = pd.DataFrame({'chr':region[0][0],'start':region[1][np.append(0,(gap+1))].values,
                                  'end':region[1][np.append(gap+1,len(region))-1].values})#+resolution
    return region_stitch

def Permutation(name):
    perm = pd.DataFrame(permutations(range(len(name)),2))
    diff = pd.DataFrame()
    for i in range(len(perm)):
        temp = pd.DataFrame(hic_score[name[perm.iloc[i,0]]]-
                            hic_score[name[perm.iloc[i,1]]])
        temp.columns = [name[perm.iloc[i,0]]+'-'+name[perm.iloc[i,1]]]
        diff = pd.concat([diff,temp],axis=1)
    return(diff)

def PageRank_Filter(hic_mtx,col,ratio):
    graph = ig.Graph.DataFrame(hic_mtx,directed=None,use_vids=False)

    graph.vs['pagerank'] = graph.pagerank(weights=graph.es[col])
    global_median = np.percentile(graph.vs['pagerank'], ratio)

    idx = np.where(graph.vs['pagerank']>global_median)[0]
    bin0 = pd.DataFrame(graph.vs['name'])[0][idx]
    
    hic_mtx = hic_mtx[(hic_mtx['reg1'].isin(bin0)) & (hic_mtx['reg2'].isin(bin0))]
    return hic_mtx
    
def Community_Detection(hic_all,weight,res):
    hic_chr = hic_all[['reg1','reg2',weight]]
#     hic_chr = PageRank_Filter(hic_chr,weight,30)
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
            hub = hub[hub['end']-hub['start']>=3*10**4]
            if len(hub)>1:
                perm = pd.DataFrame(permutations(list(range(len(hub))),2))
                hub0 = pd.concat([hub.iloc[perm[0].values,:].reset_index(drop=True),
                           hub.iloc[perm[1].values,:].reset_index(drop=True)],axis=1)
                hub0 = pd.concat([pd.concat([hub,hub],axis=1),hub0])
            if len(hub)<=1:hub0 = pd.concat([hub,hub],axis=1)
            hub0.columns = ['chr1','start1','end1','chr2','start2','end2']
            hub0 = hub0[hub0['start1']<=hub0['start2']]
            hub0['network_id'] = j
            hub_all = pd.concat([hub_all,hub0])
        else: break
    return(hub_all.reset_index(drop=True))


def extend_hichub(hub_temp):
    left_bin = list(range(int(hub_temp.iloc[:,1]),int(hub_temp.iloc[:,2]),resolution))
    right_bin = list(range(int(hub_temp.iloc[:,4]),int(hub_temp.iloc[:,5]),resolution))
    pixel = pd.DataFrame(product(left_bin,right_bin))
    pixel.columns = ['bin1','bin2']
    pixel['#chr'] = hub_temp.iloc[0,0]
    pixel = pixel[pixel['bin1']<=pixel['bin2']]
    return pixel

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


def get_hub_value_by_chrom(hub_region,chrom):
    if 'chr' in chrom: seqname = chrom
    if 'chr' not in chrom: seqname = 'chr'+chrom
    hub_temp = hub_region[hub_region.iloc[:,0] == seqname]
    hub_temp = hub_temp.reset_index(drop=True)
    for k in range(len(name)):
        hic = hicstraw.HiCFile(path + '/' + name[k] + '.hic')
        mzd = hic.getMatrixZoomData(chrom, chrom, "oe",Norm, "BP", resolution)
        hub_temp[name[k]] = 0
        for i in range(len(hub_temp)):
            hub_box = hub_temp.iloc[i,:]
            gr1 = int(hub_box[1]);gr2 = int(hub_box[2])
            gc1 = int(hub_box[4]);gc2 = int(hub_box[5])
            matrix = mzd.getRecordsAsMatrix(gr1, gr2, gc1, gc2)
            hub_temp.loc[i,name[k]] = np.mean(matrix)
    print('Finished ' + chrom)
    return hub_temp

def get_hub_value(chrom):
    return get_hub_value_by_chrom(cisland_unique,chrom)


def kmeans_cluster_hichub(hichub,n):
    hichub['idx'] = range(len(hichub))
    matrix_mean = hichub[name]
    matrix_scaled = pd.DataFrame(StandardScaler().fit_transform(np.transpose(matrix_mean)))
    matrix_scaled = np.transpose(matrix_scaled)
    matrix_scaled.columns = name
    kmeans = KMeans(n_clusters=n, random_state=0).fit(matrix_scaled)

    matrix_scaled['cluster'] = kmeans.labels_ + 1
    matrix_scaled = matrix_scaled.sort_values('cluster')

    hichub['cluster'] = kmeans.labels_ + 1
    hichub = hichub.sort_values('cluster')
    hichub = hichub.reset_index(drop=True)
    hichub['reorder'] = hichub['cluster']
    return hichub

def reorder_cluster(cisland_cluster):
    matrix_scaled = pd.DataFrame(StandardScaler().fit_transform(np.transpose(cisland_cluster[name])))
    matrix_scaled = np.transpose(matrix_scaled)
    matrix_scaled.columns = name
    matrix_scaled['cluster'] = cisland_cluster['cluster']
    group_mean = matrix_scaled.groupby('cluster')[name].mean()
    order = group_mean.sort_values(name[sort_by],ascending=True).index.tolist()

    for i in range(n_clusters):
        cisland_cluster.loc[cisland_cluster['cluster'] == order[i],'reorder'] = i + 1
    cisland_cluster = cisland_cluster.sort_values('reorder')
    cisland_cluster = cisland_cluster.reset_index(drop=True)
    return cisland_cluster


def plotheatmap(hichub):
    matrix_mean = hichub[name]
    matrix_scaled = pd.DataFrame(StandardScaler().fit_transform(np.transpose(matrix_mean)))
    matrix_scaled = np.transpose(matrix_scaled)
    matrix_scaled.columns = name
    sns.heatmap(matrix_scaled, cmap="RdYlBu_r", vmin = -2, vmax = 2.3,
                             xticklabels=True, yticklabels=False)

def iterate_merge_hub(hub_concat):
    hub_unique = merge_hub(hub_concat)
    for i in range(10):hub_unique = merge_hub(hub_unique)
    hub_unique = hub_unique.sort_values(['chr1','start1'])
    return hub_unique

def quantile_dist(hic_dist):
    hic_dist[name] = qnorm.quantile_normalize(hic_dist[name])
    return hic_dist

def wilcoxon_test(hichub_score):
    d = hichub_score.iloc[:,0] - hichub_score.iloc[:,1]
    if sum(d==0)==len(d):
        pvalue = 1
    if sum(d==0)!=len(d): 
        pvalue = wilcoxon(d).pvalue
    return pvalue

def merge_multihic(path,name,chrom):
    hic_score0 = []
    for i in name:
        hic_score0.append(implement_hicstraw(path + '/' + i +'.hic',data_type, Norm, resolution, chrom))
    hic_score = reduce(lambda  left,right: pd.merge(left,right,on=['bin1','bin2'],how='outer'), 
                       hic_score0).fillna(0)
    hic_score.columns = ['bin1','bin2']+name
    hic_score['#chr'] = chrom
    hic_score['distance'] = (hic_score['bin2'] - hic_score['bin1'])/10**4
    hic_score = hic_score[hic_score['distance']<=dist_cutoff]
    return hic_score


def ZeroInfNegBinom(a, m, psi, x):
    pmf = special.binom(x + a - 1, x) * (a / (m + a))**a * (m / (m + a))**x
    pmf[0] = (1 - psi) + pmf[0]
    pmf[1:] =  psi * pmf[1:]
    pmf /= pmf.sum()
    return pmf

def significant_score(hic_distance):
    score = hic_distance.iloc[:,0]
    score = score[score>0]
    chrom = hic_distance['#chr'].unique()[0]
    distance = hic_distance['distance'].unique()[0]
    genome_size = pd.read_csv('/home/szhu/annotation/hg38.genome',header=None,sep='\t')
    seqname = chrom
    if 'chr' not in chrom: seqname = 'chr' + chrom
    size = int(genome_size[1][genome_size[0]==seqname]/10**4)
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

    s = []; s0 = 0; expr = hic_distance.columns[0]; sig=expr+'_sig'
    significance = pd.DataFrame({expr:x})
    significance[expr+'_sig'] = 0
    for j in range(0,int(length)):
        s0 = s0 + y[j]
        significance.loc[j,sig] = -np.log10(1-s0)
    if (len(significance[sig].unique())>2):
        significance.loc[significance[sig].isna(),sig] = significance[sig].max()
        significance.loc[significance[sig]==np.inf,sig] = np.sort(significance[sig].unique())[-2]
    if (len(significance[sig].unique())<=2):
        significance.loc[significance[sig]==np.inf,sig] = 0
    significance["#chr"]=chrom
    significance["distance"]=distance
    return(significance)


def quantile_dist(hic_dist):
    hic_dist[name1] = qnorm.quantile_normalize(hic_dist[name1])
    return hic_dist

def merge_multihic(path,name,chrom):
    hic_score0 = []
    for i in name:
        hic_mtx1 = implement_hicstraw(path + '/' + i +'.hic',"oe", Norm, resolution, chrom)
        hic_mtx1.columns = ['bin1','bin2',i]
        hic_mtx1['distance'] = (hic_mtx1['bin2'] - hic_mtx1['bin1'])/10**4
        hic_mtx1 = hic_mtx1[hic_mtx1['distance']<=dist_cutoff]
        hic_mtx2 = implement_hicstraw(path + '/' + i +'.hic',"observed", Norm, resolution, chrom)
        hic_mtx2.columns = ['bin1','bin2',i+'_observed']
        hic_mtx2['distance'] = (hic_mtx2['bin2'] - hic_mtx2['bin1'])/10**4
        hic_mtx2 = hic_mtx2[hic_mtx2['distance']<=dist_cutoff]
        hic_score0.append(hic_mtx1)
        hic_score0.append(hic_mtx2)
    hic_score = reduce(lambda  left,right: pd.merge(left,right,how='outer'), hic_score0).fillna(0)
    hic_score['#chr'] = chrom
    return hic_score

def find_var_hichub_ZINB(chrom):
    combn = pd.DataFrame(combinations(range(len(name)),2))

    hic_score = merge_multihic(path,name,chrom)
#     hic_score = hic_score.groupby('distance').apply(quantile_dist)
    hic_score[name1] = hic_score[name1].round(0)

    for i in range(len(name)):
        sig_temp = hic_score[[name[i]+'_observed',"#chr","distance"]].groupby('distance').\
        apply(significant_score).droplevel('distance')
        sig_temp = sig_temp.reset_index().drop(columns=['index'])
        hic_score = hic_score.merge(sig_temp)

    hic_score['max'] = hic_score[name0].max(axis=1)
    hic_score['std'] = hic_score[name].std(axis=1)
    hic_score['mean'] = hic_score[name].mean(axis=1)
#     hic_select = hic_score[(hic_score['max']>score_thres)&
#                            (hic_score['std']/hic_score['mean']>0)]

    for i in range(len(combn)):
        idx = combn.iloc[i,:]
        hic_score['pval'+'-'+name[idx[0]]+'-'+name[idx[1]]] = abs(hic_score[name[idx[0]]] - \
                                                                  hic_score[name[idx[1]]])
    percentile = hic_score.iloc[:,-len(combn):].abs().quantile(q=0.5).max()
    hic_select = hic_score[(hic_score['max']>score_thres)&(hic_score.iloc[:,-len(combn):].max(axis=1)>=percentile)]

    hic_select['reg1'] = hic_select['#chr']+':'+hic_select['bin1'].astype(str)
    hic_select['reg2'] = hic_select['#chr']+':'+hic_select['bin2'].astype(str)

    subgraphs = Community_Detection(hic_select,'max',res_community)
    hichub = Island_Capture(subgraphs)
    hichub['idx'] = range(len(hichub))

    hichub_pixel = hichub.groupby('idx').apply(extend_hichub).reset_index().drop('level_1',axis=1)
    hichub_pixel = hichub_pixel.merge(hic_score,how='left').fillna(0)
    hichub = pd.concat([hichub,hichub_pixel.groupby('idx')[name].mean()],axis=1)

    for i in range(len(combn)):
        idx = combn.iloc[i,:]
        hichub['pval'+'-'+name[idx[0]]+'-'+name[idx[1]]] = \
        hichub_pixel.groupby('idx')[name[idx[0]],name[idx[1]]].apply(wilcoxon_test)
        hichub = hichub[(hichub.iloc[:,-len(combn):].min(axis=1)<=pval)]

    print(chrom + ' finished with ' + str(len(hichub)) + ' hichub')
    return hichub



desc="Multiple states hichub"

parser = OptionParser(description=desc)

parser.add_option("-p", "--path", action="store", type="string",
                dest="path", help="Input/output path", metavar="<file>")

parser.add_option("-n", "--names", action="store", type="string",
                dest="name", help="HiC file names separated by comma", metavar="<file>")

parser.add_option("-N", "--norm", action="store", type="string",
                dest="norm", help="HiC matrix normalization (NONE/KR)", metavar="<file>")

parser.add_option("-c", "--nclusters", action="store", type="int",
                dest="nclusters", help="Number of clusters", metavar="<int>")

parser.add_option("-r", "--resolution", action="store", type="int",
                dest="res", help="Resolution fof hic matrix", metavar="<int>")


(opt, args) = parser.parse_args()
path =  opt.path
name = opt.name.split(",")
Norm = opt.norm
n_clusters = opt.nclusters
resolution = opt.res

print("Multiple states hichub for: %s" % opt.name)
print("Resolution: %s" % resolution)

data_type = 'observed'
gap_size = 2
dist_cutoff = 200
score_thres = -np.log10(0.01)
res_community = 0.45
pval = 10**-5
sort_by = 0

name0 = [n+'_observed_sig' for n in name]
name1 = [n+'_observed' for n in name]

chrom_names = get_chrom_names(path+ '/' + name[0] +'.hic')

import click
@click.group()
@click.command()

def multiple():
    pool = multiprocessing.Pool()
    pool = multiprocessing.Pool(processes=len(chrom_names))
    outputs = pool.map(find_var_hichub_ZINB, chrom_names)
    pool.terminate()
    pool.join()
    hichub = pd.concat(outputs)

    hichub_cluster = kmeans_cluster_hichub(hichub,n_clusters)
    hichub_cluster = reorder_cluster(hichub_cluster)
    hichub_cluster.to_csv(path +'/hichub_cluster.bed' ,sep='\t',index_label=False)



if __name__ == "__main__":
    multiple()
