import os
import argparse
import multiprocessing
import random
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from functools import reduce
from itertools import permutations, combinations, product
from scipy.stats import wilcoxon, false_discovery_control
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import hicstraw
import igraph as ig
import leidenalg
import qnorm
from pybedtools import BedTool 

warnings.filterwarnings("ignore")

class HiCHub2:
    def __init__(self, args):
        self.args = args
        self.path = getattr(args, 'path', None)
        self.name = getattr(args, 'name', None)
        self.resolution = getattr(args, 'resolution', 10000)
        self.Norm = getattr(args, 'norm', 'NONE')
        self.data_type = getattr(args, 'data_type', 'oe')
        self.distance = getattr(args, 'distance', 200)
        self.score_thres = getattr(args, 'score_thres', 1.5)
        self.res_community = getattr(args, 'res_community', 0.6)
        self.gap_size = getattr(args, 'gap_size', 2)
        self.pval = getattr(args, 'pval', 0.05)
        self.qval = getattr(args, 'qval', 0.001)
        random.seed(10)

    # ----------- Original Hi-C Utility Functions -----------
    def get_chrom_names(self, hic_file):
        hic_info = hicstraw.HiCFile(hic_file)
        chrom_names = [t.name for t in hic_info.getChromosomes()]
        chrom_names = pd.Series(chrom_names)
        chrom_names = chrom_names[~chrom_names.str.contains('All|ALL|Y|M|_')]
        return chrom_names

    def implement_hicstraw(self, hic, data_type, norm, resolution, chrom):
        try:
            result = hicstraw.straw(data_type, norm, hic, chrom, chrom, 'BP', resolution)
            if result and len(result) > 0:
                bin1 = pd.DataFrame([t.binX for t in result])
                bin2 = pd.DataFrame([t.binY for t in result])
                counts = pd.DataFrame([t.counts for t in result])
                df_temp = pd.DataFrame(data={'bin1': bin1[0], 'bin2': bin2[0], 'score': counts[0]})
                return df_temp
        except: pass
        return pd.DataFrame()

    def merge_multihic(self, path, name, chrom, dt):
        hic_score0 = []
        for i in name:
            hic_temp = self.implement_hicstraw(os.path.join(path, i + '.hic'), dt, self.Norm, self.resolution, chrom)
            if not hic_temp.empty:
                hic_temp.columns = ['bin1', 'bin2', i]
                hic_score0.append(hic_temp)
        if not hic_score0: return pd.DataFrame()
        hic_score = reduce(lambda left, right: pd.merge(left, right, on=['bin2', 'bin1'], how='outer'), hic_score0).fillna(0)
        hic_score.columns = ['bin1', 'bin2'] + name
        hic_score['#chr'] = chrom
        hic_score['distance'] = (hic_score['bin2'] - hic_score['bin1']) / self.resolution
        return hic_score[hic_score['distance'] <= self.distance]

    def quantile_dist(self, hic_dist):
        hic_dist[self.name] = qnorm.quantile_normalize(hic_dist[self.name])
        return hic_dist

    def Community_Detection(self, hic_all, weight, res):
        hic_chr = hic_all[['reg1', 'reg2', weight]]
        graph = ig.Graph.DataFrame(hic_chr, directed=None, use_vids=False)
        community = leidenalg.find_partition(graph, leidenalg.CPMVertexPartition, resolution_parameter=res)
        return community.subgraphs()

    def Stitch_Region(self, region, resolution, gap_size):
        region = region.sort_values(1)
        former = region[1][range(len(region[1]) - 1)].values
        latter = region[1][range(1, len(region[1]))].values
        gap = np.where((latter - former) / resolution > gap_size)[0]
        region_stitch = pd.DataFrame({
            'chr': region[0].iloc[0],
            'start': region[1][np.append(0, (gap + 1))].values,
            'end': region[1][np.append(gap + 1, len(region)) - 1].values
        })
        return region_stitch

    def Island_Capture(self, subgraphs):
        hub_all = pd.DataFrame()
        for j in range(len(subgraphs)):
            subgraph = subgraphs[j]
            if len(subgraph.vs) >= 3:
                region = pd.DataFrame([i.split(':', 1) for i in subgraph.vs['name']])
                region[1] = region[1].astype(int)
                hub = self.Stitch_Region(region, self.resolution, self.gap_size)
                hub = hub[hub['end'] - hub['start'] >= 10000]
                if len(hub) > 1:
                    perm = pd.DataFrame(permutations(list(range(len(hub))), 2))
                    hub0 = pd.concat([hub.iloc[perm[0].values, :].reset_index(drop=True),
                                      hub.iloc[perm[1].values, :].reset_index(drop=True)], axis=1)
                    hub0 = pd.concat([pd.concat([hub, hub], axis=1), hub0])
                else: hub0 = pd.concat([hub, hub], axis=1)
                hub0.columns = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']
                hub0 = hub0[hub0['start1'] <= hub0['start2']]
                hub_all = pd.concat([hub_all, hub0])
        return hub_all.reset_index(drop=True)

    def extend_hichub(self, hub_temp):
        start1 = int(hub_temp.iloc[0, 1])
        end1   = int(hub_temp.iloc[0, 2])
        start2 = int(hub_temp.iloc[0, 4])
        end2   = int(hub_temp.iloc[0, 5])
        chrom  = hub_temp.iloc[0, 0]

        left_bin = list(range(start1, end1, self.resolution))
        right_bin = list(range(start2, end2, self.resolution))
        
        pixel = pd.DataFrame(product(left_bin, right_bin))
        pixel.columns = ['bin1', 'bin2']
        pixel['#chr'] = chrom
        pixel = pixel[pixel['bin1'] <= pixel['bin2']]
        return pixel

    def wilcoxon_test(self, hichub_score):
        d = hichub_score.iloc[:, 0] - hichub_score.iloc[:, 1]
        if sum(d == 0) == len(d): pvalue = 1
        else: pvalue = wilcoxon(d, alternative='greater').pvalue
        return pvalue

    # ----------- Merge Hub Specific Functions -----------
    def merge_hub(self, hub_concat):
        hub1 = hub_concat.iloc[:, [0, 1, 2]]
        hub2 = hub_concat.iloc[:, [3, 4, 5]]
        hub1['idx'] = list(range(len(hub1)))
        hub2['idx'] = list(range(len(hub2)))

        hub1_bedtools = BedTool.from_dataframe(hub1)
        hub2_bedtools = BedTool.from_dataframe(hub2)
        hub1_overlap = hub1_bedtools.intersect(hub1_bedtools, wao=True).to_dataframe()
        hub2_overlap = hub2_bedtools.intersect(hub2_bedtools, wao=True).to_dataframe()

        hub_intersect = hub1_overlap.merge(hub2_overlap, on=['name', 'thickEnd'])
        network = ig.Graph.DataFrame(hub_intersect[['name', 'thickEnd']])
        membership = pd.DataFrame(network.components().membership)
        hub_intersect['group'] = membership.iloc[hub_intersect['name'].values, 0].values

        region = []
        region.append(hub_intersect.iloc[:, [0, 1, 2, 16]])
        region.append(hub_intersect.iloc[:, [4, 5, 6, 16]])
        region.append(hub_intersect.iloc[:, [9, 10, 11, 16]])
        region.append(hub_intersect.iloc[:, [12, 13, 14, 16]])
        for i in range(2): region[i].columns = ['chr1', 'start1', 'end1', 'group1']
        for i in range(2, 4): region[i].columns = ['chr2', 'start2', 'end2', 'group']
        hub_net = pd.concat([pd.concat([region[0], region[1]]), pd.concat([region[2], region[3]])], axis=1)
        hub_unique = hub_net.groupby('group').agg({
            'chr1': 'first', 'start1': 'min', 'end1': 'max',
            'chr2': 'first', 'start2': 'min', 'end2': 'max'})
        return hub_unique

    def iterate_merge_hub(self, hub_concat):
        hub_unique = self.merge_hub(hub_concat)
        for i in range(5): hub_unique = self.merge_hub(hub_unique)
        hub_unique = hub_unique.sort_values(['chr1', 'start1'])
        return hub_unique

    def get_hichub_values(self, hub_to_process, chrom):
        hic_score = self.merge_multihic(self.path, self.name, chrom, "oe")
        if hic_score.empty: return pd.DataFrame()
        hic_score = hic_score.sort_values(['bin2', 'bin1'])

        hichub_chr = hub_to_process[hub_to_process['chr1'] == chrom].copy()
        if hichub_chr.empty: return pd.DataFrame()
        
        hichub_chr['idx'] = range(len(hichub_chr))
        hichub_chr = hichub_chr.reset_index(drop=True)
        hichub_pixel = hichub_chr.groupby('idx').apply(self.extend_hichub).reset_index().drop('level_1', axis=1)
        hichub_pixel = hichub_pixel.merge(hic_score, how='left').fillna(0)
        hichub_chr = pd.concat([hichub_chr, hichub_pixel.groupby('idx')[self.name].mean()], axis=1)
        print('Finished ' + chrom)
        return hichub_chr

    def kmeans_cluster_hichub(self, hichub, n, plotheatmap_flag=0):
        hichub['idx'] = range(len(hichub))
        matrix_mean = hichub[self.name]
        matrix_scaled = pd.DataFrame(StandardScaler().fit_transform(np.transpose(matrix_mean)))
        matrix_scaled = np.transpose(matrix_scaled)
        matrix_scaled.columns = self.name
        kmeans = KMeans(n_clusters=n, random_state=0).fit(matrix_scaled)

        hichub['cluster'] = kmeans.labels_ + 1
        matrix_scaled['cluster'] = kmeans.labels_ + 1

        hichub = hichub.sort_values('cluster').reset_index(drop=True)
        matrix_scaled = matrix_scaled.sort_values('cluster')

        if plotheatmap_flag == 1:
            plt.figure(figsize=(9, 12))
            sns.set(font_scale=1.8)
            sns.heatmap(matrix_scaled[self.name], cmap="bwr", vmin=-2, vmax=2.3, xticklabels=True, yticklabels=False)
            plt.savefig(getattr(self.args, 'plot_out', 'heatmap.png'))
            print(f"Heatmap saved to {getattr(self.args, 'plot_out', 'heatmap.png')}")
        return hichub

    # ----------- Command Core Logics -----------
    def find_hichub(self, chrom):
        hic = os.path.join(self.path, f"{self.name[0]}.hic")
        hic_score = self.implement_hicstraw(hic, self.data_type, self.Norm, self.resolution, chrom)
        if not hic_score.empty:
            hic_score['distance'] = (hic_score['bin2'] - hic_score['bin1']) / self.resolution
            hic_score['#chr'] = chrom
            hic_score = hic_score[(hic_score['distance'] <= self.distance)]
            hic_select = hic_score[hic_score['score'] > self.score_thres].copy()
            if hic_select.empty: return pd.DataFrame()
            hic_select['reg1'] = chrom + ':' + hic_select['bin1'].astype(str)
            hic_select['reg2'] = chrom + ':' + hic_select['bin2'].astype(str)
            subgraphs = self.Community_Detection(hic_select, 'score', self.res_community)
            hichub_chr = self.Island_Capture(subgraphs)
            if hichub_chr.empty: return pd.DataFrame()
            hichub_chr['idx'] = range(len(hichub_chr))
            hichub_pixel = hichub_chr.groupby('idx').apply(self.extend_hichub).reset_index().drop('level_1', axis=1)
            hichub_pixel = hichub_pixel.merge(hic_score, how='left').fillna(0)
            hichub_chr['score'] = hichub_pixel.groupby('idx')['score'].mean().values
            print(chrom + ' finished with ' + str(sum(hichub_chr['score'] > 1)) + ' hichubs')
            return hichub_chr[hichub_chr['score'] > 1]
        return pd.DataFrame()

    def find_diffhub(self, chrom):
        hic_score = self.merge_multihic(self.path, self.name, chrom, self.data_type)
        if hic_score.empty: return pd.DataFrame()
        hic_score = hic_score.groupby('distance').apply(self.quantile_dist)
        hic_score = hic_score.sort_values(['bin2','bin1'])
        hichub = pd.DataFrame()
        for i, j in [[0, 1], [1, 0]]:
            s1, s2 = self.name[i], self.name[j]
            hic_score['max'], hic_score['diff'] = hic_score[self.name].max(axis=1), hic_score[s1] - hic_score[s2]
            sel = hic_score[(hic_score['max'] >= self.score_thres) & (hic_score['diff'] > 0)].copy()
            if sel.empty: continue
            sel['reg1'], sel['reg2'] = chrom + ':' + sel['bin1'].astype(str), chrom + ':' + sel['bin2'].astype(str)
            subgraphs = self.Community_Detection(sel, s1, self.res_community)
            sub = self.Island_Capture(subgraphs)
            if sub.empty: continue
            sub['idx'] = range(len(sub))
            px = sub.groupby('idx').apply(self.extend_hichub).reset_index().drop('level_1', axis=1).merge(hic_score, how='left').fillna(0)
            sub['pval'] = px.groupby('idx')[[s1, s2]].apply(self.wilcoxon_test).values
            sub['qval'] = false_discovery_control(sub['pval'])
            sub = sub[sub['qval'] < self.qval].copy()
            sub['pairs'] = f"{s1}-{s2}"
            hichub = pd.concat([hichub, sub])
            print(chrom  +' ' + self.name[i] +'-' + self.name[j] + ' finished with ' + str(len(hichub)) + ' hichubs')
        return hichub

    def find_var_hichub(self, chrom):
        hic_score = self.merge_multihic(self.path, self.name, chrom, self.data_type)
        if hic_score.empty: return pd.DataFrame()
        hic_score = hic_score.groupby('distance').apply(self.quantile_dist)
        hic_score = hic_score.sort_values(['bin2','bin1'])
        hic_score['max'], hic_score['var'] = hic_score[self.name].max(axis=1), hic_score[self.name].var(axis=1)
        sel = hic_score[(hic_score['max'] > self.score_thres) & (hic_score['var'] >= 0.1)].copy()
        if sel.empty: return pd.DataFrame()
        # Corrected assignment syntax error here
        sel['reg1'] = chrom + ':' + sel['bin1'].astype(str)
        sel['reg2'] = chrom + ':' + sel['bin2'].astype(str)
        subgraphs = self.Community_Detection(sel, 'max', self.res_community)
        hubs = self.Island_Capture(subgraphs)
        if hubs.empty: return pd.DataFrame()
        hubs['idx'] = range(len(hubs))
        px = hubs.groupby('idx').apply(self.extend_hichub).reset_index().drop('level_1', axis=1).merge(hic_score, how='left').fillna(0)
        hubs = pd.concat([hubs, px.groupby('idx')[self.name].mean().reset_index(drop=True)], axis=1)
        perm = list(permutations(range(len(self.name)), 2))
        p_list = []
        for i, j in perm:
            p = px.groupby('idx')[[self.name[i], self.name[j]]].apply(self.wilcoxon_test).values
            hubs[f'pval-{self.name[i]}-{self.name[j]}'] = p
            p_list.append(p)
        hub_out = hubs[np.min(p_list, axis=0) <= self.pval]
        print(chrom + ' finished with ' + str(len(hub_out)) + ' hichub')
        return hub_out

# ----------- Main CLI -----------

def global_hub_worker(engine_instance, hub_data, chrom):
    """
    Global bridge function that can be pickled by multiprocessing.
    """
    return engine_instance.get_hichub_values(hub_data, chrom)

def main():
    parser = argparse.ArgumentParser(description="HiCHub2 CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    base = argparse.ArgumentParser(add_help=False)
    base.add_argument("--path", required=True)
    base.add_argument("--name", nargs='+', required=True)
    base.add_argument("--resolution", type=int, default=10000)
    base.add_argument("--norm", default="NONE")
    base.add_argument("--data_type", default="oe")
    base.add_argument("--distance", type=int, default=200)
    base.add_argument("--score_thres", type=float, default=2)
    base.add_argument("--res_community", type=float, default=0.6)
    base.add_argument("--gap_size", type=int, default=2)
    base.add_argument("--threads", type=int, default=23)
    base.add_argument("--output", default="hichub_results.bed")

    subparsers.add_parser("singlehub", parents=[base])
    subparsers.add_parser("pairhub", parents=[base]).add_argument("--qval", type=float, default=0.001)
    subparsers.add_parser("varhub", parents=[base]).add_argument("--pval", type=float, default=0.05)

    clust = subparsers.add_parser("cluster_hubs")
    clust.add_argument("--input", required=True)
    clust.add_argument("--name", nargs='+', required=True)
    clust.add_argument("--n_clusters", type=int, default=6)
    clust.add_argument("--plot_out", default="hichub_heatmap.png")
    clust.add_argument("--output", default="hichub_cluster.bed")

    pair_clust = subparsers.add_parser("cluster_pairhubs", parents=[base])
    pair_clust.add_argument("--input_beds", nargs='+', required=True)
    pair_clust.add_argument("--n_clusters", type=int, default=6)
    pair_clust.add_argument("--plot_out", default="hichub_pair_heatmap.png")

    args = parser.parse_args()
    engine = HiCHub2(args)

    if args.command == "cluster_hubs":
        hichub_df = pd.read_csv(args.input, sep='\t')
        clustered = engine.kmeans_cluster_hichub(hichub_df, args.n_clusters, 1)
        clustered.to_csv(args.output, sep='\t', index=False)

    elif args.command == "cluster_pairhubs":
            print('Loading and merging pairwise hichub beds...')
            hichub_list = []
            for f in args.input_beds:
                full_path = os.path.join(args.path, f)
                print(f"Reading: {full_path}")
                hichub_list.append(pd.read_csv(full_path, sep='\t'))
                
            hichub_all = pd.concat(hichub_list)
            hichub_all["idx"] = range(len(hichub_all))

            cisland_unique = engine.iterate_merge_hub(hichub_all)
            cisland_unique['idx'] = range(len(cisland_unique))
            cisland_unique = cisland_unique.reset_index(drop=True)

            print('Calculating values for merged hubs across all samples...')
            chrom_names = engine.get_chrom_names(os.path.join(args.path, args.name[0] + '.hic'))
            
            # Use a list of tuples to pass both the engine and the data
            tasks = [(engine, cisland_unique, chrom) for chrom in chrom_names]

            with multiprocessing.Pool(processes=args.threads) as pool:
                # Use starmap to unpack the tuples into the global worker function
                outputs = pool.starmap(global_hub_worker, tasks)
            
            cisland_value = pd.concat([o for o in outputs if not o.empty])
            cisland_value['max'] = cisland_value[args.name].max(axis=1)
            cisland_value = cisland_value[cisland_value['max'] > 0]

            print('Running K-means clustering on merged hubs...')
            final_cisland = engine.kmeans_cluster_hichub(cisland_value, args.n_clusters, 1)
            final_cisland.to_csv(args.output, sep='\t', index=False)
            print(f"Results saved to {args.output}")

    else:
        chroms = engine.get_chrom_names(os.path.join(args.path, f"{args.name[0]}.hic"))
        task = {"singlehub": engine.find_hichub, "pairhub": engine.find_diffhub, "varhub": engine.find_var_hichub}[args.command]
        with multiprocessing.Pool(args.threads) as pool:
            res_list = [r for r in pool.map(task, chroms) if not r.empty]
        if res_list:
            pd.concat(res_list).to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()