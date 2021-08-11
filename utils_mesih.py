import commp as cp
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from sklearn.preprocessing import normalize
from sklearn.metrics import roc_auc_score
from sklearn import metrics
from utils_pfammsa import pfammsa

colorscheme1 = ['#a93a28', '#afc8cd', '#266674', '#fb8c32', '#cbc96d',
'#60e6c1', '#d7295e', '#008ed0', '#747474']

def evaldistdca(args):
    mapfile = args[0]  # 6 K 0 K       resid to alignedIndex
    distfile = args[1] # A 6 A 7 0.123 resid
    dcafile = args[2]  # 0 1 0.123     index
    threshold = float(args[3])

    # load files
    m = cp.loadtuples(mapfile)
    dist = cp.loadtuples(distfile)
    dca = cp.loadtuples(dcafile)

    # clasify distance as contact '1' or no contact '0' based on threshold
    distScore = [ 1 if float(e) <= threshold else 0 for a,b,c,d,e in dist ]
    #dcaScore = normalize([[ c for a,b,c in dca ]], norm='max')[0]
    dcaScore = [ float(d) for r1,r2,m1,m2,i,j,M,d in dca ]
    miScore = [ float(M) for r1,r2,m1,m2,i,j,M,d in dca ]

    # create the index list to convert from 0,1... to aligned index and resid at pdb
    ind = []
    for i in range(len(m)):
        ind.append([int(m[i][2]), int(m[i][0])])

    # check if there are differences in the index sense on distfile and dcafile
    for i in range(len(dist)):
        if ind[int(dca[i][0])][1] != int(dist[i][1]) or ind[int(dca[i][1])][1] != int(dist[i][3]):
            print('Mapping error:',ind[int(dca[i][0])][1],':',
                  ind[int(dca[i][1])][1],'=',dca[i][2],',',dist[i][1],':',
                  dist[i][3],'=',dist[i][4])
            exit()

    # calculate the ROC curve
    aucScore = roc_auc_score(distScore, dcaScore)
    curve = metrics.roc_curve(distScore, dcaScore, pos_label=1)
    plt.figure(1)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.plot(curve[0], curve[1], color=colorscheme1[0], label='AUC:'+"{:.4f}".format(aucScore))
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.legend(loc='best')
    plt.savefig('roc_curve.png')
    #plt.show()

    # calculate the ROC area under the curve
    print('AUC:',"{:.4f}".format(aucScore))

# Cluster given MSA with similaritiy cutoff, and then sample each cluster until given threshold achieved 
def sampleclusters(args):
    assert len(args) == 5, 'Usage: python utils_mesih.py sampleclusters msafile scoretag{aa} similarity_cutoff{0.7} outprefix targetcnt'

    msafile = args[0]
    mapfile = args[1]
    similarity_cutoff = float(args[2])
    outprefix = args[3]
    targetcnt = int(args[4])

    # load msa and get scorebycols
    # msa2score
    pfm = pfammsa(msafile)
    cols = pfm.msareduce(['aa'],0.7,0.7)[1]
    scoremat = pfm.scorebycols('aa',cols)

    # return a list of len(rows of x), clusters[i] = 'cluster ID which i belongs'
    cp._info('cluster the msa')
    msaclusters = cp.hamming_cluster(scoremat, 1-similarity_cutoff)

    target_ids = []
    
    cp._info('sample resulting clusters')
    cluster_member_ids = [[ j for j in range(len(msaclusters)) if msaclusters[j] == i ] for i in range(1,max(msaclusters)+1)]
    cluster_member_ids.sort(key = len, reverse=True)
    #cluster_member_ids.reverse()
    print([ len(c) for c in cluster_member_ids ])
    cnt = 0
    while cnt < targetcnt:
        for c in range(len(cluster_member_ids)):
            if cnt == targetcnt: break
            if len(cluster_member_ids[c]) > 0:
                i = random.sample(cluster_member_ids[c], 1)
                target_ids += i 
                cnt += 1
                cluster_member_ids[c].remove(i[0])
        cluster_member_ids = [ c for c in cluster_member_ids if c != [] ]
    #            print(cluster_member_ids[c],i[0])
    #print(target_ids)
    #print(cluster_member_ids)

    target_cluster = [pfm.msalist[i] for i in target_ids]

    # output target cluster msa
    outmsafile = '%s_%.2f_cluster.fa' % (outprefix, similarity_cutoff)
    with open(outmsafile, 'w') as fout:
        fout.write('\n'.join(['>%s\n%s' % (s[0], s[1]) for s in target_cluster]))
    cp._info('save target cluster msa to %s' % outmsafile)

    # output target cluster score
    outscorefile = '%s_%.2f_cluster.scoremat' % (outprefix, similarity_cutoff)
    np.savetxt(outscorefile, scoremat[target_ids,:])
    cp._info('save target cluster score to %s' % outscorefile)

if __name__ == '__main__':
    cp.dispatch(__name__)
