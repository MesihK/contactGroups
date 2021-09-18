import commp as cp
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from sklearn.preprocessing import normalize
from sklearn.metrics import roc_auc_score
from sklearn import metrics
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
from utils_pfammsa import pfammsa

colorscheme1 = ['#a93a28', '#afc8cd', '#266674', '#fb8c32', '#cbc96d',
'#60e6c1', '#d7295e', '#008ed0', '#747474']

def evaldisttrans(args):
    distfile = args[0] # A 6 A 7 0.123 resid
    contacts = args[1]  # 0 1 0.123     index
    threshold = float(args[2])

    # load files
    dist = cp.loadtuplesregex(distfile, '\s{1,}')
    cnt = np.loadtxt(contacts, delimiter=',')

    # clasify distance as contact '1' or no contact '0' based on threshold
    distScore = [ 1 if float(e) <= threshold else 0 for a,b,c,d,e in dist ]

    cntScore = list()
    #I think the first token is the start token, so we ignore it
    for i in range(1,len(cnt[0])):
        for j in range(i+1,len(cnt[0])):
            cntScore.append(cnt[i][j])


    # calculate the ROC curve
    aucScore = roc_auc_score(distScore, cntScore)
    curve = metrics.roc_curve(distScore, cntScore, pos_label=1)
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

def evaldistdca(args):
    mapfile = args[0]  # 6 K 0 K       resid to alignedIndex
    distfile = args[1] # A 6 A 7 0.123 resid
    dcafile = args[2]  # 0 1 0.123     index
    threshold = float(args[3])

    # load files
    m = cp.loadtuplesregex(mapfile, '\s{1,}')
    dist = cp.loadtuplesregex(distfile, '\s{1,}')
    dca = cp.loadtuplesregex(dcafile, '\s{1,}')

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

def evalmultdistdca(args):
    distfile = args[0] # A 6 A 7 0.123 resid
    dcafile = args[1].split(',')  # 0 1 0.123     index
    labels = args[2].split(',')  
    threshold = float(args[3])
    
    colors = ['r','g','b','c','m','y','k']

    # load files
    dist = cp.loadtuplesregex(distfile, '\s{1,}')
    dca = [ cp.loadtuplesregex(i, '\s{1,}') for i in dcafile ]

    # clasify distance as contact '1' or no contact '0' based on threshold
    distScore = [ 1 if float(e) <= threshold else 0 for a,b,c,d,e in dist ]
    dcaScore = list()
    for di in dca:
        dcaScore.append([ float(c) for a,b,c in di ])

    # calculate the ROC curve
    # calculate the ROC area under the curve
    plt.figure(1)
    for i in range(len(dca)):
        aucScore = roc_auc_score(distScore, dcaScore[i])
        curve = metrics.roc_curve(distScore, dcaScore[i], pos_label=1)
        plt.plot(curve[0], curve[1], color=colors[i], label=labels[i]+':'+"{:.4f}".format(aucScore))
        print('AUC'+str(i)+':',"{:.4f}".format(aucScore))
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.legend(loc='best')
    plt.savefig('roc_curve.png')
    #plt.show()

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
    if len(pfm.msalist) < targetcnt:
        cp._info('Target count is smaller than msa length!')
        return
        
    cols = pfm.msareduce(['aa'],0.7,0.7)[1]
    scoremat = pfm.scorebycols('aa',cols)
    cp._info('msa reduce done')

    # return a list of len(rows of x), clusters[i] = 'cluster ID which i belongs'
    msaclusters = cp.hamming_cluster(scoremat, 1-similarity_cutoff)
    cp._info('cluster the msa done')

    target_ids = []
    
    # get clusters
    cluster_member_ids = [ np.where(msaclusters == i)[0].tolist() for i in range(1,max(msaclusters)+1)]
    cp._info('id done')
    # sort them according to cluster size, smallest is at top
    cluster_member_ids.sort(key = len, reverse=True)
    cp._info('sorting done')
    cnt = 0
    # sample until target count meet
    while cnt < targetcnt:
        for c in range(len(cluster_member_ids)):
            if cnt == targetcnt: break
            if len(cluster_member_ids[c]) > 0:
                i = random.sample(cluster_member_ids[c], 1)
                target_ids += i 
                cnt += 1
                cluster_member_ids[c].remove(i[0])
        cluster_member_ids = [ c for c in cluster_member_ids if c != [] ]

    target_cluster = [pfm.msalist[i] for i in target_ids]
    cp._info('sampling done')

    # output target cluster msa
    outmsafile = '%s_%.2f_cluster.fa' % (outprefix, similarity_cutoff)
    with open(outmsafile, 'w') as fout:
        fout.write('\n'.join(['>%s\n%s' % (s[0], s[1]) for s in target_cluster]))
    cp._info('save target cluster msa to %s' % outmsafile)

    # output target cluster score
    outscorefile = '%s_%.2f_cluster.scoremat' % (outprefix, similarity_cutoff)
    np.savetxt(outscorefile, scoremat[target_ids,:])
    cp._info('save target cluster score to %s' % outscorefile)

def clusterweight(args):
    assert len(args) == 1, 'Usage: python utils_mesih.py clusterweight PF0000.score'
    scorefile = args[0]

    score = np.loadtxt(scorefile, delimiter=',')
    linkage_matrix = linkage(score, "single", metric='hamming')

    def calc(linkage):
        n = len(linkage)+1
        score = n
        i = 2*n - 1
        scores = list()

        def _calc(i,score):
            if i < n:
                scores.append([i,score])
                return
            c1 = int(linkage[i-n][0])
            c2 = int(linkage[i-n][1])
            d = linkage[i-n][2]
            div = 2
            #if childs are same level as parent, then increase division
            if c1 > n and d == linkage[c1-n][2]: div += 1
            if c2 > n and d == linkage[c2-n][2]: div += 1
            score = score/div
            #if childs are same level as parent, let them divide, give score*2
            if c1 > n and d == linkage[c1-n][2]: _calc(c1,score*2)
            else: _calc(c1,score)
            if c2 > n and d == linkage[c2-n][2]: _calc(c2,score*2)
            else: _calc(c2,score)

        _calc(i-1,score)
        scores.sort()
        w = [ s for i,s in scores]
        return w

    w = calc(linkage_matrix)
    print(repr(w))


if __name__ == '__main__':
    cp.dispatch(__name__)
