import commp as cp
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from sklearn.preprocessing import normalize
from sklearn.metrics import roc_auc_score
from sklearn import metrics

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
    dcaScore = [ float(c) for a,b,c in dca ]

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
    plt.show()

    # calculate the ROC area under the curve
    print('AUC:',"{:.4f}".format(aucScore))


if __name__ == '__main__':
    cp.dispatch(__name__)
