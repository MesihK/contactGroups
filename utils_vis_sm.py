import sys
import commp as cp
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

def barplot_fgdist(arglist):
    if len(arglist) < 1:
        cp._err('Usage: python utils_vis_sm.py barplot_fgdist data.tsv outfile')
    infile = arglist[0]

    bar_width = 0.3
    opacity = 0.8
    spnum = 2
    title = ''

    data = np.loadtxt(infile, delimiter=' ')
    # be consistent with utils_pfammsa.py::wfreqcsfgdist
    xt = ['%s-%s' % (cp.aas01[i], cp.aas01[j]) for i in xrange(len(cp.aas01)) for j in xrange(i,len(cp.aas01))]
    fig ,ax = plt.subplots(spnum, figsize=(18,6), sharey=True)

    index  = np.arange(105)
    ax[0].bar(index, data[0:105], bar_width,
                alpha=opacity,
                color = '#75a8b9',
                label=title)
    plt.sca(ax[0])
    plt.xticks(index+bar_width/2, xt[0:105], rotation='vertical', fontname = "monospace")
    #plt.ylabel('Number of counts')
    plt.legend()

    ax[1].bar(index, data[105:], bar_width,
                alpha=opacity,
                color = '#75a8b9',
                label=title)
    plt.sca(ax[1])
    plt.xticks(index+bar_width/2, xt[105:], rotation='vertical', fontname = "monospace")
    #plt.ylabel('Number of counts')
    plt.legend()
    plt.title(infile)

    #plt.xlabel('Pair substitution frequency')
    plt.tight_layout()
    #plt.show()
    fig.savefig(infile+'.png')
    cp._info('save to %s.png' % infile)


def barplot_bgdist(arglist):
    if len(arglist) < 2:
        cp._err('Usage: python utils_vis_sm.py barplot_bgdist data.tsv outfile')
    infile = arglist[0]
    outfile = arglist[1]

    data = np.loadtxt(infile,delimiter=' ')
    n = len(data)
    fix, ax = plt.subplots()
    index = np.arange(n)

    bar_width = 0.3
    opacity = 0.4

    rects = plt.bar(index+bar_width / 2, data, bar_width,
                    alpha=opacity,
                    color='#5293a8',
                    label='background distribution')
    #plt.xlabel('Amino acids')
    #plt.ylabel('')
    #plt.title('Marginal ubstritution scores comparison')
    plt.ylim(0,0.2)
    xlabels = cp.aat01
    plt.xticks(index + bar_width / 2, xlabels)

    plt.legend()

    plt.tight_layout()
    #plt.show()
    plt.savefig(outfile)
    cp._info('save to %s' % outfile)

def _augmented_dendrogram(*args, **kwargs):
    ddata=dendrogram(*args, **kwargs)
    if not kwargs.get('no_plot', False):
        for i,d in zip(ddata['icoord'], ddata['dcoord']):
            x=0.5*sum(i[1:3])
            y=d[1]
            plt.plot(x,y,'ro')
            str = '%.3g' % y
            plt.annotate(str, (x,y), xytext=(0,-8), textcoords='offset points', va='top', ha='center')
    return ddata

def dendrogram_fgdist(arglist):
    if len(arglist) < 1:
	    cp._err('Usage: python utils_vis_sm.py dendrogram_fgdist cs3.210dist.vec')
    sys.setrecursionlimit(10000)

    datafile = arglist[0]
    outfile = '%s.dendrogram.png' % datafile
    x = np.loadtxt(datafile, delimiter=' ')

    linkage_matrix = linkage(x, "single", metric='correlation')
    plt.figure(figsize=(24, 15))
    #ddata = dendrogram(linkage_matrix)
    xt = ['%s-%s' % (cp.aas01[i], cp.aas01[j]) for i in xrange(len(cp.aas01)) for j in xrange(i,len(cp.aas01))]
    ddata = dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1, link_color_func=lambda x: "k")
    #ddata = _augmented_dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1)
    #plt.show()	        
    plt.savefig(outfile)
    cp._info('save to %s' % outfile)


def dendrogram_fgdist_1(arglist):
    if len(arglist) < 1:
	    cp._err('Usage: python utils_vis_sm.py dendrogram_fgdist_1 cs3.210dist_1.vec')
    sys.setrecursionlimit(10000)

    datafile = arglist[0]
    outfile = '%s.dendrogram.png' % datafile
    x = np.loadtxt(datafile, delimiter=' ')

    linkage_matrix = linkage(x, "single", metric='correlation')
    plt.figure(figsize=(24, 15))
    #ddata = dendrogram(linkage_matrix)
    xt = cp.aas01
    ddata = dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1)
    #ddata = _augmented_dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1)
    #plt.show()	        
    plt.savefig(outfile)
    cp._info('save to %s' % outfile)


def dendrogram_emboss_sm(arglist):
    if len(arglist) < 1:
	    cp._err('Usage: python utils_vis_sm.py dendrogram_emboss_sm b62.sm.dat')
    sys.setrecursionlimit(10000)

    datafile = arglist[0]
    outfile = '%s.dendrogram.png' % datafile
    x = np.loadtxt(datafile, delimiter=' ')

    linkage_matrix = linkage(x, "single", metric='correlation')
    plt.figure(figsize=(24, 15))
    xt = cp.aa201 # emboss sm order
    ddata = dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1, link_color_func=lambda x: "k")
    #ddata = _augmented_dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1)
    #plt.show()	        
    plt.savefig(outfile)
    cp._info('save to %s' % outfile)



def dendrogram_mat(arglist):
    if len(arglist) < 3:
        cp._err('Usage: python utils_vis_sm.py dendrogram_mat matfile tickfile outfile')
    matfile = arglist[0]
    tickfile = arglist[1]
    outfile = arglist[2]

    mat = np.loadtxt(matfile, delimiter= ' ')
    xt = cp.loadlines(tickfile)

    dists = squareform(mat)
    #linkage_matrix = linkage(dists, "single")
    linkage_matrix = linkage(dists, "complete")
    #plt.figure(figsize=(25, 28))
    plt.figure(figsize=(15, 15))
    ddata = dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1, link_color_func=lambda x: "k")
    plt.tight_layout()
    plt.savefig(outfile)
    #plt.show()

    '''
    dists = squareform(mat)
    linkage_matrix = linkage(dists, "single")
    dendrogram(linkage_matrix, labels=["0", "1", "2"])
    plt.title("test")
    plt.show()
    '''


# plot side-by-side barplot from the data in a two-column file
# total_num: total row number
# group_num: how many barplots to split
def bar_sbs(arglist):
    if len(arglist) < 4:
        cp._err('Usage: python utils_vis_sm.py bar_sbs data.2.vec total_num group_num xtickfile')
    datafile = arglist[0]
    totalnum = int(arglist[1])
    spnum = int(arglist[2])
    n_groups=itv= (totalnum / spnum) +1

    print n_groups

    # x-tick
    xt = cp.loadlines(arglist[3])
    #xt = ['%d' % i for i in range(1,92)]

    data = np.loadtxt(datafile, delimiter=' ')
    d1 = data[:,0]
    d2 = data[:,1]

    fig, ax = plt.subplots(spnum, figsize=(12,8), sharey=True)
    bar_width = 0.25
    opacity = 1.0
    for i in range(0, spnum):
            index = np.arange(n_groups if (i+1)*itv <totalnum else (totalnum - i*itv))
            sxt = xt[i*itv:((i+1)*itv if (i+1)*itv <totalnum else totalnum)]

            ax[i].bar(index, d1[i*itv:((i+1)*itv if (i+1)*itv < totalnum else totalnum)], bar_width,
                            alpha=opacity,
                            color='#75a8b9',
                            label='Ordered region')

            ax[i].bar(index + bar_width, d2[i*itv:((i+1)*itv if (i+1)*itv < totalnum else totalnum)], bar_width,
                            alpha=opacity,
                            color='#8d5543',
                            label='Disordered region')
            '''
            ax[i].bar(index + 2*bar_width, scsc23a[i*itv:((i+1)*itv if (i+1)*itv <210 else 210)], bar_width,
                            alpha=opacity,
                            color='#f3db81',
                            label='SCSC2')
            '''
            plt.sca(ax[i])
            plt.xticks(index + bar_width / 2, sxt, rotation='vertical')

            # plt.xlabel('xxx')
            plt.ylabel('Normalized correlation level')
            plt.legend()

    plt.tight_layout()
    plt.show()

   # fig.savefig(outname)


# generate histogram from data.list
# input: one column (float number) data
def hist(arglist):
    if len(arglist) <2:
        cp._err('Usage: python utils_vis_sm.py hist data.list bin')
    datafile = arglist[0]
    n_bins = int(arglist[1])
    x = np.loadtxt(datafile)
    n = x.shape[0]
    print n
    fig, axs = plt.subplots(tight_layout=True)
    axs.hist(x, bins=n_bins)
    plt.show()

# generate side-by-side histogram from data2.list
# input: two one-column (float number) data files
def hist_sbs(arglist):
    if len(arglist) <2:
        cp._err('Usage: python utils_vis_sm.py hist_sbs data1.list data2.list')
    datafile1 = arglist[0]
    datafile2 = arglist[1]
    data1 = np.loadtxt(datafile1)
    data2 = np.loadtxt(datafile2)
    print data1.shape
    print data2.shape
    x = data1
    y = data2
    n_bins = 50
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    axs[0].hist(x, bins=n_bins)
    axs[1].hist(y, bins=n_bins)
    plt.show()
    #fig.savefig('out.png')
    #cp._info('save to out.png')


# generate simple signal plot (using plot) from two-column (space separted columns) file
# column 1: x-tick
# column 2: values
# originally for visualizing sliding window results of IDP project
def signalplot(arglist):
    if len(arglist) < 2:
        cp._err('Usage: python utils_vis_sm.py singalplot infile')

    infile = arglist[0] 
    outfile = arglist[1]

    valuelist= []
    ticklist = []
    for line in cp.loadlines(infile):
        sarr = line.split(' ')
        ticklist.append(int(sarr[0]))
        valuelist.append(float(sarr[1]))

    print '%d records loaded.' % (len(valuelist))

    fig, ax = plt.subplots(figsize=(14,3),tight_layout=True)
    ax.plot(ticklist, valuelist, label="score")
    #ax.plot(valuelist, label="score")
    ax.grid(True)
    ax.legend()
    plt.axvline(31, linestyle='--', c='r')
    plt.xticks(rotation='vertical', fontname = "monospace")

    #plt.show()
    plt.savefig(outfile)


# batch generate plots for disordered protein
def signalplotws(arglist):
    if len(arglist) < 2:
        cp._err('Usage: python utils_vis_sm.py singalplot_ws infile outfile')

    infile = arglist[0] 
    sarr = infile.split('-')
    dstart = int(sarr[4])
    dend = int(sarr[5])
    outfile = arglist[1]

    valuelist= []
    ticklist = []
    for line in cp.loadlines(infile):
        sarr = line.split(' ')
        ticklist.append(int(sarr[0]))
        valuelist.append(float(sarr[1]))

    #print '%d records loaded.' % (len(valuelist))

    fig, ax = plt.subplots(figsize=(14,3),tight_layout=True)
    ax.plot(ticklist, valuelist, label="score")
    ax.grid(True)
    ax.legend()
    print outfile, dstart, dend
    plt.axvline(dstart, linestyle='--', c='r')
    plt.axvline(dend, linestyle='--', c='r')

    plt.xticks(rotation='vertical', fontname = "monospace")

    plt.show()
    #plt.savefig(outfile)


if __name__ == '__main__':
        cp.dispatch(__name__)
