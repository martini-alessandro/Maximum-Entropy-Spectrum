import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import os
from cosmolisa.likelihood import find_redshift
import cosmolisa.cosmology as cs
import h5py

def init_plotting():
    plt.rcParams['figure.figsize'] = (2*3.4, 3.4)
    plt.rcParams['font.size'] = 11
    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.sans-serif'] = ['Bitstream Vera Sans']
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 9
    plt.rcParams['ytick.labelsize'] = 9
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = False
    plt.rcParams['legend.loc'] = 'center left'
    plt.rcParams['axes.linewidth'] = 1
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')

def twod_kde(x,y):
    X, Y = np.mgrid[x.min()*0.9:x.max()*1.1:100j, y.min()*0.9:y.max()*1.1:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([x, y])
    kernel = gaussian_kde(values)
    return X, Y, np.reshape(kernel(positions).T, X.shape)

def FindHeightForLevel(inArr, adLevels):
    """
    compute the height of a two D function
    for a given level
    """
    # flatten the array
    oldshape = np.shape(inArr)
    adInput= np.reshape(inArr,oldshape[0]*oldshape[1])
    # GET ARRAY SPECIFICS
    nLength = np.size(adInput)

    # CREATE REVERSED SORTED LIST
    adTemp = -1.0 * adInput
    adSorted = np.sort(adTemp)
    adSorted = -1.0 * adSorted

    # CREATE NORMALISED CUMULATIVE DISTRIBUTION
    adCum = np.zeros(nLength)
    adCum[0] = adSorted[0]
    for i in range(1,nLength):
        adCum[i] = np.logaddexp(adCum[i-1], adSorted[i])
    adCum = adCum - adCum[-1]

    # FIND VALUE CLOSEST TO LEVELS
    adHeights = []
    for item in adLevels:
        idx=(np.abs(adCum-np.log(item))).argmin()
        adHeights.append(adSorted[idx])

    adHeights = np.array(adHeights)
    return np.sort(adHeights)

def format_90CR(input):
    m,l,h = np.percentile(input,[50,5,95])
    return '& $%.2f_{-%.2f}^{+%.2f}$'%(m,m-l,h-m)
    
if __name__=="__main__":
    events = ['gw150914',
              'gw151012',
              'gw151226',
              'gw170104',
              'gw170608',
              'gw170729',
              'gw170809',
              'gw170814',
              'gw170818',
              'gw170823']
    
#    for e in events:
#        b1 = np.loadtxt('/Users/wdp/Desktop/GWTC1/LALPSD/'+e+'/logB.txt',skiprows=0)
#        b2 = np.loadtxt('/Users/wdp/Desktop/GWTC1/MLCatalog/'+e+'/logB.txt',skiprows=0)
#        print(e.upper()+' '+str(b1)+' '+str(b2))
#    exit()
#    events = ['gw151226']
    import matplotlib.cm as cm
#    from matplotlib.colors  import to_rgba
    import matplotlib.lines as mlines
    
    colors = cm.rainbow(np.linspace(0, 1, len(events)))
    colors = ['royalblue','gold','olive','lime','orange','magenta','cyan','turquoise','teal','purple']
    from gwmodel.utils.utils import McQ2Masses, chi_eff
#    events = ['gw151226']
#    path = "/Users/wdp/repositories/MLGW/paper/img/GWTC-1_results/"
    path = "/Users/wdp/Desktop/GWTC1/LALPSD/"
    init_plotting()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    
    handles = []
    O = cs.CosmologicalParameters(0.68,0.31,0.69,-1.0,0.0)
    for c,e in zip(colors,events):
        try:
            print(os.path.join(path,e+'/Nested_sampler/posterior.dat'))
#            BBH_file = '/Users/wdp/Desktop/GWTC-1_sample_release'+e.upper()+'_GWTC-1.hdf5'
#            BBH = h5py.File(BBH_file, 'r')
            p = np.genfromtxt(os.path.join(path,e+'/Nested_sampler/posterior.dat'),names= True)
            z = np.array([find_redshift(O, np.exp(d)) for d in p['logdistance']])
            m1 ,m2 = McQ2Masses(p['mc']/(1+z),p['q'])
            X,Y,PDF = twod_kde(m1,m2)
            string = e.upper()+format_90CR(m1)+format_90CR(m2)+format_90CR(p['mc']/(1+z))+format_90CR(p['q'])+format_90CR(chi_eff(p['q'], p['spin1'],p['spin2']))
            print(string)
            levs = np.sort(FindHeightForLevel(np.log(PDF),[0.9]))
            ax1.contour(X,Y,np.log(PDF),levs,colors=(c,c,),linewidths=1.5)
            X,Y, PDF =  twod_kde(p['spin1'],p['spin2'])
            levs = np.sort(FindHeightForLevel(np.log(PDF),[0.9]))
            ax2.contour(X,Y,np.log(PDF),levs,colors=(c,c,),linewidths=1.5)
            handles.append(mlines.Line2D([], [], color=c, label=e.upper()))
        except:
            print(e+' posterior not found')
    fig1.legend(bbox_to_anchor=(1.05, 1.05), loc='upper left', borderaxespad=0.,fancybox=True,handles=handles)
    fig2.legend(bbox_to_anchor=(1.05, 1.05), loc='upper left', borderaxespad=0.,fancybox=True,handles=handles)
    ax1.set_xlabel('$m_1/M_\odot$')
    ax1.set_ylabel('$m_2/M_\odot$')
    ax1.set_xlim(0,80)
    ax1.set_ylim(0,50)
    ax2.set_xlabel('$s_1$')
    ax2.set_ylabel('$s_2$')
#    fig1.savefig('/Users/wdp/repositories/MLGW/paper/img/posterior_masses_source.pdf', bbox_inches='tight')
#    fig2.savefig('/Users/wdp/repositories/MLGW/paper/img/spins.pdf', bbox_inches='tight')
    fig1.savefig('/Users/wdp/Desktop/GWTC1/LALPSD/posterior_masses_source.pdf', bbox_inches='tight')
    fig2.savefig('/Users/wdp/Desktop/GWTC1/LALPSD/spins.pdf', bbox_inches='tight')
