import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import os
from cosmolisa.likelihood import find_redshift
import cosmolisa.cosmology as cs
import h5py
import matplotlib.font_manager

def init_plotting():
    plt.rcParams['figure.figsize'] = (3*3.4, 3.4)
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
    plt.rcParams['text.usetex'] = True


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
    colors = ['dodgerblue','olivedrab','goldenrod','yellowgreen','peru','magenta','mediumturquoise','aquamarine','skyblue','mediumpurple']
    line_style = ['-','--', '-','-','--','-','-','--','--','--']
    from utils import McQ2Masses, chi_eff, final_mass, final_spin
#    events = ['gw151226']
#    path = "/Users/wdp/repositories/MLGW/paper/img/GWTC-1_results/"
    path = "../img/GWTC-1_results"
    init_plotting()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(121)
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(121)
    #fig3 = plt.figure()
    ax3 = fig1.add_subplot(122)
    ax4 = fig2.add_subplot(122)
    fig1.subplots_adjust(bottom = 0.33, wspace = 0.2)  
    fig2.subplots_adjust(bottom = 0.33, wspace = 0.2)  

    handles = []
    O = cs.CosmologicalParameters(0.68,0.31,0.69,-1.0,0.0)
    for c, ls,e in zip(colors,line_style, events):
#        try:
        if 1:
            print(os.path.join(path,e+'/Nested_sampler/posterior.dat'),c)
#            BBH_file = '/Users/wdp/Desktop/GWTC-1_sample_release'+e.upper()+'_GWTC-1.hdf5'
#            BBH = h5py.File(BBH_file, 'r')
            p = np.genfromtxt(os.path.join(path,e+'/Nested_sampler/posterior.dat'),names= True)
                #ax1 - mass TEOB
            z = np.array([find_redshift(O, np.exp(d)) for d in p['logdistance']])
            m1 ,m2 = McQ2Masses(p['mc']/(1+z),p['q'])
            X,Y,PDF = twod_kde(m1,m2)
            string = e.upper()+format_90CR(m1)+format_90CR(m2)+format_90CR(p['mc']/(1+z))+format_90CR(p['q'])+format_90CR(chi_eff(p['q'], p['spin1'],p['spin2']))
            print(string)
            levs = np.sort(FindHeightForLevel(np.log(PDF),[0.9]))
            ax1.contour(X,Y,np.log(PDF),levs, linestyles = (ls,ls), colors=(c,c,),linewidths=1.5)
                #ax2 - spins TEOB
            X,Y, PDF =  twod_kde(p['spin1'],p['spin2'])
            levs = np.sort(FindHeightForLevel(np.log(PDF),[0.9]))
            ax2.contour(X,Y,np.log(PDF),levs, linestyles = (ls,ls), colors=(c,c,),linewidths=1.5)
                #ax3 - final TEOB
            mass_f = final_mass(m1,m2, p['spin1'],p['spin2'])
            spin_f = final_spin(m1,m2, p['spin1'],p['spin2'])
            X,Y, PDF =  twod_kde(mass_f,spin_f)
            levs = np.sort(FindHeightForLevel(np.log(PDF),[0.9]))
            ax3.contour(X,Y,np.log(PDF),levs, linestyles = (ls,ls), colors=(c,c,),linewidths=1.5)
            handles.append(mlines.Line2D([], [], linestyle = ls, color=c, label=e.upper()))
                #ax4 - spins SEOB
            p_SEOB = np.genfromtxt(os.path.join(path+'/SEOB',e+'/Nested_sampler/posterior.dat'),names= True)
            X,Y, PDF =  twod_kde(p_SEOB['spin1'],p_SEOB['spin2'])
            levs = np.sort(FindHeightForLevel(np.log(PDF),[0.9]))
            ax4.contour(X,Y,np.log(PDF),levs, linestyles = (ls,ls), colors=(c,c,),linewidths=1.5)
#        except:
#            print(e+' posterior not found')

    #fig1.legend(bbox_to_anchor=(.9, 0.5), loc='center left', borderaxespad=0.,fancybox=True,handles=handles)
    fig1.legend(bbox_to_anchor=(0.45, 0.2), loc='upper center', ncol=4, borderaxespad=0.,fancybox=True,handles=handles)
    fig2.legend(bbox_to_anchor=(0.45, 0.2), loc='upper center', ncol=4, borderaxespad=0.,fancybox=True,handles=handles)
    #fig2.legend(bbox_to_anchor=(.9, 0.5), loc='center left', borderaxespad=0.,fancybox=True, handles=handles)
    #fig3.legend(bbox_to_anchor=(.9, 0.5), loc='center left', borderaxespad=0.,fancybox=True,handles=handles)
    
    ax1.set_xlabel('$m_1/M_\odot$')
    ax1.set_ylabel('$m_2/M_\odot$')
    ax1.set_xlim(0,80)
    ax1.set_ylim(0,50)
    ax2.set_xlabel('$s_1$')
    ax2.set_ylabel('$s_2$')
    ax3.set_xlabel('$M_f/M_\odot$')
    ax3.set_ylabel(r"$\mathit{a_f}$")
    ax3.set_xlim(0,110)
    ax3.set_ylim(0.4,1)
    ax4.set_xlabel('$s_1$')
    ax4.set_ylabel('$s_2$')

    ax2.set_title("mlgw - TEOBResumS")
    ax4.set_title("mlgw - SEOBNRv4")

    #filling ax1
    line_width = 0.3
    ax1.fill_between([0,50],[0,50],[1000,1000], color = 'grey', zorder=100) #excluding q<1 region
    ax1.plot([0,200], [0,100], ls = '--', lw = line_width, color = 'grey', label = "2")
    ax1.plot([0,400], [0,100], ls = '-.', lw = line_width, color = 'grey', label = "4")
    ax1.plot([0,800], [0,100], ls = ':', lw = line_width, color = 'grey', label = "8")
    #ax1.legend(loc="upper left", fontsize = 8, edgecolor = 'black', facecolor = 'silver', title = "q" , framealpha = 1.)
    legend = ax1.legend(frameon = 1,loc="upper left", fontsize = 8, edgecolor = 'black', facecolor = 'silver', title = "q" , framealpha = 1., borderpad = 1.)
    legend.set_zorder(102)

#    fig1.savefig('/Users/wdp/repositories/MLGW/paper/img/posterior_masses_source.pdf', bbox_inches='tight')
#    fig2.savefig('/Users/wdp/repositories/MLGW/paper/img/spins.pdf', bbox_inches='tight')
    #fig1.savefig('../img/posterior_masses_final.pdf', bbox_inches='tight', pad_inches = .2)
    fig1.savefig('../img/posterior_masses_final.pdf', bbox_inches='tight', pad_inches = .2)
    fig2.savefig('../img/spins_TEOB_SEOB.pdf', bbox_inches='tight', pad_inches = .2)
    #fig3.savefig('../img/final_spin_mass.pdf', bbox_inches='tight', pad_inches = .2)




