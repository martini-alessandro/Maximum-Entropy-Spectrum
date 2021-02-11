import matplotlib.pyplot as plt

def init_plotting():
    """
    Initialises matplotlib *rcParams*
    
    :return: no return
    """
    
    # plotting options
    plt.rcParams['figure.figsize']             = (3.4, 3.4)
    plt.rcParams['font.size']                  = 11
    plt.rcParams['font.family']                = 'DejaVu Sans'
    plt.rcParams['font.sans-serif']            = ['Bitstream Vera Sans']
    plt.rcParams['axes.labelsize']             = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize']             = plt.rcParams['font.size']
    plt.rcParams['legend.fontsize']            = 10
    plt.rcParams['xtick.labelsize']            = 9
    plt.rcParams['ytick.labelsize']            = 9
    plt.rcParams['xtick.major.size']           = 3
    plt.rcParams['xtick.minor.size']           = 3
    plt.rcParams['xtick.major.width']          = 1
    plt.rcParams['xtick.minor.width']          = 1
    plt.rcParams['ytick.major.size']           = 3
    plt.rcParams['ytick.minor.size']           = 3
    plt.rcParams['ytick.major.width']          = 1
    plt.rcParams['ytick.minor.width']          = 1
    plt.rcParams['legend.frameon']             = False
    plt.rcParams['legend.loc']                 = 'center left'
    plt.rcParams['axes.linewidth']             = 1
    plt.rcParams['contour.negative_linestyle'] = 'solid'
    fig = plt.figure()
    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')
    return fig
