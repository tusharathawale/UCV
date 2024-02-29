import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import numpy as np

gblue = '#4486F4'
gred = '#DA483B'
gyellow = '#FFC718'
ggreen = '#1CA45C'

if __name__ == "__main__":
    # typical bar example
    fig, ax = plt.subplots(figsize=(7,4))

    # set titles
    
    ax.set_xlabel('Number of Monte Carlo samples', fontsize=20)
    ax.set_ylabel('Maximum difference', fontsize=20)

    # set tick
    N = 5
    ind = np.arange(N)    # the x locations for the groups
    width = 0.3       # the width of the bars
    ax.set_xticks(ind+0.5*width)
    ax.set_xticklabels(('1000','2000','3000','4000','5000'), fontsize='large')
    ax.tick_params(axis='x', labelsize=18)
    ax.tick_params(axis='y', labelsize=18)
    # set the value of the figure here
    # use the capsize to control the error bar
    supernova = (0.0360552,0.0340718,0.0274044,0.020167,0.0174966)
    p1 = ax.bar(ind, supernova,  width, color=[gyellow]*5)

    redsea = (0.0350144,0.0220214,0.0141039,0.0105929,0.00958919)
    p2 = ax.bar(ind+width, redsea,  width, color=[gblue]*5)

    ax.legend((p1[0], p2[0]), ('Supernova', 'Red Sea'), fontsize='large')

    plt.savefig("maximum_differences.png",bbox_inches='tight')
    #plt.savefig("exp_componentTime_sim_less_ana.pdf",bbox_inches='tight')