
import matplotlib.pyplot as plt

import numpy as np
import random

gblue = '#4486F4'
gred = '#DA483B'
gyellow = '#FFC718'
ggreen = '#1CA45C'

def beetle_results_omp_gpu():

    # uni distribution
    fig, ax = plt.subplots(figsize=(6,3.5))
    ax.set_xlabel('Different stages of computing uncertainty metrics', fontsize=13.5)
    ax.set_ylabel('Time(ms)', fontsize=13.5)
    ax.set_ylim([0,750])

    N = 3
    ind = np.arange(N)  # the x locations for the groups
    width = 0.3      # the width of the bars
    ax.set_xticks(ind + 1.5*width)
    ax.set_xticklabels(('Labeling','Downsampling','Uncertainty metrics'), fontsize=13.5)

    openmp_mean = (251.3486667,585.024,266.6013333)
    openmp_std = (0.6012839041,13.13104455,0.2597389715)
    p2 = ax.bar(ind+ width, openmp_mean,  width, color=[gred]*3, capsize=3, yerr=openmp_std, ecolor='grey')


    cuda_mean = (97.5622,538.1776667,185.4483333)
    cuda_std = (4.972551425,10.74274776,10.25276774)
    p3 = ax.bar(ind+ 2*width, cuda_mean,  width, color=[gblue]*3, capsize=3, yerr=cuda_std, ecolor='grey')
    
    ax.legend(('OpenMP', 'Cuda'), loc='upper left', ncol=1, fontsize='large')


    plt.savefig("beetle_results_uni_omp_gpu.png",bbox_inches='tight')
    plt.savefig("beetle_results_uni_omp_gpu.pdf",bbox_inches='tight')


    # ig distribution
    fig, ax = plt.subplots(figsize=(6,3.5))
    ax.set_xlabel('Different stages of computing uncertainty metrics', fontsize=13.5)
    ax.set_ylabel('Time(ms)', fontsize=13.5)
    ax.set_ylim([0,750])
    
    N = 3
    ind = np.arange(N)  # the x locations for the groups
    width = 0.3      # the width of the bars
    ax.set_xticks(ind + 1.5*width)
    ax.set_xticklabels(('Labeling','Downsampling','Uncertainty metrics'), fontsize=13.5)


    openmp_mean = (252.536,610.8803333,263.87)
    openmp_std = (0.06605300902,6.342306941,0.1704024648)
    p2 = ax.bar(ind+ width, openmp_mean,  width, color=[gred]*3, capsize=3, yerr=openmp_std, ecolor='grey')


    cuda_mean = (93.88216667,530.466,185.5746667)
    cuda_std = (7.110525455,7.525666549,10.69285389)
    p3 = ax.bar(ind+ 2*width, cuda_mean,  width, color=[gblue]*3, capsize=3, yerr=cuda_std, ecolor='grey')

    ax.legend(('OpenMP', 'Cuda'), loc='upper left', ncol=1, fontsize='large')

    plt.savefig("beetle_results_ig_omp_gpu.png",bbox_inches='tight')
    plt.savefig("beetle_results_ig_omp_gpu.pdf",bbox_inches='tight')

    # mg distribution
    # beetle_results_mg_omp_gpu split

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,figsize=(6,3.5))
    # set titles
    fig.text(0.006, 0.5, 'Time(ms)', va='center', rotation='vertical', fontsize=13.5)
    ax1.set_ylim([10000, 100000])
    ax2.set_ylim([0, 800])
    ax2.set_xlabel('Different stages of computing uncertainty metrics', fontsize=13.5)
    # delet the line at the bottom and the top of the two figures
    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax1.tick_params(axis='x', bottom=False, labeltop=False, labelbottom=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()

    d = 0.01  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    plt.subplots_adjust(hspace=.05)  # adjust distance between two subplots

    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)

    # set tick
    N = 3
    ind = np.arange(N)  # the x locations for the groups
    width = 0.3      # the width of the bars
    ax2.set_xticks(ind + 0.5*width)
    ax2.set_xticklabels(('Labeling','Downsampling','Uncertainty metrics'),fontsize=13.5)

    openmp_mean = (252.436,621.306,90478.16667)
    openmp_std = (0.2285278976,4.56341988,13.71726406)
    p1 = ax1.bar(ind, openmp_mean,  width, color=[gred]*3, yerr=openmp_std, capsize=3, ecolor='grey')

    cuda_mean = (104.003,722.8853333,17542.36667)
    cuda_std = (4.715512804,9.58137424,106.5679283)
    p2 = ax1.bar(ind+ width, cuda_mean,  width, color=[gblue]*3, yerr=cuda_std, capsize=3, ecolor='grey')

    p1 = ax2.bar(ind, openmp_mean,  width, color=[gred]*6, yerr=openmp_std, capsize=3, ecolor='grey')
    p2 = ax2.bar(ind + width, cuda_mean,  width, color=[gblue]*6, yerr=cuda_std, capsize=3, ecolor='grey')

    # add the lengend for the data by defualt
    # this can be at the centric legend
    # ax1.legend((p1[0], p2[0], p3[0], p4[0]), ('producer-responsible', 'consumer-responsible','metadata-responsible','notification-responsible'),bbox_to_anchor=(0.83, 1.0), ncol=2, fontsize='medium')

    ax1.legend((p1[0],p2[0]), ('OpenMP', 'Cuda'), loc='upper left', ncol=1, fontsize='large')

    plt.savefig("beetle_results_mg_omp_gpu_split.png",bbox_inches='tight')
    plt.savefig("beetle_results_mg_omp_gpu_split.pdf",bbox_inches='tight')

                                                                                                                                                                           
if __name__ == "__main__":
    beetle_results_omp_gpu()
