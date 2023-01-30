
import matplotlib.pyplot as plt

import numpy as np
import random

gblue = '#4486F4'
gred = '#DA483B'
gyellow = '#FFC718'
ggreen = '#1CA45C'

def beetle_results():

    # uni distribution
    fig, ax = plt.subplots(figsize=(7,4.6))
    ax.set_xlabel('Stages of computing uncertainty metrics', fontsize='large')
    ax.set_ylabel('Time(ms)', fontsize='large')

    N = 3
    ind = np.arange(N)*2.5  # the x locations for the groups
    width = 0.25       # the width of the bars
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('Labeling','Down-sampled data','Uncertainty metrics'), fontsize='large')

    serial_mean = (10554.63333,14479.36667,11067.76667)
    serial_std = (104.7227928,199.988508,31.79690761)
    p1 = ax.bar(ind, serial_mean,  width, color=[gblue]*3, capsize=3, yerr=serial_std, ecolor='grey')


    openmp_mean = (250.8206667,577.1533333,264.392)
    openmp_std = (0.347324536,13.03131188,0.3538488378)
    p2 = ax.bar(ind+ width, openmp_mean,  width, color=[gyellow]*3, capsize=3, yerr=openmp_std, ecolor='grey')


    cuda_mean = (166.4493333,630.3436667,21.16773333)
    cuda_std = (5.972815863,5.292360658,0.1159116186)
    p3 = ax.bar(ind+ 2*width, cuda_mean,  width, color=[gred]*3, capsize=3, yerr=cuda_std, ecolor='grey')

    plt.savefig("beetle_results_uni.png",bbox_inches='tight')


    # ig distribution
    fig, ax = plt.subplots(figsize=(7,4.6))
    ax.set_xlabel('Stages of computing uncertainty metrics', fontsize='large')
    ax.set_ylabel('Time(ms)', fontsize='large')

    N = 3
    ind = np.arange(N)*2.5  # the x locations for the groups
    width = 0.25       # the width of the bars
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('Labeling','Down-sampled data','Uncertainty metrics'), fontsize='large')

    serial_mean = (10468.33333,15153.13333,10741.53333)
    serial_std = (33.41800912,92.90932856,30.8230974)
    p1 = ax.bar(ind, serial_mean,  width, color=[gblue]*3, capsize=3, yerr=serial_std, ecolor='grey')


    openmp_mean = (251.4673333,587.5676667,258.049)
    openmp_std = (0.06798774399,7.164214286,0.09814784766)
    p2 = ax.bar(ind+ width, openmp_mean,  width, color=[gyellow]*3, capsize=3, yerr=openmp_std, ecolor='grey')


    cuda_mean = (169.663,666.5686667,23.2751)
    cuda_std = (2.194154963,33.53476901,0.2380756812)
    p3 = ax.bar(ind+ 2*width, cuda_mean,  width, color=[gred]*3, capsize=3, yerr=cuda_std, ecolor='grey')

    plt.savefig("beetle_results_ig.png",bbox_inches='tight')

    # mg distribution
    fig, ax = plt.subplots(figsize=(7,4.6))
    ax.set_xlabel('Stages of computing uncertainty metrics', fontsize='large')
    ax.set_ylabel('Time(ms)', fontsize='large')

    N = 3
    ind = np.arange(N)*2.5  # the x locations for the groups
    width = 0.25       # the width of the bars
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('Labeling','Down-sampled data','Uncertainty metrics'), fontsize='large')

    serial_mean = (10500.3,14943.6,3732350)
    serial_std = (0,0,0)
    p1 = ax.bar(ind, serial_mean,  width, color=[gblue]*3, capsize=3, yerr=serial_std, ecolor='grey')


    openmp_mean = (252.802,625.815,90322.3)
    openmp_std = (0,0,0)
    p2 = ax.bar(ind+ width, openmp_mean,  width, color=[gyellow]*3, capsize=3, yerr=openmp_std, ecolor='grey')


    cuda_mean = (165.696,637.963,234.827)
    cuda_std = (0,0,0)
    p3 = ax.bar(ind+ 2*width, cuda_mean,  width, color=[gred]*3, capsize=3, yerr=cuda_std, ecolor='grey')

    plt.savefig("beetle_results_mg.png",bbox_inches='tight')


def beetle_results():

    # uni distribution
    fig, ax = plt.subplots(figsize=(7,4.6))
    ax.set_xlabel('Stages of computing uncertainty metrics', fontsize='large')
    ax.set_ylabel('Time(ms)', fontsize='large')

    N = 3
    ind = np.arange(N)*2.5  # the x locations for the groups
    width = 0.25       # the width of the bars
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('Labeling','Down-sampled data','Uncertainty metrics'), fontsize='large')

    serial_mean = (10554.63333,14479.36667,11067.76667)
    serial_std = (104.7227928,199.988508,31.79690761)
    p1 = ax.bar(ind, serial_mean,  width, color=[gblue]*3, capsize=3, yerr=serial_std, ecolor='grey')


    openmp_mean = (250.8206667,577.1533333,264.392)
    openmp_std = (0.347324536,13.03131188,0.3538488378)
    p2 = ax.bar(ind+ width, openmp_mean,  width, color=[gyellow]*3, capsize=3, yerr=openmp_std, ecolor='grey')


    cuda_mean = (166.4493333,630.3436667,21.16773333)
    cuda_std = (5.972815863,5.292360658,0.1159116186)
    p3 = ax.bar(ind+ 2*width, cuda_mean,  width, color=[gred]*3, capsize=3, yerr=cuda_std, ecolor='grey')

    plt.savefig("beetle_results_uni.png",bbox_inches='tight')


    # ig distribution
    fig, ax = plt.subplots(figsize=(7,4.6))
    ax.set_xlabel('Stages of computing uncertainty metrics', fontsize='large')
    ax.set_ylabel('Time(ms)', fontsize='large')

    N = 3
    ind = np.arange(N)*2.5  # the x locations for the groups
    width = 0.25       # the width of the bars
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('Labeling','Down-sampled data','Uncertainty metrics'), fontsize='large')

    serial_mean = (10468.33333,15153.13333,10741.53333)
    serial_std = (33.41800912,92.90932856,30.8230974)
    p1 = ax.bar(ind, serial_mean,  width, color=[gblue]*3, capsize=3, yerr=serial_std, ecolor='grey')


    openmp_mean = (251.4673333,587.5676667,258.049)
    openmp_std = (0.06798774399,7.164214286,0.09814784766)
    p2 = ax.bar(ind+ width, openmp_mean,  width, color=[gyellow]*3, capsize=3, yerr=openmp_std, ecolor='grey')


    cuda_mean = (169.663,666.5686667,23.2751)
    cuda_std = (2.194154963,33.53476901,0.2380756812)
    p3 = ax.bar(ind+ 2*width, cuda_mean,  width, color=[gred]*3, capsize=3, yerr=cuda_std, ecolor='grey')

    plt.savefig("beetle_results_ig.png",bbox_inches='tight')

    # mg distribution
    fig, ax = plt.subplots(figsize=(7,4.6))
    ax.set_xlabel('Stages of computing uncertainty metrics', fontsize='large')
    ax.set_ylabel('Time(ms)', fontsize='large')

    N = 3
    ind = np.arange(N)*2.5  # the x locations for the groups
    width = 0.25       # the width of the bars
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('Labeling','Down-sampled data','Uncertainty metrics'), fontsize='large')

    serial_mean = (10500.3,14943.6,3732350)
    serial_std = (0,0,0)
    p1 = ax.bar(ind, serial_mean,  width, color=[gblue]*3, capsize=3, yerr=serial_std, ecolor='grey')


    openmp_mean = (252.802,625.815,90322.3)
    openmp_std = (0,0,0)
    p2 = ax.bar(ind+ width, openmp_mean,  width, color=[gyellow]*3, capsize=3, yerr=openmp_std, ecolor='grey')


    cuda_mean = (165.696,637.963,234.827)
    cuda_std = (0,0,0)
    p3 = ax.bar(ind+ 2*width, cuda_mean,  width, color=[gred]*3, capsize=3, yerr=cuda_std, ecolor='grey')

    plt.savefig("beetle_results_mg.png",bbox_inches='tight')

    # TODO, add figure to show total time speed up for openmp and the cuda cases
    # for different distributions

def beetle_results_omp_gpu():

    # uni distribution
    fig, ax = plt.subplots(figsize=(7,4.6))
    ax.set_xlabel('Stages of computing uncertainty metrics', fontsize='large')
    ax.set_ylabel('Time(ms)', fontsize='large')

    N = 3
    ind = np.arange(N)*2.5  # the x locations for the groups
    width = 0.25       # the width of the bars
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('Labeling','Down-sampled data','Uncertainty metrics'), fontsize='large')


    openmp_mean = (250.8206667,577.1533333,264.392)
    openmp_std = (0.347324536,13.03131188,0.3538488378)
    p2 = ax.bar(ind+ width, openmp_mean,  width, color=[gyellow]*3, capsize=3, yerr=openmp_std, ecolor='grey')


    cuda_mean = (166.4493333,630.3436667,21.16773333)
    cuda_std = (5.972815863,5.292360658,0.1159116186)
    p3 = ax.bar(ind+ 2*width, cuda_mean,  width, color=[gred]*3, capsize=3, yerr=cuda_std, ecolor='grey')

    plt.savefig("beetle_results_uni_omp_gpu.png",bbox_inches='tight')


    # ig distribution
    fig, ax = plt.subplots(figsize=(7,4.6))
    ax.set_xlabel('Stages of computing uncertainty metrics', fontsize='large')
    ax.set_ylabel('Time(ms)', fontsize='large')

    N = 3
    ind = np.arange(N)*2.5  # the x locations for the groups
    width = 0.25       # the width of the bars
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('Labeling','Down-sampled data','Uncertainty metrics'), fontsize='large')



    openmp_mean = (251.4673333,587.5676667,258.049)
    openmp_std = (0.06798774399,7.164214286,0.09814784766)
    p2 = ax.bar(ind+ width, openmp_mean,  width, color=[gyellow]*3, capsize=3, yerr=openmp_std, ecolor='grey')


    cuda_mean = (169.663,666.5686667,23.2751)
    cuda_std = (2.194154963,33.53476901,0.2380756812)
    p3 = ax.bar(ind+ 2*width, cuda_mean,  width, color=[gred]*3, capsize=3, yerr=cuda_std, ecolor='grey')

    plt.savefig("beetle_results_ig_omp_gpu.png",bbox_inches='tight')

    # mg distribution
    fig, ax = plt.subplots(figsize=(7,4.6))
    ax.set_xlabel('Stages of computing uncertainty metrics', fontsize='large')
    ax.set_ylabel('Time(ms)', fontsize='large')

    N = 3
    ind = np.arange(N)*2.5  # the x locations for the groups
    width = 0.25       # the width of the bars
    ax.set_xticks(ind + width)
    ax.set_xticklabels(('Labeling','Down-sampled data','Uncertainty metrics'), fontsize='large')


    openmp_mean = (252.802,625.815,90322.3)
    openmp_std = (0,0,0)
    p2 = ax.bar(ind+ width, openmp_mean,  width, color=[gyellow]*3, capsize=3, yerr=openmp_std, ecolor='grey')


    cuda_mean = (165.696,637.963,234.827)
    cuda_std = (0,0,0)
    p3 = ax.bar(ind+ 2*width, cuda_mean,  width, color=[gred]*3, capsize=3, yerr=cuda_std, ecolor='grey')

    plt.savefig("beetle_results_mg_omp_gpu.png",bbox_inches='tight')

def split_bar():
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,figsize=(7,4.6))
    # set titles
    fig.text(0.05, 0.5, 'Time(s)', va='center', rotation='vertical', fontsize=12)
    ax1.set_ylim([60, 100])
    ax2.set_ylim([0, 30])
    ax2.set_xlabel('The percentage of the qualified data', fontsize='large')
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

    plt.subplots_adjust(hspace=.1)  # adjust distance between two subplots

    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)

    # set tick
    N = 6
    ind = np.arange(N)    # the x locations for the groups
    width = 0.18       # the width of the bars
    ax1.set_xticks(ind + 1.5*width)
    ax1.set_xticklabels(('0%','20%','40%','60%','80%','100%'),fontsize='large')

    # set the value of the figure here
    # use the capsize to control the error bar
    tpMeans = (65.316,   67.43356667, 70.38986667, 72.1859, 74.23376667, 76.35183333 )
    tpStd = (0.6099340374,   0.04839900137,   0.5647084056 ,   0.5284408482,    0.8171366246 ,   2.759014564 )
    p1 = ax1.bar(ind, tpMeans,  width, color=[gblue]*5,yerr=tpStd, ecolor='black', capsize=3)

    tcMeans = (70.86473333,  71.28623333, 71.21433333, 71.20403333, 71.39376667, 71.64516667 )
    tcStd = (0.2115885236,   0.03775398434 ,  0.39333012,  0.8282423941,    0.560990502, 0.8876905392 )
    p2 = ax1.bar(ind + width, tcMeans,  width, color=[gred]*5, yerr=tcStd, ecolor='black', capsize=3)


    mpMeans = (78.61916667,  80.78726667, 81.15106667, 78.2954, 80.83366667, 81.21966667)
    mpStd = (1.255949467,    1.611722049, 3.626644518, 1.956503286, 2.674718677, 3.701866909)
    p3 = ax1.bar(ind + 2*width, mpMeans,  width, color=[gyellow]*5, yerr=mpStd, ecolor='black', capsize=3)

    tmMeans = (83.0089,  80.84883333, 81.65736667, 81.65436667, 81.4867, 82.08586667)
    tmStd = (1.304215316 ,   0.3882281589 ,   0.6588356725 ,   0.7596826991 ,   0.8035335774 ,   2.129914957)
    p4 = ax1.bar(ind + 3*width, tmMeans,  width, color=[ggreen]*5, yerr=tmStd, ecolor='black', capsize=3)
    

    p1 = ax2.bar(ind, tpMeans,  width, color=[gblue]*6,yerr=tpStd, ecolor='black', capsize=3)
    p2 = ax2.bar(ind + width, tcMeans,  width, color=[gred]*6, yerr=tcStd, ecolor='black', capsize=3)
    p3 = ax2.bar(ind + 2*width, mpMeans,  width, color=[gyellow]*6, yerr=mpStd, ecolor='black', capsize=3)
    p4 = ax2.bar(ind + 3*width, tmMeans,  width, color=[ggreen]*6, yerr=tmStd, ecolor='black', capsize=3)


    # add the lengend for the data by defualt
    # this can be at the centric legend
    ax1.legend((p1[0], p2[0], p3[0], p4[0]), ('producer-responsible', 'consumer-responsible','metadata-responsible','notification-responsible'),bbox_to_anchor=(0.83, 1.0), ncol=2, fontsize='medium')

    plt.savefig("exp_512_variantpercentage_sim_greater_ana.png",bbox_inches='tight')

if __name__ == "__main__":
    beetle_results()
    beetle_results_omp_gpu()
    split_bar()
