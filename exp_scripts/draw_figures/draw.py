
import matplotlib.pyplot as plt

import numpy as np
import random

gblue = '#4486F4'
gred = '#DA483B'
gyellow = '#FFC718'
ggreen = '#1CA45C'

def beetle_results_stages():
    
    # 124*208*208

    fig, ax = plt.subplots(figsize=(6,3.5))
    ax.set_xlabel('Different VTK-m backends', fontsize=13.5)
    ax.set_ylabel('Log(Time(ms))', fontsize=13.5)
    #ax.set_ylim([0,1800])
    ax.set_ylim([0,10])

    N = 4
    ind = np.arange(N)  # the x locations for the groups
    width = 0.2     # the width of the bars
    ax.set_xticks(ind + 2.5*width)
    ax.set_xticklabels(('OpenMP','Cuda','Kokkos_Cuda','Kokkos_Hip'), fontsize=12.8)

    sample_bykey_mean = (21.84766667,25.37676333,74.05266667,49.861)
    sample_bykey_stdev = (0.6148986366,0.2279637709,3.764064868,1.078393713)
    p1 = ax.bar(ind+ width, np.log(sample_bykey_mean),  width, color=[gred]*3, capsize=3, alpha=0.9)


    sample_fixed_mean = (4.135666667,15.63343333,13.832,16.732)
    sample_fixed_std = (0.03700450423,3.723371612,0.1723136675,0.1859758049)
    p2 = ax.bar(ind+ 2*width, np.log(sample_fixed_mean),  width, color=[gblue]*3, capsize=3, alpha=0.9)
    

    unmetrics_mean = (1505.236667,1720.973333,234.0573333,265.2363333)
    unmetrics_std=(2.101816675,22.16536111,3.663056966,0.2283688537)
    p3 = ax.bar(ind+ 3*width, np.log(unmetrics_mean),  width, color=[ggreen]*3, capsize=3, alpha=0.9)
   
    unmetrics_ig_mean = (11.116,18.9803,11.659,1.219)
    unmetrics_ig_std=(0.0240208243,0.009600520819,0.01907878403,0.008544003745)
    p4 = ax.bar(ind+ 4*width, np.log(unmetrics_ig_mean),  width, color=[gyellow]*3, capsize=3, alpha=0.9)
    
    ax.legend(('Subsampling (key)', 'Subsampling (fixed)', 'Metrics (MG_1k)','Metrics (IG)'), loc='upper center', ncol=2, fontsize='large')

    plt.savefig("beetle_results_small.png",bbox_inches='tight')
    plt.savefig("beetle_results_small.pdf",bbox_inches='tight')


    fig, ax = plt.subplots(figsize=(6,3.5))
    ax.set_xlabel('Different VTK-m backends', fontsize=13.5)
    ax.set_ylabel('Log(Time(ms))', fontsize=13.5)
    #ax.set_ylim([0,42000])
    ax.set_ylim([0,14])

    N = 4
    ind = np.arange(N)  # the x locations for the groups
    width = 0.2     # the width of the bars
    ax.set_xticks(ind + 2.5*width)
    ax.set_xticklabels(('OpenMP','Cuda','Kokkos_Cuda','Kokkos_HIP'), fontsize=12.8)

    sample_bykey_mean = (791.5623333,841.0503333,3432.145667,1519.451667)
    sample_bykey_stdev = (12.7493447,4.518449328,27.91769959,5.244509732)
    p1 = ax.bar(ind+ width, np.log(sample_bykey_mean),  width, color=[gred]*3, capsize=3, alpha=0.9)


    sample_fixed_mean = (163.732,384.008,468.912,160.3883333)
    sample_fixed_std = (8.390377763,15.31294848,2.82639187,10.92071867)
    p2 = ax.bar(ind+ 2*width, np.log(sample_fixed_mean),  width, color=[gblue]*3, capsize=3, alpha=0.9)
    
    # mg 
    # if use original value, it should be 0
    # log 1 is 0
    unmetrics_mean = (41385.16667,23676.76667,7748.923333,1)
    unmetrics_std=(20.65098868,71.30822767,25.25309156,1)
    p3 = ax.bar(ind+ 3*width, np.log(unmetrics_mean),  width, color=[ggreen]*3, capsize=3, alpha=0.9)

    unmetrics_ig_mean = (10842.86667,175.0313333,66.485,35.35933333)
    unmetrics_ig_std=(7.622554253,0.4154760322,0.3795009881,0.04701418226)
    p4 = ax.bar(ind+ 4*width, np.log(unmetrics_ig_mean),  width, color=[gyellow]*3, capsize=3, alpha=0.9)
    
    # label the empty position
    # plt.text(3.55, 0.1, "X")

    ax.legend(('Subsampling (key)', 'Subsampling (fixed)', 'Metrics (MG_1k)', 'Metrics (IG)'), loc='upper center', ncol=2, fontsize=11.5)

    plt.savefig("beetle_results_large.png",bbox_inches='tight')
    plt.savefig("beetle_results_large.pdf",bbox_inches='tight')

def wind_reuslts_omp_gpu():

    fig, ax = plt.subplots(figsize=(6,3.5))
    ax.set_xlabel('Number of Monte-Carlo samples', fontsize=13.5)
    ax.set_ylabel('Time(ms)', fontsize=13.5)
    ax.set_ylim([0,3000])

    N = 4
    ind = np.arange(N)  # the x locations for the groups
    width = 0.3      # the width of the bars
    ax.set_xticks(ind + 1.5*width)
    ax.set_xticklabels(('1000','2000','4000','8000'), fontsize=13.5)

    openmp_mean = (342.7656667, 727.5533333,1385.643333,2700.696667)
    openmp_std = (0.0112398102,4.148983289,5.792066413,5.944714739)
    p2 = ax.bar(ind+ width, openmp_mean,  width, color=[gred]*3, capsize=3, yerr=openmp_std, ecolor='grey')


    cuda_mean = (45.57553333,106.6176667,149.0303333,254.9003333)
    cuda_std = (0.07820206732,2.948390126,2.655132075,17.84778937)
    p3 = ax.bar(ind+ 2*width, cuda_mean,  width, color=[gblue]*3, capsize=3, yerr=cuda_std, ecolor='grey')
    
    ax.legend(('OpenMP', 'Cuda'), loc='upper left', ncol=1, fontsize='large')


    plt.savefig("wind_results_omp_gpu.png",bbox_inches='tight')
    plt.savefig("wind_results_omp_gpu.pdf",bbox_inches='tight')    

def strong_scale():
    fig, ax = plt.subplots(figsize=(6,3.5))
    # set titles
    ax.set_xlabel('The number of gpu', fontsize='large')
    ax.set_ylabel('Time(ms)', fontsize='large')
    ax.grid(axis='y')
    ax.set_axisbelow(True)

    ax.set_ylim([6,15])

    # set tick
    N = 6

    plt.xticks(range(N), ['4' , '8', '16', '32', '64', '128'], fontsize='large')
    # set the value of the figure here
    # use the capsize to control the error bar
    v1 = (2587.53,  1395.59 , 800.964 , 461.165, 294.574, 196.13)
    p1 = ax.plot(np.log2(v1), color=gblue, marker='^', label='redsea')    

    v2 = (28773.8, 16329, 8474.02, 5018, 2857.01, 1856.51)
    p2 = ax.plot(np.log2(v2), color=gred, marker='.', label='supernova')    

    ax.legend(ncol=2, fontsize='large')
    
    plt.savefig("strong_scale_redsea_supernova.png",bbox_inches='tight')
    plt.savefig("strong_scale_redsea_supernova.pdf",bbox_inches='tight')


if __name__ == "__main__":
    beetle_results_stages()
    #wind_reuslts_omp_gpu()
    #strong_scale()
