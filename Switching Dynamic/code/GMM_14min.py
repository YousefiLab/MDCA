#########  clone and install necessary packages from github repositories

import statistics
from scipy import stats as st
import scipy.io as sio
from sklearn import mixture
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

# nested function
def calc_mode_mean(lbl, t_min, t_flag, ind_flag, k_):
    mode_lbl = np.zeros((len(t_flag)-1))
    mean_lbl = np.zeros((len(t_flag)-1, k_))
    t_min_mode = np.zeros((len(t_flag)-1))
    for i_ in range(len(t_flag) - 1):
        ind_l = ind_flag[i_]
        ind_u = ind_flag[i_ + 1]

        t_i = t_min[ind_l: ind_u]
        t_min_mode[i_] = np.mean(t_i)

        lbl_i = lbl[ind_l: ind_u]
        mode_i= st.mode(lbl_i)
        if np.isscalar(mode_i):
            mode_lbl[i_] = mode_i
        else:
            mode_lbl[i_] = mode_i[0]

        len_lbl_i = lbl_i.shape[0]

        mean_lbl_i = np.zeros((k_, ))
        for i__ in range(k_):
            mean_lbl_i[i__] = np.asarray(np.where(lbl_i == i__)).shape[1] / len_lbl_i
        mean_lbl[i_, :] = mean_lbl_i

    return mode_lbl, mean_lbl, t_min_mode

def plt_mode(t_mode, mode_lbl_i, Color_, Mark_, k_, model_name):
    fig, ax = plt.subplots()
    # mode_lbl_i = np.asarray([2, 4, 1, 4, 1, 0, 2])
    for kk in range(k_):
        ind_kk = np.argwhere(mode_lbl_i == kk)

        clr = 'black'
        plt.plot(t_mode[ind_kk], mode_lbl_i[ind_kk]+1,
                 color=clr,
                 marker=Mark_[0],
                 markerfacecolor=clr,
                 markeredgecolor=clr,
                 linestyle='None',
                 markersize=20)

    alpha_ = 0.8
    # EO
    ax.axvspan(0, 2, alpha=alpha_, color='#1878b8')
    ax.axvspan(12, 14, alpha=alpha_, color='#1878b8')
    # R
    ax.axvspan(2, 4, alpha=alpha_, color='#f2f233')
    ax.axvspan(6, 8, alpha=alpha_, color='#f2f233')
    # MS
    ax.axvspan(4, 6, alpha=alpha_, color='#37a151')
    ax.axvspan(8, 10, alpha=alpha_, color='#37a151')
    # EC
    ax.axvspan(10, 12, alpha=alpha_, color='#828c85')

    str_tit = '%s - mode - K = %d' % (model_name, k_)
    plt.title(str_tit)
    plt.xlabel('Time (min)')
    plt.ylim([.5, k_+1.5])
    plt.xlim([0, 14])
    plt.ylabel('Cluster')
    plt.rcParams.update({'font.size': 20})
    # plt.legend()
    plt.show()

    return 0

def plt_mean(t_mode, mean_lbl, Color_, Mark_, k_, model_name):

    fig, ax = plt.subplots()
    for i_ in range(t_mode.shape[0]):
        mean_lbl_i = mean_lbl[i_, :]
        t_i = t_mode[i_]

        for kk in range(k_):
            plt.plot(t_i, kk,
             color=Color_[kk],
             marker=Mark_[0],
             markerfacecolor=Color_[kk],
             markeredgecolor=Color_[kk],
             linestyle='None',
             markersize=30*mean_lbl_i[kk])

    alpha_ = 0.1
    # EO
    ax.axvspan(0, 2, alpha=alpha_, color='red')
    ax.axvspan(12, 14, alpha=alpha_, color='red')
    # R
    ax.axvspan(2, 4, alpha=alpha_ + .3, color='blue')
    ax.axvspan(6, 8, alpha=alpha_ + .3, color='blue')
    # MS
    ax.axvspan(4, 6, alpha=alpha_ + .5, color='yellow')
    ax.axvspan(8, 10, alpha=alpha_ + .5, color='yellow')
    # EC
    ax.axvspan(10, 12, alpha=alpha_ + .4, color='green')

    str_tit = '%s - wheighted mean - K = %d' % (model_name, k_)
    plt.title(str_tit)
    plt.xlabel('time (min)')
    # plt.legend()
    plt.show()

    u=1

    u=1


    return 0

def plt_all(lbl, t_min, Color_, Mark_, k_, model_name):
    fig, ax = plt.subplots()

    for kk in range(k_):
        ind_kk = np.argwhere(lbl == kk)

        plt.plot(t_min[ind_kk], lbl[ind_kk], color=Color_[kk],
                 marker=Mark_[1],
                 markerfacecolor=Color_[kk],
                 markeredgecolor=Color_[kk],
                 linestyle='None',
                 label='%d' % (kk))

    alpha_ = 0.1
    # EO
    ax.axvspan(0, 2, alpha=alpha_, color='red')
    ax.axvspan(12, 14, alpha=alpha_, color='red')
    # R
    ax.axvspan(2, 4, alpha=alpha_ + .3, color='blue')
    ax.axvspan(6, 8, alpha=alpha_ + .3, color='blue')
    # MS
    ax.axvspan(4, 6, alpha=alpha_ + .5, color='yellow')
    ax.axvspan(8, 10, alpha=alpha_ + .5, color='yellow')
    # EC
    ax.axvspan(10, 12, alpha=alpha_ + .4, color='green')

    str_tit = '%s - K = %d' % (model_name, k_)
    plt.title(str_tit)
    plt.xlabel('time (min)')
    # plt.legend()
    plt.show()

    return 0

def plt_all_mode(lbl, t_min, mode_lbl_i, t_mode, Color_, Mark_, k_, model_name):
    fig, ax = plt.subplots()

    for kk in range(k_):
        ind_kk = np.argwhere(lbl == kk)

        clr = 'black'
        plt.plot(t_min[ind_kk], lbl[ind_kk] + 1,
                 color=clr,
                 marker=Mark_[1],
                 markerfacecolor=clr,
                 markeredgecolor=clr,
                 linestyle='None',
                 label='%d' % (kk))

        # mode
        ind_kk = np.argwhere(mode_lbl_i == kk)

        clr = 'black'
        plt.plot(t_mode[ind_kk], mode_lbl_i[ind_kk] + 1,
                 color=clr,
                 marker=Mark_[0],
                 markerfacecolor=clr,
                 markeredgecolor=clr,
                 linestyle='None',
                 markersize=20)


    alpha_ = 0.8
    # EO
    ax.axvspan(0, 2, alpha=alpha_, color='#1878b8')
    ax.axvspan(12, 14, alpha=alpha_, color='#1878b8')
    # R
    ax.axvspan(2, 4, alpha=alpha_, color='#f2f233')
    ax.axvspan(6, 8, alpha=alpha_, color='#f2f233')
    # MS
    ax.axvspan(4, 6, alpha=alpha_, color='#37a151')
    ax.axvspan(8, 10, alpha=alpha_, color='#37a151')
    # EC
    ax.axvspan(10, 12, alpha=alpha_, color='#828c85')

    str_tit = '%s - K = %d' % (model_name, k_)
    plt.title(str_tit)
    plt.xlabel('time (min)')
    # plt.legend()
    plt.show()

    return 0

########## code

str_data = 'ken'


mat_fname = '%s_gmm.mat'%(str_data)
mat_contents = sio.loadmat(mat_fname)

dict_da = mat_contents['dict_dat']
eVec_fr_l_u = dict_da[0, 0]['eVec_fr_l_u']
t_min = dict_da[0, 0]['t_min']
t_min = t_min.reshape((t_min.shape[1], ))
X = np.log(eVec_fr_l_u)

t_flag = [0, 2, 4, 6, 8, 10, 12, t_min[-1]]
ind_flag = np.zeros((len(t_flag), ))
cnt_ = -1
for t_ in t_flag:
    cnt_ +=1
    print(t_)
    diff_ = np.abs(t_min - t_)

    ind_flag[cnt_] = np.asarray(np.argwhere(diff_== np.min(diff_)))[0, 0]

ind_flag = ind_flag.astype(int)

K = np.arange(4, 8)

BIC_gmm = np.zeros((len(K), ))

cnt_ = -1
Color_ = ['#e00724', '#613af0', '#edab4e', '#85ede3', '#6b8067', '#4f2d66', '#0856bd', '#850d35',  '#4f2d66', '#0856bd', '#850d35',
          '#e00724', '#613af0', '#edab4e', '#85ede3', '#6b8067', '#4f2d66', '#0856bd', '#850d35',  '#4f2d66', '#0856bd', '#850d35']
Mark_ = ['o', '.']
gmm_ = 1
kmeans_ = 0
for k_ in K:

    cnt_ +=1

    if gmm_:
        gmm = mixture.GaussianMixture(n_components=k_, covariance_type="spherical")

        gmm.fit(X)

        bic_gmm = gmm.bic(X)
        BIC_gmm[cnt_] = bic_gmm
        lbl_gmm = gmm.predict(X)

        mode_lbl_gmm_i, mean_lbl_gmm_i, t_mode =  calc_mode_mean(lbl_gmm, t_min, t_flag, ind_flag, k_)

        model_name = 'gmm'
        # plot mode
        plt_mode(t_mode, mode_lbl_gmm_i, Color_, Mark_, k_, 'gmm')

        # plot mean
        # plt_mean(t_mode, mean_lbl_gmm_i, Color_, Mark_, k_)

        # plot all points
        # plt_all(lbl_gmm, t_min, Color_, Mark_, k_ ,'gmm')

        plt_all_mode(lbl_gmm, t_min, mode_lbl_gmm_i, t_mode, Color_, Mark_, k_, model_name)

        u=1
    elif kmeans_:
        kmeans = KMeans(n_clusters=k_, random_state=0).fit(X)

        lbl_kmeans = kmeans.labels_

        model_name = 'kmeans'
        



    u=1

if gmm_:
    plt.figure()
    plt.plot(K, BIC_gmm, marker=Mark_[0], linestyle = 'None')
    plt.title("BIC")
    plt.xlabel('number of clusters')
    plt.show()

    ind_k = np.asarray(np.where(BIC_gmm == np.min(BIC_gmm)))[0, 0]

    K_best = K[ind_k]

    print('Best number of cluster is = %d'%K_best)

u=1


u=1