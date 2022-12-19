#########  clone and install necessary packages from github repositories

import statistics
from scipy import stats as st
import scipy.io as sio
from sklearn import mixture
import numpy as np
import matplotlib.pyplot as plt

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

def plt_mode(t_mode, mode_lbl_i, Color_, Mark_, k_):
    fig, ax = plt.subplots()

    for kk in range(k_):
        ind_kk = np.argwhere(mode_lbl_i == kk)

        plt.plot(t_mode[ind_kk], mode_lbl_i[ind_kk],
                 color=Color_[kk],
                 marker=Mark_[0],
                 markerfacecolor=Color_[kk],
                 markeredgecolor=Color_[kk],
                 linestyle='None',
                 markersize=20)

    alpha_ = 0.1
    ax.axvspan(0, 1, alpha=alpha_, color='red')
    ax.axvspan(7, 8, alpha=alpha_, color='red')
    ax.axvspan(1, 2, alpha=alpha_ + .3, color='blue')
    ax.axvspan(5, 6, alpha=alpha_ + .3, color='blue')
    ax.axvspan(2, 3, alpha=alpha_ + .5, color='yellow')
    ax.axvspan(4, 5, alpha=alpha_ + .5, color='yellow')
    ax.axvspan(3, 4, alpha=alpha_ + .4, color='green')
    ax.axvspan(6, 7, alpha=alpha_ + .4, color='green')

    str_tit = 'mode - K = %d' % (k_)
    plt.title(str_tit)
    plt.xlabel('time (min)')
    plt.ylim([-.5, k_+.5])
    # plt.legend()
    plt.show()

    return 0

def plt_mean(t_mode, mean_lbl, Color_, Mark_, k_):

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
    ax.axvspan(0, 1, alpha=alpha_, color='red')
    ax.axvspan(7, 8, alpha=alpha_, color='red')
    ax.axvspan(1, 2, alpha=alpha_ + .3, color='blue')
    ax.axvspan(5, 6, alpha=alpha_ + .3, color='blue')
    ax.axvspan(2, 3, alpha=alpha_ + .5, color='yellow')
    ax.axvspan(4, 5, alpha=alpha_ + .5, color='yellow')
    ax.axvspan(3, 4, alpha=alpha_ + .4, color='green')
    ax.axvspan(6, 7, alpha=alpha_ + .4, color='green')

    str_tit = 'wheighted mean - K = %d' % (k_)
    plt.title(str_tit)
    plt.xlabel('time (min)')
    # plt.legend()
    plt.show()

    u=1

    u=1


    return 0

def plt_all(lbl, t_min, Color_, Mark_, k_):
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
    ax.axvspan(0, 1, alpha=alpha_, color='red')
    ax.axvspan(7, 8, alpha=alpha_, color='red')
    ax.axvspan(1, 2, alpha=alpha_ + .3, color='blue')
    ax.axvspan(5, 6, alpha=alpha_ + .3, color='blue')
    ax.axvspan(2, 3, alpha=alpha_ + .5, color='yellow')
    ax.axvspan(4, 5, alpha=alpha_ + .5, color='yellow')
    ax.axvspan(3, 4, alpha=alpha_ + .4, color='green')
    ax.axvspan(6, 7, alpha=alpha_ + .4, color='green')

    str_tit = 'K = %d' % (k_)
    plt.title(str_tit)
    plt.xlabel('time (min)')
    # plt.legend()
    plt.show()

    return 0

########## code

# str_data = '8_12_2min_11973'
# str_data = 'ali_reading_6_16_8529'
# str_data = 'ali_reading_6_16_eye_8526'
# str_data = 'reza_3'
# str_data = 'ali_0001'
str_data = 'ali_11_30_22'

# str_data = '8_12_2min_4113'
# str_data = '8_12_2min_4699'
# str_data = '8_12_2min_4707'
# str_data = 'A1'
# str_data = 'A2'


mat_fname = '%s_gmm.mat'%(str_data)
mat_contents = sio.loadmat(mat_fname)

dict_da = mat_contents['dict_dat']
eVec_fr_l_u = dict_da[0, 0]['eVec_fr_l_u']
t_min = dict_da[0, 0]['t_min']
t_min = t_min.reshape((t_min.shape[1], ))
X = np.log(eVec_fr_l_u)

t_flag = [0, 1, 2, 3, 4, 5, 6, 7, 8]
ind_flag = np.zeros((len(t_flag), ))
for t_ in t_flag:
    diff_ = np.abs(t_min - t_)
    ind_flag[t_] = np.asarray(np.argwhere(diff_== np.min(diff_)))[0, 0]

ind_flag = ind_flag.astype(int)

K = [2, 3, 4, 5, 6, 7, 8]

BIC = np.zeros((len(K), ))
LBL = np.zeros((len(K), X.shape[0]))

cnt_ = -1
Color_ = ['#e00724', '#613af0', '#edab4e', '#85ede3', '#6b8067', '#4f2d66', '#0856bd', '#850d35']
Mark_ = ['o', '.']
for k_ in K:

    cnt_ +=1
    gmm = mixture.GaussianMixture(n_components=k_, covariance_type="full")

    gmm.fit(X)

    bic = gmm.bic(X)
    BIC[cnt_] = bic
    lbl = gmm.predict(X)
    LBL[cnt_, :] = lbl

    mode_lbl_i, mean_lbl_i, t_mode =  calc_mode_mean(lbl, t_min, t_flag, ind_flag, k_)

    # plot mode
    plt_mode(t_mode, mode_lbl_i, Color_, Mark_, k_)

    # plot mean
    plt_mean(t_mode, mean_lbl_i, Color_, Mark_, k_)

    # plot all points
    plt_all(lbl, t_min, Color_, Mark_, k_)

    u=1

plt.figure()
plt.plot(K, BIC, marker=Mark_[0], linestyle = 'None')
plt.title("BIC")
plt.xlabel('number of clusters')
plt.show()

ind_k = np.asarray(np.where(BIC == np.min(BIC)))[0, 0]

K_best = K[ind_k]

print('Best number of cluster is = %d'%K_best)
u=1


u=1