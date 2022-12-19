import mne
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import savemat

## reading data

data_raw = mne.io.read_raw_edf('reza_record_2_raw.edf')
data_raw = data_raw.get_data()

l_data = data_raw.shape[0]
i_not = [8, 16, 17, 20, 21, 24]
data_rem_mean = np.zeros((19, data_raw.shape[1]))

cntr = -1
for i in range(l_data):
    A1_ = data_raw[i, :]
    A1_ = A1_ - np.mean(A1_)


    if not i_not.count(i):
        cntr += 1
        print(cntr)
        data_rem_mean[cntr, :] = A1_

    plt.figure()
    plt.plot(A1_)
    plt.title('A1 - %d'%(i+1))
    plt.show()
    u=1

mdic = {"data": data_rem_mean}

savemat('reza_data.mat', mdic)
u=1