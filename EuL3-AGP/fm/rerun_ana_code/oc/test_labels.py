import numpy as np

data = np.load('t1/binned_IE.npz')
mean_IE = data['mean_IE']
resids = data['resids']

neg_cutoff = -40
pos_cutoff = 5

for i, res in enumerate(resids):
    min_val = np.nanmin(mean_IE[i])
    max_val = np.nanmax(mean_IE[i])
    if min_val < neg_cutoff:
        print(f"Res {res} neg peak: {min_val:.1f}")
    if max_val > pos_cutoff:
        print(f"Res {res} pos peak: {max_val:.1f}")
