import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import warnings

print("Loading data...")
data = np.load('t1/binned_IE.npz')
mean_IE = data['mean_IE']
resids = data['resids']
bins = data['bins']

print("Data loaded. Processing...")
min_res = resids.min()
max_res = resids.max()
full_res_range = np.arange(min_res, max_res + 1)
mapped_mean_IE = np.full((len(full_res_range), mean_IE.shape[1]), np.nan)
for idx, res in enumerate(resids):
    row_idx = res - min_res
    mapped_mean_IE[row_idx, :] = mean_IE[idx, :]

extent = [bins[0], bins[-1], min_res - 0.5, max_res + 0.5]
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    vmin = np.nanmin(mapped_mean_IE)
vmax = 25.0
norm = TwoSlopeNorm(vcenter=0, vmin=vmin if vmin < 0 else -0.1, vmax=vmax)

print("Plotting...")
plt.figure(figsize=(10, 8))
im = plt.imshow(mapped_mean_IE, aspect='auto', origin='lower', extent=extent, cmap='bwr', norm=norm, interpolation='nearest')

neg_cutoff = -40
pos_cutoff = 10
labels_to_plot = []

with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=RuntimeWarning)
    for idx, res in enumerate(resids):
        row_IE = mean_IE[idx, :]
        min_val = np.nanmin(row_IE)
        if min_val < neg_cutoff:
            min_idx = np.nanargmin(row_IE)
            peak_d1 = (bins[min_idx] + bins[min_idx+1]) / 2.0
            labels_to_plot.append((peak_d1, res, f"RES{res}", min_val))

print("Annotating...")
for peak_d1, res, text, val in labels_to_plot:
    plt.plot(peak_d1, res, '*', color='yellow', markersize=8, markeredgecolor='black', markeredgewidth=0.5)
    plt.annotate(text, xy=(peak_d1, res), xytext=(peak_d1 + 0.15, res + 2),
                 color='black', fontsize=12, fontweight='bold',
                 bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec="none"),
                 arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color="black", alpha=0.7))

print("Saving...")
plt.savefig('test_anno.png')
print("Done.")
