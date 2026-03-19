import pandas as pd
import numpy as np

# Let's read the saved dataframe to see what columns are in frame_weights.csv
df_calc = pd.read_csv('t1/frame_weights.csv')

# And read PLUMED test output
with open('weights.dat', 'r') as f:
    fields = f.readline().strip().split()[2:]
df_plumed = pd.read_csv('weights.dat', sep=r'\s+', comment='#', names=fields)

print("Our frame_weights.csv columns:", df_calc.columns.tolist())
print("PLUMED weights.dat columns:", df_plumed.columns.tolist())

# Our weight was: w = exp(bias/kt - max(bias/kt))
kt = 2.494339
# PLUMED's "as" output from REWEIGHT_METAD is exactly bias/kt
plumed_as = df_plumed['as'].values

# The actual weight used in HISTOGRAM LOGWEIGHTS=as is exp(as)
# To avoid overflow, PLUMED subtracts max(as) before exponentiating internally
# which is exactly what we did!

# Let's compare our w to exp(plumed_as - max(plumed_as))
min_len = min(len(df_calc), len(plumed_as))
our_w = df_calc['weight'].values[:min_len]

plumed_w_normalized = np.exp(plumed_as[:min_len] - np.max(plumed_as))

diff = np.abs(our_w - plumed_w_normalized)

print(f"\nMax difference between our normalized weights and PLUMED shifted exp(as): {np.max(diff):.6e}")

print("\nSample check (first 5 frames):")
for i in range(5):
    print(f"Frame {i}: our_w = {our_w[i]:.6e}, plumed_w (normalized) = {plumed_w_normalized[i]:.6e}")

if np.max(diff) < 1e-5:
    print("\nSUCCESS: Calculated weights exactly match PLUMED REWEIGHT_METAD logic.")
else:
    print("\nWARNING: Weights do not match.")
