import pandas as pd
import numpy as np

# Load our calculated weights
df_calc = pd.read_csv('t1/frame_weights.csv')

# Load PLUMED weights
with open('weights.dat', 'r') as f:
    fields = f.readline().strip().split()[2:]
df_plumed = pd.read_csv('weights.dat', sep=r'\s+', comment='#', names=fields)

# Align frames
min_frames = min(len(df_calc), len(df_plumed))
calc_w = df_calc['weight'].values[:min_frames]
plumed_w = df_plumed['as'].values[:min_frames]

# PLUMED's REWEIGHT_METAD returns as = bias/kT
# While we calculated w = exp(bias/kT - max(bias/kT))
# Let's check if our exponents match PLUMED's "as" output
kt = 2.494339
metad_bias = df_calc['metad.bias'] if 'metad.bias' in df_calc.columns else df_plumed['m.bias'].values[:min_frames]
expected_as = metad_bias / kt

diff_as = np.abs(expected_as - plumed_w)
max_diff = np.max(diff_as)

print(f"Max difference between calculated bias/kT and PLUMED 'as': {max_diff:.6e}")
print(f"Mean PLUMED 'as': {np.mean(plumed_w):.6f}")

# To match exactly what REWEIGHT_METAD provides to HISTOGRAM (which takes exponent internally for LOGWEIGHTS)
# HISTOGRAM with LOGWEIGHTS=as computes weight = exp(as)
# Since exp(as) can overflow, PLUMED subtracts the max before exponentiating
# Let's verify our final weights are proportional to what PLUMED would use
plumed_logweights = plumed_w
plumed_shifted = plumed_logweights - np.max(plumed_logweights)
plumed_final_weights = np.exp(plumed_shifted)

diff_w = np.abs(calc_w - plumed_final_weights)
print(f"Max difference between our normalized weights and PLUMED shifted exp(as): {np.max(diff_w):.6e}")

# Print a few samples to show they match
print("\nSample check (first 5 frames):")
for i in range(5):
    print(f"Frame {i}: calc_w = {calc_w[i]:.6e}, plumed_w (normalized) = {plumed_final_weights[i]:.6e}")

if np.max(diff_w) < 1e-5:
    print("\nSUCCESS: Calculated weights exactly match PLUMED REWEIGHT_METAD logic.")
else:
    print("\nWARNING: Weights do not match.")
