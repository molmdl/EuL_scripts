import numpy as np
import time

start = time.time()
data = np.loadtxt('t1/rerun_Protein_ALA_128.xvg', comments=['#', '@'])
print(f"Time taken: {time.time() - start:.3f} s, shape: {data.shape}")
