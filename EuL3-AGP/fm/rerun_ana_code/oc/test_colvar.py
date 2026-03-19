import pandas as pd
import numpy as np

with open('t1/trj.COLVAR.o') as f:
    header = f.readline().strip().split()[2:] # Skip "#! FIELDS"
    
df = pd.read_csv('t1/trj.COLVAR.o', sep=r'\s+', comment='#', names=header)
print(df.head())
print("min d1:", df['d1'].min(), "max d1:", df['d1'].max())
