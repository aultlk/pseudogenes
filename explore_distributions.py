'''
   This script groups RPM distributions by locus, locus bin, and rearrangement type. It 
   then gathers statistics and structures it to be visually appealing for quick analysis.
'''


from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

ROOT_DIR = Path(__file__).resolve().parent.parent
APPENDED_META = f'{ROOT_DIR}/appended_meta.csv'


meta = pd.read_csv(APPENDED_META)

counts = (meta
          .groupby(['locus', 'locus_bin', 'rearrangement'])
          [['RPM']]
          .describe()
          .stack(level=0)
          .loc[:, ['count', 'min', 'mean', 'max']]
          .astype(int))
#          .applymap(lambda x: '{:,}'.format(x)))

counts.index = counts.index.droplevel(3)

counts.to_csv(f'{ROOT_DIR}/distribution_stats.csv')
