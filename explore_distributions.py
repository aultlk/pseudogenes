'''
   This script groups RPM distributions by locus, locus bin, and rearrangement type. It 
   then gathers statistics and structures it to be visually appealing for quick analysis.
'''

import pandas as pd
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent
APPENDED_META = f'{ROOT_DIR}/appended_meta.csv'


meta = pd.read_csv(APPENDED_META)
counts = (meta
          .groupby(['locus', 'locus_bin', 'rearrangement'])
          [['RPM']]
          .describe()
          .stack(level=0)
          [['count', 'min', 'mean', 'max']]
          .applymap(lambda x : "{:,}".format(x)))

counts.to_csv(f'{ROOT_DIR}/distribution_stats.csv')
