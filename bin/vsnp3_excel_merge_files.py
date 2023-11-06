#!/usr/bin/env python

import glob
import pandas as pd
import time
from datetime import datetime

import locale
# Set the locale to United States
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

'''
working directory to contain Excel files to simply merge based on header to a single file
'''

ts = time.time()
st = datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')

files = glob.glob("*.xlsx")
df_list = []
for fff in files:
    df_list.append(pd.read_excel(fff, index_col='sample'))

result = pd.concat(df_list, sort=False)
result.to_excel(f'combined_excelworksheets-{st}.xlsx')
