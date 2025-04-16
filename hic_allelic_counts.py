import numpy as np
import pandas as pd
import math


import sys
import os

CC1=str(sys.argv[1])
CC2=str(sys.argv[2])
CC1_m=str(sys.argv[3])
CC2_m=str(sys.argv[4])
CC1_lift=str(sys.argv[5])
CC2_lift=str(sys.argv[6])
mm10=str(sys.argv[7])
outputf=str(sys.argv[8])

cc1_uniq=pd.read_csv(CC1, header=None,sep="\t")
cc2_uniq=pd.read_csv(CC2, header=None,sep="\t")
cc1_multi=pd.read_csv(CC1_m, header=None,sep="\t")
cc2_multi=pd.read_csv(CC2_m, header=None,sep="\t")
cc1_liftover=pd.read_csv(CC1_lift, header=None,sep="\t")
cc2_liftover=pd.read_csv(CC2_lift, header=None,sep="\t")


mm10=pd.read_csv(mm10, header=0,sep="\t")
mm10.columns=['chr1','mm10_x1','mm10_x2','chr2','mm10_y1','mm10_y2','loop']
lst = pd.Series(['chr' + str(i) for i in range(1,20)])
mm10=mm10[(mm10.chr1.isin(lst))].drop_duplicates().reset_index(drop=True)

col =['chr1','CC1_x1','CC1_y1','uniq_mat_count']
cola =['chr1','CC2_x1','CC2_y1','uniq_pat_count']
col1 =['chr1','CC1_x1','CC1_y1','multi_mat_count']
col1a =['chr1','CC2_x1','CC2_y1','multi_pat_count']

cc1_uniq.columns=col
cc2_uniq.columns=cola

cc1_multi.columns=col1
cc2_multi.columns=col1a

cc1_liftover.columns = ['chr1','CC1_x1','CC1_x2','chr2','CC1_y1','CC1_y2','loop']
cc2_liftover.columns = ['chr1','CC2_x1','CC2_x2','chr2','CC2_y1','CC2_y2','loop']

merged_cc1_uniq = cc1_liftover.merge(cc1_uniq, how='inner', on=['chr1','CC1_x1','CC1_y1'])
merged_cc2_uniq = cc2_liftover.merge(cc2_uniq, how='inner', on=['chr1','CC2_x1','CC2_y1'])
merged_cc1_multi = cc1_liftover.merge(cc1_multi, how='inner', on=['chr1','CC1_x1','CC1_y1'])
merged_cc2_multi = cc2_liftover.merge(cc2_multi, how='inner', on=['chr1','CC2_x1','CC2_y1'])

merged_cc1_uniq_mm10=pd.merge(merged_cc1_uniq,mm10, how='right', on=['chr1','loop']).fillna(0)
merged_cc2_uniq_mm10=pd.merge(merged_cc2_uniq,mm10, how='right', on=['chr1','loop']).fillna(0)
merged_cc1_multi_mm10=pd.merge(merged_cc1_multi,mm10, how='right', on=['chr1','loop']).fillna(0)
merged_cc2_multi_mm10=pd.merge(merged_cc2_multi,mm10, how='right', on=['chr1','loop']).fillna(0)

merged_cc1=pd.merge(merged_cc1_uniq_mm10,merged_cc1_multi_mm10, how='inner', on=['chr1','loop','mm10_x1','mm10_x2','mm10_y1','mm10_y2','CC1_x1','CC1_x2','CC1_y1','CC1_y2']).fillna(0)
merged_cc2=pd.merge(merged_cc2_uniq_mm10,merged_cc2_multi_mm10, how='inner', on=['chr1','loop','mm10_x1','mm10_x2','mm10_y1','mm10_y2','CC2_x1','CC2_x2','CC2_y1','CC2_y2']).fillna(0)

all_merged=pd.merge(merged_cc2,merged_cc1, how='inner', on=['chr1','loop','mm10_x1','mm10_x2','mm10_y1','mm10_y2']).fillna(0)
all_merged['total_counts']=all_merged['uniq_mat_count']+all_merged['uniq_pat_count']+((all_merged['multi_mat_count']+all_merged['multi_pat_count'])/2)

all_merged=all_merged[['chr1','mm10_x1','mm10_x2','mm10_y1','mm10_y2','CC1_x1','CC1_x2',
                      'CC1_y1','CC1_y2','CC2_x1', 'CC2_x2','CC2_y1', 'CC2_y2','loop', 'uniq_mat_count',
                      'uniq_pat_count','multi_mat_count','multi_pat_count','total_counts']]

all_merged.to_csv(outputf, sep='\t', index=False, header=True)
