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
mm10_loc=[]
mm10_loc=mm10[['chr1','x1','x2','y1','y2']]
mm10_loc.columns=['chr','mm10_x1','mm10_x2','mm10_y1','mm10_y2']
lst = pd.Series(['chr' + str(i) for i in range(1,20)])
mm10_loc=mm10_loc[(mm10_loc.chr.isin(lst))].drop_duplicates().reset_index(drop=True)


col =['chr','CC1_x1','CC1_y1','uniq_mat_count']
cola =['chr','CC2_x1','CC2_y1','uniq_pat_count']
col1 =['chr','CC1_x1','CC1_y1','multi_mat_count']
col1a =['chr','CC2_x1','CC2_y1','multi_pat_count']

cc1_uniq.columns=col
cc2_uniq.columns=cola

cc1_multi.columns=col1
cc2_multi.columns=col1a

cc1_liftover.columns = ['chr','mm10_start','mm10_end','start','end']
cc2_liftover.columns = ['chr','mm10_start','mm10_end','start','end']

cc1_start=(cc1_liftover['start'] / 10000).apply(np.floor)*10000
cc1_end=(cc1_liftover['end'] / 10000).apply(np.floor)*10000
cc1_mm10_start=(cc1_liftover['mm10_start'] / 10000).apply(np.floor)*10000
cc1_mm10_end=(cc1_liftover['mm10_end'] / 10000).apply(np.floor)*10000

cc1_liftover_bin=pd.DataFrame(columns=['chr','mm10_x1','mm10_x2','mm10_y1','mm10_y2','CC1_x1','CC1_x2','CC1_y1','CC1_y2'])

cc2_start=(cc2_liftover['start'] / 10000).apply(np.floor)*10000
cc2_end=(cc2_liftover['end'] / 10000).apply(np.floor)*10000
cc2_mm10_start=(cc2_liftover['mm10_start'] / 10000).apply(np.floor)*10000
cc2_mm10_end=(cc2_liftover['mm10_end'] / 10000).apply(np.floor)*10000
cc2_liftover_bin=pd.DataFrame(columns=['chr','mm10_x1','mm10_x2','mm10_y1','mm10_y2','CC2_x1','CC2_x2','CC2_y1','CC2_y2'])

cc1_liftover_bin['chr']=cc1_liftover['chr']
cc1_liftover_bin['CC1_x1']=cc1_start
cc1_liftover_bin['CC1_x2']=cc1_start+10000
cc1_liftover_bin['CC1_y1']=cc1_end
cc1_liftover_bin['CC1_y2']=cc1_end+10000
cc1_liftover_bin['mm10_x1']=cc1_mm10_start
cc1_liftover_bin['mm10_x2']=cc1_mm10_start+10000
cc1_liftover_bin['mm10_y1']=cc1_mm10_end
cc1_liftover_bin['mm10_y2']=cc1_mm10_end+10000

cc2_liftover_bin['chr']=cc2_liftover['chr']
cc2_liftover_bin['CC2_x1']=cc2_start
cc2_liftover_bin['CC2_x2']=cc2_start+10000
cc2_liftover_bin['CC2_y1']=cc2_end
cc2_liftover_bin['CC2_y2']=cc2_end+10000
cc2_liftover_bin['mm10_x1']=cc2_mm10_start
cc2_liftover_bin['mm10_x2']=cc2_mm10_start+10000
cc2_liftover_bin['mm10_y1']=cc2_mm10_end
cc2_liftover_bin['mm10_y2']=cc2_mm10_end+10000


merged_cc1_uniq = cc1_liftover_bin.merge(cc1_uniq, how='inner', on=['chr','CC1_x1','CC1_y1'])
merged_cc2_uniq = cc2_liftover_bin.merge(cc2_uniq, how='inner', on=['chr','CC2_x1','CC2_y1'])
merged_cc1_multi = cc1_liftover_bin.merge(cc1_multi, how='inner', on=['chr','CC1_x1','CC1_y1'])
merged_cc2_multi = cc2_liftover_bin.merge(cc2_multi, how='inner', on=['chr','CC2_x1','CC2_y1'])

merged_cc1=pd.merge(merged_cc1_uniq,merged_cc1_multi, how='outer', on=['chr','mm10_x1','mm10_x2','mm10_y1','mm10_y2','CC1_x1','CC1_x2','CC1_y1','CC1_y2']).fillna(0)
merged_cc2=pd.merge(merged_cc2_uniq,merged_cc2_multi, how='outer', on=['chr','mm10_x1','mm10_x2','mm10_y1','mm10_y2','CC2_x1','CC2_x2','CC2_y1','CC2_y2']).fillna(0)

all_merged=pd.merge(merged_cc1.drop_duplicates(['chr','mm10_x1','mm10_x2','mm10_y1','mm10_y2']),merged_cc2.drop_duplicates(['chr','mm10_x1','mm10_x2','mm10_y1','mm10_y2']), how='outer', on=['chr','mm10_x1','mm10_x2','mm10_y1','mm10_y2']).fillna(0)
mm10_merged=mm10_loc.merge(all_merged.drop_duplicates(['chr','mm10_x1','mm10_x2','mm10_y1','mm10_y2']), how='left', on=['chr','mm10_x1','mm10_x2','mm10_y1','mm10_y2']).fillna(0)
mm10_merged['total_counts']=mm10_merged['uniq_mat_count']+mm10_merged['uniq_pat_count']+((mm10_merged['multi_mat_count']+mm10_merged['multi_pat_count'])/2)
mm10_merged[['mm10_x1','mm10_x2','mm10_y1','mm10_y2','CC1_x1','CC1_x2','CC1_y1','CC1_y2','CC2_x1','CC2_x2','CC2_y1','CC2_y2',
                    'uniq_mat_count','uniq_pat_count','multi_mat_count','multi_pat_count','total_counts']]=mm10_merged[['mm10_x1','mm10_x2','mm10_y1','mm10_y2','CC1_x1','CC1_x2','CC1_y1','CC1_y2','CC2_x1','CC2_x2','CC2_y1','CC2_y2',
                    'uniq_mat_count','uniq_pat_count','multi_mat_count','multi_pat_count','total_counts']].astype(int)
mm10_merged=mm10_merged[['chr','mm10_x1','mm10_x2','mm10_y1','mm10_y2','CC1_x1','CC1_x2','CC1_y1','CC1_y2','CC2_x1','CC2_x2','CC2_y1','CC2_y2',
                    'uniq_mat_count','uniq_pat_count','multi_mat_count','multi_pat_count','total_counts']]
mm10_merged.to_csv(outputf, sep='\t', index=False, header=True)
