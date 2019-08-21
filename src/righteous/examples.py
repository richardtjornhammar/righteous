"""
Copyright 2019 RICHARD TJÖRNHAMMAR

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
import pandas as pd
import numpy as np
import sys,os
from righteous.convert import *
from righteous.quantification import *

def run_mofa():
    if bMOFA_data :
        # os.system('mkdir ../data')
        # os.system('wget https://github.com/bioFAM/MOFAdata/blob/master/data/CLL_data.RData')
        # os.system('mv CLL_data.RData ../data/.')
        df_dict = convert_rdata_to_dataframe( filename = '../data/CLL_data.RData' )
        pruned_df = df_dict[2].T.dropna().T
        journal_df = df_dict[3].loc[['IGHV'],:]
        mask = [ v>=0 for v in journal_df.loc['IGHV'].values ]
        journal_df = journal_df.iloc[ :,mask ]
        use = list(set(pruned_df.columns)&set(journal_df.columns))
        analyte_df =  pruned_df .loc[ :,use ].apply( lambda x:np.log2(1+x) )
        journal_df = journal_df .loc[ :,use ]
        qdf = quantify ( analyte_df , journal_df , 'anova~C(IGHV)' , '~/Work/data/naming_and_annotations/reactome/reactome.gmt' )

