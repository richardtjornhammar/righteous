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

def run_mofa_data():
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


def show_results ( qdf , sort_on , show_info=['description'] , nitems=15 ) :
    for label in sort_on :
        print ( qdf.sort_values(label)[[*show_info,*sort_on]][:nitems] )


#examples

def run_tcga_brca():
    #
    # READ THE TCGA DATA
    results = []
    analyte_df =                pd.read_csv( './data/analyte_TCGABRCA.csv', '\t' , index_col=0 )
    journal_df = prune_journal( pd.read_csv( './data/journal_TCGABRCA.csv', '\t' , index_col=0 ) )
    journal_df = journal_df .loc[ :,analyte_df.columns.values ].copy()
    #
    # ONLY KEEP WELL DETERMINED CASES
    test = ' C(HER2pos) + C(ERpos) + C(PRpos)' 
    bSel = [ not ('ND' in erp or 'ND' in her2p or 'ND' in prp) for ( erp,her2p,prp ) in zip ( 
                  journal_df.loc['ERposLabel'] , journal_df.loc['HER2posLabel'] , journal_df.loc['PRposLabel'] 
            ) ]
    #
    adf = analyte_df.iloc[:,bSel].copy().apply( lambda x:np.log2(x+1.) )
    jdf = journal_df.iloc[:,bSel].copy()
    #
    # CONDUCT THE RIGHTEOUS QUANTIFICATION
    qdf = quantify ( adf , jdf ,
                    'anova ~' + test , './data/vanilla_reactome.gmt',
                     test_type = 'random' )
    qdf.sort_values(test.replace(' ','').split('+')[0]+',q').to_csv('./tcga_her2p_erp_prp_results.tsv','\t')
    results.append(qdf)

    print ( 'SHOWING:: P ~ ' , test )
    print ( qdf.sort_values('C(HER2pos),q')['description'][:5] )
    print ( qdf.sort_values('C(ERpos),q' ) ['description'][:5] )
    print ( qdf.sort_values('C(PRpos),q')['description'][:5] )
    #
    test = 'C(TripleNegative)' 
    bSel = [ not ('ND' in lab ) for lab in journal_df.loc['TripleNegativeLabel'] ]
    #
    adf = analyte_df.iloc[:,bSel].copy().apply( lambda x:np.log2(x+1.) )
    jdf = journal_df.iloc[:,bSel].copy()
    #
    # CONDUCT THE RIGHTEOUS QUANTIFICATION
    qdf = quantify ( adf , jdf ,
               'anova ~' + test , './data/reactome.gmt' ,
                test_type = 'random' )
    qdf.sort_values(test.replace(' ','').split('+')[0]+',q').to_csv('./tcga_tripleneg_results.tsv','\t')
    results.append(qdf)

    print ( 'SHOWING:: P ~ ' , test )
    print ( qdf.sort_values(test+',q')['description'][:5] )
    test = 'C(HER2pos)+C(PRpos)+C(HER2pos):C(PRpos)' 
    bSel = [ not ('ND' in her2p or 'ND' in prp) for ( her2p,prp ) in zip ( 
                       journal_df.loc['HER2posLabel'] , journal_df.loc['PRposLabel'] 
             ) ]
    #
    adf = analyte_df.iloc[:,bSel].copy().apply( lambda x:np.log2(x+1.) )
    jdf = journal_df.iloc[:,bSel].copy()
    #
    # CONDUCT THE RIGHTEOUS QUANTIFICATION
    qdf = quantify ( adf , jdf ,
                    'anova ~' + test , './data/reactome.gmt' ,
                     test_type = 'random' )
    qdf.sort_values(test.replace(' ','').split('+')[0]+',q').to_csv('./tcga_her2pvsprp_prp_results.tsv','\t')
    results.append(qdf)
    print ( 'SHOWING:: P ~ ' , test )
    print ( qdf.sort_values('C(HER2pos):C(PRpos)'+',q')['description'][:5] )
    return ( results )

def run_CLL():
    #
    analyte_df =                pd.read_csv('./data/analyte_MOFACLL.csv', '\t' , index_col=0 )
    journal_df = prune_journal( pd.read_csv('./data/journal_MOFACLL.csv', '\t' , index_col=0 ) )
    journal_df = journal_df.loc[:,analyte_df.columns.values].copy()
    #
    test = 'Ibrutinib4'
    qdf = quantify ( analyte_df , journal_df , 'anova~'+test ,
                         'data/vanilla_reactome.gmt',
                         test_type='random' )
    qdf.sort_values(test+',q').to_csv('./cll_ighv_results.tsv','\t')
    return(qdf)


def run_MS():
    analyte_df =                 pd.read_csv ( './data/analyte_MS.csv'  , '\t' , index_col=0 )
    journal_df = prune_journal ( pd.read_csv ( './data/journal_MS.csv', '\t' , index_col=0 ) )
    journal_df = journal_df .loc[ :,analyte_df.columns.values ].copy( )
    #
    # ONLY KEEP WELL DETERMINED CASES
    bSel = [ not ('ND' in stat ) for ( stat ) in journal_df.loc['StatusLabel'] ]
    test = 'C(Status)'
    #
    adf = drop_duplicate_indices( analyte_df.iloc[:,bSel].copy() )
    jdf = journal_df.iloc[:,bSel].copy()
    #
    # CONDUCT THE RIGHTEOUS QUANTIFICATION
    qdf = quantify ( adf , jdf ,
                     'anova ~' + test , './data/vanilla_reactome.gmt' ,
                      test_type = 'random' )
    show_results( qdf,sort_on=['C(Status),q'],show_info=['description','C(Status),p'] )
    return ( qdf )


if __name__ == '__main__' :
    #
    bMS = True
    if bMS :
        qdf = run_MS()

    bTCGA = True
    if bTCGA :
        run_tcga_brca()

    bMOFA = True
    if bMOFA :
        run_CLL()

