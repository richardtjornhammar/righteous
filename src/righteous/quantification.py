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
import numpy as np
import pandas as pd

def qvalues( pvalues, pi0=1 ):
    m = float(len(pvalues))
    assert(m>0)
    pvalues.sort()
    num_p, p_sum, qs = 0, 0.0, []
    for p in pvalues:
        num_p += 1
        p_sum += p
        q = pi0*p*(m/float(num_p))
        qs.append((q,p))
    qs.reverse()
    old_q = 1.0
    for ix in range(len(qs)):
        q = min(old_q,qs[ix][0])
        old_q  = q
        qs[ix] = (q,qs[ix][1])
    qs.reverse()
    return qs

from scipy import stats
from statsmodels.stats.anova import anova_lm as anova
import statsmodels.api as sm
import patsy
def anova_test ( formula, group_expression_df, journal_df, test_type = 'random' ) :
    type_d = { 'paired':1 , 'random':2 , 'fixed':1 }
    formula = formula.replace(' ','')
    tmp_df = pd.concat([ journal_df, group_expression_df ])
    gname = tmp_df.index.tolist()[-1]
    formula_l = formula.split('~')
    rename = { gname:formula_l[0] }
    tmp_df.rename( index=rename, inplace=True )
    tdf = tmp_df.T.iloc[ :,[ col in formula for col in tmp_df.T.columns] ].apply( pd.to_numeric )
    y, X = patsy.dmatrices( formula, tdf, return_type='dataframe')
    model = sm.OLS(endog=y,exog=X).fit()
    model .model.data.design_info = X.design_info
    table = sm.stats.anova_lm(model,typ=type_d[test_type])
    return table.iloc[ [(idx in formula) for idx in table.index],-1]

def parse_test ( statistical_formula, group_expression_df , journal_df , test_type = 'random' ) :
    result = anova_test( statistical_formula, group_expression_df , journal_df , test_type=test_type )
    return ( result )

def prune_journal ( journal_df , remove_units_on = '_' ) :
    journal_df = journal_df.loc[ [ 'label' in idx.lower() or '[' in idx for idx in journal_df.index.values] , : ].copy()
    bSel = [ ('label' in idx.lower() ) for idx in journal_df.index.values]
    bool_dict = { False:0 , True:1 , 'False':0 , 'True':1 }
    str_journal = journal_df.iloc[ bSel ]
    journal_df = journal_df.replace({'ND':np.nan})
    nmr_journal = journal_df.iloc[ [ not b for b in bSel ] ].replace(bool_dict).apply( pd.to_numeric )
    if not remove_units_on is None :
        nmr_journal.index = [ idx.split(remove_units_on)[0] for idx in nmr_journal.index ]
    journal_df = pd.concat( [nmr_journal,str_journal] )
    return( journal_df )

class RCA( object ) :
    def __init__(self):
        self.components_ = None
        self.F_ = None
        self.U_ , self.S_, self.V_ = None,None,None
        self.evr_ = None
        self.var_ = None

    def fit_transform(self,X):
        Xc = X - np.mean( X , 0 )
        u, s, v = np.linalg.svd( Xc, full_matrices=False )
        S = np.diag( s )
        self.F_ = np.dot(u,S)
        self.var_ = s ** 2 / Xc.shape[0]
        self.explained_variance_ratio_ = self.var_/self.var_.sum()
        self.U_, self.S_, self.V_ = u,s,v
        self.components_ = self.V_
        return(self.F_)

def quantify ( analyte_df , journal_df , formula , grouping_file , synonyms = None ,
                delimiter = '\t' , test_type = 'random' ,
                split_id = None , skip_line_char = '#' 
              ) :
    dimred = RCA()
    statistical_formula = formula
    if not split_id is None :
        nidx = [ idx.split(split_id)[-1].replace(' ','') for idx in analyte_df.index.values ]
        analyte_df.index = nidx
    sidx = set( analyte_df.index.values ) ; nidx=len(sidx)
    eval_df = None
    with open ( grouping_file ) as input:
        for line in input:
            if line[0] == skip_line_char :
                continue
            vline = line.replace('\n','').split(delimiter)
            gid,gdesc,analytes_ = vline[0],vline[1],vline[2:]
            if not synonyms is None :
                [ analytes_.append(synonyms[a]) for a in analytes_ if a in synonyms ]
            try :
                group = analyte_df.loc[[a for a in analytes_ if a in sidx] ].dropna( axis=0, how='any', thresh=analyte_df.shape[1]/2 ).drop_duplicates()
            except KeyError as e :
                continue
            L_ = len( group ); str_analytes=','.join(group.index.values)
            if L_>0 :
                Xnew = dimred.fit_transform(group.T.values)
                group_expression_df = pd.DataFrame([Xnew.T[0]],columns=analyte_df.columns.values,index=[gid])
                rdf = pd.DataFrame( parse_test( statistical_formula, group_expression_df , journal_df , test_type=test_type )).T
                rdf .columns = [ col+',p' if (not ',s' in col) else col+',s' for col in rdf.columns ]
                rdf['description'] = gdesc+','+str(L_)
                rdf['analytes'] = str_analytes
                rdf.index = [ gid ] ; ndf = pd.concat([rdf.T,group_expression_df.T]).T
                if eval_df is None :
                    eval_df = ndf
                else :
                    eval_df = pd.concat([eval_df,ndf])
    edf = eval_df.T
    for col in eval_df.columns :
        if ',p' in col :
            q = [q_[0] for q_ in qvalues(eval_df.loc[:,col].values)]; l=col.split(',')[0]+',q'
            edf.loc[l] = q
    return ( edf.T )

def add_spearmanr( analyte_results_df, journal_df, what='M') :
    if what in set( journal_df.index.values ) :
        from scipy.stats import spearmanr
        K = []
        patients = [ c for c in analyte_results_df.columns if '_' in c ]
        for idx in analyte_results_df.index :
            y = journal_df.loc[what,patients].values
            x = analyte_results_df.loc[[idx],patients].values[0] # IF DUPLICATE GET FIRST
            k = spearmanr( x,y )
            K .append( k )
        analyte_results_df['Spearman'] = K
    return ( analyte_results_df )

if __name__ == '__main__' :
    #
    pvs = [0.00001,0.01,0.0002,0.00005,0.01,0.1,0.2,0.4,0.5,0.6,0.7,0.8,0.9,0.99,0.0114,0.15,0.23,0.20]
    print ( [q for q in qvalues(pvs) ] )

    """
    path_ = './'
    analyte_file  = path_ + 'fine.txt'
    journal_file  = path_ + 'coarse.txt'
    grouping_file = path_ + 'groups.gmt'
    analyte_df = pd.read_csv(analyte_file,'\t' , index_col=0 )
    journal_df = prune_journal( pd.read_csv(journal_file,'\t', index_col=0 ) )
    print ( quantify( analyte_df, journal_df, 'Group ~ Var + C(Cat) ', grouping_file ) )
    """
