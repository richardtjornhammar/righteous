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
import sys

def drop_duplicate_indices( df ):
    df_ = df.loc[~df.index.duplicated(keep='first')]
    return df_

def frac_rank_stats( X , method='ordinal' ):
    X_ = rankdata( X , method=method )/len(X)
    return(X_)

def make_group_analytes_unique( grouping_file , delimiter='\t' ):
    uniqe_grouping_file = '/'.join(grouping_file.split('/')[:-1]) + '/unique_'+grouping_file.split('/')[-1]
    with open( uniqe_grouping_file , 'w' ) as of:
        with open( grouping_file ) as input :
            for line in input :
                vline = line.replace('\n','').split(delimiter)
                gid, gdesc, analytes_ = vline[0], vline[1], list(set(vline[2:]))
                nvec = [gid,gdesc] ; [ nvec.append(a) for a in analytes_ ]
                print ( delimiter.join(nvec) , file = of )

def read_conversions(file_name) :
    gene2ens = {} ; non_unique = []
    with open( file_name , 'r' ) as infile:
        if sys.version_info[0] < 3:
            infile.next()
        else :
            next(infile)
        for line in infile:
            words = line.strip().split('\t')
            if len(words)==2 :
                ens, gene = words
                if gene in gene2ens:
                    gene2ens[gene].append(ens)
                else :
                    gene2ens[gene] = [ens]
    return gene2ens

def read_gene_ensemble_conversion(file_name):
    gene2ens = {} ; non_unique = []
    with open(file_name,'r') as infile:
        if sys.version_info[0] < 3:
            infile.next()
        else:
            next(infile)
        for line in infile:
            words = line.strip().split('\t')
            if len(words)==2 :
                ens, gene = words
                if gene in gene2ens:
                    non_unique.append((gene,ens,gene2ens[gene]))
                else :
                    gene2ens[gene] = ens
        if len(non_unique)>0:
            print(' WARNING ' )
            print( 'FOUND ', len(non_unique), ' NON UNIQUE ENTRIES' )
    return gene2ens

def create_synonyms( convert_file , unique_mapping=False ):
    # CREATE SYNONYMS
    ens2sym , sym2ens = {} , {}
    if unique_mapping:
        sym2ens = read_gene_ensemble_conversion( convert_file )
        ens2sym = { v:k for k,v in sym2ens.items() }
    else :
        sym2ens_list = read_conversions( convert_file )
        ens2sym_list = {}
        for s,L in sym2ens_list.items() :
            for e in L:
                if e in ens2sym_list:
                    ens2sym_list.append(s)
                else:
                    ens2sym_list[e] = [s]
        ens2sym = ens2sym_list
        sym2ens = sym2ens_list
    return ( ens2sym , sym2ens )

def convert_rdata_to_dataframe ( filename ) :
    #
    from rpy2.robjects import r as R
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    import rpy2.robjects as ro
    #
    print ( 'WARNING THIS PROGRAM NEED VALUE ERROR CHECKING' )
    rd_ = R.load( filename )
    if 'matrix' in str( type( R[rd_[0]] ) ).lower() :
        column_names = [ R[rd_[0]].colnames ]
        index_names  = [ R[rd_[0]].rownames ]
    else :
        column_names = [ [r for r in _rd_.colnames] for _rd_ in R[rd_[0]]]
        index_names  = [ [r for r in _rd_.rownames] for _rd_ in R[rd_[0]]]
    #
    pandas2ri.activate()
    #
    # SMALL HELPER FUNCTION THAT TRANSFORMS A RDATA OBJECT INTO
    # A PANDAS DATAFRAME. CURRENTLY THERE IS NO VALUE ERROR CHECKING
    #
    rd = R.load( filename )
    raw_df_l = []
    if 'ndarray' in str( type( R[rd[0]] ) ).lower() :
        [ raw_df_l.append( R[rd[0]] ) ]
    else :
        [ raw_df_l.append( rdf ) for rdf in ro.vectors.DataFrame(R[rd[0]]) ]
    full_df_dict = {} ; i_ = 0
    for raw_df,colnames,rownames in zip( raw_df_l,column_names,index_names ) :
        pdf = pd.DataFrame( raw_df , columns=colnames , index=rownames ) 
        full_df_dict[i_] = pdf
        i_ = i_ + 1
    pandas2ri.deactivate()
    return ( full_df_dict )

import os
if __name__ == '__main__' :
    print(' ')
