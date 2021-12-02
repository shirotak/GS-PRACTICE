#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from sklearn.neighbors import KNeighborsClassifier as KN
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.linear_model import LogisticRegression as LR
import joblib

# specify library directory
script_dir=os.path.dirname(os.path.abspath(__file__))
lib_dir=script_dir+"/data/"
out_dir=script_dir+"/data/"

def main():
    # import tcga data & labels
    with open(lib_dir+"TCGA_7181_cluster_numbers.tsv") as f:
        cluster_numbers=f.read().splitlines()
    with open(lib_dir+"TCGA_7181_cluster_colors.tsv") as f:
        cluster_colors=f.read().splitlines()    
    with open(lib_dir+"TCGA_7181_cluster_annotations.tsv") as f:
        cluster_anno=f.read().splitlines()
    df_mutsig=pd.read_csv(lib_dir+"TCGA_7181_df_mutsig_log10.tsv",sep="\t",index_col=0)
    X=df_mutsig
    Y=cluster_anno

    # make four classifiers
    ## KNeighborsClassifier 
    dict_kn=dict( n_neighbors=5,
        weights='distance',
        algorithm='auto',
        leaf_size=30,
        p=2,
        metric='minkowski',
        metric_params=None,
        n_jobs=None)
    kn=KN(**dict_kn)
    print("KNN parameters")
    print(["%s : %s"%(k,v) for k,v in dict_kn.items()] )

    ## SVC
    dict_svc=dict( C=1.0, 
        kernel='rbf',
        degree=3, 
        gamma=0.1, 
        coef0=0.0, 
        shrinking=True,
        probability=False, 
        tol=0.001, 
        cache_size=200, 
        class_weight=None, 
        max_iter=-1,
        decision_function_shape='ovr', 
        break_ties=False,
        random_state=None)
    svc=SVC(**dict_svc)
    print("SVC parameters")
    print(["%s : %s"%(k,v) for k,v in dict_svc.items()] )

    ## RandomForestClassifier
    dict_rfc=dict(n_estimators=100,
        criterion='gini',
        max_depth= None,
        min_samples_split= 2,
        min_samples_leaf= 1,
        min_weight_fraction_leaf= 0.0,
        max_features= 'auto',
        max_leaf_nodes= None,
        min_impurity_decrease= 0.0,
        bootstrap= True,
        oob_score= False,
        n_jobs= None,
        random_state= None,
        class_weight= None,
        ccp_alpha= 0.0,
        max_samples= None)
    rfc=RFC(**dict_rfc)
    print("RFC parameters")
    print(["%s : %s"%(k,v) for k,v in dict_rfc.items()] )

    ## LogisticRegression
    dict_lr=dict(C=1,
        class_weight=None,
        dual=False,
        fit_intercept=True,
        intercept_scaling=1,
        l1_ratio=None,
        max_iter=2000,
        multi_class='auto',
        n_jobs=None,
        penalty='l2',
        random_state=None,
        solver='lbfgs',
        tol=0.0001,
        verbose=0,
        warm_start=False)
    lr=LR(**dict_lr)
    print("LR parameters")
    print(["%s : %s"%(k,v) for k,v in dict_lr.items()] )

    ## write out
    print("Write out four classifiers...")
    names=["KNN","SVC","RFC","LRC"]
    clfs=[kn,svc,rfc,lr]
    for clf,name in zip(clfs,names):
        clf.fit(X,Y)
        with open(out_dir+name+".joblib", mode="wb") as fw:
            joblib.dump(clf,fw) 

    # UMAP Projection
    import umap
    dict_umap=dict(n_neighbors=10,
        min_dist=0.35,
        n_components=2,
        random_state=1,
        metric="euclidean")
    print( "UMAP parameters")
    print(["%s : %s"%(k,v) for k,v in dict_umap.items()] )
    fit = umap.UMAP(**dict_umap)
    u= fit.fit(df_mutsig)
    print("Write out umap projector...")
    joblib.dump(u,out_dir+"UMAP_projector.joblib")
    print("Finished")

if __name__ == "__main__":
    main()

