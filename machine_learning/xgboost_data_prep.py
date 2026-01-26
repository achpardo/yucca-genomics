#!/usr/bin/env python3

__author__ = "Anna Pardo"

# import modules
import argparse
import os
import pandas as pd
import numpy as np
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from functools import reduce

# define functions

def iscomplete(gt,mdonly):
    gtmd = mdonly[mdonly["genotype"]==gt]
    gtcount = gtmd.groupby(["treat","ZT"]).count().reset_index()[["treat","ZT","sample_name"]].rename(columns={"sample_name":"nreps"})
    if len(gtcount.index)==12:
        return "complete"
    elif len(gtcount.index)<12:
        return "incomplete"

def import_tpm(file):
    mdtpm = pd.read_csv(file,sep="\t",header="infer")
    mdtpm.drop(["time","species"],axis=1,inplace=True)
    mdtpm["genotype"] = mdtpm["genotype"].astype(str)
    zttpm = mdtpm[mdtpm["ZT"].isin([1.0,5.0,9.0,13.0,17.0,21.0])]
    mdonly = zttpm[["sample_name","genotype","treat","ZT"]]
    zttpm = zttpm[~zttpm["genotype"].isin(["15","52","55","70","53","16","20","50"])]
    return zttpm,mdonly

def egtnum(gid,comptpm):
    tpmsub = comptpm[["genotype",gid]]
    maxtpm = tpmsub.groupby("genotype").max().reset_index()
    # count the number of genotypes with max>0
    max0 = maxtpm[maxtpm[gid]>0]
    return len(max0)

def variance_threshold_selector(data):
    selector = VarianceThreshold()
    selector.fit(data)
    return data[data.columns[selector.get_support(indices=True)]]

def genemeans(gid,vttpm_log):
    sublog = vttpm_log[["genotype","treat_ZT",gid]]
    sublog.groupby(["genotype","treat_ZT"]).mean().reset_index()
    slp = sublog.pivot_table(index="genotype",columns="treat_ZT",values=gid).reset_index()
    for c in slp.columns:
        if c!="genotype":
            slp.rename(columns={c:gid+"_"+c},inplace=True)
    return slp

def main():
    parser = argparse.ArgumentParser(description="Parse args")
    # command line arguments
    parser.add_argument("--tpmfile","-t",type=str,help="full path to input TPM")
    parser.add_argument("--outfile","-o",type=str,help="full path to output file")
    args = parser.parse_args()
    tpmfile = str(args.tpmfile)
    outfile = str(args.outfile)

    zttpm,mdonly = import_tpm(tpmfile)

    dropgenes = []
    for c in zttpm.columns:
        if c.startswith("Y"):
            if egtnum(c,zttpm)<=2:
                dropgenes.append(c)

    filtcomp = zttpm.drop(dropgenes,axis=1)
    filttpm = filtcomp.set_index(["sample_name","genotype","treat","ZT"])
    vttpm = variance_threshold_selector(filttpm)
    vttpm_log = vttpm.apply(lambda x: np.log2(x+1))
    vttpm_log.reset_index(inplace=True)
    vttpm_log.drop("sample_name",axis=1,inplace=True)
    vttpm_log["treat_ZT"] = vttpm_log["treat"]+"_"+vttpm_log["ZT"].astype(int).astype(str)
    vttpm_log.drop(["treat","ZT"],axis=1,inplace=True)

    dflist = []
    for c in vttpm_log.columns:
        if c.startswith("Y"):
            dflist.append(genemeans(c,vttpm_log))

    allfeatures = reduce(lambda x,y: pd.merge(x,y,on="genotype"),dflist)

    allfeatures.to_csv(outfile,sep="\t",header=True,index=False)

if __name__ == "__main__":
    main()
