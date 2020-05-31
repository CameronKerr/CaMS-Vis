import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def compliment(str):
    return str.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()

def data_input(vcf):
    f = pd.read_table(vcf)
    f["ALT"] = f["ALT"].str[0]
    
    h = 1
    for i in f["ALT"]:
        if h == (len(f["ALT"])-1):
            break
        if f["ALT"][h] != "<":
            f["ALT"][h] = f["REF"][h-1] + f["ALT"][h] + f["REF"][h+1]
        h += 1
    g = "global"
    g = f[f["ALT"]!= "<"]
    g.drop(g.columns.difference(['REF','ALT']), 1, inplace=True)
    
def mut_catalog(data):
    reference = g["REF"].tolist()
    substitution = g["ALT"].tolist()
    basic_mutation_type = ['CA','CG','CT','TA','TC','TG']
    types = []
    codon_type = []

    for i in range(len(reference)):
        if reference[i] == 'A' or reference[i] == 'G':
            reference[i] = compliment(reference[i])
            substitution[i] = compliment(substitution[i])
        types.append(reference[i]+substitution[i][1])
        codon_type.append(substitution[i][0]+'X'+substitution[i][2])

    d = {'Reference':reference, 'Substitution':substitution, 'Type':types, 'Codon':codon_type}
    df1 = pd.DataFrame(data=d)
    
    N_rows = len(set(codon_type))
    d = {'CA':np.zeros(N_rows),'CT':np.zeros(N_rows),'CG':np.zeros(N_rows),'TA':np.zeros(N_rows),'TC':np.zeros(N_rows),'TG':np.zeros(N_rows)}
    df2 = "global"
    df2 = pd.DataFrame(data=d)
    df2.index = set(codon_type)

    for typ in df2.columns:
        for cod in df2.index:
            lis=[]
            h=0
            for sub in df1["Codon"]:
                if sub==cod and df1["Type"][h]==typ:
                    lis.append(sub)
                h+=1
            df2[typ][cod] = len(lis)
    df2.rename(columns={"CA":"C>A", "CT":"C>T", "CG":"C>G", "TA":"T>A", "TC":"T>C", "TG":"T>G"})

def mutsig_plt():
    axes = df2.plot.bar(rot=90, subplots=True, legend=None, ylim=(0,1000), 
                    layout=(1,6), figsize=(20,3), sharey=True,
                   fontsize=10)
    plt.tight_layout()
    plt.savefig('mut_sig.png')

def main(vcf):
    data_input(vcf)
    mut_catalog(g)
    mutsig_plt()

if __name__ == "__main__":
    vcf = input("Enter your bcf in text file format:")
    main(vcf)    

