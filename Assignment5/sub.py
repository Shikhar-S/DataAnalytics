import pandas as pd
df=pd.read_csv('../data/Raw_Data_GeneSpring.csv', sep='\t')


def get_value(category,idx):
    #category def
    #1 male non smoker
    #2 male smoker
    #3 women non smoker
    #4 women smoker
    columns=[]
    if category==1:
        columns=df.columns[1:13]
    elif category==2:
        columns=df.columns[13:25]
    elif category==3:
        columns=df.columns[25:37]
    elif category==4:
        columns=df.columns[37:49]
    return df.iloc[idx][columns].values

import numpy as np
import scipy.stats as stats

def setup_models():
    A=np.zeros((48,4))
    for i in range(4):
        for j in range(12):
            A[j+(i*12),i]=1
            
    A_null=np.zeros((48,4))
    for i in range(4):
        for j in range(12):
            if i==0:
                A_null[j+(i*12),0]=1
                A_null[j+(i*12),1]=1
            elif i==1:
                A_null[j+(i*12),0]=1
                A_null[j+(i*12),2]=1
            elif i==2:
                A_null[j+(i*12),3]=1
                A_null[j+(i*12),1]=1
            elif i==3:
                A_null[j+(i*12),3]=1
                A_null[j+(i*12),2]=1
    return (A,A_null)

def setup_1way_model():
    B=np.zeros((24,2))
    for i in range(2):
        for j in range(12):
            B[j+i*12,i]=1
    B_null=np.zeros((24,1))+1
    return (B,B_null)

def get_sum_sq(A,h):
    k=h.shape[0]
    P=np.matmul(np.transpose(A),A)
    S=np.linalg.pinv(P)
    Inner=np.eye(k)-np.matmul(np.matmul(A,S),np.transpose(A))
    sum_sq=np.matmul(np.matmul(np.transpose(h),Inner),h)
    return sum_sq

def generate_p_value(h,A,A_null):
    denom=get_sum_sq(A,h)
    if denom<1e-8:
        return 1
    A_rank=np.linalg.matrix_rank(A)
    A_null_rank=np.linalg.matrix_rank(A_null)
    F=((get_sum_sq(A_null,h)/denom)-1)*((h.shape[0]-A_rank)/(A_rank-A_null_rank)) #(48-4/4-3)
    p=1-stats.f.cdf(F,A_rank-A_null_rank,h.shape[0]-A_rank)
    return p

def compute_gene_type(i):
    #1- down in women
    #2- up in women
    #3- down in men
    #4- up in men
    data=[]
    for category in range(1,5):
        data.append(get_value(category,i))
    means=[np.mean(x) for x in data]
    delta_men=means[0]-means[1]
    delta_women=means[2]-means[3]
    ans=[-1,-1] #-1 denotes no list
    if significant_difference(i,male=True):
        ans[0]=2 if delta_men>0 else 3
    if significant_difference(i,male=False):
        ans[1]=0 if delta_women>0 else 1
    return ans

def significant_difference(i,male=True):
    B,B_null=setup_1way_model()
    h=None
    if male:
        h=df.iloc[i][df.columns[1:25]].values
    else:
        h=df.iloc[i][df.columns[25:49]].values
    p_value_anova=generate_p_value(h,B,B_null)
    if p_value_anova<0.05:
        return True

print('Computing p values, please wait.')

P=[]
A,A_null=setup_models()
for idx,row in df.iterrows():
    h=row[df.columns[1:49]].values
    if idx%10000==0:
    	print('On row,',idx)
    P.append(generate_p_value(h,A,A_null))

import matplotlib.pyplot as plt
plt.hist(P, bins = 20)
plt.show()

print('Since histogram is clustered towards 1, no better estimate for n_0 is justifiable than n.')

shortlisted_genes=[]
for i,p in enumerate(P):
    if(p<0.05): #as instructed on forum
        shortlisted_genes.append(i)

print('shortlisted_genes count:',len(shortlisted_genes))


DNA_REPAIR=["ABL1","ALKBH1","APEX1","APTX","ASF1A","ATM","ATP23","ATR","ATRX","ATXN3","BLM","BRCA1","BRCA2","BTG2","CCNO","CDKN2D","CEBPG","CIB1","CSNK1D","CSNK1E","DDB1","DDB2","ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERCC6","ERCC8","EXO1","FANCA","FANCC","FANCG","FEN1","GADD45A","GADD45G","GTF2H1","GTF2H4","HMGB1","HMGB1P10","HMGB2","HUS1","IGHMBP2","KAT5","LIG1","LIG3","LIG4","MLH1","MMS19","MNAT1","MPG","MRE11","MSH2","MSH3","MSH5","MSH6","MUTYH","NBN","NHEJ1","NTHL1","OGG1","PARP1","PARP3","PMS1","PMS2","PMS2P1","PNKP","POLA1","POLD1","POLE","POLE2","POLG","POLH","POLI","POLL","POLQ","PRKCG","RAD1","RAD17","RAD21","RAD23A","RAD23B","RAD50","RAD51","RAD51B","RAD51C","RAD52","RAD54B","RAD54L","RAD9A","RBBP8","RECQL","RECQL4","RECQL5","REV1","RFC3","RPA1","RPAIN","RUVBL2","SETX","SMC1A","SMUG1","SOD1","SUMO1","TDG","TNP1","TP53","TP73","TREX2","UBE2A","UBE2B","UBE2N","UBE2V1","UBE2V2","UNG","UPF1","UVRAG","VCP","WRNIP1","XAB2","XPC","XRCC2","XRCC3","XRCC4","XRCC6"]
FREE_RADICAL=["ADPRHL2","APOA4","ATP7A","BMP7","CCS","CD36","DHFR","DHFRP1","ERCC6","FANCC","FBLN5","GCH1","GLRX2","MIR21","MPO","MT3","NFE2L2","NOS3","NQO1","PARK7","PRDX1","PRDX2","RGN","SOD1","SOD2","SOD3","SZT2","TNF","TXNRD2","UCP2","UCP3"]
CYTOTOXICITY=["ARAF","BID","BRAF","CASP3","CD244","CD247","CD48","CHP1","CHP2","CSF2","FAS","FASLG","FCER1G","FCGR3A","FCGR3B","FYN","GRB2","GZMB","HCST","HLA-A","HLA-B","HLA-C","HLA-E","HLA-G","HRAS","ICAM1","ICAM2","IFNA1","IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA2","IFNA21","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8","IFNAR1","IFNAR2","IFNB1","IFNG","IFNGR1","IFNGR2","ITGAL","ITGB2","KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL4","KIR2DL5A","KIR2DS1","KIR2DS3","KIR2DS4","KIR2DS5","KIR3DL1","KIR3DL2","KLRC1","KLRC2","KLRC3","KLRD1","KLRK1","KRAS","LAT","LCK","LCP2","MAP2K1","MAP2K2","MAPK1","MAPK3","MICA","MICB","NCR1","NCR2","NCR3","NFAT5","NFATC1","NFATC2","NFATC3","NFATC4","NRAS","PAK1","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2","PIK3R3","PIK3R5","PLCG1","PLCG2","PPP3CA","PPP3CB","PPP3CC","PPP3R1","PPP3R2","PRF1","PRKCA","PRKCB","PRKCG","PTK2B","PTPN11","PTPN6","RAC1","RAC2","RAC3","RAET1E","RAET1G","RAET1L","RAF1","SH2D1A","SH2D1B","SH3BP2","SHC1","SHC2","SHC3","SHC4","SOS1","SOS2","SYK","TNF","TNFRSF10A","TNFRSF10B","TNFRSF10C","TNFRSF10D","TNFSF10","TYROBP","ULBP1","ULBP2","ULBP3","VAV1","VAV2","VAV3","ZAP70"]
XENOBIOTIC=["AADAC","ACAA1","ACSL1","ACSM1","ACSM2B","ACY1","ACY3","AHR","AHRR","AIP","AKR1C1","AKR7A2","AKR7A3","AKR7L","ALDH3A1","AOC1","AOC2","AOC3","ARNT","ARNT2","AS3MT","BCHE","BPHL","CBR3","CES1","CES2","CES3","CMBL","CRYZ","CYB5B","CYB5R3","CYP1A1","CYP1A2","CYP1B1","CYP26A1","CYP26B1","CYP2A13","CYP2A6","CYP2A7","CYP2B6","CYP2C18","CYP2C19","CYP2C8","CYP2C9","CYP2D6","CYP2D7","CYP2E1","CYP2F1","CYP2G1P","CYP2J2","CYP2R1","CYP2S1","CYP2U1","CYP2W1","CYP3A4","CYP3A5","CYP3A7","CYP46A1","DPEP1","EPHX1","EPHX2","FMO1","FMO2","FMO3","GGT1","GLYAT","GRIN1","GSTA4","GSTM1","GSTM2","GSTM3","GSTM4","GSTO1","GSTO2","GSTP1","HNF4A","HSP90AB1","LPO","MARC1","MARC2","MGST1","MGST2","MGST3","N6AMT1","NAT1","NAT2","NQO1","NQO2","NR1I2","PON3","POR","PTGES3","PTGS1","RORA","RORC","S100A12","STAR","SULT1A1","SULT1A2","SULT1A3","SULT1B1","UGT1A1","UGT1A10","UGT1A3","UGT1A4","UGT1A5","UGT1A6","UGT1A7","UGT1A8","UGT1A9","UGT2B11","UGT2B15","UGT2B28"]

dna_repair=[0 for i in range(4)]
free_radical=[0 for i in range(4)]
cytotoxicity=[0 for i in range(4)]
xenobiotic=[0 for i in range(4)]


for i in shortlisted_genes:
    Gene=df.iloc[i]['GeneSymbol']
    z=compute_gene_type(i)
    flag=False
    
    if Gene in DNA_REPAIR:
        if z[0]!=-1:
            dna_repair[z[0]]+=1
        if z[1]!=-1:
            dna_repair[z[1]]+=1
        flag=True
    if Gene in FREE_RADICAL:
        if z[0]!=-1:
            free_radical[z[0]]+=1
        if z[1]!=-1:
            free_radical[z[1]]+=1
        flag=True
    if Gene in CYTOTOXICITY:
        if z[0]!=-1:
            cytotoxicity[z[0]]+=1
        if z[1]!=-1:
            cytotoxicity[z[1]]+=1
        flag=True
    if Gene in XENOBIOTIC:
        if z[0]!=-1:
            xenobiotic[z[0]]+=1
        if z[1]!=-1:
            xenobiotic[z[1]]+=1
        flag=True
    # if flag:
    #     print(z)

print('Intersection Counts- Down in women, Up in women, Down in men, Up in women')
print("DNA Repair",dna_repair)
print("Free radical response",free_radical)
print("Natural killer cell cytotoxicity",cytotoxicity)
print("Xenobiotic Metabolism",xenobiotic)