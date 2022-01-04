import numpy as np
import pandas as pd
import regex
import re
from Bio import SeqIO
import Bio

###CHANGE PAM AND WINDOW INFO HERE###
PAM = 'NNNNCN'
windowstart = 6
windowend = 13
###CHANGE PAM AND WINDOW INFO HERE###

def RC(seq):
	encoder = {'A':'T','T':'A','C':'G','G':'C','N':'N','R':'Y','Y':'R', 'M':'K', 'K':'M', 'S':'S', 'W':'W', 'H':'D', 'B':'V', 'V':'B', 'D':'H'}
	rc = []
	for n in reversed(seq):
		rc.append(encoder[n])
	return ''.join(rc)

def create_PAM(pam):
	encoder = {'A':'A','T':'T','G':'G','C':'C','R':'[A|G]','Y':'[C|T]','N':'[A|T|C|G]','M':'[A|C]','K':'[G|T]','S':'[C|G]','W':'[A|T]','H':'[A|C|T]','B':'[C|G|T]','V':'[A|C|G]','D':'[A|G|T]'}
	enc_pam = {'f':[],'r':[]}
	rc_pam = RC(pam)
	for n,m in zip(pam, rc_pam):
		enc_pam['f'].append(encoder[n])
		enc_pam['r'].append(encoder[m])
	enc_pam['f'] = ''.join(enc_pam['f'])
	enc_pam['r'] = ''.join(enc_pam['r'])
	return enc_pam

enc_pam = create_PAM(PAM)


ClinVar=pd.read_csv('Clinvar.csv', sep='\t',encoding = "ISO-8859-1", low_memory=False)
Phenotypes=pd.read_csv('DiseaseNames.csv', sep='\t')
PhenotypeDict=dict(zip(Phenotypes.CUI, Phenotypes.name))

#open flanking sequence fasta files for all Y-type pathogenic human SNPs (includes both C>T and T>C ref>variant)
#downloaded as fasta file from flanks.py
handle = open("Rfasta.txt", 'r')
flanks={}
#save as a dictionary keyed on rsID as an Integer with values being 25nt of flanking sequence on each side of the SNP
for record in SeqIO.parse(handle, "fasta") :
    flanks[int(record.id)]=regex.findall('.{25}[^A,T,C,G].{25}', str(record.seq))
    #flanks[int(record.id)]= str(record.seq)
handle.close()
#merge flanking sequences to the CtoT frame on rsID
F=pd.DataFrame({'RS# (dbSNP)': list(flanks.keys()), 'Flanks': [x for x in flanks.values()]})
CtoT=F.merge(ClinVar, left_on='RS# (dbSNP)', right_on='RS# (dbSNP)', how='left')
# clinvar may refer to the opposite strand that was used in dbSNP; 
# we want to allow clinvar reference alleles A and T with alternate alleles G and C respectively
# we do not want to allow reference alleles G and C with alternate alleles A and T respectively; these Y-type SNPs must be removed
CtoT=CtoT[(CtoT.ReferenceAlleleVCF=='A') | (CtoT.ReferenceAlleleVCF=='T')].drop_duplicates('RS# (dbSNP)')

#define window limits and the length of the pam including all N residues
windowlen=windowend-windowstart+1
lenpam=len(PAM)
#define a positional preference dictionary

CtoT['gRNAs']=None
CtoT['gRNAall']=None
for i in range(len(CtoT)):
    if type(CtoT.iloc[i].Flanks)==list and CtoT.iloc[i].Flanks!=[]:
        test=CtoT.iloc[i].Flanks[0]
        # define a potential gRNA spacer for each window positioning
        gRNAoptions=[test[(25+windowstart+j-23-lenpam):(25+windowstart+j)] for j in range(windowlen)]
        gRNA=[(gRNAoptions[k],[23+lenpam-x.start() for x in re.finditer('T',gRNAoptions[k]) if windowstart-1<20+lenpam-x.start()<windowend+1]) for k in range(len(gRNAoptions)) if regex.match(enc_pam['r'], gRNAoptions[k][:lenpam])]
        gRNAsingleC=[]
        for g,c in gRNA:
            if g[23+lenpam-windowstart-windowlen+1:23+lenpam-windowstart+1].count('T')==0:
                gRNAsingleC.append(g)
            #elif all([p<priority[x] for x in c]):
                #gRNAsingleC.append(g)
        CtoT.gRNAs.iloc[i]=gRNAsingleC
        CtoT.gRNAall.iloc[i]=[g for g,c in gRNA]

#merge in phenotypes based upon MedGen IDs; remove redundant columns
CtoT=CtoT[['RS# (dbSNP)','GeneSymbol','Name', 'PhenotypeIDS', 'Origin', 'ReviewStatus', 'NumberSubmitters', 'LastEvaluated', 'gRNAs', 'gRNAall']]
ids=[re.findall('MedGen:C.{7}', x) for x in CtoT.PhenotypeIDS.values]
CtoT['Phenotypes']=[[PhenotypeDict[y.lstrip('MedGen:')] for y in x if y.lstrip('MedGen:') in PhenotypeDict.keys()] for x in ids]
CtoT.drop('PhenotypeIDS', inplace=True, axis=1)

#open flanking sequence fasta files for all R-type pathogenic human SNPs (includes both G>A and A>G ref>variant)
#downloaded as fasta file from flanks.py
handle = open("Yfasta.txt", 'r') 
flanks={}
#save as a dictionary keyed on rsID as an Integer with values being 25nt of flanking sequence on each side of the SNP
for record in SeqIO.parse(handle, "fasta") :
    flanks[int(record.id)]=regex.findall('.{25}[^A,T,C,G].{25}', str(record.seq))
    #flanks[int(record.id)]= str(record.seq)
handle.close()
#merge flanking sequences to the AtoG frame on rsID
F=pd.DataFrame({'RS# (dbSNP)': list(flanks.keys()), 'Flanks': [x for x in flanks.values()]})
AtoG=F.merge(ClinVar, left_on='RS# (dbSNP)', right_on='RS# (dbSNP)', how='left')
# clinvar may refer to the opposite strand that was used in dbSNP; 
# we want to allow clinvar reference alleles G and C with alternate alleles A and T respectively
# we do not want to allow reference alleles A and T with alternate alleles G and C respectively; these R-type SNPs must be removed
AtoG=AtoG[(AtoG.ReferenceAlleleVCF=='A') | (AtoG.ReferenceAlleleVCF=='T')].drop_duplicates('RS# (dbSNP)')


AtoG['gRNAs']=None
AtoG['gRNAall']=None
for i in range(len(AtoG)):
    if type(AtoG.iloc[i].Flanks)==list and AtoG.iloc[i].Flanks!=[]:
        test=AtoG.iloc[i].Flanks[0]
        gRNAoptions=[test[(26-windowstart-j):(26-windowstart-j+lenpam+23)] for j in range(windowlen)]
        #if there is an appropriate PAM placed for a given gRNA spacer
        #save tuple of gRNA spacer, and the position of off-target Gs in the window
        gRNA=[(gRNAoptions[k],[x.start()+1 for x in re.finditer('G',gRNAoptions[k]) if windowstart-1<x.start()+1<windowend+1]) for k in range(len(gRNAoptions)) if regex.match(enc_pam['f'], gRNAoptions[k][-lenpam:])]
        gRNAsingleG=[]
        for g,c in gRNA:
            #if the target G is the only G in the window save this as a single G site
            if g[windowstart-1:windowend].count('G')==0:
                gRNAsingleG.append(g)
            #OPTIONAL uncomment the ELIF statement if you are interest in filtered based upon position of off-target G
            #if the target G is expected to be edited more efficiently than the off-target Gs, also save as a single G Site
            #elif all([p<priority[x] for x in c]):
                #gRNAsingleG.append(g)        
        AtoG.gRNAs.iloc[i]=gRNAsingleG
        AtoG.gRNAall.iloc[i]=[g for g,c in gRNA]

AtoG=AtoG[['RS# (dbSNP)','GeneSymbol','Name', 'PhenotypeIDS', 'Origin', 'ReviewStatus', 'NumberSubmitters', 'LastEvaluated', 'gRNAs', 'gRNAall']]
ids=[re.findall('MedGen:C.{7}', x) for x in AtoG.PhenotypeIDS.values]
AtoG['Phenotypes']=[[PhenotypeDict[y.lstrip('MedGen:')] for y in x if y.lstrip('MedGen:') in PhenotypeDict.keys()] for x in ids]
AtoG.drop('PhenotypeIDS', inplace=True, axis=1)

CtoT.to_csv('pathogenic_CtoT_all.csv')
AtoG.to_csv('pathogenic_AtoG_all.csv')

pathogenic_CtoT_hasPAM=CtoT[[type(x)==list and x!=[] for x in CtoT.gRNAall]]
pathogenic_AtoG_hasPAM=AtoG[[type(x)==list and x!=[] for x in AtoG.gRNAall]]

pathogenic_AtoG_hasPAM.to_csv('pathogenic_AtoG_hasPAM.csv')
pathogenic_CtoT_hasPAM.to_csv('pathogenic_CtoT_hasPAM.csv')

pathogenic_CtoT_SingleC=CtoT[[type(x)==list and x!=[] for x in CtoT.gRNAs]]
pathogenic_AtoG_SingleG=AtoG[[type(x)==list and x!=[] for x in AtoG.gRNAs]]

pathogenic_AtoG_SingleG.to_csv('pathogenic_AtoG_SingleG.csv')
pathogenic_CtoT_SingleC.to_csv('pathogenic_CtoT_SingleC.csv')

with open("Summary.txt", "w") as text_file:
    text_file.write("singleC %s \n" % (len(pathogenic_CtoT_SingleC)+len(pathogenic_AtoG_SingleG)))
    text_file.write("hasPAM %s \n" % (len(pathogenic_CtoT_hasPAM)+len(pathogenic_AtoG_hasPAM)))
    text_file.write("Pathogenic SNPs that can be targeted with CBE %s" % (len(CtoT)+len(AtoG)))
    
    
