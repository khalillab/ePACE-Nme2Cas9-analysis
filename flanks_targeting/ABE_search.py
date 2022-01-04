import numpy as np
import pandas as pd
import regex
import re
from Bio import SeqIO
import Bio

###CHANGE PAM AND WINDOW INFO HERE###
PAM = 'NNNNNCN'
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


ClinVar=pd.read_csv('Clinvar.csv',sep='\t', encoding = "ISO-8859-1", low_memory=False)
Phenotypes=pd.read_csv('DiseaseNames.csv', sep='\t')
PhenotypeDict=dict(zip(Phenotypes.CUI, Phenotypes.name))

#open flanking sequence fasta files for all Y-type pathogenic human SNPs (includes both C>T and T>C ref>variant)
#downloaded as fasta file from flanks.py script
handle = open("Yfasta.txt", 'r')
flanks={}
#save as a dictionary keyed on rsID as an Integer with values being 25nt of flanking sequence on each side of the SNP
for record in SeqIO.parse(handle, "fasta") :
    flanks[int(record.id)]=regex.findall('.{25}[^A,T,C,G].{25}', str(record.seq))
handle.close()
#merge flanking sequences to the TtoC frame on rsID
F=pd.DataFrame({'RS# (dbSNP)': list(flanks.keys()), 'Flanks': [x for x in flanks.values()]})
TtoC=F.merge(ClinVar, left_on='RS# (dbSNP)', right_on='RS# (dbSNP)', how='left')
# clinvar may refer to the opposite strand that was used in dbSNP; 
# we want to allow clinvar reference alleles G and C with alternate alleles A and T respectively
# we do not want to allow reference alleles A and T with alternate alleles G and C respectively; these Y-type SNPs must be removed
TtoC=TtoC[(TtoC.ReferenceAlleleVCF=='G') | (TtoC.ReferenceAlleleVCF=='C')].drop_duplicates('RS# (dbSNP)')

#define window limits and the length of the pam including all N residues
windowlen=windowend-windowstart+1
lenpam=len(PAM)
#define a positional preference dictionary

TtoC['gRNAs']=None
TtoC['gRNAall']=None
for i in range(len(TtoC)):
    if type(TtoC.iloc[i].Flanks)==list and TtoC.iloc[i].Flanks!=[]:
        test=TtoC.iloc[i].Flanks[0]
        # define a potential gRNA spacer for each window positioning
        gRNAoptions=[test[(25+windowstart+j-23-lenpam):(25+windowstart+j)] for j in range(windowlen)]
        gRNA=[(gRNAoptions[k],[23+lenpam-x.start() for x in re.finditer('C',gRNAoptions[k]) if windowstart-1<23+lenpam-x.start()<windowend+1]) for k in range(len(gRNAoptions)) if regex.match(enc_pam['r'], gRNAoptions[k][:lenpam])]
        gRNAsingleA=[]
        for g,c in gRNA:
            if g[23+lenpam-windowstart-windowlen+1:23+lenpam-windowstart+1].count('C')==0:
                gRNAsingleA.append(g)
            #elif all([p<priority[x] for x in c]):
                #gRNAsingleA.append(g)
        TtoC.gRNAs.iloc[i]=gRNAsingleA
        TtoC.gRNAall.iloc[i]=[g for g,c in gRNA]

#merge in phenotypes based upon MedGen IDs; remove redundant columns
TtoC=TtoC[['RS# (dbSNP)','GeneSymbol','Name', 'PhenotypeIDS', 'Origin', 'ReviewStatus', 'NumberSubmitters', 'LastEvaluated', 'gRNAs', 'gRNAall']]
ids=[re.findall('MedGen:C.{7}', x) for x in TtoC.PhenotypeIDS.values]
TtoC['Phenotypes']=[[PhenotypeDict[y.lstrip('MedGen:')] for y in x if y.lstrip('MedGen:') in PhenotypeDict.keys()] for x in ids]
TtoC.drop('PhenotypeIDS', inplace=True, axis=1)

#open flanking sequence fasta files for all R-type pathogenic human SNPs (includes both G>A and A>G ref>variant)
#downloaded as fasta file from flanks.py
handle = open("Rfasta.txt", 'r') 
flanks={}
#save as a dictionary keyed on rsID as an Integer with values being 25nt of flanking sequence on each side of the SNP
for record in SeqIO.parse(handle, "fasta") :
    flanks[int(record.id)]=regex.findall('.{25}[^A,T,C,G].{25}', str(record.seq))
handle.close()
#merge flanking sequences to the AtoG frame on rsID
F=pd.DataFrame({'RS# (dbSNP)': list(flanks.keys()), 'Flanks': [x for x in flanks.values()]})
GtoA=F.merge(ClinVar, left_on='RS# (dbSNP)', right_on='RS# (dbSNP)', how='left')
# clinvar may refer to the opposite strand that was used in dbSNP; 
# we want to allow clinvar reference alleles A and T with alternate alleles G and C respectively
# we do not want to allow reference alleles G and C with alternate alleles A and T respectively; these R-type SNPs must be removed
GtoA=GtoA[(GtoA.ReferenceAlleleVCF=='G') | (GtoA.ReferenceAlleleVCF=='C')].drop_duplicates('RS# (dbSNP)')


GtoA['gRNAs']=None
GtoA['gRNAall']=None
for i in range(len(GtoA)):
    if type(GtoA.iloc[i].Flanks)==list and GtoA.iloc[i].Flanks!=[]:
        test=GtoA.iloc[i].Flanks[0]
        gRNAoptions=[test[(26-windowstart-j):(26-windowstart-j+lenpam+23)] for j in range(windowlen)]
        #if there is an appropriate PAM placed for a given gRNA spacer
        #save tuple of gRNA spacer, and the position of off-target As in the window
        gRNA=[(gRNAoptions[k],[x.start()+1 for x in re.finditer('G',gRNAoptions[k]) if windowstart-1<x.start()+1<windowend+1]) for k in range(len(gRNAoptions)) if regex.match(enc_pam['f'], gRNAoptions[k][-lenpam:])]
        gRNAsingleA=[]
        for g,c in gRNA:
            #if the target A is the only A in the window save this as a single A site
            if g[windowstart-1:windowend].count('G')==0:
                gRNAsingleA.append(g)
            #OPTIONAL uncomment the ELIF statement if you are interest in filtered based upon position of off-target A
            #if the target A is expected to be edited more efficiently than the off-target As, also save as a single A Site
            #elif all([p<priority[x] for x in c]):
                #gRNAsingleA.append(g)        
        GtoA.gRNAs.iloc[i]=gRNAsingleA
        GtoA.gRNAall.iloc[i]=[g for g,c in gRNA]

GtoA=GtoA[['RS# (dbSNP)','GeneSymbol','Name', 'PhenotypeIDS', 'Origin', 'ReviewStatus', 'NumberSubmitters', 'LastEvaluated', 'gRNAs', 'gRNAall']]
ids=[re.findall('MedGen:C.{7}', x) for x in GtoA.PhenotypeIDS.values]
GtoA['Phenotypes']=[[PhenotypeDict[y.lstrip('MedGen:')] for y in x if y.lstrip('MedGen:') in PhenotypeDict.keys()] for x in ids]
GtoA.drop('PhenotypeIDS', inplace=True, axis=1)

TtoC.to_csv('pathogenic_TtoC_all.csv')
GtoA.to_csv('pathogenic_GtoA_all.csv')

pathogenic_TtoC_hasPAM=TtoC[[type(x)==list and x!=[] for x in TtoC.gRNAall]]
pathogenic_GtoA_hasPAM=GtoA[[type(x)==list and x!=[] for x in GtoA.gRNAall]]

pathogenic_GtoA_hasPAM.to_csv('pathogenic_GtoA_hasPAM.csv')
pathogenic_TtoC_hasPAM.to_csv('pathogenic_TtoC_hasPAM.csv')

pathogenic_TtoC_SingleA=TtoC[[type(x)==list and x!=[] for x in TtoC.gRNAs]]
pathogenic_GtoA_SingleA=GtoA[[type(x)==list and x!=[] for x in GtoA.gRNAs]]

pathogenic_GtoA_SingleA.to_csv('pathogenic_GtoA_SingleA.csv')
pathogenic_TtoC_SingleA.to_csv('pathogenic_TtoC_SingleA.csv')

with open("Summary.txt", "w") as text_file:
    text_file.write("singleA %s \n" % (len(pathogenic_TtoC_SingleA)+len(pathogenic_GtoA_SingleA)))
    text_file.write("hasPAM %s \n" % (len(pathogenic_TtoC_hasPAM)+len(pathogenic_GtoA_hasPAM)))
    text_file.write("Pathogenic SNPs that can be targeted with ABE %s" % (len(TtoC)+len(GtoA)))
    
    
