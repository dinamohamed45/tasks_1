#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().system('pip install biopython')


# In[12]:


get_ipython().system('pip install pyopenms')


# In[30]:


from pyopenms import *
seq = AASequence.fromString("DFPIANGER")
prefix = seq.getPrefix(4) 
suffix = seq.getSuffix(5) 
concat = seq + seq 
print("sequence:", seq)
print("prefix:", prefix)
print("suffix:", suffix)
print("concatenated:", concat)
mfull = seq.getMonoWeight() 
mprecursor = seq.getMonoWeight(Residue.ResidueType.Full, 2)

mz = seq.getMonoWeight(Residue.ResidueType.Full, 2) / 2.0 
mz = seq.getMZ(2) 
print()
print("monoisotopic mass of peptide [M] is", mfull)
print("monoisotopic mass of peptide precursor [M+2H]2+ is", mprecursor)
print("monoisotopic m/z of [M+2H]2+ is", mz)


# In[16]:


seq = AASequence.fromString("DFPIANGER")

print("The peptide", str(seq), "consists of following amino acid:")

for aa in seq:

    print(aa.getName(), ":", aa.getMonoWeight())


# In[17]:


seq = AASequence.fromString("C[143]PKCK(Label:13C(6)15N(2))CR")

if seq.hasNTerminalModification():
    print("N-Term Modification: ", seq.getNTerminalModification().getFullId())
    
if seq.hasCTerminalModification():
    print("C-Term Modification: ", seq.getCTerminalModification().getFullId())
    
for aa in seq:
    if (aa.isModified()):
        
        print(aa.getName(), ":", aa.getMonoWeight(), ":", aa.getModificationName())
    else:
        print(aa.getName(), ":", aa.getMonoWeight())


# # Molecular formula
# 

# In[18]:


seq = AASequence.fromString("DFPIANGER")
seq_formula = seq.getFormula()
print("Peptide", seq, "has molecular formula", seq_formula)


# # Isotope patterns
# 

# In[20]:


coarse_isotopes = seq_formula.getIsotopeDistribution( CoarseIsotopePatternGenerator(6) )

for iso in coarse_isotopes.getContainer():
    print ("Isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")


# In[21]:


fine_isotopes = seq_formula.getIsotopeDistribution( FineIsotopePatternGenerator(0.01) ) 
for iso in fine_isotopes.getContainer():
    print ("Isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")


# In[32]:



suffix = seq.getSuffix(3)
print("="*35)
print("y3 ion sequence:", suffix)
y3_formula = suffix.getFormula(Residue.ResidueType.YIon, 2) 

suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0 
suffix.getMonoWeight(Residue.ResidueType.XIon, 2) / 2.0 
suffix.getMonoWeight(Residue.ResidueType.BIon, 2) / 2.0 

print("y3 mz:", suffix.getMonoWeight(Residue.ResidueType.YIon, 2) / 2.0 )
print("y3 molecular formula:", y3_formula)


# In[42]:


bsa = FASTAEntry() 
bsa.sequence = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGE"

bsa.description = "BSA Bovine Albumin (partial sequence)"

bsa.identifier = "BSA"

alb = FASTAEntry()
alb.sequence = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGE"
alb.description = "ALB Human Albumin (partial sequence)"
alb.identifier = "ALB"

entries = [bsa, alb]

f = FASTAFile()
f.store("example.fasta", entries)
    
entries = []
f = FASTAFile()
f.load("example.fasta", entries)
print( len(entries) )
for e in entries:
 print (e.identifier, e.sequence)
    


# In[39]:





# In[ ]:





# In[ ]:




