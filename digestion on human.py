#!/usr/bin/env python
# coding: utf-8

# In[4]:


from pyopenms import *
dig = ProteaseDigestion()
dig.getEnzymeName()
bsa = "".join([l.strip() for l in open("human.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
result = []
dig.digest(bsa, result)
peptides = [AASequence.fromString(s.toString()) for s in result]
for s in result:
    print(s.toString())


# In[ ]:




