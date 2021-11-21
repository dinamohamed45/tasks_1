#!/usr/bin/env python
# coding: utf-8

# In[3]:


from pyopenms import *
from urllib.request import urlretrieve

dig = ProteaseDigestion()
dig.getEnzymeName() 
bsa = "".join([l.strip() for l in open("yeast.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
result = []
dig.digest(bsa, result)
for s in result:
    print(s.toString())


# In[ ]:




