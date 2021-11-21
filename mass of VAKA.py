#!/usr/bin/env python
# coding: utf-8

# In[33]:


s=0
seq = AASequence.fromString("VAKA")

for aa in seq:
    s+=aa.getMonoWeight()
s
   


# In[34]:


seq = AASequence.fromString("VAKA")
t = seq.getMonoWeight()
t


# In[36]:


t==s


# In[ ]:





# In[ ]:





# In[ ]:




