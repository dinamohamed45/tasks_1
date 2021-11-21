#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pyopenms
help(pyopenms.Constants)
print ("Avogadro's number is ", pyopenms.Constants.AVOGADRO )


# In[4]:


from pyopenms import *
edb = ElementDB()

edb.hasElement("O")
edb.hasElement("S")

oxygen = edb.getElement("O")
print(oxygen.getName())
print(oxygen.getSymbol())
print(oxygen.getMonoWeight())
print(oxygen.getAverageWeight())

sulfur = edb.getElement("S")
print(sulfur.getName())
print(sulfur.getSymbol())
print(sulfur.getMonoWeight())
print(sulfur.getAverageWeight())
isotopes = sulfur.getIsotopeDistribution()

print ("one mole of oxygen weighs", 2 * oxygen.getAverageWeight(), "grams")
print ("one mole of 16O2 weighs", 2* oxygen.getMonoWeight(), "grams")


# In[5]:


edb = ElementDB()
oxygen_isoDist = {"mass": [], "abundance": []}
sulfur_isoDist = {"mass": [], "abundance": []}

oxygen = edb.getElement("O")
isotopes = oxygen.getIsotopeDistribution()

for iso in isotopes.getContainer():
    print ("Oxygen isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")
    oxygen_isoDist["mass"].append(iso.getMZ())
    oxygen_isoDist["abundance"].append((iso.getIntensity() * 100))

sulfur = edb.getElement("S")
isotopes = sulfur.getIsotopeDistribution()

for iso in isotopes.getContainer():
    print ("Sulfur isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")
    sulfur_isoDist["mass"].append(iso.getMZ())
    sulfur_isoDist["abundance"].append((iso.getIntensity() * 100))


# In[11]:


import math
from matplotlib import pyplot as plt
def adjustText(x1, y1, x2, y2):
    if y1 > y2:
        plt.annotate('%0.3f' % (y2), xy=(x2, y2), xytext=(x2+0.5,y2+9),
                     textcoords='data',
                     arrowprops=dict(arrowstyle="->", color='r', lw=0.5),
                     horizontalalignment='right', verticalalignment='top')
    else:
        plt.annotate('%0.3f' % (y1), xy=(x1, y1), xytext=(x1+0.5,y1+9),
                     textcoords='data',
                     arrowprops=dict(arrowstyle="->", color='r', lw=0.5),
                     horizontalalignment='right', verticalalignment='top')


def plotDistribution(distribution):
    n = len(distribution["mass"])
    for i in range(0, n):
        plt.vlines(x=distribution["mass"][i], ymin=0, ymax=distribution["abundance"][i])
        if int(distribution["mass"][i - 1]) == int(distribution["mass"][i])                 and i != 0:
            adjustText(distribution["mass"][i - 1], distribution["abundance"][i - 1],
                       distribution["mass"][i], distribution["abundance"][i])
        else:
            plt.text(x=distribution["mass"][i],
                     y=(distribution["abundance"][i] + 2),
                     s='%0.3f' % (distribution["abundance"][i]), va='center',
                     ha='center')
    plt.ylim([0, 110])
    plt.xticks(range(math.ceil(distribution["mass"][0]) - 2,
                     math.ceil(distribution["mass"][-1]) + 2))


plt.figure(figsize=(10,7))

plt.subplot(1,2,1)
plt.title("isotopic distribution of oxygen")
plotDistribution(oxygen_isoDist)
plt.xlabel("atomic mass (u)")
plt.ylabel("relative abundance (%)")

plt.subplot(1,2,2)
plt.title("isotopic distribution of sulfur")
plotDistribution(sulfur_isoDist)
plt.xlabel("atomic mass (u)")
plt.ylabel("relative abundance (%)")

plt.show()


# In[12]:


edb = ElementDB()
isotopes = edb.getElement("C").getIsotopeDistribution().getContainer()
carbon_isotope_difference = isotopes[1].getMZ() - isotopes[0].getMZ()
isotopes = edb.getElement("N").getIsotopeDistribution().getContainer()
nitrogen_isotope_difference = isotopes[1].getMZ() - isotopes[0].getMZ()

print ("mass difference between 12C and 13C:", carbon_isotope_difference)
print ("mass difference between 14N and N15:", nitrogen_isotope_difference)
print ("relative deviation:", 100 *(carbon_isotope_difference - nitrogen_isotope_difference)/carbon_isotope_difference, "%")


# In[18]:


from pyopenms.Constants import *

he = ElementDB().getElement("He")
isotopes = he.getIsotopeDistribution()

mass_sum = 2*PROTON_MASS_U + 2*ELECTRON_MASS_U + 2*NEUTRON_MASS_U
he4 = isotopes.getContainer()[1].getMZ()
print ("sum of masses of 2 protons, neutrons and electrons:", mass_sum)

print ("mass of He4:", he4)

print ("difference between the two masses:", 100*(mass_sum - he4)/mass_sum, "%")


# In[19]:


methanol = EmpiricalFormula("CH3OH")
water = EmpiricalFormula("H2O")
ethanol = EmpiricalFormula("CH2") + methanol
print("ethanol chemical formula:", ethanol.toString())
print("ethanol composition :", ethanol.getElementalComposition())
print("ethanol has", ethanol.getElementalComposition()[b"H"], "hydrogen atoms")


# In[22]:


methanol = EmpiricalFormula("CH3OH")
ethanol = EmpiricalFormula("CH2") + methanol

methanol_isoDist = {"mass": [], "abundance": []}
ethanol_isoDist = {"mass": [], "abundance": []}

print("coarse isotope distribution:")
isotopes = ethanol.getIsotopeDistribution( CoarseIsotopePatternGenerator(4) )
prob_sum = sum([iso.getIntensity() for iso in isotopes.getContainer()])
print("this covers", prob_sum, "probability")
for iso in isotopes.getContainer():
    print ("isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, " %")
    methanol_isoDist["mass"].append(iso.getMZ())
    methanol_isoDist["abundance"].append((iso.getIntensity() * 100))

print("fine isotope distribution:")
isotopes = ethanol.getIsotopeDistribution( FineIsotopePatternGenerator(1e-3) )
prob_sum = sum([iso.getIntensity() for iso in isotopes.getContainer()])
print("this covers", prob_sum, "probability")
for iso in isotopes.getContainer():
    print ("isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "% ")
    ethanol_isoDist["mass"].append(iso.getMZ())
    ethanol_isoDist["abundance"].append((iso.getIntensity() * 100))


# In[23]:


plt.figure(figsize=(10,7))

plt.subplot(1,2,1)
plt.title("isotopic distribution of methanol")
plotDistribution(methanol_isoDist)
plt.xlabel("Atomic mass (u)")
plt.ylabel("Relative abundance (%)")

plt.subplot(1,2,2)
plt.title("isotopic distribution of ethanol")
plotDistribution(ethanol_isoDist)
plt.xlabel("atomic mass (u)")
plt.ylabel("relative abundance (%)")

plt.savefig("methanol_ethanol_isoDistribution.png")


# In[24]:


methanol = EmpiricalFormula("CH3OH")
ethanol = EmpiricalFormula("CH2") + methanol

print("fine isotope distribution:")
isotopes = ethanol.getIsotopeDistribution( FineIsotopePatternGenerator(1e-6) )
prob_sum = sum([iso.getIntensity() for iso in isotopes.getContainer()])
print("this covers", prob_sum, "probability")
for iso in isotopes.getContainer():
    print ("isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "% ")


# In[25]:


isotopes = ethanol.getIsotopeDistribution( CoarseIsotopePatternGenerator(5, True) )
for iso in isotopes.getContainer():
    print ("Isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")


# In[26]:


lys = ResidueDB().getResidue("Lysine")
print(lys.getName())
print(lys.getThreeLetterCode())
print(lys.getOneLetterCode())
print(lys.getAverageWeight())
print(lys.getMonoWeight())
print(lys.getPka())
print(lys.getFormula().toString())


# In[28]:


ox = ModificationsDB().getModification("Oxidation")
print(ox.getUniModAccession())
print(ox.getUniModRecordId())
print(ox.getDiffMonoMass())
print(ox.getId())
print(ox.getFullId())
print(ox.getFullName())
print(ox.getDiffFormula())


# In[29]:


isotopes = ox.getDiffFormula().getIsotopeDistribution(CoarseIsotopePatternGenerator(5))
for iso in isotopes.getContainer():
    print (iso.getMZ(), ":", iso.getIntensity())


# In[32]:


uridine = RibonucleotideDB().getRibonucleotide(b"U")
print(uridine.getName())
print(uridine.getCode())
print(uridine.getAvgMass())
print(uridine.getMonoMass())
print(uridine.getFormula().toString())
print(uridine.isModified())
methyladenosine = RibonucleotideDB().getRibonucleotide(b"m1A")
print(methyladenosine.getName())
print(methyladenosine.isModified())


# In[ ]:





# In[ ]:




