from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rBRICS_public import *
import rdkit

#POPC
popc=Chem.MolFromSmiles("CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\CCCCCCCC")
pbonds=FindrBRICSBonds(popc)
ppieces3=BreakrBRICSBonds(popc,pbonds)
pfrags3=Chem.GetMolFrags(ppieces3,asMols=True)
Draw.MolsToGridImage(pfrags3,subImgSize=(700,200))

#Lanostane
lano=Chem.MolFromSmiles('C[C@H](CCCC(C)C)[C@@]1([H])CC[C@@]2(C)[C@]3([H])CC[C@@]4([H])C(C)(C)CCC[C@]4(C)[C@@]3([H])CC[C@@]21C')
lbonds=FindrBRICSBonds(lano)
lpieces3=BreakrBRICSBonds(lano,lbonds)
lfrags3=Chem.GetMolFrags(lpieces3,asMols=True)
Draw.MolsToGridImage(lfrags3,subImgSize=(700,200))

#Ritonavir
rito=Chem.MolFromSmiles("CC(C)C1=NC(=CS1)CN(C)C(=O)NC(C(C)C)C(=O)NC(CC2=CC=CC=C2)CC(C(CC3=CC=CC=C3)NC(=O)OCC4=CN=CS4)O")
pbonds=FindrBRICSBonds(rito)
ppieces3=BreakrBRICSBonds(rito,pbonds)
pfrags3=Chem.GetMolFrags(ppieces3,asMols=True)
Draw.MolsToGridImage(pfrags3,subImgSize=(700,200))

#P135530454, this one keeps aromatic rings
hcv=Chem.MolFromSmiles("NC1=CC=C2N(CCC3CC3)C(=O)C(=C(O)C2=C1)C1=NS(=O)(=O)C2=C(N1)C=CC=C2")
pbonds=FindrBRICSBonds(hcv)
ppieces3=BreakrBRICSBonds(hcv,pbonds)
pfrags3=Chem.GetMolFrags(ppieces3,asMols=True)
Draw.MolsToGridImage(pfrags3,subImgSize=(700,300))

#P135530454, this one breaks aromatic rings
hcv=Chem.MolFromSmiles("NC1=CC=C2N(CCC3CC3)C(=O)C(=C(O)C2=C1)C1=NS(=O)(=O)C2=C(N1)C=CC=C2")
Chem.Kekulize(hcv,clearAromaticFlags=True)
#AllChem.Compute2DCoords(hcv)
pbonds=FindrBRICSBonds(hcv)
ppieces3=BreakrBRICSBonds(hcv,pbonds)
pfrags3=Chem.GetMolFrags(ppieces3,asMols=True)
Draw.MolsToGridImage(pfrags3,subImgSize=(700,300))

#reBRICS, which breaks carbon chain polymers
ccc=Chem.MolFromSmiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
testbonds=FindrBRICSBonds(ccc)
testpieces=BreakrBRICSBonds(ccc,testbonds)
testfrags=Chem.GetMolFrags(testpieces,asMols=True)
newfrags=reBRICS(testfrags)
Draw.MolsToGridImage(newfrags,subImgSize=(700,200))
Draw.MolsToGridImage(testfrags,subImgSize=(700,300))
