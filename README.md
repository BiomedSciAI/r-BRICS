# r-BRICS – a revised BRICS module that breaks ring structures and carbon chains
![toc](https://github.com/BiomedSciAI/r-BRICS/blob/main/toc.jpg)

Figure 1. An example of using r-BRICS using lanostane.

## Reference
L. Zhang, V. Rao, W. Cornell, r-BRICS – a revised BRICS module that breaks ring structures and carbon chains (submitted).

## Pre-requisites
```
rdkit >= 2022.09.5
```

## Python environment
We recommend to use miniconda with python 3.6+.

## Example
First, rdkit and r-BRICS need to be imported. Simply put ```reBRICS_public.py``` in your workfolder to use it (or python library folder).
```
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from reBRICS_public import *
import rdkit
```
(1) POPC, as shown in publication
```
popc=Chem.MolFromSmiles("CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\CCCCCCCC")
Draw.MolToImage(popc,size=(800,300))
bonds=BRICS.FindrBRICSBonds(popc)
pieces=BRICS.BreakrBRICSBonds(popc,bonds)
frags=Chem.GetMolFrags(pieces,asMols=True)
Draw.MolsToGridImage(frags,subImgSize=(700,200))
```
(2) Lanostane, as shown in publication
```
lano=Chem.MolFromSmiles('C[C@H](CCCC(C)C)[C@@]1([H])CC[C@@]2(C)[C@]3([H])CC[C@@]4([H])C(C)(C)CCC[C@]4(C)[C@@]3([H])CC[C@@]21C')
bonds=FindrBRICSBonds(lano)
pieces=BreakrBRICSBonds(lano,bonds)
frags=Chem.GetMolFrags(pieces,asMols=True)
Draw.MolsToGridImage(frags,subImgSize=(700,200))
```
(3) n-hexapentacontane with ```reBRICS```, as shown in publication
```
ccc=Chem.MolFromSmiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
### Regular r-BRICS, without reBRICS
bonds=FindrBRICSBonds(ccc)
pieces=BreakrBRICSBonds(ccc,testbonds)
frags=Chem.GetMolFrags(pieces,asMols=True)
Draw.MolsToGridImage(frags,subImgSize=(700,300))
### With reBRICS
newfrags=reBRICS(frags) ### Note that input for reBRICS is a list of fragments from regular r-BRICS
Draw.MolsToGridImage(newfrags,subImgSize=(700,200))
```
(4) Breaking aromatic fused rings
* Regular r-BRICS does not break aromatic bonds. This might be preferrable in some situations.
```
p135530454=Chem.MolFromSmiles("NC1=CC=C2N(CCC3CC3)C(=O)C(=C(O)C2=C1)C1=NS(=O)(=O)C2=C(N1)C=CC=C2")
bonds=FindrBRICSBonds(p135530454)
pieces=BreakrBRICSBonds(p135530454,pbonds)
frags=Chem.GetMolFrags(pieces,asMols=True)
Draw.MolsToGridImage(frags,subImgSize=(700,200))
```
* In order to break aromatic bonds, all bonds must be kekulized.
```
p135530454=Chem.MolFromSmiles("NC1=CC=C2N(CCC3CC3)C(=O)C(=C(O)C2=C1)C1=NS(=O)(=O)C2=C(N1)C=CC=C2")
Chem.Kekulize(p135530454,clearAromaticFlags=True)
bonds=FindrBRICSBonds(p135530454)
pieces=BreakrBRICSBonds(p135530454,pbonds)
frags=Chem.GetMolFrags(pieces,asMols=True)
Draw.MolsToGridImage(frags,subImgSize=(700,200))
```
(5) Other examples shown in the publication
```
ritonavir=Chem.MolFromSmiles("CC(C)C1=NC(=CS1)CN(C)C(=O)NC(C(C)C)C(=O)NC(CC2=CC=CC=C2)CC(C(CC3=CC=CC=C3)NC(=O)OCC4=CN=CS4)O")
```
