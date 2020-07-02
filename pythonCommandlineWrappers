# Learning to use commandline wrappers in Python

### Kalina Grabb

### June 2020



## Goal

To use commandline wrappers in Biopython in order to access teh `musclecommandline` tools.

Documentation from Biopython is here: https://biopython.org/DIST/docs/api/Bio.Align.Applications._Muscle.MuscleCommandline-class.html

Biopython muscle source code is here: https://biopython.org/DIST/docs/api/Bio.Align.Applications._Muscle-pysrc.html

Biopython cookbook for muscle is here: http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec92



## Script

This script is written in Python program. I have also tried running a Python script in bash.

The source file is: `bioCoralSeqKwNoxD.fas`



```python
#!/anaconda3/bin/python

"""
Section 1
Following biopython pipeline in order to perform commandline muscle for msa
"""
# %% Import functions
import Bio 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
import Bio.Align.Applications
from Bio.Align.Applications import MuscleCommandline



# %% import file names
inputFilename = "bioCoralSeqKwNoxD.fas"

import os
path = os.getcwd()
#path = "/Users/kgrabb/Documents/2018.05CoralLarvae/Genomes/Poseidon/blastResults/v2"
print(path)
inputFile = path+"/"+inputFilename
print(inputFile)
viewFile = pd.read_csv(inputFile)
print(viewFile.head(5))

# %% Run Muscle to obtain a multiple sequence alignment - this DID NOT WORK
#instead, took muscle file that was run on Poseidon and proceeded
from Bio.Align.Applications import MuscleCommandline
musFileName = "bioCoralSeqMuscleOut.aln"
musFile = path+"/"+musFileName
muscle_exe = "//anaconda3/lib/python3.7/site-packages/Bio/Application/_Phyml.py"
cmdline = MuscleCommandline(muscle_exe, input = inputFile, out=musFile, clw=True)
print(cmdline)
cmdline()

#assert os.path.isfile(muscle_exe), "Muscle executable missing"
stdout, stderr = cmdline()

# %% test to run clustalW to try and learn commandline
import os
from Bio.Align.Applications import ClustalwCommandline
clusFileName = "bioCoralSeqMuscleOut.aln"
clusFile = path+"/"+clusFileName

clustalw_exe = r"//anaconda3/lib/python3.7/site-packages/Bio/Application/__init__.py"
cline = ClustalwCommandline(clustalw_exe, infile=clusFile)
print(cline)
#cline()

assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
stdout, stderr = cline()

```



## Problems

I cannot access the environment correctly to execute a commandline wrapper.