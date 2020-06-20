#!/anaconda3/bin/python

"""
Created on May 26,2020
Kalina Grabb

This script it to process the data of coral sequences blasted against NOX
Data from blastp_v2_all_py.csv
Blast data from coral online and NOX from Kawahara et al., 2007
"""


"""
Section 1 - create heatmaps for the data to display the blast results
"""
# %%
import pandas as pd 
import matplotlib 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np 
import seaborn as sb
import math



# %% import data   
# Chose which file to import by changing the filename, coralSeqUniq, and qSeqUniq
filename = "/blastp_v2_OR_py.csv"
coralSeqUniq = "/coralSeqUniqOR"
qSeqUniq = "/qSeqUniqOR"

#File import
import os
path = os.getcwd()
print(path)
print(path+filename)
data = pd.read_csv(path+filename, header=0)
print(data.head(5))

#Names import for coral and nox
coralU = pd.read_csv(path+coralSeqUniq, names=["coral"])
print(coralU.head(5))
qU =pd.read_csv(path+qSeqUniq, names=["nox"])
print(qU.head(5))
namesU = pd.concat([coralU, qU], axis=1)
print(namesU.head(5))

# %% name variables to be graphed in heatmap
bits = data.bitScore.values
nox = data.qSequence.values
coral = data.coralSequence.values
perID = data.pID.values


# %% Create an array for the bit score with coral and nox on each axis
"""
Use I and J to do a loop and label the unique values into a new array. Try and graph this
if i = nox
and j = coral
insert bit score into cell i, j
"""

corU = namesU.coral
noxU = namesU.nox
lenN = len(noxU)
lenC = len(corU)
bitArray = np.empty([lenC, lenN])
pidArray = np.empty([lenC, lenN])
evalArray = np.empty([lenC, lenN])
for j in range(lenN):
    for i in range(lenC):
            corSame = np.where(data.coralSequence==corU[i])[0]
            noxSame = np.where(data.qSequence==noxU[j])[0]
            matchSame = (set(corSame) & set(noxSame))
            if len(matchSame)!=0:
                rowSame = list(matchSame)[0]
                bitSame = data.bitScore[rowSame]
                pidSame = data.pID[rowSame]
                evalSame = data.evalue[rowSame]
            else:
                bitSame = float("nan")
                pidSame = float("nan")
                evalSame = float("nan")
            bitArray[i][j] = bitSame
            pidArray[i][j] = pidSame
            evalArray[i][j] = evalSame
# two = np.arange(4*5).reshape(4,5)
#any(namesU.coral=='KXJ16767.1')
#namesU.coral.isin(['KXJ16767.1'])
#np.where(res==True)
#np.where(corU=='KXJ16767.1')[0]


# %% Create heatmap for bitscore
#Heat map of the bit score
ax=sb.heatmap(bitArray, xticklabels=qU.nox, yticklabels=coralU.coral, linewidths=.05, 
    cmap="PuBuGn", vmin = 30, vmax = 650)


ax.set_title(f"Bit score \n {filename}", fontsize=15)
ax.set_ylabel("Coral Sequences", fontsize=12)
ax.set_xlabel("NOX Sequences", fontsize=12)
ax.tick_params(labelsize=2)
plt.show()




#plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
#

# %% Heat map of pID
ax=sb.heatmap(pidArray, xticklabels=qU.nox, yticklabels=coralU.coral, linewidths=.05, 
    cmap="PuBuGn", vmin=50, vmax=58)


ax.set_title(f"Percent ID \n {filename}", fontsize=15)
ax.set_ylabel("Coral Sequences", fontsize=12)
ax.set_xlabel("NOX Sequences", fontsize=12)
ax.tick_params(labelsize=2)
plt.show()


# %% Heat map of evalue
logArray = np.log10(evalArray)
minArray = np.nanmin(evalArray)
maxArray = np.nanmax(evalArray)
ax=sb.heatmap(logArray, xticklabels=qU.nox, yticklabels=coralU.coral, linewidths=.05, 
    cmap="PuBuGn", vmin = -180, vmax = -70)
#, norm=colors.LogNorm(vmin=0, vmax=1e-97)

ax.set_title(f"log(e value) \n {filename}", fontsize=15)
ax.set_ylabel("Coral Sequences", fontsize=12)
ax.set_xlabel("NOX Sequences", fontsize=12)
ax.tick_params(labelsize=2)
plt.show()




"""
Section 2
Phylogenetic tree using coral sequences that show most simlartiy with NOX
Using Biopython, pieced together wtih muscle on Poseidon and different commands
This one uses a fasta file to build a tree
"""
# %%
import Bio 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure



from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
#from Bio.Phylo.Applications import PhymlCommandline
#from Bio.Phylo.PAML import codeml
#from Bio.Phylo.PhyloXML import Phylogeny

#import networkx
#from networkx.drawing import nx_agraph
#networkx.graphviz_layout = nx_agraph.graphviz_layout

# %% File import
msaFile = "/muscleAlignmentCoralSeq.msa"

import os
path = os.getcwd()
#path = "/Users/kgrabb/Documents/2018.05CoralLarvae/Genomes/Poseidon/blastResults/v2"
print(path)
print(path+msaFile)
inputFile = path+msaFile
viewFile = pd.read_csv(inputFile)
print(viewFile.head(5))


alignment = AlignIO.read(inputFile, "fasta")
print(alignment)

#Calculate distance matrix
calculator = DistanceCalculator('identity')
dm=calculator.get_distance(alignment)

#print("Distance Matrix")
#print(dm)

#Construct phylogenetic tree with UPGMA algorithm
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)



# %% plot the tree
#print(tree)

figFileName = "treeFig.png"
figFile = path+"/"+figFileName

SMALL_SIZE = 3
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# plt.rcParams.update({'font.size':2})
# plt.rcParams.update({'lines.linewidth':1})
tree.ladderize()
Phylo.draw(tree, branch_labels=None)
plt.savefig(figFile)
#plt.figure(figsize=(12, 12), dpi=100)

#Phylo.draw_ascii(tree)



# %% Plot tree using a def function
figFileName = "treeFig2.png"
figFile = path+"/graphs"+figFileName



def plotTree(treedata, treeFigure):
    SMALL_SIZE = 7
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    fig = plt.figure(figsize=(15,20), dpi=100)
    axes=fig.add_subplot(1,1,1)

    treedata.root.color = "gray"
    #treedata.clade({"name": "adi2mcaRNA8168_R0"}).color = "blue"
    Phylo.draw(treedata, axes=axes, branch_labels=None)
    plt.savefig(treeFigure, dpi=100)
    return

plotTree(tree, figFile)

# %% Plot tree using a def function and importing xml file from Poseidon processing
from Bio import Phylo
#from Bio.Phylo.Consesus import *

figFileName = "treeFig2.png"
figFile = path+"/graphs"+figFileName
treeFileName = "treeFileCoral.xml"
treeFile = path+"/"+treeFileName
treeFilePhylo = Phylo.read(treeFile, 'phyloxml')

nox2File = path+"/treeCoralSeqOnlyNox2"
nox2Names = pd.read_csv(nox2File)

def plotTree(treedata, treeFigure, nox2):
    SMALL_SIZE = 8
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 40
    titleName = "All NOX-like"
    colorRoot = "gray"

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    fig = plt.figure(figsize=(15,20), dpi=100)
    axes=fig.add_subplot(1,1,1)

    treedata.root.color = colorRoot
    #treedata.common_ancestor({"name":"KXJ24335.1"}).color = "blue"
    #treedata.common_ancestor("KXJ24335.1").color = "blue"
    treedata.common_ancestor(nox2).color = "blue"
    
    #for clade in tree.find_clades(terminal=True):
     #   clade.color=Phylo.PhyloXML.BranhColor(0,0,250)
     #   Phylo.write(treedata, "colored.xml", "phyloxml")

    plt.title(titleName, fontsize=BIGGER_SIZE)
    #treedata.clade({"name": "adi2mcaRNA8168_R0"}).color = "blue"
    Phylo.draw(treedata, axes=axes, branch_labels=None)
    plt.savefig(treeFigure, dpi=100)
    return 


returned = plotTree(treeFilePhylo, figFile, nox2Names)


# %%
"""
Section 3
Following biopython pipeline in order to perform commandline muscle for msa
uses AlignIO to convert the fasta file and create an alignment and PhyML. 
then Phylo tree to build/draw a tree
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
inputFilename = "bioCoralSeqRecord"

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

# %% convert alignment file from fasta from muscle in Poseidon to phyli-relaxed
msaFile = "/muscleAlignmentCoralSeq.msa"
outFile = "/muscleAlignmentCoralSeq.phy"

import os
path = os.getcwd()
#path = "/Users/kgrabb/Documents/2018.05CoralLarvae/Genomes/Poseidon/blastResults/v2"
print(path)
print(path+msaFile)
inputFile = path+msaFile
outputFile = path+outFile
viewFile = pd.read_csv(inputFile)
print(viewFile.head(5))

AlignIO.convert(inputFile, "fasta", outputFile, "phylip-relaxed")



# %% use PhyML. feed in phy alignment with the command line wrapper
from Bio.Phylo.Applications import PhymlCommandline
cmdline = PhymlCommandline(input=outputFile, datatype="aa", model="WAG", alpha="e", bootstrap=100)
out_log, err_log = cmdline()


# %%

