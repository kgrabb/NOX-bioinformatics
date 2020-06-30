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
# %% Import all programs
import Bio 
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib.lines import Line2D as ln
import pandas as pd 

from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import os
#from Bio.Phylo.Applications import PhymlCommandline
#from Bio.Phylo.PAML import codeml
#from Bio.Phylo.PhyloXML import Phylogeny

#import networkx
#from networkx.drawing import nx_agraph
#networkx.graphviz_layout = nx_agraph.graphviz_layout

# %% File import
msaFile = "/muscleAlignmentCoralSeq.msa"


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
from Bio.Phylo import PhyloXML as PX
#from Bio.Phylo.Consesus import *

path = os.getcwd()
figFileName = "treeFig2.png"
figFile = path+"/graphs"+figFileName
treeFileName = "treeFileCoralUPGMAD.xml" #change this filename
treeFile = path+"/"+treeFileName
treeFilePhylo = Phylo.read(treeFile, 'phyloxml')

nox1File = path+"/treeCoralSeqOnlyNox1"
nox1Names = pd.read_csv(nox1File, squeeze=True)
nox2File = path+"/treeCoralSeqOnlyNox2"
nox2Names = pd.read_csv(nox2File, squeeze=True)
nox3File = path+"/treeCoralSeqOnlyNox3"
nox3Names = pd.read_csv(nox3File, squeeze=True)
nox4File = path+"/treeCoralSeqOnlyNox4"
nox4Names = pd.read_csv(nox4File, squeeze=True)
nox5File = path+"/treeCoralSeqOnlyNox5"
nox5Names = pd.read_csv(nox5File, squeeze=True)
noxDFile = path+"/treeCoralSeqOnlyNoxD"
noxDNames = pd.read_csv(noxDFile, squeeze=True)

def plotTree(treedata, treeFigure, nox1, nox2, nox3, nox4, nox5, noxD):
    SMALL_SIZE = 8
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 40
    titleName = f"NOX2-like Protein Sequences in Coral \n Bootstrap 100, UPMGA \n {treeFileName}"
    colorRoot = "gray"

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=15)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    fig = plt.figure(figsize=(15,20), dpi=100)
    axes=fig.add_subplot(1,1,1)

    treedata.root.color = colorRoot
    treedata.root.width = 2
    for Clade in treedata.find_clades():
        if Clade.is_terminal():
            Clade.color = PX.BranchColor.from_name('black')
            for i in nox1:
                if i==Clade.name: 
                    Clade.color = PX.BranchColor.from_hex('#7fc97f') #green
            for i in nox2:
                if i==Clade.name: 
                    Clade.color = PX.BranchColor.from_hex('#f0027f') #pink
            for i in nox3Names:
                if i==Clade.name: 
                    Clade.color = PX.BranchColor.from_hex('#ffff99') #yellow
            for i in nox4:
                if i==Clade.name: 
                    Clade.color = PX.BranchColor.from_hex('#386cb0') #blue
            for i in nox5:
                if i==Clade.name: 
                    Clade.color = PX.BranchColor.from_hex('#fdc086') #peach
            for i in noxD:
                if i==Clade.name: 
                    Clade.color = PX.BranchColor.from_hex('#beaed4') #lavender
        else:
            Clade.color = PX.BranchColor.from_name('gray')


    plt.title(titleName, fontsize=BIGGER_SIZE)
    legendLines = [ln([0], [0], color=('#7fc97f'), lw=2), 
        ln([0], [0], color=('#f0027f'), lw=2), ln([0], [0], color=('#ffff99'), lw=2),
        ln([0], [0], color=('#386cb0'), lw=2), ln([0], [0], color=('#fdc086'), lw=2),
        ln([0], [0], color=('#beaed4'), lw=2)]
    plt.legend(legendLines, ['NOX1', 'NOX2', 'NOX3', 'NOX4', 'NOX5', 'DUOX'], loc="lower left")
    Phylo.draw(treedata, axes=axes, branch_labels=None)
    plt.savefig(treeFigure, dpi=100)
    return 

plotTree(treeFilePhylo, figFile, nox1Names, nox2Names, nox3Names, nox4Names, nox5Names, noxDNames)


# %% Plot bootstrapped trees using a def function and importing xml file from Poseidon processing
# this one has the colored tree and is used most recently!
#This code requries a set of bootstrapped trees to be input
#import the packages specific to this section
from Bio import Phylo
from Bio.Phylo import PhyloXML as PX

#Import tree file. the filename can be changed to draw different trees
path = os.getcwd()
figFileName = "treeFig2.png" #this does not work, currently
figFile = path+"/graphs"+figFileName
treeFileName = "treeFileCoralNox_NJ100.xml" #Change filename here
treeFile = path+"/"+treeFileName
treeFilePhylo = Phylo.parse(treeFile, 'phyloxml')

#Import files with names only, no sequences for each color desired
nox1File = path+"/treeCoralSeqOnlyNox1"
nox1Names = pd.read_csv(nox1File, squeeze=True)
nox2File = path+"/treeCoralSeqOnlyNox2"
nox2Names = pd.read_csv(nox2File, squeeze=True)
nox3File = path+"/treeCoralSeqOnlyNox3"
nox3Names = pd.read_csv(nox3File, squeeze=True)
nox4File = path+"/treeCoralSeqOnlyNox4"
nox4Names = pd.read_csv(nox4File, squeeze=True)
nox5File = path+"/treeCoralSeqOnlyNox5"
nox5Names = pd.read_csv(nox5File, squeeze=True)
noxDFile = path+"/treeCoralSeqOnlyNoxD"
noxDNames = pd.read_csv(noxDFile, squeeze=True)

#def function to plot trees
def plotTree(trees, treeFigure, nox1, nox2, nox3, nox4, nox5, noxD):
    #Title and root color
    titleName = f"NOX Protein Sequences \n Bootstrap 100, Neighbor Joining \n {treeFileName}"
    colorRoot = "gray"
    widthRoot = 2

    #Set sizes of fonts
    SMALL_SIZE = 8
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 30
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=15)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    #Set figure size and dpi
    fig = plt.figure(figsize=(15,20), dpi=100)
    axes=fig.add_subplot(1,1,1)

    #For each tree within the bootstrapped trees, do the following
    for treedata in trees:
        #Set color and width of root. variables are changed above
        treedata.root.color = colorRoot
        treedata.root.width = widthRoot

        #Comment in and out below for Root options
        #can use either inner 66 or midpoint for data that has both human duox1 and duox2
        #Inner66 for upgma and Inner1 for NJ for Hum files

        #treedata.root_with_outgroup('Inner65')
        #treedata.root_at_midpoint()

        #for each clade, change color based on imported names in file
        #colors are in hex format. 
        #Can also use BranchColor.from_name('name') or BranchColor('RGB')
        for Clade in treedata.find_clades():
            if Clade.is_terminal():
                Clade.color = PX.BranchColor.from_name('black')
                for i in nox1:
                    if i==Clade.name: 
                        Clade.color = PX.BranchColor.from_hex('#7fc97f') #green
                for i in nox2:
                    if i==Clade.name: 
                        Clade.color = PX.BranchColor.from_hex('#f0027f') #pink
                for i in nox3Names:
                    if i==Clade.name: 
                        Clade.color = PX.BranchColor.from_hex('#ffff99') #yellow
                for i in nox4:
                    if i==Clade.name: 
                        Clade.color = PX.BranchColor.from_hex('#386cb0') #blue
                for i in nox5:
                    if i==Clade.name: 
                        Clade.color = PX.BranchColor.from_hex('#fdc086') #peach
                for i in noxD:
                    if i==Clade.name: 
                        Clade.color = PX.BranchColor.from_hex('#beaed4') #lavender
            else:
                Clade.color = PX.BranchColor.from_name('gray')

    #title and legend
    plt.title(titleName, fontsize=BIGGER_SIZE)
    legendLines = [ln([0], [0], color=('#7fc97f'), lw=2), 
        ln([0], [0], color=('#f0027f'), lw=2), ln([0], [0], color=('#ffff99'), lw=2),
        ln([0], [0], color=('#386cb0'), lw=2), ln([0], [0], color=('#fdc086'), lw=2),
        ln([0], [0], color=('#beaed4'), lw=2)]
    #plt.legend(legendLines, ['NOX1', 'NOX2', 'NOX3', 'NOX4', 'NOX5', 'DUOX'], loc="upper right")
    
    #ladderize tree - this places the closest node on top
    treedata.ladderize()

    #Draw tree
    Phylo.draw(treedata, axes=axes, branch_labels=None)

    #Save figure, this does not work currently
    plt.savefig(treeFigure, dpi=100)
    return 

#Call function to plot tree with tree data, figure file, and names for the colors
plotTree(treeFilePhylo, figFile, nox1Names, nox2Names, nox3Names, nox4Names, nox5Names, noxDNames)


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
