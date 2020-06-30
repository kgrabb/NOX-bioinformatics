# Processing pipeline for NOX proteins

#### Kalina Grabb

#### June 2020

*This pipeline uses Python and bash code on the HPC Poseidon in order to search for NOX-like proteins in coral genes. and then display the NOX-like genes in coral through heatmaps and phylogenetic trees*

---

## Outline

1. [Coral sequence query](#Coral-sequence-query)
   1. [Visualize query results](#Visualize-query-results)
2. [Blast analysis with NOX proteins](#Blast-analysis-with-NOX-proteins) 
   1. [Blast script](#Blast-script)
   2. [Sorting blast results](#Sorting-blast-results)
   3. [Visualizing blast results](#Visualizing-blast-results)
3. [Tree Construction](#Tree-Construction)
   1. [Create multiple sequence alignments](#create-multiple-sequence-alignments)
   2. [Create a tree file](#Create-a-tree-file)
   3. [Visualizing trees](#Visualizing-trees)
4. [Pfam and Hmm](#Pfam-and-Hmm)

---

## Coral sequence query

 In this section, the annotated coral protein sequences from reefgenomics.org were downloaded into a file on Poseidon. Basic commandline bash was then used to query the annotated genomes and determine how often the keywords occured.

Files from reefgenomics.org were obtained using the `wget` code as follows. This then places a file in the current directory on Poseidon

``` bash
wget http://comparative.reefgenomics.org/faa/Coral/Porites_lobata_peptides_100.final.clstr.faa
```



In the same folder on Poseidon with all of the annotated protein files, run the `grep` command in order to search for the number of occurances that each key word. This example searches all `.tsv` files for the term `superoxide` and returns the lines where `superoxide` is found into the files `superoxideHitstsv_v2`, which is located up one directory and in the `ReefGenomeQueryAll` folder. The `-i` option for `grep` ignores case for matching.

``` bash
grep -i superoxide *.tsv > ../ReefGenomeQueryAll/superoxideHitstsv_v2
```



The next step was to create a file that queried each coral species for a set of keywords. In the output file, it will display the keyword, `$j`, then the filename, `$i`,  and then the number of lines where the keyword appeared, as specified by the `-c` option in `grep`.

``` bash
for j in NOX NADPH SUPEROXIDE SOD DUOX CATALASE B-245; 
	do echo "$j" >> ../HitsCounts/HitsCountsReefGenomics; 
	for i in *; 
		do echo "$i" >> ../HitsCounts/HitsCountsReefGenomics; 
		grep -ci $j $i >> ../HitsCounts/HitsCountsReefGenomics; 
	done; 
done
```



#### Visualize query results

Using Matlab, the reef genome query results were plotted into a bar graph. This is in the script "genomeSummaryQueries.m".



---

## Blast analysis with NOX proteins

#### Blast script

**Preparation of NOX and Coral sequences**

This portion of analysis was completed by taking a reference database of NOX protein sequences from Kawahara et al., 2007 and querying the coral sequences. The coral sequences were obtained from different sources on the internet.

First, the reference NOX proteins and coral sequences were prepared by combining all NOX sequences into one file and all coral sequences into a different file. The NOX sequences were also edited by hand to make the the fasta files in the format of:

> \>header_without_spaces: description 
>
> Amino_Acid_Sequence_wihtout_breaks

NOX proteins and coral sequences were combined using bash in Poseidon with the following commands.

``` bash
#while in folder with all nox files, copy nox database
cat noxSeqKawahara2007 > ../../blastSlrm/v2Script/noxSeqv2 
less ../../blastSlrm/v2Script/noxSeqv2 

#while in folder with all coral files (excluding COT files), copy all files
cat * >> ../../blastSlrm/v2Script/coralSeqv2
less ../../blastSlrm/v2Script/coralSeqv2 
```



**Blastp the coral sequences for NOX-like proteins**

Next, blastp was run on Poseidon. There are two ways, first, blastp was run through a slurm script. This took a few seconds. The slurm script was as follows:

```bash
#!/bin/bash

#SBATCH --partition=compute                             # Queue selection
#SBATCH --job-name=blastp_nox_coral_kw_all_v1           # Job name
#SBATCH --mail-type=ALL                                 # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kgrabb@whoi.edu                     # Where to send mail
#SBATCH --ntasks=1                                      # Run on a single CPU
#SBATCH --mem=5gb                                       # Job memory request
#SBATCH --time=00:45:00                                 # Time limit hrs:min:sec
#SBATCH --output=serial_job_%j.log                      # Standard output/error
 
pwd; hostname; date
 
echo "Running blastp for Kawahara NOX and all coral sequences on a single CPU core"

module load python3             # Load the python module
module load bio                 # Load bioinformatics package
module load blast               # Load blast

echo "Running graph script on a single CPU core"
 
makeblastdb -in coralSeqCat -dbtype prot
blastp -query noxSeqCat -db coralSeqCat -out blastp_nox_coral_v1_taboutput.txt -outfmt 6

date
```

The script was run by typing `sbatch "scriptname"` in the `bash` command line on Poseidon.



The second way to do it was throught the `bash` command line:

```bash
module load python3             # Load the python module
module load bio                 # Load bioinformatics package
module load blast               # Load blast

makeblastdb -in coralSeqCat -dbtype prot
blastp -query noxSeqCat -db coralSeqCat -out blastp_nox_coral_v1_taboutput.txt -outfmt 6
```

`outfit 6` formats the output for tabular form.



#### Sort blastp results

**Sort by cutoff values**

Blastp results were sorted by percent ID (column 3), e-value (column 11), and bit score (column 12). The `sort` command sorts alphabetically. These results were piped into a new file.

```bash
awk '{if ($11<=1e-70) print}' blastp_nox_coral_v1_taboutput.txt | sort -nrk 11,11 >> ../v1Analysis/blastp_v1_evalue_sort.txt
awk '{if ($3>=50) print}' blastp_nox_coral_v1_taboutput.txt | sort -nrk 3,3 >> ../v1Analysis/blastp_v1_perID_sort.txt
awk '{if ($12>=300) print}' blastp_nox_coral_v1_taboutput.txt | sort -nrk 12,12 > ../v1Analysis/blastp_v1_bitscore_sort.txt
```

Analysis can be combined with the `&&` command.

``` bash
awk '{if ($11<=1e-70 && $3>=50 && $12>=300) print}' blastp_nox_coral_v1_taboutput.txt | sort -nrk 11,11 >> ../v1Analysis/blastp_v1_evalue_sort.txt
```



**Edit headers and FASTA files**

To remove the carrot in front of the header, either can be used in Poseidon commandline. Output files can be provided to pipe `|` the files into. The `$0` prints the whole line and the first line prints just the headers. `sed` removes the `>`

``` bash
awk -e '/^>.*/ {print $0}' coralSeqCopy >> coralSeqHeaders
awk '$1 ~/^>/' noxSeqKawaharaEdited | sed -e 's/>//'  > noxSeqKawaharaList 
cat ../v2Script/coralSeqv2 | sed 's/>//' >test_coralseq.txt
```

To remove the `:` after the header:

```bash
awk -e '/^>.*/ {print $0}' coralSeqCopy >> coralSeqHeaders
```

To copy the unique circumstances, use this to parse out the headers.

```bash
cut -f1 blastp_v2_bit_pID_ev_sort.txt | sed -e 's/://' | sort -u > qSeqUniq
cut -f2 blastp_v2_bit_pID_ev_sort.txt | sort -u > coralSeqUniq
```

To print line numbers

```bash 
awk '{print NR,$0}' matches_bit_pID_ev > matches_bit_pID_ev.txt
```

To copy specific columns

```bash
cut -f1,2,3,11,12 matches_bit_pID_ev.txt > selected_matches.tmp
```

To cut and sort headers

```bash 
sort selected_matches.tmp -k 3 > selected_matches_sorted.tmp
cut -f2 selected_matches_sorted.tmp > gene_matches.tmp
```

To choose specific headers that are in a category, type the category into a column and call the column, in this case column`$15`. The category in this example is `NOX2`.

``` bash
awk '{if ($15=="NOX2") print}' coralSeqSelectTreeNox > coralSeqSelectTreeNoxOnly2
cut -f3 coralSeqSelectTreeNoxOnly2 > treeCoralSeqOnlyNox2
```



#### Visualizing Blast Results

To visualize general patterns of the results from Blast, heatmaps were created. Heat maps show the relationship between coral species (y-axis) and NOX sequences (x-axis) and the given Blast result (color gradient associated iwth e-value, bit score, or percent ID). To do this, the parsed blast results were imported into Python and a Python script in the desktop program was run as follows.

``` python
# %% Import all programs needed
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

#File import for blast results
import os
path = os.getcwd()
print(path)
print(path+filename)
data = pd.read_csv(path+filename, header=0)
print(data.head(5))

#Names import for coral and nox labels on y- and x-axes
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

# %% Create heatmap for bitscore. for percent ID and e-value, it is done the same way.
#Heat map of the bit score
ax=sb.heatmap(bitArray, xticklabels=qU.nox, yticklabels=coralU.coral, linewidths=.05, 
    cmap="PuBuGn", vmin = 30, vmax = 650)

ax.set_title(f"Bit score \n {filename}", fontsize=15)
ax.set_ylabel("Coral Sequences", fontsize=12)
ax.set_xlabel("NOX Sequences", fontsize=12)
ax.tick_params(labelsize=2)
plt.show()
```





---

## Tree Construction

#### Create multiple sequence alignments

**Run a Python script on Poseidon**

Use the SeqIO package in Biopython to parse out specific sequences for the tree. This portion of the code will be run on Poseidon. The script will be written in Python language, but run by commandline on Poseidon. Here is a brief outline of how to run a Python script on Poseidon:

1. Open Poseidon and logon
2. In the Poseidon `bash` commandline, load Python. 
   1. Type `module load python3/3.6.5`. Python should load, it may take a few minutes, but if no error returns, then it worked.
   2. Check that python is working by typing `python`
   3. The version of Python should show, along with a few lines of information. Python mode is indicated by `>>>` on the commandline. 
   4. To exit Python mode, type `exit()`
3. Create a Python script. 
   1. in `bash` commandline, type `nano {scriptname}` . In the file, type your script and save. Note: unless you put directories in the filenames, it will access and create all files in the current directory where the script is and where the script is run
4. Run the Python scritp. 
   1. In the`bash` commandline (exit Python mode first), type `python {scriptname}`. Note: if the script is in a different directory than where the command is being run, be sure to include the filepath in `{filename}`

**Script to parse out sequences**

Use the `SeqIO` package in Biopython to parse out specific sequences. Start with the coral sequence database `coralSeqCopy` with sequencese in fasta format. Create a file with a list of headers for the specific sequences that you want `treeCoralSeqOnly`. the output file will be `bioCoralSeqRecord`. Write the following script in Python and store it on Poseidon. Use the command `python {scriptname}` to run it, as described above. This script is called `biopython_script1_v2`. 

``` python
import Bio
from Bio import SeqIO
import pandas as pd
inputFile="coralSeqCopy"
outputFile="bioCoralSeqRecord"
coralSeqFile="treeCoralSeqOnly"

seqs = pd.read_csv(coralSeqFile)

count=0
total=0
outputHandle = (open(outputFile,"w"))
for record in SeqIO.parse(inputFile, "fasta"):
        total = total+1
        # if seqs['coralSequence'].str.contains(record.id).any():
        # if seqs['coralSequence'].str.lower().any() == record.id:
        if seqs['coralSequence'].str.contains(record.id+"$").any():     
                count = count +1
                SeqIO.write(record,outputHandle,"fasta")
outputHandle.close()
print(str(count) + " records selected out of " + str(total))

```



**Multiple Sequence Alignment**

Use `muscle` to create a multiple sequence alignement. The input is a `fasta` file and the output is a multiple sequence alignment. This is completed on Poseidon `bash` commandline:

```bash
module load bio
module load muscle
muscle -in bioCoralSeqRecord -out muscleAlignmentCoralSeq.msa
```

To create new multiple sequence alignments with give sequences, add new sequences to the input file. This can be done with `cat` in Poseidon commandline.

``` bash
cat bioCoralSeqRecord bioNoxSeqRecordHajjar.fas > bioCoralNoxSeqRecord
```



#### Create a tree file

Use the multiple sequence alignment to create a tree file. This is completed through a series of Python scripts that are run on Poseidon, as described above. These scripts use the Biopython module.

*Convert msa to phylip-relaxed*

Use the `AlignIO` package in Biopython to convert the multiple sequence alignment format to a phylip-relaxed so that it can be used in the next script. This script is called `biopython_align_convert` and is run in the Poseidon commandline.

```python
import Bio
from Bio import AlignIO
inFile = "muscleAlignmentCoralSeq.msa"
outFile = "muscleAlignmentCoralSeq.phy"

AlignIO.convert(inFile, "fasta", outFile, "phylip-relaxed")
```



**Create tree files**

Use `AlignIO` to to create tree files from the `.phy` file created in the step above. The output file is a `phyloxml` , `.xml`. This format allows the Biopython package to color the nodes and therefore is advised to use. To change to `nexus` or `newick`, change the `phyloxml`. The calculator default is `identity`. The constructor can process neighbor joining, `nj`, or `upgma`. `bootstrap_trees` creates the number of trees indicated so that the tree can then be bootstraped when it is visualized. The following script is called `biopython_align_tree-v2` and is run in the Poseidon commandline.

```python
import Bio
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.Consensus import *

inputFile = "muscleAlignmentCoralSeq.phy"
outputFile = "treeFileCoralNJ.xml"

alignment = AlignIO.read(inputFile, "phylip")

calculator = DistanceCalculator('identity')
constructor = DistanceTreeConstructor(calculator, 'nj')
trees = bootstrap_trees(alignment, 100, constructor)

print(trees)

Phylo.write(trees, outputFile, 'phyloxml')
```



#### Visualizing trees

Using the `.xml` tree file, the tree can be visualized in Python. To do this, the `.xml` file was pulled down from Poseidon onto the local computer and the following code was written in a Python program. This file is called `blastp_v2_analysis_2.py`. 

```python
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
treeFileName = "treeFileNoxHajjar_NJ100.xml" #Change filename here
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

```



---

## Pfam and Hmm

**Building Hmm from Pfam** 

Proteins have domains that define the functional regions of the protein. Pfam is a database of the conserved domains. The Pfam consists of multiple sequence alignments and the probablistic representation is displated in a hidden Markov model (HMM). These alignments are made from seed alignments, which are a smaller subset of alignments known to belong to the domain. Hmm profiles can be used to search against a larger database to find other homologous sequences.

Here, a few Pfams were available for NOX online at pfam.xam.org. Pfam `.seed` files were downloaded to Poseidon using `wget`. A slurm script was written on Poseidon to first build a `hmm` profile and then perform `hmmscan` and `hmmsearch` against coral sequences for three different pfam seeds. 

The slurm script was run by typing `sbatch "scriptname"` in the `bash` command line on Poseidon. The script was as follows:

```bash
#!/bin/bash

#SBATCH --partition=compute                             # Queue selection
#SBATCH --job-name=hmmer_coral_nox1_dom_pfam_v2           # Job name
#SBATCH --mail-type=ALL                                 # Mail events (BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=kgrabb@whoi.edu                     # Where to send mail
#SBATCH --ntasks=1                                      # Run on a single CPU
#SBATCH --mem=50gb                                       # Job memory request
#SBATCH --time=09:00:00                                 # Time limit hrs:min:sec
#SBATCH --output=hmmer_coral_nox1_dom_pfam_v2_%j.log                      # Standard output/error
 
pwd; hostname; date
 
echo "Running hmmer on pfams from website with nox1, fad, nad, and ferric reductase against coral seq from online on a single CPU core"

module load python3             # Load the python module
module load bio                 # Load bioinformatics package
module load hmmer               # Load hmmer

echo "Running hmmer for PF01794"

hmmbuild PF01794.hmm PF01794.seed
hmmpress PF01794.hmm
hmmscan --domtblout coral_PF01794_pfam_scan_domtblout PF01794.hmm coralSeqv2
hmmsearch --tblout coral_PF01794_pfam_tblout PF01794.hmm coralSeqv2

echo "Running hmmer for PF08022"

hmmbuild PF08022.hmm PF08022.seed
hmmpress PF08022.hmm
hmmscan --domtblout coral_PF08022_pfam_scan_domtblout PF08022.hmm coralSeqv2
hmmsearch --tblout coral_PF08022_pfam_tblout PF08022.hmm coralSeqv2

echo "Running hmmer for PF08030"

hmmbuild PF08030.hmm PF08030.seed
hmmpress PF08030.hmm
hmmscan --domtblout coral_PF08030_pfam_scan_domtblout PF08030.hmm coralSeqv2
hmmsearch --tblout coral_PF08030_pfam_tblout PF08030.hmm coralSeqv2
 

date

```

This script took ~1hr to run, however 9hr was requested.



**Processing pfam search and scan results**

Pfam search and scan results appear similar to blast results, yet they tend to be more strict on the cut offs. To analyze the results briefly, the sequence headers for the sequences that appears in the results were parsed out. The following code was run on Poseidon commandline. The input files were the direct results from the script above.

``` bash
sed '/^#/d' coral_PF01794_pfam_search |sort -u | awk '{print $1}'  > coralSeqUniqPF01794Search 
sed '/^#/d' coral_PF01794_pfam_scan |awk '{print $4}' | sort -u  > coralSeqUniqPF01794Scan
```

