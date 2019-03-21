<font size=20>__MetaMiner 1.0  Manual__</font>

->Liu Cao, Alexey Gurevich, Hosein Mohimani<-

* [About MetaMiner](#sec_intro)
* [Input and output](#sec_inout)
    * [Input file format](#sec_input)
    * [Output reports](#sec_output)  
* [Pipeline](#sec_pipline)
* [Installation](#sec_installation)
    * [Download MetaMiner](#sec_install_metaminer)
    * [Download SPAdes](#sec_install_spades)
    * [Download antiSMASH](#sec_install_antismash)
    * [Download BOA](#sec_install_antismash)
* [How to run MetaMiner](#sec_howto)
    * [Run MetaMiner with raw nucleotides sequence files (.fasta)](#sec_howto_fasta)
    * [Run MetaMiner with raw reads files (.fastq)](#sec_howto_fastq)
    * [Run MetaMiner with antiSMASH results (.final.gbk)](#sec_howto_antismash)
    * [Run MetaMiner with BOA results (.fasta)](#sec_howto_boa)
    * [Important MetaMiner parameters](#sec_para)
* [Citation](#sec_citation)
* [Feedback and bug reports](#sec_feedback)

<a name="sec_intro"></a>
# About MetaMiner

MetaMiner is a metabologenomic pipeline which integrates metabolomic (tandem mass spectra) and genomic data to identify novel **Ri**bosomally synthesized and **P**ost-translationally modified **P**eptides (RiPPs) and the biosynthetic gene clusters encoding them.

MetaMiner is developed in collaboration of [Carnegie Mellon University](http://mohimanilab.cbd.cmu.edu) (PA, USA), [Saint Petersburg State University](http://cab.spbu.ru) (Russia) and [University of California San Diego](http://cseweb.ucsd.edu/~ppevzner/) (CA, USA) under the Apache 2.0 License. The latest version will be updated in the **N**atural **P**roduct **D**iscovery **tool**kit **NPDtools** at <https://github.com/ablab/npdtools>.

<a name="sec_inout"></a>
# Input and output

MetaMiner takes paired metabolomic and genomic data as input, and output a report summaring all the RiPPs that MetaMiner detected.   

## Input files format
For metabolomic data, MetaMiner work with liquid chromatography–tandem mass spectrometry data (LS-MS/MS). Spectra files must be centroided and in an open spectrum format (MGF, mzXML, mzML or mzData). MetaMiner supports MGF (Mascot Generic Format) and uses msconvert utility from the ProteoWizard package to convert spectra in other formats to MGF.

For genomic data, MetaMiner uses either raw nucleotide sequences or specific genome mining tools' output:

* raw nucleotide sequences `.fasta` format (a high-quality reference or a draft assembly) 
* *antiSMASH*'s `.final.gbk` file
* *BOA*'s protein `.txt` file

For users who only have DNA short read files (`.fastq`), you can first assemble reads with *SPAdes* or *metaSPAdes*. We provide a brief tutorial about how to assemble DNA short reads to nucleotide sequences using *SPAdes* in section [Run MetaMiner with raw read files](#sec_howto_fastq).

## Output reports

All the detected RiPPs are reported in plain text tab-separated value files (`.tsv`). 
Each file starts with a header line containing column descriptions. 
The rest lines represent compound–spectrum matches which include information about both the corresponding mass spectrum and the compound. The columns in the report includes :
-  `SpecFile` (filepath of the spectra file)
-  `Scan` (scan number of the identified spectrum inside the spectra file)
-  `SpectrumMass` (mass of the spectrum in Daltons)
-  `Retention` (retention time of the spectrum in seconds)
-  `Charge` (charge of the spectrum)
-  `Score` (score of the compound–spectrum match)
-  `P-Value` (statistical significance of the compound–spectrum match)
-  `FDR` (estimated FDR at the corresponding P-Value level)
-  `PeptideMass` (mass of the compound in Daltons)
-  `SeqFile` (filepath of the genome sequence file)
-  `Class` (class of the identified RiPP compound)
-  `FragmentSeq` (raw initial sequence of the identified compound)
-  `ModifiedSeq` (the sequence of the identified compound with all applied modifications) 

<a name="sec_pipeline"></a>
# Pipline

MetaMiner pipeline is as follows:

![alt text](https://github.com/mohimanilab/MetaMiner/blob/master/Figure1_combined_v3.png "MetaMiner pipeline")

Please refer the MetaMiner paper for details.

<a name="sec_installation"></a>
# Installation

MetaMiner requires pre-installation of Python 2.7 on a 64-bit Linux system or macOS.  It also requires GNU sed to be present in the PATH environment variable as sed (this is always true for Linux systems but may require additional configurations on macOS since GNU sed is usually installed there as gsed). 

There is no need for installation. Users can directly download and run the binaries.

<a name="sec_install_MetaMiner"></a>
##  Download MetaMiner 

For Linux users, to download [MetaMiner Linux binaries](https://github.com/ablab/npdtools/releases/download/npdtools-2.3.0/NPDtools-2.3.0-Linux.tar.gz) and extract them, go to the directory in which you wish NPDtools to be installed and run:

``` bash
wget https://github.com/ablab/npdtools/releases/download/npdtools-2.3.0/NPDtools-2.3.0-Linux.tar.gz
tar -xzf NPDtools-2.3.0-Linux.tar.gz
cd NPDtools-2.3.0-Linux
```

For mac users, to download [MetaMiner macOS binaries](https://github.com/ablab/npdtools/releases/download/npdtools-2.3.0/NPDtools-2.3.0-Darwin.tar.gz) and extract them, go to the directory in which you wish NPDtools to be installed and run:

``` bash
curl -L https://github.com/ablab/npdtools/releases/download/npdtools-2.3.0/NPDtools-2.3.0-Darwin.tar.gz -o NPDtools-2.3.0-Darwin.tar.gz 
tar -xzf NPDtools-2.3.0-Darwin.tar.gz
cd NPDtools-2.3.0-Darwin
```

<a name="sec_install_spades"></a>
## Download SPAdes binaries

SPAdes will help assembly DNA short reads (`.fastq` file) to nucleotide sequences. For users who have already had `.fasta` or `.final.gbk` file as input, there is no need to download SPAdes.

For Linux users, to download SPAdes:

```bash
wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz
tar -xzf SPAdes-3.13.0-Linux.tar.gz
cd SPAdes-3.13.0-Linux/
```

For macOS binaries:
```bash
curl http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Darwin.tar.gz -o SPAdes-3.13.0-Darwin.tar.gz
tar -zxf SPAdes-3.13.0-Darwin.tar.gz
cd SPAdes-3.13.0-Darwin/
```

<a name="sec_install_antismash"></a>
## Download antiSMASH
After installing conda with Python 2.7, users can easily install antiSMASH 4 as follows:

```bash
# add bioconda channel and its dependencies to conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# install antiSMASH and create antiSMASH environment
conda create -n antismash antismash
source activate antismash
download-antismash-databases
source deactivate antismash
```

<a name="sec_install_boa"></a>
## Download BOA

After installing conda with python 2.7, users can create an environment for BOA and install its dependencies.

```bash
# create BOA environment in conda and install its dependencies.
conda create -n boa python=2.7
source activate boa
conda install python-biopython python-matplotlib python-panda python-numpy nltk clustalw cd-hit hmmer
conda install -c bioconda bx-python
source deactivate boa
```

Then users can directly download BOA from its [github](https://github.com/idoerg/BOA) repository.

<a name="sec_howto"></a>
# How to run MetaMiner?

<a name="sec_howto_fasta"></a>
## Run MetaMiner with raw nucleotides sequence files (.fasta)

For users who have already obtained the assembled genome/metagenome, a sample run of MetaMiner may look like this:
```bash
python metaminer.py test_data/msms/ -s test_data/genome/fasta/ -o metaminer_outdir
```
In this case, all spectra files in `test_data/msms/` will be searched against 
all sequence files in `test_data/fasta/`. In this particular case, it is a search of `test_data/msms/AmfS.mgf` spectrum
against `test_data/fasta/AmfS.fasta` genome fragment. The search mode (considered RiPP class) is 'lantibiotic' (by default).
The identification results will be saved in `metaminer_outdir`. 
The search is performed with all default parameters, 
see the [corresponding subsection](#sec_para) for the default values and available options.

If the run is finished correctly, you will see identification of a lantibiotic with "TGSQVSLLVCEYSSLSVVLCTP" original sequence 
and "T-18GS-18QVS-18LLVCEYS-18SLSVVLCTP" sequence after modifications in `metaminer_outdir/significant_matches.tsv`.
The modifications "T-18" and "S-18" correspond to dehydrobutyrine and dehydroalanine, respectively.
These sequences correspond to AmfS peptide, you may read more about it in [Ueda et al, 2002](https://www.ncbi.nlm.nih.gov/pubmed/11844785).

<a name="sec_howto_fastq"></a>
## Run MetaMiner with raw reads files (.fastq)

For users who have raw DNA short read files (`.fastq`), it is convinient to assemble DNA short reads into nucleitides sequences with *SPAdes* or 'metaSPAdes'. An example of `SPAdes` is as follows:

```bash
# assemble paired-end short reads using SPAdes
python /path/to/SPAdes/bin/spades.py -1 SRR3309439_R1.fastq -2 SRR3309439_R2.gz -o spades_outdir
# take contigs.fasta as input and run MetaMiner
python metaminer.py test_data/msms/ -s spades_outdir/contigs.fasta -o metaminer_outdir
```

where `reads1.fastq.gz` and `reads2.fastq.gz` are paired-end reads file, and `spades_outdir` is the directory saving the genome assembly output. Users can use either `output_dir/contigs.fasta` or `output_dir/scaffolds.fasta` as input of 'metaminer'

<a name="sec_howto_antismash"></a>
## Run MetaMiner with antiSMASH results (.final.gbk)
For users who have obtained antiSMASH results, you can directly apply MetaMiner to the gene bank file (`.final.gbk`). An example is as follows:

```bash
# take contigs.fasta as input and run MetaMiner
python metaminer.py test_data/msms/ -s test_data/antismash/ --antismash -o metaminer_outdir
```

where `--antismash` specifies that `metaminer.py` will only look for `.final.gbk` files in the sequence directory `test_data/antismash` indicated by `-s`. The `.final.gbk` file in the test data folder is generated from `contigs.fasta` by `antiSMASH`. While `MetaMiner` successfully detect AmfS using the `contigs.fasta ` file, it fails with antiSMASH result as input. 

<a name="sec_howto_boa"></a>
## Run MetaMiner with BOA results (.txt)

For users who have obtained BOA results, you can apply MetaMiner to the gene annotation file (`.txt`), which contains a list of annotated genes with the peptide sequence. An example is as follows:

```bash
# extract peptide sequences from BOA result file test.annotated.txt
python boa2fasta.py test.annotated.txt test_data/boa/
# apply metaminer to the .fasta files, each of which contains a peptide sequence
python metaminer.py test_data/msms/ -s test_data/boa/ --boa -o metaminer_outdir
```

where `--boa` specifies that `metaminer.py` will only look for protein sequence files`.fasta` in the sequence directory `test_data/boa/` indicated by `-s`. 

<a name="sec_para"></a>
## Important MetaMiner parameters

Use the following command to see available parameters and details of the options:

```bash
python /path/to/metaminer/metaminer.py
```

Here are a list of frequently used parameters of MetaMiner:

`-s <path>` (or `--sequence <path>`)  
    Path to a sequnce file or to a directory with multiple sequence files inside.
    In the latter case, NPDtools recursively walks through the directory and picks up all files 
    with appropriate extensions (by default: `.fna`, `.fasta`, or `.fa`; case insensitive). 
    You can specify an unlimited number of input sequence files/directories, 
    they will be processed independently (see also `--correspondence` option below).
    *By default, we assume that sequences are raw nucleotide input but this can be modified by specific options (see below)*.
    **At least one sequence file is required**.

`-C <filepath>` (or `--correspondence <filepath>`)  
    Path to a file describing correspondence between sequence and spectra files. 
    The file should be tab-separated and has two columns listing basenames of spectra and sequence filepaths.
    If not provided, the all-vs-all analysis will be performed.

`-a` (or `--antismash`)  
    Sequence files are antiSMASH output (`.final.gbk`). If not specified, the input files are expected to 
    be raw genome nucleotide sequences in FASTA format (see also `--boa` option). Tested with antiSMASH v.2 output.
       
`--boa`                   
    Sequence files are BOA output (protein `.fasta`). If not specified, the input files are expected to 
    be raw genome nucleotide sequences in FASTA format (see also `--antismash` option).
    
`-c <class>` (or `--class <class>`)  
    Class of RiPPs to look for. Valid choices are: 'formylated',
    'glycocin', 'lantibiotic', 'lap', 'lassopeptide', 'linaridin',
    'proteusin', 'cyanobactin', and 'methanobactin'. You can also specify 'all' to try all classes one by one.
    *The default value is 'lantibiotic'*.  


<a name="sec_citation"></a>
# Citation
Please cite Cao et al, 
*MetaMiner: A Peptidogenomics Approach for the Discovery of Ribosomally Synthesized and Post-translationally Modified Peptides,* 
Submitted, *2019* (initial version of the paper is [available on bioRxiv](https://www.biorxiv.org/content/early/2017/12/03/227504)).

<a name="sec_feedback"></a>
# Feedback and bug reports
Your comments, bug reports, and suggestions are very welcomed. 
They will help us to further improve NPDtools.
You can leave them at [our GitHub repository tracker](https://github.com/ablab/npdtools/issues) 
or sent them via support e-mail: <npdtools.support@cab.spbu.ru>.


