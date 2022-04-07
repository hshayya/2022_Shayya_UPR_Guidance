This repository contains code associated with the following study:

**ER stress transforms stochastic olfactory receptor gene choice into stereotypic axon guidance programs**

Hani J. Shayya<sup>1,2,3</sup>, Jerome K. Kahiapo<sup>1,3</sup>, Rachel Duffié<sup>1</sup>, Katherine S. Lehmann<sup>4</sup>, Lisa Bashkirova<sup>1</sup>, Kevin Monahan<sup>1</sup>, Ryan P. Dalton<sup>5</sup>, Joanna Gao<sup>6</sup>, Song Jiao<sup>4</sup>, Ira Schieren<sup>1</sup>, Leonardo Belluscio<sup>4</sup>, and Stavros Lomvardas<sup>1,7,8</sup>

<sup>1</sup>Mortimer B. Zuckerman Mind, Brain and Behavior Institute, Columbia University, New York, NY 10027, USA

<sup>2</sup>Medical Scientist Training Program, Vagelos College of Physicians and Surgeons, Columbia University, New York, NY, 10032, USA

<sup>3</sup>Integrated Program in Cellular, Molecular, and Biomedical Studies, Columbia University Irving Medical Center, Vagelos College of Physicians and Surgeons, Columbia University, New York, NY, 10032, USA

<sup>4</sup>Developmental Plasticity Section, National Institute of Neurological Disorders and Stroke, National Institutes of Health, Bethesda, MD, 20892 USA

<sup>5</sup>The Miller Institute for Basic Research in Science, University of California Berkeley, Berkeley, CA 94720 USA

<sup>6</sup>Barnard College, New York, NY, 10025 USA 

<sup>7</sup>Department of Biochemistry and Molecular Biophysics, Columbia University Irving Medical Center, Vagelos College of Physicians and Surgeons, Columbia University, New York, NY, 10032, USA

<sup>8</sup>Department of Neuroscience, Columbia University Irving Medical Center, Vagelos College of Physicians and Surgeons, Columbia University, New York, NY, 10032, USA

**Raw Data:**

GEO: GSE198886

Zenodo: (available after publication)

**Contents:**

***annot***: Create custom GTFs/STAR & Salmon indices used in the study

***align***: Alignment

***bk_rnaseq***: Bulk RNA-seq analysis of OSN differentiation, stress score computation in WT & swap mice, Omp-Cre Ddit3 and M28-iCre/iGFP Ddit3 experiments

***or_cds_clus***: Alignment and clustering of OR protein sequences, visualization

***scrna***: scRNA-seq clustering and DE analysis high/low stress OSNs

***aracneviper***: Prep & run ARACNE-AP network reverse-engineering, run & viz msVIPER high/low stress OSNs

***flowcyt***: Analysis/Viz for OR-iGFP; Atf5-Rep experiments & Omp-GFP;Atf5-Rep;Perk Het experiments

***imageproc_helpers***: Scripts for useful image processing tasks (max projections, montages, downsampling)

***imageproc_viz***: Example scripts showing how figures in the paper were generated (4x OB montages for quantiative/qualitative comparisons, 20x Montages, OE montages w/ red/gree ROIs delineated)

***imageproc_blinding***: Example scripts showing how images of 20x OB glomeruli were blinded prior to intermixed/compartmentalized/adjacent calls, script for blinding OE IF images prior to counting tdtom and GFP cells in M28/M71 Perk mice

***imageproc_gloms_pearson***: Example scripts showing how red/green pixel intensities were extracted in blinded fashion from 20x images of glomeruli & pearson correlations computed to quantify fiber overlap ~ genotype.

***imageproc_OE_IF_quant***: Example scripts showing how cells were extracted from OE IF images, M28 signals measured and local differences in M28 levels GFP - tdtom cells computed.

**Notes**

1) For genomics and flow cytometry folders, see corresponding folders in Zenodo dataset for useful intermediate & output files. The posted code contains commented-out lines that read the intermediate files and generate the final figures in the paper, which will save some time on re-doing the anlaysis. 
2) Please note that some aspects of the genomics (ie. alignment) introduce randomness that will cause the included intermediate files to differ (VERY) slightly from the files you might get re-running the analysis from scratch. This is expected and will not affect the conclusions of the work.
3) Image processing scripts are generally meant to serve as examples to show how the analysis was done. Useful raw, intermediate, and output files are included on Zenodo where possible, but the size of the dataset and complexity of the file-structure here make it so that the scripts provided here will not work "out-of-the-box". Please contact Hani if you'd like anything that isn't posted and we'll figure out a way to get it to you.

***Software Info***

**Genomics**: [sessioninfo_linux.txt](https://github.com/hshayya/2022_Shayya_UPR_Guidance/files/8438755/sessioninfo_linux.txt)

Salmon v.0.13.1

STAR v2.5.3a

samtools 1.7 (using htslib 1.7-6-g6d2bfb7)

wigToBigWig v 4

cellranger v2.1.1-3.1.0 

**Flow Cytometry and Image Processing**:
