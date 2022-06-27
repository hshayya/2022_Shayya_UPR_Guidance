OB IF Analysis
*This folder covers two experiments: 
(1) M28 IF in p5 Mor28 Ddit3 WT vs. cKO OBs
(2) Nrp2 IF in p28 Mor28 Ddit3 cKO OBs

*The computation here is trivial (subtracting means, only done for Nrp2 experiment). 
*The complexity comes from taking huge stacks of OB sections imaged at 20x, finding the glomeruli of interest, and extracting pixel values for analysis. 

*Below is a summary of the analysis, and what has been provided on Github/Zenodo as examples. 
*Contact hani (hjs2148 at columbia dot edu) if you have questions or need something that isn't here.

Additional Notes:
(filename): indicates file on github (/2022_Shayya_UPR_Guidance/imageproc_OB_IF/filename)
[filename]: indicates file on Zenodo (/imageproc_OB_IF/filename)

Imaging/Analysis Pipeline: 
1. OB sections imaged using Nikon Ti2E microscope (running NIS Elements) at 4x magnification. No autofocus. 
  *In initial experiments, DAPI channel used -> goal is to mark X,Y coordinates for each section for step 2
  *In follow up experiments, used GFP/tdtomato channels -> marked only X,Y coordinates for sections containing glomeruli for step 2. Required more "hands-on" time but less overall scope time, smaller final dataset. 
  
2. Sections of interest (found in 1) marked in NIS Elements and re-imaged at 20x (autofocus used to find Z) -> saved to large .nd2 files. 
  *These .nd2 files can be somewhat complex: contain images from up to 4 slides. Each slice is an entire OB section
  *A given .nd2 file could contain sections from >1 animal since multiple animals were sometimes imaged on the same day. Could also have one animal spread out over >1 .nd2 file (pending #s of slide, scope time etc.).
  *Every slide imaged contained two rows of OB sections, one row of which was flipped (this helped with the mechanics of the sectioning). So ~1/2 of the OBs in these large .nd2 files were upside down!
  
3. (process_nd2_to_slices.py) Initial Image Processing of large .nd2 files. 
  *During imaging setup, I recorded which OBs would need to be flipped (see above) and which OBs map to which animals. 
  *Used these annotations & jython to parse large .nd2 files -> downsample each slice, flip if needed, name based on slice # in animal and save in folder for that animal.
  *Example script on Github is for a simple case: 3 .nd2 files with slices from the same animal. 
  **[Example_Slice_Lookup_Table.tsv] example of intermediate lookup table created by this script, used in later steps
      *slice_order: slice # in .nd2 file 
      *out_order: slice # in the final dataset for this animal
      
4. Using the downsampled dataset from 3, marked Rois for LLateral, RLateral, LMedial, RMedial glomeruli. Rois were marked for each section that the glomerulus was present. IJ>Measure tool used to save X/Y/Slice coordinates as csv
  ** [Example_Glomeruli.csv] example of the output from this step
  
5. (5_annotate_rois_for_extraction.R): annotates the X/Y/Slice measurements from [Example_Glomeruli.csv] w/ glomerular id (LLateral, RLateral...) and incorporates information from [Example_Slice_Lookup_Table.tsv] to map back to original (full-res) .nd2 file.
  ** [Example_Glomeruli_Annotated.csv]: example of output from this step
  
6. (6_extract_glomerular_rois.py): steps 1-5 are done independently for all animals analyzed. Now, everything is combined. This script reads the [Example_Glomeruli_Annotated.csv] table generated for all the animals and extracts each Roi as a small .tif file into a new directory. 
  ** All the extracted images are provided on Zenodo
      ***[M28_Ddit3_Nrp2IFs_ExtractedImages/*tif]: extracted images for Nrp2 staining in 4wk M28 Ddit3 cKO mice
      ***[M28_Ddit3_M28IFs_ExtractedImages/*tif]: extracted images for the M28 staining in p5 M28 Ddit3 WT and cKO mice
      
7. (7_analyzeOBs_make_montage.py): Iterates through all the extracted images, allows user to select GFP/tdtomato Rois, then extracts those pixels as custom .pxintz files (done only for Nrp2 IFs) and creates a montage of images for each animal
  ** All the Rois are provided on Zenodo
      ***[M28_Ddit3_Nrp2IFs_Rois/*zip]: Rois for Nrp2 staining in 4wk M28 Ddit3 cKO mice
      ***[M28_Ddit3_M28IFs_Rois/*zip]: Rois for M28 staining in p5 M28 Ddit3 WT and cKO mice
  ** PxintZ Outs are Provided for Nrp2 experiment
      ***[M28_Ddit3_Nrp2IFs_PxintZs/*pxintz]: custom pxintz files for Nrp2 staining in 4wk M28 Ddit3 cKO mice
  ** Montages are provided on Zenodo showing all glomeruli imaged
      ***[M28_Ddit3_Nrp2IFs_Montages/*jpg]: montages showing the full dataset of Nrp2 staining in the 4wk M28 Ddit3 cKO mice
      ***[M28_Ddit3_M28IFs_Montages/*jpg]: montages summarizing the full dataset of M28 staining in the p5 M28 Ddit3 cKO and WT OBs.
      
8. (8_quant_nrp2_ifs.R): Analyzes the 4wk M28 Ddit3 cKO Nrp2 IF pxintz files, taking sections where there is both a gfp and tdtomato glomerulus and comparing the Nrp2 levels. See methods. 

9. (9_montage_for_paper_nrp2.py): Prepares a montage of representative images for the p28 Mor28 Ddit3 cKO Nrp2 IF experiment

10. (10_montage_for_paper_m28.py): Prepares a montage of representative images for the p5 Mor28 Ddit3 WT vs. cKO M28 IF experiment. 