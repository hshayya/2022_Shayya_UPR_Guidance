from __future__ import division
from ij import ImagePlus, ImageStack, CompositeImage, IJ, io
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from loci.formats import ImageReader, MetadataTools
import itertools
import re
from collections import OrderedDict, defaultdict
import os
from ij.plugin.frame import RoiManager
from ij.gui import WaitForUserDialog, Toolbar, NonBlockingGenericDialog
from ij.plugin.frame import RoiManager
import random

#Example code to display Max projections and Z-stacks of 20x glomerular images in a blinded fashion, allow user to select Roi and write red/green pixel values to a custom .pxintz format
#Red/Green pixel values are then used to compute per-glomerulus pearson correlations to quantify overlap as a function of genotype.
#Briefly: recursive search to find new 20x glomerular images -> blinding -> display max projection and allow user to select ROI -> load Z stack and allow user to adjust ROI for each slice, saving px values -> pxintz out

def write_roi_from_maxproj(max_path, roi_out):
	'''Load a Max projection and allow user to select Roi to write to disk

	max_path: (str) path to .tif file to be loaded. 
	roi_out: (str) path where .roi file should be saved

	Returns: (null)
	'''
	#Load the Image and Display as a Composite. 
	imp = CompositeImage(ImagePlus(max_path)) #Expects a Tiff File
	imp.setDisplayMode(1)
	imp.show()

	#Brighten
	for i in range(1,imp.getDimensions()[2] + 1):
		imp.setC(i)
		IJ.run(imp, "Enhance Contrast", "saturated=0.35"); #adjusts LUT not pixel data

	#Allow User to Select ROI
	IJ.setTool(Toolbar.POLYGON)
	WaitForUserDialog('Select the ROI').show()
	my_roi = imp.getRoi()

	#Write to RoiManager and Save to Disk
	rm = RoiManager.getInstance()
	if not rm:
		rm = RoiManager()
	rm.runCommand("reset")
	
	rm.addRoi(my_roi)
	rm.runCommand("Save", roi_out)

	#Clean Up 
	rm.runCommand("reset")
	imp.close()

def get_roi_from_file(roi_path, index_ = 0):
	'''Loads an ROI from File (using the ROI Manager) and returns it as an Roi object
	roi_path: (str) path to .roi file on disk
	index_: (int) index of roi to return, useful in multi_roi arrays

	Returns: (ij.gui.Roi) 
	'''
	#Get a Blank ROIManager
	rm = RoiManager.getInstance()
	if not rm:
		rm = RoiManager()
	rm.runCommand("reset")
	
	#Import ROI to Manager and extract it
	rm.runCommand("Open", roi_path)
	my_roi = rm.getRoi(index_)

	#Clean up
	rm.reset()
	rm.close()

	return my_roi

def convert_zstack_to_pxintz(fpath, roi_obj, fout_full_path):
	'''Worker to create a .pxintz file from a Virtual Stack (Multi-C, Multi-Z) .nd2 File and Base ROI
	fpath: (str) path to .nd2 file on disk. 
	roi_path: (ij.gui.Roi) Roi object to show on initial Z stack. User can refine for each Z stack.
	fout_full_path: (str) full path for .pxintz file to be saved to, INCLUDE the final .pxintz
	
	Returns (null)
	'''
	#Import the Image as a Virtual Stack and Show it as a Composite
	opts = ImporterOptions()
	opts.setVirtual(True)
	opts.setId(fpath)
	
	load_image = BF.openImagePlus(opts)[0]
	load_image.setDisplayMode(1)
	load_image.show()
	
	#Get C and Z Dimensions of Stack
	n_channels_ = load_image.getDimensions()[2]
	n_z_ = load_image.getDimensions()[3]

	#Define fout
	#fname = re.search('.*(?=\\.)', os.path.split(fpath)[-1]).group(0) #no longer doing it this way, just have the user pass this directly
	fout_path = fout_full_path

	#Set up variables for first iteration of loop
	iter_ = 0 #if we are first line of output file (0) or not
	roi_obj_iter = roi_obj.clone() #this will update with each iteration, start with what was input

	#Loop over Z slices and write to file
	with file(fout_path,'w') as fout:
		for z_ in range(1, n_z_ + 1):
			print 'Measuring Z: ' + str(z_) + '/' + str(n_z_)
			
			#Adjust LUT for Slice of Interest
			load_image.setZ(z_)
			load_image.killRoi() #If present, otherwise it will scale based on the ROI!
			labs = []
			for i in range(1,n_channels_ + 1):
				load_image.setC(i)
				IJ.run(load_image, "Enhance Contrast", "saturated=0.35"); #adjusts LUT not pixel data
				labs.append(re.search('(?<=rgb\[255\]=).*?(?=,)',load_image.getProcessor().getLut().toString()).group(0))
			
			#Put the loaded ROI on the image
			load_image.setRoi(roi_obj_iter) 
			
			#Allow user to refine ROI or skip this Z slice
			gd = NonBlockingGenericDialog('Roi Delineation')
			gd.enableYesNoCancel('Proceed','Skip')
			gd.addMessage('Select the Roi, then click proceed. If you want to skip this slice, click skip.')
			gd.hideCancelButton()
			gd.showDialog()
			
			#If dialog 'x'-d break, if slice skipped, skip
			if (gd.wasCanceled()):
				print 'Cancelled on z: ' + str(z_) 
				break 
			elif not gd.wasOKed():
				continue
			
			#Otherwise, get Roi after refinement
			my_roi = load_image.getRoi()

			#Write the Slice index 
			if iter_ == 0: #first z written 
				fout.write('#' + str(z_))
				iter_+=1
			else: #subsequently need line break
				fout.write('\n#' + str(z_))
			
			#Write X,Y 
			X,Y = zip(*[(x.getX(),x.getY()) for x in my_roi])
			fout.write('\n>X\n')
			fout.write(','.join([str(x) for x in X]))
			fout.write('\n>Y\n')
			fout.write(','.join([str(x) for x in Y]))
			
			#Write the pixel values for each channel
			for c_ in range(1, n_channels_ + 1):
				load_image.setC(c_)
				proc = load_image.getProcessor()
				px_vals = [proc.get(int(pt.getX()), int(pt.getY())) for pt in my_roi]
				fout.write('\n>' + labs[c_ - 1] + '\n')
				fout.write(','.join([str(x) for x in px_vals]))

			#Save roi from this slice to display on the next iteration
			roi_obj_iter = my_roi.clone()
			
			print '...done!'
		#
	
	#Close the Image
	load_image.close()
	print 'Pxintz File Complete'

#Find all 20x glomerular Z-stacks in base_
base_ = '/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Experiments'
gloms_ = [os.path.join(base_, dir_, f_) for dir_ in next(os.walk(base_))[1] for f_ in os.listdir(os.path.join(base_, dir_)) if re.search('(Ctrl|WT|Exp).*20[Xx].*nd2$', f_) is not None]

#Set up the Write Directories. 
roi_write_dir = '/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Experiments/Analysis_Glomeruli/Zstacks_Confocal_Quant/maxproj_rois'
pxintz_write_dir = '/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Experiments/Analysis_Glomeruli/Zstacks_Confocal_Quant/pxintz_outs'
#use pxintz_write_dir to determine which glomeruli have already been analyzed (will skip these, only analyze new glomeruli)

#Set up Random File that will be used to blind analysis (no peeking!).
f_link = os.path.join(pxintz_write_dir, 'file.tif')

#Find which glomeruli are new
existing_gloms = [i[:-7] for i in os.listdir(pxintz_write_dir)]
new_gloms = [i for i in gloms_ if os.path.split(i)[1][:-4] not in existing_gloms]
new_gloms = [i for i in new_gloms if re.search('4wk',os.path.split(i)[1]) is None] #needed to remove 4wk images (this example is for p5 data)

#Print new glomeruli, useful to check you've got the right stuff
for i in new_gloms:
	print os.path.split(i)[1]

#Shuffle so that we can't bias with the order
random.shuffle(new_gloms)

for (i,fpath) in enumerate(new_gloms):
	split_path = os.path.split(fpath)
	#print 'Working on ' + split_path[1] + '...iteration ' + str(i+1) + '/' + str(len(new_gloms))
	print 'Working on ' + '...iteration ' + str(i+1) + '/' + str(len(new_gloms)) #stay blinded

	#Find the Max Proj Image and Have User Write ROI
	max_path = os.path.join(split_path[0], 'MAX_' + split_path[1][:-4] + '.tif') #always use this convention mapping Zstack <-> MaxProj. See MaxProj code.
	os.symlink(max_path, f_link) #blind the path, no peeking!
	out_roi = os.path.join(roi_write_dir, split_path[1][:-4] + '.roi')
	write_roi_from_maxproj(f_link, out_roi) 
	os.remove(f_link)

	#Write the Pxintz file
	init_roi = get_roi_from_file(out_roi)
	os.symlink(fpath, f_link[:-4]+'.nd2') #blind the path
	f_out_label = re.search('.*(?=\\.)', os.path.split(fpath)[-1]).group(0)
	convert_zstack_to_pxintz(f_link[:-4]+'.nd2', init_roi, fout_full_path = os.path.join(pxintz_write_dir, f_out_label + '.pxintz'))
	os.remove(f_link[:-4]+'.nd2')
	#
