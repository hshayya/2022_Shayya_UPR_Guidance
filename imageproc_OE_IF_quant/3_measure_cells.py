import os
import re
import math
import random
from ij import ImagePlus, CompositeImage
from ij import IJ, gui
from ij import io
from ij.plugin import ImagesToStack, MontageMaker
from ij.gui import NonBlockingGenericDialog, Overlay, TextRoi
from ij.plugin.frame import RoiManager
from java.awt import Color, Font
from ij.process import ColorProcessor, ImageProcessor

def parse_imp(fpath, rotate = False, autolevel = False):
	if fpath[-3:] == 'tif':
		parsed = CompositeImage(ImagePlus(fpath))
	elif fpath[-3:] == 'nd2':
		opts = ImporterOptions()
		opts.setId(fpath)
		parsed = BF.openImagePlus(opts)[0] #creates array
	else:
		return None

	parsed.setDisplayMode(1)
	if rotate:
		IJ.run(parsed, "Rotate 90 Degrees Right", "")
	if autolevel:
		autolevel_all_channels(parsed, 0.35)
	return parsed

def segment_cell(imp, channel_to_segment, channel_to_measure):
	'''Segment an image of a single cell
	imp: (ImagePlus) to segment. 
	channel_to_segment: (int) channel to use for segmentation. 1-based (uses .setC() method)
	channel_to_measure: (int) channel to measure mean for segmented area. 1-based (uses .setC() method)

	Returns:
	{dict}. 'area':segmented area in cal units^2, 'mean':mean for indicated channel in segmented area, 'roi':(ij.gui.Roi) with segmented area
	'''
	#Intialize RoiManager
	rm = RoiManager.getInstance()
	if not rm:
		rm = RoiManager()
		rm.runCommand("reset")

	#Get Calibration Information
	cal_table = imp.getCalibration()
	um2_per_px2 = cal_table.pixelHeight * cal_table.pixelWidth
	px_center = [i/2 for i in imp.getDimensions()[0:2]] #in pixel

	#Extract correct channel for thresholding
	channel_of_interest = ImagePlus('threshold',imp.getStack().getProcessor(channel_to_segment))
	channel_of_interest.setCalibration(cal_table)

	#Auto-Threshold and Mask 
	IJ.setAutoThreshold(channel_of_interest, "Huang dark");
	IJ.run(channel_of_interest, "Convert to Mask", "");
	IJ.run(channel_of_interest, "Watershed", "");

	#Find Particles
	IJ.run(channel_of_interest, "Analyze Particles...", "size=20-150 clear include add"); #size in um^2.
	#channel_of_interest.show() #was helpful to debug this part/single images
	#imp.show() #was helpful to debug this part/single images (refine list to have 1 image, then run code to this point)
	
	#Check nROI Found
	#print imp.getTitle() + '....' + str(rm.getCount()) + ' rois' #was useful for debugging- find cases with >1 Roi and refine approach to pick best
	if rm.getCount() == 0:
		raise ValueError('0 Rois found for ' + imp.getTitle() + '...')
	
	#Get Areas and Centroids for All Rois
	px_arrays = [[pt for pt in roi] for roi in (rm.getRoi(i) for i in range(rm.getCount()))] #list of lists, each contining java.awt.points for a found roi
	areas = [len(i) * um2_per_px2 for i in px_arrays] #Area ROI = npixels, convert to um^2
	centroids = [[(max(i) + min(i))/2 for i in zip(*[(pt.getX(), pt.getY()) for pt in ls])]for ls in px_arrays] #in px, leave it this way
	dist_to_center = [math.sqrt((x-px_center[0])**2 + (y-px_center[1])**2) for (x,y) in (a for a in centroids)] #in px
	
	#Select the Roi Closest to center of Image
	closest_roi = min(range(len(dist_to_center)), key=dist_to_center.__getitem__)

	#Make Measurements of Deduced Roi
	imp.setC(channel_to_measure)
	proc = imp.getProcessor()
	or_vals = [proc.get(pt.x, pt.y) for pt in px_arrays[closest_roi]]
	mean_val = sum(or_vals)/len(or_vals)
	
	#Return Roi and Measurements
	output_dict = {'area':areas[closest_roi], 'mean':mean_val, 'roi':rm.getRoi(closest_roi)}
	return output_dict
	
#Find All Extracted Files
search_dir = '/path/to/dir/where/extracted/cell/images/saved'
f_ls = [os.path.join(root, f) for (root, dirs, files) in os.walk(search_dir) for f in files if f[-3:] == 'tif']

#Subset (only used for testing)
#f_ls = random.sample(f_ls,100)

#Compute Cell Type and Set Lookup Dict
lookup_dict = {'gfp':3,'tdtom':2}
cell_type = [re.search('gfp|tdtom', i).group() for i in f_ls]

#Make Measurements
fout_path = '/path/to/output.tsv'
with open(fout_path,'w') as f_out:
	f_out.write('\t'.join(['cell','area','mean']))
	for e,(f_in, cell_) in enumerate(zip(f_ls, cell_type)):
		print 'Index ' + str(e) + '/' + str(len(f_ls))
		imp = parse_imp(f_in)
		segmented_out = segment_cell(imp,lookup_dict[cell_],1)
		f_out.write('\n'+'\t'.join([imp.getTitle()[:-4],str(segmented_out['area']), str(segmented_out['mean'])]))
