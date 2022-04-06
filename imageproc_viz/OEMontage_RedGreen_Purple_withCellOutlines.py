import xml.etree.ElementTree as ET
import csv
import os
import re
import math
from loci.plugins.in import ImporterOptions
from loci.plugins import BF
from ij.plugin import ImagesToStack, MontageMaker
from ij import ImagePlus, CompositeImage
from ij import IJ, gui
from ij import io
from ij.plugin import ImagesToStack, MontageMaker
from ij.gui import NonBlockingGenericDialog, Overlay, TextRoi, PolygonRoi, Roi
from ij.plugin.frame import RoiManager
from java.awt import Color, Font
from ij.process import ColorProcessor, ImageProcessor

#Example code to make montages like Fig S4D-E, S5D or S7A
#For each genotype of interest, one panel merged GFP/tdtomato, one panel A647 channel w/ red/green overlays of the relevent cells
#ROIs were extracted from images of the OE, Cells were counted manually using CellCounter plugin and those files saved prior to running this code.

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

def autolevel_all_channels_v2(imp, sat = 0.35):
	'''Autolevel a composite Image, possibily with a different saturation for each channel

	--Input--
	imp (ImagePlus): to autolevel
	sat (int, float or list of int/floats): saturation to pass to Enhance contrast. If a single int/float all channels autoleveled to that. If a list, the saturation for each channel

	--Output--
	Auto-leveled (ImagePlus)
	'''
	imp.setDisplayMode(1)
	if isinstance(sat, list):
		for c,sat_ in zip(range(imp.getDimensions()[2]), sat):
			sat_str = "saturated=" + str(sat_)
			imp.setC(c+1) #1-based...
			IJ.run(imp, "Enhance Contrast", sat_str)
	else:
		sat_str = "saturated=" + str(sat)
		for c in range(imp.getDimensions()[2]):
			imp.setC(c+1) #1-based...
			IJ.run(imp, "Enhance Contrast", sat_str)
	return imp

def compute_spanning_lut_v2(ls_imps, sat_frac = 0.35):
	'''Places list of Composite Images on same scale by first autoleveling all channels for all images, then computing min(lut), max(lut) for each channel across all images
	Unlike v1, can use different saturations for each scale
	
	Input:
	ls_imps: (list) of (CompositeImages). Assumes have the same #s of channels in all cases. 
	sat_frac: (int, float or list of int/floats) passed to autolevel_all_channels_v2. See that function for details.

	Return:
	(array of LUTs) containing the spanning LUT for the images  
	'''
	autoscaled_imps = [autolevel_all_channels_v2(i, sat_frac) for i in ls_imps]
	luts = [i.getLuts() for i in autoscaled_imps]
	lut_lims = [(min([z.min for z in i]), max([z.max for z in i])) for i in zip(*luts)]
	output_lut = ls_imps[0].getLuts()
	for e,i in enumerate(output_lut):
		output_lut[e].min, output_lut[e].max = lut_lims[e]

	return output_lut

def force_lut_update(imp):
	'''
	Setting Luts via imp.setLuts() will update the Luts but will not always update the display. 
	This helper function ensures that the displayed image is correctly using the lut in imp.getLuts()
	
	imp: (CompositeImage) to be analyzed, in-place

	Returns: Null
	'''
	lut = imp.getLuts()
	for i in range(imp.getDimensions()[2]):
		imp.setC(i+1) #1-based
		imp.getProcessor().setLut(lut[i])
	imp.updateAndDraw()

def parse_cellcounter_to_dict(fpath):
	'''Parse Cell-Counter Xml file to Dictionary
	Inputs:
	fpath (str) path to xml file on disk

	Values:
	(dict). Keys 'x_cal', 'y_cal' = (float) calibrations in each axis. 
	Keys '1'-'8' = (lists) of tuples containing cell positions in the form (x,y)
	'''
	tree = ET.parse(fpath)
	
	cells_dict = {}
	cells_dict['x_cal'] = float(tree.find('./Image_Properties/X_Calibration').text)
	cells_dict['y_cal'] = float(tree.find('./Image_Properties/Y_Calibration').text)
	
	rt = tree.find('Marker_Data') #re-root the tree
	for type_ in rt.iter('Marker_Type'):
		cells = []
		for marker_ in type_.iter('Marker'):
			cells.append((int(marker_[0].text), int(marker_[1].text)))
		#
		cells_dict[type_.find('Type').text] = cells

	return cells_dict

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

def find_all_single_cell_rois(path_to_imp, path_to_counter, counter_to_color, rotate_bool = False, rotate_angle = 0, px_size_cell_box=60):
	'''Outline single cells on a larger image.
	--Input--
	path_to_imp (str): disk location of image to parse
	path_to_counter (str): disk location of counter file to parse
	counter_to_color (dict): of form {counter#:{'color':Color.xxx, 'channel':channel_#}} mapping cell counter type to color of outline & channel to use for autosegmentation
	rotate_bool (bool): should the parsed image be rotated before starting analysis?
	rotate_angle (str): angle to rotate (careful!! expects string). Only used if rotate_bool = T
	px_size_cell_box (int): pixels for one side of square tile to crop around each cell. Default works well for 20x Spinning Disk

	--Return--
	(dict) of the form {'imp':Parsed_ImagePlus, 'cell_counter_#':Overlay ... for all cell counter types}
	'''
	#Parse Imp and Counter
	imp = parse_imp(path_to_imp)
	if rotate_bool:
		IJ.run(full_image, "Rotate... ", "angle=" + rotate_angle + " grid=1 interpolation=Bilinear stack")
	
	counter_ = parse_cellcounter_to_dict(path_to_counter)
	
	#Crop tiles around single cells -> store images and tile locations in dictionary
	cells_out = {k:{k_inner:list() for k_inner in ['imps','locs']} for k in counter_to_color.iterkeys()}
	for cell_type in counter_to_color.iterkeys():
		for i in range(len(counter_[cell_type])):
			#Set Roi and Crop
			x_topleft = int(counter_[cell_type][i][0] - px_size_cell_box/2)
			y_topleft = int(counter_[cell_type][i][1] - px_size_cell_box/2)
			imp.setRoi(x_topleft, y_topleft, int(px_size_cell_box), int(px_size_cell_box))
			out_image = imp.crop('stack')
	
			#Image
			cells_out[cell_type]['imps'].append(out_image)
			cells_out[cell_type]['locs'].append((x_topleft, y_topleft))
			imp.setRoi(None)
	
	#Outline cell for each tile -> convert Rois back to full-image coordinates
	final_out_dict = {}
	for cell_type,val in cells_out.iteritems():
		out_overly = Overlay()
		for e,(singlecell_imp, roi_locs) in enumerate(zip(val['imps'], val['locs'])):
			print 'Segmenting Cell ' + str(e+1) + '/' + str(len(val['imps']))
			try:
				segmented_roi = segment_cell(singlecell_imp, counter_to_color[cell_type]['channel'], 1)['roi']
				segmented_poly = segmented_roi.getFloatPolygon()
				new_xs = [x + roi_locs[0] for x in segmented_poly.xpoints]
				new_ys = [y + roi_locs[1] for y in segmented_poly.ypoints]
				new_roi = PolygonRoi(new_xs, new_ys, Roi.POLYGON)
				out_overly.add(new_roi)
			except ValueError:
				print 'Segmentation Failed: Type' + cell_type + ' cell ' + str(e+1)  
				singlecell_imp.show()
				continue
		
		out_overly.setStrokeColor(counter_to_color[cell_type]['color'])
		out_overly.setStrokeWidth(1)
		out_overly.setFillColor(None)
		final_out_dict[cell_type] = out_overly

	#Returns
	final_out_dict['imp'] = imp
	return final_out_dict

def add_text_roi(imp,txt,x,y,ft_size,color_, just_=0):
	'''Helper Function to add text roi to image and flatten, using some typical defaults
	--Input--
	imp: (ImagePlus) Not modified by function.
	txt: (str) 
	x,y,ft_size: (int)
	color_: (java.awt.Color)
	just_: (int) (0-left, 1-center, 2-right) justification

	--Returns--
	(ImagePlus) with text burned in. 
	'''
	overlay_ = Overlay()
	out = imp.duplicate() #prevents making changes to the original imp
	txtroi = TextRoi(x, y, txt,Font("SansSerif", Font.PLAIN, ft_size))
	txtroi.setJustification(just_)
	overlay_.add(txtroi)
	overlay_.setFillColor(None)
	overlay_.setStrokeColor(color_)
	overlay_.setStrokeWidth(1)
	out.setOverlay(overlay_)
	return out.flatten()

#Script
base_path = '/Users/hani/Desktop' #root path
montages = ['p5_M71_Perk','p5_M28_Perk'] #subdirs for each genotype
imp_order = ['WT','Ctrl','Exp'] #image order to use for each genotype
counter_lookup = {'6':{'color':Color.yellow,'channel':2},'7':{'color':Color.red,'channel':2},'8':{'color':Color.green,'channel':3}} #mapping btwn cell counter #'s & colors on overlay
settings = [
[[(256,2044),(85,855),(83,453)],[(256,3406),(84,1035),(83,461)]], #M71 Ctrl, Exp Settings B&C max/min for each channel
[[(89,1353),(85,645),(84,332)],[(89,3529),(85,640),(84,274)],[(89,1893),(85,522),(84,332)]] #M28 Wt, Ctrl,Exp B&C max/min for each channel
] #NB no quantiative comparisons being made, no need to force same LUT on all images.

for f_, s_ in zip(montages, settings):
	#Find Images/Counters
	fs = [os.path.join(base_path, f_, i) for i in os.listdir(os.path.join(base_path, f_)) if re.search('^ROI',i) is not None]
	cs = [os.path.join(os.path.split(i)[0], 'CellCounter_' + os.path.split(i)[1][:-3] + 'xml') for i in fs]

	#Sort Images/Counters
	fs.sort(key = lambda x: imp_order.index(re.search('WT|Ctrl|Exp',os.path.split(x)[1]).group()))
	cs.sort(key = lambda x: imp_order.index(re.search('WT|Ctrl|Exp',os.path.split(x)[1]).group()))
	
	#Find and Outline the Cells
	parsed_output = [find_all_single_cell_rois(path_to_imp = f, path_to_counter = c, counter_to_color = counter_lookup, rotate_bool = False, rotate_angle = 0, px_size_cell_box=60) for f,c in zip(fs, cs)]
	
	#Create Panels
	out_monts = []
	for x,lims in zip(parsed_output, s_):
		#Get the Analyzed Image
		imp = autolevel_all_channels_v2(x['imp'], sat = [0.5,0.5,0.5,10])

		#Set the Luts based on Manual Annotations
		lut = imp.getLuts()
		for ch, l_ in zip(lut[:3], lims):
			ch.min, ch.max = l_
		imp.setLuts(lut)
		force_lut_update(imp)

		#Make the Montage
		for cs_ in ['0110','1000']:
			imp.setActiveChannels(cs_)
			imp.setOverlay(None)
			if cs_ == '0110':
				out_monts.append(imp.flatten())
			else:
				out = imp.flatten()
				for k,v in x.iteritems():
					if k != 'imp':
						v.setStrokeWidth(3)
						out.setOverlay(v)
						out = out.flatten()
				out_monts.append(out.flatten())
				
	#Add Scale bar to last one
	IJ.run(out_monts[-1], "Scale Bar...", "width=50 height=5 font=0 color=White background=None location=[Lower Right] hide overlay")
	out_monts[-1] = out_monts[-1].flatten()
	imp_stack = ImagesToStack.run(out_monts)
	montage_ = MontageMaker()
	montage_.makeMontage2(imp_stack, 2,len(fs), 1, 1, len(fs)*2, 1, 2, False).show()
