from __future__ import division
from ij import ImagePlus, CompositeImage, IJ, gui
from ij.plugin import ImagesToStack, MontageMaker
import os
import re
from collections import defaultdict
import string, math, random
import itertools
import csv
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from ij.gui import NonBlockingGenericDialog, Roi, Overlay
from java.awt import Color
from ij.plugin.frame import RoiManager

#This code was used to make 4x Montages for M71 and Mor28 Perk Experiments 
#We observe a severe reduction of tdtom+ fibers at 4wks in Perk cKO mice 
#To allow these cross-timepoint comparisons, a single LUTs was computed spanning both p5 and 4wk data.
#The single LUT was applied to ALL WT/Het/cKO mice AT BOTH AGES for the figures.

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

#Must level p5 and 4wk data together to see the effect here 
dirs_ = {'p5':'/path/to/p5/data',
		 '4wk':'/path/to/4wk/data'} #each directory contains WT/Het/cKO images for the appropriate timepoint and an OR (M71 or Mor28)

files = {k:[os.path.join(v,i) for i in os.listdir(v) if re.search('4xTiled\-1.tif$',i) is not None] for k,v in dirs_.iteritems()}

#Compute Spanning LUT across p5 and 4wk Images
lut_ = compute_spanning_lut_v2([parse_imp(f) for v in files.itervalues() for f in v], sat_frac = 0.75)
print lut_

#Initialize RoiManager
rm = RoiManager.getInstance()
if not rm:
	rm = RoiManager()

for k in ['p5','4wk']:
	#Sort Files
	f_ls = files[k]
	f_ls.sort(key = lambda x: ['WT','Ctrl','Exp'].index(re.search('WT|Ctrl|Exp',x).group()))

	#Load RoiManager
	rm.reset()
	rm.runCommand('Open',os.path.join(dirs_[k],'4x_Roiset.zip')) #Contains traced OB outlines for images at the given age, index 0 = WT, 1=Het, 2=cKO.

	#Prep the Images
	montages_per_animal = []
	for e,i in enumerate(f_ls):
		#Set up Image
		imp = parse_imp(i)
		imp.setLuts(lut_)
		force_lut_update(imp)

		#Add Overlay
		overlay = Overlay(rm.getRoi(e))
		overlay.setStrokeColor(Color.white)
		overlay.setStrokeWidth(5)
		imp.setOverlay(overlay)

		#Add Individual Images
		out_ls = []
		for channels in ['11','01','10']:
			imp.setActiveChannels(channels)
			out_ls.append(imp.flatten())

		#Make Montage
		out_stack = ImagesToStack.run(out_ls)
		montage_ = MontageMaker()
		out = montage_.makeMontage2(out_stack, len(out_ls), 1, 1, 1, len(out_ls), 1, 2, False)
		out.setTitle('')
		montages_per_animal.append(out)

	#Make the Montages
	IJ.run(montages_per_animal[-1], "Scale Bar...", "width=500 height=8 font=28 color=White background=None location=[Lower Right] hide overlay");
	montages_per_animal[-1] = montages_per_animal[-1].flatten()
	
	montages_per_animal_stack = ImagesToStack.run(montages_per_animal)
	final_out = montage_.makeMontage2(montages_per_animal_stack, 1, len(montages_per_animal), 1, 1, len(montages_per_animal), 1, 2, False)
	final_out.setTitle(k)
	final_out.show()