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
from ij.gui import NonBlockingGenericDialog, Toolbar, Overlay, TextRoi
from java.awt import Color, Font
from ij.plugin import ImagesToStack, MontageMaker
from ij.process import ColorProcessor
from collections import Counter

def show_generic_dialog(title,text):
	gd = NonBlockingGenericDialog(title)
	gd.addMessage(text)
	gd.showDialog()

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

def keyfun_anim_code(x):
	'''Keyfunction to code .*(WT|Ctrl|Exp)([0-9]+).* pattern for sorting WT..., Ctrl..., Exp...

	Input:
	x: (str) to code

	Output:
	(tuple) of form (A,B) where A is 0=WT, 1=Ctrl, 2=Exp, B is animal #
	'''
	ord_ = ['WT','Ctrl','Exp']
	extracted_ = re.search('(WT|Ctrl|Exp)([0-9]+)',x).groups()
	return (ord_.index(extracted_[0]), extracted_[1])

#Basic Locations of things
roi_base = '/Volumes/Backup_HJS/Papers/UPR_AxonGuidance_2021/Zenodo/imageproc_OB_IF/M28_Ddit3_Nrp2IFs_Rois'
imp_base = '/Volumes/Backup_HJS/Papers/UPR_AxonGuidance_2021/Zenodo/imageproc_OB_IF/M28_Ddit3_Nrp2IFs_ExtractedImages'

#Sections of Interest
#sections_of_interest_base = ['7-30-21_4250_4wkMor28_Exp1_LLateral_16','7-30-21_4250_4wkMor28_Exp1_LMedial_25','7-30-21_4250_4wkMor28_Exp1_RMedial_29','9-16-21_4390_4wkMor28_Exp3_LLateral_2']
sections_of_interest_base = ['7-30-21_4250_4wkMor28_Exp1_LLateral_16','7-30-21_4250_4wkMor28_Exp1_LMedial_25']

#Parse Imps and Annotate
imps_of_interest = [parse_imp(os.path.join(imp_base, i + '.tif')) for i in sections_of_interest_base]

#Spanning LUT
spanning_lut = compute_spanning_lut_v2(imps_of_interest, sat_frac = [15,0.35,0.35,0.35])

#Set up for the montage
channels = ['0110','1000'] #panels in mini-montage
color_dict = {'gfp':Color.green, 'tdt':Color.red}
ft_size = 50

out = []
for i in range(len(imps_of_interest)):
	imp = imps_of_interest[i]
	f_base = sections_of_interest_base[i]
	
	#Set LUT
	imp.setDisplayMode(1)
	autolevel_all_channels_v2(imp,sat = [15,0.35,0.35,0.35])
	#imp.setLuts(spanning_lut) #not needed- comparing within NOT between panels
	#force_lut_update(imp)

	#Create the Individual Montages
	indiv_mont = []
	for (e,i) in enumerate(channels):
		imp.setActiveChannels(i)
		imp.setOverlay(None) #clear any overlays (ie. scale bar)
		
		#Add Overlays to the A647-Only Image
		if e == 1:
			for color_ in color_dict.iterkeys():
				roi_path = os.path.join(roi_base, f_base + '_' + color_ + '.zip')
				if os.path.exists(roi_path):
					#Open RoiManager and parse Rois
					rm = RoiManager.getInstance()
					if not rm:
						rm = RoiManager()
					rm.reset()
					rm.runCommand("Open", roi_path)

					#Convert to Overlay
					overlays = Overlay()
					for index_ in range(rm.getCount()):
						overlays.add(rm.getRoi(index_))
					
					overlays.setStrokeColor(color_dict[color_])
					overlays.setStrokeWidth(5)
					imp.setOverlay(overlays)
					imp = imp.flatten() #overwrite with flattened imp

		#Flatten the Imp
		indiv_mont.append(imp.flatten())
	
	#Make Individual Montage GFP/tdtom:A647
	indiv_mont_stack = ImagesToStack.run(indiv_mont)
	montage_ = MontageMaker()
	indiv_mont_out = montage_.makeMontage2(indiv_mont_stack, 1, len(channels), 1, 1, len(indiv_mont), 1, 2, False)
	out.append(indiv_mont_out)

IJ.run(out[-1], "Scale Bar...", "width=50 height=8 font=28 color=White background=None location=[Lower Right] hide overlay")
out[-1] = out[-1].flatten()

out_stack = ImagesToStack.run(out)
montage_ = MontageMaker()
out = montage_.makeMontage2(out_stack, len(out), 1, 1, 1, len(out), 1, 2, False)
out.show()
