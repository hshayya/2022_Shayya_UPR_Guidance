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
		autolevel_all_channels_v2(parsed, 0.35)
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

####################################################################
####                 Find and Organize All Files                ####
####################################################################

#Find all the files
dir_ = '/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Nrp2IFs/Analysis_Glomeruli/Extracted_Images'
files_ls = [i for i in os.listdir(dir_) if i[-3:] == 'tif']

regex = '([0-9]+\\-[0-9]+\\-[0-9]+_(?:(?:B[0-9]+_[LR][1-5])|(?:[0-9]+_4wkMor28))_(?:WT|Ctrl|Exp)[0-9]).*([LR](?:Lateral|Medial))'
annotations = [re.search(regex,i).groups() for i in files_ls]

#Create Empty Dictionary of Form Animal:{Glomerulus:[list of images]}
unique_annotations = [set(i) for i in zip(*annotations)]

anim_dict = {}
for anim in unique_annotations[0]:
	anim_dict[anim] = {}
	for glom in unique_annotations[1]:
		anim_dict[anim][glom] = list()

#Populate dictionary and sort
for (f, (anim, glom)) in zip(files_ls, annotations):
	anim_dict[anim][glom].append(f)

for (a,d) in anim_dict.iteritems():
	for k,v in d.iteritems():
		v.sort(key = lambda x: int(re.search('[0-9]+(?=.tif)', x).group()))

#Record the first image in each glomerulus (LLat...) which will get the lab
first_images = []
for a,d in anim_dict.iteritems():
	for k,v in d.iteritems():
		try:
			first_images.append(v[0])
		except IndexError:
			continue

#Max Images/location 
max_n_images = {a:max([len(c) for c in b.itervalues()]) for a,b in anim_dict.iteritems()}

##############################################################################
####                 Parse Images and Compute Spanning LUT                ####
##############################################################################

#Read in parsed files (right now ignoring top level key: animal, just going through positions).
parsed_files = {}
for k,v in anim_dict.iteritems():
	parsed_files[k] = []
	for pos in ['LLateral','LMedial','RLateral','RMedial']:
		for e in range(max_n_images[k]):
			try: 
				parsed_files[k].append(parse_imp(os.path.join(dir_, v[pos][e])))
			except IndexError:
				parsed_files[k].append('blank')

#Compute spanning lut on all non-blank images.
spanning_lut = {k:compute_spanning_lut_v2([i for i in v if i != 'blank'], sat_frac = 0.35) for k,v in parsed_files.iteritems()}

####################################################################
####               Delimit Rois GFP/tdtomato                    ####
####################################################################

#Write Red & Green Rois for all non-blank images. Will only run on files that aren't there already!
roi_out_base = '/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Nrp2IFs/Analysis_Glomeruli/Extracted_Images/roi_annotations'
for k,v in parsed_files.iteritems():
	for imp in v:
		#If blank, continue
		if imp == 'blank':
			continue
	
		#If exists, continue
		out_path_base = os.path.join(roi_out_base, imp.getTitle()[:-4])
		if os.path.exists(out_path_base + '_gfp.zip') or os.path.exists(out_path_base + '_tdt.zip'):
			continue
	
		#Set Luts, Show Image and Set up RoiManager
		imp.setLuts(spanning_lut[k])
		force_lut_update(imp)
		imp.show()
		rm = RoiManager.getInstance()
		if not rm:
			rm = RoiManager()
		rm.reset()
	
		#Select GFP Roi and Save
		IJ.setTool(Toolbar.POLYGON)
		show_generic_dialog('GFP Roi Selection','Add GFP Roi(s) to RoiManager and click ok')
		rm.runCommand("Save", out_path_base + '_gfp.zip')
		rm.reset()
		imp.killRoi()
		
		#Select tdt Roi and Save
		show_generic_dialog('Tdtomato Roi Selection','Add Tdtomato Roi(s) to RoiManager and click ok')
		rm.runCommand("Save", out_path_base + '_tdt.zip')
		rm.reset()
		imp.killRoi()
	
		imp.hide()

#############################################################
####                Write PxintZ Files                   ####
#############################################################
pxintz_out_base = '/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Nrp2IFs/Analysis_Glomeruli/Extracted_Images/pxintz_outs'
for k,v in parsed_files.iteritems():
	for imp in v:
		if imp == 'blank':
			continue
		
		roi_files = [os.path.join(roi_out_base, i) for i in os.listdir(roi_out_base) if re.search('.*(?=_(gfp|tdt).zip$)', i).group() == imp.getTitle()[:-4]]
		for roi_f in roi_files:
			pxintz_out_f = os.path.join(pxintz_out_base, os.path.split(roi_f)[1][:-4] + '.pxintz')
			if os.path.exists(pxintz_out_f):
				continue
			else:
				#Load Roi File
				rm = RoiManager.getInstance()
				if not rm:
					rm = RoiManager()
				rm.reset()
				rm.runCommand("Open", roi_f)
	
				#Setup for Pxintz
				n_channels = imp.getDimensions()[2]
				channel_labs = [re.search("(?<=rgb\[255\]=).*?(?=,)", i.toString()).group(0) for i in imp.getLuts()]
				with file(pxintz_out_f,'w') as fout:
					for iter_ in range(rm.getCount()):
						if iter_ == 0: #first z written 
							fout.write('#' + str(iter_))
						else: #subsequently need line break
							fout.write('\n#' + str(iter_))
						my_roi = rm.getRoi(iter_)
	
						#Write X,Y
						X,Y = zip(*[(x.getX(),x.getY()) for x in my_roi])
						fout.write('\n>X\n')
						fout.write(','.join([str(x) for x in X]))
						fout.write('\n>Y\n')
						fout.write(','.join([str(x) for x in Y]))
	
						#Write pixel values for each channel
						for c_ in range(1, n_channels + 1):
							imp.setC(c_)
							proc = imp.getProcessor()
							px_vals = [proc.get(int(pt.getX()), int(pt.getY())) for pt in my_roi]
							fout.write('\n>' + channel_labs[c_ - 1] + '\n')
							fout.write(','.join([str(x) for x in px_vals]))

######################################################
####               Make Montage                   ####
######################################################

#Find first non-blank image and use it to create blank image
for i in parsed_files.values()[0]:
	if i=='blank':
		pass
	else:
		dims_ = i.getDimensions()
		break #once found one, stop loop

blank_image = IJ.createImage(' ',"16-bit", *dims_)
blank_image.setDisplayMode(1)

#Fill in blanks
for k,v in parsed_files.iteritems():
	for (i,elem) in enumerate(v):
		if elem == 'blank':
			parsed_files[k][i] = blank_image

#For all images, apply spanning lut, create mini-montage and save output for final montage
channels = ['0110','1000'] #panels in mini-montage
color_dict = {'gfp':Color.green, 'tdt':Color.red}
ft_size = 50

final_imps = {}
iter_order = sorted(parsed_files.keys(), key = keyfun_anim_code)
for k in iter_order:
	v = parsed_files[k]
	out = []
	for iter_, imp in enumerate(v):
		#Store the name, helps avoid "flat_" issues later
		imp_name = imp.getTitle()
		
		#Ensure mode and LUT are correct- Use AUTO LUT in this case
		imp.setDisplayMode(1)
		imp = autolevel_all_channels_v2(imp, sat = [15,0.35,0.35,0.35])
		#imp.setLuts(spanning_lut[k]) #apply spanning lut again, just in case something changed
		#force_lut_update(imp)
		
		#Create the Individual Montages
		indiv_mont = []
		for (e,i) in enumerate(channels):
			imp.setActiveChannels(i)
			imp.setOverlay(None) #clear any overlays (ie. scale bar)
			
			#Add Scale Bar if 1st Image, GFP/tdtom
			if iter_ == 0 and e == 0: 
				IJ.run(imp, "Scale Bar...", "width=50 height=8 font=28 color=White background=None location=[Lower Left] show overlay")
			
			#Add Overlays to the A647-Only Image
			elif e == 1:
				for color_ in color_dict.iterkeys():
					roi_path = os.path.join(roi_out_base, imp_name[:-4] + '_' + color_ + '.zip')
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
						overlays.setStrokeWidth(2)
						imp.setOverlay(overlays)
						imp = imp.flatten() #overwrite with flattened imp
	
			#Flatten the Imp
			indiv_mont.append(imp.flatten())
	
		#Make Individual Montage GFP/tdtom:A647
		indiv_mont_stack = ImagesToStack.run(indiv_mont)
		montage_ = MontageMaker()
		indiv_mont_out = montage_.makeMontage2(indiv_mont_stack, 1, len(channels), 1, 1, len(indiv_mont), 1, 2, False)
		indiv_mont_out.setTitle(imp_name)
	
		#Annotate the Individual Montages
		if indiv_mont_out.getTitle() != ' ':
			panel_overlay = Overlay()
			slice_lab = TextRoi(indiv_mont_out.getDimensions()[0]- 10,10, re.search("(?<=_)[0-9]+(?=.tif)",imp_name).group(), Font("SansSerif", Font.PLAIN, ft_size))
			slice_lab.setJustification(2)
			panel_overlay.add(slice_lab)
			if imp_name in first_images:
				panel_overlay.add(TextRoi(10,10, re.search("[LR](Lateral|Medial)",imp_name).group(), Font("SansSerif", Font.PLAIN, ft_size)))
		
			panel_overlay.setStrokeColor(Color.white)
			panel_overlay.setStrokeWidth(20)
		
			indiv_mont_out.setOverlay(panel_overlay)
			indiv_mont_out = indiv_mont_out.flatten()
		
		out.append(indiv_mont_out)
	
	out_stack = ImagesToStack.run(out)
	montage_ = MontageMaker()
	out = montage_.makeMontage2(out_stack, max_n_images[k], 4, 1, 1, len(out), 1, 2, False)
	out.setTitle(k)
	out.show()
	final_imps[k] = out

#Put all the Montages Together...
#x_dims, y_dims = zip(*[final_imps[i].getDimensions()[:2] for i in iter_order])

#spacer_px = 10
#final_proc = ColorProcessor(sum(x_dims) + spacer_px*(len(x_dims)-1), y_dims[0])

#for e,i in enumerate(iter_order):
#	print "Adding Pixels from " + i + '...' + str(e+1) + '/' + str(len(iter_order))
#	#Add Pixels
#	if e == 0:
#		x_offset = 0
#	else:
#		x_offset = sum(x_dims[:e]) + spacer_px*(e-1)
#	imp = final_imps[i]
#	proc = imp.getProcessor()
#	for x in range(imp.getDimensions()[0]):
#		for y in range(imp.getDimensions()[1]):
#			final_proc.putPixel(x + x_offset, y, proc.get(x,y))
#
#final_imp = ImagePlus('Summary Montage',final_proc)
#final_imp.show()

