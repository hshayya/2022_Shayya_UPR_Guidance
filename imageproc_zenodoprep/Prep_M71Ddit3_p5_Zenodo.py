from __future__ import division
import shutil
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
from ij.gui import NonBlockingGenericDialog, Roi, Overlay, TextRoi
from java.awt import Color, Font

def autolevel_all_channels(imp, sat = 0.35):
	imp.setDisplayMode(1)
	sat_str = "saturated=" + str(sat)
	for c in range(imp.getDimensions()[2]):
		imp.setC(c+1) #1-based...
		IJ.run(imp, "Enhance Contrast", sat_str)
	return imp

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

#User Inputs
dir_ = '/Volumes/Image_Data/M71iCreiGFP_Ddit3cKO_Experiments' 
last_csv = '/Volumes/Image_Data/M71iCreiGFP_Ddit3cKO_Experiments/Analysis_Glomeruli/Blinded_Assignments/LookupTable_8_20_21_Output.tsv'
regex = '^MAX.*(WT|Ctrl|Exp).*20x.*.tif$'
file_copy_dir = '/Volumes/Backup_HJS/Papers/UPR_AxonGuidance_2021/Zenodo/Images_20xGlomeruli/p5_M71_Ddit3'
fix_mappings = {'WT':'WT','Ctrl':'Het','Exp':'cKO'} #on my filesystem WT = WT, Ctrl = Het, Exp = cKO. Update names so images posted to Zenodo use same names as we do in paper.

#Restrict search to top level subdirectories
subdirs = [os.path.join(dir_, d) for d in next(os.walk(dir_))[1]]
f_list = [os.path.join(subdir, f) for subdir in subdirs for f in os.listdir(subdir) if re.match(regex, f) is not None and re.search('4wk|tdtomonly',f) is None]

#Group files by animal
files_dict = defaultdict(list)
for f in f_list:
	animal = re.search('(WT|Ctrl|Exp)[0-9]+', os.path.split(f)[1]).group()
	files_dict[animal].append(f)

#Order the files within each animal, inserting blanks in the correct places
glom_order = ['LLateral','RLateral']

files_out_dict = {k:['blank']*2 for k in files_dict.iterkeys()}
for k,ls in files_dict.iteritems():
	for f in ls:
		for i,regex in enumerate(glom_order):
			if re.search(regex, os.path.split(f)[1]) is not None:
				files_out_dict[k][i] = f

#Load annations to color images!
reader = csv.DictReader(open(last_csv,'r'),delimiter = '\t') #Load the Last OUTPUT CSV File
parsed_entries = [i for i in reader] 
parsed_entries_dict = {i['animal']+'_'+i['glomerulus']:i for i in parsed_entries}

color_mapping = {"compartmentalized":Color.blue, "adjacent/accessory":Color.red, "intermixed":Color.white,'no_coalescence':Color.green}
border_width = 20

for gt in ['WT','Ctrl','Exp']:
	out_ls = []
	configs = []
	top_lab = []
	bot_lab = []
	#Loop through animals in order
	anims = sorted([a for a in files_out_dict.iterkeys() if re.search(gt,a) is not None], key = lambda x: int(re.search('[0-9]+',x).group()))
	for anim_ in anims:
		#Correct anim_ for labels (WT = WT, Ctrl = Het, Exp = cKO)
		anim_for_lab = fix_mappings[gt] + re.search('[0-9]+',anim_).group()
		
		#Loop through gloms for that animal, in order
		for lab,fpath_ in zip(['LLateral','RLateral'], files_out_dict[anim_]):
			#Default is to have a blank image
			imp = 'blank'
			config_write = 'intermixed'
			top_lab_write = ''
			bot_lab_write = ''
			#Find Configuration
			try:
				config = parsed_entries_dict[anim_ + '_' + lab]['config']
				if config == 'none': 
					pass #just leave blank, do not include in montage
				else:
					#only case where not a blank image, overwrite defaults
					imp = parse_imp(fpath_, rotate = False, autolevel = True)
					imp = imp.flatten()
					config_write = config
					top_lab_write = anim_for_lab 
					bot_lab_write = lab

					#Write Files that were analyzed
					#shutil.copy(fpath_, os.path.join(file_copy_dir, anim_for_lab + '_' + lab + '.tif'))
			except KeyError: 
				print "No annotation found for " + anim_ + '_' + lab + "...skipped!"
				pass #if blank, do not include in montage
			
			#Add to lists
			out_ls.append(imp)
			configs.append(config_write)
			top_lab.append(top_lab_write)
			bot_lab.append(bot_lab_write)
	
	#Find the first non "blank" image and use it to create the blank image
	for i in out_ls:
		if i=='blank':
			pass
		else:
			dims_ = [i.getDimensions()[0], i.getDimensions()[1], 1, i.getBitDepth()]
			break #once found one, stop loop
	
	blank_image = IJ.createImage(' ', *dims_)
	blank_image = blank_image.flatten()
	blank_image.setTitle(' ')
	
	#Replace "blanks" with blank image
	for (i,elem) in enumerate(out_ls):
		if elem == 'blank':
			out_ls[i] = blank_image

	#Make Stack (will adjust sizes to all the same size!)
	out_show = ImagesToStack.run(out_ls)

	#Add the labels/annotation
	for slice_,c,t,b in zip(range(out_show.getDimensions()[3]), configs, top_lab, bot_lab):
		#Grab Resized Images
		out_show.setSlice(slice_+1)
		imp = out_show.crop() #just pulls that slice
		#Add Roi
		imp.setOverlay(Roi(0,0,imp.getDimensions()[0]-1, imp.getDimensions()[1]-1), color_mapping[c], border_width, None)
		imp = imp.flatten()
		#Add Labs
		ovrly = Overlay()
		ovrly.add(TextRoi(50,50, top_lab[slice_], Font("SansSerif", Font.PLAIN, 150)))
		ovrly.add(TextRoi(50,imp.getDimensions()[1]-200, bot_lab[slice_], Font("SansSerif", Font.PLAIN, 150)))
		ovrly.setStrokeColor(Color.white)
		ovrly.setStrokeWidth(20)
		imp.setOverlay(ovrly)
		out_ls[slice_] = imp.flatten()

	#Re-make Stack
	out_show = ImagesToStack.run(out_ls)
	
	#Make Montage
	montage_ = MontageMaker()
	out = montage_.makeMontage2(out_show, 2, int(len(out_ls)/2), 1, 1, len(out_ls), 1, 0, False)
	out.setTitle(gt)
	out.show()