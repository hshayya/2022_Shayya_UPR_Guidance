import xml.etree.ElementTree as ET
import csv
import os
import re
from ij import IJ
from loci.plugins.in import ImporterOptions
from loci.plugins import BF
from ij.plugin import ImagesToStack
from ij import io

#Takes a set of CellCounter.xml files (generated from downsampled .tif image) & associated .nd2 images of OE sections and extracts the counted cells -> tiny images for downstream analysis
#Counter.xml files are in downsampled space (see downsampling script), so conversion is made -> full-res coordinates before tiles surrounding cells are extracted

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

def extract_cells(input_dict, rotate_angle, px_size_cell_box, types_of_interest, dir_out_base):
	'''
	---Input---
	(input_dict): (tuple) of form section_name, dict{'fullres':/path/to/fullres, 'counter':{dict_gen_by parse_cellcounter_to_dict}}. 
	(rotate_angle): (str!!) angle to rotate section (270 = rotate left, 90 = rotate right)
	(px_size_cell_box): (int) must be EVEN number, size in pixels of one side of roi to extract
	(types_of_interest): (dict) of form keys = strs of cell types in input_dict[1]['counter'], vals = labels
	(dir_out_base): base output directory

	--Output--
	If /dir_out_base/section_name unless is already present, does nothing. 
	Otherwise creates directory, extracts cells to /dir_out_base/section_name/section_name{_label}{_n}.tif
	No return value. 
	'''
	anim, vals = input_dict
	
	#Make Output Folder
	dir_out_section = os.path.join(dir_out_base, anim)
	try:
		os.makedirs(dir_out_section)
	except OSError:
		print 'Section ' + anim + ' Already Analyzed'
		return
	
	#Load the Fullres image
	opts = ImporterOptions()
	opts.setId(vals['fullres'])
	full_image = BF.openImagePlus(opts)[0]
	
	#Rotate
	IJ.run(full_image, "Rotate... ", "angle=" + rotate_angle + " grid=1 interpolation=Bilinear stack")
	
	#Extract Cal Table
	cal_table = full_image.getCalibration()
	
	#Loop through Cells and Extract Rois
	for cell_type, cell_label in types_of_interest.iteritems():
		print 'Working on cell_type ' + cell_label
		for i in range(len(vals['counter'][cell_type])):
			print 'Iteration ' + str(i+1) + '/' + str(len(vals['counter'][cell_type]))
			
			#Convert Px Downsampled -> Px Full Res
			x_full_px = cal_table.getRawX(vals['counter'][cell_type][i][0] * vals['counter']['x_cal']) #in px
			y_full_px = cal_table.getRawY(vals['counter'][cell_type][i][1] * vals['counter']['y_cal']) #in px
	
			#Set Roi and Crop
			full_image.setRoi(int(x_full_px - px_size_cell_box/2), int(y_full_px - px_size_cell_box/2), int(px_size_cell_box), int(px_size_cell_box))
			out_image = full_image.crop('stack')
	
			#Save Output
			out_title = '_'.join([anim, cell_label, str(i)])
			out_image.setTitle(out_title)
			io.FileSaver(out_image).saveAsTiff(os.path.join(dir_out_section, out_title + '.tif'))
			print '...done!'

	full_image.close()
	return

#Load Xml Files
xml_locs = ['/path/to/cellcounter/xmls']
xml_files = [os.path.join(base_, f) for base_ in xml_locs for f in os.listdir(base_) if f[-3:] == 'xml' and f[0] != '.']

#Work through each xml file
for e,xml_ in enumerate(xml_files):
	print 'Working on file: ' + os.path.split(xml_)[1] + '...' + str(e+1) + '/' + str(len(xml_files))
	
	#Find the orig .nd2 file
	orig_f_name = re.search('(?<=CellCounter_).*(?=\\-Downsampled)', os.path.split(xml_)[1]).group() + '.nd2'
	search_dir = '/'.join(os.path.split(xml_)[0].split('/')[:-1]) #.nd2 files are in ../subfolder relative to xmls. Go up one level
	files_found = [os.path.join(root, f) for (root, dirs, files) in os.walk(search_dir) for f in files if f == orig_f_name] #then recursive search down for the right .nd2 file

	#Just make sure you found the right file...works!
	if len(files_found) == 1:
		fullres_image = files_found[0]
	else:
		print "Could not find fullres image."
		raise ValueError('Found 0 or >1 matching file')

	#Extract Cells
	input_item = (re.search('(?<=_).*',orig_f_name[:-4]).group(), {'fullres':fullres_image, 'counter':parse_cellcounter_to_dict(xml_)})
	extract_cells(input_dict=input_item, rotate_angle='270', px_size_cell_box=60, types_of_interest={'7':'tdtom','8':'gfp'}, dir_out_base='/path/to/write/extracted/images')

#Makes a new folder in the dir_out_base for each section analyzed. Images are of the form SlideName_{gfp|tdtom}_cell#.tif