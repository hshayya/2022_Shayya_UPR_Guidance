import xml.etree.ElementTree as ET
import csv
import os
import re
from ij import IJ
from loci.plugins.in import ImporterOptions
from loci.plugins import BF
from ij.plugin import ImagesToStack
from ij import io

#Records metadata (x,y location) for cells that were extracted with 1_find_extract_cells.py
#metadata will be used in subsequent analysis to cluster cells from similar locations on the section -> semi-quantiative, local, analysis

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

#Load Xml Files
xml_locs = ['/path/to/xml/files'] #same as used in find_extract_cells
xml_files = [os.path.join(base_, f) for base_ in xml_locs for f in os.listdir(base_) if f[-3:] == 'xml' and f[0] != '.']

#Work through each xml file
f_out_path = '/path/to/annotation/out.tsv' 

with open(f_out_path,'w') as fout:
	fout.write('\t'.join(['cell','x_um','y_um']))
	for e,xml_ in enumerate(xml_files):
		print 'Working on file: ' + os.path.split(xml_)[1] + '...' + str(e+1) + '/' + str(len(xml_files))
		
		#Find the orig .nd2 file, copied from find_extract_cells.py, see that code for more details.
		orig_f_name = re.search('(?<=CellCounter_).*(?=\\-Downsampled)', os.path.split(xml_)[1]).group() + '.nd2'
		search_dir = '/'.join(os.path.split(xml_)[0].split('/')[:-1])
		files_found = [os.path.join(root, f) for (root, dirs, files) in os.walk(search_dir) for f in files if f == orig_f_name]
		
		if len(files_found) == 1:
			fullres_image = files_found[0]
		else:
			print "Could not find fullres image."
			raise ValueError('Found 0 or >1 matching file')
	
		#Generate the original inputs that were passed to extract_cells
		input_item = (re.search('(?<=_).*',orig_f_name[:-4]).group(), {'fullres':fullres_image, 'counter':parse_cellcounter_to_dict(xml_)})
		input_dict = input_item
		types_of_interest={'7':'tdtom','8':'gfp'}
		
		#Copied from the "Extract Cells", recovering positional info and writing to disk instead of extracting cell -> small image.
		anim, vals = input_dict
		#Loop through Cells and Annotate.
		for cell_type, cell_label in types_of_interest.iteritems():
			print 'Working on cell_type ' + cell_label
			for i in range(len(vals['counter'][cell_type])):
				print 'Iteration ' + str(i+1) + '/' + str(len(vals['counter'][cell_type]))
				
				#Convert Px Downsampled -> Px Full Res
				x_full_px = vals['counter'][cell_type][i][0] * vals['counter']['x_cal'] #in um
				y_full_px = vals['counter'][cell_type][i][1] * vals['counter']['y_cal'] #in um
		
				#Write Information
				out_title = '_'.join([anim, cell_label, str(i)])
				fout.write('\n' + '\t'.join([out_title, str(x_full_px), str(y_full_px)]))

#Final tsv of form cell_label,x,y.