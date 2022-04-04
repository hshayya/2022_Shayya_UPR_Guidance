from __future__ import division
from ij import io
from ij import plugin
from ij import IJ
from ij import gui
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from loci.formats import ImageReader
import csv
import re
from glob import glob
import os

#Convert folder of .nd2 Z stacks -> .tif max projections. 

#Worker Function 
def run_and_save_zproj(fpath, show_max = True, save_output = False, use_virtual = False):
	'''Load a file from a given image path. If there are multiple Z timepoints, run a max projection and display or save output. Otherwise do nothing. 

	fpath: (str) path to file of interest
	show_max: (bool) Should the max projection be displayed in FIJI?
	save_output: (bool) Should the max projection be saved to disk? 

	Return: None
	If show_max, image is shown. If save_output, a new file named ('MAX_' + filename) is saved in the same folder as input.
	'''
	#Load the Image
	opts = ImporterOptions()
	opts.setId(fpath)
	
	if use_virtual:
		opts.setVirtual(True) 
	
	imps = BF.openImagePlus(opts)[0] #creates array
	
	#Check the Z dimension and Max project/save 
	if imps.getDimensions()[3] == 1:
		print 'Stack z dimension = 1, Program Complete'
		return
	else: 
		imp_max = plugin.ZProjector().run(imps, 'max')
		if show_max:
			imp_max.show()
		if save_output:
			file_split = os.path.split(fpath)
			newname = 'MAX_' + re.search('(.*)(?=\\.[A-Za-z0-9]+$)', file_split[-1]).group(0) + '.tif'
			fout_path = os.path.join(file_split[0], newname)
			io.FileSaver(imp_max).saveAsTiff(fout_path)
		print 'Ran Z projection as directed... done!'
		return

#Run the Script
my_glob = '/path/to/folder/*nd2'
ram_limit = 4444*10**6


for (i,file_) in enumerate(glob(my_glob)):
	print 'Working on' + file_ + '...' + str(i+1) + '/' + str(len(glob(my_glob)))
	
	#Compute Output Name & Check if this Exists
	file_split = os.path.split(file_)
	newname = 'MAX_' + re.search('(.*)(?=\\.[A-Za-z0-9]+$)', file_split[-1]).group(0) + '.tif'
	fout_path = os.path.join(file_split[0], newname)
	
	if os.path.exists(fout_path):
		print 'Output at ' + fout_path + ' already exists, skipped...'
		continue

	#Run the Analysis
	fsize = os.stat(file_).st_size
	if fsize >= ram_limit:
		print 'File Size of ' + str(fsize/10**9) + 'GB exceeds RAM Limit of ' + str(ram_limit/10**9) + 'GB. \nUsing Virtual Stack...'
		run_and_save_zproj(file_, show_max = False, save_output = True, use_virtual = True)
	else:
		run_and_save_zproj(file_, show_max = False, save_output = True)
	#
