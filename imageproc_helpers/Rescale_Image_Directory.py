from __future__ import division
from ij import io
from ij import plugin
from ij import IJ
from ij import gui
from loci.plugins import BF
from ij import ImagePlus, ImageStack, CompositeImage
from ij.process import ImageProcessor
from glob import glob
import re
import os

#Rescale folder of tiled .nd2 images -> downsampled .tif images of specific pixel dimensions
#Useful for downsampling tiled images of the whole OE to speed up loading etc. when counting cells.

def rescale_compositeImage(imp, dims_resize):
	''' 
	Rescale a CompositeImage to specified pixel dimensions recomputing scale and retaining colorization 
	
	Input:
	imp: CompositeImage to be rescaled
	dims_resize: int of len = 2 giving pixel dimensions to resize to
	
	Output: rescaled CompositeImage
	'''
	
	#Get calibration from CompositeImage 
	cal_table = imp.getCalibration()
	micron_width = imp.getDimensions()[0] * cal_table.pixelWidth
	micron_height = imp.getDimensions()[1] * cal_table.pixelHeight
	
	#Get Luts from the CompositeImage
	luts = imp.getLuts()
	
	#Convert to ImageStack and Interate through ImageProcessors to Resize
	out = imp.getImageStack()
	
	out_stack = ImageStack()
	for i in range(1, imp.getNChannels()+1):
		out_stack.addSlice(out.getProcessor(i).resize(*dims_resize))
		
	#Convert Output -> ImagePlus -> CompositeImage and colorize
	new_title = imp.getTitle() + '-Downsampled'
	out_image = CompositeImage(ImagePlus(new_title,out_stack))
	out_image.setLuts(luts)
	out_image.show()
	
	#Set up new calibration
	cal_table.pixelWidth = micron_width/dims_resize[0]
	cal_table.pixelHeight = micron_height/dims_resize[1]
	out_image.setCalibration(cal_table)
	
	return out_image

#Set these params. 
dir_ = '/path/to/input/folder/*.nd2'
out_dir = '/path/to/output/folder'
dims = [5000,5000] #px x, px y

for (i,file_) in enumerate(glob(dir_)):
	print 'Working on' + file_ + '...' + str(i+1) + '/' + str(len(glob(dir_)))

	#Run the Script
	load_image = BF.openImagePlus(file_)[0]
	out_image = rescale_compositeImage(load_image, dims)

	#Figure out where to Save
	file_name = re.split('/', file_)[-1]
	file_name = re.search('.*(?=.nd2)', file_name).group(0)
	out_path = os.path.join(out_dir, file_name + '-Downsampled.tif')

	#Save 
	print 'Saving file to ' + out_path
	io.FileSaver(out_image).saveAsTiff(out_path)
	IJ.run('Close All')
	print '...done!'
