from __future__ import division
from ij import ImagePlus, ImageStack, CompositeImage, IJ, io
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from loci.formats import ImageReader, MetadataTools
import itertools
import re
from collections import OrderedDict, defaultdict
import os

def rescale_compositeImage(imp, dims_resize):
	''' 
	Rescale a CompositeImage to specified pixel dimensions recomputing scale and retaining colorization 
	
	Input:
	imp: CompositeImage to be rescaled
	dims_resize: list of ints, len = 2 giving pixel dimensions to resize to
	
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

#Set up the organization for this dataset in Python
datasets = {'a':'5-19-22_4wkM28_Ddit3_Exp1_Set1.nd2', 'b':'5-19-22_4wkM28_Ddit3_Exp1_Set2.nd2','c':'5-19-22_4wkM28_Ddit3_Exp1_Set3.nd2'}

#(anim,slicestart,sliceend(+1), flip_bool)
dset_mapping = {'a':[('exp',1,8,True)],'b':[('exp',1,16,True)], 'c':[('exp',1,12,True)]}
output = {'exp':'7-30-21_4250_4wkMor28_Exp1'}

#Create Lists with Information for All Slices, as Imaged
slice_names = list() 
image_file = list()
output_file = list()
flip_bool = list()
for (file_, ls_) in dset_mapping.iteritems():
	for tup in ls_:
		lab, start, stop, flip = tup
		for slice_ in range(start,stop):
			slice_names.append('Slice' + str(slice_))
			image_file.append(datasets[file_])
			output_file.append(output[lab])
			flip_bool.append(flip)

slice_order = [i for i in itertools.chain.from_iterable([range(7), range(15), range(11)])] #order of slices in the files
out_order = [str(index).zfill(4) for index in itertools.chain.from_iterable([range(33)])] #order to be named with

#Check that everything is the proper length and save that length
lens_ = set([len(i) for i in [slice_names, image_file, output_file, slice_order, out_order, flip_bool]])
if len(lens_) == 1:
	lens_ = list(lens_)[0]
	print 'Length Checks Passed...'
else:
	raise Exception('Computed Lists of Unequal Length')

#Check everything is correctly ordered and write to file
#for ls in itertools.izip(slice_names, image_file, output_file, slice_order, out_order, flip_bool):
#	print '\t'.join([str(i) for i in ls])

#with file('/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Nrp2IFs/5-19-22_4wkOBSect_Exp1/Slice_Lookup_Table.tsv', 'w') as fout:
#	fout.write('slice_names\timage_file\toutput_file\tslice_order\tout_order\tflip_bool\n')
#	for ls in itertools.izip(slice_names, image_file, output_file, slice_order, out_order, flip_bool):
#		fout.write('\t'.join([str(i) for i in ls]) + '\n')
#	#

#By image, set everything up
run_ = defaultdict(list)
for ls in itertools.izip(slice_names, image_file, output_file, slice_order, out_order, flip_bool):
	run_[ls[1]].append([i for (elem, i) in enumerate(ls) if elem != 1])

#print run_.keys()
#for i in run_.items()[0]:
#	print(i)

#Run the Analysis
ds_frac = 0.5
iter_ = 1
base_path = '/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Nrp2IFs/5-19-22_4wkOBSect_Exp1'
for (k,v) in run_.iteritems():
	fpath = os.path.join(base_path, k)
	
	#Parse Metadata & Set Scaler Options (Latter uses only the first series, assume all same size...)
	print('Parsing Metadata...' + fpath)
	reader = ImageReader()
	reader.setId(fpath)
	
	series_ct = reader.getSeriesCount()
	print(str(series_ct) + ' series')
	
	orig_x = reader.getSizeX()
	orig_y = reader.getSizeY()
	print('X: ' + str(orig_x) + '\n' + 'Y: ' + str(orig_y))
	out_px = [int(round(i*ds_frac)) for i in [orig_x, orig_y]]


	#Set up BF Importer Conserved Options
	opts = ImporterOptions()
	#opts.setVirtual(True) #Way slow- do NOT use if reading one slice at a time like this.
	opts.setId(fpath)
	
	#Iterate through the Slices
	for slice_ in v:
		#Naming Stuff
		name_, fout_, idx_, name_num_, flip_ = slice_
		out_path = os.path.join(base_path, fout_, name_num_ + '_' + name_ + '.tif')
		print('Iteration ' + str(iter_) + '/' + str(lens_) + '...')
		print('File: ' + fpath + '\nOutput:' + out_path + '\nFlip:' + str(flip_))

		#Load the Right Slice and Rescale
		opts.setSeriesOn(idx_, True)
		load_image = BF.openImagePlus(opts)[0]
		out_image = rescale_compositeImage(load_image, out_px)	
		out_image.setTitle(name_num_ + '_' + name_)

		#Flip if Needed
		if flip_:
			#IJ.run(out_image, "Flip Vertically", "stack") #Should be "ROTATE 180 Degrees" NOT flip... this mixes up L/R
			IJ.run(out_image, "Rotate... ", "angle=180 grid=1 interpolation=Bilinear stack")
		#
		
		#Save and Close
		io.FileSaver(out_image).saveAsTiff(out_path)
		out_image.changes = False
		out_image.close()
		opts.clearSeries()

		#Adjust iter
		iter_+=1
	#
