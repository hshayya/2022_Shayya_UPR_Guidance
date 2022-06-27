from __future__ import division
import csv
import itertools
from collections import defaultdict
import os
from ij import ImagePlus, ImageStack, CompositeImage, IJ, io
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from loci.formats import ImageReader

def parse_csv_to_dict(path_, delim = ','):
	'''
	Parse csv to dictionary where each key is a column name of the dataframe and each value is a list of values for that column

	Input:
	path_ (str): file location on disk
	
	Output: dictionary with structure above. Not ordered.
	'''
	with open(path_, 'r') as fin:
		iterator_ = csv.reader(fin, delimiter = delim)
		names_ = next(iterator_)
		out = zip(*[i for i in iterator_])
		out = {k:v for (k,v) in itertools.izip(names_, out)}
		#
	
	return out

def view_dict_as_dataframe(dict_):
	'''Convert a dictionary of the form column names: list of values to a list of tab-separated strings useful for printing or writing to disk.
	'''
	out_ls = []
	out_ls.append('\t'.join(dict_.keys()))
	for i in zip(*dict_.values()):
		out_ls.append('\t'.join([str(v) for v in i]))
	#
	return out_ls

def extract_tile_from_multislice_stack(fpath, slice_, flip_, x_um, y_um, box_side_pixel, out_title, force_calibration_20x = True):
	'''Load a multi-series .nd2 image and extract a square ROI centered on point x_um, y_um. 
	
	fpath: (str) disk location of .nd2 image to load
	slice_: (int) series to extract tile from. 0-indexed
	flip_: (str) ['True' or 'False']. Should the image be rotated 180 degrees before tile extraction?
	x_um: (float) x-coordinate for glomerlus in microns (will be converted to pixel internally)
	y_um: (float) y-coordinate for glomerulus in microns 
	box_side_pixel: (int) side length for square ROI to create around glomerulus
	out_title: (str) title of output ImagePlus
	force_calibration_20x: (bool) Should a default x/y calibration of 0.325150970771092 be used regardless of metadata in .nd2 file? Helpful in cases where metadata recorded wrong!

	Return: ImagePlus contained cropped region
	'''
	#Set up Importer Options and Load Image
	opts = ImporterOptions()
	opts.setId(fpath)
	opts.setSeriesOn(slice_, True)
	load_image = BF.openImagePlus(opts)[0]

	#Flip if necessary
	if flip_.lower() == 'true':
		IJ.run(load_image, "Rotate... ", "angle=180 grid=1 interpolation=Bilinear stack")
	elif flip_.lower() != 'false':
		raise Exception('flip_ must be True or False')
	
	#Load Calibration of High-Res Image and convert X,Y from um (roi manager) -> pixel
	cal_table = load_image.getCalibration()
	px_x = int(x_um/cal_table.pixelWidth) #Relies on the fact that X,Y in um, calibration properly calculated between downsampled image and original
	px_y = int(y_um/cal_table.pixelHeight)

	if force_calibration_20x:
		px_x = int(x_um/0.325150970771092)
		px_y = int(y_um/0.325150970771092)

	#Draw Roi and Crop
	final_x = px_x - int(box_side_pixel/2)
	final_y = px_y - int(box_side_pixel/2)
	
	load_image.setRoi(final_x,final_y, box_side_pixel, box_side_pixel)
	out = load_image.crop('stack')
	out.setTitle(out_title)
	
	return out

#Script
csv_base_path = '/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Nrp2IFs/Analysis_Glomeruli/Rois_To_Extract'
f_ls = [os.path.join(csv_base_path, i) for i in os.listdir(csv_base_path) if i[0] != '.']
out_path = '/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Nrp2IFs/Analysis_Glomeruli/Extracted_Images'

for f in f_ls:
	print 'Working on Animal: ' + os.path.split(f)[1][:-4]
	run_dict = parse_csv_to_dict(f)

	for i in range(len(run_dict['output_file'])):
		print 'Iteration ' + str(i+1) + '/' + str(len(run_dict['output_file']))
		out_title = '_'.join([run_dict['output_file'][i], run_dict['Glomerulus'][i], str(run_dict['Slice'][i])])
		if os.path.exists(os.path.join(out_path, out_title + '.tif')):
			print 'Found Output: ' + out_title + '. Skipping...'
			continue
		output_image = extract_tile_from_multislice_stack(fpath = run_dict['image_file'][i], slice_ = int(run_dict['slice_order'][i]), flip_ = run_dict['flip_bool'][i], x_um = float(run_dict['X'][i]), y_um=float(run_dict['Y'][i]), box_side_pixel = 500, out_title = out_title, force_calibration_20x = True)
		out_file = os.path.join(out_path, out_title + '.tif')
		io.FileSaver(output_image).saveAsTiff(out_file)
		print '...done!'

#NB: there was one file (5-19-22_4wkM28_Ddit3_Exp1_Set2.nd2) where calibration is recorded incorrectly! This was causing huge issues, so I introduced the force_calibration_20x to fix it! 
#Probaby the issue happened w/ the Set 2 image b/c that one had to be "finished" mid-imaging b/c my reservation time ran out! All other images are perfect, the calibration at 20x on the Spinning disk is fixed regardless of file size etc. so the code here should be safe anyways!
