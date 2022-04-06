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
from ij.gui import NonBlockingGenericDialog

#Framework used in blinded anlaysis of OB glomerular configurations AFTER running _InitialBatch code for a given OR/Flox Allele combo
#Same idea as _InitialBatch, but only analyzes NEW images that were not already analyzed! This makes it so do not have to analyze ALL images EVERY TIME want to make the figure!
#Gives random names to all the OB glomerular images, presents MaxProjections & Z-stacks to experimenter -> records blinded call intermixed/compartmentalized...
#NB: analysis was always performed in batches with approx identical #'s of WT/Het/cKO images!

#User Inputs
dir_ = '/path/to/root/dir/where/all/images/stored'
regex = '^MAX.*(WT|Ctrl|Exp).*20x.*.tif$'
ft_size = 50

last_csv = '/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Experiments/Analysis_Glomeruli/Blinded_Assignment/LookupTable_05_13_21_Output.tsv' #Output file from the last batch that was run (or initial batch)
date_append = '2-14-22' #Today's date. Used to make the blinded annotations unique btwn batches.

#Restrict search to top level subdirectories
subdirs = [os.path.join(dir_, d) for d in next(os.walk(dir_))[1]]
f_list = [os.path.join(subdir, f) for subdir in subdirs for f in os.listdir(subdir) if re.match(regex, f) is not None and re.search('4wk|tdtomonly',f) is None]

#Group files by animal
files_dict = defaultdict(list)
for f in f_list:
	animal = re.search('(WT|Ctrl|Exp)[0-9]+', os.path.split(f)[1]).group()
	files_dict[animal].append(f)

#Order the files within each animal, inserting blanks in the correct places
glom_order = ['LLateral','L(Medial|me)','RLateral','R(Medial|me)']

files_out_dict = {k:['blank']*4 for k in files_dict.iterkeys()}
for k,ls in files_dict.iteritems():
	for f in ls:
		for i,regex in enumerate(glom_order):
			if re.search(regex, os.path.split(f)[1]) is not None:
				files_out_dict[k][i] = f

#Need same number of animals across WT/Ctrl if doing two stacks. 
nums = set([re.search('[0-9]+', i).group() for i in files_out_dict.iterkeys()])
nums = sorted([int(i) for i in list(nums)])

#Grab animals in order and append (already ordered) images (with blanks) to final list. Add four blanks if that animal is missing
animal_order = [g+str(i) for i in nums for g in ['WT','Ctrl','Exp']]
fout_final = []
for animal in animal_order:
	try:
		ls = files_out_dict[animal]
		for f in ls:
			fout_final.append(f)
	except KeyError:
		for f in ['blank']*4:
			fout_final.append(f)

#Get the Animals/Labels for All Files
glom_labels = ['LLateral','LMedial','RLateral','RMedial']*int((len(fout_final)/4))
animal_labels = [f for i in animal_order for f in [i]*4]

#Load the Last OUTPUT CSV File
reader = csv.DictReader(open(last_csv,'r'),delimiter = '\t') 
parsed_entries = [i for i in reader]
parsed_entries_dict = {i['animal']+'_'+i['glomerulus']:i for i in parsed_entries}

#Need to get random codes for new images
new_images = [i for i in fout_final if i not in (f['file'] for f in parsed_entries) and i != 'blank'] #blank can lead to weird things here

#Helpful to run up to here to just ensure that it found all the new images you expect in your batch...
for i in new_images:
	print os.path.split(i)[1]

rep_len = int(math.ceil(len(new_images)/26))
alphanum = [f*(i+1) + '_' + date_append for i in range(rep_len) for f in list(string.ascii_lowercase)][:len(new_images)] #in order
random.Random(4).shuffle(alphanum)

#Now iterate through everything in fout_final that isn't blank and either add the existing mappings or create new ones
new_input = os.path.join(os.path.split(last_csv)[0], 'LookupTable_' + date_append + '.tsv')
headers = ['animal','glomerulus','file','random_code']

#Uncomment this to write the new input file including the new batch
#with open(new_input,'w') as fout:
#	fout.write('\t'.join(headers) + '\n')
#	iter_ = 0
#	for (e, z) in enumerate(itertools.izip(animal_labels, glom_labels, fout_final)):
#		if z[2] == 'blank':
#			continue #only images that have an image!
#		elif z[2] in new_images:
#			fout.write('\t'.join(list(z) + [alphanum[iter_]]) + '\n')
#			iter_+=1
#		else:
#			fout.write('\t'.join([parsed_entries_dict[z[0]+'_'+z[1]][k] for k in headers]) + '\n')
#		#
#	#

#Read in the New Lookup Table Line by Line 
new_reader = csv.DictReader(open(new_input,'r'),delimiter = '\t') 
new_parsed_entries = [i for i in new_reader]
new_parsed_entries.sort(key = lambda x: x['random_code'])

#Go line by line. If not new -> copy old result. If new, analyze.
out_final_path = os.path.join(os.path.split(new_input)[0], os.path.split(new_input)[1][:-4] + '_Output.tsv')
with open(out_final_path,'w') as fout:
	fout.write('\t'.join(headers) + '\t' + 'config' + '\n')
	for e,entry in enumerate(new_parsed_entries):
		print "Working on line " + str(e+1) + '/' + str(len(new_parsed_entries))
		if entry['file'] not in new_images:
			fout.write('\t'.join([parsed_entries_dict[entry['animal'] + '_' + entry['glomerulus']][x] for x in headers + ['config']]) + '\n')
		else: 
			#Load Max Projection
			imp = CompositeImage(ImagePlus(entry['file']))
			imp.setTitle(entry['random_code'])
			for c in range(imp.getDimensions()[2]):
				imp.setC(c+1) #1-based...
				IJ.run(imp, "Enhance Contrast", "saturated=0.5")
			
			imp.setDisplayMode(1)
			imp.show()
			
			#Load Nd2 File
			rt, fname = os.path.split(entry['file'])
			full_path = os.path.join(rt, re.search('(?<=MAX_).*(?=.tif)', fname).group(0) + '.nd2')
			symlink_path = os.path.join(os.path.split(last_csv)[0], entry['random_code'] + '.nd2') #create a symlink with correct filename (otherwise will read that into FIJI and affect blinding!)
			os.symlink(full_path, symlink_path)
			
			opts = ImporterOptions()
			opts.setId(symlink_path)
			opts.setVirtual(True) 
	
			imp = BF.openImagePlus(opts)[0]
			imp.setTitle(entry['random_code'] + '_full')
			imp.setZ(int(imp.getDimensions()[3]/2))
			for c in range(imp.getDimensions()[2]):
				imp.setC(c+1) #1-based...
				IJ.run(imp, "Enhance Contrast", "saturated=0.5")
			
			imp.setDisplayMode(1)
			imp.show()
	
			#Show Dialog
			gd = NonBlockingGenericDialog('Image Coding')
			choices = ['none','intermixed','compartmentalized','adjacent/accessory','no_coalescence']
			gd.addChoice('Code',choices,'none')
			gd.showDialog()
	
			out = choices[gd.getNextChoiceIndex()]
			
			#Write out
			fout.write('\t'.join([entry[i] for i in headers]) + '\t' + out + '\n')

			IJ.run("Close All", "");
			os.remove(symlink_path) #remove the symlink
			
		
		