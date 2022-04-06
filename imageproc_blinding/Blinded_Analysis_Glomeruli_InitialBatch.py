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

#Framework used in blinded anlaysis of OB glomerular configurations
#Essentially gives random names to all the OB glomerular images, presents MaxProjections & Z-stacks to experimenter -> records blinded call intermixed/compartmentalized...
#Analysis was done as the images were collected, so code here covers only the first batch of glomeruli to analyze of a given OR/Flox Allele
#For subsequent batches, see _AdditionalBatch code. 
#NB: analysis was always performed in batches with approx identical #'s of WT/Het/cKO images!

#User Inputs
dir_ = '/path/to/root/dir/with/all/images'
regex = '^MAX.*(WT|Ctrl|Exp).*20x.*.tif$'
ft_size = 50

#Restrict search to top level subdirectories
subdirs = [os.path.join(dir_, d) for d in next(os.walk(dir_))[1]]
f_list = [os.path.join(subdir, f) for subdir in subdirs for f in os.listdir(subdir) if re.match(regex, f) is not None and re.search('4wk',f) is None]

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

#Write Glomeruli in order with a random code
glom_labels = ['LLateral','LMedial','RLateral','RMedial']*int((len(fout_final)/4))
animal_labels = [f for i in animal_order for f in [i]*4]
rep_len = int(math.ceil(len(fout_final)/26))
alphanum = [f*(i+1) for i in range(rep_len) for f in list(string.ascii_lowercase)][:len(fout_final)]
random.Random(4).shuffle(alphanum)
print alphanum

#Uncomment this to write the lookup table to disk (required for subsequent steps of code)
#with open('/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Experiments/Analysis_Glomeruli/Blinded_Assignment/LookupTable_4-2-21.tsv','w') as fout:
#	fout.write('\t'.join(['animal','glomerulus','file','random_code']) + '\n')
#	for z in itertools.izip(animal_labels, glom_labels, fout_final, alphanum):
#		fout.write('\t'.join(list(z)) + '\n')
#	#

#Read in the Lookup Table Line by Line
reader = csv.DictReader(open('/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Experiments/Analysis_Glomeruli/Blinded_Assignment/LookupTable_4-2-21.tsv','r'),delimiter = '\t') 
parsed_entries = [i for i in reader]
parsed_entries.sort(key = lambda x: x['random_code'])
dict_order = parsed_entries[1].keys()

with open('/Volumes/Image_Data/M28iCreiGFP_Ddit3cKO_Experiments/Analysis_Glomeruli/Blinded_Assignment/LookupTable_4-2-21_Output.tsv','w') as fout:
	fout.write('\t'.join(dict_order) + '\t' + 'config' + '\n')
	for e,entry in enumerate(parsed_entries):
		print "Working on line " + str(e+1) + '/' + str(len(parsed_entries))
		if entry['file'] == 'blank':
			continue
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
			full_path = os.path.join(rt, re.search('(?<=MAX_).*(?=.tif)', fname).group(0) + '.nd2') #see MaxProj code. Mapping MaxProj <-> Z stack always works this way.
			
			opts = ImporterOptions()
			opts.setId(full_path)
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
			fout.write('\t'.join([entry[i] for i in dict_order]) + '\t' + out + '\n')

			IJ.run("Close All", "");
			
		
		