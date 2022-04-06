import csv
from ij import ImagePlus, CompositeImage, IJ, gui
from ij.plugin import ImagesToStack

#Prepare Stack of 20x Images for a given OR WT/Ctrl/cKO
#Auto-levels each panel (biology of interest here is overlap & correlation red/green, not absolute levels).

reader = csv.DictReader(open('/path/to/blinded/annotation/out','r'), delimiter = '\t') #used blinded annotation output tsv to select images
#ensured that fractions of intermixed/compartmentalized etc. on final montage ~= the observed frequencies in the blinded annotations for that OR/gt combo.

#Parse the dictionary
reader = [i for i in reader]

slides_of_interest = ['b','k','l','hh','o','e_05_13_21','jj','z','i_05_13_21']
random_codes = []
imps = []

for elem in reader:
	if elem['random_code'] in slides_of_interest:
		random_codes.append(elem['random_code'])
		
		imp = CompositeImage(ImagePlus(elem['file']))
		#Stretch Histogram for Each Channel
		for c in range(imp.getDimensions()[2]):
			imp.setC(c+1) #1-based...
			IJ.run(imp, "Enhance Contrast", "saturated=0.35")
		
		#Flatten to RGB
		title = imp.getTitle()
		imp.setDisplayMode(1)
		out_ = imp.flatten()
		out_.setTitle(title)
		imps.append(out_)

order = [slides_of_interest.index(i) for i in random_codes]
final_imps = [x for _, x in sorted(zip(order, imps))]

#Add scale bar to last image
IJ.run(final_imps[len(final_imps)-1], "Scale Bar...", "width=100 height=8 font=28 color=White background=None location=[Lower Right] hide overlay");
final_imps[len(final_imps)-1] = final_imps[len(final_imps)-1].flatten()

out_show = ImagesToStack.run(final_imps)
out_show.show() #made the montage manually rather than programatically for these. (See Image -> Stacks -> Make Montage)