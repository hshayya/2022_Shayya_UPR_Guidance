from ij.plugin.frame import RoiManager
from ij import IJ
from java.awt import Color
from ij.plugin import ImagesToStack, MontageMaker
import re

#Create montage (merged, c1, c2, c3... for all c) from Composite Image, overlaying an ROI from RoiManager in each panel and adding 100um scale bar

which_roi = 2 #0-indexed, Roi from manager to pull for overlay

rm = RoiManager.getInstance()
imp = IJ.getImage()
imp.setOverlay(rm.getRoi(which_roi), Color.white, 5, None)
imp.setDisplayMode(1)

n_channel = imp.getDimensions()[2]
ls = list([['1']*n_channel])
for i in range(n_channel):
	x = ['0']*n_channel
	x[i] = '1'
	ls.append(x)

out = []
for e,i in enumerate(ls):
	imp.setActiveChannels(''.join(i))
	if e == len(ls)-1:
		IJ.run(imp, "Scale Bar...", "width=100 height=8 font=28 color=White background=None location=[Lower Right] hide overlay");
	out.append(imp.flatten())

out_stack = ImagesToStack.run(out)
#out_stack.setOverlay(rm.getRoi(0), Color.white, 5, None)

montage_ = MontageMaker()
out = montage_.makeMontage2(out_stack, len(out), 1, 1, 1, len(out), 1, 2, False)
out.setTitle(re.search('.*(?=\\.)',imp.getTitle()).group(0) + '_Montage')
out.show()