import os
import pydicom
import sys
import time
#import matplotlib.pyplot as plt
import numpy as np
from operator import attrgetter, itemgetter
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


"""
We assume RL-AP-IS orientation:
e.g. +X goes to the left, +Y goes to the posterior, and +Z goes to the superior
"""


# constants
MAX = sys.maxint
MIN = -sys.maxint - 1
CT_IMG_UID = 'CT Image Storage'

# read inputs
arg = sys.argv[1]
if arg == '-h' or arg == '--help' or arg == '-help':
	print 'Generates a CT where the voxels of a specified contour are overriden with a specific HU value'
	print 'args: contour_name HU_override_value CT_dir_path RT_contour_file_path out_dir'
	exit()

contour = arg
overr_hu = sys.argv[2]
path_to_ct_dir = sys.argv[3]
path_to_ct_dir = os.path.join(path_to_ct_dir, '')

rt_contour_file = sys.argv[4]
out_dir = sys.argv[5]
out_dir = os.path.join(out_dir, '')

# auxiliary functions
def group(lst, n):
  for i in range(0, len(lst), n):
    val = lst[i:i+n]
    if len(val) == n:
      yield list(val)

def reshapeContourData(contour_data):
	res = list( group(contour_data, 3) )
	for i in range(0,len(res)):
		res[i] = [float(j) for j in res[i]]
		del res[i][-1]
	res = np.array(res)
	return res

def findContourSlice(z, contour_slices):
	for cs in contour_slices:
		cs_z = cs.ContourData[-1]
		if z == cs_z:
			return cs
	return None

# read CT and store some info
lstFilesDCM = []
for dirName, subdirList, fileList in os.walk(path_to_ct_dir):
    for filename in fileList:
        if ".dcm" in filename and filename.lower().startswith('ct'):
            lstFilesDCM.append(os.path.join(dirName,filename))

ds = pydicom.read_file(lstFilesDCM[0])
pixdim = (int(ds.Rows), int(ds.Columns), len(lstFilesDCM))
pixspac = (float(ds.PixelSpacing[0]), float(ds.PixelSpacing[1]), float(ds.SliceThickness))
hu_intercept = ds.RescaleIntercept
overr_hu = float(overr_hu) - float(hu_intercept)
xyoffset = (float(ds.ImagePositionPatient[0]), float(ds.ImagePositionPatient[1]))


# read RT file
rt = pydicom.read_file(rt_contour_file)
rt_contour_idx = -1
found_contours = []
for i in range(0, len(rt.StructureSetROISequence)):
	cname = rt.StructureSetROISequence[i].ROIName
	found_contours.append(cname)
	if cname == contour:
		rt_contour_idx = i

if rt_contour_idx == -1:
	print contour, 'not found in RT DICOM. The following were found:'
	for cname in found_contours:
		print cname
	exit()

contour_slices = sorted(rt.ROIContourSequence[rt_contour_idx].ContourSequence, key=lambda x: x.ContourData[-1])
lowest_contour_slice_z = contour_slices[0].ContourData[-1]
highest_contour_slice_z = contour_slices[-1].ContourData[-1]

full_CT = []
for i in range(0,len(lstFilesDCM)):
	sce = pydicom.read_file(lstFilesDCM[i])
	if not hasattr(sce, 'SliceLocation'):
		continue	
	full_CT.append(sce)

full_CT = np.asarray(full_CT)
full_CT = sorted(full_CT, key=attrgetter("SliceLocation"))


for i in range(0,len(full_CT)):

	cts = full_CT[i]
	rz = cts.SliceLocation

	if rz < lowest_contour_slice_z or rz > highest_contour_slice_z:
		continue

	max_xy_contour = [MIN,MIN]
	min_xy_contour = [MAX,MAX]

	contour_slice = findContourSlice(rz, contour_slices)
	if contour_slice == None:
		continue

	pp_contour_data = reshapeContourData(contour_slice.ContourData)

	max_xy_contour = np.amax(pp_contour_data, axis=0)
	min_xy_contour = np.amin(pp_contour_data, axis=0)
	
	shapely_pp_contour = Polygon( list(tuple(map(tuple, pp_contour_data))) )
	
	x_voxels = cts.pixel_array.shape[0]
	y_voxels = cts.pixel_array.shape[1]

	new_data = np.copy(cts.pixel_array)

	for x in range(0, x_voxels):
		rx = (x)*pixspac[0] + xyoffset[0]
		if rx < min_xy_contour[0] or rx > max_xy_contour[0]:
			continue
		for y in range(0, y_voxels):
			ry = (y)*pixspac[1] + xyoffset[1]
			if ry < min_xy_contour[1] or ry > max_xy_contour[1]:
				continue
			p = Point(rx,ry)
			if shapely_pp_contour.contains(p):
				# y is inverted with x because of how data is stored in DICOM
				new_data[y][x] = overr_hu


	cts.PixelData = new_data.tobytes()
	#plt.imshow(cts.pixel_array, cmap=plt.cm.bone)
	#plt.show(block=False)
	#time.sleep(0.5)
	#plt.close()

# save the DICOM file
for i in range(0,len(full_CT)):
	cts = full_CT[i]
	cts.save_as(out_dir+'CT_'+str(i)+'.dcm')