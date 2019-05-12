import os
import pydicom
import sys
import time
#import matplotlib.pyplot as plt
import numpy as np
from operator import attrgetter, itemgetter
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import vtk
import scipy
import scipy.ndimage
import cv2

"""
We assume RL-AP-IS orientation:
e.g. +X goes to the left, +Y goes to the posterior, and +Z goes to the superior
"""


##### CONSTANTS #####
MAX = sys.maxint
MIN = -sys.maxint - 1
CT_IMG_UID = 'CT Image Storage'

##### INPUTS #####
arg = sys.argv[1]
if arg == '-h' or arg == '--help' or arg == '-help':
	print 'Implants the organ of donor into the body of receiver. An RTSTRUCT DICOM file is returned.'
	print 'Args: refpoint_receiver_RL,AP,IS refpoint_donor_RL,AP,IS struct_displacement_RL,AP,IS CT_receiver_dir CT_donor_dir RTSTRUCT_receiver RTSTRUCT_donor organ_name out_dir'
	print 'E.g.: -12.0,203.0,1020.0 33.0,220.0,110.0 55.0,-10.0,25.5 path/to/receiver/CTdir/ path/to/donor/CTdir/ path/to/receiver/RTSTRUCTFile path/to/donor/RTSTRUCTFile Liver path/to/outputdir/'
	exit()

coords_ref_rec = sys.argv[1]
coords_ref_donor = sys.argv[2]
coords_struct_displacement = sys.argv[3]
ctdir_rec = sys.argv[4]
ctdir_rec = os.path.join(ctdir_rec, '')
ctdir_donor = sys.argv[5]
ctdir_donor = os.path.join(ctdir_donor, '')
rt_receiver = sys.argv[6]
rt_donor = sys.argv[7]
organ_name = sys.argv[8]
out_dir = sys.argv[9]
out_dir = os.path.join(out_dir, '')

coords_ref_rec = np.array(coords_ref_rec.split(',')).astype(np.float)
coords_ref_donor = np.array(coords_ref_donor.split(',')).astype(np.float)
coords_struct_displacement = np.array(coords_struct_displacement.split(',')).astype(np.float)



##### FUNCTIONS #####
def computeCenterOfMass(polydata):
	centerFilter = vtk.vtkCenterOfMass()
	centerFilter.SetInputData(polydata)
	centerFilter.SetUseScalarsAsWeights(False)
	centerFilter.Update()
	com = np.array( centerFilter.GetCenter() )
	return com


def readCTDicom(dirName):
	lstFilesDCM = []
	for dirName, subdirList, fileList in os.walk(dirName):
	    for filename in fileList:
	        if ".dcm" in filename.lower()[-4:] and "ct" in filename.lower()[:2]:
	            lstFilesDCM.append(os.path.join(dirName,filename))
	full_CT = []
	for i in range(0,len(lstFilesDCM)):
		sce = pydicom.read_file(lstFilesDCM[i])
		if not hasattr(sce, 'SliceLocation'):
			continue
		full_CT.append(sce)

	full_CT = np.asarray(full_CT)
	full_CT = sorted(full_CT, key=attrgetter("SliceLocation")) # sort along IS, increasing
	return full_CT


def findRTContour(rt, organname):
	rt_contour_idx = -1
	found_contours = []
	for i in range(0, len(rt.StructureSetROISequence)):
		cname = rt.StructureSetROISequence[i].ROIName
		found_contours.append(cname)
		if cname == organname:
			rt_contour_idx = i
			break

	if rt_contour_idx == -1:
		print organname, "not found in RTSTRUCT"
		print "The following organs were found:"
		for cname in found_contours:
			print cname
		exit()

	# sort along IS, increasing
	contour_slices = sorted(rt.ROIContourSequence[rt_contour_idx].ContourSequence, key=lambda x: x.ContourData[-1]) 
	return contour_slices

def group(lst, n):
  for i in range(0, len(lst), n):
    val = lst[i:i+n]
    if len(val) == n:
      yield list(val)

def convertContourDataToNPArray(contour_data):
	res = list( group(contour_data, 3) )
	for i in range(0,len(res)):
		res[i] = [float(j) for j in res[i]]
	res = np.array(res)
	return res

def convertContourSequenceToNPArray(seq):
	res = []
	for cs in seq:
		res.append(convertContourDataToNPArray(cs))
	return res


def replaceRTDicomContour(rt_to_replace, rt_replacing, organ_name, translation):
	rt_toreplace_contour_idx = -1
	rt_replacing_contour_idx = -1
	found_contours_toreplace = []
	found_contours_replacing = []
	for i in range(0, len(rt_to_replace.StructureSetROISequence)):
		cname = rt_to_replace.StructureSetROISequence[i].ROIName
		found_contours_toreplace.append(cname)
		if cname == organ_name:
			rt_toreplace_contour_idx = i
			break

	if rt_toreplace_contour_idx == -1:
		print organ_name, "not found in RTSTRUCT of the receiver. The following were found:"
		for cname in found_contours_toreplace:
			print cname
		exit()

	for i in range(0, len(rt_replacing.StructureSetROISequence)):
		cname = rt_replacing.StructureSetROISequence[i].ROIName
		found_contours_replacing.append(cname)
		if cname == organ_name:
			rt_replacing_contour_idx = i
			break

	if rt_replacing_contour_idx == -1:
		print organ_name, "not found in RTSTRUCT of the donor. The following were found:"
		for cname in found_contours_replacing:
			print cname
		exit()

	replacing_sequence = rt_replacing.ROIContourSequence[rt_replacing_contour_idx].ContourSequence
	j = 0
	for r in replacing_sequence:
		data = r.ContourData
		new_data = []
		for d in data:
			j = j % 3
			new_data.append(  str(float(d) + translation[j])  )
			j = j+1
		r.ContourData = new_data

	rt_to_replace.ROIContourSequence[rt_toreplace_contour_idx].ContourSequence = replacing_sequence
	return rt_to_replace

def convert_scan_to_image(ct):
	image = np.stack([s.pixel_array for s in ct])
	image = image.astype(np.int16)   
	return np.array(image, dtype=np.int16)

def resample(scan, new_spacing=[1,1,1]):
	image = convert_scan_to_image(scan)
	# Determine current pixel spacing
	spacing = map(float, [scan[0].SliceThickness]) + map(float, scan[0].PixelSpacing)
	spacing = np.array(list(spacing))

	resize_factor = spacing / new_spacing
	new_real_shape = image.shape * resize_factor
	new_shape = np.round(new_real_shape)
	real_resize_factor = new_shape / image.shape
	new_spacing = spacing / real_resize_factor

	image = scipy.ndimage.interpolation.zoom(image, real_resize_factor)

	return image, new_spacing

def find_nearest(array,value):
	idx = (np.abs(array-value)).argmin()
	return array[idx]

def find_nearest_by_arg(array,value,attr):
	temp = []
	for item in array:
		temp.append(item.attrgetter(attr));
	idx = (np.abs(temp-value)).argmin()
	return array[idx]

def get_nearest_slice_idx(z, slices, spacing, offset):
	# assumes slices are sorted by increasing IS
	dist = MAX
	best = -1
	for i in range(0,len(slices)):
		zs = i*spacing + offset
		cdist = np.abs(z - zs)
		if cdist < dist:
			dist = cdist
			best = i
		else:
			break
	
	return best


def get_nearest_slice_idx_from_CT(z, ct):
	# assumes CT is sorted by increasing IS
	dist = MAX
	best = -1
	for i in range(0,len(ct)):
		zct = float( ct[i].SliceLocation )
		cdist = np.abs(z - zct)
		if cdist < dist:
			dist = cdist
			best = i
		else:
			break

	return best


def get_nearest_voxel_idcs(xy, slic, spacings, xyoffset):
	dists = [MAX,MAX]
	bests = [float('nan'),float('nan')]

	xlen = slic.shape[0]
	ylen = slic.shape[1]

	for x in range(0,xlen):
		rx = x*spacings[0] + xyoffset[0]
		d = np.abs(rx - xy[0])
		if d < dists[0]:
			dists[0] = d
			bests[0] = x

	for y in range(0,ylen):
		ry = y*spacings[1] + xyoffset[1]
		d =  np.abs(ry - xy[1])
		if d < dists[1]:
			dists[1] = d
			bests[1] = y
	
	return np.array( bests )	


def check_organ_nicely_inside( ct_rec, ct_rec_img, left_or_right, contour_donor_np_translated):
	print 'check_organ_nicely_inside organ not implemented'
	exit()
	"""
	if left_or_right != 'left' and left_or_right != 'right':
		print 'invalid value for left or right control parameter:', left_or_right
		exit()
	
	max_x = MIN
	min_x = MAX
	for cont_slic in contour_donor_np_translated :
		x = np.amax(cont_slic, axis=0)
		if (x > max_x):
			max_x = x
			max_x_slic = cont_slic
		x = np.amin(cont_slic, axis=0)
		if (x < min_x):
			min_x = x
			min_x_slic = cont_slic

	if left_or_right == 'left':
		z = min_x_slic[0][-1]	
		rec_slic = find		
	elif left_or_right == 'right':
		z = max_x_slic[0][-1]

	rec_slic_idx = get_nearest_slice_idx( z, ct_rec_img, ct_rec[0].SliceThickness, ct_rec[0].SliceLocation )
	ct_rec_slic = ct_rec_img[rec_slic_idx]

	# now check if the organ is inside the body / does not overlap with bones
	
	# if it is not good, compute delta x to put it inside
	# while( not good enough) if left, x - delta x; if right, x + delta x (keep track of cumulative translation)
	#theta_prime = np.arccos( x_prime / com_organ_body_distance )
	#z_prime = com_organ_body_distance * np.sin(theta_prime)
	
	# translate everything by x_prime, y_prime
	#for cont_slic in contour_donor_np_translated:
	#	cont_slic = cont_slic + [x_prime - x, y_prime - y, z_prime - z]

	return True
	"""


def reposition_organ():
	print 'reposition_organ not implemented'
	exit()

def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
	# Code by aubricus: https://gist.github.com/aubricus/f91fb55dc6ba5557fbab06119420dd6a
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '#' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()



##### PROGRAM #####

# Read data
print 'Reading CT receiver... ',
ct_rec = readCTDicom(ctdir_rec)
if len(ct_rec) > 0:
	print 'OK'
else:
	print 'CT receiver reading error'
	exit()

"""
You can use the following snippet to assess the orientation and change the code:
"""
"""for dd in dir(ct_rec[0]):
	if dd.startswith('Image') or dd.startswith('Anatomical')  or dd.startswith('Patient'):
		try:
			print dd, ':', getattr(ct_rec[0], dd)
		except:
			continue
"""


print 'Reading CT donor... ',
ct_donor = readCTDicom(ctdir_donor)
if len(ct_donor) > 0:
	print 'OK'
else:
	print 'CT donor reading error'
	exit()


xyoffset_donor = (float(ct_donor[0].ImagePositionPatient[0]), float(ct_donor[0].ImagePositionPatient[1]))
zoffset_donor = float(ct_donor[0].SliceLocation)
pixspac_donor = (float(ct_donor[0].PixelSpacing[0]), float(ct_donor[0].PixelSpacing[1]), float(ct_donor[0].SliceThickness))
xyoffset_rec = (float(ct_rec[0].ImagePositionPatient[0]), float(ct_rec[0].ImagePositionPatient[1]))

zoffset_rec = float(ct_rec[0].SliceLocation)
pixspac_rec = (float(ct_rec[0].PixelSpacing[0]), float(ct_rec[0].PixelSpacing[1]), float(ct_rec[0].SliceThickness))


# Read contour to trasplant and transform its coordinates to be relative to its center of mass
print 'Reading RTSTRUCT receiver... ',
rt_receiver_dicom = pydicom.read_file(rt_receiver)
if len(rt_receiver_dicom) == 0:
	print 'RTSTRUCT receiver reading error'
	exit()
print 'Looking for',organ_name,'... ',
contour_receiver = findRTContour(rt_receiver_dicom, organ_name)
print 'OK'


print 'Reading RTSTRUCT donor... ',
rt_donor_dicom = pydicom.read_file(rt_donor)
if len(rt_donor_dicom) == 0:
	print 'RTSTRUCT receiver reading error'
	exit()
print 'Looking for',organ_name,'... ',
contour_donor = findRTContour(rt_donor_dicom, organ_name)
print 'OK'


# Create VTK representation of organ, determine Center of Mass, then translate it
print 'Translating',organ_name,'of donor...',

vertices = vtk.vtkCellArray()
points = vtk.vtkPoints()

contour_donor_np = []
for cs in contour_donor:
	conv = convertContourDataToNPArray(cs.ContourData)
	contour_donor_np.append(conv)		
	for cc in conv:
		id = points.InsertNextPoint(cc)
		vertices.InsertNextCell(1)
		vertices.InsertCellPoint(id)	
contour_donor_np = np.array(contour_donor_np)

vtkPolyStruct = vtk.vtkPolyData()
vtkPolyStruct.SetPoints(points)
vtkPolyStruct.SetVerts(vertices)
com_contour_donor = computeCenterOfMass(vtkPolyStruct)

contour_donor_np_translated = []
for cs in contour_donor_np:
	contour_donor_np_translated.append( cs - com_contour_donor + coords_ref_rec + coords_struct_displacement )
contour_donor_np_translated = np.array(contour_donor_np_translated)


print 'OK'


# Prepare CT images
print 'Preparing to transplant...',
ct_donor_img = convert_scan_to_image(ct_donor) 
 # adjust for different HU values intercepts
intercept_diff = ct_rec[0].RescaleIntercept - ct_donor[0].RescaleIntercept
ct_donor_img = ct_donor_img - intercept_diff
ct_rec_img = convert_scan_to_image(ct_rec)

# Insert organ donor into body receiver
vals = []

x_voxels_num = ct_rec_img[0].shape[0]
y_voxels_num = ct_rec_img[0].shape[1]

print 'OK'

# for test purpose, painting the donor organ (no rescaling)

print ' > Transplanting',organ_name,'of donor into receiver'

prev_slic_idx = -1
prev_slic_contour_application = []
for i_o_slic in range(0,len(contour_donor_np_translated)):
	
	# pick respective slice for receiver
	o_don_slic = contour_donor_np[i_o_slic]
	o_don_slic_translated = contour_donor_np_translated[i_o_slic]
	z_o_don_slic_translated = o_don_slic_translated[0][-1]
	ct_slic_idx = get_nearest_slice_idx_from_CT( z_o_don_slic_translated, ct_rec )
	ct_rec_slic = ct_rec_img[ct_slic_idx]

	# propagate previous contour if the CT of the receiver is oversampled compared to the contour
	if prev_slic_idx != -1 and ct_slic_idx - prev_slic_idx > 1:
		for i_prev_slic in range(prev_slic_idx, ct_slic_idx):
			prev_slic = ct_rec_img[i_prev_slic]
			# propagate previous to this intermediate slices
			for item in prev_slic_contour_application:
				prev_slic[item[0]][item[1]] = item[2]
				vals.append(item[2])
			ct_rec[i_prev_slic].PixelData = prev_slic.tobytes()
		prev_slic_contour_application = []
			
	# pick respective slice for donor
	ct_don_slic = ct_donor_img[  get_nearest_slice_idx_from_CT( o_don_slic[0][-1], ct_donor ) ]
	
	max_xy_contour = np.amax(o_don_slic_translated, axis=0)
	min_xy_contour = np.amin(o_don_slic_translated, axis=0)

	# Prepare voxel-indices version of the organ slice for the receiver CT
	pix_o_don_slic_translated = []
	pix_o_don_slic = []
	
	for i in range(0, len(o_don_slic_translated[:,:-1])):
		xy = o_don_slic[:,:-1][i]
		xy_t = o_don_slic_translated[:,:-1][i]
		pix = np.round( (xy - xyoffset_donor) / pixspac_donor[:-1] ).astype(int)
		pix_t = np.round( (xy_t - xyoffset_rec) / pixspac_rec[:-1] ).astype(int)
		# the following inversions between x and y is due to how pixeldata is stored
		pix_o_don_slic.append( [pix[1],pix[0]] )	
		pix_o_don_slic_translated.append([pix_t[1],pix_t[0]])

	pix_o_don_slic = np.array(pix_o_don_slic)
	pix_o_don_slic_translated = np.array(pix_o_don_slic_translated)
	
	minpix_o_don_slic = np.min(pix_o_don_slic, axis=0)
	maxpix_o_don_slic = np.max(pix_o_don_slic, axis=0)
	minpix_o_don_slic_translated = np.min(pix_o_don_slic_translated, axis=0)
	maxpix_o_don_slic_translated = np.max(pix_o_don_slic_translated, axis=0)

	poly_o_don_slic = Polygon( list(tuple(map(tuple, pix_o_don_slic ))) )
	poly_o_don_slic_translated = Polygon( list(tuple(map(tuple, pix_o_don_slic_translated ))) )

	# Create matrix of HU values that need to be written
	hu_vals = np.ndarray( ( maxpix_o_don_slic[0] - minpix_o_don_slic[0], maxpix_o_don_slic[1] - minpix_o_don_slic[1]) )
	hu_vals[:,:] = np.nan
	for i, x in enumerate( range( minpix_o_don_slic[0], maxpix_o_don_slic[0] )):
		for j, y in enumerate( range( minpix_o_don_slic[1], maxpix_o_don_slic[1] )):
			if poly_o_don_slic.contains(Point(x,y)):
				hu_vals[i][j] = ct_don_slic[x][y]
				vals.append(ct_don_slic[x][y])
	medianHU = np.median(vals)
	where_are_NaNs = np.isnan(hu_vals)
	hu_vals[where_are_NaNs] = medianHU
	hu_vals = cv2.resize(hu_vals, dsize=(maxpix_o_don_slic_translated[1]-minpix_o_don_slic_translated[1], maxpix_o_don_slic_translated[0]-minpix_o_don_slic_translated[0]))

	# Overwrite HU values of voxels
	for i, x in enumerate( range( minpix_o_don_slic_translated[0] , maxpix_o_don_slic_translated[0]) ):	
		for j, y in enumerate( range(minpix_o_don_slic_translated[1] , maxpix_o_don_slic_translated[1]) ):
			if poly_o_don_slic_translated.contains(Point(x,y)):
				hu_to_write = hu_vals[i][j]
				if np.isnan(hu_to_write):
					hu_to_write = medianHU
				ct_rec_slic[x][y] = hu_to_write
				prev_slic_contour_application.append( (x,y,hu_to_write) )
			
	# Overwrite CT slice
	ct_rec[ct_slic_idx].PixelData = ct_rec_slic.tobytes()
	prev_slic_idx = ct_slic_idx

	print_progress(i_o_slic + 1, len(contour_donor_np_translated))

#print 'HU values of the inserted organ -- Median:', np.median(vals), '25th perc:',np.percentile(vals,25),'75th perc:',np.percentile(vals,75)

print 'Saving CT and RTSTRUCT files in',out_dir,'...',
for i in range(0,len(ct_rec)):
	cts = ct_rec[i]
	cts.save_as(out_dir+'CT_'+str(i)+'.dcm')

rt_receiver_dicom = pydicom.read_file(rt_receiver)

last_translation = ( coords_ref_rec + coords_struct_displacement - com_contour_donor )

new_rt_receiver = replaceRTDicomContour( rt_receiver_dicom, rt_donor_dicom, organ_name, last_translation )
new_rt_receiver.save_as(out_dir+'RT_replaced.dcm')


contour_donor = findRTContour(new_rt_receiver, organ_name)
print 'OK'

"""
The following code can be used to check that the 
center of mass of the transplant organ went where intended
"""

"""
print 'Checking C.O.M.',organ_name,'after transplant...',
# create vtk representation of organ to then translate it 
# so that it has center of mass in origin
vertices = vtk.vtkCellArray()
points = vtk.vtkPoints()

contour_donor_np = []
for cs in contour_donor:
	conv = convertContourDataToNPArray(cs.ContourData)
	contour_donor_np.append(conv)		
	for cc in conv:
		id = points.InsertNextPoint(cc)
		vertices.InsertNextCell(1)
		vertices.InsertCellPoint(id)	
contour_donor_np = np.array(contour_donor_np)

vtkPolyStruct = vtk.vtkPolyData()
vtkPolyStruct.SetPoints(points)
vtkPolyStruct.SetVerts(vertices)
com_contour_donor = computeCenterOfMass(vtkPolyStruct)

print 'obtained:',np.round(com_contour_donor,2), 'vs expected:', coords_ref_rec + coords_struct_displacement

print 'OK'
"""