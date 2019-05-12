import numpy as np
from lookup_table_sDSC import lookup as LOOK
from xxx_to_xxx import XXX_to_XXX
import os



PATS = [000]

OARs = ['Liver','Spleen']

DIR_CTs = 'xxx/'
DIR_CTs_XXX = 'xxx/'
HU_softtissue = 78	# like NCI phantoms

rt_struct_file_prefix = 'rs'

for pat in PATS:

	os.system('rm -r temp/*')
	os.system('rm -r temp_resection/*')
	os.system('rm -r temp_transplant/*')

	if pat <= 0:
		continue

	body_donor = LOOK['DSC']['Body'][pat]
	print 'Getting body for pat',pat,'from pat', body_donor

	ct_dir = DIR_CTs + 'XXX_'+str(body_donor)
	if body_donor > 000:
		ct_dir = DIR_CTs_XXX + XXX_to_XXX['xxx_'+str(body_donor)]

	rtstruct_file = [x for x in os.listdir(ct_dir) if rt_struct_file_prefix in x.lower()]

	if len(rtstruct_file) == 0:
		rtstruct_file = [x for x in os.listdir(ct_dir) if 'rtstruct' in x.lower()]
	if len(rtstruct_file) == 0:
		print 'Could not find RTSTRUCT for', pat
		continue

	rtstruct_file = ct_dir+'/'+rtstruct_file[0]


	some_donor_OAR_does_not_exist = False
	for oar in OARs:
		# before proceeding, assess that the OAR of the donor exists
		donor = LOOK['DSC'][oar][pat]

		don_dir = DIR_CTs + 'XXX_'+str(donor)
		if donor > 000:
			don_dir = DIR_CTs_XXX + XXX_to_XXX['XXX_'+str(donor)]


		rtstruct_donor_file = [x for x in os.listdir(don_dir) if rt_struct_file_prefix in x.lower()]
		if len(rtstruct_donor_file) == 0:
			rtstruct_donor_file = [x for x in os.listdir(don_dir) if 'rtstruct' in x.lower()]

		if len(rtstruct_donor_file) == 0:
			print 'Could not find RTSTRUCT for', donor
			some_donor_OAR_does_not_exist = True
			break

	if some_donor_OAR_does_not_exist:
		continue
	
	os.system('cp -r '+ct_dir+'/* temp')

	for oar in OARs:
		print 'Resection of',oar
		os.system('python ResectNTransplant/resect.py '+oar+' '+str(HU_softtissue)+' temp/ "'+rtstruct_file+'" temp_resection/')
		os.system('rm temp/*')
		os.system('cp -r temp_resection/* temp/')

	os.system('cp -r temp_resection/* temp')
	os.system('cp '+rtstruct_file+' temp/RT_replaced.dcm')
	
	for oar in OARs:
		donor = LOOK['DSC'][oar][pat]
		print 'Transplant of',oar,'using as donor pat',donor

		# reference point in the receiver CT scan
		rp_rec = ( -LOOK['X']['Ref'][body_donor], -LOOK['Y']['Ref'][body_donor], LOOK['Z']['Ref'][body_donor] )
		rp_don = ( -LOOK['X']['Ref'][donor], -LOOK['Y']['Ref'][donor], LOOK['Z']['Ref'][donor] )
		displ = ( -LOOK['X'][oar][pat], -LOOK['Y'][oar][pat], LOOK['Z'][oar][pat] )

		don_dir = DIR_CTs + 'XXX_'+str(donor)
		if donor > 000:
			don_dir = DIR_CTs_XXX + XXX_to_XXX['XXX_'+str(donor)]
		rtstruct_donor_file = [x for x in os.listdir(don_dir) if rt_struct_file_prefix in x.lower()]
		if len(rtstruct_donor_file) == 0:
			rtstruct_donor_file = [x for x in os.listdir(don_dir) if 'rtstruct' in x.lower()]

		if len(rtstruct_donor_file) == 0:
			print 'Could not find RTSTRUCT for', donor
			continue
		rtstruct_donor_file = don_dir+'/'+rtstruct_donor_file[0]

		os.system('python ResectNTransplant/transplant.py '+
			str(rp_rec[0])+','+str(rp_rec[1])+','+str(rp_rec[2])+' '+
			str(rp_don[0])+','+str(rp_don[1])+','+str(rp_don[2])+' '+
			str(displ[0])+','+str(displ[1])+','+str(displ[2])+' '+
			'temp/' +' '+don_dir+' temp/RT_replaced.dcm '+rtstruct_donor_file+' '+
			oar + ' temp_transplant/'
		)
		os.system('rm -r temp/*')
		os.system('cp -r temp_transplant/* temp/')

	os.system('mkdir '+str(pat))
	os.system('cp -r temp_transplant/* '+str(pat))

	os.system('rm -r temp/*')
	os.system('rm -r temp_resection/*')
	os.system('rm -r temp_transplant/*')
