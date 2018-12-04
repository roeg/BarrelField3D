#!/usr/bin/python

#######################################
# prepare .hoc files from NeuroConv   #
# for use with NeuroMap:              #
# Pia and WM contours are checked and #
# (if necessary) corrected for equal  #
# spacing and integer z-values        #
#                                     #
# Written by Robert Egger             #
#######################################

import sys

def read_hoc_file(hoc_filename):
	read_pia_z = 0
	read_wm_z = 0
	pia_z_list = []
	wm_z_list = []
	hoc_file = open(hoc_filename, 'r')
	for line in hoc_file:
		if line and 'create' in line:
			if 'alpha' in line:
				read_pia_z = 1
			else:
				read_pia_z = 0
			if 'WM' in line:
				read_wm_z = 1
			else:
				read_wm_z = 0
		if line and 'pt3dadd' in line:
			z_val = line.split(',')[2]
			z_val = float(z_val)
			if read_pia_z and not pia_z_list.count(z_val):
				pia_z_list.append(z_val)
				read_pia_z = 0
			elif read_pia_z and pia_z_list.count(z_val):
				err_str = 'Multiple pia contours at z = %f in file %s' % (z_val, hoc_filename)
				raise RuntimeError(err_str)
			if read_wm_z and not wm_z_list.count(z_val):
				wm_z_list.append(z_val)
				read_wm_z = 0
			elif read_wm_z and wm_z_list.count(z_val):
				err_str = 'Multiple WM contours at z = %f in file %s' % (z_val, hoc_filename)
				raise RuntimeError(err_str)
		if line and 'EOF' in line:
			hoc_file.close()
			return pia_z_list, wm_z_list
	hoc_file.close()
	return pia_z_list, wm_z_list

def write_hoc_file(hoc_filename, pia_z_map, wm_z_map):
	read_pia_z = 0
	read_wm_z = 0
	pia_z_list = []
	wm_z_list = []
	out_filename = hoc_filename[:-4]
	out_filename += '_NeuroMap_ready.hoc'
	hoc_file = open(hoc_filename, 'r')
	out_file = open(out_filename, 'w')
	out_file.write('/***************************************/\n')
	out_file.write('/* hoc file from NeuroConv prepared    */\n')
	out_file.write('/* for use with NeuroMap:              */\n')
	out_file.write('/* Pia and WM contours are checked and */\n')
	out_file.write('/* corrected for equal spacing and     */\n')
	out_file.write('/* integer z-values                    */\n')
	out_file.write('/***************************************/\n')
	out_file.write('/*Pia z mapping:                       */\n')
	pia_key_list = pia_z_map.keys()
	pia_key_list.sort()
	pia_key_list.reverse()
	for pia_key in pia_key_list:
		tmp_str = '/* '
		tmp_str += str(pia_key)
		tmp_str += '\t->\t'
		tmp_str += str(pia_z_map[pia_key])
		tmp_str += '\t*/\n'
		out_file.write(tmp_str)
	if len(wm_z_map):
		out_file.write('/***************************************/\n')
		out_file.write('/*WM z mapping:                        */\n')
		wm_key_list = wm_z_map.keys()
		wm_key_list.sort()
		wm_key_list.reverse()
		for wm_key in wm_key_list:
			tmp_str = '/* '
			tmp_str += str(wm_key)
			tmp_str += '\t->\t'
			tmp_str += str(wm_z_map[wm_key])
			tmp_str += '\t*/\n'
			out_file.write(tmp_str)
	out_file.write('/***************************************/\n')
	for line in hoc_file:
		if line and '/*' in line and not 'EOF' in line:
			out_file.write(line)
			continue
		if line and 'create' in line:
			if 'alpha' in line:
				read_pia_z = 1
			else:
				read_pia_z = 0
			if 'WM' in line:
				read_wm_z = 1
			else:
				read_wm_z = 0
			out_file.write(line)
			continue
		if line and 'pt3dadd' in line:
			line_parts = line.split(',')
			z_val = line_parts[2]
			z_val = float(z_val)
			if read_pia_z:
				z_val = pia_z_map[z_val]
			if read_wm_z:
				z_val = wm_z_map[z_val]
			out_file.write(line_parts[0]+',')
			out_file.write(line_parts[1]+',')
			out_file.write(str(z_val)+',')
			out_file.write(line_parts[3])
			continue
		if line and 'EOF' in line:
			break
		out_file.write(line)
	hoc_file.close()
	out_file.close()

def get_z_val_map(z_val_list):
	z_map = {}
	z_map_list = []
	z_val_list.sort()
	z_val_list.reverse()
	#print z_val_list
	avg_dz = 0.0
	for ii in range(0, len(z_val_list)-1):
		avg_dz += z_val_list[ii] - z_val_list[ii+1]
	if len(z_val_list) > 1:
		avg_dz /= (len(z_val_list)-1)
	avg_dz = 10*round(avg_dz/10)
	for ii in range(len(z_val_list)):
		z_map_list.append(round(z_val_list[0])-ii*avg_dz)
	z_map = dict(zip(z_val_list, z_map_list))
	return z_map

#def fix_neuroconv_file(filename):
#	pia_z_list, wm_z_list = read_hoc_file(filename)
#	pia_z_map = get_z_val_map(pia_z_list)
#	wm_z_map = get_z_val_map(wm_z_list)
#	write_hoc_file(filename, pia_z_map, wm_z_map)

def is_identity(map):
	for z_before in map.keys():
		z_after = map[z_before]
		if z_after != z_before:
			return False
	return True

if __name__ == '__main__':
	hoc_file = ''
	if len(sys.argv) > 1:
		hoc_file = sys.argv[1]
	else:
		hoc_file = raw_input('Enter .hoc filename:')
	pia_z_list, wm_z_list = read_hoc_file(hoc_file)
	pia_z_map = get_z_val_map(pia_z_list)
	wm_z_map = get_z_val_map(wm_z_list)
	#print pia_z_map
	#print wm_z_map
	if not is_identity(pia_z_map) or not is_identity(wm_z_map):
		write_hoc_file(hoc_file, pia_z_map, wm_z_map)
