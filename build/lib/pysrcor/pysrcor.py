#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from collections import defaultdict
import astropy.units as u
from astropy.coordinates import SkyCoord
import math
import sys

def match(ra1, de1, ra2, de2, rad, opt = 2, silent=False):
	'''
	This function will perform a match between two catalogs

	Inputs:
		ra1, de1, ra2, de2: coordiantes in two catalogs, should be array-like unit: deg
		rad: matching radius, unit: arcsec

	Keywords:
		opt: same as srcor.pro in IDL
				0, source in catalog 2 may correspond to >1 sources in catalog 1
				1, forced 1-1 match
				2, forced 1-1 match, but the program will run two times. At the first run is to obtain					 the median coordinate offset between the two catalogs. Then we apply the offset to
				   one catalog to eliminate the systematic difference. The second run is based on the
				   updated coordinates. (default)
		silent: True, do not display results for opt=0 or 1
				False, display results (default)

	Outputs:
		id1, indexes of matched sources in catalog I
        id2, indexes of matched sources in catalog II
		dis, distances of the pairs, unit: arcsec
	'''

	# Format unification
	ra1 = np.array(ra1)
	de1 = np.array(de1)
	ra2 = np.array(ra2)
	de2 = np.array(de2)

	# Check if all RA, DEC are float
	if ( np.where(np.isfinite(ra1))[0].size != ra1.size ) | \
	   ( np.where(np.isfinite(de1))[0].size != de1.size ) | \
	   ( np.where(np.isfinite(ra2))[0].size != ra2.size ) | \
	   ( np.where(np.isfinite(de2))[0].size != de2.size ) : \
		sys.exit('Error: non-float in ra,dec')

	if ra1.shape != de1.shape:
		raise ValueError('ra1 and dec1 do not match!')
	if ra2.shape != de2.shape:
		raise ValueError('ra2 and dec2 do not match!')

	# Initial match, note that id_2 may contain repetitive elements
	if opt == 0:
		id1, id2, dis = _ini_match(ra1, de1, ra2, de2, rad/3600.)
		if not(silent):
			print('Multi-one match: ', len(id1), ' sources')
		return id1, id2, dis*3600.

	# First match, only the nearest neighbor is left when multi-to-one occurs
	id1, id2, dis = _1st_match(ra1, de1, ra2, de2, rad/3600.)
	if opt == 1:
		if not(silent):
			print('Forced one-one match: ', len(id1), ' sources')
		return id1, id2, dis*3600.
	print('First match: ', len(id1), ' sources')

	# Second match, elliminate the offset
	print( 'RA offset:', \
		   np.median(ra1[id1] - ra2[id2]) * 3600 \
           * math.cos(np.deg2rad(np.median(de1))), 'arcsec', \
		   'DEC offset:', \
		   np.median(de1[id1] - de2[id2])*3600, 'arcsec')
	ra2_new = ra2 + np.median(ra1[id1] - ra2[id2])
	de2_new = de2 + np.median(de1[id1] - de2[id2])
	id1, id2, dis = _1st_match(ra1, de1, ra2_new, de2_new, rad/3600.)
	print('Second match: ', len(id1), ' sources')

	return id1, id2, dis*3600


def _1st_match(ra1, de1, ra2, de2, rad):
	'''
	This function performs the forced 1-1 match. Similar to srcor.pro (opt=1) in IDL.

	Inputs:
		ra1, de1, ra2, de2: coordiantes in two catalogs, unit: deg
		rad: matching radius, unit: deg

	Outputs:
		id1, indexes of matched sources in catalog I
        id2, indexes of matched sources in catalog II
		dis, distances of the pairs, unit: deg
    '''

	# Initial match, might produce multi-one match
	id1, id2, dis = _ini_match(ra1, de1, ra2, de2, rad)

	# Find the repetitive ones  (routine)
	rep = defaultdict(list)
	for i, item in enumerate(id2):
		rep[item].append(i)
	rep = {k:v for k,v in rep.items() if len(v) > 1}
	# note that k is indexes in catalog II, v is indexes in the matching array

	# For them, reserve the nearest one only
	kill = []
	for dum, idx in rep.items():
		dis_rep = dis[idx]	# Their distances
		best = dis_rep.argmin()  # Pick up the nearest one
		kill.append(np.delete(idx, best)) # Creat the to-be-delete list

	buf = np.arange(0)	# Transform format from list to numpy array
	while kill != []:
		buf = np.append(buf, kill.pop())
	kill = buf

	id1 = np.delete(id1, kill)	# Remove them
	id2 = np.delete(id2, kill)
	dis = np.delete(dis, kill)

	return id1, id2, dis


def _ini_match(ra1, de1, ra2, de2, rad):
	'''
	This function is modified from pyspherematch.py, available at
	https://gist.github.com/eteq/4599814
	The main change is, turn to match_coordinates_sky (astropy) from KDTree (scipy).
	Because the second is based on 'straight-line' distance rather than
	'great-circle' distance.

	Inputs:
		ra1, de1, ra2, de2: coordiantes in two catalogs, unit: deg
		rad: matching radius, unit: deg

	Outputs:
		id1, indexes of matched sources in catalog I
        id2, indexes of matched sources in catalog II
		dis, distances of the pairs, unit: deg
	'''

	cor1 = SkyCoord(ra=ra1*u.deg, dec=de1*u.deg, frame='icrs')	 # Change to ICRS system
	cor2 = SkyCoord(ra=ra2*u.deg, dec=de2*u.deg, frame='icrs')

	id2, dis, buf = cor1.match_to_catalog_sky(cor2)	 # Match
	dis = np.array(dis)	# Drop the unit
	id1 = np.arange(ra1.size)
	#print dis, rad
	if id2.size<=1:
		id2 = np.array([id2])
	# I changed id2 = np.array(id2) -> id2 = np.array([id2])
	# So it works for the case when id2 is a single value
	#											12/14/2015
	#												 Guang
	else:
		id2 = np.array(id2)

	msk = dis < rad	 # Filter out those with distance > rad
	id1 = id1[msk]
	id2 = id2[msk]
	dis = dis[msk]

	return id1, id2, dis
