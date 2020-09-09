#!/usr/bin/env python3
"""Extract cities from OSM.

It required the xml extracts:
for F in central-america-latest.osm.pbf australia-oceania-latest.osm.pbf south-america-latest.osm.pbf africa-latest.osm.pbf asia-latest.osm.pbf north-america-latest.osm.pbf europe-latest.osm.pbf; do  osmosis --read-pbf ${F} --tf accept-nodes place=city,town,village --tf reject-relations --tf reject-ways  --lp --write-xml ${F/.osm.pbf/.xml}; done

Then you can run as:
for F in *xml; do nohup python parse_osm.py --xml $F --outdir /var/tmp/ekt248fg6h4jh8m90d8bn6d/osm/parsed/ >> ${F/.xml/.log} & sleep 1;  done
"""

import argparse
import time
import os
from os.path import join as pjoin
import inspect
import sys
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from myutils import info, create_readme
import xml.etree.ElementTree as ET
import fiona
import pandas as pd
from shapely.geometry import shape, Point
from scipy.spatial import cKDTree

def load_countries_shp(countriesshp):
    countries = fiona.open(countriesshp)

    centers = np.zeros((30000, 2), dtype=float)
    inds = []

    for i, country in enumerate(countries):
        coords = np.array(country['geometry']['coordinates'])
        
        if len(coords) == 1:
            centers[len(inds), :] = np.mean(np.array(coords[0]), axis=0)
            inds.append(i)
        else:
            for cs in coords:
                centers[len(inds), :] = np.mean(np.array(cs[0]), axis=0)
                inds.append(i)
        
    centers = centers[:len(inds), :]
    tree = cKDTree(centers)
    
    ret = {'countries': countries,
            'inds': inds,
            'tree': tree}
    return ret

##########################################################
def find_country(lat, lon, countriesstruct):
    tree = countriesstruct['tree']

    _, lstids = tree.query([[lon, lat]], k=10)
    inds = np.unique(np.array(countriesstruct['inds'])[lstids[0]])
    
    for country in np.array(countriesstruct['countries'])[inds]:
        point = Point(lon, lat)
        if point.within(shape(country['geometry'])):
            return country['properties']['CNTRY_NAME']
    return ''

##########################################################
def main():
    info(inspect.stack()[0][3] + '()')
    t0 = time.time()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--xml', required=True, help='OSM in xml format')
    parser.add_argument('--outdir', default='/tmp/out/', help='Output directory')
    args = parser.parse_args()

    if not os.path.isdir(args.outdir): os.mkdir(args.outdir)
    readmepath = create_readme(sys.argv, args.outdir)

    tree = ET.parse(args.xml)
    root = tree.getroot()
    child = root[0]

    tagattribs = {
            'place': 'classification',
            'is_in:country': 'country',
            'is_in:state': 'state',
            'name': 'name',
            }

    keys = ['osmid', 'classification', 'name', 'state', 'country',  
        'lat', 'lon', 'timestamp',]

    countriesshp = 'countries/99bfd9e7-bb42-4728-87b5-07f8c8ac631c2020328-1-1vef4ev.lu5nk.shp'
    countries = load_countries_shp(countriesshp)
    
    MAX = 350000
    cities = np.empty((MAX, len(keys)), dtype='object')
    cities[:] = ''

    for i, child in enumerate(root[1:]):
        if child.tag != 'node':
            raise Exception('Something wrong with the xml format.')

        info('i:{}'.format(i))
        cities[i, keys.index('osmid')] = child.attrib['id']
        cities[i, keys.index('lat')] = child.attrib['lat']
        cities[i, keys.index('lon')] = child.attrib['lon']
        cities[i, keys.index('timestamp')] = child.attrib['timestamp']

        for tag in child.iter('tag'):
            if 'k' not in tag.attrib: continue
            attrib = tag.attrib['k']

            if attrib == 'is_in:continent':
                continent = tag.attrib['v']
                continue

            for k, v in tagattribs.items():
                if attrib == k:
                    cities[i, keys.index(v)] = tag.attrib['v']
                    break

        if cities[i, keys.index('country')] == '':
            cities[i, keys.index('country')] = \
                    find_country(float(cities[i, keys.index('lat')]),
                            float(cities[i, keys.index('lon')]), countries)
            
    cities = cities[:i+1]
    df = pd.DataFrame(cities, columns=keys)
    df['continent'] = continent
    filename = os.path.basename(args.xml)
    df.to_csv(pjoin(args.outdir, filename.replace('.xml', '.csv')),
            sep='\t', index=False)

    info('Elapsed time:{}'.format(time.time()-t0))
    info('Output generated in {}'.format(args.outdir))

##########################################################
if __name__ == "__main__":
    main()
