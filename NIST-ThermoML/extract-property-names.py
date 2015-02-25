#!/usr/bin/env python

"""
Extract property names from ThermoML

"""

import re
import glob

property_counts = dict()
substances = dict()

for filename in glob.glob('*/*.xml'):
    try:
        # Read file contents.
        file = open(filename, 'r')
        contents = file.read()
        file.close()
        # Look for ePropName tags.
        property_names = re.findall(u'<ePropName>(.+)</ePropName>', contents)
        substance_names = re.findall(u'<sCommonName>(.+)</sCommonName>', contents)
        for property_name in property_names:
            if property_name in property_counts:
                property_counts[property_name] += 1
                for substance_name in substance_names:
                    substances[property_name].add(substance_name)
            else:
                property_counts[property_name] = 1
                substances[property_name] = set()
    except Exception as e:
        pass

property_names = sorted(list(property_counts.keys()))
for property_name in property_names:
    print "%5d %s" % (property_counts[property_name], property_name)
#    substance_list = sorted(list(substances[property_name]))
#    for substance in substance_list:
#        print "        %s" % substance


