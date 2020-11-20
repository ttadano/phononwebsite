#!/usr/bin/env python
# Copyright (c) 2017, Henrique Miranda
# All rights reserved.
#
# This file is part of the phononwebsite project
#
""" Read phonon dispersion from anaddb """
import argparse
import sys
from phononweb.alamodephonon import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Read anaddb netCDF file and write a .json file to use in the phononwebsite')
    parser.add_argument('filename', help='netCDF filename')
    parser.add_argument('name', help='name of the material')
    parser.add_argument('-r', '--reps', help='number of default repetitions')
    parser.add_argument('-l', '--labels', help='string with the labels of the k-points. Eg. \"GMKG\" ')
    parser.add_argument('-w', '--writeonly', help='only write json file (do not open the web browser)',
                        action="store_true")
    parser.add_argument('-br', '--band_range', help='specify the range of bands to save in json',
                        type=str, default=None, metavar='start end')

    # check if enough arguments are present
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    q = AlamodePhonon(args.filename, args.name, band_range = args.band_range)
    if args.labels: q.set_labels(args.labels)
    if args.reps:   q.set_repetitions(args.reps)

    # diplsay information
    print(q)

    if args.writeonly:
        q.write_json()
    else:
        q.write_json()
        q.open_json()
