#!/usr/bin/env python2

import argparse
import subprocess


def main():

    parser=argparse.ArgumentParser(description="Generate symlinks for raw files from metadata GEO table")
    parser.add_argument('-i', help='input metadata table (ex. GEO_Houda_2019_07_14.txt)', type=str, required=True)
    args=parser.parse_args()

    with open(args.i) as FH:
        lines = FH.readlines()
        raw = False
        skipheader = False
        for i, line in enumerate(lines):
            splitline = line.strip().split('\t')
            if raw and splitline[0] == '':
                raw = False
            if raw and skipheader:
                if splitline[0][:4] == 'run2':
                    fpath = splitline[6]
                    f = fpath.split('/')[-1]
                    subprocess.call('ln -s ' + splitline[6] + ' run2-' + f, shell=True)
                else:
                    subprocess.call('ln -s ' + splitline[6], shell=True)
            if splitline[0] == 'RAW FILES':
                raw = True
            if splitline[0] == 'file name' and raw:
                skipheader = True
            
if __name__ == '__main__':
    main()
