"sort_gnashyfile(filename)"

import argparse
import pandas as pd

def sort_gnashyfile(gnashyfilename):
    gnashy_frame = pd.read_csv(gnashyfilename,delimiter='\t',header=None,names=['chr','pos','reads'])
    gnashy_frame = gnashy_frame.sort_values(['chr','pos'])
    gnashy_frame.to_csv(gnashyfilename,sep='\t',header=False,index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='sort_gnashyfile.py')
    parser.add_argument('-i','--infilename',help='gnashy filename (full path)',required=True)
    args = parser.parse_args()
    sort_gnashyfile(args.infilename)
