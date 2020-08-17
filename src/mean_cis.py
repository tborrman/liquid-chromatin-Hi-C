
from __future__ import print_function

import numpy as np
import scipy as sp
import pdb
import h5py
import sys
import argparse
import time
import gzip
import re
import uuid
import math
import os

def main():

    parser=argparse.ArgumentParser(description='Extract mean cis interactions per bin from hdf5 file.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-in',help='interaction matrix hdf5 file',dest='infile',type=str,required=True)
    
    parser.add_argument('-chrs',help='subset of chromosomes to extract (default=all)',dest='rel_chrs',nargs='+',type=str)
    parser.add_argument('-b',help='block size for extracting (default=hdf chunk size)',dest='blocksize',type=int)
    parser.add_argument('-p',help='output precision (# of digits)',dest='precision',type=int,default=4)
    
    args=parser.parse_args()

    infile=args.infile
    # this argument doesn't work
    rel_chrs=args.rel_chrs
    blocksize=args.blocksize
    precision=args.precision
    
    print("\n",end="")
    
    infile_name=os.path.basename(infile)
    infile_name=re.sub(".hdf5", "", infile_name)
    inhdf=h5py.File(infile,'r')

    # build chr dict    
    chrs=inhdf['chrs'][:]
    chr_dict={}
    for i,c in enumerate(chrs):
        chr_dict[c]=i
    
    print("available chromosomes\n",end="")
    for i,c in enumerate(chrs):
        print("\t",c,"\n",end="")
    
    print("\n",end="")
   
    print("selected chromosomes\n",end="")
    cis_mode=0
    if rel_chrs!=None:
        # check relevant chromosomes    
        for c in rel_chrs:
            print("\t",c,"\n",end="")
            if not c in chr_dict:
                sys.exit('specificed chr '+c+'not found in file!')
        # ensure rel_chrs are sorted same as the HDF
        rel_chrs=sorted(rel_chrs,key=lambda x:chr_dict[x])
        cis_mode=1
    else:
        rel_chrs=chrs
        print("\tall\n",end="")
        cis_mode=0
        
    print("\n",end="")
    
    format_func=("{:0."+str(precision)+"f}").format
    
    bin_positions=inhdf['bin_positions'][:]
    genome=inhdf.attrs['genome'][:]
    chr_bin_range=inhdf['chr_bin_range'][:]
    n=inhdf['interactions'].shape[0]
    bin_mask=np.zeros(n,dtype=bool)
    chr_bin_mask = {}
    
   
    for i,c in enumerate(rel_chrs):
        #build bin mask
        c_ind=chr_dict[c]
        
        tmp_chr_bin_mask=np.zeros(n,dtype=bool)
        r=chr_bin_range[c_ind]
        tmp_chr_bin_mask[r[0]:r[1]+1]=True
        bin_mask[r[0]:r[1]+1]=True
        
        chr_bin_mask[c]=tmp_chr_bin_mask
          
    if blocksize==None:
        blocksize=inhdf['interactions'].chunks[0]

    print("matrix is",n,"x",n)
    print("blocksize is",blocksize)
    
    print("")
    print("calculating mean cis ...")
    
    
    meancis_bedGraphFile=infile_name+".meancis.bedGraph"
    meancis_fh=open(meancis_bedGraphFile, "w")
    
    k=0
    for c in rel_chrs:
        c_ind=chr_dict[c]
        c_start,c_end=inhdf['chr_bin_range'][c_ind]
        
        # redefine chr bounds by bin mask 
        tmp_bin_mask=bin_mask[c_start:c_end+1]
       
        c_end=max(np.nonzero(tmp_bin_mask))[-1]+c_start
        c_start=max(np.nonzero(tmp_bin_mask))[0]+c_start
              
        for i in xrange(c_start,c_end+1,blocksize):
            b=min(c_end-i+1,blocksize)
            current_block=inhdf['interactions'][i:i+b,:][:,bin_mask] # hdf to memory
            #print(inhdf['interactions'][i:i+b,:][:,bin_mask].shape)
            # current block is interaction matrix with row starting on chrom but number of rows = size of blocksize
            # number of cols  = number of bins in chrom
            
            for j in xrange(current_block.shape[0]):
                
                cis_mask=chr_bin_mask[c]
                cis_mean=np.nanmean(current_block[j,cis_mask])
                                    
                chr_id=deGroupChr(chrs[bin_positions[k][0]])
                bin_start=bin_positions[k][1]
                bin_end=bin_positions[k][2]
                
                # progress bar
                if(k != 0):
                    sys.stdout.write('\r')
                pc=((float(k)/float(np.count_nonzero(bin_mask)-1))*100)
                sys.stdout.write("\t"+str(k)+" / "+str((np.count_nonzero(bin_mask)-1)-1)+" ["+str("{0:.2f}".format(pc))+"%] complete")
                sys.stdout.flush()
                
                if(cis_mean != 0):
                    print(str(chr_id)+"\t"+str(bin_start)+"\t"+str(bin_end)+"\t"+str("{0:.{1}f}".format(cis_mean,precision)),end="\n",file=meancis_fh)
                else:
                    print(str(chr_id)+"\t"+str(bin_start)+"\t"+str(bin_end)+"\t"+str(np.nan),end="\n",file=meancis_fh)
                
                k+=1
                
    meancis_fh.close()
    print("\ndone")
    
    print("")

def getSmallUniqueString():  
    tmp_uniq=str(uuid.uuid4())
    tmp_uniq=tmp_uniq.split('-')[-1]
    return(tmp_uniq)
    
def bin2header(bin,genome,chrs,index=getSmallUniqueString()):
    #name|assembly|chr:start-end
    header=str(index)+'|'+genome+'|'+str(chrs[bin[0]])+':'+str(bin[1])+'-'+str(bin[2])
    return(header)

def deGroupChr(chr_id):
    return(chr_id.split('-')[0])
    
def deGroupHeader(header,extractBy="liteChr",index=getSmallUniqueString()):
    m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
    if m==None:
        sys.exit('error: incorrect input format!')
                
    bin_id,genome,chr_id,bin_start,bin_end=m.groups()
    chr_id=chr_id.split('-')[0]

    header=str(bin_id)+'|'+genome+'|'+str(chr_id)+':'+str(bin_start)+'-'+str(bin_end)
    
    return(header)
    
def flip_intervals(a,b):
    """flip intervals, to ensure a < b
    """
    
    return(b,a)
    
def is_overlap(a, b):
    """test to for overlap between two intervals.
    """
    
    if(a[0] > a[1]):
        sys.exit('\nerror: incorrectly formated interval! start '+str(a[0])+' > end '+str(a[1])+'!\n\t'+str(a)+' '+str(b)+'\n')
    if(b[0] > b[1]):
        sys.exit('\nerror: incorrectly formated interval! start '+str(b[0])+' > end '+str(b[1])+'!\n\t'+str(a)+' '+str(b)+'\n')
    
    if a[0] < b[0] and a[1] > b[1]:
        return((b[1]-b[0])+1)
    
    if b[0] < a[0] and b[1] > a[1]:   
        return((a[1]-a[0])+1)
        
    if b[0] < a[0]:
        a,b=flip_intervals(a,b)
           
    return max(0, ( min(a[1],b[1]) - max(a[0],b[0]) ) ) 
    
if __name__=="__main__":
    main()

   