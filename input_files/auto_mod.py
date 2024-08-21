import csv
import os
import time
import math
import multiprocessing
import logging
import pandas as pd
from multiprocessing import Process, current_process
from multiprocessing import log_to_stderr,get_logger
from pysigproc import SigprocFile
import argparse

def write_path(filpath,cfile):
    with open(cfile,'r+') as csv_file:
     string=[]
     for line in csv_file:
      string.append(filpath+','+line)
     csv_file.seek(0)
     for i in string:
      csv_file.write(i)
    csv_file.close()

class Modfil(SigprocFile):
  def __init__(self, fp=None, copy_hdr=None, dm=None, tcand=0, width=0, label=-1, snr=0, min_samp=256, device=0,kill_mask=None):
    SigprocFile.__init__(self, fp, copy_hdr)
  
  def getchans(self):
    nchans = getattr(self, 'fch1')
    frequency_chans = float(nchans)
    print(f'number of channels = {nchans}')
    return frequency_chans


if __name__ == '__main__':
  
  parser = argparse.ArgumentParser()
 #logging info
  log_to_stderr()
  logger = get_logger()
  logger.setLevel(logging.INFO)
  
 #optional arguments
  parser.add_argument('--pout', help='Output file directory for candidate h5', type=str,default = '/home/guest/suyash/default_cand/')
  parser.add_argument('--pin', help='input csv file for candmaker', type=str,default='/home/guest/suyash/fetch/cands.csv')
  parser.add_argument('--fpath', help='filterbank path', type=str)
  parser.add_argument('--clean',help='deletes contents of default',default = 0)
  args =parser.parse_args()

  fil=str(args.fpath)
  in_path= str(args.pin)
  out_path= str(args.pout)
  clean = int(args.clean)
  
  if (clean == 1):
    print("cleaning folder")
    os.system(f"find {out_path} -type f -delete")
 
  write_path(fil,in_path)
 
  start_time = time.time()
  names=[]
  csvfile = open(in_path, 'r').readlines()
  print(len(csvfile))
  filindex = 1
  for i in range(len(csvfile)):
   filename = str(out_path)+str(filindex)+ '.csv'
   #names.append(str(filename))
   #open(str(filename) , 'w+').writelines(csvfile[i])
   row_csv = str(csvfile[i]).split(',')
   if len(row_csv) == 5:
    filep= row_csv[0]
    #filep= fil
    snr =  float(row_csv[1])
    start = float(row_csv[2])
    dm = float(row_csv[3])  
    width = float(row_csv[4])
    print(filep)
    print(dm)
    print(width)
    mod = Modfil(filep, snr=snr, width=2, dm=dm, tcand=start) 
    fch1 = float(mod.getchans())
    dec = 1
    if(fch1 < 500):
      if(dm > 1000) and (dm < 1500):
       dec = 2
       width = width -1
       print(f"decimate file {filep}")
       newfilep=(filep.rsplit('/')[-1]).split('.fil')[0]+str('dec')+str('.fil')
       newfilename = str(out_path)+str(newfilep)
       os.system(f"decimate {filep} -c 1 -t {dec} > {newfilename}")
       print(f"new filterbank file is {newfilename}")


       
      if(dm >= 1500):
        dec = 2
        chan = 2
        width = width -1
        newfilep=(filep.rsplit('/')[-1]).split('.fil')[0]+str('dec')+str('.fil')
        newfilename = str(out_path)+str(newfilep)
        os.system(f"decimate {filep} -c {chan} -t {dec} > {newfilename}")
        print(f"decimating file {filep}")
        print(f"new filterbank file is {newfilename}")
    if dec == 1:
     newline = filep + str(",") + str(row_csv[1]) + str(",") + str(row_csv[2]) + str(",") + str(row_csv[3]) + str(",") +str(width)
    else :
     newline = newfilename + str(",") + str(row_csv[1]) + str(",") + str(row_csv[2]) + str(",") + str(row_csv[3]) + str(",") +str(width) 
    names.append(str(filename))
    open(str(filename) , 'w+').writelines(newline)
    filindex += 1



    
  for mod_path in names:
     os.system(f"candmaker.py  --nproc 50 --frequency_size 256 --time_size 256 --cand_param_file {mod_path}  --plot --fout {out_path} ")


 # os.system(f"candmaker.py  --nproc 10 --frequency_size 256 --time_size 256 --cand_param_file {in_path}  --plot --fout {out_path} --src1 {s1} --src2 {s2}")
  #os.system(f"candmaker.py  --nproc 50 --frequency_size 256 --time_size 256 --cand_param_file {in_path}  --plot --fout {out_path} ")
  os.system(f"predict.py --data_dir {out_path} --model a")
  print("displaying results")
  end_time = time.time()
  total_time = end_time - start_time
  logging.info(f'End to End time: {total_time}')
  print (f'End to End time: {total_time}')
  if (os.path.exists(f"{out_path}auto_output_level0")==False): 
   os.system(f"mkdir {out_path}auto_output_level0")
  if (os.path.exists(f"{out_path}auto_output_level1")==False): #check for folder
   os.system(f"mkdir {out_path}auto_output_level1")
  df=pd.read_csv(f'{out_path}results_a.csv')
  for ind in df.index:
    if(df['label'][ind]==1):
        fil=df['candidate'][ind]
        ret=fil.split('/')[-1]
        target=ret.split('.h5')[0]
        img= target+'.png'
        os.system(f"mv {out_path}/{ret} {out_path}auto_output_level1/")
        os.system(f"mv {out_path}/{img} {out_path}auto_output_level1/")
    if(df['label'][ind]==0):
        fil=df['candidate'][ind]
        ret=fil.split('/')[-1]
        target=ret.split('.h5')[0]
        img= target+'.png'
        os.system(f"mv {out_path}/{ret} {out_path}auto_output_level0/")
        os.system(f"mv {out_path}/{img} {out_path}auto_output_level0/")
       # print(img)
       # print(ret)

