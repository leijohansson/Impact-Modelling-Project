# -*- coding: utf-8 -*-
"""
Code to unzip ERA and S2S data
"""

import os 
import zipfile

if __name__ == '__main__':
    #unzipping ERA5 data
    zipdir = 'ZipData'
    filenames = os.listdir(zipdir)
    
    if 'dataall' not in os.listdir():
        os.mkdir('dataall')
            
    for file in filenames:
        with zipfile.ZipFile(f'{zipdir}/{file}', 'r') as zip_ref:
            zip_ref.extractall(f'./dataall')
    
    #unzipping S2S data        
    for name in ['ECMF', 'LFPW']:
        zipdir = f'S2SData/{name}/Zip'
        subdir = f'S2SData/{name}'
        filenames = os.listdir(zipdir)
        if 'Data' not in os.listdir(subdir):
            os.mkdir(f'{subdir}/Data')
        for file in filenames:
            with zipfile.ZipFile(f'{zipdir}/{file}', 'r') as zip_ref:
                zip_ref.extractall(f'{subdir}/Data')

