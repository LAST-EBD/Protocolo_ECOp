
# coding: utf-8
#%matplotlib inline

import os, shutil, re, time, subprocess, pandas, rasterio, sys, stat 
import numpy as np
import gdal, gdalconst
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
from datetime import datetime, date
from scipy import ndimage
from scipy.stats import linregress
from pymasker import LandsatMasker
from pymasker import LandsatConfidence



class NLandsat(object):
    
       
    
    '''This python class is made to be an automated alternative to Landsat Normalization Process of 
    Remote Sensing and GIS Laboratory (Donana Biological Station (CSIC)). 
    This process consists of three steps: 
    1. Geometric Correction, which basically apply a common extent to all the Landsat Scenes (keeping crs WGS84 h29)
    2. Radiometric Correction, Getting Surface Reflectance from DNs. This process is made by applying the same 
    algoritm that MIRAMON Software has implemented in Corrad Module (Dark Object Subtraction, Kl + DTM)
    3. Normalization, this step applies a normalization of the reflectivity values through a linear regression over
    some predefined Pseudo Invariant Areas (PIAs)
    
    Originally, this process uses Fmask, but to avoid Matlab in the VLab, we use pymkaser library, and get the water
    mask with the same thresholds that Fmask uses(pymasker uses Landsat BQA band, which doesn't have a water mask)
    
    The process requires some folders on the same level (/ ori, / rad, / nor and / data).

    In the data folder there must be some files required for the process:
    
        1) Landsat bands of referency scene (Landsat 7, 20020718l7etm202_34)
        2) Shapefile with the limits of Donana National Park
        3) Digital Terrain Model of the projected Extent
        4) PIAs raster mask'''
    
    
    
    def __init__(self, ruta, umbral=50, hist=1000):
        
        
        '''We made the class instance with the current scene to process, we just need to include the path to ori,
        and all the other variables are getting from there. Default parameters are threshold confidence for Fmask 
        and number of bins for the histogram in the kl step'''
        
    
        self.ruta_escena = ruta
        self.ori = os.path.split(ruta)[0]
        self.escena = os.path.split(ruta)[1]
        self.raiz = os.path.split(self.ori)[0]
        self.rad = os.path.join(self.raiz, 'rad')
        self.nor = os.path.join(self.raiz, 'nor')
        self.data = os.path.join(self.raiz, 'data')
        self.temp = os.path.join(self.raiz, 'temp')
        self.umbral = umbral
        self.hist = hist
        if 'l7etm' in self.escena:
            self.sat = 'L7'
        elif 'l8oli' in self.escena:
            self.sat = 'L8'
        elif 'l5tm' in self.escena:
            self.sat = 'L5'
        else:
            print(' No reconozco el satelite')
        print(self.sat, self.escena)
        self.kl = {}
        self.equilibrado = os.path.join(self.data, 'Equilibrada.tif')
        self.noequilibrado = os.path.join(self.data, 'NoEquilibrada.tif')
        self.parametrosnor = {}
        self.iter = 1
        self.cloud_mask = 'None'
        
        self.mtl = {}
        for i in os.listdir(self.ruta_escena):
            #print i
            if i.endswith('MTL.txt'):
                mtl = os.path.join(self.ruta_escena,i)
                
                f = open(mtl, 'r') 

                #Dict
                for line in f.readlines(): 
                    if "=" in line: 
                        l = line.split("=") 
                        self.mtl[l[0].strip()] = l[1].strip()             
                        
        #Quicklook download 
        #print('We are going to save the quicklook')
        #self.quicklook = os.path.join(self.ruta_escena, self.mtl['LANDSAT_SCENE_ID'] + '.jpg')

        #if self.sat == 'L7':
            #sensor = 'etm'
        #elif self.sat == 'L5' or self.sat == 'L4':
            #sensor = 'tm'

        #qcklk = open(self.quicklook,'wb')
        #if self.sat == 'L8':
            #s = 'https://earthexplorer.usgs.gov/browse/landsat_8_c1/{}/202/034/{}.jpg'.format(self.escena[:4], self.mtl['LANDSAT_PRODUCT_ID'][1:-1])
            #print(s)
        #else:
            #s = 'https://earthexplorer.usgs.gov/browse/landsat_{}_c1/{}/202/034/{}.jpg'.format(sensor, self.escena[:4], self.mtl['LANDSAT_PRODUCT_ID'][1:-1])
            #print(s)
            
        #We save the miage quicklook to /ori
        #u2=urlopen(s)
        #junk=u2.read()
        #qcklk.write(junk)
        #qcklk.close()
        
        
        
    def get_water(self):
    
        '''This method generate the Cloud, shadow and water mask alternative to Fmask'''
        
        #We are saving the watermask needed to calculate kl values in temp (same values that are used in Fmask)
        print('starting get water')
        outfile = os.path.join(self.temp, self.escena + '_watermask.tif')
        
        for i in os.listdir(self.ruta_escena):
            
            
            
            if i.endswith('.TIF') and not 'BQA' in i:
                
                banda = i.split('_')[-1].split('.')[0]
                
                if banda not in ['B10', 'B11', 'B6', '1', '2']:
                    REFMULT = float(self.mtl['REFLECTANCE_MULT_BAND_'+ banda[1:]])
                    REFADD = float(self.mtl['REFLECTANCE_ADD_BAND_'+ banda[1:]])

                if self.sat == 'L8' and banda == 'B4':
                    with rasterio.open(os.path.join(self.ruta_escena, i)) as src:
                        B4 = src.read()
                        RED = REFMULT * B4 + REFADD

                elif self.sat == 'L8' and banda == 'B5':
                    with rasterio.open(os.path.join(self.ruta_escena, i)) as src:
                        B5 = src.read()
                        NIR = REFMULT * B5 + REFADD 

                elif self.sat != 'L8' and banda == 'B3':
                    with rasterio.open(os.path.join(self.ruta_escena, i)) as src:
                        B3 = src.read()
                        RED = REFMULT * B3 + REFADD

                elif self.sat != 'L8' and banda == 'B4':
                    with rasterio.open(os.path.join(self.ruta_escena, i)) as src:
                        B4 = src.read()
                        NIR = REFMULT * B4 + REFADD

            else: continue

        num = NIR.astype(float)-RED.astype(float)
        den = NIR+RED
        NDVI = np.true_divide(num, den)
        #Hay que hacer un ndvi y un toa del nir temporal!
        WATER = np.logical_or(np.logical_and(NDVI < 0.01, NIR < 0.11), np.logical_and(NDVI < 0.1, NIR < 0.05))
        water = np.where((num == 0), 0, WATER)
        profile = src.meta
        profile.update(nodata=0)
        profile.update(dtype=rasterio.ubyte)

        with rasterio.open(outfile, 'w', **profile) as dst:
            dst.write(water.astype(rasterio.ubyte))
                    
            
    def fmask(self):
            
            '''-----\n
            This method applies Fmask algorithm or its equivalent with the BQA band (default for Ecopotential)'''
            
            os.chdir(self.ruta_escena)
                
            print('Starting Cloud Mask')
            
            try:
                
                print('Starting Fmask')
                t = time.time()
                #Last value is the cloud confidence value 
                a = os.system('removethisisFmaskinstalled/usr/GERS/Fmask_4_0/application/run_Fmask_4_0.sh /usr/local/MATLAB/MATLAB_Runtime/v93 3 3 1 {}'.format(self.umbral))
                a
                if a == 0:
                    self.cloud_mask = 'Fmask'
                    print('Cloud Mask (Fmask) generated in ' + str(t-time.time()) + ' seconds')
                    
                else:
                    
                    print('Starting Cloud Mask with BQA band')
                    outfile = os.path.join(os.path.join(self.ruta_escena, self.escena + '_Fmask.tif'))
                    
                    for i in os.listdir(self.ruta_escena):
                        if i.endswith('BQA.TIF'):
                            bqa = os.path.join(self.ruta_escena, i)
                            masker = LandsatMasker(bqa, collection=1)
                            conf = LandsatConfidence.high

                            cloud = masker.get_cloud_mask(conf)
                            cirrus = masker.get_cirrus_mask(conf)
                            shadow = masker.get_cloud_shadow_mask(conf)

                            MASCARA = cloud + cirrus + shadow                    
                            masker.save_tif(MASCARA, outfile)
                            
                    self.cloud_mask = 'BQA'
                    print('Cloud Mask (BQA) generated in ' + str(t-time.time()) + ' seconds')
                    print('Calling get_water')
                    self.get_water()
                                           
            except Exception as e:
                
                print("Unexpected error:", type(e), e)
                
                
    def get_cloud_pn(self):
        
        '''-----\n
        This method clips the fmask with the shp of the National Park, to obtain the cloud coverage 
        in Donana National Park in the next step'''
        
        shape = os.path.join(self.data, 'Limites_PN_Donana.shp')
        crop = "-crop_to_cutline"
                    
        for i in os.listdir(self.ruta_escena):
            if i.endswith('Fmask4.tif') or i.endswith('Fmask.tif'):
                cloud = os.path.join(self.ruta_escena, i)
                
        #we use Gdalwarp for the masks
        cmd = ["gdalwarp", "-dstnodata" , "0" , "-cutline", ]
        path_masks = os.path.join(self.ruta_escena, 'masks')
        os.makedirs(path_masks, exist_ok=True)

        
        salida = os.path.join(path_masks, 'cloud_PN.TIF')
        cmd.insert(4, shape)
        cmd.insert(5, crop)
        cmd.insert(6, cloud)
        cmd.insert(7, salida)

        proc = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout,stderr=proc.communicate()
        exit_code=proc.wait()

        if exit_code: 
            raise RuntimeError(stderr)
                          
        
        ds = gdal.Open(salida)
        cloud = np.array(ds.GetRasterBand(1).ReadAsArray())
        
        if self.cloud_mask == 'BQA':
            mask = (cloud >= 1)
        elif self.cloud_mask == 'Fmask': 
            mask = (cloud == 2) | (cloud == 4)
        else:
            print('There is no Cloud Mask... Bad thing')
        
        cloud_msk = cloud[mask]
        print(cloud_msk.size)
        clouds = float(cloud_msk.size*900)
        PN = 534158729.313 
        pn_cover = round((clouds/PN)*100, 2)
        ds = None
        cloud = None
        cloud_msk = None
        clouds = None
               
            
        print("The percentage of clouds in the National Park is " + str(pn_cover))
        
        
    def remove_masks(self):
        
        '''-----\n
        This method removes the folder in which we have saved the masks used to obtain the kl and 
        the cloud cover in the National Park'''
        
        path_masks = os.path.join(self.ruta_escena, 'masks')
        for i in os.listdir(path_masks):
            
            name = os.path.join(path_masks, i)
            os.chmod(name, stat.S_IWRITE)
            os.remove(name)

        shutil.rmtree(path_masks)
        
        
    def projwin(self):
        
        '''This method provided the same extent for every scene'''
        
        path_rad = os.path.join(self.rad, self.escena)
        os.makedirs(path_rad, exist_ok=True)
        
        if self.sat == 'L8':
            
            for i in os.listdir(self.ruta_escena):
                if re.search('B[2-7]', i):

                    ins = os.path.join(self.ruta_escena, i)
                    out = os.path.join(path_rad, i)

                    cmd = "gdal_translate -projwin  623385.0 4266315.0 867615.0 4034685.0 {} {}".format(ins, out)
                    print(cmd)
                    os.system(cmd)
                    
                elif re.search('Fmask4', i):

                    ins = os.path.join(self.ruta_escena, i)
                    out = os.path.join(path_rad, i)

                    cmd = "gdal_translate -projwin  623385.0 4266315.0 867615.0 4034685.0 -a_nodata 255 {} {}".format(ins, out)
                    print(cmd)
                    os.system(cmd)
                    
                elif re.search('Fmask', i):

                    ins = os.path.join(self.ruta_escena, i)
                    out = os.path.join(path_rad, i)

                    cmd = "gdal_translate -projwin  623385.0 4266315.0 867615.0 4034685.0 -a_nodata 255 {} {}".format(ins, out)
                    print(cmd)
                    os.system(cmd)
                    
                    #In this case we also need to proj water_mask
                    wmask = [i for i in os.listdir(self.temp) if i.endswith('watermask.tif')]
                    ins = os.path.join(self.temp, wmask[0])
                    out = os.path.join(self.temp, wmask[0][:-4] + '_proj.tif')
                    
                    cmd = "gdal_translate -projwin  623385.0 4266315.0 867615.0 4034685.0 -a_nodata 255 {} {}".format(ins, out)
                    print(cmd)
                    os.system(cmd)                    

                    
                else: continue
                    
        else:
            
            for i in os.listdir(self.ruta_escena):
                if (re.search('B[1-7]', i) and not 'B6' in i) or re.search('Fmask4', i):

                    ins = os.path.join(self.ruta_escena, i)
                    out = os.path.join(path_rad, i)

                    cmd = "gdal_translate -projwin  623385.0 4266315.0 867615.0 4034685.0 {} {}".format(ins, out)
                    print(cmd)
                    os.system(cmd)
                    
                elif re.search('Fmask', i):

                    ins = os.path.join(self.ruta_escena, i)
                    out = os.path.join(path_rad, i)

                    cmd = "gdal_translate -projwin  623385.0 4266315.0 867615.0 4034685.0 -a_nodata 255 {} {}".format(ins, out)
                    print(cmd)
                    os.system(cmd)
                    
                    #In this case we also need to proj water_mask
                    wmask = [i for i in os.listdir(self.temp) if i.endswith('watermask.tif')]
                    ins = os.path.join(self.temp, wmask[0])
                    out = os.path.join(self.temp, wmask[0][:-4] + '_proj.tif')
                    
                    cmd = "gdal_translate -projwin  623385.0 4266315.0 867615.0 4034685.0 -a_nodata 255 {} {}".format(ins, out)
                    print(cmd)
                    os.system(cmd)                    

                    
                else: continue
                
        print('Projections to the common extent finished')
                
    def get_kl_csw(self):
        
        '''This method obtains the Kl for each band. It does so by looking for the minimum values within 
        the zones classified as water and orographic shadow, as long as the orographic shadow is not covered 
        by clouds or cloud shadow. The quality of the mask is very important, that is why the scenes that 
        can not be done with Fmask would have to check the value of kl.'''
    
        #We started deleting files from temp, but this are not needed for the VLab or we can do it at the end of the process
        
        #for i in os.listdir(self.temp):
            #arz = os.path.join(self.temp, i)
            #os.remove(arz)

        
        t = time.time()
             
        dtm = os.path.join(self.data, 'dtm_extent_l8.tif') #WGS84 h29 (Path 202, Row 034; Donana National Park Scene)
        azimuth = self.mtl['SUN_AZIMUTH']
        elevation = self.mtl['SUN_ELEVATION']
        
        #Hillshade generation
        salida = os.path.join(self.temp, 'hillshade.img')
        cmd = ["gdaldem", "hillshade", "-az", "-alt", "-of", "GTIFF"]
        cmd.append(dtm)
        cmd.append(salida)
        cmd.insert(3, str(azimuth))
        cmd.insert(5, str(elevation))
        proc = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        stdout,stderr=proc.communicate()
        exit_code=proc.wait()

        if exit_code: 
            raise RuntimeError(stderr)
        else:
            print(stdout)
            print('Hillshade generado')
        
        ruta = os.path.join(self.rad, self.escena)
        print('ruta:', ruta)
        
        for i in os.listdir(ruta):
            
            if 'Fmask' in i: 
                rs = os.path.join(ruta, i)
                Fmask = gdal.Open(rs).ReadAsArray()
                print('Fmask: min, max, size: ', Fmask.min(), Fmask.max(), Fmask.size)
        
        for i in os.listdir(self.temp):
            
            if i.endswith('shade.img'):
                
                hill = os.path.join(self.temp, i)
                Hillshade = gdal.Open(hill).ReadAsArray()
                print('Hillshade: min, max, size: ', Hillshade.min(), Hillshade.max(), Hillshade.size)       

        #We want the pixels of each band that are within the water value (1) or with any category defined (0) 
        #but in orographic shadow, so we will only look for kl values within water or shadows. Avoiding problems
        #with clouds and clouds shadows. We also use an "inner_buffer" to avoid problems with edges pixels in the 
        #NoData value (this could be a problem in Landsat 5 and Landsat 7)
        
        
        #Now we apply the mask and make the histograms
        if self.sat == 'L8':
            
            bandas = ['B2', 'B3', 'B4','B5', 'B6', 'B7']
            
        else:
            
            bandas = ['B1', 'B2', 'B3', 'B4','B5', 'B7']
            
            
        lista_kl = []
        for i in os.listdir(ruta):
            banda = i[-6:-4]
            if banda in bandas:
                raster = os.path.join(ruta, i)
                print('Raster:', raster)
                #data = gdal.Open(raster).ReadAsArray()
                data = rasterio.open(raster).read()
                #anadimos la distincion entre Fmask y BQA
                if self.cloud_mask == 'Fmask':
                    
                    if self.sat == 'L8':
                        print('using Fmask with Landsat 8')
                        data2 = data[(((Fmask==1) | (((Fmask==0) & (data != 0)) & (Hillshade<(np.percentile(Hillshade, 20))))))]
                    else:
                        print('using Fmask with Landsat', self.sat)
                        #Abrimos la mascara para evitar problemas con los valores bajos en los bordes de Fmask
                        inner = os.path.join(self.data, 'intern_buffer.tif')
                        Inner = gdal.Open(inner).ReadAsArray()
                        data2 = data[((Inner == 1) & ((Fmask==1) | (((Fmask==0) & (data != 0)) & (Hillshade<(np.percentile(Hillshade, 20))))))]
                
                
                elif self.cloud_mask == 'BQA':
                    
                    #Pick the water mask generated to complete BQA info
                    for i in os.listdir(self.temp):
                        if i.endswith('_proj.tif'):
                            water = os.path.join(self.temp, i)
                    Water = rasterio.open(water).read()    
                    
                    if self.sat == 'L8':
                            print('using BQA with Landsat 8')
                            data2 = data[(data != 0) & (((Fmask==0) & (data != 0)) & (Hillshade<(np.percentile(Hillshade, 20))) | (Water==1))]
                    else:
                        print('using BQA with Landsat', self.sat)
                        #Abrimos la mascara para evitar problemas con los valores bajos en los bordes de Fmask
                        inner = os.path.join(self.data, 'intern_buffer.tif')
                        Inner = gdal.Open(inner).ReadAsArray()
                        data2 = data[(data != 0) & (Inner == 1) & (((Fmask==0) & (data != 0)) & (Hillshade<(np.percentile(Hillshade, 20))) | (Water==1))]
                                
                                
                        #mask = np.copy(data)
                        #mask[mask != 0] == 1
                        #print('Ahora viene el erode')
                        #erode = ndimage.grey_erosion(data, size=(5,5,1))
                        #print(type(erode), erode.size, erode.shape)
                        
                        #data2 = data[((erode != 0) & ((Fmask==1) | (((Fmask==0) & (data != 0)) & (Hillshade<(np.percentile(Hillshade, 20))))))]
                    
                    #Aqui iria una alternativa a Fmask si fallara

                print('data 2 obtained')

                lista = sorted(data2.tolist())
                print(sorted(lista)[:10])
                self.kl[banda] = data2.min() #np.mean(lista10)data2.min()
                data3 = data2[:self.hist]

                print('data3: ', data3.min(), data3.max())

                df = pandas.DataFrame(lista[:10000])
                plt.figure(); df.hist(figsize=(10,8), bins = 20, cumulative=False, color="Red"); 
                plt.title(self.escena + '_gr_' + banda, fontsize = 18)
                plt.xlabel("Pixel Value", fontsize=16)  
                plt.ylabel("Count", fontsize=16)
                path_rad = os.path.join(self.rad, self.escena)
                os.makedirs(path_rad, exist_ok=True)
                name = os.path.join(path_rad, self.escena + '_gr_'+ banda.lower() + '.png')
                plt.savefig(name)

        plt.close('all')
        
        print('Histograms generated')
        print('kl_values:', sorted(self.kl.items()))
            
                       
    def get_radiance(self):
                  
        '''-----\n
        This method generates the radiance values of the image, calculating them from *MTL.txt coefficients'''   
                  
        path_rad = os.path.join(self.rad, self.escena)
        
        #if self.sat == 'L8':
        
        for i in os.listdir(path_rad):

            if re.search('B[1-7].TIF', i):

                banda = i.split('_')[-1].split('.')[0]
                print(banda)


                RADMULT = float(self.mtl['RADIANCE_MULT_BAND_'+ banda[1:]])
                RADADD = float(self.mtl['RADIANCE_ADD_BAND_'+ banda[1:]])


                with rasterio.open(os.path.join(path_rad, i)) as src:

                    B = src.read()
                    RAD = RADMULT * B + RADADD

                    profile = src.meta
                    profile.update(dtype=rasterio.float32)

                    outfile = os.path.join(path_rad, banda + '_rad.tif')
                    print(outfile)                

                    with rasterio.open(outfile, 'w', **profile) as dst:
                        dst.write(RAD.astype(rasterio.float32))
                            
                            
    def corrad(self):
                  
        '''-----\n
        This method generates the equivalent to the Radiometric Correction performed by MiraMon, paper:
        Pons X, Solé-Sugrañes L (1994) "A Simple Radiometric Correction Model to Improve Automatic Mapping of 
        Vegetation from Multispectral Satellite Data." Remote Sensing of Environment, 48:191-204.'''      
        
        #Parameters obtained from MIRAMON file "m_atmos_rad (2014 version)"
        bandnames = {'B1': 'blue', 'B2': 'green', 'B3': 'red', 'B4': 'nir', 'B5': 'swir1', 'B7': 'swir2'}
        essunl8 = {'B2':2081.34, 'B3': 1776.50, 'B4': 1552.82, 'B5': 1113.74, 'B6': 171.061, 'B7': 105.38}
        essunl7 = {'B1':1997, 'B2': 1812, 'B3': 1533, 'B4': 1039, 'B5': 230.8, 'B7': 84.9}
        taul8 = {'B2': 0.4, 'B3': 0.34, 'B4': 0.29, 'B5': 0.21, 'B6': 0.11, 'B7': 0.08}
        taul7 = {'B1': 0.5, 'B2': 0.3, 'B3': 0.25, 'B4': 0.2, 'B5': 0.125, 'B7': 0.075}

        if self.sat == 'L8':
            essun, tau = essunl8, taul8
        else:
            essun, tau = essunl7, taul7
            
            
        path_rad = os.path.join(self.rad, self.escena)
        
        for i in os.listdir(path_rad):
            
            if i.endswith('_rad.tif'):
                banda = i.split('_')[0]
                
                if banda in tau.keys():

                    RADMULT = float(self.mtl['RADIANCE_MULT_BAND_'+ banda[1:]])
                    RADADD = float(self.mtl['RADIANCE_ADD_BAND_'+ banda[1:]])
                    print(banda)
                    with rasterio.open(os.path.join(path_rad, i)) as src:
                        B = src.read()

                    klradiancia = (RADMULT * self.kl[banda]) +  RADADD
                    print(klradiancia)
                    radkl = B - klradiancia
                    NUM = np.pi * radkl * np.power(float(self.mtl['EARTH_SUN_DISTANCE']), 2) 

                    t1 = np.power(np.e, float(tau[banda]) / np.cos(np.deg2rad(float(self.mtl['SUN_AZIMUTH'])))) #FALTA COMPROBAR T1 Y T2 d['SUN_AZIMUTH']
                    t2 = np.power(np.e, float(tau[banda]))
                    DEN = np.cos(np.deg2rad(90 - float(self.mtl['SUN_ELEVATION']))) * essun[banda] * t1 * t2 

                    SR = np.divide(NUM, DEN)
                    SR = np.around(SR*10000)
                    SR = np.where(SR>10000, 0, SR)
                    SR = np.where(SR<0, 0, SR)

                    profile = src.meta
                    profile.update(dtype=rasterio.uint16)

                    outfile = os.path.join(path_rad, self.escena + '_gr2_' + banda + '.tif')
                    print(outfile)                

                    with rasterio.open(outfile, 'w', **profile) as dst:
                        dst.write(SR.astype(rasterio.uint16))

    
    
    def clean_rad(self):
        
        '''-----\n
        This method removes unnecessary files from the rad folder and renames Fmask taking it to /nor'''    
        
        
        path_rad = os.path.join(self.rad, self.escena)
        path_nor = os.path.join(self.nor, self.escena)
        os.makedirs(path_nor, exist_ok=True)

        for i in os.listdir(path_rad):
            if not re.search('^[0-9]', i):
                #if not 'Fmask' in i:
                os.remove(os.path.join(path_rad, i))
            elif 'Fmask' in i:
                print('moving Fmask to nor')
                old = os.path.join(path_rad, i)
                new = os.path.join(path_nor, self.escena + '_Fmask.tif')
                os.rename(old, new)
                    
                    
    
    def normalize(self):
        
        '''-----\n
        This method controls the flow of normalization, if the coefficients are not obtained 
        (R> 0.85 and N_Pixels> = 10,  is going to the next level, until you get those values 
        or until you reach the last step)'''
        
        path_rad = os.path.join(self.rad, self.escena)
                
        bandasl8 = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7']
        bandasl7 = ['B1', 'B2', 'B3', 'B4', 'B5', 'B7']
        
        if self.sat == 'L8':
            print('landsat 8\n')
            lstbandas = bandasl8
        else:
            print('landsat', self.sat, '\n')
            lstbandas = bandasl7
        
        #Start the normalization process
        for i in os.listdir(path_rad):
                    
            if re.search('B[1-7].tif$', i):
                
                banda = os.path.join(path_rad, i)
                banda_num = i[-6:-4]
                
                print(banda, 'from normalize')
                #First call to nor1
                self.iter = 1
                self.nor1(banda, self.noequilibrado)
                
                #Run over the different thresholds for the normalization 
                if banda_num not in self.parametrosnor.keys():
                    
                    self.iter += 1
                    print('Iteration', self.iter)
                    self.nor1(banda, self.noequilibrado, coef = 2)
                    if banda_num not in self.parametrosnor.keys():
                        self.iter += 1
                        print('Iteration', self.iter)
                        self.nor1(banda, self.equilibrado)
                        if banda_num not in self.parametrosnor.keys():
                            self.iter += 1
                            print('Iteration', self.iter)
                            self.nor1(banda, self.equilibrado, coef = 2)
                            if banda_num not in self.parametrosnor.keys():
                                self.iter += 1
                                print('Iteration', self.iter)
                                self.nor1(banda, self.noequilibrado, coef = 3,)
                                if banda_num not in self.parametrosnor.keys():
                                    self.iter += 1
                                    print('Iteration', self.iter)
                                    self.nor1(banda, self.equilibrado, coef = 3)
                                else:
                                    print('Can not normalize band', banda_num)
                                    
            #coefficients saved to a file, originally saved in a MongoDB, but it seems that we can't do it in the VLab
            path_nor = os.path.join(self.nor, self.escena)
            #os.makedirs(path_nor, exist_ok=True)
            arc = os.path.join(path_nor, 'coeficientes.txt')
            f = open(arc, 'w')
            for i in sorted(self.parametrosnor.items()):
                f.write(str(i)+'\n')
            f.close()  
            
        
                
    def nor1(self, banda, mascara, coef = 1):
        
        '''-----\n
        This method gets the necessary coefficients to carry out the normalization,'''

        print('starting nor1')
        
        #path to reference scene bands
        path_b1 = os.path.join(self.data, '20020718l7etm202_34_ref_B1.tif')
        path_b2 = os.path.join(self.data, '20020718l7etm202_34_ref_B2.tif')
        path_b3 = os.path.join(self.data, '20020718l7etm202_34_ref_B3.tif')
        path_b4 = os.path.join(self.data, '20020718l7etm202_34_ref_B4.tif')
        path_b5 = os.path.join(self.data, '20020718l7etm202_34_ref_B5.tif')
        path_b7 = os.path.join(self.data, '20020718l7etm202_34_ref_B7.tif')
        
        dnorbandasl8 = {'B2': path_b1, 'B3': path_b2, 'B4': path_b3, 'B5': path_b4, 'B6': path_b5, 'B7': path_b7}
        dnorbandasl7 = {'B1': path_b1, 'B2': path_b2, 'B3': path_b3, 'B4': path_b4, 'B5': path_b5, 'B7': path_b7}
        
        if self.sat == 'L8':
            dnorbandas = dnorbandasl8
        else:
            dnorbandas = dnorbandasl7
        
        path_nor = os.path.join(self.nor, self.escena)
        
        
        for i in os.listdir(path_nor):
            if 'Fmask' in i:
                mask_nubes = os.path.join(path_nor, i)
        print('Cloud: ', mask_nubes)
        
        if mascara == self.noequilibrado:
            poly_inv_tipo = os.path.join(self.data, 'NoEquilibrada.tif')
        else:
            poly_inv_tipo = os.path.join(self.data, 'Equilibrada.tif')

        print('mask: ', mascara)
                            
        with rasterio.open(mask_nubes) as nubes:
            CLOUD = nubes.read()
                
        #Open PIAs raster
        with rasterio.open(poly_inv_tipo) as pias:
            PIAS = pias.read()

        banda_num = banda[-6:-4]
        print(banda_num)
        if banda_num in dnorbandas.keys():
            with rasterio.open(banda) as current:
                CURRENT = current.read()
                print('Current: ', banda, 'Shape:', CURRENT.shape)
            #Dict ensures that we cross current band with correct reference band
            with rasterio.open(dnorbandas[banda_num]) as ref:
                REF = ref.read()
                print('Reference: ', dnorbandas[banda_num], 'Shape:', REF.shape)
            
            #Bands of current scene and reference scene like arrays
            #THIS IS FOR ECOPOTENTIAL WITHOUT FMASK, SO CLOUD MASK IS GOING TO BE ALWAYS BQA
            REF2 = REF[((CURRENT != 0) & (PIAS != 0)) & (CLOUD == 0)]
            BANDA2 = CURRENT[((CURRENT != 0) & (PIAS != 0)) & (CLOUD == 0)]
            PIAS2 = PIAS[((CURRENT != 0) & (PIAS != 0)) & (CLOUD == 0)]
            
            #Apply first regression
            First_slope, First_intercept, r_value, p_value, std_err = linregress(BANDA2,REF2)
            print ('\n++++++++++++++++++++++++++++++++++')
            print('slope: '+ str(First_slope), 'intercept:', First_intercept, 'r', r_value, 'N:', PIAS2.size)
            print ('++++++++++++++++++++++++++++++++++\n')
                        
            esperado = BANDA2 * First_slope + First_intercept
            residuo = REF2 - esperado
            
            print('RESIDUO STD:', residuo.std())
            print('RESIDUO STD_DDOF:', residuo.std(ddof=1))
            std = residuo.std() * coef
            print('STD:', std, 'COEF:', coef)
                        
            #calculating remnant of the regression

            mask_current_PIA_NoData_STD = np.ma.masked_where(abs(residuo)>=std, BANDA2)
            mask_ref_PIA_NoData_STD = np.ma.masked_where(abs(residuo)>=std,REF2)
            mask_pias_PIA_NoData_STD = np.ma.masked_where(abs(residuo)>=std,PIAS2)

            current_PIA_NoData_STD = np.ma.compressed(mask_current_PIA_NoData_STD)
            ref_PIA_NoData_STD = np.ma.compressed(mask_ref_PIA_NoData_STD)
            pias_PIA_NoData_STD = np.ma.compressed(mask_pias_PIA_NoData_STD)
                       
            
            #remnant of the regression masked, applying second regression
            slope, intercept, r_value, p_value, std_err = linregress(current_PIA_NoData_STD,ref_PIA_NoData_STD)
            print ('\n++++++++++++++++++++++++++++++++++')
            print ('slope: '+ str(slope), 'intercept:', intercept, 'r', r_value, 'N:', len(ref_PIA_NoData_STD))
            print ('++++++++++++++++++++++++++++++++++\n')
            
            
            #Check pixel count for every PIA type
            values = {}
            values_str = {1: 'Mar', 2: 'Embalses', 3: 'Pinar', 
                          4: 'Urbano-1', 5: 'Urbano-2', 6: 'Aeropuertos', 7: 'Arena', 8: 'Pastizales', 9: 'Mineria'}
            
            print('Vamos a sacar el count de cada zona (dict)')
            for i in range(1,8):

                mask_pia_= np.ma.masked_where(pias_PIA_NoData_STD != i, pias_PIA_NoData_STD)
                PIA = np.ma.compressed(mask_pia_)
                a = PIA.tolist()
                values[values_str[i]] = len(a)
                print('Values_dict:', values)
            
            #PIAs keys to string
            print(banda_num)
            #Regression equation applied
            if r_value > 0.85 and min(values.values()) >= 10:
                self.parametrosnor[banda_num]= {'Parametros':{'slope': slope, 'intercept': intercept, 'std': std,
                        'r': r_value, 'N': len(ref_PIA_NoData_STD), 'iter': self.iter}, 'Tipo_Area': values}
                
                print('nor1 parameters: ', self.parametrosnor)
                print('\nor2 starting with band:', banda[-6:-4], '\n')
                #Hemos calculado la regresion con las bandas recortadas con Rois_extent
                #Ahora vamos a pasar las bandas de rad (completas) para aplicar la ecuacion de regresion
                path_rad = os.path.join(self.rad, self.escena)
                print('Ruta Rad:', path_rad)
                for r in os.listdir(path_rad):
                    if banda[-6:-4] in r and r.endswith('.tif'):
                        print('band:', r)
                        raster = os.path.join(path_rad, r)
                        print('Band to normalize:', raster)
                        self.nor2l8(raster, slope, intercept)
                        print('\nNormalization of ', banda_num, ' finished.\n')
                 
                        fig = plt.figure(figsize=(15,10))
                        ax1 = fig.add_subplot(121)
                        ax2 = fig.add_subplot(122)
                        ax1.set_ylim((0, 10000))
                        ax1.set_xlim((0, 10000))
                        ax2.set_ylim((0, 10000))
                        ax2.set_xlim((0, 10000))

                        sns.regplot(BANDA2, REF2, color='g', ax=ax1,
                         line_kws={'color': 'grey', 'label':"y={0:.5f}x+{1:.5f}".format(First_slope,First_intercept)}).set_title('Regresion PIAs')

                        sns.regplot(current_PIA_NoData_STD, ref_PIA_NoData_STD, color='b', ax=ax2,
                         line_kws={'color': 'grey', 'label':"y={0:.5f}x+{1:.5f}".format(slope,intercept)}).set_title('Regresion PIAs-STD')

                        #Legend
                        ax1.legend()
                        ax2.legend()

                        title_ = os.path.split(banda)[1][:-4] + '. Iter: ' + str(self.iter)
                        fig.suptitle(title_, fontsize=15, weight='bold')
                        
                        plt.savefig(os.path.join(path_nor, os.path.split(banda)[1][:-4])+'.png')
                        plt.show()
                            
            else:
                pass
                                       
                    
    def nor2l8(self, banda, slope, intercept):
    
        '''-----\n
        This method applies the equation of the regression line to each band (when it's possible)'''
        
        
        print('now we are in nor2!')
        path_rad = os.path.join(self.rad, self.escena)
        path_nor = os.path.join(self.nor, self.escena)
        
        banda_num = banda[-6:-4]
        outFile = os.path.join(path_nor, self.escena + '_grn2_' + banda_num + '.tif')
        print('Outfile', outFile)
        
        #NoData Ref
        for i in os.listdir(path_rad):
            
            if 'B5' in i:
                ref = os.path.join(path_rad, i)
        
        with rasterio.open(ref) as src:
            ref_rs = src.read()
        
        with rasterio.open(banda) as src:

            rs = src.read()
            rs = rs*slope+intercept

            nd = (ref_rs == 0)
            min_msk =  (rs < 0)             
            max_msk = (rs>=10001)

            rs[min_msk] = 0
            rs[max_msk] = 10000

            rs = np.around(rs)
            rs[nd] = 0

            profile = src.meta
            profile.update(dtype=rasterio.uint16)

            with rasterio.open(outFile, 'w', **profile) as dst:
                dst.write(rs.astype(rasterio.uint16))
       
    
                      
    def run(self):
        
        t0 = time.time()
        self.fmask()
        self.get_cloud_pn()
        self.remove_masks()
        self.projwin()
        self.get_kl_csw()
        self.get_radiance()
        self.corrad()
        self.clean_rad()
        self.normalize()
        print('Process finsihed in', abs(t0-time.time()), 'seconds')
