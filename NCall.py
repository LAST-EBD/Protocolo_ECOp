
import os, tarfile, time

os.chdir('/home/diego/Documentos/GitHub/Protocolo_ECOp')

from NProtocolo import NLandsat
from NProductos import Product

sats = {'LT05': 'l5tm', 'LE07': 'l7etm', 'LC08': 'l8oli'}
base_path = '/home/diego/ECOprotocolo'

def run_ECOp(scene_tar, data_tar):
   
    
    t0 = time.time()
    
    sat = os.path.split(scene_tar)[1].split('_')[0]
    serie = os.path.split(scene_tar)[1].split('_')[1]
    scene = os.path.split(scene_tar)[1].split('_')[3] + sats[sat] + '202_34'
    
    if sat in sats.keys() or serie == 'L1TP':
                
        print(sat, serie, scene)
        
        #Untar scene in /ori
        print('uncompressing scene')
        ecop = os.path.join(base_path, os.path.join('ori', scene))
        os.makedirs(ecop, exist_ok=True)
        os.chdir(ecop)
        tar = tarfile.open(scene_tar)
        tar.extractall()
        tar.close()
        
        #Untar data in /data
        print('uncompressing data')
        os.chdir(base_path)
        tar = tarfile.open(data_tar)
        tar.extractall()
        tar.close()
        
        #Create rad, nor and temp path before run the process
        print('creating paths')
        os.makedirs(os.path.join(base_path, 'rad'), exist_ok=True)
        os.makedirs(os.path.join(base_path, 'nor'), exist_ok=True)
        os.makedirs(os.path.join(base_path, 'temp'), exist_ok=True)
        os.makedirs(os.path.join(base_path, 'pro'), exist_ok=True)
        
        #Run the process
        print('starting the process')
        escena = NLandsat(ecop)
        escena.run()
        producto = Product(os.path.join(escena.nor, escena.escena))
        producto.ndvi()
        producto.turbidity(producto.flood())
        
        print('process finished in', time.time() - t0, 'segundos')

        
    else: 
    
       print('This process only works with L1TP Landsats 5-TM, 7-ETM and 8-OLI')
        
        
    #with zipfile.ZipFile(scene_zip, 'r') as inp_ref:
        #inp_ref.extractall(os.path.join(ecop,'/ori/'))
        
run_ECOp('/home/diego/Descargas/LE07_L1TP_202034_20180730_20180825_01_T1.tar.gz', '/home/diego/protocolo/data.tar.gz')

