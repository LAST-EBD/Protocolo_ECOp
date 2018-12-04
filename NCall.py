
import os, tarfile, time

#os.chdir('/home/diego/Documentos/GitHub/Protocolo_ECOp')
base_path = os.getcwd()

from NProtocolo import NLandsat
from NProductos import Product

sats = {'LT05': 'l5tm', 'LE07': 'l7etm', 'LC08': 'l8oli'}


def run_ECOp(scene_tar, data_tar):
   
    print('There we go!')
    t0 = time.time()
    
    sat = os.path.split(scene_tar)[1].split('_')[0]
    serie = os.path.split(scene_tar)[1].split('_')[1]
    scene = os.path.split(scene_tar)[1].split('_')[3] + sats[sat] + '202_34'
    
    if sat in sats.keys() or serie == 'L1TP':
                
        print(sat, serie, scene)
        
        #Untar scene in /ori
        print('uncompressing scene')
        ecop = os.path.join('ori', scene)
        os.makedirs(ecop, exist_ok=True)
        os.chdir(ecop)
        tar = tarfile.open(scene_tar)
        tar.extractall()
        tar.close()
        
        #Untar data in /data
        print('uncompressing data')
        os.chdir('../../')
        print(os.getcwd())
        tar = tarfile.open(data_tar)
        tar.extractall()
        tar.close()
        
        #Create rad, nor and temp path before run the process
        print('creating paths')
        os.makedirs('rad', exist_ok=True)
        os.makedirs('nor', exist_ok=True)
        os.makedirs('temp', exist_ok=True)
        os.makedirs('pro', exist_ok=True)
        
        #Run the process
        print('starting the process')
        escena = NLandsat(os.path.join(base_path, ecop))
        escena.run()
        producto = Product(os.path.join(escena.nor, escena.escena))
        producto.ndvi()
        producto.turbidity(producto.flood())
        
        #Comnpress the output files (ndvi, flood mask and water turbidity in the output folder)
        print('compressing the products')
        out_compress = tarfile.open(producto.pro_esc + '.tar.gz',  mode='w:gz')
        for i in os.listdir(producto.pro_esc):
            out_compress.add(os.path.join(producto.pro_esc, i))
            
        
        print('process finished in', time.time() - t0, 'segundos')

        
    else: 
    
       print('This process only works with L1TP Landsats 5-TM, 7-ETM and 8-OLI')
        
        
    
#Run the code taking like i/o the ones defined in the VLab/iodescription.json
if __name__ == '__main__':
    run_ECOp(os.path.join(base_path,'Data/LC08_L1TP_202034_20180924_20180929_01_T1.tar.gz'), os.path.join(base_path, 'Data/data.tar.gz'))
