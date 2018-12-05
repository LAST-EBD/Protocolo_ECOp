
import os, tarfile, time

from NProtocolo import NLandsat
from NProductos import Product

sats = {'LT05': 'l5tm', 'LE07': 'l7etm', 'LC08': 'l8oli'}


def run_ECOp(scene_tar, data_tar):
   
    print('There we go!')
    
    os.chdir('../')
    base_path = os.getcwd()
    
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
        
        #We don't need to compress each scene product folder beacuse we are going to compress the whole nor folder later
        
        #Compress the output files (ndvi, flood mask and water turbidity in the output folder)
        #print('compressing the products')
        #out_compress = tarfile.open(producto.pro_esc + '.tar.gz',  mode='w:gz')
        #for i in os.listdir(producto.pro_esc):
            #out_compress.add(os.path.join(producto.pro_esc, i))
            
        
        print('process finished in', time.time() - t0, 'segundos')

        
    else: 
    
       print('This process only works with L1TP Landsats 5-TM, 7-ETM and 8-OLI')
        
           
    
#Run the code taking like i/o the ones defined in the VLab/iodescription.json
if __name__ == '__main__':
    
    t = time.time()
    
    print('Empezamos...')
    os.chdir('Data')
    path = os.getcwd()
    data_tar = os.path.join(path, 'scene.tar.gz')
    tar = tarfile.open(data_tar)
    tar.extractall()
    tar.close()
    
    scenes = os.path.join(path, 'scene')
    for i in os.listdir(scenes):
        try:
            scene = os.path.join(scenes, i)
            print(scene)
            run_ECOp(scene, os.path.join(path, 'data.tar.gz'))
        except Exception as e:
            print(e)
            continue
        
    print('Moving nor to scene_products')
    os.chdir('../../')
    print('We are in', os.getcwd())
    dst = 'Data/scene_products'
    os.rename('pro', dst)
    
    print('Compressing scene_products')
    out_compress = tarfile.open(dst + '.tar.gz',  mode='w:gz')
    for i in os.listdir(dst):
        out_compress.add(os.path.join(dst, i))
    
    print('Process finished in', time.time() - t, 'segundos')
