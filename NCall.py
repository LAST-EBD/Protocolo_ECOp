
import os
os.chdir('/home/diego/Documentos/GitHub/Protocolo_ECOp')


from NProtocolo import NLandsat
from NProductos import Product

landsat = '/home/diego/protocolo/ori/20050920l5tm202_34' #Here it must be the input 1 of the VLab
escena = NLandsat(landsat)
escena.run()
producto = Product(os.path.join(escena.nor, escena.escena))
producto.ndvi()
#producto.flood()
producto.turbidity(producto.flood())
