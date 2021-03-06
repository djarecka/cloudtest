from cffi import FFI
import numpy as np
import pdb

ffi = FFI()

keepalive = []

def zero_cffiarray(nz):
     tab_np = np.zeros((nz))
     keepalive.append(tab_np)
     tab  = ffi.cast("double*", tab_np.__array_interface__['data'][0])
     return tab

def morrison_2momNc_simplewarm(nx, ny, nz, qc_np, qr_np, nc_np, nr_np,
                               temp_np, qv_np, press_np, dz_np, dt_in):
                               
     ffi.cdef("void c_init();", override=True)
     
     ffi.cdef("void c_morr(double QC3DTEN[], double QI3DTEN[], double QNI3DTEN[], double QR3DTEN[], double NC3DTEN[], double NI3DTEN[], double NS3DTEN[], double NR3DTEN[], double QC3D[], double QI3D[], double QNI3D[], double QR3D[], double NI3D[], double NS3D[], double NR3D[], double NC3D[], double T3DTEN[], double QV3DTEN[], double T3D[], double QV3D[], double PRES[], double DZQ[], double W3D[], double WVAR[], double PRECRT, double SNOWRT, double EFFC[], double EFFI[], double EFFS[], double EFFR[], double DT, int IMS, int IME, int JMS, int JME, int KMS, int KME, int ITS, int ITE, int JTS, int JTE, int KTS, int KTE, double QG3DTEN[], double NG3DTEN[], double QG3D[], double NG3D[], double EFFG[], double qrcu1d[], double qscu1d[], double qicu1d[], double QGSTEN[], double QRSTEN[], double QISTEN[], double QNISTEN[], double QCSTEN[]);", override=True)

     lib = ffi.dlopen('libmorrison_2momNc.so')

     QC3DTEN, QI3DTEN, QNI3DTEN, QR3DTEN = (zero_cffiarray(nz) for i in range(4))
     NC3DTEN, NI3DTEN, NS3DTEN, NR3DTEN  = (zero_cffiarray(nz) for i in range(4))
     QC3D = ffi.cast("double*", qc_np.__array_interface__['data'][0])
     QI3D, QNI3D =  (zero_cffiarray(nz) for i in range(2))
     QR3D = ffi.cast("double*", qr_np.__array_interface__['data'][0])
     NI3D, NS3D =  (zero_cffiarray(nz) for i in range(2))
     NR3D = ffi.cast("double*", nr_np.__array_interface__['data'][0])
     NC3D = ffi.cast("double*", nc_np.__array_interface__['data'][0])
     T3DTEN, QV3DTEN = (zero_cffiarray(nz) for i in range(2))
     T3D = ffi.cast("double*", temp_np.__array_interface__['data'][0])
     QV3D = ffi.cast("double*", qv_np.__array_interface__['data'][0])
     PRES = ffi.cast("double*", press_np.__array_interface__['data'][0])
     DZQ = ffi.cast("double*", dz_np.__array_interface__['data'][0])
     W3D, WVAR = (zero_cffiarray(nz) for i in range(2))
     PRECRT = 0.
     SNOWRT = 0.
     EFFC, EFFI, EFFS, EFFR = (zero_cffiarray(nz) for i in range(4))
     DT = dt_in
     [IMS, IME, ITS, ITE] = [1, nx] * 2
     [JMS, JME, JTS, JTE] = [1, ny] * 2
     [KMS, KME, KTS, KTE] = [1, nz] * 2
     QG3DTEN, NG3DTEN, QG3D, NG3D, EFFG = (zero_cffiarray(nz) for i in range(5))
     qrcu1d, qscu1d, qicu1d = (zero_cffiarray(nz) for i in range(3))
     QGSTEN, QRSTEN, QISTEN, QNISTEN, QCSTEN = (zero_cffiarray(nz) for i in range(5))

     print "cffi", QGSTEN, QRSTEN, QISTEN, QNISTEN, QCSTEN 
     lib.c_init()
     lib.c_morr(QC3DTEN, QI3DTEN, QNI3DTEN, QR3DTEN,       
                                 NC3DTEN, NI3DTEN, NS3DTEN, NR3DTEN,              
                                 QC3D, QI3D, QNI3D, QR3D, NI3D, NS3D, NR3D, NC3D, 
                                 T3DTEN, QV3DTEN, T3D, QV3D,                      
                                 PRES, DZQ, W3D, WVAR, PRECRT, SNOWRT,            
                                 EFFC, EFFI, EFFS, EFFR,                          
                                 DT,                                              
                                 IMS, IME, JMS, JME, KMS, KME,                    
                                 ITS, ITE, JTS, JTE, KTS, KTE,                     
                                 QG3DTEN, NG3DTEN, QG3D, NG3D, EFFG,              
                                 qrcu1d, qscu1d, qicu1d,                          
                                 QGSTEN, QRSTEN, QISTEN, QNISTEN, QCSTEN           
                                 )
     
