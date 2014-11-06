module morrison_2mom_wrap

  use iso_c_binding, only: c_int, c_double
  use module_mp_morr_two_moment

  implicit none

contains

  subroutine c_init(hail_opt) bind(c)
    integer(c_int), intent(in), value::  hail_opt 
    integer(c_int) :: n, i
    call MORR_TWO_MOMENT_INIT(hail_opt)
  end subroutine c_init
    

  subroutine c_morr(QC3DTEN, QI3DTEN, QNI3DTEN, QR3DTEN,         &
                    NI3DTEN, NS3DTEN, NR3DTEN,                   &
                    QC3D, QI3D, QNI3D, QR3D, NI3D, NS3D, NR3D,   &
                    T3DTEN, QV3DTEN, T3D, QV3D,                  &
                    PRES, DZQ, W3D, WVAR, PRECRT, SNOWRT,        &
                    SNOWPRT, GRPLPRT,                            & !new       
                    EFFC, EFFI, EFFS, EFFR,                      &
                    DT,                                          &
                    IMS, IME, JMS, JME, KMS, KME,                &
                    ITS, ITE, JTS, JTE, KTS, KTE,                &
                    QG3DTEN, NG3DTEN, QG3D, NG3D, EFFG,          &
                    qrcu1d, qscu1d, qicu1d,                      &
                    QGSTEN, QRSTEN, QISTEN, QNISTEN, QCSTEN,     &
                    nc3d, nc3dten, iinum,                        & ! wrf-chem
                    c2prec, CSED, ISED, SSED, GSED, RSED         & !wrf-chem 
                    ) bind(c)

    real(c_double) :: QC3DTEN(KTS:KTE), QI3DTEN(KTS:KTE), QNI3DTEN(KTS:KTE),    &
                      QR3DTEN(KTS:KTE), NI3DTEN(KTS:KTE),  NS3DTEN(KTS:KTE),    &
                      NR3DTEN(KTS:KTE), QC3D(KTS:KTE), QI3D(KTS:KTE),           &
                      QNI3D(KTS:KTE), QR3D(KTS:KTE), NI3D(KTS:KTE),             &
                      NS3D(KTS:KTE), NR3D(KTS:KTE), T3DTEN(KTS:KTE),            & 
                      QV3DTEN(KTS:KTE), T3D(KTS:KTE), QV3D(KTS:KTE)
    real(c_double) :: PRES(KTS:KTE), DZQ(KTS:KTE), W3D(KTS:KTE), WVAR(KTS:KTE)  
    real(c_double), value:: PRECRT, SNOWRT, SNOWPRT, GRPLPRT
    real(c_double) :: EFFC(KTS:KTE), EFFI(KTS:KTE), EFFS(KTS:KTE), EFFR(KTS:KTE)
    real(c_double), value:: DT
    integer(c_int), intent(in), value:: IMS,IME, JMS,JME, KMS,KME,              &
                                        ITS,ITE, JTS,JTE, KTS,KTE
    real(c_double) :: QG3DTEN(KTS:KTE), NG3DTEN(KTS:KTE), QG3D(KTS:KTE),        &
                      NG3D(KTS:KTE), EFFG(KTS:KTE),                             &
                      qrcu1d(KTS:KTE), qscu1d(KTS:KTE), qicu1d(KTS:KTE),        &
                      QGSTEN(KTS:KTE), QRSTEN(KTS:KTE), QISTEN(KTS:KTE),        &
                      QNISTEN(KTS:KTE), QCSTEN(KTS:KTE),                        &
                      nc3d(KTS:KTE), nc3dten(KTS:KTE)          !wrf-chem
    integer(c_int), intent(in), value:: iinum                  !wrf-chem
    real(c_double) :: c2prec(KTS:KTE), CSED(KTS:KTE), ISED(KTS:KTE),            &
                      SSED(KTS:KTE), GSED(KTS:KTE), RSED(KTS:KTE)      !wrf-chem


    call MORR_TWO_MOMENT_MICRO(QC3DTEN, QI3DTEN, QNI3DTEN, QR3DTEN,         &
                               NI3DTEN, NS3DTEN, NR3DTEN,                   &
                               QC3D, QI3D, QNI3D, QR3D, NI3D, NS3D, NR3D,   &
                               T3DTEN, QV3DTEN, T3D, QV3D,                  &
                               PRES, DZQ, W3D, WVAR, PRECRT, SNOWRT,        &
                               SNOWPRT, GRPLPRT,                            & !new
                               EFFC, EFFI, EFFS, EFFR,                      &
                               DT,                                          &
                               IMS, IME, JMS, JME, KMS, KME,                &
                               ITS, ITE, JTS, JTE, KTS, KTE,                & 
                               QG3DTEN, NG3DTEN, QG3D, NG3D, EFFG,          & 
                               qrcu1d, qscu1d, qicu1d,                      &
                               QGSTEN, QRSTEN, QISTEN, QNISTEN, QCSTEN,     &
                               nc3d, nc3dten, iinum,                        & ! wrf-chem
                               c2prec, CSED, ISED, SSED, GSED, RSED         & !wrf-chem
                               )
  end subroutine c_morr

end module morrison_2mom_wrap
