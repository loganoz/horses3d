module Physics
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
   use Physics_NS
#elif defined(NAVIERSTOKES) && (SPALARTALMARAS)
   use Physics_NSSA
#elif defined(INCNS)
   use Physics_iNS
#elif defined(MULTIPHASE)
   use Physics_MU
#elif defined(SCALAR)
   use Physics_SLR
#elif defined(SCALAR_INS_V04)
   use Physics_SLR_INS_V04
#endif

#if defined(CAHNHILLIARD)
   use Physics_CH
#endif
   implicit none

end module Physics