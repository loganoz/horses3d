module Physics
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
   use Physics_NS
#elif defined(NAVIERSTOKES) && (SPALARTALMARAS)
   use Physics_NSSA
#elif defined(INCNS)
   use Physics_iNS
#elif defined(MULTIPHASE)
   use Physics_MU
#elif defined(ACOUSTIC)
   use Physics_CAA
#endif
#if defined(CAHNHILLIARD)
   use Physics_CH
#endif
   implicit none

end module Physics
