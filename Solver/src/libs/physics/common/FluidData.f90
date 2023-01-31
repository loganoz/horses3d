module FluidData
#if defined(NAVIERSTOKES) && (!(SPALARTALMARAS))
   use FluidData_NS
#elif defined(INCNS)
   use FluidData_iNS
#elif defined(MULTIPHASE)
   use FluidData_MU
#elif defined(SPALARTALMARAS)
   use FluidData_NSSA
#endif
#if defined(CAHNHILLIARD)
   use FluidData_CH
#endif
   implicit none

end module FluidData