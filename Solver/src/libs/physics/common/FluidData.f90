module FluidData
#if defined(SPALARTALMARAS)
   use FluidData_NSSA
#elif defined(NAVIERSTOKES)
   use FluidData_NS
#elif defined(INCNS)
   use FluidData_iNS
#elif defined(MULTIPHASE)
   use FluidData_MU
#endif
#if defined(CAHNHILLIARD)
   use FluidData_CH
#endif
   implicit none

end module FluidData