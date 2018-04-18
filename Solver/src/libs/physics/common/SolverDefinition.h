#if (defined(NAVIERSTOKES) && defined(CAHNHILLIARD))
#define SOLVER NSCH_SOLVER

#elif defined(NAVIERSTOKES)
#define SOLVER NS_SOLVER

#else
#define SOLVER CH_SOLVER

#endif
