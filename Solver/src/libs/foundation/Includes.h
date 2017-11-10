
#define POW2(x) ((x)*(x))
#define POW3(x) (POW2((x))*(x))

#define errorMessage(UNIT) write(UNIT,'(A,A,A,I0,A)')   "Error in file ", __FILE__ , ", in line " , __LINE__ ,"."
#define stopMessage(UNIT)  write(UNIT,'(A,A,A,I0,A)') "Stopped in file ", __FILE__ , ", in line " , __LINE__ ,"."
#define safedeallocate(x) if (allocated(x)) deallocate(x)
