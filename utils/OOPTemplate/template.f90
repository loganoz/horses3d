!
!//////////////////////////////////////////////////////
!
module ClassTemplate
!  use ... statements
   use SMConstants, only: RP
   implicit none
!
!  *****************************
!  Default everything to private
!  *****************************
!
   private
!
!  ****************
!  Public variables
!  ****************
!
   public I1, I2
!
!  ******************
!  Public definitions
!  ******************
!
   public Class_t, SUB1_f
!
!  ****************************
!  Static variables definitions
!  ****************************
!
   integer, parameter   :: I1 = 1, I2 = 2
!
!  ****************
!  Class definition
!  ****************
!
   type Class_t
      logical                    :: constructed = .false.
      integer                    :: integer_property
      real(kind=RP)              :: real_property
      real(kind=RP), allocatable :: array(:)
      procedure(SUB1_f), pointer :: Variable_procedure => NULL()         ! Useful to have different procedures in run-time (e.g. interior or boundary condition)
      contains
         procedure         :: Destruct          => Class_Destruct
         procedure         :: Fixed_procedure   => Class_Fixed_procedure  ! Fixed_procedure == Class_Fixed_procedure always
         procedure, nopass :: Static_procedure  => Class_Static_procedure ! Procedure that does not depend on the class (i.e. no "self" input argument)
         generic           :: assignment(=)     => Class_Assign           ! Custom (=) operator
         procedure         :: Class_Assign
         generic, public   :: Generic_procedure => Class_generic1, &      ! Procedure that can be "Class_generic1" or "Class_generic" depending
                                                   Class_generic2         !    on input arguments.
         procedure, private :: Class_generic1              ! (Declarations of Class_generic#)
         procedure, private :: Class_generic2              ! (Declarations of Class_generic#)
   end type Class_t
!
!  *******************************************************************
!  Traditionally, constructors are exported with the name of the class
!  *******************************************************************
!
   interface Class_t
      module procedure ConstructClass1, ConstructClass2
   end interface Class_t
!
!  *******************
!  Function prototypes
!  *******************
!
   abstract interface
      subroutine SUB1_f(self,x,y)
         import Class_t
         implicit none
         class(Class_t),   intent(in)  :: self
         real(kind=RP),    intent(in)  :: x
         real(kind=RP),    intent(out) :: y
      end subroutine SUB1_f
   end interface
!
!  ========
   contains
!  ========
!
!/////////////////////////////////////////////////////////
!
!        Class constructor
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      function ConstructClass1()
!
!        ************************************************************
!        This is a default constructor (i.e. without input arguments)
!        ************************************************************
!
         implicit none
         type(Class_t)  :: ConstructClass1

         ConstructClass1 % constructed      = .true.
         ConstructClass1 % integer_property = 0
         ConstructClass1 % real_property    = 0.0
         ConstructClass1 % Variable_procedure => Default_procedure

      end function ConstructClass1
   
      function ConstructClass2(integer_property_, real_property_, Procedure_)
!
!        *************************************************
!        This is a constructor where we specify everything
!        *************************************************
!
         implicit none
         integer, intent(in)          :: integer_property_
         real(kind=RP),    intent(in) :: real_property_
         procedure(SUB1_f)            :: Procedure_
         type(Class_t)                :: ConstructClass2

         ConstructClass2 % constructed      = .true.
         ConstructClass2 % integer_property = integer_property_
         ConstructClass2 % real_property    = real_property_

         allocate(ConstructClass2 % array(ConstructClass2 % integer_property))
         ConstructClass2 % array = real_property_

         ConstructClass2 % Variable_procedure => Procedure_
         
      end function ConstructClass2
!
!/////////////////////////////////////////////////////////
!
!        Class destructors
!        -----------------
!
!/////////////////////////////////////////////////////////
!
      subroutine Class_Destruct(self)
         implicit none
         class(Class_t)    :: self

         if (.not. self % constructed) return
      
         self % integer_property = 0
         self % real_property    = 0.0

         if(allocated(self % array)) deallocate(self % array)

         self % Variable_procedure => NULL()

         self % constructed      = .false.

      end subroutine Class_Destruct
!
!/////////////////////////////////////////////////////////
!  
!        Suitable subroutines for Variable_procedure
!        ---------------------------------------------
!
!/////////////////////////////////////////////////////////
!
      subroutine Default_procedure(self,x,y)
         implicit none
         class(Class_t), intent(in) :: self
         real(kind=RP), intent(in)           :: x
         real(kind=RP), intent(out)          :: y

         y = x + self % real_property

      end subroutine Default_procedure
!
!/////////////////////////////////////////////////////////
!
!        Overloading operators
!        ---------------------
!
!/////////////////////////////////////////////////////////
!
      subroutine Class_Assign(to, from)
!
!        ********************************************
!        Custom subroutine to perform class2 = class1
!        ********************************************
!
         implicit none
         class(Class_t), intent(out)   :: to
         class(Class_t), intent(in)    :: from

         to % integer_property = from % integer_property
         to % real_property    = from % real_property

         if ( allocated(to % array)) deallocate(to % array)
         if ( allocated(from % array)) allocate(to % array(size(from % array)))
         to % array = from % array

         to % Variable_procedure => from % Variable_procedure
         to % constructed = .true.

      end subroutine Class_Assign
!
!/////////////////////////////////////////////////////////
!
!        Additional procedures
!        ---------------------
!
!/////////////////////////////////////////////////////////
!
      subroutine Class_Fixed_procedure(self, arg1)
         implicit none
         class(Class_t),   intent(inout) :: self
         integer,          intent(in)    :: arg1

         self % array = arg1

      end subroutine Class_Fixed_procedure

      subroutine Class_Static_procedure(x,y)
!
!        ********************************************************
!        This procedure does not depend on the class (e.g. y=x^2)
!        ********************************************************
!
         implicit none  
         real(kind=RP), intent(in)  :: x
         real(kind=RP), intent(out) :: y

         y = x*x

      end subroutine Class_Static_procedure

      subroutine Class_generic1(self)
         implicit none
         class(Class_t) :: self       

      end subroutine Class_generic1

      subroutine Class_generic2(self, arg1)
         implicit none
         class(Class_t) :: self
         integer        :: arg1

      end subroutine Class_generic2
   
end module ClassTemplate
