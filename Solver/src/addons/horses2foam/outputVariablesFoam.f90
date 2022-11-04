!
! ////////////////////////////////////////////////////////////////////
!   HORSES3D to Foam Result - outputVariablesFoam Module
!
!      This module generate outputVariables
!
!////////////////////////////////////////////////////////////////////////////////////////

module outputVariablesFoam
!
!//////////////////////////////////////////////////////////////////////////////////
!
!        This module selects user-defined output variables for the 
!     tecplot file. The user may add extra variables, this can be done by
!     following the steps:
!
!           * Adding a new ID for the variable NEWVAR_V (following the given order).
!           * Adding a new key for the variable, NEWVARKey.
!           * Adding the new key to the keys array.
!           * Adding the "select case" procedure that computes the output
!                 variable.
!
!//////////////////////////////////////////////////////////////////////////////////
!
   use SMConstants
   use PhysicsStorage
   use Headers
   use storageFoam
   use createFoamMeshFile

   private
   public   no_of_outputVariables, getNoOfCommas
   public   getOutputVariables, ComputeOutputVariables, getOutputVariablesLabel

   type resultVariables_t
      character(len=LINE_LENGTH) :: Name
      integer      , allocatable :: Code (:)
      integer                    :: Dimension(7)
      real(kind=RP), allocatable :: outputVars(:,:,:,:,:)
   end type resultVariables_t
   integer, parameter   :: STR_VAR_LEN = 16
!
!  ***************************
!  Variables without gradients
!  ***************************
!
   enum, bind(C)
      enumerator :: Q_V = 1, RHO_V, U_V, V_V, W_V
      enumerator :: P_V, P0_V, T_V, Mach_V, S_V, Vabs_V
      enumerator :: Vvec_V, Ht_V, RHOU_V, RHOV_V
      enumerator :: RHOW_V, RHOE_V, C_V, Nxi_V, Neta_V
      enumerator :: Nzeta_V, Nav_V, N_V
      enumerator :: Xi_V, Eta_V, Zeta_V, ThreeAxes_V, Axes_V, eID_V, MPIRANK_V
      enumerator :: GRADV_V, UX_V, VX_V, WX_V
      enumerator :: UY_V, VY_V, WY_V, UZ_V, VZ_V, WZ_V
      enumerator :: CX_V, CY_V, CZ_V
      enumerator :: OMEGA_V, OMEGAX_V, OMEGAY_V, OMEGAZ_V
      enumerator :: OMEGAABS_V, QCRIT_V
      enumerator :: LASTVARIABLE
   end enum

   integer, parameter   :: NO_OF_VARIABLES = LASTVARIABLE-1
   integer, parameter   :: NO_OF_INVISCID_VARIABLES = MPIRANK_V

   character(len=STR_VAR_LEN), parameter  :: QKey          = "Q"
   character(len=STR_VAR_LEN), parameter  :: RHOKey        = "rho"
   character(len=STR_VAR_LEN), parameter  :: UKey          = "u"
   character(len=STR_VAR_LEN), parameter  :: VKey          = "v"
   character(len=STR_VAR_LEN), parameter  :: WKey          = "w"
   character(len=STR_VAR_LEN), parameter  :: PKey          = "p"
   character(len=STR_VAR_LEN), parameter  :: P0Key         = "pt"
   character(len=STR_VAR_LEN), parameter  :: TKey          = "T"
   character(len=STR_VAR_LEN), parameter  :: MachKey       = "Mach"
   character(len=STR_VAR_LEN), parameter  :: SKey          = "s"
   character(len=STR_VAR_LEN), parameter  :: VabsKey       = "Vabs"
   character(len=STR_VAR_LEN), parameter  :: VvecKey       = "V"
   character(len=STR_VAR_LEN), parameter  :: HtKey         = "Ht"
   character(len=STR_VAR_LEN), parameter  :: RHOUKey       = "rhou"
   character(len=STR_VAR_LEN), parameter  :: RHOVKey       = "rhov"
   character(len=STR_VAR_LEN), parameter  :: RHOWKey       = "rhow"
   character(len=STR_VAR_LEN), parameter  :: RHOEKey       = "rhoe"
   character(len=STR_VAR_LEN), parameter  :: cKey          = "c"
   character(len=STR_VAR_LEN), parameter  :: NxiKey        = "Nxi"
   character(len=STR_VAR_LEN), parameter  :: NetaKey       = "Neta"
   character(len=STR_VAR_LEN), parameter  :: NzetaKey      = "Nzeta"
   character(len=STR_VAR_LEN), parameter  :: NavKey        = "Nav"
   character(len=STR_VAR_LEN), parameter  :: NKey          = "N"
   character(len=STR_VAR_LEN), parameter  :: XiKey   = "Ax_Xi"
   character(len=STR_VAR_LEN), parameter  :: EtaKey  = "Ax_Eta"
   character(len=STR_VAR_LEN), parameter  :: ZetaKey = "Ax_Zeta"
   character(len=STR_VAR_LEN), parameter  :: ThreeAxesKey = "ThreeAxes"
   character(len=STR_VAR_LEN), parameter  :: AxesKey = "Axes"   
   character(len=STR_VAR_LEN), parameter  :: eIDKey  = "eID"   
   character(len=STR_VAR_LEN), parameter  :: mpiRankKey    = "mpi_rank"
   character(len=STR_VAR_LEN), parameter  :: gradVKey      = "gradV"
   character(len=STR_VAR_LEN), parameter  :: uxKey         = "u_x"
   character(len=STR_VAR_LEN), parameter  :: vxKey         = "v_x"
   character(len=STR_VAR_LEN), parameter  :: wxKey         = "w_x"
   character(len=STR_VAR_LEN), parameter  :: uyKey         = "u_y"
   character(len=STR_VAR_LEN), parameter  :: vyKey         = "v_y"
   character(len=STR_VAR_LEN), parameter  :: wyKey         = "w_y"
   character(len=STR_VAR_LEN), parameter  :: uzKey         = "u_z"
   character(len=STR_VAR_LEN), parameter  :: vzKey         = "v_z"
   character(len=STR_VAR_LEN), parameter  :: wzKey         = "w_z"
   character(len=STR_VAR_LEN), parameter  :: cxKey         = "c_x"
   character(len=STR_VAR_LEN), parameter  :: cyKey         = "c_y"
   character(len=STR_VAR_LEN), parameter  :: czKey         = "c_z"
   character(len=STR_VAR_LEN), parameter  :: omegaKey      = "omega"
   character(len=STR_VAR_LEN), parameter  :: omegaxKey     = "omega_x"
   character(len=STR_VAR_LEN), parameter  :: omegayKey     = "omega_y"
   character(len=STR_VAR_LEN), parameter  :: omegazKey     = "omega_z"
   character(len=STR_VAR_LEN), parameter  :: omegaAbsKey   = "omega_abs"
   character(len=STR_VAR_LEN), parameter  :: QCriterionKey = "Qcrit"
   
   character(len=STR_VAR_LEN), dimension(NO_OF_VARIABLES), parameter  :: variableNames = (/ QKey, RHOKey, UKey, VKey, WKey, &
                                                                            PKey, P0Key, TKey, MachKey, SKey, VabsKey, &
                                                                            VvecKey, HtKey, RHOUKey, RHOVKey, RHOWKey, &
                                                                            RHOEKey, cKey, NxiKey, NetaKey, NzetaKey, NavKey, &
                                                                            NKey, &
                                                                            XiKey, EtaKey, ZetaKey, ThreeAxesKey, AxesKey, &
                                                                            eIDKey, mpiRankKey, &
                                                                            gradVKey, uxKey, vxKey, wxKey, &
                                                                            uyKey, vyKey, wyKey, uzKey, vzKey, wzKey, &
                                                                            cxKey, cyKey, czKey, &
                                                                            omegaKey, omegaxKey, omegayKey, omegazKey, &
                                                                            omegaAbsKey, QCriterionKey/)
                                                               
   integer                :: no_of_outputVariables
   integer, allocatable   :: outputVariableNames(:)
   logical                :: outScale

   contains
!
!/////////////////////////////////////////////////////////////////////////////
!
        SUBROUTINE constructResultVariables(nVariables,codeVariables,mesh)
            IMPLICIT NONE
            INTEGER                    ,INTENT(IN)     :: nVariables
            INTEGER                    ,INTENT(IN)     :: codeVariables(nVariables)
            TYPE(Mesh_t)               ,INTENT(IN)     :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
            TYPE(resultVariables_t) ::self
            INTEGER                 :: i,j,k,nV
            INTEGER                 :: fid
            INTEGER                 :: eID
            INTEGER                 :: N(3)
            INTEGER     ,ALLOCATABLE:: output(:)
            character(len=LINE_LENGTH) :: formatout

            N(:)=Mesh % elements(1) % Nout(:)

            DO nV=1,nVariables
                j=outputVariablesForVariable(codeVariables(nV))
                ALLOCATE (output(j), self % Code(j))
                CALL OutputVariablesForPreliminarVariable(codeVariables(nV), output)
                self % Code(:) = output
                self % Name = variableNames(codeVariables(nV))
                self % Dimension(:) = 0 ! Assign 0 or no Unit (need modification)
                ALLOCATE(self % outputVars(1:size(mesh % elements), j, 0:N(1), 0:N(2), 0:N(3)))
                self % outputVars(:,:,:,:,:)=0
                DEALLOCATE(output)
!
!
                DO eID = 1, size(mesh % elements)
                    ASSOCIATE ( e => mesh % elements(eID) )
                    N = e % Nout
                    CALL ComputeOutputVariables(eID, e % Nout, e, self, mesh % refs, mesh % hasGradients)
                    END ASSOCIATE
                END DO
!
!
                CALL createFileResultHeader (self % Name, size(self % code), fid)

                formatout=getFormat(size(self % code))

                WRITE(fid,'(A)')'dimensions      [0 0 0 0 0 0 0];'
                WRITE(fid,*)
                IF (size(self % code).EQ.1) THEN
                    WRITE(fid,'(A)')'internalField   nonuniform List<scalar> '
                ELSE IF (size(self % code).EQ.3) THEN
                    WRITE(fid,'(A)')'internalField   nonuniform List<vector> '
                ELSE
                    WRITE(fid,'(A)')'internalField   nonuniform List<tensor> '
                END IF
                WRITE(fid,*)
                WRITE(fid,'(I0)')size(mesh % elements)*(N(1)+1)*(N(2)+1)*(N(3)+1)
                WRITE(fid,'(A,I10,I10,I10,I10,A)')'('
                DO eID = 1, size(mesh % elements)
                    do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                    IF (size(self % code).EQ.1) THEN
                        write(fid,'(E18.10)')self % outputVars(eID,1,i,j,k)
                    ELSE
                        write(fid,trim(formatout))'( ', self % outputVars(eID,:,i,j,k), ' )'
                    END IF
                    end do                ; end do                ; end do
                END DO
                WRITE(fid,'(A)')')'
                WRITE(fid,'(A)')';'
                CLOSE(fid)
                DEALLOCATE (self % Code, self % outputVars)
            END DO

        END SUBROUTINE constructResultVariables
!
!/////////////////////////////////////////////////////////////////////////////
!
      subroutine getOutputVariables(flag_name, mesh)
         implicit none
         character(len=*), intent(in)  :: flag_name
         TYPE(Mesh_t)               ,INTENT(IN)     :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         logical                       :: flagPresent
         integer                       :: pos, pos2, i, preliminarNoOfVariables
         integer, allocatable          :: preliminarVariables(:)
         character(len=LINE_LENGTH)    :: flag
         character(len=STR_VAR_LEN)    :: inputVar

		 write(STD_OUT,'(10X,A,A)') "Creating Output File:"
         write(STD_OUT,'(10X,A,A)') "--------------------"
		 
         flagPresent = .false.
         outScale    = .true.

         pos = LEN(flag_name)

         if ( pos .ne. 0 ) then
            flagPresent = .true.
         end if
!
!           Also, check if the dimensionless version is requested
!           -----------------------------------------------------
            pos = index(trim(flag_name),"--dimensionless")
            if ( pos .ne. 0 ) then
               outScale = .false.
            end if

            flag=flag_name
!
!        ***********************************************************
!              Read the variables. They are first loaded into
!           a "preliminar" variables, prior to introduce them in
!           the real outputVariables. This is because some of
!           the output variables lead to multiple variables (e.g.
!           Q, V or N).
!        ***********************************************************
!
         if ( .not. flagPresent ) then
!
!           Default: export Q
!           -------
            preliminarNoOfVariables = 1
            allocate( preliminarVariables(preliminarNoOfVariables) )
            preliminarVariables(1) = Q_V

         else
!
!           Prepare to read the variable names
!           ----------------------------------
            preliminarNoOfVariables = getNoOfCommas(trim(flag)) + 1
            allocate( preliminarVariables(preliminarNoOfVariables) )

            if ( preliminarNoOfVariables .eq. 1 ) then
               read(flag(1:len_trim(flag)),*) inputVar
               preliminarVariables(1) = outputVariableForName(ADJUSTL(trim(inputVar)))

            else
               pos=0
               do i = 1, preliminarNoOfVariables-1
                  pos2 = index(trim(flag(pos+1:)),",") + pos
                  read(flag(pos+1:pos2),*) inputVar
                  preliminarVariables(i) = outputVariableForName(ADJUSTL(trim(inputVar)))

                  pos = pos2
               end do

               pos = index(trim(flag),",",BACK=.true.)
               preliminarVariables(preliminarNoOfVariables) = outputVariableForName(ADJUSTL(flag(pos+1:)))

            end if
         end if
!
!        *******************************************************
!        Convert the preliminary variables into output variables
!        *******************************************************
!
         no_of_outputVariables = 0

         do i = 1, preliminarNoOfVariables
            no_of_outputVariables = no_of_outputVariables + outputVariablesForVariable(preliminarVariables(i))
         end do

         if (hasMPIranks)  no_of_outputVariables = no_of_outputVariables + 1

         allocate( outputVariableNames(no_of_outputVariables) )

         pos = 1
         do i = 1, preliminarNoOfVariables
            pos2 = pos + outputVariablesForVariable(preliminarVariables(i)) - 1

            call outputVariablesForPreliminarVariable(preliminarVariables(i), outputVariableNames(pos:pos2) )

            pos = pos + outputVariablesForVariable(preliminarVariables(i))

         end do

         if (hasMPIranks) outputVariableNames(pos) = MPIRANK_V
!
!        **************************
!        Construct Result Variables
!        **************************
!
         CALL constructResultVariables (preliminarNoOfVariables,preliminarVariables,mesh)
!
!        *****************************
!        Describe the output variables
!        *****************************
!		 
		 write(STD_OUT,'(20X,A,A)') "Selected Output Variables:"

         do i = 1, no_of_outputVariables
            write(STD_OUT,'(30X,A,A20)') "->", trim(variableNames(outputVariableNames(i)))
         end do

         if ( outScale ) then
            write(STD_OUT,'(30X,A,A)') "-> Variables are exported with dimensions."

         else
            write(STD_OUT,'(30X,A,A)') "-> Dimensionless mode."

         end if
		 
		 deallocate(outputVariableNames)

      end subroutine getOutputVariables

!
!////////////////////////////////////////////////////////////////////////
!   This subroutine compute the output variable and store it in the ...
!
!
      subroutine ComputeOutputVariables(eID, N, e, resultVariables, refs, hasGradients)
         use SolutionFile
         use storageFoam
         implicit none
         integer, intent(in)                    :: eID
         integer, intent(in)                    :: N(3)
         type(resultVariables_t)        ,INTENT(INOUT)  :: resultVariables
         type(Element_t), intent(in)           :: e
         real(kind=RP), intent(in)              :: refs(NO_OF_SAVED_REFS)
         logical,       intent(in)              :: hasGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: var2, i, j, k
         real(kind=RP) :: Sym, Asym


         do var2 = 1, size(resultVariables % Code )
            if ( hasGradients .or. (resultVariables % Code(var2) .le. NO_OF_INVISCID_VARIABLES ) ) then
                        associate ( Q   => e % Qout, &
                           U_x => e % U_xout, & 
                           U_y => e % U_yout, & 
                           U_z => e % U_zout, &
                        output => resultVariables)
               select case (resultVariables % Code(var2))
   
               case(RHO_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(RHO_REF) * output % outputVars(eID,var2,:,:,:)
   
               case(U_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = Q(IRHOU,i,j,k) / Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(V_REF) * output % outputVars(eID,var2,:,:,:)
    
               case(V_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = Q(IRHOV,i,j,k) / Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(V_REF) * output % outputVars(eID,var2,:,:,:)
    
               case(W_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = Q(IRHOW,i,j,k) / Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(V_REF) * output % outputVars(eID,var2,:,:,:)
   
               case(P_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = (refs(GAMMA_REF) - 1.0_RP)*(Q(IRHOE,i,j,k) - 0.5_RP*&
                                       ( (Q(IRHOU,i,j,k)**2) + (Q(IRHOV,i,j,k)**2) + (Q(IRHOW,i,j,k)**2)) /Q(IRHO,i,j,k))
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(RHO_REF) * (refs(V_REF)**2) &
                                    * output % outputVars(eID,var2,:,:,:)
									
			   case(P0_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = ( (Q(IRHOU,i,j,k)**2) + (Q(IRHOV,i,j,k)**2) + (Q(IRHOW,i,j,k)**2)) &
                                                             /(Q(IRHO,i,j,k)**2)     ! Vabs**2
                     output % outputVars(eID,var2,i,j,k) =  output % outputVars(eID,var2,i,j,k) / ( refs(GAMMA_REF)*  &
                                        (refs(GAMMA_REF)-1.0_RP)*(Q(IRHOE,i,j,k) /Q(IRHO,i,j,k)-0.5_RP * &
                                            output % outputVars(eID,var2,i,j,k)) )  ! Mach Number**2
					 output % outputVars(eID,var2,i,j,k) = (1+((refs(GAMMA_REF)-1.0_RP)*0.5_RP)*output % outputVars(eID,var2,i,j,k))** &
										(refs(GAMMA_REF)/(refs(GAMMA_REF)-1.0_RP))
											
                     output % outputVars(eID,var2,i,j,k) = (refs(GAMMA_REF) - 1.0_RP)*(Q(IRHOE,i,j,k) - 0.5_RP*&
                                       ( (Q(IRHOU,i,j,k)**2) + (Q(IRHOV,i,j,k)**2) + (Q(IRHOW,i,j,k)**2)) /Q(IRHO,i,j,k)) &
									    * output % outputVars(eID,var2,i,j,k) 
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(RHO_REF) * (refs(V_REF)**2) &
                                    * output % outputVars(eID,var2,:,:,:)
   
               case(T_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = (refs(GAMMA_REF) - 1.0_RP) * refs(GAMMA_REF) * (refs(MACH_REF)**2) &
                              / Q(IRHO,i,j,k) * (Q(IRHOE,i,j,k) - 0.5_RP*( (Q(IRHOU,i,j,k)**2) + (Q(IRHOV,i,j,k)**2) + &
                              (Q(IRHOW,i,j,k)**2))  /Q(IRHO,i,j,k))
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(T_REF) * output % outputVars(eID,var2,:,:,:)
   
               case(MACH_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = ( (Q(IRHOU,i,j,k)**2) + (Q(IRHOV,i,j,k)**2) + (Q(IRHOW,i,j,k)**2)) &
                                                             /(Q(IRHO,i,j,k)**2)     ! Vabs**2
                     output % outputVars(eID,var2,i,j,k) = sqrt( output % outputVars(eID,var2,i,j,k) / ( refs(GAMMA_REF)*  &
                                        (refs(GAMMA_REF)-1.0_RP)*(Q(IRHOE,i,j,k) /Q(IRHO,i,j,k)-0.5_RP * &
                                            output % outputVars(eID,var2,i,j,k)) ) )
                  end do         ; end do         ; end do
   
               case(S_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = refs(GAMMA_REF) * (refs(GAMMA_REF) - 1.0_RP)*(Q(IRHOE,i,j,k) - 0.5_RP * &
                                       ( (Q(IRHOU,i,j,k)**2) + (Q(IRHOV,i,j,k)**2) + (Q(IRHOW,i,j,k)**2)) /Q(IRHO,i,j,k)) &
                                       / (Q(IRHO,i,j,k)**refs(GAMMA_REF))
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(RHO_REF)*(refs(V_REF)**2)/refs(RHO_REF)** &
                                                        refs(GAMMA_REF) * output % outputVars(eID,var2,:,:,:)
      
               case(Vabs_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = sqrt((Q(IRHOU,i,j,k)**2) + (Q(IRHOV,i,j,k)**2) + &
                                                        (Q(IRHOW,i,j,k)**2))/Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(V_REF) * output % outputVars(eID,var2,:,:,:)
   
               case(Ht_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = refs(GAMMA_REF)*Q(IRHOE,i,j,k) - 0.5_RP*(refs(GAMMA_REF)-1.0_RP)*&
                                       ( (Q(IRHOU,i,j,k)**2) + (Q(IRHOV,i,j,k)**2) + (Q(IRHOW,i,j,k)**2)) /Q(IRHO,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(RHO_REF) * (refs(V_REF)**2) * output % &
                                    outputVars(eID,var2,:,:,:)
   
               case(RHOU_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = Q(IRHOU,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(RHO_REF) * refs(V_REF) * output % &
                                    outputVars(eID,var2,:,:,:)
   
               case(RHOV_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = Q(IRHOV,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(RHO_REF) * refs(V_REF) * output % &
                                    outputVars(eID,var2,:,:,:)
         
               case(RHOW_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = Q(IRHOW,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(RHO_REF) * refs(V_REF) * output % &
                                    outputVars(eID,var2,:,:,:)
   
               case(RHOE_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = Q(IRHOE,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = refs(RHO_REF) * (refs(V_REF)**2) * output % &
                                    outputVars(eID,var2,:,:,:)
                  
               case(Nxi_V)
                  output % outputVars(eID,var2,:,:,:) = e % Nsol(1)
                  
               case(Neta_V)
                  output % outputVars(eID,var2,:,:,:) = e % Nsol(2)
                  
               case(Nzeta_V)
                  output % outputVars(eID,var2,:,:,:) = e % Nsol(3)
                  
               case(Nav_V)
                  output % outputVars(eID,var2,:,:,:) = sum(e % Nsol)/real(NDIM)
                  
               case(Xi_V)
                  output % outputVars(eID,var2,:,:,:) = 0
                  output % outputVars(eID,var2,:,0,0) = 3
                  output % outputVars(eID,var2,0,:,0) = 0
                  output % outputVars(eID,var2,0,0,:) = 0
                  
               case(Eta_V)
                  output % outputVars(eID,var2,:,:,:) = 0
                  output % outputVars(eID,var2,:,0,0) = 0
                  output % outputVars(eID,var2,0,:,0) = 3
                  output % outputVars(eID,var2,0,0,:) = 0
                  
               case(Zeta_V)
                  output % outputVars(eID,var2,:,:,:) = 0
                  output % outputVars(eID,var2,:,0,0) = 0
                  output % outputVars(eID,var2,0,:,0) = 0
                  output % outputVars(eID,var2,0,0,:) = 3
                  
               case(ThreeAxes_V)
                  output % outputVars(eID,var2,:,:,:) = 0
                  output % outputVars(eID,var2,:,0,0) = 3
                  output % outputVars(eID,var2,0,:,0) = 1.5
                  output % outputVars(eID,var2,0,0,:) = 1.5
               
               case(eID_V)
                  output % outputVars(eID,var2,:,:,:) = e % eID
                  
               case(MPIRANK_V)
                  output % outputVars(eID,var2,:,:,:) = e % mpi_rank
!
!
!              ******************
!              Gradient variables   
!              ******************
!
               case(UX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = U_x(1,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)**2) * output % outputVars(eID,var2,:,:,:)
               
               case(VX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = U_x(2,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)**2) * output % outputVars(eID,var2,:,:,:)
               
               case(WX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = U_x(3,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)**2) * output % outputVars(eID,var2,:,:,:)
               
               case(UY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = U_y(1,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)**2) * output % outputVars(eID,var2,:,:,:)
               
               case(VY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = U_y(2,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)**2) * output % outputVars(eID,var2,:,:,:)
               
               case(WY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = U_y(3,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)**2) * output % outputVars(eID,var2,:,:,:)
               
               case(UZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = U_z(1,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)**2) * output % outputVars(eID,var2,:,:,:)
               
               case(VZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = U_z(2,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)**2) * output % outputVars(eID,var2,:,:,:)
               
               case(WZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = U_z(3,i,j,k)
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)**2) * output % outputVars(eID,var2,:,:,:)

               case(OMEGAX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = ( U_y(3,i,j,k) - U_z(2,i,j,k) )
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)) * output % outputVars(eID,var2,:,:,:)

               case(OMEGAY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = ( U_z(1,i,j,k) - U_x(3,i,j,k) )
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)) * output % outputVars(eID,var2,:,:,:)

               case(OMEGAZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = ( U_x(2,i,j,k) - U_y(1,i,j,k) )
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)) * output % outputVars(eID,var2,:,:,:)

               case(OMEGAABS_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = sqrt(  ( U_y(3,i,j,k) - U_z(2,i,j,k) )**2 &
                                              + ( U_z(1,i,j,k) - U_x(3,i,j,k) )**2 &
                                              + ( U_x(2,i,j,k) - U_y(1,i,j,k) )**2 )
                  end do         ; end do         ; end do
                  if ( outScale ) output % outputVars(eID,var2,:,:,:) = (refs(V_REF)) * output % outputVars(eID,var2,:,:,:)

               case(QCRIT_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     Sym =   ( U_x(1,i,j,k) )**2 + ( U_y(2,i,j,k) )**2 + ( U_z(3,i,j,k) )**2  &
                           + 2.0_RP *( ( 0.5_RP * (U_x(2,i,j,k) + U_y(1,i,j,k)) )**2 +          &
                                       ( 0.5_RP * (U_x(3,i,j,k) + U_z(1,i,j,k)) )**2 +          &
                                       ( 0.5_RP * (U_y(3,i,j,k) + U_z(2,i,j,k)) )**2 )

                     Asym =   2.0_RP *( ( 0.5_RP * (U_x(2,i,j,k) - U_y(1,i,j,k)) )**2 +        &
                                        ( 0.5_RP * (U_x(3,i,j,k) - U_z(1,i,j,k)) )**2 +        &
                                        ( 0.5_RP * (U_y(3,i,j,k) - U_z(2,i,j,k)) )**2 )

                     output % outputVars(eID,var2,i,j,k) = 0.5_RP*( Asym - Sym )
                  end do            ; end do            ; end do

               case(C_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = Q(size(Q,1),i,j,k)
                  end do         ; end do         ; end do
               
               case(CX_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = U_x(size(U_x,1),i,j,k)
                  end do         ; end do         ; end do
               
               case(CY_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = U_y(size(U_y,1),i,j,k)
                  end do         ; end do         ; end do
               
               case(CZ_V)
                  do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
                     output % outputVars(eID,var2,i,j,k) = U_z(size(U_z,1),i,j,k)
                  end do         ; end do         ; end do
               
               end select
               end associate
   
            else
               ASSOCIATE(output => resultVariables)
               output % outputVars(eID,var2,:,:,:) = 0.0_RP
               END ASSOCIATE

            end if
         end do
      end subroutine ComputeOutputVariables

      character(len=1024) function getOutputVariablesLabel()
         implicit none
!
!        ---------------
!        Local variables         
!        ---------------
!
         integer  :: iVar

         getOutputVariablesLabel = ""
         do iVar = 1, no_of_outputVariables
            write(getOutputVariablesLabel,'(A,A,A,A)') trim(getOutputVariablesLabel), ',"', &
            trim(variableNames(outputVariableNames(iVar))),'"'
         end do

      end function getOutputVariablesLabel

      integer function outputVariablesForVariable(iVar)
!
!        ************************************************
!           This subroutine specifies if a variable
!           implies more than one variable, e.g.,
!           Q = [rho,rhou,rhov,rhow,rhoe], or
!           V = [u,v,w]
!        ************************************************
!
         implicit none
         integer,    intent(in)     :: iVar
   
         select case(iVar)
   
         case(Q_V)
            outputVariablesForVariable = NVARS

         case(Vvec_V)
            outputVariablesForVariable = 3

         case(gradV_V)
            outputVariablesForVariable = 9

         case(omega_V)
            outputVariablesForVariable = 3
            
         case(N_V)
            outputVariablesForVariable = 4
            
         case(Axes_V)
            outputVariablesForVariable = 4

         case default
            outputVariablesForVariable = 1
      
         end select

      end function outputVariablesForVariable

      subroutine OutputVariablesForPreliminarVariable(iVar, output)
         implicit none
         integer, intent(in)     :: iVar
         integer, intent(out)    :: output(:)

         select case(iVar)

         case(Q_V)
            if ( NVARS .eq. 5 ) then
               output = (/RHO_V, RHOU_V, RHOV_V, RHOW_V, RHOE_V/)
      
            elseif ( NVARS .eq. 6 ) then
               output = (/RHO_V, RHOU_V, RHOV_V, RHOW_V, RHOE_V, C_V/)

            elseif ( NVARS .eq. 1 ) then
               output = (/C_V/)

            end if

         case(Vvec_V)
            output = (/U_V, V_V, W_V/)

         case(gradV_V)
            output = (/UX_V, VX_V, WX_V, UY_V, VY_V, WY_V, UZ_V, VZ_V, WZ_V/)

         case(omega_V)
            output = (/OMEGAX_V, OMEGAY_V, OMEGAZ_V/)
            
         case(N_V)
            output = (/Nxi_V, Neta_V, Nzeta_V, Nav_V/)
         
         case(Axes_V)
            output = (/Xi_V, Eta_V, Zeta_V, ThreeAxes_V/)
            
         case default
            output = iVar

         end select

      end subroutine OutputVariablesForPreliminarVariable

      integer function outputVariableForName(variableName)
         implicit none
         character(len=*),    intent(in)     :: variableName

         do outputVariableForName = 1, NO_OF_VARIABLES
            if ( trim(variableNames(outputVariableForName)) .eq. trim(variableName) ) return
         end do
!
!        Return "-1" if not found
!        ------------------------
         outputVariableForName = -1
      
      end function outputVariableForName

      integer function getNoOfCommas(input_line)
         implicit none
         character(len=*), intent(in)     :: input_line
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i

         getNoOfCommas = 0
         do i = 1, len_trim(input_line)
            if ( input_line(i:i) .eq. "," ) getNoOfCommas = getNoOfCommas + 1 
         end do

      end function getNoOfCommas
!
!////////////////////////////////////////////////////////////////////////
!   This subroutine compute the output variable and store it in the ...
!
!
      character(len=LINE_LENGTH) function getFormat(nData)
         implicit none
         INTEGER            ,INTENT(IN) :: nData

         INTEGER        :: i

         getFormat = '(A,'

         DO i=1,nData
            write(getFormat,'(A,A)')trim(getFormat),'E18.10,'
         END DO
         write(getFormat,'(A,A)')trim(getFormat),'A)'

      end function getFormat

end module outputVariablesFoam
