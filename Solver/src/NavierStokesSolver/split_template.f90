CAT(subroutine TimeDerivative_VolumetricContribution_Split_,SUFFIX)(mesh)
    use HexMeshClass
    use NodalStorageClass, only: NodalStorage
    use ElementClass
    use DGIntegrals
    use RiemannSolvers_NS
    implicit none    
    type(HexMesh), intent (inout) :: mesh

    ! ---------------
    ! Local variables
    ! ---------------
    integer                                :: i, j, k,l,eq, eID
    real(kind=RP)                          :: Flux(1:NCONS, 1:NDIM), Q_acc

    !$acc parallel present(mesh) vector_length(128) num_gangs(9750) async(1)
    !$acc loop gang
    do eID = 1 , size(mesh % elements)

       !$acc loop vector collapse(3) private(Flux)
       do k = 0, mesh % elements(eID) % Nxyz(3)  
          do j = 0, mesh % elements(eID) % Nxyz(2)  
             do i = 0, mesh % elements(eID) % Nxyz(1)
                call ViscousFlux_STATE( NCONS, NGRAD,  mesh % elements(eID) % storage % Q(:,i,j,k), mesh % elements(eID) % storage % U_x(:,i,j,k), & 
                                        mesh % elements(eID) % storage % U_y(:,i,j,k) , mesh % elements(eID) % storage % U_z(:,i,j,k), &
                                        mesh % elements(eID) % storage % mu_ns(1,i,j,k), 0.0_RP, &
                                        mesh % elements(eID) % storage % mu_ns(2,i,j,k), Flux)
                !$acc loop seq
                do eq = 1, NCONS
          
                  mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k,IX)  = - Flux(eq,IX) * mesh % elements(eID) % geom % jGradXi(IX,i,j,k)  &
                                                                                     - Flux(eq,IY) * mesh % elements(eID) % geom % jGradXi(IY,i,j,k)  &
                                                                                     - Flux(eq,IZ) * mesh % elements(eID) % geom % jGradXi(IZ,i,j,k)

                  mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k,IY) = - Flux(eq,IX) * mesh % elements(eID) % geom % jGradEta(IX,i,j,k)  &
                                                                                     - Flux(eq,IY) * mesh % elements(eID) % geom % jGradEta(IY,i,j,k)  &
                                                                                     - Flux(eq,IZ) * mesh % elements(eID) % geom % jGradEta(IZ,i,j,k)
                
                  mesh % elements(eID) % storage % contravariantFlux(eq,i,j,k,IZ) = - Flux(eq,IX) * mesh % elements(eID) % geom % jGradZeta(IX,i,j,k)  &
                                                                                     - Flux(eq,IY) * mesh % elements(eID) % geom % jGradZeta(IY,i,j,k)  &
                                                                                     - Flux(eq,IZ) * mesh % elements(eID) % geom % jGradZeta(IZ,i,j,k)
                end do
             end do
          end do
       end do

       !$acc loop vector collapse(4)
       do k = 0, mesh % elements(eID) % Nxyz(3)
          do j = 0, mesh % elements(eID) % Nxyz(2)
             do i = 0, mesh % elements(eID) % Nxyz(1)
                do eq = 1, NCONS
                   Q_acc = 0.0_RP
                   !$acc loop seq 
                   do l = 0, mesh % elements(eID) % Nxyz(1)
                      Q_acc = Q_acc + NodalStorage(mesh % elements(eID) % Nxyz(1)) % hatD(i,l) * mesh % elements(eID) % storage % contravariantFlux(eq,l,j,k,IX) + &
                                      NodalStorage(mesh % elements(eID) % Nxyz(2)) % hatD(j,l) * mesh % elements(eID) % storage % contravariantFlux(eq,i,l,k,IY) + &
                                      NodalStorage(mesh % elements(eID) % Nxyz(3)) % hatD(k,l) * mesh % elements(eID) % storage % contravariantFlux(eq,i,j,l,IZ)
                   end do
                   mesh % elements(eID) % storage % QDot(eq,i,j,k) = Q_acc
                end do
             end do
          end do
       end do

       ! *************************************
       ! Compute interior contravariant fluxes
       ! *************************************
       ! Compute inviscid contravariant flux
       ! -----------------------------------
       !$acc loop vector collapse(3) private(Flux)
       do k = 0, mesh % elements(eID) % Nxyz(3)  
          do j = 0, mesh % elements(eID) % Nxyz(2)  
             do i = 0, mesh % elements(eID) % Nxyz(1)
                !$acc loop seq
                do l = 0, mesh % elements(eID) % Nxyz(1)
                   call TWO_POINT_FLUX_FUNC( mesh % elements(eID) % storage % Q(:,i,j,k),  mesh % elements(eID) % storage % Q(:,l,j,k), mesh % elements(eID) % geom % jGradXi(:,i,j,k),  mesh % elements(eID) % geom % jGradXi(:,l,j,k), Flux(:,IX))
                   call TWO_POINT_FLUX_FUNC( mesh % elements(eID) % storage % Q(:,i,j,k),  mesh % elements(eID) % storage % Q(:,i,l,k), mesh % elements(eID) % geom % jGradEta(:,i,j,k), mesh % elements(eID) % geom % jGradEta(:,i,l,k), Flux(:,IY))
                   call TWO_POINT_FLUX_FUNC( mesh % elements(eID) % storage % Q(:,i,j,k),  mesh % elements(eID) % storage % Q(:,i,j,l), mesh % elements(eID) % geom % jGradZeta(:,i,j,k), mesh % elements(eID) % geom % jGradZeta(:,i,j,l), Flux(:,IZ))
                   !$acc loop seq
                   do eq = 1, NCONS
                      mesh % elements(eID) % storage % QDot(eq,i,j,k) = mesh % elements(eID) % storage % QDot(eq,i,j,k) &
                                                                      - NodalStorage(mesh % elements(eID) % Nxyz(1)) % sharpD(i,l) *  Flux(eq,IX) &
                                                                      - NodalStorage(mesh % elements(eID) % Nxyz(2)) % sharpD(j,l) *  Flux(eq,IY) &
                                                                      - NodalStorage(mesh % elements(eID) % Nxyz(3)) % sharpD(k,l) *  Flux(eq,IZ)
                   end do
                end do
             end do
          end do
       end do
    enddo
    !$acc end parallel loop
CAT(end subroutine TimeDerivative_VolumetricContribution_Split_,SUFFIX)