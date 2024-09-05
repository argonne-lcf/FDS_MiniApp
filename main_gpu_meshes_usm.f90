PROGRAM LOOP3D
   ! This program computes the velocity flux on a list of meshes,
   ! using the "point-to-mesh" feature,
   ! with unified shared memory allocations.
   USE OMP_LIB
   IMPLICIT NONE
   
   ! Miscellaneous declarations
   INTEGER :: ISTEP
   INTEGER :: NM
   INTEGER :: N_MESHES
   INTEGER :: LOWER_MESH_INDEX, UPPER_MESH_INDEX
   INTEGER, PARAMETER :: IBAR = 256, JBAR = 256, KBAR =256
   INTEGER, PARAMETER :: IBP1 = IBAR+1, JBP1 = JBAR+1, KBP1 = KBAR+1
   INTEGER, PARAMETER :: NUM_TIME_STEPS = 1
   INTEGER, PARAMETER :: EB = SELECTED_REAL_KIND(12)
   INTEGER, PARAMETER :: NEDGE = 12
   REAL(EB), PARAMETER :: FOTH = 4.0_EB/3.0_EB
   
   REAL(EB), POINTER, DIMENSION(:,:,:) :: DP, RHOP, UU, VV, WW, OMY, OMX, OMZ, TXZ, TXY, TYZ
   REAL(EB), POINTER, DIMENSION(:) :: X, Y, Z, RDXN, RDYN, RDZN, RDX, RDY, RDZ, RHO_0
   REAL(EB), POINTER, DIMENSION(:,:,:) :: D, U, V, W, RHO, MU, FVX, FVY, FVZ, WORK1, WORK2, WORK3, WORK4, WORK5, WORK6
   INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX
                                          
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: GX, GY, GZ
   REAL(EB) :: MUX,MUY,MUZ,UP,UM,VP,VM,WP,WM,VTRM,OMXP,OMXM,OMYP,OMYM,OMZP,OMZM,TXYP,TXYM,TXZP,TXZM,TYZP,TYZM, &
               DTXYDY,DTXZDZ,DTYZDZ,DTXYDX,DTXZDX,DTYZDY, &
               DUDX,DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY, &
               VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX, &
               RRHO,TXXP,TXXM,TYYP,TYYM,TZZP,TZZM,DTXXDX,DTYYDY,DTZZDZ,T_NOW,T_END
   INTEGER :: I,J,K,IEXP,IEXM,IEYP,IEYM,IEZP,IEZM,IC,IE,MAX_EDGE,NT
   CHARACTER(LEN=50) :: FILENAME
   
   TYPE CELL_TYPE
      INTEGER :: EDGE_INDEX(NEDGE)=0
   END TYPE CELL_TYPE
   TYPE(CELL_TYPE), ALLOCATABLE, DIMENSION(:) :: CELL
   
   TYPE EDGE_TYPE
      REAL(EB) :: OMEGA(-2:2)=-1.E6_EB
      REAL(EB) :: TAU(-2:2)=-1.E6_EB
   END TYPE EDGE_TYPE
   TYPE(EDGE_TYPE), ALLOCATABLE, DIMENSION(:) :: EDGE

   TYPE MESH_TYPE
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: U
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: V
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: W
       REAL(EB), ALLOCATABLE, DIMENSION(:) :: X
       REAL(EB), ALLOCATABLE, DIMENSION(:) :: Y
       REAL(EB), ALLOCATABLE, DIMENSION(:) :: Z
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: D
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: RHO
       REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHO_0
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: WORK1
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: WORK2
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: WORK3
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: WORK4
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: WORK5
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: WORK6
       REAL(EB), ALLOCATABLE, DIMENSION(:) :: RDXN
       REAL(EB), ALLOCATABLE, DIMENSION(:) :: RDYN
       REAL(EB), ALLOCATABLE, DIMENSION(:) :: RDZN
       REAL(EB), ALLOCATABLE, DIMENSION(:) :: RDX
       REAL(EB), ALLOCATABLE, DIMENSION(:) :: RDY
       REAL(EB), ALLOCATABLE, DIMENSION(:) :: RDZ
       INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: CELL_INDEX
       !TYPE(CELL_TYPE), ALLOCATABLE, DIMENSION(:) :: CELL
       !TYPE(EDGE_TYPE), ALLOCATABLE, DIMENSION(:) :: EDGE
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: MU
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FVX
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FVY
       REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FVZ
   END TYPE MESH_TYPE
   
   TYPE (MESH_TYPE), POINTER :: M
   TYPE (MESH_TYPE), ALLOCATABLE, TARGET, DIMENSION(:) :: MESHES

   ! Write out Starting:
   !$OMP PARALLEL
   !$OMP MASTER
   !$ NT = OMP_GET_NUM_THREADS()
   !$OMP END MASTER
   !$OMP BARRIER
   !$OMP END PARALLEL

   N_MESHES = 2
   LOWER_MESH_INDEX = 1
   UPPER_MESH_INDEX = 2
   ALLOCATE(MESHES(N_MESHES))

   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
       CALL ALLOCATE_MESHES(NM)
   END DO
   
   WRITE(FILENAME,'(A,I3,A,I3,A,I3,A)') 'loop3d_',IBAR,'GRID_',NT,'THR_',NUM_TIME_STEPS,'STEPS_GPU_V0.txt'
   WRITE(*,*) 'Starting Loop3D, out file: ',TRIM(FILENAME)
   OPEN(UNIT=10,FILE=TRIM(FILENAME),STATUS='UNKNOWN')
   WRITE(10,*) 'Number of devices=',OMP_GET_NUM_DEVICES()
   WRITE(10,*) 'Starting Loop3D'
   WRITE(10,*) 'IBAR=',IBAR,' JBAR=',JBAR,' KBAR=',KBAR,' OMP_NUM_THREADS=',NT
   
   ! gravity (not part of the mesh)
   ALLOCATE(GX(0:IBAR), GY(0:IBAR), GZ(0:IBAR))
   
   ! Initialize meshes:
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
       CALL INITIALIZE_MESHES(NM)
   END DO

   T_NOW = OMP_GET_WTIME()
   SIM_LOOP: DO ISTEP = 1, NUM_TIME_STEPS
       DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
         CALL LOOP3D_OMP_GPU(NM)
       END DO
   END DO SIM_LOOP
   T_END = OMP_GET_WTIME()
   
   
   WRITE(10,*) 'Time=',T_END-T_NOW
   DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
       WRITE(10,*) 'mean FVX (MESH ',NM,')=',SUM(MESHES(NM)%FVX(1:IBAR,1:JBAR,1:KBAR))/(IBAR*JBAR*KBAR)
       WRITE(10,*) 'mean FVY (MESH ',NM,')=',SUM(MESHES(NM)%FVY(1:IBAR,1:JBAR,1:KBAR))/(IBAR*JBAR*KBAR)
       WRITE(10,*) 'mean FVZ (MESH ',NM,')=',SUM(MESHES(NM)%FVZ(1:IBAR,1:JBAR,1:KBAR))/(IBAR*JBAR*KBAR)
   END DO
   WRITE(10,*) 'Ending Loop3D'
   CLOSE(10)
   WRITE(*,*) 'Loop3D done.'
   
   CONTAINS
   
   SUBROUTINE POINT_TO_MESH(NM)
       INTEGER, INTENT(IN) :: NM
       M => MESHES(NM)
       X => M%X
       Y => M%Y
       Z => M%Z
       U => M%U
       V => M%V
       W => M%W
       RDXN => M%RDXN
       RDYN => M%RDYN
       RDZN => M%RDZN
       RDX => M%RDX
       RDY => M%RDY
       RDZ => M%RDZ
       RHO => M%RHO
       RHO_0 => M%RHO_0
       D => M%D
       MU => M%MU
       FVX => M%FVX
       FVY => M%FVY
       FVZ => M%FVZ
       WORK1 => M%WORK1
       WORK2 => M%WORK2
       WORK3 => M%WORK3
       WORK4 => M%WORK4
       WORK5 => M%WORK5
       WORK6 => M%WORK6
       CELL_INDEX => M%CELL_INDEX

   END SUBROUTINE POINT_TO_MESH

   SUBROUTINE ALLOCATE_MESHES(NM)
       INTEGER, INTENT(IN) :: NM
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%X)
       ALLOCATE(MESHES(NM)%X(0:IBAR))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%Y)
       ALLOCATE(MESHES(NM)%Y(0:JBAR))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%Z)
       ALLOCATE(MESHES(NM)%Z(0:KBAR))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%RDXN)
       ALLOCATE(MESHES(NM)%RDXN(0:IBAR))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%RDYN)
       ALLOCATE(MESHES(NM)%RDYN(0:JBAR))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%RDZN)
       ALLOCATE(MESHES(NM)%RDZN(0:KBAR))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%RDX)
       ALLOCATE(MESHES(NM)%RDX(0:IBAR))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%RDY)
       ALLOCATE(MESHES(NM)%RDY(0:JBAR))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%RDZ)
       ALLOCATE(MESHES(NM)%RDZ(0:KBAR))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%U)
       ALLOCATE(MESHES(NM)%U(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%V)
       ALLOCATE(MESHES(NM)%V(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%W)
       ALLOCATE(MESHES(NM)%W(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%MU)
       ALLOCATE(MESHES(NM)%MU(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%D)
       ALLOCATE(MESHES(NM)%D(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%RHO)
       ALLOCATE(MESHES(NM)%RHO(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%WORK1)
       ALLOCATE(MESHES(NM)%WORK1(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%WORK2)
       ALLOCATE(MESHES(NM)%WORK2(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%WORK3)
       ALLOCATE(MESHES(NM)%WORK3(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%WORK4)
       ALLOCATE(MESHES(NM)%WORK4(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%WORK5)
       ALLOCATE(MESHES(NM)%WORK5(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%WORK6)
       ALLOCATE(MESHES(NM)%WORK6(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%RHO_0)
       ALLOCATE(MESHES(NM)%RHO_0(0:KBAR))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%FVX)
       ALLOCATE(MESHES(NM)%FVX(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%FVY)
       ALLOCATE(MESHES(NM)%FVY(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%FVZ)
       ALLOCATE(MESHES(NM)%FVZ(0:IBP1,0:JBP1,0:KBP1))
       !$OMP ALLOCATORS ALLOCATE(OMP_TARGET_SHARED_MEM_ALLOC:MESHES(NM)%CELL_INDEX)
       ALLOCATE(MESHES(NM)%CELL_INDEX(0:IBP1,0:JBP1,0:KBP1))
       !? Can we define the data type we want for the allocator
   END SUBROUTINE ALLOCATE_MESHES

   SUBROUTINE INITIALIZE_MESHES(NM)
       INTEGER, INTENT(IN) :: NM
       CALL POINT_TO_MESH(NM)

       DO I=0,IBAR
           X(I)    = REAL(I,EB)
           RDXN(I) = 1.0_EB
           RDX(I)  = 1.0_EB
        ENDDO
        DO J=0,JBAR
           Y(J)    = REAL(J,EB)
           RDYN(J) = 1.0_EB
           RDY(J)  = 1.0_EB
        ENDDO
        DO K=0,KBAR
           Z(K)    = REAL(K,EB)
           RDZN(K) = 1.0_EB
           RDZ(K)  = 1.0_EB
        ENDDO
        ! Cell Index, CELL and EDGE:
        CELL_INDEX = 0
        IC = 0
        DO K=1,KBAR
           DO J=1,JBAR
              DO I=1,IBAR
                 IF( .NOT. (ANY( K==(/1,KBAR/) ) .OR. ANY( J==(/1,JBAR/) ) .OR. ANY( I==(/1,IBAR/) )) ) CYCLE
                 IC = IC + 1
                 CELL_INDEX(I,J,K) = IC
              ENDDO
           ENDDO
        ENDDO
        IF(.NOT.ALLOCATED(CELL)) ALLOCATE(CELL(0:IC))
        IC = 0
        MAX_EDGE=-1
        DO K=1,KBAR
           DO J=1,JBAR
              DO I=1,IBAR
                 IF( .NOT. (ANY( K==(/1,KBAR/) ) .OR. ANY( J==(/1,JBAR/) ) .OR. ANY( I==(/1,IBAR/) )) ) CYCLE
                 IC = IC + 1
                 CELL(IC)%EDGE_INDEX(1)  = CELL_INDEX(I  ,J  ,K  ) + 1
                 CELL(IC)%EDGE_INDEX(2)  = CELL_INDEX(I+1,J  ,K  ) + 1
                 CELL(IC)%EDGE_INDEX(3)  = CELL_INDEX(I  ,J  ,K  ) + 2
                 CELL(IC)%EDGE_INDEX(4)  = CELL_INDEX(I  ,J+1,K  ) + 2
                 CELL(IC)%EDGE_INDEX(5)  = CELL_INDEX(I  ,J  ,K  ) + 3
                 CELL(IC)%EDGE_INDEX(6)  = CELL_INDEX(I  ,J  ,K-1) + 3
                 CELL(IC)%EDGE_INDEX(7)  = CELL_INDEX(I  ,J  ,K  ) + 4
                 CELL(IC)%EDGE_INDEX(8)  = CELL_INDEX(I  ,J+1,K  ) + 4
                 CELL(IC)%EDGE_INDEX(9)  = CELL_INDEX(I  ,J  ,K  ) + 5
                 CELL(IC)%EDGE_INDEX(10) = CELL_INDEX(I  ,J  ,K-1) + 5
                 CELL(IC)%EDGE_INDEX(11) = CELL_INDEX(I  ,J  ,K  ) + 6
                 CELL(IC)%EDGE_INDEX(12) = CELL_INDEX(I+1,J  ,K  ) + 6
                 DO IE=1,NEDGE
                    MAX_EDGE = MAX(MAX_EDGE,CELL(IC)%EDGE_INDEX(IE))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        IF(.NOT.ALLOCATED(EDGE)) ALLOCATE(EDGE(0:MAX_EDGE))
        DO IE=1,MAX_EDGE,2
           EDGE(IE)%OMEGA = 1.5E-4_EB
           EDGE(IE)%TAU   = 2.5E-4_EB
        ENDDO
        RHO   = 1.19_EB
        RHO_0 = 1.19_EB
        D     = 0.0015_EB
        MU    = 0.0019_EB
        WORK1 = 0.0_EB; WORK2 = 0.0_EB; WORK3 = 0.0_EB; WORK4 = 0.0_EB; WORK5 = 0.0_EB; WORK6 = 0.0_EB
        GX(:) = 0.0_EB; GY(:) = 0.0_EB; GZ(:) = 1.0_EB
        ! U, V, W:
        DO K=0,KBAR
           DO J=0,JBAR
              DO I=0,IBAR
                 ! Some Trig functions:
                 U(I,J,K) =  SIN(X(I))*COS(Y(J))*COS(Z(K))
                 V(I,J,K) = -COS(X(I))*SIN(Y(J))*COS(Z(K))
                 W(I,J,K) =  COS(X(I))*COS(Y(J))*SIN(Z(K))
              ENDDO
           ENDDO
        ENDDO
        ! Compute Tau OMG:

        UU => U
        VV => V
        WW => W
        DP => D
        RHOP => RHO
        TXY => WORK1
        TXZ => WORK2
        TYZ => WORK3
        OMX => WORK4
        OMY => WORK5
        OMZ => WORK6
        DO K=0,KBAR
           DO J=0,JBAR
              DO I=0,IBAR
                 DUDY = RDYN(J)*(UU(I,J+1,K)-UU(I,J,K))
                 DVDX = RDXN(I)*(VV(I+1,J,K)-VV(I,J,K))
                 DUDZ = RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
                 DWDX = RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
                 DVDZ = RDZN(K)*(VV(I,J,K+1)-VV(I,J,K))
                 DWDY = RDYN(J)*(WW(I,J+1,K)-WW(I,J,K))
                 OMX(I,J,K) = DWDY - DVDZ
                 OMY(I,J,K) = DUDZ - DWDX
                 OMZ(I,J,K) = DVDX - DUDY
                 MUX = 0.25_EB*(MU(I,J+1,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I,J+1,K+1))
                 MUY = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I+1,J,K+1))
                 MUZ = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J+1,K)+MU(I+1,J+1,K))
                 TXY(I,J,K) = MUZ*(DVDX + DUDY)
                 TXZ(I,J,K) = MUY*(DUDZ + DWDX)
                 TYZ(I,J,K) = MUX*(DVDZ + DWDY)
              ENDDO
           ENDDO
        ENDDO
   END SUBROUTINE INITIALIZE_MESHES

   SUBROUTINE LOOP3D_OMP_GPU(NM)
      INTEGER, INTENT(IN) :: NM
      CALL POINT_TO_MESH(NM)
      ! Compute x-direction flux term FVX
      !$OMP TARGET DATA USE_DEVICE_ADDR(UU,VV,WW,DP,RHOP,TXY,TXZ,TYZ,OMX,OMY,OMZ) &
      !$OMP USE_DEVICE_ADDR(X,Y,Z,RDXN,RDYN,RDZN,RDX,RDY,RDZ,MU,RHO_0) &
      !$OMP USE_DEVICE_ADDR(FVX,FVY,FVZ,CELL_INDEX)

      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO NUM_THREADS(128) COLLAPSE(3) 

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=0,IBAR
               WP    = WW(I,J,K)   + WW(I+1,J,K)
               WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
               VP    = VV(I,J,K)   + VV(I+1,J,K)
               VM    = VV(I,J-1,K) + VV(I+1,J-1,K)
               OMYP  = OMY(I,J,K)
               OMYM  = OMY(I,J,K-1)
               OMZP  = OMZ(I,J,K)
               OMZM  = OMZ(I,J-1,K)
               TXZP  = TXZ(I,J,K)
               TXZM  = TXZ(I,J,K-1)
               TXYP  = TXY(I,J,K)
               TXYM  = TXY(I,J-1,K)
               IC    = CELL_INDEX(I,J,K)
               IEYP  = CELL(IC)%EDGE_INDEX(8)
               IEYM  = CELL(IC)%EDGE_INDEX(6)
               IEZP  = CELL(IC)%EDGE_INDEX(12)
               IEZM  = CELL(IC)%EDGE_INDEX(10)
               IF (EDGE(IEYP)%OMEGA(-1)>-1.E5_EB) THEN
                  OMYP = EDGE(IEYP)%OMEGA(-1)
                  TXZP = EDGE(IEYP)%TAU(-1)
               ENDIF
               IF (EDGE(IEYM)%OMEGA( 1)>-1.E5_EB) THEN
                  OMYM = EDGE(IEYM)%OMEGA( 1)
                  TXZM = EDGE(IEYM)%TAU( 1)
               ENDIF
               IF (EDGE(IEZP)%OMEGA(-2)>-1.E5_EB) THEN
                  OMZP = EDGE(IEZP)%OMEGA(-2)
                  TXYP = EDGE(IEZP)%TAU(-2)
               ENDIF
               IF (EDGE(IEZM)%OMEGA( 2)>-1.E5_EB) THEN
                  OMZM = EDGE(IEZM)%OMEGA( 2)
                  TXYM = EDGE(IEZM)%TAU( 2)
               ENDIF
               WOMY  = WP*OMYP + WM*OMYM
               VOMZ  = VP*OMZP + VM*OMZM
               RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
               DVDY  = (VV(I+1,J,K)-VV(I+1,J-1,K))*RDY(J)
               DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*RDZ(K)
               TXXP  = MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2._EB*(DVDY+DWDZ) )
               DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
               DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
               TXXM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DVDY+DWDZ) )
               DTXXDX= RDXN(I)*(TXXP-TXXM)
               DTXYDY= RDY(J) *(TXYP-TXYM)
               DTXZDZ= RDZ(K) *(TXZP-TXZM)
               VTRM  = DTXXDX + DTXYDY + DTXZDZ
               FVX(I,J,K) = 0.25_EB*(WOMY - VOMZ) - GX(I) + RRHO*(GX(I)*RHO_0(K) - VTRM)
            ENDDO
         ENDDO
      ENDDO
   
      ! Compute y-direction flux term FVY
      
      !$OMP TEAMS DISTRIBUTE PARALLEL DO NUM_THREADS(128) COLLAPSE(3)

      DO K=1,KBAR
         DO J=0,JBAR
            DO I=1,IBAR
               UP    = UU(I,J,K)   + UU(I,J+1,K)
               UM    = UU(I-1,J,K) + UU(I-1,J+1,K)
               WP    = WW(I,J,K)   + WW(I,J+1,K)
               WM    = WW(I,J,K-1) + WW(I,J+1,K-1)
               OMXP  = OMX(I,J,K)
               OMXM  = OMX(I,J,K-1)
               OMZP  = OMZ(I,J,K)
               OMZM  = OMZ(I-1,J,K)
               TYZP  = TYZ(I,J,K)
               TYZM  = TYZ(I,J,K-1)
               TXYP  = TXY(I,J,K)
               TXYM  = TXY(I-1,J,K)
               IC    = CELL_INDEX(I,J,K)
               IEXP  = CELL(IC)%EDGE_INDEX(4)
               IEXM  = CELL(IC)%EDGE_INDEX(2)
               IEZP  = CELL(IC)%EDGE_INDEX(12)
               IEZM  = CELL(IC)%EDGE_INDEX(11)
               IF (EDGE(IEXP)%OMEGA(-2)>-1.E5_EB) THEN
                  OMXP = EDGE(IEXP)%OMEGA(-2)
                  TYZP = EDGE(IEXP)%TAU(-2)
               ENDIF
               IF (EDGE(IEXM)%OMEGA( 2)>-1.E5_EB) THEN
                  OMXM = EDGE(IEXM)%OMEGA( 2)
                  TYZM = EDGE(IEXM)%TAU( 2)
               ENDIF
               IF (EDGE(IEZP)%OMEGA(-1)>-1.E5_EB) THEN
                  OMZP = EDGE(IEZP)%OMEGA(-1)
                  TXYP = EDGE(IEZP)%TAU(-1)
               ENDIF
               IF (EDGE(IEZM)%OMEGA( 1)>-1.E5_EB) THEN
                  OMZM = EDGE(IEZM)%OMEGA( 1)
                  TXYM = EDGE(IEZM)%TAU( 1)
               ENDIF
               WOMX  = WP*OMXP + WM*OMXM
               UOMZ  = UP*OMZP + UM*OMZM
               RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
               DUDX  = (UU(I,J+1,K)-UU(I-1,J+1,K))*RDX(I)
               DWDZ  = (WW(I,J+1,K)-WW(I,J+1,K-1))*RDZ(K)
               TYYP  = MU(I,J+1,K)*( FOTH*DP(I,J+1,K) - 2._EB*(DUDX+DWDZ) )
               DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
               DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
               TYYM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DWDZ) )
               DTXYDX= RDX(I) *(TXYP-TXYM)
               DTYYDY= RDYN(J)*(TYYP-TYYM)
               DTYZDZ= RDZ(K) *(TYZP-TYZM)
               VTRM  = DTXYDX + DTYYDY + DTYZDZ
               FVY(I,J,K) = 0.25_EB*(UOMZ - WOMX) - GY(I) + RRHO*(GY(I)*RHO_0(K) - VTRM)
            ENDDO
         ENDDO
      ENDDO
   
      ! Compute z-direction flux term FVZ
      
      !$OMP TEAMS DISTRIBUTE PARALLEL DO NUM_THREADS(128) COLLAPSE(3)

      DO K=0,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               UP    = UU(I,J,K)   + UU(I,J,K+1)
               UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
               VP    = VV(I,J,K)   + VV(I,J,K+1)
               VM    = VV(I,J-1,K) + VV(I,J-1,K+1)
               OMYP  = OMY(I,J,K)
               OMYM  = OMY(I-1,J,K)
               OMXP  = OMX(I,J,K)
               OMXM  = OMX(I,J-1,K)
               TXZP  = TXZ(I,J,K)
               TXZM  = TXZ(I-1,J,K)
               TYZP  = TYZ(I,J,K)
               TYZM  = TYZ(I,J-1,K)
               IC    = CELL_INDEX(I,J,K)
               IEXP  = CELL(IC)%EDGE_INDEX(4)
               IEXM  = CELL(IC)%EDGE_INDEX(3)
               IEYP  = CELL(IC)%EDGE_INDEX(8)
               IEYM  = CELL(IC)%EDGE_INDEX(7)
               IF (EDGE(IEXP)%OMEGA(-1)>-1.E5_EB) THEN
                  OMXP = EDGE(IEXP)%OMEGA(-1)
                  TYZP = EDGE(IEXP)%TAU(-1)
               ENDIF
               IF (EDGE(IEXM)%OMEGA( 1)>-1.E5_EB) THEN
                  OMXM = EDGE(IEXM)%OMEGA( 1)
                  TYZM = EDGE(IEXM)%TAU( 1)
               ENDIF
               IF (EDGE(IEYP)%OMEGA(-2)>-1.E5_EB) THEN
                  OMYP = EDGE(IEYP)%OMEGA(-2)
                  TXZP = EDGE(IEYP)%TAU(-2)
               ENDIF
               IF (EDGE(IEYM)%OMEGA( 2)>-1.E5_EB) THEN
                  OMYM = EDGE(IEYM)%OMEGA( 2)
                  TXZM = EDGE(IEYM)%TAU( 2)
               ENDIF
               UOMY  = UP*OMYP + UM*OMYM
               VOMX  = VP*OMXP + VM*OMXM
               RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
               DUDX  = (UU(I,J,K+1)-UU(I-1,J,K+1))*RDX(I)
               DVDY  = (VV(I,J,K+1)-VV(I,J-1,K+1))*RDY(J)
               TZZP  = MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2._EB*(DUDX+DVDY) )
               DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
               DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
               TZZM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DVDY) )
               DTXZDX= RDX(I) *(TXZP-TXZM)
               DTYZDY= RDY(J) *(TYZP-TYZM)
               DTZZDZ= RDZN(K)*(TZZP-TZZM)
               VTRM  = DTXZDX + DTYZDY + DTZZDZ
               FVZ(I,J,K) = 0.25_EB*(VOMX - UOMY) - GZ(I) + RRHO*(GZ(I)*0.5_EB*(RHO_0(K)+RHO_0(K+1)) - VTRM)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END TARGET DATA
      END SUBROUTINE LOOP3D_OMP_GPU
   
   END PROGRAM LOOP3D
   
   
   
   
   
