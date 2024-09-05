PROGRAM LOOP3D
    ! This program computes the velocity flux over a single mesh,
    ! using a custom mapper to define how the mesh variables are to be
    ! automatically transferred to the device.
    USE OMP_LIB
    IMPLICIT NONE
    
    ! Miscellaneous declarations
    INTEGER :: ISTEP
    INTEGER, PARAMETER :: NUM_TIME_STEPS = 100
    INTEGER, PARAMETER :: EB = SELECTED_REAL_KIND(12)
    INTEGER, PARAMETER :: NEDGE = 12
    REAL(EB), PARAMETER :: FOTH = 4.0_EB/3.0_EB
    
    REAL(EB), ALLOCATABLE, DIMENSION(:) :: GX, GY, GZ
    REAL(EB) :: MUX,MUY,MUZ,UP,UM,VP,VM,WP,WM,VTRM,OMXP,OMXM,OMYP,OMYM,OMZP,OMZM,TXYP,TXYM,TXZP,TXZM,TYZP,TYZM, &
                DTXYDY,DTXZDZ,DTYZDZ,DTXYDX,DTXZDX,DTYZDY, &
                DUDX,DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY, &
                VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX, &
                RRHO,TXXP,TXXM,TYYP,TYYM,TZZP,TZZM,DTXXDX,DTYYDY,DTZZDZ,T_NOW,T_END
    INTEGER :: I,J,K,IEXP,IEXM,IEYP,IEYM,IEZP,IEZM,IC,IC1,IC2,IE,MAX_EDGE,NT
    CHARACTER(LEN=50) :: FILENAME
    
    TYPE CELL_TYPE
       INTEGER :: EDGE_INDEX(NEDGE)=0
    END TYPE CELL_TYPE
    
    TYPE EDGE_TYPE
       REAL(EB) :: OMEGA(-2:2)=-1.E6_EB
       REAL(EB) :: TAU(-2:2)=-1.E6_EB
    END TYPE EDGE_TYPE

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
        TYPE(CELL_TYPE), ALLOCATABLE, DIMENSION(:) :: CELL
        TYPE(EDGE_TYPE), ALLOCATABLE, DIMENSION(:) :: EDGE
        REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: MU
        REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FVX
        REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FVY
        REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FVZ
        INTEGER :: IBAR
        INTEGER :: JBAR
        INTEGER :: KBAR
        INTEGER :: IBP1
        INTEGER :: JBP1
        INTEGER :: KBP1
    END TYPE MESH_TYPE
    
    ! target allocatable array of meshes
    TYPE (MESH_TYPE) :: M
    !TYPE (MESH_TYPE), ALLOCATABLE, TARGET, DIMENSION(:) :: MESHES

    !$OMP DECLARE MAPPER (MESH_MAP : MESH_TYPE :: M) &
    !$OMP MAP(TOFROM: M%X(0:M%IBAR)) &
    !$OMP MAP(TOFROM: M%Y(0:M%JBAR)) &
    !$OMP MAP(TOFROM: M%Z(0:M%KBAR)) &
    !$OMP MAP(TOFROM: M%RDXN(0:M%IBAR)) &
    !$OMP MAP(TOFROM: M%RDYN(0:M%JBAR)) &
    !$OMP MAP(TOFROM: M%RDZN(0:M%KBAR)) &
    !$OMP MAP(TOFROM: M%RDX(0:M%IBAR)) &
    !$OMP MAP(TOFROM: M%RDY(0:M%JBAR)) &
    !$OMP MAP(TOFROM: M%RDZ(0:M%KBAR)) &
    !$OMP MAP(TOFROM: M%U(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%V(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%W(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%MU(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%D(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%RHO(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%WORK1(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%WORK2(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%WORK3(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%WORK4(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%WORK5(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%WORK6(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%RHO_0(0:M%KBAR)) &
    !$OMP MAP(TOFROM: M%FVX(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%FVY(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%FVZ(0:M%IBP1,0:M%JBP1,0:M%KBP1)) &
    !$OMP MAP(TOFROM: M%CELL_INDEX(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    
    ! Write out Starting:
    !$OMP PARALLEL
    !$OMP MASTER
    !$ NT = OMP_GET_NUM_THREADS()
    !$OMP END MASTER
    !$OMP BARRIER
    !$OMP END PARALLEL

    ! Initialize the number of mesh cells
    M%IBAR = 256
    M%JBAR = 256
    M%KBAR = 256
    M%IBP1 = M%IBAR+1
    M%JBP1 = M%JBAR+1
    M%KBP1 = M%KBAR+1
    
    WRITE(FILENAME,'(A,I3,A,I3,A,I3,A)') 'loop3d_',M%IBAR,'GRID_',NT,'THR_',NUM_TIME_STEPS,'STEPS_GPU_V0.txt'
    WRITE(*,*) 'Starting Loop3D, out file: ',TRIM(FILENAME)
    OPEN(UNIT=10,FILE=TRIM(FILENAME),STATUS='UNKNOWN')
    WRITE(10,*) 'Number of devices=',OMP_GET_NUM_DEVICES()
    WRITE(10,*) 'Starting Loop3D'
    WRITE(10,*) 'IBAR=',M%IBAR,' JBAR=',M%JBAR,' KBAR=',M%KBAR,' OMP_NUM_THREADS=',NT
    
    ! gravity (not part of the mesh)
    ALLOCATE(GX(0:M%IBAR), GY(0:M%IBAR), GZ(0:M%IBAR))

    ! pointer to mesh 1, point to mesh 2, rinse and repeat allocation, just like in FDS. meshes(NM).LOWER_MESH_INDEX to UPPER_MESH_INDEX
    ! Allocate mesh vars in CPU:
    ALLOCATE(M%X(0:M%IBAR))
    ALLOCATE(M%Y(0:M%JBAR))
    ALLOCATE(M%Z(0:M%KBAR))
    ALLOCATE(M%RDXN(0:M%IBAR))
    ALLOCATE(M%RDYN(0:M%JBAR))
    ALLOCATE(M%RDZN(0:M%KBAR))
    ALLOCATE(M%RDX(0:M%IBAR))
    ALLOCATE(M%RDY(0:M%JBAR))
    ALLOCATE(M%RDZ(0:M%KBAR))
    ALLOCATE(M%U(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%V(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%W(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%MU(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%D(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%RHO(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%WORK1(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%WORK2(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%WORK3(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%WORK4(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%WORK5(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%WORK6(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%RHO_0(0:M%KBAR))
    ALLOCATE(M%FVX(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%FVY(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%FVZ(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    ALLOCATE(M%CELL_INDEX(0:M%IBP1,0:M%JBP1,0:M%KBP1))
    
    ! Initialize:
    DO I=0,M%IBAR
       M%X(I)    = REAL(I,EB)
       M%RDXN(I) = 1.0_EB
       M%RDX(I)  = 1.0_EB
    ENDDO
    DO J=0,M%JBAR
       M%Y(J)    = REAL(J,EB)
       M%RDYN(J) = 1.0_EB
       M%RDY(J)  = 1.0_EB
    ENDDO
    DO K=0,M%KBAR
       M%Z(K)    = REAL(K,EB)
       M%RDZN(K) = 1.0_EB
       M%RDZ(K)  = 1.0_EB
    ENDDO
    
    ! Cell Index, CELL and EDGE:
    M%CELL_INDEX = 0
    IC = 0
    DO K=1,M%KBAR
       DO J=1,M%JBAR
          DO I=1,M%IBAR
             IF( .NOT. (ANY( K==(/1,M%KBAR/) ) .OR. ANY( J==(/1,M%JBAR/) ) .OR. ANY( I==(/1,M%IBAR/) )) ) CYCLE
             IC = IC + 1
             M%CELL_INDEX(I,J,K) = IC
          ENDDO
       ENDDO
    ENDDO
    ALLOCATE(M%CELL(0:IC))
    
    IC = 0
    MAX_EDGE=-1
    DO K=1,M%KBAR
       DO J=1,M%JBAR
          DO I=1,M%IBAR
             IF( .NOT. (ANY( K==(/1,M%KBAR/) ) .OR. ANY( J==(/1,M%JBAR/) ) .OR. ANY( I==(/1,M%IBAR/) )) ) CYCLE
             IC = IC + 1
             M%CELL(IC)%EDGE_INDEX(1)  = M%CELL_INDEX(I  ,J  ,K  ) + 1
             M%CELL(IC)%EDGE_INDEX(2)  = M%CELL_INDEX(I+1,J  ,K  ) + 1
             M%CELL(IC)%EDGE_INDEX(3)  = M%CELL_INDEX(I  ,J  ,K  ) + 2
             M%CELL(IC)%EDGE_INDEX(4)  = M%CELL_INDEX(I  ,J+1,K  ) + 2
             M%CELL(IC)%EDGE_INDEX(5)  = M%CELL_INDEX(I  ,J  ,K  ) + 3
             M%CELL(IC)%EDGE_INDEX(6)  = M%CELL_INDEX(I  ,J  ,K-1) + 3
             M%CELL(IC)%EDGE_INDEX(7)  = M%CELL_INDEX(I  ,J  ,K  ) + 4
             M%CELL(IC)%EDGE_INDEX(8)  = M%CELL_INDEX(I  ,J+1,K  ) + 4
             M%CELL(IC)%EDGE_INDEX(9)  = M%CELL_INDEX(I  ,J  ,K  ) + 5
             M%CELL(IC)%EDGE_INDEX(10) = M%CELL_INDEX(I  ,J  ,K-1) + 5
             M%CELL(IC)%EDGE_INDEX(11) = M%CELL_INDEX(I  ,J  ,K  ) + 6
             M%CELL(IC)%EDGE_INDEX(12) = M%CELL_INDEX(I+1,J  ,K  ) + 6
             DO IE=1,NEDGE
                MAX_EDGE = MAX(MAX_EDGE,M%CELL(IC)%EDGE_INDEX(IE))
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ALLOCATE(M%EDGE(0:MAX_EDGE))
    DO IE=1,MAX_EDGE,2
       M%EDGE(IE)%OMEGA = 1.5E-4_EB
       M%EDGE(IE)%TAU   = 2.5E-4_EB
    ENDDO
    
    M%RHO   = 1.19_EB
    M%RHO_0 = 1.19_EB
    M%D     = 0.0015_EB
    M%MU    = 0.0019_EB
    M%WORK1 = 0.0_EB; M%WORK2 = 0.0_EB; M%WORK3 = 0.0_EB; M%WORK4 = 0.0_EB; M%WORK5 = 0.0_EB; M%WORK6 = 0.0_EB
    GX(:) = 0.0_EB; GY(:) = 0.0_EB; GZ(:) = 1.0_EB
    
    ! U, V, W:
    DO K=0,M%KBAR
       DO J=0,M%JBAR
          DO I=0,M%IBAR
             ! Some Trig functions:
             M%U(I,J,K) =  SIN(M%X(I))*COS(M%Y(J))*COS(M%Z(K))
             M%V(I,J,K) = -COS(M%X(I))*SIN(M%Y(J))*COS(M%Z(K))
             M%W(I,J,K) =  COS(M%X(I))*COS(M%Y(J))*SIN(M%Z(K))
          ENDDO
       ENDDO
    ENDDO
    
    ! Compute Tau OMG:
    DO K=0,M%KBAR
       DO J=0,M%JBAR
          DO I=0,M%IBAR
             DUDY = M%RDYN(J)*(M%U(I,J+1,K)-M%U(I,J,K))
             DVDX = M%RDXN(I)*(M%V(I+1,J,K)-M%V(I,J,K))
             DUDZ = M%RDZN(K)*(M%U(I,J,K+1)-M%U(I,J,K))
             DWDX = M%RDXN(I)*(M%W(I+1,J,K)-M%W(I,J,K))
             DVDZ = M%RDZN(K)*(M%V(I,J,K+1)-M%V(I,J,K))
             DWDY = M%RDYN(J)*(M%W(I,J+1,K)-M%W(I,J,K))
             M%WORK4(I,J,K) = DWDY - DVDZ
             M%WORK5(I,J,K) = DUDZ - DWDX
             M%WORK6(I,J,K) = DVDX - DUDY
             MUX = 0.25_EB*(M%MU(I,J+1,K)+M%MU(I,J,K)+M%MU(I,J,K+1)+M%MU(I,J+1,K+1))
             MUY = 0.25_EB*(M%MU(I+1,J,K)+M%MU(I,J,K)+M%MU(I,J,K+1)+M%MU(I+1,J,K+1))
             MUZ = 0.25_EB*(M%MU(I+1,J,K)+M%MU(I,J,K)+M%MU(I,J+1,K)+M%MU(I+1,J+1,K))
             M%WORK1(I,J,K) = MUZ*(DVDX + DUDY)
             M%WORK2(I,J,K) = MUY*(DUDZ + DWDX)
             M%WORK3(I,J,K) = MUX*(DVDZ + DWDY)
          ENDDO
       ENDDO
    ENDDO
    
    
    T_NOW = OMP_GET_WTIME()
    SIM_LOOP: DO ISTEP = 1, NUM_TIME_STEPS
       CALL LOOP3D_OMP_GPU()
    END DO SIM_LOOP
    T_END = OMP_GET_WTIME()
    
    
    WRITE(10,*) 'Time=',T_END-T_NOW
    WRITE(10,*) 'mean FVX =',SUM(M%FVX(1:M%IBAR,1:M%JBAR,1:M%KBAR))/(M%IBAR*M%JBAR*M%KBAR)
    WRITE(10,*) 'mean FVY =',SUM(M%FVY(1:M%IBAR,1:M%JBAR,1:M%KBAR))/(M%IBAR*M%JBAR*M%KBAR)
    WRITE(10,*) 'mean FVZ =',SUM(M%FVZ(1:M%IBAR,1:M%JBAR,1:M%KBAR))/(M%IBAR*M%JBAR*M%KBAR)
    WRITE(10,*) 'Ending Loop3D'
    CLOSE(10)
    WRITE(*,*) 'Loop3D done.'
    
    CONTAINS
    
    SUBROUTINE LOOP3D_OMP_GPU()
       
       ! Compute x-direction flux term FVX
       !$OMP TARGET DATA

       !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO NUM_THREADS(128) COLLAPSE(3)
       DO K=1,M%KBAR
          DO J=1,M%JBAR
             DO I=0,M%IBAR
                WP    = M%W(I,J,K)   + M%W(I+1,J,K)
                WM    = M%W(I,J,K-1) + M%W(I+1,J,K-1)
                VP    = M%V(I,J,K)   + M%V(I+1,J,K)
                VM    = M%V(I,J-1,K) + M%V(I+1,J-1,K)
                OMYP  = M%WORK5(I,J,K)
                OMYM  = M%WORK5(I,J,K-1)
                OMZP  = M%WORK6(I,J,K)
                OMZM  = M%WORK6(I,J-1,K)
                TXZP  = M%WORK2(I,J,K)
                TXZM  = M%WORK2(I,J,K-1)
                TXYP  = M%WORK1(I,J,K)
                TXYM  = M%WORK1(I,J-1,K)
                IC    = M%CELL_INDEX(I,J,K)
                IEYP  = M%CELL(IC)%EDGE_INDEX(8)
                IEYM  = M%CELL(IC)%EDGE_INDEX(6)
                IEZP  = M%CELL(IC)%EDGE_INDEX(12)
                IEZM  = M%CELL(IC)%EDGE_INDEX(10)
                IF (M%EDGE(IEYP)%OMEGA(-1)>-1.E5_EB) THEN
                   OMYP = M%EDGE(IEYP)%OMEGA(-1)
                   TXZP = M%EDGE(IEYP)%TAU(-1)
                ENDIF
                IF (M%EDGE(IEYM)%OMEGA( 1)>-1.E5_EB) THEN
                   OMYM = M%EDGE(IEYM)%OMEGA( 1)
                   TXZM = M%EDGE(IEYM)%TAU( 1)
                ENDIF
                IF (M%EDGE(IEZP)%OMEGA(-2)>-1.E5_EB) THEN
                   OMZP = M%EDGE(IEZP)%OMEGA(-2)
                   TXYP = M%EDGE(IEZP)%TAU(-2)
                ENDIF
                IF (M%EDGE(IEZM)%OMEGA( 2)>-1.E5_EB) THEN
                   OMZM = M%EDGE(IEZM)%OMEGA( 2)
                   TXYM = M%EDGE(IEZM)%TAU( 2)
                ENDIF
                WOMY  = WP*OMYP + WM*OMYM
                VOMZ  = VP*OMZP + VM*OMZM
                RRHO  = 2._EB/(M%RHO(I,J,K)+M%RHO(I+1,J,K))
                DVDY  = (M%V(I+1,J,K)-M%V(I+1,J-1,K))*M%RDY(J)
                DWDZ  = (M%W(I+1,J,K)-M%W(I+1,J,K-1))*M%RDZ(K)
                TXXP  = M%MU(I+1,J,K)*( FOTH*M%D(I+1,J,K) - 2._EB*(DVDY+DWDZ) )
                DVDY  = (M%V(I,J,K)-M%V(I,J-1,K))*M%RDY(J)
                DWDZ  = (M%W(I,J,K)-M%W(I,J,K-1))*M%RDZ(K)
                TXXM  = M%MU(I,J,K)  *( FOTH*M%D(I,J,K)   - 2._EB*(DVDY+DWDZ) )
                DTXXDX= M%RDXN(I)*(TXXP-TXXM)
                DTXYDY= M%RDY(J) *(TXYP-TXYM)
                DTXZDZ= M%RDZ(K) *(TXZP-TXZM)
                VTRM  = DTXXDX + DTXYDY + DTXZDZ
                M%FVX(I,J,K) = 0.25_EB*(WOMY - VOMZ) - GX(I) + RRHO*(GX(I)*M%RHO_0(K) - VTRM)
             ENDDO
          ENDDO
       ENDDO
    
       ! Compute y-direction flux term FVY
       
       !$OMP TEAMS DISTRIBUTE PARALLEL DO NUM_THREADS(128) COLLAPSE(3)
       DO K=1,M%KBAR
          DO J=0,M%JBAR
             DO I=1,M%IBAR
                UP    = M%U(I,J,K)   + M%U(I,J+1,K)
                UM    = M%U(I-1,J,K) + M%U(I-1,J+1,K)
                WP    = M%W(I,J,K)   + M%W(I,J+1,K)
                WM    = M%W(I,J,K-1) + M%W(I,J+1,K-1)
                OMXP  = M%WORK4(I,J,K)
                OMXM  = M%WORK4(I,J,K-1)
                OMZP  = M%WORK6(I,J,K)
                OMZM  = M%WORK6(I-1,J,K)
                TYZP  = M%WORK3(I,J,K)
                TYZM  = M%WORK3(I,J,K-1)
                TXYP  = M%WORK1(I,J,K)
                TXYM  = M%WORK1(I-1,J,K)
                IC    = M%CELL_INDEX(I,J,K)
                IEXP  = M%CELL(IC)%EDGE_INDEX(4)
                IEXM  = M%CELL(IC)%EDGE_INDEX(2)
                IEZP  = M%CELL(IC)%EDGE_INDEX(12)
                IEZM  = M%CELL(IC)%EDGE_INDEX(11)
                IF (M%EDGE(IEXP)%OMEGA(-2)>-1.E5_EB) THEN
                   OMXP = M%EDGE(IEXP)%OMEGA(-2)
                   TYZP = M%EDGE(IEXP)%TAU(-2)
                ENDIF
                IF (M%EDGE(IEXM)%OMEGA( 2)>-1.E5_EB) THEN
                   OMXM = M%EDGE(IEXM)%OMEGA( 2)
                   TYZM = M%EDGE(IEXM)%TAU( 2)
                ENDIF
                IF (M%EDGE(IEZP)%OMEGA(-1)>-1.E5_EB) THEN
                   OMZP = M%EDGE(IEZP)%OMEGA(-1)
                   TXYP = M%EDGE(IEZP)%TAU(-1)
                ENDIF
                IF (M%EDGE(IEZM)%OMEGA( 1)>-1.E5_EB) THEN
                   OMZM = M%EDGE(IEZM)%OMEGA( 1)
                   TXYM = M%EDGE(IEZM)%TAU( 1)
                ENDIF
                WOMX  = WP*OMXP + WM*OMXM
                UOMZ  = UP*OMZP + UM*OMZM
                RRHO  = 2._EB/(M%RHO(I,J,K)+M%RHO(I,J+1,K))
                DUDX  = (M%U(I,J+1,K)-M%U(I-1,J+1,K))*M%RDX(I)
                DWDZ  = (M%W(I,J+1,K)-M%W(I,J+1,K-1))*M%RDZ(K)
                TYYP  = M%MU(I,J+1,K)*( FOTH*M%D(I,J+1,K) - 2._EB*(DUDX+DWDZ) )
                DUDX  = (M%U(I,J,K)-M%U(I-1,J,K))*M%RDX(I)
                DWDZ  = (M%W(I,J,K)-M%W(I,J,K-1))*M%RDZ(K)
                TYYM  = M%MU(I,J,K)  *( FOTH*M%D(I,J,K)   - 2._EB*(DUDX+DWDZ) )
                DTXYDX= M%RDX(I) *(TXYP-TXYM)
                DTYYDY= M%RDYN(J)*(TYYP-TYYM)
                DTYZDZ= M%RDZ(K) *(TYZP-TYZM)
                VTRM  = DTXYDX + DTYYDY + DTYZDZ
                M%FVY(I,J,K) = 0.25_EB*(UOMZ - WOMX) - GY(I) + RRHO*(GY(I)*M%RHO_0(K) - VTRM)
             ENDDO
          ENDDO
       ENDDO
    
       ! Compute z-direction flux term FVZ
       
       !$OMP TEAMS DISTRIBUTE PARALLEL DO NUM_THREADS(128) COLLAPSE(3)
       DO K=0,M%KBAR
          DO J=1,M%JBAR
             DO I=1,M%IBAR
                UP    = M%U(I,J,K)   + M%U(I,J,K+1)
                UM    = M%U(I-1,J,K) + M%U(I-1,J,K+1)
                VP    = M%V(I,J,K)   + M%V(I,J,K+1)
                VM    = M%V(I,J-1,K) + M%V(I,J-1,K+1)
                OMYP  = M%WORK5(I,J,K)
                OMYM  = M%WORK5(I-1,J,K)
                OMXP  = M%WORK4(I,J,K)
                OMXM  = M%WORK4(I,J-1,K)
                TXZP  = M%WORK2(I,J,K)
                TXZM  = M%WORK2(I-1,J,K)
                TYZP  = M%WORK3(I,J,K)
                TYZM  = M%WORK3(I,J-1,K)
                IC    = M%CELL_INDEX(I,J,K)
                IEXP  = M%CELL(IC)%EDGE_INDEX(4)
                IEXM  = M%CELL(IC)%EDGE_INDEX(3)
                IEYP  = M%CELL(IC)%EDGE_INDEX(8)
                IEYM  = M%CELL(IC)%EDGE_INDEX(7)
                IF (M%EDGE(IEXP)%OMEGA(-1)>-1.E5_EB) THEN
                   OMXP = M%EDGE(IEXP)%OMEGA(-1)
                   TYZP = M%EDGE(IEXP)%TAU(-1)
                ENDIF
                IF (M%EDGE(IEXM)%OMEGA( 1)>-1.E5_EB) THEN
                   OMXM = M%EDGE(IEXM)%OMEGA( 1)
                   TYZM = M%EDGE(IEXM)%TAU( 1)
                ENDIF
                IF (M%EDGE(IEYP)%OMEGA(-2)>-1.E5_EB) THEN
                   OMYP = M%EDGE(IEYP)%OMEGA(-2)
                   TXZP = M%EDGE(IEYP)%TAU(-2)
                ENDIF
                IF (M%EDGE(IEYM)%OMEGA( 2)>-1.E5_EB) THEN
                   OMYM = M%EDGE(IEYM)%OMEGA( 2)
                   TXZM = M%EDGE(IEYM)%TAU( 2)
                ENDIF
                UOMY  = UP*OMYP + UM*OMYM
                VOMX  = VP*OMXP + VM*OMXM
                RRHO  = 2._EB/(M%RHO(I,J,K)+M%RHO(I,J,K+1))
                DUDX  = (M%U(I,J,K+1)-M%U(I-1,J,K+1))*M%RDX(I)
                DVDY  = (M%V(I,J,K+1)-M%V(I,J-1,K+1))*M%RDY(J)
                TZZP  = M%MU(I,J,K+1)*( FOTH*M%D(I,J,K+1) - 2._EB*(DUDX+DVDY) )
                DUDX  = (M%U(I,J,K)-M%U(I-1,J,K))*M%RDX(I)
                DVDY  = (M%V(I,J,K)-M%V(I,J-1,K))*M%RDY(J)
                TZZM  = M%MU(I,J,K)  *( FOTH*M%D(I,J,K)   - 2._EB*(DUDX+DVDY) )
                DTXZDX= M%RDX(I) *(TXZP-TXZM)
                DTYZDY= M%RDY(J) *(TYZP-TYZM)
                DTZZDZ= M%RDZN(K)*(TZZP-TZZM)
                VTRM  = DTXZDX + DTYZDY + DTZZDZ
                M%FVZ(I,J,K) = 0.25_EB*(VOMX - UOMY) - GZ(I) + RRHO*(GZ(I)*0.5_EB*(M%RHO_0(K)+M%RHO_0(K+1)) - VTRM)
             ENDDO
          ENDDO
       ENDDO
       !$OMP END TARGET DATA
       END SUBROUTINE LOOP3D_OMP_GPU
    
    END PROGRAM LOOP3D
    
    
    
    
    
