PROGRAM WSMTD_rw_2D
!---------------------------------------------------------------------------------------------------------------------!
!FORTRAN PROGRAM WRITTEN TO COMPUTE UNBIASED DISTRIBUTION (1D AND 2D)FROM TASS SIMULATION ALONG USER DEFINED          !
!COLLECTIVE COORDINATES..                                                                                             !
!ORIGINAL CODE WRITTEN BY Shalini Awasthi (ashalini@iitk.ac.in)                                                       !
!MODIFIED BY Rahul Verma (vrahul@iitk.ac.in)                                                                          !
!                                                                                                                     !
!kt0 = System Temeprature ; kt = Extended CV Temperature ; bias_fact = Biased Factor for MTD ; ct = ct Factor         !
!t_min = Minimum MD steps ; t_max = Maximum MD steps ; narg = Argumets ; ncv = Numer of CV ; v_baias = Total MTD bias !
!cv_mtd = m = MTD CV index ; u_cv = u = Umbrella CV index ; t_cv = t = Temeprature CV index                           !
!UCV = U cv argument ; MTD = MTD cv argument ; Prob_nD = Dimension of Unbiased Probability                            !
!CV_num = Probability is the Dimension of ; pfrqMD = Print Frequency argument ; w_cv Print Frequency in cvmdck File   !
!dtMTD  = Print Frequency argument for MTD bias added ; w_hill = Print Frequency in colvar File                       !
!gridmin = Minimum Grid Size = gridmax = Maximum Grid Size ; griddif = Grid Difference                                !
!width = Hill Width of Gaussian Bias in MTD ; ht = Hill Height of Gaussian Bias in MTD ; ht = MTD CV Displacement     !
!kb = Boltzman Constant in A.U. ; prob = 1D Probability ; prob_2D = 2D Probability                                    !
!---------------------------------------------------------------------------------------------------------------------!
USE GetSteps
USE GetFileName
USE Input_file
USE Error_msg
USE MTD_Unbais
USE MTD_Potential
USE US_Prob
USE US_MTD
USE US_TEMP
USE MeanForce
!USE wham_code
USE B_Spline
IMPLICIT NONE
REAL*8              :: dummy11,den,kt0,kt,bias_fact
REAL*8, ALLOCATABLE :: prob(:),fes(:),pcons(:),kcons(:)
REAL*8, ALLOCATABLE :: dummy(:,:,:),cv(:,:,:),prob_2D(:,:),prob_mtd(:,:,:)
REAL*8, ALLOCATABLE :: gridmin(:),gridmax(:),griddif(:),vbias(:,:),ct(:,:),norm(:)
INTEGER,ALLOCATABLE :: nbin(:),indx(:),t(:),t_cv(:)
INTEGER :: md_steps,mtd_steps,dummy1,i,j,t_min,t_max,narg
INTEGER :: i_md,ncv,w_hill,w_cv,prob_nD,cv_mtd,cv_us,cv_num(3)
INTEGER :: ii,jj,kk,u,m,ir,nr,ios
LOGICAL :: pmf,probT,spline,inpgrid,read_ct,read_vbias,max_step
CHARACTER*5   :: mtd,tool
CHARACTER*10  :: code_name
CHARACTER*120 :: arg 
CHARACTER*50  :: filename_loc
CHARACTER(LEN=50),ALLOCATABLE :: filename(:),filename_mtd(:,:)
REAL*8, PARAMETER :: kb = 1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: au_to_kcal = 627.51 
REAL*8, PARAMETER :: kj_to_kcal = 0.239006 
!-----------------------------------------------!
md_steps    = 9999999  ; mtd_steps  = 999999    !
kt0         = 300.D0   ; kt         = 300.D0    !
t_min       = 1        ; t_max      = 7000      !
w_hill      = 0        ; w_cv       = 0         !
pmf         = .FALSE.  ; inpgrid    = .FALSE.   !
read_ct     = .FALSE.  ; read_vbias = .FALSE.   !
probT       = .FALSE.  ; spline     = .FALSE.   !
narg        = IARGC()  ; bias_fact  = 1500.D0   !
max_step = .FALSE.                              !
!-----------------------------------------------!
PRINT*,"!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!+!"
PRINT*,"ARGUMENTS IN THE RUN FILE ARE CASE SENSITIVE"
!----------------------------------------------------!
DO i=1,narg
  CALL GETARG(i,arg)
  IF(INDEX(arg,'-T0').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)kt0
  ELSEIF(INDEX(arg,'-T').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)kt
  ELSEIF(INDEX(arg,'-prog_name').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)code_name
     IF(code_name .ne. 'PLUMED') THEN
     IF(code_name .ne. 'CPMD') THEN
     PRINT*,"ERROR : ONLY WORKS WITH CPMD or PLUMED"
     STOP ; ENDIF ; ENDIF
  ELSEIF(INDEX(arg,'-bias_fact').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)bias_fact
  ELSE IF(INDEX(arg,'-tmin').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)t_min
  ELSE IF(INDEX(arg,'-tmax').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)t_max
     max_step = .TRUE.
   ELSEIF(INDEX(arg, '-ncv') .NE. 0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)ncv
   ELSEIF(INDEX(arg, '-UCV') .NE. 0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)cv_us
   ELSEIF(INDEX(arg, '-MTD') .NE. 0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)mtd
   ELSEIF(INDEX(arg, '-MCV') .NE. 0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)cv_mtd
  ELSE IF(INDEX(arg,'-tool').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)tool
       IF (tool .eq. 'pmf') pmf=.TRUE.
       IF (tool .eq. 'probT') probT=.TRUE.
   ELSEIF(INDEX(arg, '-Prob_nD') .NE. 0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)prob_nD
       IF (tool .eq. 'pmf') STOP "ERROR: PLESE CHECK RUN FILE (may be 'probT' variables are active)"
       IF (prob_nD .gt. 3) STOP  "MORE THAN 3 DIMENSION IS NOT IMPLIMENTED"
   ELSEIF(INDEX(arg, '-CV_num') .NE. 0)THEN
         DO ii = 1,prob_nD
           CALL GETARG(i+ii,arg)
           READ(arg,*)cv_num(ii)
         ENDDO
  ELSE IF(INDEX(arg,'-grid').NE.0)THEN

ALLOCATE(gridmin(ncv),gridmax(ncv),griddif(ncv)) ; ALLOCATE(indx(ncv))

j = 0 
DO ii = 1,ncv
j = j + 1
      CALL GETARG(i+j,arg)
      READ(arg,*)gridmin(ii)
      CALL GETARG(i+j+1,arg)
      READ(arg,*)gridmax(ii)
      CALL GETARG(i+j+2,arg)
      READ(arg,*)griddif(ii)
j = j + 2 
!100 FORMAT (I4,2X,3F16.4)
!WRITE(*,100)ii,gridmin(ii),gridmax(ii),griddif(ii)
ENDDO
      inpgrid=.true.
  ELSE IF(INDEX(arg,'-pfrqMD').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)w_cv
  ELSE IF(INDEX(arg,'-dtMTD').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)w_hill
  ELSE IF(INDEX(arg,'-nr').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)nr
  ELSE IF(INDEX(arg,'-interpolate').NE.0)THEN
        spline=.TRUE.
  ELSE IF(INDEX(arg,'-read_ct').NE.0)THEN
        read_ct=.TRUE.
  ELSE IF(INDEX(arg,'-read_vbias').NE.0)THEN
        read_vbias=.TRUE.
  END IF
END DO
!=======================================================================!
IF(code_name .eq. 'CPMD')   OPEN(11,FILE='cvmdck_mtd',STATUS='unknown')
IF(code_name .eq. 'PLUMED') OPEN(11,FILE='COLVAR',STATUS='unknown')
!=======================================================================!
IF (mtd .eq. 'y') THEN
IF(code_name .eq. 'CPMD')   OPEN(12,FILE='parvar_mtd',STATUS='unknown')
IF (code_name .eq. 'CPMD')  OPEN(13,FILE='colvar_mtd',STATUS='unknown')
IF (code_name .eq. 'PLUMED')OPEN(13,FILE='HILLS',STATUS='unknown')
ENDIF
!=======================================================================!
IF (mtd .eq. 'y' .and. cv_mtd .eq. 0) STOP "***ERROR !! PLEASE SPECIFY METADYNAMICS CV COORDINATE INDEX"
WRITE(*,'(A)')'!------------------------------------------------------------------------------------'
   WRITE(*,'(A,2X,A)')'! MOLECULAR DYNAMICS PROGRAM USED FOR SAMPLING : ',code_name
   WRITE(*,'(A,2X,I2)')'! No. of Umbrella Windows                      = ',nr
!IF(probT)   WRITE(*,'(A,I10)')'! No: of MD  steps                            =',md_steps
   WRITE(*,'(A,I9)')'! No. of min MD  steps                         =',t_min
IF (max_step) WRITE(*,'(A,I9)')'! No. of max MD  steps                         =',t_max
   WRITE(*,'(A,2X,I3)')'! No. of CV                                    =',ncv
   WRITE(*,'(A,2X,I3)')'! Umbrella CV Column                           =',cv_us
IF (mtd .eq. 'n') THEN
   WRITE(*,'(A,4X,A)')'! Metadynamics is Enabled                      =','NO'
ELSEIF(mtd .eq. 'y') THEN
   WRITE(*,'(A,5X,A)')'! Metadynamics is Enabled                      =','YES'
   WRITE(*,'(A,2X,I5)')'! Metadynamics CV Column                       =',cv_mtd
ENDIF
   WRITE(*,'(A,F10.2)')'! Physical System Temperature T0 (K)           =',kt0
   WRITE(*,'(A,F11.2)')'! CV Temperature T (K)                         =',kt
IF (mtd .eq. 'y')  WRITE(*,'(A,F11.2)')'! Bias Factor (K)                              =',bias_fact
   WRITE(*,'(A,I6)')'! Print Freq. in cvmdck_mtd file               =',w_cv
   IF(mtd.eq.'y') WRITE(*,'(A,I7)')'! Freq. of Hill Update                         =',w_hill
IF (pmf) THEN
   WRITE(*,'(A,2X,A)')'! ####### FREE ENERGY WILL BE CALCULATED USING MEAN FORCE METHOD #######'
ELSEIF (probT)THEN
   WRITE(*,'(A,I6)')'! Dimension of Probability                     =',prob_nD
   WRITE(*,'(A,3I7)')'! Probability Index                            =',cv_num(1:prob_nD)
   WRITE(*,'(A,2X,A)')'!####### PROBABLITIES WILL BE GENERATED THROUGH UNBIASING #######'
ENDIF
WRITE(*,'(A)')'!------------------------------------------------------------------------------------'
ii = cv_num(1); jj = cv_num(2) ; kk = cv_num(3)
u  = cv_us    ; m  = cv_mtd  !! #### m = METADYNAMICS CV INDEX !! #### u = UMBRELLA CV INDEX
!===========================================================================================================!
IF(probT) THEN
WRITE(*,'(A85)')'=========================================================================================='! 
IF (ii .ne. u .and. jj .ne. u .and. kk .ne. u) THEN
WRITE(*,'(10X,A)') "!!SORRY!! PLEASE SPECIFY UMBRELLA CV INDEX PROPERLY" 
WRITE(*,'(A85)')'=========================================================================================='!
STOP ; ENDIF
IF (mtd .eq. 'y' .and. cv_mtd .eq. cv_us) THEN
!IF (ii .ne. m .and. jj .ne. m .and. kk .ne. m) THEN 
WRITE(*,'(10X,A)')"METADYNAMICS ENABLED, CAN'T UNBIAS PROBABILITY WITHOUT UNBIASING 'MTD'"
WRITE(*,'(A85)')'=========================================================================================='!
STOP
ENDIF ; ENDIF 
!-----------------------------------------------------------------------------------------------------!
ALLOCATE(filename(nr))          ; ALLOCATE(filename_mtd(2,nr))
ALLOCATE(cv(nr,ncv,md_steps))   ; ALLOCATE(dummy(nr,ncv,md_steps))
ALLOCATE(nbin(ncv))             ; ALLOCATE(vbias(nr,md_steps))
ALLOCATE(ct(nr,mtd_steps))      ; ALLOCATE(t(ncv))
ALLOCATE(kcons(nr))             ; ALLOCATE(pcons(nr))
ALLOCATE(t_cv(ncv))             ; ALLOCATE(norm(nr))
!-----------------------------------------------------------------------------------------------------!
OPEN(10,FILE='input.inp',STATUS='old',IOSTAT=ios,ACTION='read')
IF (ios .ne. 0) STOP "ERROR : FILE input.inp doesn't exist..!"

101 FORMAT (I10,10F16.6)
102 FORMAT (F10.6,10F16.6)

DO ir = 1,nr
!-----------------------------------------------------------------------------------------------------!
   IF(code_name .eq. 'CPMD') THEN
     READ(10,*)pcons(ir),kcons(ir)
      kcons(ir)=kcons(ir)*au_to_kcal
      READ(10,'(A)')filename(ir)
     IF(mtd .eq. 'y')THEN
       READ(10,'(A)')filename_mtd(1,ir)
       READ(10,'(A)')filename_mtd(2,ir)
     ENDIF
CALL get_filename('cv.dat_',filename_loc,ir)
OPEN(14,FILE=filename_loc,STATUS='unknown')
OPEN(11,FILE=filename(ir),STATUS='old',IOSTAT=ios)
CALL get_steps(11,md_steps)
DO i_md=1,md_steps
   READ(11,*)dummy1,dummy1,(dummy(ir,j,i_md), j=1,ncv),(cv(ir,j,i_md) ,j=1,ncv)
!   WRITE(*,*)dummy1,(dummy(j,i_md), j=1,ncv),(cv(j,i_md) ,j=1,ncv)
!-------------------------------------------------------------------------
DO j = 1,ncv
IF ( cv(ir,j,i_md) .lt. gridmin(j)) CALL cv_error (j)
IF ( cv(ir,j,i_md) .gt. gridmax(j)) CALL cv_error (j)
ENDDO
!-------------------------------------------------------------------------
   WRITE(14,101)dummy1,(cv(ir,j,i_md) ,j=1,ncv)
END DO
WRITE(6,'(A,I2,3X,A,I10)')'! No. of MD STEPS in umb:',ir,'=',md_steps
CLOSE(11);CLOSE(14)
!-------------------------------------------------------------------------
   ELSEIF (code_name .eq. 'PLUMED') THEN
     READ(10,*)pcons(ir),kcons(ir)
     kcons(ir)=kcons(ir)*kj_to_kcal
     READ(10,'(A)')filename(ir)
   IF(mtd .eq. 'y') THEN
     READ(10,'(A)')filename_mtd(1,ir)
   ENDIF

CALL get_filename('cv.dat_',filename_loc,ir)
OPEN(14,FILE=filename_loc,STATUS='unknown')
OPEN(11,FILE=filename(ir),STATUS='old',IOSTAT=ios)
CALL get_steps(11,md_steps)
DO i_md=1,md_steps
!   READ(11,*)dummy11,dummy11,cv(ir,u,i_md),dummy11,cv(ir,m,i_md)
   READ(11,*)dummy11,(cv(ir,j,i_md) ,j=1,ncv)
!-------------------------------------------------------------------------
DO j = 1,ncv
IF ( cv(ir,j,i_md) .lt. gridmin(j)) CALL cv_error (j)
IF ( cv(ir,j,i_md) .gt. gridmax(j)) CALL cv_error (j)
ENDDO
!-------------------------------------------------------------------------
!   WRITE(*,*)dummy11,(cv(j,i_md) ,j=1,ncv)
   WRITE(14,102)dummy11,(cv(ir,j,i_md) ,j=1,ncv)
END DO
WRITE(6,'(A,I2,3X,A,I10)')'! No. of MD STEPS in umb:',ir,'=',md_steps
CLOSE(11);CLOSE(14)
ENDIF
ENDDO
CLOSE(10)
!-----------------------------------------------------------------------------------------------------!
WRITE(*,'(A85)')'=========================================================================================='!

  DO i = 1,ncv 
   nbin(i) = NINT((gridmax(i)-gridmin(i))/griddif(i)) + 1
  ENDDO

j = 0
WRITE(*,'(9X,4A9)')'GRIDMIN','GRIDMAX','GRIDBIN','GRIDSIZE'
WRITE(*,'(A10,3F8.4,I10)')'US   COORD:', gridmin(u),gridmax(u),griddif(u),nbin(u)
IF (mtd .eq. 'y') THEN
WRITE(*,'(A10,3F8.4,I10)')'MTD  COORD:', gridmin(m),gridmax(m),griddif(m),nbin(m)
ENDIF
DO i = 1,ncv ; IF (i .ne. u .and. i .ne. m) THEN
WRITE(*,'(A10,3F8.4,2I10)')'TASS COORD:', gridmin(i),gridmax(i),griddif(i),nbin(i)
j = j + 1 ; t_cv(j) = i
ENDIF ; ENDDO
WRITE(*,'(A85)')'=========================================================================================='!
t(1:j) = t_cv(1:j) !; PRINT*,t(1:j)    !# t_cv TASS CV INDEX
!n1 = nbin(1) ;n2 = nbin(2) !;n3 = nbin(3) !;n4 = nbin(4) 
!---------------------------------------------------------------------------------------------------------------------------!
ALLOCATE(prob((nbin(u))))
IF (mtd .ne. 'y') m = t(1)
ALLOCATE(prob_2D(nbin(u),nbin(m)))
ALLOCATE(prob_mtd(nr,nbin(u),nbin(m)))
prob = 0.d0 ; prob_2D = 0.d0 ; prob_mtd = 0.d0
IF (jj .eq. 0 .and. kk .eq. 0) jj = 1 ; kk = 1
!ALLOCATE(prob_3D(nbin(ii),nbin(jj),nbin(kk)))
!---------------------------------------------------------------------------------------------------------------------------!
IF(mtd .eq. 'y') THEN

CALL mtd_unbiased(au_to_kcal,bias_fact,kt,kb,md_steps,mtd_steps,w_cv,w_hill &
           & ,ncv,t_min,t_max,gridmin,gridmax,griddif,vbias,ct,nbin,m,mtd,code_name,nr)

CALL mtd_pot(md_steps,mtd_steps,w_cv,w_hill,t_min,t_max,gridmin,gridmax,griddif,vbias, &
           & ct,m,u,ncv,kt,nbin,cv,den,prob_mtd,norm,ir,nr,mtd,code_name)
ENDIF
IF (pmf) THEN
ALLOCATE (fes(nr))
  CALL mean_force(max_step,u,ncv,cv,nr,kt,nbin,t_min,t_max,pcons,kcons &
                 & ,fes,mtd,w_cv,w_hill,ct,vbias,code_name)
  IF(spline) CALL bspline(nr,pcons,fes)
DEALLOCATE(pcons,fes)
ELSE
!---------------------------------------------------------------------------------------------------------------------------!
!Computing 1-Dimensional Probability along Umbrella Coordinate
IF (Prob_nD .eq. 1) THEN
CALL oneD_prob(ii,nr,u,m,mtd,max_step,t_min,t_max,md_steps,den,prob,prob_mtd,ncv,cv,nbin &
               & ,gridmin,griddif,norm,code_name)
DEALLOCATE(prob)
!---------------------------------------------------------------------------------------------------------------------------!
!Computing 2-Dimensional Probability along Umbrella/MTD (if MTD is enabled) or US/TEMP
ELSEIF (Prob_nd .eq. 2) THEN
IF (mtd .eq. 'y') THEN
CALL twoD_prob(ii,jj,u,m,nr,kt,w_cv,w_hill,ct,vbias,max_step,t_min,t_max,md_steps,den,prob_mtd, &
              & ncv,cv,nbin,gridmin,gridmax,griddif,norm,mtd,code_name)

ELSEIF (mtd .eq. 'n') THEN
CALL twoD_temp_prob(jj,u,nr,max_step,t_min,t_max,md_steps,prob_2D,ncv,cv,nbin, &
              & gridmin,griddif,mtd,code_name)
!---------------------------------------------------------------------------------------------------------------------------!
DEALLOCATE(prob_2D) ; DEALLOCATE(prob_mtd)
DEALLOCATE(ct,cv,vbias,nbin,indx)
ENDIF ; ENDIF ; ENDIF
!===========================================================================================================================!
END PROGRAM WSMTD_rw_2D
