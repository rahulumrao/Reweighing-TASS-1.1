!--------------------------------------------------------------------------------
! THIS SUBROUTINE COMPUTE'S MEAN FORCE AT EACH UMBRELLA WINDOW AND THEN 
! INTEGRATE THAT MEAN FORCE TO CALCULATE FREE ENERGY (1D) ALONG UMBRELLA CV
! u --> UMBRELLA CV INDEX ; nr --> NUMBER OF UMBRELLA WINDOWS
! pcons --> POSITION OF UMBRELL ; kcons --> Kappa VALUE IN EACH UMBRELLA WINDOW
! dfds = delF/delS ; fes --> FREE ENERGY ALONG UMBRELLA CV
! WRITTEN BY : Rahul Verma (vrahul@iitk.ac.in)
!--------------------------------------------------------------------------------
MODULE MeanForce
USE Input_file        
USE GetSteps
USE GetFileName
USE MTD_Unbais
CONTAINS
SUBROUTINE mean_force(max_step,u,ncv,cv,nr,kt,nbin,t_min,t_max,pcons,kcons &
                 & ,fes,mtd,w_cv,w_hill,ct,vbias,code_name)
IMPLICIT NONE
INTEGER                 :: i,i_md,i_mtd,t_min,t_max,nr
INTEGER                 :: ncv,w_cv,w_hill,md_steps
INTEGER                 :: nbin(*),u,vt_max
REAL*8                  :: diff_s,den,num,kt,dum
!REAL*8                  :: gridmin(*),gridmax(*),griddif(*)
REAL*8                  :: vbias(nr,*),ct(nr,*),cv(nr,ncv,*)
REAL*8                  :: pcons(nr),kcons(nr)
REAL*8,ALLOCATABLE      :: dummy(:,:,:),prob(:)
REAL*8,ALLOCATABLE      :: dfds(:,:),av_dfds(:),norm(:),fes(:)
REAL*8, PARAMETER       :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER       :: au_to_kcal = 627.51
REAL*8, PARAMETER       :: kj_to_kcal = 0.239006

LOGICAL                 :: max_step
!LOGICAL                 :: read_ct,read_vbias
CHARACTER*5             :: mtd
!CHARACTER(LEN=50)       :: filename_loc
CHARACTER(LEN=10)       :: code_name
CHARACTER(LEN=50)       :: filename(nr),filename_mtd(2,nr)

IF(nbin(1) .eq. 0 ) STOP "ERROR : NUMBER OF BINS CAN NOT BE ZERO"

CALL file_input(nr,code_name,mtd,pcons,kcons,filename,filename_mtd)

101 FORMAT (A8,1X,A14,4X,A20,1X,A20,20X,A10,1X,A8,1X,A8)
102 FORMAT (I4,5X,F8.2,11X,F10.2,10X,A30,10X,I8,2X,I8,2X,I8)
write(*,101)'#replica','umbrella_mean','umbrella_k(kcal/mol)','CV_VAL_file',"MD Steps","t_min","t_max"
!------------------------------------------------------------------------------------------------------!
ALLOCATE(dfds(nr,9999999))
ALLOCATE(av_dfds(nr))
!------------------------------------------------------------------------------------------------------!

vt_max = t_max
DO i = 1,nr
OPEN(11,FILE=filename(i),status='old')
CALL get_steps(11,md_steps)
CALL max_t (max_step,t_min,t_max,vt_max,md_steps)   ! get t_max

WRITE(*,102)i, pcons(i), kcons(i), filename(i),md_steps,t_min,t_max
!------------------------------------------------------------------------------------------------------!
ALLOCATE(dummy(nr,ncv,md_steps))

DO i_md=1,md_steps
  diff_s = cv(i,u,i_md) - pcons(i)
  dfds(i,i_md) = -diff_s*kcons(i)                 ! calculating df/ds
ENDDO
!
IF (mtd .eq. 'y') THEN
den = 0.d0 ; num = 0.d0 ; dum = 0.d0

DO i_md = t_min,t_max
  i_mtd = ((i_md-1)*w_cv/w_hill) + 1
    dum = vbias(i,i_md) - ct(i,i_mtd)
    num = num + dfds(i,i_md)*DEXP(dum/kt)
    den = den + DEXP(dum/kt)
ENDDO
ELSEIF (mtd .eq. 'n') THEN
den = 0.d0 ; num = 0.d0 ; dum = 0.d0
   DO i_md = t_min,t_max
     num = num + dfds(i,i_md)
     den = den + DEXP(dum/kt)
   ENDDO
ENDIF
   av_dfds(i) = num/den                            ! average of df/ds along every umbrella
DEALLOCATE(dummy)
ENDDO

PRINT*,"av_dfds is computed"

OPEN(12,FILE='av_dfds.dat')
DO i = 1,nr
WRITE(12,*)pcons(i),av_dfds(i)
ENDDO
!------------------------------------------------------------------------------------------------------!
ALLOCATE(prob(nbin(1)))
ALLOCATE(norm(nr))

OPEN(13,FILE="free_energy.dat")
!
fes = 0.d0 ; num = 0.d0
!Integration using TrapeZoidal Rule 
DO i = 1,nr-1
   dum = pcons(i+1) - pcons(i)
   num = num + dum*(av_dfds(i+1) + av_dfds(i))
   fes(i+1) = num*0.5d0
   WRITE(13,*)pcons(i+1),fes(i+1)
ENDDO
PRINT*,"Free Energy is Computed"
!
DEALLOCATE(av_dfds,dfds)
!------------------------------------------------------------------------------------------------------!
CLOSE(11) ; CLOSE(12) ; CLOSE(13)
END SUBROUTINE 
END MODULE MeanForce

