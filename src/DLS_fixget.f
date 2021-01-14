!      SUBROUTINE MATRIX_FIX(key,key_RD,keyEL,keySUB
!     &                     ,key_org,key_fx,key_grid1
!     &                     ,KR,R,RD,KN1,grid1,KM,ANGLE
!     &                     ,WAVEL,KRE,KIM,ARE,AIM,pomin,pomax
!     &                     ,NRATN,RATIO
!     &                     ,distname_O,distname_F,distname_N, NDP)
      SUBROUTINE MATRIX_FIX(KN1, grid1
      , WAVEL, KRE, KIM, ARE, AIM
      , RATIO
      , NDP)
!c ** matrix_fixget_S.f previous file name 30/06/2011
!c ** ANGLE SPLINE interpolation version
!c ** rootdir is defined in "mo_par_DLS.f90"
!c ** 12/04/03 f22 interpolation is logarithmic
!c ** 05/05/03 this version can be used to retrieve an aspect ratio
!c ** 13/10/03 IF(ANGLE<=40.)ln(f33)-interpol.
!c **          IF(ANGLE<=50.)ln(f44)-interpol.
!c ** 16/07/18 IF(ANGLE<=20.)ln(f33)-interpol.
!c **          IF(ANGLE<=20.)ln(f44)-interpol.

!c **************************************************************** c
!c **   Subroutine gets original and calculates fixed kernel     ** c
!c **   matrices for given                                       ** c
!c **   aspect ratio distribution and scattering angles          ** c
!c **                                                            ** c
!c **************************************************************** c
!c **                                                            ** c
!c ** INPUT:                                                     ** c
!c **                                                            ** c
!c **   key  = 1 - create fixed kernels (for fixed axis          ** c
!c **              ratio distr.) and save them                   ** c
!c **              into 'Rke...fix...' files and calculate       ** c
!c **              opt.characteristics                           ** c
!c **          2 - read fixed kernels from 'Rke...fix...'  files ** c
!c **          3 - create fixed kernels but don't save them      ** c
!c **          4 - don't create fixed kernels, calculate         ** c
!c **              opt.characteristics from original kernels     ** c
!c **   key_RD =1 - volume mixture of spheroids                  ** c
!c **           2 - surface area  mixture of spheroids           ** c
!c **   keyEL=  1 - calculate F11                                ** c
!c **           2 - F11,F12                                      ** c
!c **           3 - F11,F12,F22                                  ** c
!c **           4 - F11,F12,F22,F33                              ** c
!c **           5 - F11,F12,F22,F33,F34                          ** c
!c **           6 - F11,F12,F22,F33,F34,F44                      ** c
!c **   key_org=0 - read original kernels simultaneously make    ** c
!c **               angle interpolation and calculate opt.char.  ** c
!c **           1 -  -"-, save new kernels in                    ** c
!c **                /NEWORG/ directory, STOP                    ** c
!c **   key_fx works when key=1                                  ** c
!c **        =0 -  create fixed kernels (for fixed axis          ** c
!c **              ratio distr.) and save them                   ** c
!c **              into 'Rke...fix...' files and calculate       ** c
!c **              opt.characteristics                           ** c
!c **         1 - save fixed kernels with original kernel format ** c
!c **             in order to be used as input kernels;          ** c
!c **             'Rke...fix' kernels have to be renamed and moved* c
!c **             into directory 'dir_name'(see 'matrix_fixget.f')* c
!c **             The files can not be used if key=2.            ** c
!c **                                                            ** c
!c **   key_grid1 read grid radii and scat.angles which were used** c
!c **             for kernel look up tables or fixed kernels     ** c
!c **          =0 - 'grid1.txt'                                  ** c
!c **           1 - 'grid1_fix.txt'                              ** c
!c **   KR  - number of aspect ratios                            ** c
!c **   R(KR)  - grid ratios                                     ** c
!c **   RD(KR) - aspect ratio distribution for grid aspect ratios** c
!c **   KM   - number of scattering angles for fixed kernels     ** c
!c **   ANGLE(KM)- scattering angles for fixed kernels           ** c
!c **   distname_O - original kernel directory name              ** c
!c **   distname_F - .fix kernel directory name                  ** c
!c **   distname_N - new original kernel directory               ** c
!c **                                      name (key_org=1)      ** c
!c **                                                            ** c
!c ** OUTPUT:                                                    ** c
!c **                                                            ** c
!c **   UF...(KMpar,KN1par,KIMpar,KREpar) - kernels for          ** c
!c **           given aspect ratio distribution                  ** c
!c **                                                            ** c
!c **   UFEA(2,KN1par,KIMpar,KREpar) - extinction and absorption ** c
!c **                                                    kernels ** c
!c **              1 - extinction                                ** c
!c **              2 - absorption                                ** c
!c **   KN1 - number of grid radii from original or fixed kernels** c
!c **   grid1(KN1) - grid radii                                  ** c
!c **   WAVEL - wavelength from original or fixed kernels        ** c
!c **   KRE   - number of real parts of refr.ind.                ** c
!c **   KIM   - number of imaginary parts of refr.ind.           ** c
!c **   ARE(KRE) - real parts of refr.ind.                       ** c
!c **   AIM(KIM) - imaginary parts of refr.ind.                  ** c
!c **************************************************************** c
!c **************************************************************** c
!c
!c **   key_grid1 read grid radii and scat.angles which were used** c
!c **             for kernel look up tables or fixed kernels     ** c
!c **          =0 - 'grid1.txt'                                  ** c
!c **           1 - 'grid1_fix.txt'                              ** c
!c

      use, intrinsic :: IEEE_ARITHMETIC ! fd
      use alloc
      use alloc1
      USE mo_par_DLS
      use mo_DLS, only: key, key_RD, keyEL, keySUB
      , key_org, key_fx, key_grid1
      , KR, R, RD, KM, ANGLE
      , pomin, pomax
      , distname_O, distname_F, distname_N
      , NRATN, comm_name
!      use mod_os
      use mo_intrpl_linear
      use mo_intrpl_spline

      integer          ::  NEL(0:6)
      character(len=2) ::  NELC(0:6)

      real rmin1, rmax1, rgrid1, rgrid2
      real RATIO(KR1par), RRATN
      dimension grid1(KN1par)
      , ANGLE1(KM1par)
      , ANGLE2(KM1par)
      , RAB1(KM1par)
      , RAB2(KM1par, KN1par, KIMpar, KREpar)
      dimension RDc(KR), X(2), Y(2)
      , ARE(KREpar), AIM(KIMpar)
      real WAVEL
!      real LINEAR                !AH original
      integer NDP
!cl     &, LINEAR_LN
      CHARACTER(255) name, dir_name_O, dir_name_F, dir_name_N
      CHARACTER(255) full_name
      CHARACTER(255) comm_name1

!c      save NDP
      INTEGER                      :: key_spln, NNEL, KKEL
      double precision, dimension(KM1par)    :: XXS2, YYS2
      double precision                       :: XARG, YFIT
      double precision, dimension(KM1par + 4)  :: KS1, CS1
!c ----- Timer ------
      !real*4, external :: etime,dtime

      real ::  tarray(2), T_RFM, T_CFM, T_CFM0
      real   ::  UTEST, tiny = 1e-29
      real :: xx1, xx2
      real :: ang_f33, ang_f44

      ang_f33 = 20.0
      ang_f44 = 20.0

!cl    write(*,*)'distname_O=',TRIM(distname_O)
!cl          write(*,*)'distname_F=',TRIM(distname_F)
!cl          write(*,*)'distname_N=',TRIM(distname_N)

      !dir_name_O=TRIM(rootdir)//TRIM(distname_O)//'/'
      !dir_name_F=TRIM(rootdir)//TRIM(distname_F)//'/'
      !dir_name_N=TRIM(rootdir)//TRIM(distname_N)//'/'

      dir_name_O = TRIM(distname_O)//'/'
      dir_name_F = TRIM(distname_F)//'/'
      dir_name_N = TRIM(distname_N)//'/'

      PI2 = 2.*ACOS(-1.)
      NEL = (/0, 11, 12, 22, 33, 34, 44/)
      NELC = (/'00', '11', '12', '22', '33', '34', '44'/)
      T_RFM = 0.
      T_CFM = 0.
      write (*, *) 'key_grid1=', key_grid1

      if (key_grid1 .eq. 0) then
         full_name = TRIM(dir_name_O)//'grid1.txt'
      else
         full_name = TRIM(dir_name_F)//'grid1_fix.txt'
      end if ! key_grid1
      write (*, *) 'file=', TRIM(full_name)
      OPEN (10, file=full_name, status='unknown')
      READ (10, *) KN1, XWL

      if (KN1 .gt. KN1par) then
         write (*, *) 'in GET_MATRIX: KN1=', KN1, ' KN1par=', KN1par
         stop ' !!! KN1.ne.KN1par'
      end if
      DO I = 1, KN1
         READ (10, *) grid1(I)
      END DO ! I
      rgrid1 = grid1(1)
      rgrid2 = grid1(KN1)
      READ (10, *) KM1
      if (KM1 .gt. KM1par) then
         write (*, *) 'in GET_MATRIX: KM1=', KM1, ' KM1par=', KM1par
         stop ' !!! KM1.ne.KM1par'
      end if
      DO J = 1, KM1
         READ (10, *) ANGLE1(J)
      END DO ! J
      CLOSE (10)
      if ((key .eq. 4) .and. (KM .ne. KM1)) then
         write (*, *) 'if key=4, in input.dat KM=', KM,
         ' should be equal to KM1=', KM1, ' in grid1.txt'
         stop 'STOP in MATRIX_FIX (matrix_fixget.f)'
      end if ! key
!c
!c ** Redefine grids for fixed or NEWORG kernels
!c
      NDPP = 0
      if (key .eq. 1 .or. key_org .eq. 1) then
         rgmin = pomin*XWL/pi2
         rgmax = pomax*XWL/pi2
!tl        write(*,*) 'pomin=',pomin,' XWL=',XWL,' pi2=',pi2,rgmin
!tl        stop
!c
         xx1 = rgmin
         xx2 = grid1(1)
         if (xx1 .lt. xx2) then
            write (*, *) 'XWL=', XWL,&
             &                        ' IS wave length in kernels equal to XWL?'
            write (*, *) 'check input.dat: rgmin=', rgmin, ' < ',
            'grid1(1)=', grid1(1)
            STOP 'STOP in MARTRIX_FIX (matrix_fixget.f)'
         end if
         xx1 = rgmax
         xx2 = grid1(KN1)
         if (xx1 .gt. xx2) then
            write (*, *) 'XWL=', XWL,&
             &                        ' IS wave length in kernels equal to XWL?'
            write (*, *) 'check input.dat: rgmax=', rgmax, ' > ',
            'grid1(KN1)=', grid1(KN1)
            STOP 'STOP in MARTRIX_FIX (matrix_fixget.f)'
         end if

         ind = 0
         do i = 1, KN1 - 1
!c        write(*,*) 'i=',i,' rgmin=',rgmin,' grid1(i+1)=',grid1(i+1)
            if (rgmin .le. grid1(i + 1)) then
!c        xr1=grid1(i)
               nn1 = i
               ind = ind + 1
            end if
            if (ind .eq. 1) EXIT
         end do ! i
!c
         ind = 0
         do i = KN1, 2, -1
!c        write(*,*) 'i=',i,' rgmax=',rgmax,' grid1(i-1)=',grid1(i-1)
            if (rgmax .ge. grid1(i - 1)) then
!c        xr2=grid1(i)
               nn2 = i
               ind = ind + 1
            end if
            if (ind .eq. 1) EXIT
         end do ! i
!c
         if (rgmin .eq. grid1(1)) then
!c        xr1=grid1(1)
            nn1 = 1
         end if ! rgmin
         if (rgmax .eq. grid1(KN1)) then
!c        xr2=grid1(KN1)
            nn2 = KN1
         end if ! rgmax
         nn3 = nn2 - nn1 + 1
!cl        write(*,*) 'nn1=',nn1,' nn2=',nn2,' nn3=',nn3
!cl        write(*,*) 'xr1=',xr1,' xr2=',xr2,' KN1=',KN1

!cl        STOP 'TEST STOP'
!c
!      write(*,*) 'before grid1_fix.txt'
         if (key_org .eq. 0) then
            full_name = TRIM(dir_name_F)//'grid1_fix.txt'
         else
            full_name = TRIM(dir_name_N)//'grid1.txt'
         end if
         write (*, *) 'file=', TRIM(full_name)
         open (77, file=full_name, status='replace')
         WRITE (77, '(i4,f8.3)') nn3, XWL
         DO I = nn1, nn2
            WRITE (77, '(e15.7)') grid1(I)
         END DO ! I
         WRITE (77, *) KM
         DO J = 1, KM
            WRITE (77, '(f6.2)') ANGLE(J)
         END DO ! J
         close (77)
         WRITE (*, *)
         WRITE (*, *) '  ATTENTION:'
         WRITE (*, *) '    New input file grid1.txt (key=1,key_org=1) or'
         WRITE (*, *) '    grid1_fix.txt (key=1,key_org=0) has been created'
         WRITE (*, *) '    for your further calculations !!!'
         WRITE (*, *)
      end if ! key&key_org
!			этот кусок оставляем
			
!c
!c ** Define directories
!c
!c ** for AERONET
!cl      NNEL=2
!cl      if(key_f11.eq.0) NNEL=7
      NNEL = keyEL + 1
      KKEL = 6 - keyEL

      IF (key .EQ. 2) THEN

         open (10, file=TRIM(TRIM(dir_name_F)//"Rkernel1_fix_00.txt"),
         status = "old")
         IF (keyEL .gt. 0) THEN
            open (11, file=TRIM(TRIM(dir_name_F)//"Rkernel1_fix_11.txt"),
            status = "old")
            if (keyEL .gt. 1)
            open (12, file=TRIM(dir_name_F)//"Rkernel1_fix_12.txt",
            status = "old")
            if (keyEL .gt. 2)
            open (13, file=TRIM(dir_name_F)//"Rkernel1_fix_22.txt",
            status = "old")
            if (keyEL .gt. 3)
            open (14, file=TRIM(dir_name_F)//"Rkernel1_fix_33.txt",
            status = "old")
            if (keyEL .gt. 4)
            open (15, file=TRIM(dir_name_F)//"Rkernel1_fix_34.txt",
            status = "old")
            if (keyEL .gt. 5)
            open (16, file=TRIM(dir_name_F)//"Rkernel1_fix_44.txt",
            status = "old")
         END IF ! keyEL
         open (20, file='CHECK.dat', status='unknown')

         DO II = 1, NNEL
            WRITE (20, *) 'INPUT CHECK: Nelement=', NEL(II - 1)
            READ (10 + II - 1, *) key_RD1
            IF (key_RD .ne. key_RD1) then
               WRITE (20, 14) II, key_RD, key_RD1
               WRITE (20, *) 'STOP: key_RD.ne.keyRD1 in Rke...fix'
               WRITE (*, 14) II, key_RD, key_RD1
               WRITE (*, *) 'STOP: key_RD.ne.keyRD1 in Rke...fix'
               STOP
            END IF ! key_RD

            READ (10 + II - 1, *) rmin, rmax
            rmin1 = rmin
            rmax1 = rmax
            IF (rmin1 .ne. rgrid1 .or. rmax1 .ne. rgrid2) THEN
               WRITE (20, *) 'grid1(1)=', rgrid1, ' rmin=', rmin1,
               ' rmax=', rmax1, ' grid1(KN1)=', rgrid2
               WRITE (*, *) 'grid1(1)=', rgrid1, ' rmin=', rmin1,
               ' rmax=', rmax1, ' grid1(KN1)=', rgrid2
               WRITE (*, *) 'Compare1: grid1(1) in grid1.dat',
               ' with RMIN in Rke...fix'
               WRITE (*, *) 'Compare2: grid1(KN1) in grid1.dat',
               ' with RMAX in Rke...fix'
               STOP 'STOP: grid1(1)/grid(KN1).ne.rmin/rmax'
            END IF
            READ (10 + II - 1, *) kk
            WRITE (20, *) 'KR=', KR, ' KR=', kk
            READ (10 + II - 1, *)
            DO I = 1, KR
               READ (10 + II - 1, *) xx, yy
               WRITE (20, *) R(I), RD(I), ' R,RD', xx, yy, ' R,RD'
            END DO ! I
            READ (10 + II - 1, *) kk
            WRITE (20, *) 'KN1=', KN1, ' KN1=', -kk
            IF (KN1 .ne. -kk) THEN
               WRITE (*, *) 'in grid1.dat KN1=', KN1,
               ' .ne. KN1=', -kk, ' in Rke...fix'
               STOP 'STOP in matrix_fixget'
            END IF
         END DO ! II

         DO II = 2, NNEL
            WRITE (20, *) 'INPUT CHECK: Nelement=', NEL(II - 1)
            READ (10 + II - 1, *) kk
            WRITE (20, *) 'KM=', KM, ' KM=', kk
            READ (10 + II - 1, *) ANGLE2(:KM)
            WRITE (20, *) 'SCATTERING ANGLES2:'
            WRITE (20, 15) ANGLE2(:KM)
            WRITE (20, *) 'SCATTERING ANGLES:'
            WRITE (20, 15) ANGLE(:KM)
         END DO ! II

         DO II = 1, NNEL
            WRITE (20, *) 'INPUT CHECK: Nelement=', NEL(II - 1)
            READ (10 + II - 1, *) xx, yy
            WRITE (20, *) 'ARE(1),ARE(KRE)', xx, yy
            READ (10 + II - 1, *) xx, yy
            WRITE (20, *) 'AIM(1),AIM(KIM)', xx, yy
            READ (10 + II - 1, *) KRE, KIM
            WRITE (20, *) 'KRE,KIM', KRE, KIM
            IF (KIM .lt. 0) KIM = -KIM
         END DO ! II
         WRITE (*, *) 'Fixed kernel matrices have been read :'
         DO IRE = 1, KRE
            DO IIM = 1, KIM
               DO II = 1, NNEL
                  if (key_fx .eq. 1) then
                     READ (10 + II - 1, *) kk, aa
                  else ! key_fx=0
                     READ (10 + II - 1, *) kk
                  end if ! key_fx
                  READ (10 + II - 1, *) WAVEL, ARE(IRE), AIM(IIM)
                  IF (AIM(IIM) .LT. 0) AIM(IIM) = -AIM(IIM)
                  !cl        WRITE(20,*) 'WAVEL,ARE,AIM,NEL',
                  !cl     &                WAVEL,ARE(IRE),-AIM(IIM),NEL(II-1)
               END DO ! II

               READ (10, *)
               READ (10, 11) UFEA(1, :KN1, IIM, IRE)
               READ (10, *)
               READ (10, 11) UFEA(2, :KN1, IIM, IRE)

               if (keyEL .gt. 0) then
                  DO I = 1, KN1
                     READ (11, 11) UF11(1:KM, I, IIM, IRE)
                  END DO ! I
               end if
               if (keyEL .gt. 1) then
                  DO I = 1, KN1
                     READ (12, 11) UF12(1:KM, I, IIM, IRE)
                  END DO ! I
               end if
               if (keyEL .gt. 2) then
                  DO I = 1, KN1
                     READ (13, 11) UF22(1:KM, I, IIM, IRE)
                  END DO ! I
               end if

               if (keyEL .gt. 3) then
                  DO I = 1, KN1
                     READ (14, 11) UF33(1:KM, I, IIM, IRE)
                  END DO ! I
               end if
               if (keyEL .gt. 4) then
                  DO I = 1, KN1
                     READ (15, 11) UF34(1:KM, I, IIM, IRE)
                     do j = 1, KM
                        UTEST = UF34(j, I, IIM, IRE)
                        if ((abs(UTEST) .lt. tiny) .and. (abs(UTEST) .gt. 0.0))
                        UF34(j, I, IIM, IRE) = 0.0
                     end do ! j

                  END DO ! I
               end if
               if (keyEL .gt. 5) then
                  DO I = 1, KN1
                     READ (16, 11) UF44(1:KM, I, IIM, IRE)

                  END DO ! I
               end if

            END DO ! IIM
            WRITE (*, 16) WAVEL, ARE(IRE)
         END DO ! IRE
         close (10)
         if (keyEL .gt. 0) close (11)
         if (keyEL .gt. 1) close (12)
         if (keyEL .gt. 2) close (13)
         if (keyEL .gt. 3) close (14)
         if (keyEL .gt. 4) close (15)
         if (keyEL .gt. 5) close (16)
         close (20)

         IF (keyEL .gt. 0) THEN
            UF11(1:KM, 1:KN1, 1:KIM, 1:KRE) =
            LOG(UF11(1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL)
            if (keyEL .gt. 1)
            UF12(1:KM, 1:KN1, 1:KIM, 1:KRE) =
            UF12(1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL
            if (keyEL .gt. 2)
            UF22(1:KM, 1:KN1, 1:KIM, 1:KRE) =
            LOG(UF22(1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL)
            if (keyEL .gt. 3) then
               DO J = 1, KM
                  UF33(J, 1:KN1, 1:KIM, 1:KRE) =
                  UF33(J, 1:KN1, 1:KIM, 1:KRE)*WAVEL
                  IF (ANGLE(J) .LE. ang_f33)
                  UF33(J, 1:KN1, 1:KIM, 1:KRE) =
                  LOG(UF33(J, 1:KN1, 1:KIM, 1:KRE))
               END DO ! J
            end if
            if (keyEL .gt. 4)
            UF34(1:KM, 1:KN1, 1:KIM, 1:KRE) =
            UF34(1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL
            if (keyEL .gt. 5) then
               DO J = 1, KM
                  UF44(J, 1:KN1, 1:KIM, 1:KRE) =
                  UF44(J, 1:KN1, 1:KIM, 1:KRE)*WAVEL
                  IF (ANGLE(J) .LE. ang_f44)
                  UF44(J, 1:KN1, 1:KIM, 1:KRE) =
                  LOG(UF44(J, 1:KN1, 1:KIM, 1:KRE))
               END DO ! J
            end if
         END IF ! keyEL

13       FORMAT(3E12.4, I4)
14       FORMAT('II=', i2, ' key_RD=', i2, ' key_RD1=', i2)
16       FORMAT(12x, 'wl=', f5.2, 2x, 'n=', f8.5)
!cl      WRITE(*,*) 'Fixed kernel matrices have been read'
         IF (key_RD .EQ. 1) WRITE (*, *) 'Volume mixture of spheroids'
         IF (key_RD .EQ. 2) WRITE (*, *) 'Surface area mixture of spheroids'

61       format('  Read fixed matr. ........ ', f8.3, ' min.')
      END IF ! key=2

C 			!========================== этот блок и ниже возможно надо удалить
C       IF (key .NE. 2) THEN
C          IF (key .EQ. 1) THEN
C !c
C !c *** ALLOCATE and INITIALIZE ARRAYS
C !c
C             ALLOCATE (UEA(KR1par, 2, KN1par, KIMpar, KREpar), stat=ierr)
C             if (ierr /= 0) stop 'Can not allocate UEA array'
C             UEA = 0.
C             IF (keyEL .gt. 0) THEN
C                ALLOCATE (U11(KR1par, KMpar, KN1par, KIMpar, KREpar), stat=ierr)
C                if (ierr /= 0) stop 'Can not allocate U11 array'
C                U11 = 0.
C                if (keyEL .gt. 1) then
C                   ALLOCATE (U12(KR1par, KMpar, KN1par, KIMpar, KREpar), stat=ierr)
C                   if (ierr /= 0) stop 'Can not allocate U12 array'
C                   U12 = 0.
C                end if
C                if (keyEL .gt. 2) then
C                   ALLOCATE (U22(KR1par, KMpar, KN1par, KIMpar, KREpar), stat=ierr)
C                   if (ierr /= 0) stop 'Can not allocate U22 array'
C                   U22 = 0.
C                end if
C                if (keyEL .gt. 3) then
C                   ALLOCATE (U33(KR1par, KMpar, KN1par, KIMpar, KREpar), stat=ierr)
C                   if (ierr /= 0) stop 'Can not allocate U33 array'
C                   U33 = 0.
C                end if
C                if (keyEL .gt. 4) then
C                   ALLOCATE (U34(KR1par, KMpar, KN1par, KIMpar, KREpar), stat=ierr)
C                   if (ierr /= 0) stop 'Can not allocate U34 array'
C                   U34 = 0.
C                end if
C                if (keyEL .gt. 5) then
C                   ALLOCATE (U44(KR1par, KMpar, KN1par, KIMpar, KREpar), stat=ierr)
C                   if (ierr /= 0) stop 'Can not allocate U44 array'
C                   U44 = 0.
C                end if
C             END IF ! keyEL=0
C          END IF ! key=1
C          IF (NDP .EQ. 0) THEN
C             T_CFM = 0.
C !c
C !c ** READ ORIGINAL kernels
C !c
C !c      OPEN (10,file=TRIM(rootdir)//'name.dat',status='old')
C             !OPEN (10,file='name.dat',status='old')
C             !READ(10,*) NRATN
C             IF (NRATN .gt. KR1par)
C             STOP ' in GET_MATRIX 1: NRATN.gt.KR1par !!!'
C
C             DO IRATN = 1, NRATN
C !c **
C !c ** READ U11, U12, U22, U33, U34, U44 matrices
C !c **
C                DO KEL = 1, NNEL - 1
C                   !READ(10,*) name
C !cl        OPEN(11,FILE=name,status='old')
C
C                   comm_name1 = trim(comm_name(IRATN))//'_'//NELC(KEL)//'.txt'
C                   full_name = TRIM(dir_name_O)//trim(comm_name1)
C                   write (*, *) trim(full_name)
C                   OPEN (11, FILE=trim(full_name), status='old')
C                   if (key_org .eq. 1) then
C                      full_name = TRIM(dir_name_N)//trim(comm_name1)
C !cl        full_name=TRIM(name)
C                      OPEN (21, FILE=trim(full_name), status='unknown')
C                      write (*, *) trim(full_name)
C                   end if ! key_org
C                   READ (11, *) rmin, rmax, RATIO(iratn)
C                   READ (11, *) KNN
C                   IF (KNN .LT. 0) KNN = -KNN
C                   READ (11, *) KMM
C                   IF (KNN .ne. KN1) STOP ' in GET_MATRIX 1: KNN.ne.KN1 !!!'
C                   IF (KMM .ne. KM1) STOP ' in GET_MATRIX 1: KMM.ne.KM1 !!!'
C                   READ (11, *) ANGLE2(1:KMM)
C !cl        WRITE(*,*)  'ANGLE2 in GET_MATRIX NEL=',NEL(KEL)
C !cl        WRITE(*,11) ANGLE2(:KM1)
C !cl        WRITE(*,*)  'ANGLE1  in GET_MATRIX NEL=',NEL(KEL)
C !cl        WRITE(*,11) ANGLE1(:KM1)
C
C                   WRITE (*, *) 'READ matrix U', NEL(KEL)
C                   READ (11, *) RREMIN, RREMAX
C !cl        WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'
C                   READ (11, *) RIMMIN, RIMMAX
C !cl        WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'
C                   READ (11, *) KRE, KIM
C !cl        WRITE(*,*) KRE, KIM,' KRE, KIM'
C                   if (key_org .eq. 1) then
C                      WRITE (21, *) grid1(nn1), grid1(nn2), RATIO(iratn),
C                      '  rmin, rmax, RATIO'
C                      WRITE (21, *) - nn3, ' number of intervals'
C                      WRITE (21, *) KM, '  number of angles'
C                      WRITE (21, 27) ANGLE(1:KM)
C
C                      WRITE (*, *) 'WRITE matrix U', NEL(KEL)
C                      WRITE (21, *) RREMIN, RREMAX, ' real refr. indices'
C !cl        WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'
C                      WRITE (21, *) RIMMIN, RIMMAX, ' imag refr. indices'
C !cl        WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'
C                      WRITE (21, *) KRE, KIM, ' number of intervals for opt. const'
C !cl        WRITE(*,*) KRE, KIM,' KRE, KIM'
C                   end if ! key_org=1
C
C                   IF (KRE .LT. 0) KRE = -KRE
C                   IF (KIM .LT. 0) KIM = -KIM
C                   DO IRE = 1, KRE
C                      DO IIM = 1, KIM
C                         READ (11, *) nn, xx
C !cl        WRITE(*,*) nn,xx,' KEL,RATIO'
C                         READ (11, *) WAVEL, ARE(IRE), AIM(IIM)
C                         AIM(IIM) = -AIM(IIM)
C
C                         if (key_org .eq. 1) then
C                            if (NDPP .eq. 0) then
C                               if (XWL .ne. WAVEL) then
C                                  write (*, *) 'in grid1.dat: XWL=', XWL, '.NE.',
C                                  ' in Rke...: WAVEL=', WAVEL
C                                  STOP 'STOP in MATRIX_FIX (matrix_fixget.f)'
C                               end if ! XWL
C                               NDPP = 1
C                            end if ! NDPP
C                            WRITE (21, 21) nn, xx
C !cl        WRITE(*,*) nn,xx,' KEL,RATIO'
C                            WRITE (21, 20) WAVEL, ARE(IRE), -AIM(IIM)
C                         end if ! key_org=1
C
C ! key .lt. 4
C                         IF (key .lt. 4) THEN
C                            IF (KEL .EQ. 1) then
C                               DO I = 1, KN1
C                                  READ (11, *) RAB1(1:KM1)
C                                  RAB2(1:KM1, I, IIM, IRE) = RAB1(1:KM1)
C                                  XXS2(1:KM1) = ANGLE1(1:KM1)
C                                  YYS2(1:KM1) = LOG(RAB1(1:KM1))
C                                  key_spln = 0
C                                  DO J = 1, KM
C                                     XARG = ANGLE(J)
C                                     CALL intrpl_spline(KM1, XXS2(1:KM1), YYS2(1:KM1),&
C                                     &XARG, YFIT, key_spln, KS1(1:KM1 + 4), CS1(1:KM1 + 4))
C                                     U11(IRATN, J, I, IIM, IRE) = EXP(YFIT)
C                                     key_spln = 1
C                                  END DO ! J
C                               END DO ! I
C                               if (key_org .eq. 1) then
C                                  DO I = nn1, nn2
C                                     WRITE (21, 27) U11(IRATN, 1:KM, I, IIM, IRE)
C                                  END DO ! I
C                               end if ! key_org=1
C                            ELSE IF (KEL .EQ. 2) THEN
C                               DO I = 1, KN1
C                                  READ (11, *) RAB1(1:KM1)
C                                  if (ANGLE1(1) .eq. 0.) RAB1(1) = 0.
C                                  if (ANGLE1(KM1) .eq. 180.) RAB1(KM1) = 0.
C                                  XXS2(1:KM1) = ANGLE1(1:KM1)
C                                  YYS2(1:KM1) = RAB1(1:KM1)/RAB2(1:KM1, I, IIM, IRE)
C                                  key_spln = 0
C                                  DO J = 1, KM
C                                     XARG = ANGLE(J)
C                                     CALL intrpl_spline(KM1, XXS2(1:KM1),&
C                                           &        YYS2(1:KM1), XARG, YFIT, key_spln,&
C                                                      & KS1(1:KM1 + 4), CS1(1:KM1 + 4))
C                                     if (J .eq. 1 .and. ANGLE1(J) .eq. 0.) YFIT = 0.
C                                     if (J .eq. KM1 .and. ANGLE1(KM1) .eq. 180.) YFIT = 0.
C                                     U12(IRATN, J, I, IIM, IRE) =&
C                                 & YFIT*U11(IRATN, J, I, IIM, IRE)
C                                     key_spln = 1
C                                  END DO ! J
C                               END DO ! I
C                               if (key_org .eq. 1) then
C                                  DO I = nn1, nn2
C                                     WRITE (21, 27) U12(IRATN, 1:KM, I, IIM, IRE)
C                                  END DO ! I
C                               end if ! key_org=1
C                            ELSE IF (KEL .EQ. 3) THEN
C                               DO I = 1, KN1
C                                  READ (11, *) RAB1(1:KM1)
C                                  XXS2(1:KM1) = ANGLE1(1:KM1)
C                                  YYS2(1:KM1) = LOG(RAB1(1:KM1))
C                                  key_spln = 0
C                                  DO J = 1, KM
C                                     XARG = ANGLE(J)
C                                     CALL intrpl_spline(KM1, XXS2(1:KM1),
C                                     YYS2(1:KM1)
C                                     , XARG, YFIT, key_spln, KS1(1:KM1 + 4),
C                                     CS1(1:KM1 + 4))
C                                     U22(IRATN, J, I, IIM, IRE) = EXP(YFIT)
C                                     key_spln = 1
C                                  END DO ! J
C                               END DO ! I
C                               if (key_org .eq. 1) then
C                                  DO I = nn1, nn2
C                                     WRITE (21, 27) U22(IRATN, 1:KM, I, IIM, IRE)
C                                  END DO ! I
C                               end if ! key_org=1
C                            ELSE IF (KEL .EQ. 4) THEN
C                               DO I = 1, KN1
C                                  READ (11, *) RAB1(1:KM1)
C                                  XXS2(1:KM1) = ANGLE1(1:KM1)
C                                  YYS2(1:KM1) = RAB1(1:KM1)/RAB2(1:KM1, I, IIM, IRE)
C                                  key_spln = 0
C                                  DO J = 1, KM
C                                     XARG = ANGLE(J)
C                                     CALL intrpl_spline(KM1, XXS2(1:KM1), YYS2(1:KM1)
C                                     , XARG, YFIT, key_spln, KS1(1:KM1 + 4), CS1(1:KM1 + 4))
C                                     U33(IRATN, J, I, IIM, IRE) = YFIT*U11(IRATN, J, I, IIM, IRE)
C                                     key_spln = 1
C                                  END DO ! J
C                               END DO ! I
C                               if (key_org .eq. 1) then
C                                  DO I = nn1, nn2
C                                     WRITE (21, 27) U33(IRATN, 1:KM, I, IIM, IRE)
C                                  END DO ! I
C                               end if ! key_org=1
C                            ELSE IF (KEL .EQ. 5) THEN
C                               DO I = 1, KN1
C                                  READ (11, *) RAB1(1:KM1)
C                                  do j = 1, KM1
C                                     UTEST = RAB1(j)
C                                     if ((abs(UTEST) .lt. tiny) .and. (abs(UTEST) .gt. 0.0))
C                                     RAB1(j) = 0.0
C                                  end do ! j
C                                  if (ANGLE1(1) .eq. 0.) RAB1(1) = 0.
C                                  if (ANGLE1(KM1) .eq. 180.) RAB1(KM1) = 0.
C                                  XXS2(1:KM1) = ANGLE1(1:KM1)
C                                  YYS2(1:KM1) = RAB1(1:KM1)/RAB2(1:KM1, I, IIM, IRE)
C                                  key_spln = 0
C                                  DO J = 1, KM
C                                     XARG = ANGLE(J)
C                                     CALL intrpl_spline(KM1, XXS2(1:KM1), YYS2(1:KM1)
C                                     , XARG, YFIT, key_spln, KS1(1:KM1 + 4), CS1(1:KM1 + 4))
C                                     if (J .eq. 1 .and. ANGLE1(J) .eq. 0.) YFIT = 0.
C                                     if (J .eq. KM1 .and. ANGLE1(J) .eq. 180.) YFIT = 0.
C                                     if (J .eq. -180) then
C                                        write (*, *) 'U34,YFIT,U11 - ',
C                                        U34(IRATN, J, I, IIM, IRE), YFIT, U11(IRATN, J, I, IIM, IRE)
C                                     end if
C                                     if ((abs(YFIT) .lt. tiny) .and. (abs(YFIT) .gt. 0.0))
C                                     YFIT = 0.0
C                                     U34(IRATN, J, I, IIM, IRE) = YFIT*U11(IRATN, J, I, IIM, IRE)
C                                     if ((abs(U34(IRATN, J, I, IIM, IRE)) .lt. tiny) .and.
C                                     (abs(U34(IRATN, J, I, IIM, IRE)) .gt. 0.0))
C                                     U34(IRATN, J, I, IIM, IRE) = 0.0
C                                     key_spln = 1
C                                  END DO ! J
C                               END DO ! I
C                               if (key_org .eq. 1) then
C                                  DO I = nn1, nn2
C                                     WRITE (21, 27) U34(IRATN, 1:KM, I, IIM, IRE)
C                                  END DO ! I
C                               end if ! key_org=1
C                            ELSE IF (KEL .EQ. 6) THEN
C                               DO I = 1, KN1
C                                  READ (11, *) RAB1(1:KM1)
C                                  XXS2(1:KM1) = ANGLE1(1:KM1)
C                                  YYS2(1:KM1) = RAB1(1:KM1)/RAB2(1:KM1, I, IIM, IRE)
C                                  key_spln = 0
C                                  DO J = 1, KM
C                                     XARG = ANGLE(J)
C                                     CALL intrpl_spline(KM1, XXS2(1:KM1), YYS2(1:KM1)
C                                     , XARG, YFIT, key_spln, KS1(1:KM1 + 4), CS1(1:KM1 + 4))
C                                     U44(IRATN, J, I, IIM, IRE) = YFIT*U11(IRATN, J, I, IIM, IRE)
C                                     key_spln = 1
C                                  END DO ! J
C                               END DO ! I
C                               if (key_org .eq. 1) then
C                                  DO I = nn1, nn2
C                                     WRITE (21, 27) U44(IRATN, 1:KM, I, IIM, IRE)
C                                  END DO ! I
C                               end if ! key_org=1
C                            END IF ! KEL
C ! key .eq. 4
C                         ELSE ! key=4
C                            IF (KEL .EQ. 1) then
C                               DO I = 1, KN1
C                                  READ (11, *) U11(IRATN, 1:KM1, I, IIM, IRE)
C                               END DO ! I
C                            ELSE IF (KEL .EQ. 2) THEN
C                               DO I = 1, KN1
C                                  READ (11, *) U12(IRATN, 1:KM1, I, IIM, IRE)
C                               END DO ! I
C                               if (ANGLE1(1) .eq. 0.) U12(IRATN, 1, 1:KN1, IIM, IRE) = 0.
C                               if (ANGLE1(KM1) .eq. 180.) U12(IRATN, KM1, 1:KN1, IIM, IRE) = 0.
C                            ELSE IF (KEL .EQ. 3) THEN
C                               DO I = 1, KN1
C                                  READ (11, *) U22(IRATN, 1:KM1, I, IIM, IRE)
C                               END DO ! I
C                            ELSE IF (KEL .EQ. 4) THEN
C                               DO I = 1, KN1
C                                  READ (11, *) U33(IRATN, 1:KM1, I, IIM, IRE)
C                               END DO ! I
C                            ELSE IF (KEL .EQ. 5) THEN
C                               DO I = 1, KN1
C                                  READ (11, *) U34(IRATN, 1:KM1, I, IIM, IRE)
C                                  do j = 1, KM1
C                                     UTEST = U34(IRATN, j, I, IIM, IRE)
C                                     if ((abs(UTEST) .lt. tiny) .and. (abs(UTEST) .gt. 0.0))
C                                     U34(IRATN, j, I, IIM, IRE) = 0.0
C                                  end do ! j
C                               END DO ! I
C                               if (ANGLE1(1) .eq. 0.) U34(IRATN, 1, 1:KN1, IIM, IRE) = 0.
C                               if (ANGLE1(KM1) .eq. 180.) U34(IRATN, KM1, 1:KN1, IIM, IRE) = 0.
C                            ELSE IF (KEL .EQ. 6) THEN
C                               DO I = 1, KN1
C                                  READ (11, *) U44(IRATN, 1:KM1, I, IIM, IRE)
C                               END DO ! I
C                            END IF ! KEL
C                         END IF ! key
C                      END DO ! IIM
C                   END DO ! IRE
C                   CLOSE (11)
C                   if (key_org .eq. 1) CLOSE (21)
C                END DO ! KEL
C !c **
C !c ** READ UEA MATRIX
C !c **
C                !if(keyEL.lt.6) then
C                !do KEL=1,KKEL
C                !READ(10,*) name
C                !enddo ! KEL
C                !endif ! key_EL
C                !READ(10,*) name
C
C                comm_name1 = trim(comm_name(IRATN))//'_'//NELC(0)//'.txt'
C                full_name = TRIM(dir_name_O)//trim(comm_name1)
C                write (*, *) trim(full_name)
C                OPEN (11, FILE=trim(full_name), status='old')
C
C                if (key_org .eq. 1) then
C                   full_name = TRIM(dir_name_N)//trim(comm_name1)
C                   OPEN (21, FILE=trim(full_name), status='unknown')
C                   write (*, *) trim(full_name)
C                end if ! key_org
C
C                READ (11, *) rmin, rmax, RATIO(iratn)
C !cl      WRITE(*,*) rmin,rmax,RATIO(iratn),' rmin,rmax,ratio'
C                READ (11, *) KNN
C                IF (KNN .LT. 0) KNN = -KNN
C                IF (KNN .ne. KN1) THEN
C                   WRITE (*, *) KNN, KN1, ' KNN,KN1'
C                   STOP ' in GET_MATRIX 2: KNN.ne.KN1 !!!'
C                END IF
C                WRITE (*, *) 'READ matrix U', NEL(0)
C                READ (11, *) RREMIN, RREMAX
C !cl      WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'
C                READ (11, *) RIMMIN, RIMMAX
C !cl      WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'
C                READ (11, *) KRE, KIM
C                if (key_org .eq. 1) then
C                   WRITE (21, *) grid1(nn1), grid1(nn2), RATIO(iratn),
C                   '  rmin, rmax, RATIO'
C                   WRITE (21, *) - nn3, ' number of intervals'
C                   WRITE (*, *) 'WRITE matrix U', NEL(KEL)
C                   WRITE (21, *) RREMIN, RREMAX, ' real refr. indices'
C !cl        WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'
C                   WRITE (21, *) RIMMIN, RIMMAX, ' imag refr. indices'
C !cl        WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'
C                   WRITE (21, *) KRE, KIM, ' number of intervals for opt. const'
C !cl        WRITE(*,*) KRE, KIM,' KRE, KIM'
C                end if ! key_org=1
C
C                IF (KRE .LT. 0) KRE = -KRE
C                IF (KIM .LT. 0) KIM = -KIM
C                DO IRE = 1, KRE
C                   DO IIM = 1, KIM
C                      READ (11, *) nn, xx
C !cl      WRITE(*,*) nn,xx,' KEL,RATIO'
C                      READ (11, *) WAVEL1, ARE(IRE), AIM(IIM)
C                      IF (keyEL .EQ. 0) THEN
C                         WAVEL = WAVEL1
C                      ELSE
C                         IF (WAVEL .ne. WAVEL1) THEN
C                            WRITE (*, *) WAVEL, WAVEL1, ' WAVEL,WAVEL1'
C                            STOP 'in GET_MATRIX: WAVEL.ne.WAVEL1'
C                         END IF
C                      END IF ! keyEL=0
C                      AIM(IIM) = -AIM(IIM)
C                      if (key_org .eq. 1) then
C                         WRITE (21, 21) nn, xx
C !cl        WRITE(*,*) nn,xx,' KEL,RATIO'
C                         WRITE (21, 20) WAVEL, ARE(IRE), -AIM(IIM)
C                      end if ! key_org=1
C
C                      READ (11, *)
C                      READ (11, *) (UEA(IRATN, 1, I, IIM, IRE), I=1, KN1)
C                      if (key_org .eq. 1) then
C                         WRITE (21, *) ' EXTINCTION (1/km, for d()/dlnr m3/m3*km):'
C                         WRITE (21, 27) (UEA(IRATN, 1, I, IIM, IRE), I=nn1, nn2)
C                      end if ! key_org=1
C                      READ (11, *)
C                      READ (11, *) (UEA(IRATN, 2, I, IIM, IRE), I=1, KN1)
C                      if (key_org .eq. 1) then
C                         WRITE (21, *) ' ABSORPTION (1/km, for dv/dlnr m3/m3*km):'
C                         WRITE (21, 27) (UEA(IRATN, 2, I, IIM, IRE), I=nn1, nn2)
C                      end if ! key_org=1
C                   END DO ! IIM
C                END DO ! IRE
C                CLOSE (11)
C !cl        endif ! key_org=0
C             END DO ! IRATN
C             IF (NRATN .lt. KR .and. (R(1) .lt. RATIO(1) .or. R(KR) .lt. RATIO(NRATN)))
C             STOP ' in GET_MATRIX: NRATN.lt.KR'
C             !check input.dat and name.dat)!!!
C             !CLOSE(10)
C
C 27          FORMAT(7E16.7)
C 20          FORMAT(F14.5, 2E16.7, '  wavel,rreal,rimag')
C 21          FORMAT(I3, F14.5, '  element,ratn')
C
C             if (key_org .eq. 1)
C             STOP 'STOP: key_org=1, new kernels have been saved'
C
C             NDP = 1
C
C          END IF ! NDP
C
C !cl ** TEST for R=1 (sphere) write kernel elements
C
C !cl      write(*,*) 'EXT, ABS : n=',ARE(13),' k=',AIM(8)
C !cl      do i=17,38
C !cl      write(*,'(3e17.5)') grid1(i),UEA(1,1,I,8,13),UEA(1,2,I,8,13)
C !cl        enddo
C
C          IF (key .EQ. 4) then
C             goto 555
C             write (iu_output, '(a,i0)') 'test read kernels, key = ', key
C             write (iu_output, '(es12.4,a,es12.4)')
C             minval(UEA(1, 1, 1:KN1, 1:KIM, 1:KRE)),
C             ' <= UEXT <=', maxval(UEA(1, 1, 1:KN1, 1:KIM, 1:KRE))
C             write (iu_output, '(es12.4,a,es12.4)')
C             minval(UEA(1, 2, 1:KN1, 1:KIM, 1:KRE)),
C             ' <= UABS <=', maxval(UEA(1, 2, 1:KN1, 1:KIM, 1:KRE))
C
C             if (keyEL .gt. 0)
C             write (iu_output, '(es12.4,a,es12.4)')
C             minval(U11(1, 1:KM, 1:KN1, 1:KIM, 1:KRE)),
C             ' <= U11 <=', maxval(U11(1, 1:KM, 1:KN1, 1:KIM, 1:KRE))
C             if (keyEL .gt. 1)
C             write (iu_output, '(es12.4,a,es12.4)')
C             minval(U12(1, 1:KM, 1:KN1, 1:KIM, 1:KRE)),
C             ' <= U12 <=', maxval(U12(1, 1:KM, 1:KN1, 1:KIM, 1:KRE))
C             if (keyEL .gt. 2)
C             write (iu_output, '(es12.4,a,es12.4)')
C             minval(U22(1, 1:KM, 1:KN1, 1:KIM, 1:KRE)),
C             ' <= U22 <=', maxval(U22(1, 1:KM, 1:KN1, 1:KIM, 1:KRE))
C             if (keyEL .gt. 3)
C             write (iu_output, '(es12.4,a,es12.4)')
C             minval(U33(1, 1:KM, 1:KN1, 1:KIM, 1:KRE)),
C             ' <= U33 <=', maxval(U33(1, 1:KM, 1:KN1, 1:KIM, 1:KRE))
C             if (keyEL .gt. 4)
C             write (iu_output, '(es12.4,a,es12.4)')
C             minval(U34(1, 1:KM, 1:KN1, 1:KIM, 1:KRE)),
C             ' <= U34 <=', maxval(U34(1, 1:KM, 1:KN1, 1:KIM, 1:KRE))
C             if (keyEL .gt. 5)
C             write (iu_output, '(es12.4,a,es12.4)')
C             minval(U44(1, 1:KM, 1:KN1, 1:KIM, 1:KRE)),
C             ' <= U44 <=', maxval(U44(1, 1:KM, 1:KN1, 1:KIM, 1:KRE))
C 555         continue
C
C !c ** LOG(U...)
C !c
C             if (keyEL .gt. 0)
C             U11(1:KR, 1:KM, 1:KN1, 1:KIM, 1:KRE) =
C             LOG(U11(1:KR, 1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL)
C             !write(*,*) 'after u11'
C             if (keyEL .gt. 1)
C             U12(1:KR, 1:KM, 1:KN1, 1:KIM, 1:KRE) =
C             U12(1:KR, 1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL
C             !write(*,*) 'i am here after u12'
C             if (keyEL .gt. 2)
C             U22(1:KR, 1:KM, 1:KN1, 1:KIM, 1:KRE) =
C             LOG(U22(1:KR, 1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL)
C             !write(*,*) 'after u22'
C             if (keyEL .gt. 3) then
C             DO J = 1, KM
C                U33(1:KR, J, 1:KN1, 1:KIM, 1:KRE) =
C                U33(1:KR, J, 1:KN1, 1:KIM, 1:KRE)*WAVEL
C                IF (ANGLE(J) .LE. ang_f33) THEN
C !         do l=1,kr
C !         do n=1,kre
C !         do k=1,kim
C !         do i=1,kn1
C !         if(U33(l,J,i,k,n) .le. 0.0)
C !     &     write(*,*) l,n,k,i,J,ANGLE(J),U33(l,J,i,k,n)/WAVEL,
C !     &     '  l,n,k,i,J,ANGLE(J),U33(l,J,i,k,n)'
C !         enddo
C !         enddo
C !         enddo
C !         enddo
C
C                   U33(1:KR, J, 1:KN1, 1:KIM, 1:KRE) =
C                   LOG(U33(1:KR, J, 1:KN1, 1:KIM, 1:KRE))
C                END IF
C             END DO ! J
C             do l = 1, KR
C             do j = 1, KM
C             do i = 1, KN1
C             do k = 1, KIM
C             do n = 1, KRE
C !fd
C                if (.not. ieee_is_normal(U33(l, j, i, k, n))) then
C                   write (*, '(a,5i5,a,e16.6)') 'l,j,i,k,n: ', l, j, i, k, n,
C                   '  U33=', U33(l, j, i, k, n)
C                end if
C             end do
C             end do
C             end do
C             end do
C             end do
C             end if
C             !write(*,*) 'after u33'
C
C             if (keyEL .gt. 4) then
C                U34(1:KR, 1:KM, 1:KN1, 1:KIM, 1:KRE) =
C                U34(1:KR, 1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL
C                do l = 1, KR
C                do j = 1, KM
C                do i = 1, KN1
C                do k = 1, KIM
C                do n = 1, KRE
C !fd
C                   if (.not. ieee_is_normal(U34(l, j, i, k, n))) then
C                      write (*, '(a,5i5,a,e16.6)') 'l,j,i,k,n: ', l, j, i, k, n,
C                      '  U34=', U34(l, j, i, k, n)
C                   end if
C                end do
C                end do
C                end do
C                end do
C                end do
C             end if
C             !write(*,*) 'after u34'
C
C             if (keyEL .gt. 5) then
C             DO J = 1, KM
C                U44(1:KR, J, 1:KN1, 1:KIM, 1:KRE) =
C                U44(1:KR, J, 1:KN1, 1:KIM, 1:KRE)*WAVEL
C                IF (ANGLE(J) .LE. ang_f44) THEN
C                   do l = 1, kr
C                   do n = 1, kre
C                   do k = 1, kim
C                   do i = 1, kn1
C                      if (U44(l, J, i, k, n) .le. 0.0)
C                      write (*, *) l, n, k, i, J, ANGLE(J), U44(l, J, i, k, n)/WAVEL,
C                      '  l,n,k,i,J,ANGLE(J),U44(l,J,i,k,n)'
C                   end do
C                   end do
C                   end do
C                   end do
C                   U44(1:KR, J, 1:KN1, 1:KIM, 1:KRE) =
C                   LOG(U44(1:KR, J, 1:KN1, 1:KIM, 1:KRE))
C !cl       write(*,*) 'U44 before intrpl orgnl:',U44(1,j,1,1,1)
C                END IF ! ANGLE(J) .LE. ang_f44
C             END DO ! J
C             end if
C !      write(*,*) 'after u44'
C             RETURN
C          END IF ! key .EQ. 4
C          if (keySUB .eq. 0) then
C             !T_CFM0=dtime(tarray)!+++
C             !T_CFM0=tarray(1)+tarray(2)
C          end if
C
C          IF (key_RD .eq. 2) then
C ! ** RECALCULATE ASPECT RATIO DISTRIBUTION (RDc()=SAREA/VOLUME)
C ! ** RDc()=RD()/RDc(); sumRD=sum(RDc())
C !
C ! ** OBLATE
C             do IR = 1, KR
C             if (R(IR) .lt. 1.) then
C                E = SQRT(1.-R(IR)*R(IR))
C                xa1 = LOG((1.+E)/(1.-E))/E
C                RDc(IR) = 1.5*(R(IR)**(-2./3.) +
C 0              .5*xa1*R(IR)**(4./3.))
C !c ** PROLATE
C             elseif (R(IR) .gt. 1.) then
C                E = SQRT(1.-1./R(IR)/R(IR))
C                xa2 = ASIN(E)/E
C                RDc(IR) = 1.5*(R(IR)**(-2./3.) +
C                xa2*R(IR)**(1./3.))
C !c ** SPHERE
C             elseif (R(IR) .eq. 1.) then
C                RDc(IR) = 3.
C             end if ! R()
C !c ** WRITE ASPECT RATIO DISTRIBUTION
C !c          write(*,*) 'R=',R(IR),' B=',RDc(IR),
C !c     &  ' 1/B=',1./RDc(IR),' RD=',RD(IR)
C
C             end do ! IR
C             RDc(:KR) = RD(:KR)/RDc(:KR)
C          END IF ! key_RD
C
C          IF (key_RD .eq. 1) RDc(:KR) = RD(:KR)
C          sumRD = sum(RDc(:KR))
C          if (keySUB .eq. 0) then
C             write (*, *)
C             do IR = 1, KR
C                write (*, '(''R='',f8.4,4x,'' RDc='',e13.5)') R(IR), RDc(IR)
C             end do ! IR
C             write (*, '(''sumRD='',e13.5)') sumRD
C          end if
C
C !c ** ALL SCATTERING MATRIX ELEMENTS
C          DO IRE = 1, KRE
C          DO IIM = 1, KIM
C !c
C !c **  SCATTERING MATRIX ELEMENTS
C !c
C             DO I = 1, KN1
C             IF (keyEL .gt. 0) THEN
C             DO J = 1, KM
C
C                DO IR = 1, KR
C                IF (NRATN .EQ. 1) then
C                   RRATN = RATIO(1)
C                   IF (RRATN .NE. R(IR)) THEN
C                      WRITE (*, *) 'R=', R(IR), ' .NE. RATIO=', RATIO(1)
C                      STOP 'in subroutine MATRIX_FIX 1'
C !c        WRITE(*,*) 'R has been changed',
C !c     &                  R(IR),' => ',RATIO(1)
C !c        R(IR)=RRATN
C                   END IF
C                ELSE
C                   RRATN = R(IR)
C                END IF
C                IF (NRATN .NE. 1) THEN
C                DO IRATN = 1, NRATN - 1
C                IF (RRATN .GE. RATIO(IRATN) .AND. RRATN .LE. RATIO(IRATN + 1)) THEN
C                   L0 = IRATN
C                   L1 = IRATN + 1
C                END IF
C                END DO
C                IF (RRATN .LE. RATIO(1)) THEN
C                   IF (RRATN .LT. RATIO(1)) THEN
C                      WRITE (*, *) 'R=', R(IR), ' is out of the range:',
C                      RATIO(1), '< R <', RATIO(NRATN)
C                      STOP 'in subroutine MATRIX_FIX 2'
C !c        WRITE(*,*) 'R has been changed',R(IR),' => ',RATIO(1)
C                   END IF
C                   L0 = 1
C                   L1 = 2
C                   R(IR) = RATIO(1)
C                   RRATN = RATIO(1)
C                END IF
C                IF (RRATN .GE. RATIO(NRATN)) THEN
C                   IF (RRATN .GT. RATIO(NRATN)) THEN
C                      WRITE (*, *) 'R=', R(IR), ' is out of the range:',
C                      RATIO(1), '< R <', RATIO(NRATN)
C                      STOP 'in subroutine MATRIX_FIX 3'
C !c        WRITE(*,*) 'R has been changed',R(IR),
C !c     &                            ' => ',RATIO(NRATN)
C                   END IF
C                   L0 = NRATN - 1
C                   L1 = NRATN
C                   R(IR) = RATIO(NRATN)
C                   RRATN = RATIO(NRATN)
C                END IF
C                ELSE
C                L0 = 1
C                L1 = 1
C                END IF
C
C                RRATN1 = RRATN
C
C                X(1) = RATIO(L0)
C                X(2) = RATIO(L1)
C !c
C !c ** U11 ->    UF11
C !c
C                if (keyEL .gt. 0) then
C                   Y(1) = U11(L0, J, I, IIM, IRE)
C                   Y(2) = U11(L1, J, I, IIM, IRE)
C                   UF11(J, I, IIM, IRE) = UF11(J, I, IIM, IRE) +
C                   LINEAR(X, Y, 2, RRATN1)*RDc(IR)
C                end if
C !      write(*,*) 'after UF11',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
C !c
C !c ** U12 ->    UF12
C !c
C                if (keyEL .gt. 1) then
C                   Y(1) = U12(L0, J, I, IIM, IRE)
C                   Y(2) = U12(L1, J, I, IIM, IRE)
C                   UF12(J, I, IIM, IRE) = UF12(J, I, IIM, IRE) +
C                   LINEAR(X, Y, 2, RRATN1)*RDc(IR)
C                end if
C !      write(*,*) 'after UF12',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
C !c
C !c ** U22 ->    UF22
C !c
C                if (keyEL .gt. 2) then
C                   Y(1) = U22(L0, J, I, IIM, IRE)
C                   Y(2) = U22(L1, J, I, IIM, IRE)
C                   UF22(J, I, IIM, IRE) = UF22(J, I, IIM, IRE) +
C                   LINEAR(X, Y, 2, RRATN1)*RDc(IR)
C                end if
C !      write(*,*) 'after UF22',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
C !c
C !c ** U33 ->    UF33
C !c
C                if (keyEL .gt. 3) then
C                   Y(1) = U33(L0, J, I, IIM, IRE)
C                   Y(2) = U33(L1, J, I, IIM, IRE)
C                   UF33(J, I, IIM, IRE) = UF33(J, I, IIM, IRE) +
C                   LINEAR(X, Y, 2, RRATN1)*RDc(IR)
C                end if
C !      write(*,*) 'after UF33',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
C !c
C !c ** U34 ->    UF34
C !c
C                if (keyEL .gt. 4) then
C                   Y(1) = U34(L0, J, I, IIM, IRE)
C                   Y(2) = U34(L1, J, I, IIM, IRE)
C !      UTEST = LINEAR(X,Y,2,RRATN1)
C !      if(ir .eq. 2 .and. j .eq. 180 .and. i .eq. 5 .and.
C !     &   iim .eq. 1 .and. ire .eq. 1) then
C !        write(*,*) X(1),X(2),Y(1),Y(2),RRATN1,RDc(IR),UTEST,
C !     &  '  X(1),X(2),Y(1),Y(2),RRATN1,RDc(IR),UTEST'
C !      endif
C                   UF34(J, I, IIM, IRE) = UF34(J, I, IIM, IRE) +
C                   LINEAR(X, Y, 2, RRATN1)*RDc(IR)
C                end if
C !      write(*,*) 'after UF34',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
C !c
C !c ** U44 ->    UF44
C !c
C                if (keyEL .gt. 5) then
C                   Y(1) = U44(L0, J, I, IIM, IRE)
C                   Y(2) = U44(L1, J, I, IIM, IRE)
C                   UF44(J, I, IIM, IRE) = UF44(J, I, IIM, IRE) +
C                   LINEAR(X, Y, 2, RRATN1)*RDc(IR)
C
C                end if
C !      write(*,*) 'after UF44',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
C                END DO ! IR RD
C                if (keyEL .gt. 0) UF11(J, I, IIM, IRE) = UF11(J, I, IIM, IRE)/sumRD
C                if (keyEL .gt. 1) UF12(J, I, IIM, IRE) = UF12(J, I, IIM, IRE)/sumRD
C                if (keyEL .gt. 2) UF22(J, I, IIM, IRE) = UF22(J, I, IIM, IRE)/sumRD
C                if (keyEL .gt. 3) UF33(J, I, IIM, IRE) = UF33(J, I, IIM, IRE)/sumRD
C                if (keyEL .gt. 4) UF34(J, I, IIM, IRE) = UF34(J, I, IIM, IRE)/sumRD
C                if (keyEL .gt. 5) UF44(J, I, IIM, IRE) = UF44(J, I, IIM, IRE)/sumRD
C
C             END DO   ! J KM
C             END IF ! keyEL
C
C !c
C !c ** EXTINCTION & ABSORPTION
C !c
C             DO J = 1, 2
C             DO IR = 1, KR
C             IF (NRATN .EQ. 1) then
C                RRATN = RATIO(NRATN)
C             ELSE
C                RRATN = R(IR)
C             END IF
C             IF (NRATN .NE. 1) THEN
C             DO IRATN = 1, NRATN - 1
C             IF (RRATN .GE. RATIO(IRATN) .AND. RRATN .LE. RATIO(IRATN + 1)) THEN
C                L0 = IRATN
C                L1 = IRATN + 1
C             END IF
C             END DO
C
C             IF (RRATN .GE. RATIO(NRATN)) THEN
C                L0 = NRATN - 1
C                L1 = NRATN
C                R(IR) = RATIO(NRATN)
C                RRATN = RATIO(NRATN)
C             END IF
C             IF (RRATN .LE. RATIO(1)) THEN
C                L0 = 1
C                L1 = 2
C                R(IR) = RATIO(1)
C                RRATN = RATIO(1)
C             END IF
C             ELSE
C             L0 = 1
C             L1 = 1
C             END IF
C
C             RRATN1 = RRATN
C
C             X(1) = RATIO(L0)
C             X(2) = RATIO(L1)
C
C             Y(1) = UEA(L0, J, I, IIM, IRE)
C             Y(2) = UEA(L1, J, I, IIM, IRE)
C             UFEA(J, I, IIM, IRE) = UFEA(J, I, IIM, IRE) +
C             LINEAR(X, Y, 2, RRATN1)*RDc(IR)
C
C !      write(*,*) 'after UFEA',IR,J,I,IIM,IRE,'  IR,J,I,IIM,IRE'
C             END DO ! IR KR
C
C             UFEA(J, I, IIM, IRE) = UFEA(J, I, IIM, IRE)/sumRD
C
C             END DO ! J 2
C
C             END DO ! I KN1
C
C          END DO ! IIM
C          END DO ! IRE
C
C          !goto 554
C          write (iu_output, '(a)') 'in kernels_tb_shd_equal'
C          write (iu_output, '(es12.4,a,es12.4)')
C          minval(UFEA(1, 1:KN1, 1:KIM, 1:KRE)),
C          ' <= UE  <=', maxval(UFEA(1, 1:KN1, 1:KIM, 1:KRE))
C          write (iu_output, '(es12.4,a,es12.4)')
C          minval(UFEA(2, 1:KN1, 1:KIM, 1:KRE)),
C          ' <= UA  <=', maxval(UFEA(2, 1:KN1, 1:KIM, 1:KRE))
C
C 554      continue
C
C          !T_CFM0=dtime(tarray)
C          !T_CFM0=tarray(1)+tarray(2)
C          !T_CFM=T_CFM+T_CFM0 !+++
C          IF (key .EQ. 1) THEN
C !c
C !c *** SAVE FIXED matrices
C !c
C
C             open (10, file=TRIM(dir_name_F)//'Rkernel1_fix_00.txt',
C             status = 'unknown')
C             IF (keyEL .gt. 0) THEN
C                open (11, file=TRIM(dir_name_F)//'Rkernel1_fix_11.txt',
C                status = 'unknown')
C                if (keyEL .gt. 1)
C                open (12, file=TRIM(dir_name_F)//'Rkernel1_fix_12.txt',
C                status = 'unknown')
C                if (keyEL .gt. 2)
C                open (13, file=TRIM(dir_name_F)//'Rkernel1_fix_22.txt',
C                status = 'unknown')
C                if (keyEL .gt. 3)
C                open (14, file=TRIM(dir_name_F)//'Rkernel1_fix_33.txt',
C                status = 'unknown')
C                if (keyEL .gt. 4)
C                open (15, file=TRIM(dir_name_F)//'Rkernel1_fix_34.txt',
C                status = 'unknown')
C                if (keyEL .gt. 5)
C                open (16, file=TRIM(dir_name_F)//'Rkernel1_fix_44.txt',
C                status = 'unknown')
C             END IF
C             DO II = 1, NNEL
C             if (key_fx .eq. 1) then
C                WRITE (10 + II - 1, 10) grid1(nn1), grid1(nn2), R(KR)
C             else ! key_fx=0
C                WRITE (10 + II - 1, *) key_RD,
C                '  key_RD=1-volume mixture, 2-surface area mixture'
C                WRITE (10 + II - 1, 10) grid1(nn1), grid1(nn2)
C                WRITE (10 + II - 1, *) KR, ' a number of grid aspect ratios'
C                WRITE (10 + II - 1, *) 'aspect ratio distribution'
C                DO I = 1, KR
C                   WRITE (10 + II - 1, 11) R(I), RDc(I)
C                END DO ! I
C             end if ! key_fx
C             WRITE (10 + II - 1, *) - nn3, ' a number of grid radii'
C             END DO ! II
C             DO II = 2, NNEL
C                WRITE (10 + II - 1, *) KM, ' a number of scattering angles'
C                WRITE (10 + II - 1, 15) ANGLE(1:KM)
C             END DO ! II
C             DO II = 1, NNEL
C                WRITE (10 + II - 1, '(2e13.5,'' real  refr. index'')')
C                ARE(1), ARE(KRE)
C                WRITE (10 + II - 1, '(2e13.5,'' imag refr. index'')')
C                -AIM(1), -AIM(KIM)
C                WRITE (10 + II - 1, *) KRE, -KIM, ' number of intervals for opt. const'
C             END DO ! II
C
C             DO IRE = 1, KRE
C             DO IIM = 1, KIM
C             DO II = 1, NNEL
C             if (key_fx .eq. 1) then
C                WRITE (10 + II - 1, '(i3,E13.5,'' element number,ratn'')')
C                NEL(II - 1), R(KR)
C             else ! key_fx=0
C                WRITE (10 + II - 1, '(i3,E13.5,'' element number'')')
C                NEL(II - 1)
C             end if ! key_fx
C             WRITE (10 + II - 1, '(E11.4,2E15.7,'' wavel,rreal,rimag'')')
C             WAVEL, ARE(IRE), -AIM(IIM)
C             END DO ! II
C
C             IF (keyEL .gt. 0) THEN
C                DO I = nn1, nn2
C                   WRITE (11, 11) UF11(1:KM, I, IIM, IRE)
C                END DO ! I
C                if (keyEL .gt. 1) then
C                DO I = nn1, nn2
C                   WRITE (12, 11) UF12(1:KM, I, IIM, IRE)
C                END DO ! I
C                end if
C                if (keyEL .gt. 2) then
C                DO I = nn1, nn2
C                   WRITE (13, 11) UF22(1:KM, I, IIM, IRE)
C                END DO
C                end if
C                if (keyEL .gt. 3) then
C                DO I = nn1, nn2
C                   WRITE (14, 11) UF33(1:KM, I, IIM, IRE)
C                END DO ! I
C                end if
C                if (keyEL .gt. 4) then
C                   DO I = nn1, nn2
C                      WRITE (15, 11) UF34(1:KM, I, IIM, IRE)
C                   END DO ! I
C                end if
C                if (keyEL .gt. 5) then
C                   DO I = nn1, nn2
C                      WRITE (16, 11) UF44(1:KM, I, IIM, IRE)
C                   END DO ! I
C                end if
C             END IF
C
C             WRITE (10, *) 'EXTINCTION (1/km, for d()/dlnr m3/m3*km)'
C             WRITE (10, 11) UFEA(1, nn1:nn2, IIM, IRE)
C             WRITE (10, *) 'ABSORPTION (1/km, for d()/dlnr m3/m3*km)'
C             WRITE (10, 11) UFEA(2, nn1:nn2, IIM, IRE)
C             END DO ! IRE
C             END DO ! IIM
C             close (10)
C             IF (keyEL .gt. 0) THEN
C                close (11)
C                if (keyEL .gt. 1) close (12)
C                if (keyEL .gt. 2) close (13)
C                if (keyEL .gt. 3) close (14)
C                if (keyEL .gt. 4) close (15)
C                if (keyEL .gt. 5) close (16)
C             END IF ! keyEL
C
C             WRITE (*, *) 'Fixed kernels have been calculated and saved'
C             IF (key_RD .EQ. 1) WRITE (*, *) 'Volume mixture of spheroids'
C             IF (key_RD .EQ. 2) WRITE (*, *) 'Surface area mixture of spheroids'
C !c
C !c *** DEALLOCATE ARRAYS
C !c
C             IF (keyEL .gt. 0) THEN
C                DEALLOCATE (U11, stat=ierr)
C                if (ierr /= 0) stop 'Can not deallocate U11 array'
C                DEALLOCATE (UEA, stat=ierr)
C                if (ierr /= 0) stop 'Can not deallocate UEA array'
C                if (keyEL .gt. 1) then
C                   DEALLOCATE (U12, stat=ierr)
C                   if (ierr /= 0) stop 'Can not deallocate U12 array'
C                end if
C                if (keyEL .gt. 2) then
C                   DEALLOCATE (U22, stat=ierr)
C                   if (ierr /= 0) stop 'Can not deallocate U22 array'
C                end if
C                if (keyEL .gt. 3) then
C                   DEALLOCATE (U33, stat=ierr)
C                   if (ierr /= 0) stop 'Can not deallocate U33 array'
C                end if
C                if (keyEL .gt. 4) then
C                   DEALLOCATE (U34, stat=ierr)
C                   if (ierr /= 0) stop 'Can not deallocate U34 array'
C                end if
C                if (keyEL .gt. 5) then
C                   DEALLOCATE (U44, stat=ierr)
C                   if (ierr /= 0) stop 'Can not deallocate U44 array'
C                end if
C             END IF ! keyEL
C          END IF ! key=1
C
C          IF (keyEL .gt. 0) THEN
C             UF11(1:KM, 1:KN1, 1:KIM, 1:KRE) =
C             LOG(UF11(1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL)
C             if (keyEL .gt. 1)
C             UF12(1:KM, 1:KN1, 1:KIM, 1:KRE) =
C             UF12(1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL
C             if (keyEL .gt. 2)
C             UF22(1:KM, 1:KN1, 1:KIM, 1:KRE) =
C             LOG(UF22(1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL)
C             if (keyEL .gt. 3) then
C             DO J = 1, KM
C                UF33(J, 1:KN1, 1:KIM, 1:KRE) =
C                UF33(J, 1:KN1, 1:KIM, 1:KRE)*WAVEL
C                IF (ANGLE(J) .LE. ang_f33)
C                UF33(J, 1:KN1, 1:KIM, 1:KRE) =
C                LOG(UF33(J, 1:KN1, 1:KIM, 1:KRE))
C             END DO ! J
C             end if
C             if (keyEL .gt. 4)
C             UF34(1:KM, 1:KN1, 1:KIM, 1:KRE) =
C             UF34(1:KM, 1:KN1, 1:KIM, 1:KRE)*WAVEL
C             if (keyEL .gt. 5) then
C             DO J = 1, KM
C                UF44(J, 1:KN1, 1:KIM, 1:KRE) =
C                UF44(J, 1:KN1, 1:KIM, 1:KRE)*WAVEL
C                IF (ANGLE(J) .LE. ang_f44)
C                UF44(J, 1:KN1, 1:KIM, 1:KRE) =
C                LOG(UF44(J, 1:KN1, 1:KIM, 1:KRE))
C             END DO ! J
C             end if
C          END IF ! keyEL
C
C          !if(keySUB .eq. 0) then
C          !write(6,*)
C          !write(6,*)
C          !&        '------------------ T I M I N G ------------------'
C          !WRITE(*,62) T_CFM/60.
C          !endif
C 62       format('  Calcul. fixed kernels ... ', f8.3, ' min.')
C
C       END IF ! key.NE.2
C 			!!! Удялять по эту строку


10    FORMAT(3E15.7, ' rmin, rmax, RATIO')
11    FORMAT(7E15.7)
12    FORMAT(3E15.7, I4, '    wavelength, n, k, NEL')
15    FORMAT(7F12.2)

      RETURN
   END SUBROUTINE MATRIX_FIX

!c ****************************************************************
