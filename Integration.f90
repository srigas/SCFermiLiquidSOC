program SC_SOC_INT
    implicit none

    integer :: I, J, TOTPHISTEPS, IO, TOTRSTEPS, NUMKPOLES, NINTRAPOLE, NFINAL, TOTPTS, POLECHECKER, &
    &PTCHECKER, INNERVAL, NEND
    real*8 :: HBAR, M, PI, KFERMI, DELTA, ARASHBA, BETA, MU, ENERGY, ETA, KEND, Q, KSTEP, QMID, R, PHI, PHISTEP, &
    &KVAL, POLESCALE, DENOM, KFINAL, XFACTOR, DFACT, DAMP
    complex*16 :: CI, GREENK(4,4), INTEGRAL(4,4), OUTERINTPREV(4,4), OUTERINTMID(4,4), OUTERINTNEW(4,4), &
    &INNERINTPREV(4,4), INNERINTMID(4,4), INNERINTNEW(4,4), EXPA, EXPB, EXACTCALCMAT(4,4)
    integer, allocatable, dimension(:) :: NUMPTS
    real*8, allocatable, dimension(:) :: KPOLES, MAINKVALS, MIDKVALS
    complex*16, allocatable, dimension(:,:,:) :: MAINGREENS, MIDGREENS

    ! Constants
    HBAR = 1.D0
    M = 0.5D0
    CI = (0.D0,1.D0)
    PI = 4.D0*atan(1.D0)

    ! Variables read from config.txt
    open(1, file = 'config.txt', action = 'read')
    read(1,*) KFERMI
    read(1,*) DELTA
    read(1,*) ARASHBA
    read(1,*) BETA
    read(1,*) ENERGY
    close(1)
    MU = (HBAR*KFERMI)**2/(2.D0*M)

    ! z = Ε + iη
    ETA = 0.001D0

    ! To introduce damping, set DFACT higher than zero
    DFACT = 0.1D0

    ! Counts the number of k-poles produced by the Julia code
    ! Note that the k-poles are already arranged from smallest to largest
    NUMKPOLES = 0
    open(2, file = 'kpoles.dat', action = 'read')
    do
        read(2,*,iostat=io)
        if (io/=0) exit
        NUMKPOLES = NUMKPOLES + 1
    end do
    close(2)

    allocate(KPOLES(NUMKPOLES))
    open(2, file = 'kpoles.dat', action = 'read')
    do i = 1,NUMKPOLES
        read(2,*) KPOLES(i)
    end do
    close(2)


    ! The integration is split into two main parts.
    ! The first has to do with the k-values near the poles, while the second deals with the k-values away
    ! from the final pole. For the first part, the usual scheme is implemented, using a double Simpson
    ! method, as will be mentioned further down. For the second part, the integrand will be approximated
    ! as a second-order polynomial, and will therefore be exactly calculated (within this approximation).

    ! First we calculate all the Green functions to be considered.
    ! -----------------------------------------------------------------------------------

    ! This part calculates the number of points to be considered at the first part.
    allocate(NUMPTS(NUMKPOLES))
    ! NUMPTS(1) is the number of steps to be taken from 0 till the first pole.
    NUMPTS(1) = 50000
    ! NINTRAPOLE is the number of steps to be taken in the vicinity (k_p - η, k_p + η).
    NINTRAPOLE = 80000

    POLECHECKER = 2
    DENOM = KPOLES(1)-4.D0*ETA
    DO WHILE (POLECHECKER < NUMKPOLES+1)
        POLESCALE = (KPOLES(POLECHECKER)-KPOLES(POLECHECKER-1)-8.D0*ETA)/DENOM
        NUMPTS(POLECHECKER) = INT(POLESCALE*NUMPTS(1))
        POLECHECKER = POLECHECKER + 1
    END DO

    ! First part finishes at the point 2*k_p(final), i.e. 2 times the k-value of the final pole.
    POLESCALE = (KPOLES(NUMKPOLES)-4.D0*ETA)/DENOM
    NFINAL = INT(POLESCALE*NUMPTS(1))

    ! Second part starts at the point 2*k_p(final) and ends at KFINAL
    KFINAL = 10000
    NEND = 100000

    TOTPTS = SUM(NUMPTS) + NINTRAPOLE*NUMKPOLES + NFINAL + NEND

    allocate(MAINKVALS(TOTPTS+1))
    allocate(MIDKVALS(TOTPTS))
    allocate(MAINGREENS(4,4,TOTPTS+1))
    allocate(MIDGREENS(4,4,TOTPTS))

    Q = 0.D0
    PTCHECKER = 1
    POLECHECKER = 1
    DO WHILE (POLECHECKER < NUMKPOLES+1)

        IF (POLECHECKER == 1) THEN
            KSTEP = (KPOLES(POLECHECKER)-4.D0*ETA)/NUMPTS(1)
        ELSE
            KSTEP = (KPOLES(POLECHECKER)-KPOLES(POLECHECKER-1)-8.D0*ETA)/NUMPTS(POLECHECKER)
        ENDIF

        DO I = 1, NUMPTS(POLECHECKER)
            MAINKVALS(PTCHECKER) = Q
            CALL GREENKBUILD(HBAR,M,MU,DELTA,ARASHBA,BETA,ENERGY,ETA,Q,GREENK)
            MAINGREENS(:,:,PTCHECKER) = GREENK
            QMID = Q + 0.5D0*KSTEP
            CALL GREENKBUILD(HBAR,M,MU,DELTA,ARASHBA,BETA,ENERGY,ETA,QMID,GREENK)
            MIDGREENS(:,:,PTCHECKER) = GREENK
            MIDKVALS(PTCHECKER) = QMID
            Q = Q + KSTEP
            PTCHECKER = PTCHECKER + 1
        END DO

        KSTEP = 8.D0*ETA/NINTRAPOLE

        DO I = 1, NINTRAPOLE
            MAINKVALS(PTCHECKER) = Q
            CALL GREENKBUILD(HBAR,M,MU,DELTA,ARASHBA,BETA,ENERGY,ETA,Q,GREENK)
            MAINGREENS(:,:,PTCHECKER) = GREENK
            QMID = Q + 0.5D0*KSTEP
            CALL GREENKBUILD(HBAR,M,MU,DELTA,ARASHBA,BETA,ENERGY,ETA,QMID,GREENK)
            MIDGREENS(:,:,PTCHECKER) = GREENK
            MIDKVALS(PTCHECKER) = QMID
            Q = Q + KSTEP
            PTCHECKER = PTCHECKER + 1
        END DO

        POLECHECKER = POLECHECKER + 1

    END DO

    ! Part from k_p(final) + 4η to 2k_p(final)
    KSTEP = (KPOLES(NUMKPOLES)-4.D0*ETA)/NFINAL

    DO I = 1, NFINAL
        MAINKVALS(PTCHECKER) = Q
        CALL GREENKBUILD(HBAR,M,MU,DELTA,ARASHBA,BETA,ENERGY,ETA,Q,GREENK)
        MAINGREENS(:,:,PTCHECKER) = GREENK
        QMID = Q + 0.5D0*KSTEP
        CALL GREENKBUILD(HBAR,M,MU,DELTA,ARASHBA,BETA,ENERGY,ETA,QMID,GREENK)
        MIDGREENS(:,:,PTCHECKER) = GREENK
        MIDKVALS(PTCHECKER) = QMID
        Q = Q + KSTEP
        PTCHECKER = PTCHECKER + 1
    END DO

    ! Part from 2k_p(final) to KFINAL
    KSTEP = (KFINAL-2.D0*KPOLES(NUMKPOLES))/NEND

    DO I = 1, NEND
        MAINKVALS(PTCHECKER) = Q
        CALL GREENKBUILD(HBAR,M,MU,DELTA,ARASHBA,BETA,ENERGY,ETA,Q,GREENK)
        MAINGREENS(:,:,PTCHECKER) = GREENK
        QMID = Q + 0.5D0*KSTEP
        CALL GREENKBUILD(HBAR,M,MU,DELTA,ARASHBA,BETA,ENERGY,ETA,QMID,GREENK)
        MIDGREENS(:,:,PTCHECKER) = GREENK
        MIDKVALS(PTCHECKER) = QMID
        Q = Q + KSTEP
        PTCHECKER = PTCHECKER + 1
    END DO

    CALL GREENKBUILD(HBAR,M,MU,DELTA,ARASHBA,BETA,ENERGY,ETA,Q,GREENK)
    MAINGREENS(:,:,PTCHECKER) = GREENK
    MAINKVALS(PTCHECKER) = Q

    ! Now, the arrays MAINGREENS, MIDGREENS, MAINKVALS and MIDKVALS
    ! contain all the required info for the k-integrations to be performed.
    ! -----------------------------------------------------------------------------------    

    TOTRSTEPS = 100 ! Values for r
    TOTPHISTEPS = 101 ! Steps for the polar angle
    PHISTEP = 2.D0*PI/TOTPHISTEPS ! Step per φ iteration
    R = 0.D0 ! Initial value for r

    open (1, file = 'Green.dat', action = 'write')

    ! This begins the r-loop
    DO I = 1, TOTRSTEPS

        R = R + 1.D0/(20.D0*KFERMI) ! Step for R is equal to 0.1/k_F
        INTEGRAL(:,:) = (0.D0,0.D0)

        ! In general, a double-integral must be calculated.
        ! In any case, the outer sum corresponds to the φ-integration
        ! while the inner sum corresponds to the k-integration. For each φ-iteration,
        ! three k-iterations must be performed: one evaluating Int(φ), one evaluating
        ! Int(φ+dφ) and one evaluating Int(φ+0.5dφ), so that the Simpson method can be applied
        ! for the φ-integration.

        ! For the first part, the same method is applied within the inner sum, where 
        ! Int(φ,k), Int(φ,k+dk) and Int(φ,k+0.5dk) must be calculated, in order to apply the
        ! Simpson method for the k-integration.

        ! For the second part, the k-integrand is approximated by its Taylor second order
        ! expansion, so the integral of the approximate integrand can be analytically
        ! evaluated.

        ! Calculation of the first Int(φ)
        OUTERINTPREV = (0.D0,0.D0)

        PHI = 0.D0 ! First value of φ

        Q = MAINKVALS(1) ! First value of k
        INNERINTPREV = MAINGREENS(:,:,1) ! G(k) for first k
        CALL GTRANSWEXP(INNERINTPREV,PHI,PI,CI,Q,R,DFACT) ! I(φ,k) = k*U(φ)G(k)U+(φ)*exp(ikxcosφ)
        ! This begins the k loop for φ
        ! First part, until 2k_p(final)
        DO INNERVAL = 1, TOTPTS-NEND

            QMID = MIDKVALS(INNERVAL) ! k + 0.5dk
            INNERINTMID = MIDGREENS(:,:,INNERVAL) ! G(k+0.5dk)
            CALL GTRANSWEXP(INNERINTMID,PHI,PI,CI,QMID,R,DFACT) ! I(φ,k+0.5dk)
            Q = MAINKVALS(INNERVAL+1) ! k + dk
            INNERINTNEW = MAINGREENS(:,:,INNERVAL+1) ! G(k+dk)
            CALL GTRANSWEXP(INNERINTNEW,PHI,PI,CI,Q,R,DFACT) ! I(φ,k+0.5dk)

            KSTEP = Q-MAINKVALS(INNERVAL)
            OUTERINTPREV = OUTERINTPREV + (KSTEP/6.D0)*(INNERINTPREV + 4.D0*INNERINTMID + INNERINTNEW)
            INNERINTPREV = INNERINTNEW

        END DO

        ! At this point the k value is equal to MAINKVALS(TOTPTS-NEND+1). From now on, the second part begins.
        ! INNERINTPREV = f(a), INNERINTMID = f(k_0), INNERINTNEW = f(b), where the first a is 2*k_p(final) and
        ! the last b is KFINAL.
        INNERINTPREV = MAINGREENS(:,:,TOTPTS-NEND+1)
        CALL GTRANSWOEXP(INNERINTPREV,PHI,PI,CI,Q) ! I(φ,k) = k*U(φ)G(k)U+(φ), without the exponent
        
        XFACTOR = R*COS(PHI)
        DAMP = R*DFACT*(SIN(PHI)**10)
        EXPA = EXP(CI*Q*XFACTOR)*EXP(-DAMP*Q)

        DO INNERVAL = TOTPTS-NEND+1,TOTPTS

            QMID = MIDKVALS(INNERVAL) ! k_0
            INNERINTMID = MIDGREENS(:,:,INNERVAL)
            CALL GTRANSWOEXP(INNERINTMID,PHI,PI,CI,QMID) ! f(k_0)
            Q = MAINKVALS(INNERVAL+1) ! b
            INNERINTNEW = MAINGREENS(:,:,INNERVAL+1)
            CALL GTRANSWOEXP(INNERINTNEW,PHI,PI,CI,Q) ! f(b)

            KSTEP = Q-MAINKVALS(INNERVAL)
            EXPB = EXP(CI*Q*XFACTOR)*EXP(-DAMP*Q)
            CALL EXACTCALCINT(CI,INNERINTPREV,INNERINTMID,INNERINTNEW,XFACTOR,KSTEP,EXPA,EXPB,EXACTCALCMAT,DAMP)

            OUTERINTPREV = OUTERINTPREV + EXACTCALCMAT
            EXPA = EXPB

        END DO

        ! This begins the φ-loop
        DO J = 1, TOTPHISTEPS

            ! Calculation of Int(φ + 0.5dφ)
            OUTERINTMID = (0.D0,0.D0)

            PHI = PHI + 0.5D0*PHISTEP ! Corresponds to φ + 0.5dφ

            Q = MAINKVALS(1)
            INNERINTPREV = MAINGREENS(:,:,1)
            CALL GTRANSWEXP(INNERINTPREV,PHI,PI,CI,Q,R,DFACT)

            ! This begins the k loop for φ+0.5dφ
            ! First part
            DO INNERVAL = 1, TOTPTS-NEND

                QMID = MIDKVALS(INNERVAL) ! k + 0.5dk
                INNERINTMID = MIDGREENS(:,:,INNERVAL) ! G(k+0.5dk)
                CALL GTRANSWEXP(INNERINTMID,PHI,PI,CI,QMID,R,DFACT) ! I(φ+0.5dφ,k+0.5dk)
                Q = MAINKVALS(INNERVAL+1) ! k + dk
                INNERINTNEW = MAINGREENS(:,:,INNERVAL+1) ! G(k+dk)
                CALL GTRANSWEXP(INNERINTNEW,PHI,PI,CI,Q,R,DFACT) ! I(φ+0.5dφ,k+0.5dk)

                KSTEP = Q-MAINKVALS(INNERVAL)
                OUTERINTMID = OUTERINTMID + (KSTEP/6.D0)*(INNERINTPREV + 4.D0*INNERINTMID + INNERINTNEW)
                INNERINTPREV = INNERINTNEW

            END DO

            ! Second part
            INNERINTPREV = MAINGREENS(:,:,TOTPTS-NEND+1)
            CALL GTRANSWOEXP(INNERINTPREV,PHI,PI,CI,Q)
            
            XFACTOR = R*COS(PHI)
            DAMP = R*DFACT*(SIN(PHI)**10)
            EXPA = EXP(CI*Q*XFACTOR)*EXP(-DAMP*Q)

            DO INNERVAL = TOTPTS-NEND+1,TOTPTS

                QMID = MIDKVALS(INNERVAL) ! k_0
                INNERINTMID = MIDGREENS(:,:,INNERVAL)
                CALL GTRANSWOEXP(INNERINTMID,PHI,PI,CI,QMID) ! f(k_0)
                Q = MAINKVALS(INNERVAL+1) ! b
                INNERINTNEW = MAINGREENS(:,:,INNERVAL+1)
                CALL GTRANSWOEXP(INNERINTNEW,PHI,PI,CI,Q) ! f(b)

                KSTEP = Q-MAINKVALS(INNERVAL)
                EXPB = EXP(CI*Q*XFACTOR)*EXP(-DAMP*Q)
                CALL EXACTCALCINT(CI,INNERINTPREV,INNERINTMID,INNERINTNEW,XFACTOR,KSTEP,EXPA,EXPB,EXACTCALCMAT,DAMP)

                OUTERINTMID = OUTERINTMID + EXACTCALCMAT
                EXPA = EXPB

            END DO

            ! Calculation of Int(φ + dφ)
            OUTERINTNEW = (0.D0,0.D0)

            PHI = PHI + 0.5D0*PHISTEP ! Corresponds to φ + dφ

            Q = MAINKVALS(1)
            INNERINTPREV = MAINGREENS(:,:,1)
            CALL GTRANSWEXP(INNERINTPREV,PHI,PI,CI,Q,R,DFACT)

            ! This begins the k loop for φ+dφ
            ! First part
            DO INNERVAL = 1, TOTPTS-NEND

                QMID = MIDKVALS(INNERVAL) ! k + 0.5dk
                INNERINTMID = MIDGREENS(:,:,INNERVAL) ! G(k+0.5dk)
                CALL GTRANSWEXP(INNERINTMID,PHI,PI,CI,QMID,R,DFACT) ! I(φ+dφ,k+0.5dk)
                Q = MAINKVALS(INNERVAL+1) ! k + dk
                INNERINTNEW = MAINGREENS(:,:,INNERVAL+1) ! G(k+dk)
                CALL GTRANSWEXP(INNERINTNEW,PHI,PI,CI,Q,R,DFACT) ! I(φ+dφ,k+0.5dk)

                KSTEP = Q-MAINKVALS(INNERVAL)
                OUTERINTNEW = OUTERINTNEW + (KSTEP/6.D0)*(INNERINTPREV + 4.D0*INNERINTMID + INNERINTNEW)
                INNERINTPREV = INNERINTNEW

            END DO

            ! Second part
            INNERINTPREV = MAINGREENS(:,:,TOTPTS-NEND+1)
            CALL GTRANSWOEXP(INNERINTPREV,PHI,PI,CI,Q)
            
            XFACTOR = R*COS(PHI)
            DAMP = R*DFACT*(SIN(PHI)**10)
            EXPA = EXP(CI*Q*XFACTOR)*EXP(-DAMP*Q)

            DO INNERVAL = TOTPTS-NEND+1,TOTPTS

                QMID = MIDKVALS(INNERVAL) ! k_0
                INNERINTMID = MIDGREENS(:,:,INNERVAL)
                CALL GTRANSWOEXP(INNERINTMID,PHI,PI,CI,QMID) ! f(k_0)
                Q = MAINKVALS(INNERVAL+1) ! b
                INNERINTNEW = MAINGREENS(:,:,INNERVAL+1)
                CALL GTRANSWOEXP(INNERINTNEW,PHI,PI,CI,Q) ! f(b)

                KSTEP = Q-MAINKVALS(INNERVAL)
                EXPB = EXP(CI*Q*XFACTOR)*EXP(-DAMP*Q)
                CALL EXACTCALCINT(CI,INNERINTPREV,INNERINTMID,INNERINTNEW,XFACTOR,KSTEP,EXPA,EXPB,EXACTCALCMAT,DAMP)

                OUTERINTNEW = OUTERINTNEW + EXACTCALCMAT
                EXPA = EXPB

            END DO

            ! Simpson's rule for the φ-integral
            INTEGRAL = INTEGRAL + (PHISTEP/6.0)*(OUTERINTPREV + 4.0*OUTERINTMID + OUTERINTNEW)

            OUTERINTPREV = OUTERINTNEW

        END DO

        ! At the moment, only the 1-1 and 2-2 elements are printed
        write (1,'(F17.8, A, F17.8, A, F17.8, A, F17.8, A, F17.8, A, F17.8, A, F17.8, A, F17.8, A, F17.8)') R, ',', &
        &REAL(INTEGRAL(1,1)), ',', AIMAG(INTEGRAL(1,1)), ',', REAL(INTEGRAL(2,2)), ',', AIMAG(INTEGRAL(2,2)), ',', &
        &REAL(INTEGRAL(3,3)), ',', AIMAG(INTEGRAL(3,3)), ',', REAL(INTEGRAL(4,4)), ',', AIMAG(INTEGRAL(4,4))

        print *, 'Finished running for r =', R

    END DO

    close(1)

    !------------------------------------------------------------------------------------------------------------------

    contains

    subroutine GREENKBUILD(HBAR,M,MU,DELTA,ARASHBA,BETA,ENERGY,ETA,Q,GREENK)
        implicit none

        real*8, intent(in) :: HBAR, M, ENERGY, MU, DELTA, ARASHBA, BETA, ETA, Q
        integer :: IPIV(4), INFO
        real*8 :: SINGLEEN, WORK(16)
        complex*16 :: ZETA, GREENK(4,4)

        SINGLEEN = ((HBAR**2)/(2.D0*M))*(Q**2) - MU
        ZETA = DCMPLX(ENERGY,ETA)

        GREENK(:,:) = (0.D0,0.D0) ! Initialized

        ! Setup of z-H(k)
        GREENK(1,1) = ZETA - SINGLEEN - BETA
        GREENK(1,2) = ARASHBA*Q
        GREENK(1,4) = -DELTA

        GREENK(2,1) = GREENK(1,2)
        GREENK(2,2) = ZETA - SINGLEEN + BETA
        GREENK(2,3) = DELTA

        GREENK(3,2) = GREENK(2,3)
        GREENK(3,3) = ZETA + SINGLEEN + BETA
        GREENK(3,4) = GREENK(1,2)

        GREENK(4,1) = GREENK(1,4)
        GREENK(4,3) = GREENK(1,2)
        GREENK(4,4) = ZETA + SINGLEEN - BETA
        
        ! Compute an LU factorization of GREENK
        CALL ZGETRF(4,4,GREENK,4,IPIV,INFO)
        ! Invert GREENK, i.e. get g(k;z)
        CALL ZGETRI(4,GREENK,4,IPIV,WORK,16,INFO)

    end subroutine GREENKBUILD

    subroutine GTRANSWEXP(GREENK,PHI,PI,CI,Q,R,DFACT)
        implicit none
        real*8, intent(in) :: PHI, PI, Q, R, DFACT
        complex*16, intent(in) :: CI
        real*8 :: GAMMA, SINFACT, FOURCONST
        complex*16 :: GREENK(4,4), MINEXP, PLUSEXP, EXPON

        GAMMA = PHI - 0.5D0*PI

        FOURCONST = 1.D0/((2.D0*PI)**2)

        SINFACT = SIN(PHI)**10

        MINEXP = EXP(-CI*GAMMA)
        PLUSEXP = EXP(CI*GAMMA)

        GREENK(1,2) = GREENK(1,2)*MINEXP
        GREENK(1,3) = GREENK(1,3)*MINEXP
        GREENK(2,1) = GREENK(2,1)*PLUSEXP
        GREENK(2,4) = GREENK(2,4)*PLUSEXP
        GREENK(3,1) = GREENK(3,1)*PLUSEXP
        GREENK(3,4) = GREENK(3,4)*PLUSEXP
        GREENK(4,2) = GREENK(4,2)*MINEXP
        GREENK(4,3) = GREENK(4,3)*MINEXP

        EXPON = EXP(CI*Q*R*COS(PHI))*EXP(-Q*R*DFACT*SINFACT)

        ! The Fourier integrand is g(k,φ;z)*exp(ikrcosφ)*k
        GREENK = GREENK*EXPON*Q*FOURCONST

    end subroutine GTRANSWEXP

    subroutine EXACTCALCINT(CI,FOFA,FOFKO,FOFB,XFACTOR,KSTEP,EXPA,EXPB,EXACTCALCMAT,DAMP)
        implicit none
        real*8, intent(in) :: KSTEP, XFACTOR, DAMP
        complex*16, intent(in) :: CI, FOFA(4,4), FOFB(4,4), FOFKO(4,4), EXPA, EXPB
        complex*16 :: EXACTCALCMAT(4,4), DENOM, FIRSTTERM(4,4), SECONDTERM(4,4)

        EXACTCALCMAT = (0.D0,0.D0)
        DENOM = (CI*XFACTOR)-DAMP

        FIRSTTERM = FOFB/DENOM - (3.D0*FOFB + FOFA - 4.D0*FOFKO)/(KSTEP*DENOM*DENOM) + &
        &4.D0*(FOFB+FOFA-2.D0*FOFKO)/(DENOM*DENOM*DENOM*(KSTEP**2))

        SECONDTERM = FOFA/DENOM + (3.D0*FOFA + FOFB - 4.D0*FOFKO)/(KSTEP*DENOM*DENOM) + &
        &4.D0*(FOFB+FOFA-2.D0*FOFKO)/(DENOM*DENOM*DENOM*(KSTEP**2))

        EXACTCALCMAT = EXPB*FIRSTTERM - EXPA*SECONDTERM

    end subroutine EXACTCALCINT

    subroutine GTRANSWOEXP(GREENK,PHI,PI,CI,Q)
        implicit none
        real*8, intent(in) :: PHI, PI, Q
        complex*16, intent(in) :: CI
        real*8 :: GAMMA, FOURCONST
        complex*16 :: GREENK(4,4), MINEXP, PLUSEXP

        GAMMA = PHI - 0.5D0*PI

        FOURCONST = 1.D0/((2.D0*PI)**2)

        MINEXP = EXP(-CI*GAMMA)
        PLUSEXP = EXP(CI*GAMMA)

        GREENK(1,2) = GREENK(1,2)*MINEXP
        GREENK(1,3) = GREENK(1,3)*MINEXP
        GREENK(2,1) = GREENK(2,1)*PLUSEXP
        GREENK(2,4) = GREENK(2,4)*PLUSEXP
        GREENK(3,1) = GREENK(3,1)*PLUSEXP
        GREENK(3,4) = GREENK(3,4)*PLUSEXP
        GREENK(4,2) = GREENK(4,2)*MINEXP
        GREENK(4,3) = GREENK(4,3)*MINEXP

        ! The function we need to calculate is now g(k,φ;z)*k
        GREENK = GREENK*Q*FOURCONST

    end subroutine GTRANSWOEXP

end program SC_SOC_INT