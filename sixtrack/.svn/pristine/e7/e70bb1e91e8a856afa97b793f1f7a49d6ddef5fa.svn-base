module mod_beam
  implicit none

  ! DIST input block present?
  logical, public :: dist_enable = .false.

  ! dump a detailed echo of the read particle distribution
  logical, public :: dist_echo = .false.

  integer, parameter :: outlun = 54

  character(len=20) :: beam_type

  character(len=15) :: beam_filename

  double precision :: beam_emitx0, &
                      beam_emity0, &
                      beam_enom, &
                      beam_nex, &
                      beam_ney, &
                      beam_dex, &
                      beam_dey, &
                      beam_enerror, &
                      beam_bunchsize, &
                      beam_nr, &
                      beam_ndr

  contains

  !----------------------------------------------------------------------------
  subroutine beam_init(npart, maxn, alphax, alphay, betax, betay, &
                      enom, pnom, c, x, y, xp, yp, s, pc)
    integer :: npart, maxn
    double precision :: alphax, &
                        alphay, &
                        betax, &
                        betay, &
                        enom, &
                        pnom, &
                        c
    double precision :: x(maxn),  &
                        y(maxn),  &
                        xp(maxn), &
                        yp(maxn), &
                        s(maxn),  &
                        pc(maxn)

    integer :: j

    !if(beam_type.eq.0) then
      call beam_readdis(beam_filename, npart, maxn, enom, pnom, c, x, y, xp, yp, s, pc)
    !else if(beam_type.eq.4) then
    !  call beam_makedis_radial(npart, alphax, alphay, betax, betay, &
    !            beam_emitx0, beam_emity0, beam_enom, beam_nr, beam_ndr, &
    !            x, xp, y, yp, pc, s)
    !else if(beam_type.eq.3) then
    !  call beam_makedis_de(npart, alphax, alphay, betax, betay, &
    !            beam_emitx0, beam_emity0, beam_enom, beam_nex, beam_ney, &
    !            beam_dex, beam_dey, beam_enerror, beam_bunchsize, &
    !            x, xp, y, yp, pc, s)
    !else if(beam_type.eq.2) then
    !  call beam_makedis_st(npart, alphax, alphay, betax, betay, &
    !            beam_emitx0, beam_emity0, beam_enom, beam_nex, beam_ney, &
    !            beam_dex, beam_dey, x, xp, y, yp, pc, s)
    !else 
    !  call beam_makedis(npart, alphax, alphay, betax, betay, &
    !            beam_emitx0, beam_emity0, beam_enom, beam_nex, beam_ney, &
    !            beam_dex, beam_dey, x, xp, y, yp, pc, s)
    !end if

  end subroutine beam_init

  !----------------------------------------------------------------------------
  subroutine beam_makedis(np, maxn, alphax, alphay, betax, betay, &
    emitx0, emity0, enom, nex, ney, dex, dey, &
    x, xp, y, yp, p, s)
    !
    !  Generate distribution
    !

    integer :: np, maxn
    double precision :: alphax,   &
                        alphay,   &
                        betax,    &
                        betay,    &
                        emitx0,   &
                        emity0,   &
                        enom,     &
                        nex,      &
                        dex,      &
                        ney,      &
                        dey,      &
                        x(maxn),  &
                        xp(maxn), &
                        y(maxn),  &
                        yp(maxn), &
                        p(maxn),  &
                        s(maxn)

    double precision :: pi
    double precision :: emitx,    &
                        gammax,   &
                        emity,    &
                        gammay,   &
                        xsigmax,  &
                        ysigmay

    integer :: j
    !---------------------------------------------------------------------
    !+ Generate particle distribution
    !+ Generate random distribution, assuming optical parameters at IP1


    !+ Calculate the gammas

    pi=4d0*atan(1d0)
    gammax = (1d0+alphax**2)/betax
    gammay = (1d0+alphay**2)/betay

    !+ Number of points and generate distribution

    write(*,*)
    write(*,*) 'Generation of particle distribution Version 1'
    write(*,*)
    write(*,*) 'This routine generates particles in phase space'
    write(*,*) 'X/XP and Y/YP ellipses, as defined in the input'
    write(*,*) 'parameters. Distribution is flat in the band.'
    write(*,*) 'X and Y are fully uncorrelated.'
    write(*,*)

    write(outlun,*)
    write(outlun,*) 'Generation of particle distribution Version 1'
    write(outlun,*)
    write(outlun,*) 'This routine generates particles in phase space'
    write(outlun,*) 'X/XP and Y/YP ellipses, as defined in the input'
    write(outlun,*) 'parameters. Distribution is flat in the band.'
    write(outlun,*) 'X and Y are fully uncorrelated.'
    write(outlun,*)
    write(outlun,*) 'INFO>  Number of particles   = ', np
    write(outlun,*) 'INFO>  Av number of x sigmas = ', nex
    write(outlun,*) 'INFO>  +- spread in x sigmas = ', dex
    write(outlun,*) 'INFO>  Av number of y sigmas = ', ney
    write(outlun,*) 'INFO>  +- spread in y sigmas = ', dey
    write(outlun,*) 'INFO>  Nominal beam energy   = ', enom
    write(outlun,*) 'INFO>  Sigma_x0 = ', sqrt(betax*emitx0)
    write(outlun,*) 'INFO>  Sigma_y0 = ', sqrt(betay*emity0)
    write(outlun,*) 'INFO>  Beta x   = ', betax
    write(outlun,*) 'INFO>  Beta y   = ', betay
    write(outlun,*) 'INFO>  Alpha x  = ', alphax
    write(outlun,*) 'INFO>  Alpha y  = ', alphay
    write(outlun,*)

    do j = 1, np

      emitx = emitx0*(nex + ((2d0*dble(rand()-0.5))*dex) )**2
      xsigmax = sqrt(betax*emitx)
      x(j) = xsigmax * sin((2d0*pi)*dble(rand()))

      if (rand().gt.0.5) then
        xp(j) = +1d0*sqrt(emitx/betax-x(j)**2/betax**2)- &
                  (alphax*x(j))/betax
      else
        xp(j) = -1d0*sqrt(emitx/betax-x(j)**2/betax**2)- &
                  (alphax*x(j))/betax
      endif

      emity = emity0*(ney + ((2d0*dble(rand()-0.5))*dey) )**2
      ysigmay = sqrt(betay*emity)
      y(j) = ysigmay * sin((2d0*pi)*dble(rand()))

      if (rand().gt.0.5) then
        yp(j) = +1d0*sqrt(emity/betay-y(j)**2/betay**2)- &
                  (alphay*y(j))/betay
      else
        yp(j) = -1d0*sqrt(emity/betay-y(j)**2/betay**2)- &
                  (alphay*y(j))/betay
      endif

      p(j) = enom
      s(j) = 0d0

    end do

  end subroutine beam_makedis

  subroutine beam_makedis_de(np, maxn, alphax, alphay, betax,     &
    betay, emitx0, emity0, enom, nex, ney, dex, dey, &
    enerror, bunchlen, &
    x, xp, y, yp, p, s)
    !
    !  Generate distribution
    !

    integer :: np, maxn
    double precision :: alphax,   &
                        alphay,   &
                        betax,    &
                        betay,    &
                        emitx0,   &
                        emity0,   &
                        enom,     &
                        nex,      &
                        dex,      &
                        ney,      &
                        dey,      &
                        x(maxn),  &
                        xp(maxn), &
                        y(maxn),  &
                        yp(maxn), &
                        p(maxn),  &
                        s(maxn)

    double precision :: pi
    double precision :: emitx,    &
                        gammax,   &
                        emity,    &
                        gammay,   &
                        xsigmax,  &
                        ysigmay,  &
                        enerror,  &
                        bunchlen

    integer :: j

!     Uses the old routine 'MAKEDIS' for the halo plane and adds the
!     transverse beam size in the other plane (matched distributions
!     are generated starting from thetwiss functions).
!     If 'nex' and 'ney' are BOTH set to zero, nominal bunches
!     centred in the aperture centre are generated. (SR, 08-05-2005)

      double precision :: iix, iiy, phix, phiy
      double precision :: en_error, bunch_len
!
      double precision :: long_cut
      double precision :: a_st, b_st

!
!-----------------------------------------------------------------------
!++  Generate particle distribution
!
!
!++  Generate random distribution, assuming optical parameters at IP1
!++  Calculate the gammas

      pi=4d0*atan(1d0)
      gammax = (1d0+alphax**2)/betax
      gammay = (1d0+alphay**2)/betay

      en_error = enerror
      bunch_len = bunchlen

      write(outlun,*) "Generation of bunch with dp/p and length:"
      write(outlun,*) "  RMS bunch length  = ", bunch_len
      write(outlun,*) "  RMS energy spread = ", en_error

      do j=1, np

         if ((nex.gt.0d0).and.(ney.eq.0d0)) then

            emitx = emitx0*(nex+((2d0*dble(rand()-0.5))*dex))**2
            xsigmax = sqrt(betax*emitx)
            x(j) = xsigmax * sin((2d0*pi)*dble(rand()))

            if (rand().gt.0.5) then
              xp(j) = sqrt(emitx/betax-x(j)**2/betax**2) - (alphax*x(j))/betax
            else
              xp(j) = -1d0*sqrt(emitx/betax-x(j)**2/betax**2)- (alphax*x(j))/betax
            endif

            phiy = (2d0*pi)*dble(rand())
            iiy = (-1d0*emity0) * log( dble(rand()) )

            y(j) = sqrt((2d0*iiy)*betay) * cos(phiy)
            yp(j) = (-1d0*sqrt((2d0*iiy)/betay)) * (sin(phiy) + alphay * cos(phiy))

         elseif ( nex.eq.0d0.and.ney.gt.0d0 ) then

            emity = emity0*(ney+((2d0*dble(rand()-0.5))*dey))**2
            ysigmay = sqrt(betay*emity)
            y(j) = ysigmay * sin((2d0*pi)*dble(rand()))

            if (rand().gt.0.5) then
              yp(j) = sqrt(emity/betay-y(j)**2/betay**2) - (alphay*y(j))/betay
            else
              yp(j) = -1d0*sqrt(emity/betay-y(j)**2/betay**2) - (alphay*y(j))/betay
            endif

            phix = (2d0*pi)*dble(rand())
            iix = (-1d0*emitx0) * log( dble(rand()) )

            x(j) = sqrt((2d0*iix)*betax) * cos(phix)
            xp(j) = (-1d0*sqrt((2d0*iix)/betax)) * (sin(phix) + alphax * cos(phix))

         elseif ( nex.eq.0d0.and.ney.eq.0d0 ) then

            phix = (2d0*pi)*dble(rand())
            iix = (-1d0*emitx0) * log( dble(rand()) )

            x(j) = sqrt((2d0*iix)*betax) * cos(phix)
            xp(j) = (-1d0*sqrt((2d0*iix)/betax)) * (sin(phix) + alphax * cos(phix))
            phiy = (2d0*pi)*dble(rand())
            iiy = (-1d0*emity0) * log( dble(rand()) )

            y(j) = sqrt((2d0*iiy)*betay) * cos(phiy)
            yp(j) = (-1d0*sqrt((2d0*iiy)/betay)) * (sin(phiy) + alphay * cos(phiy))
         else
            write(*,*) "Error - beam parameters not correctly set!"
         endif

      end do

      ! SR, 11-08-2005 For longitudinal phase-space, add a cut at 2 sigma

      long_cut = 2
      j = 1

      ! 1st: generate npnumbers within the chose cut
      do while (j.le.np)
         a_st = ran_gauss(5d0)
         b_st = ran_gauss(5d0)
         do while ((a_st**2+b_st**2).gt.long_cut**2)
            a_st = ran_gauss(5d0)
            b_st = ran_gauss(5d0)
         enddo
         s(j) = a_st
         p(j) = b_st
         j = j + 1
      enddo

      ! 2nd: give the correct values
      do j=1,np
         p(j) = enom * (1d0 + p(j) * en_error)
         s(j) = bunch_len * s(j)
      enddo

      return
  end subroutine beam_makedis_de

  subroutine beam_makedis_st(np, maxn, alphax, alphay, betax, betay, &
    emitx0, emity0, enom, nex, ney, dex, dey,        &
    x, xp, y, yp, p, s)

!     Uses the old routine 'MAKEDIS' for the halo plane and adds the
!     transverse beam size in the other plane (matched distrubutions
!     are generated starting from thetwiss functions).
!     If 'nex' and 'ney' are BOTH set to zero, nominal bunches
!     centred in the aperture centre are generated. (SR, 08-05-2005)

    integer :: np, maxn
    double precision :: alphax,   &
                        alphay,   &
                        betax,    &
                        betay,    &
                        emitx0,   &
                        emity0,   &
                        enom,     &
                        nex,      &
                        dex,      &
                        ney,      &
                        dey,      &
                        x(maxn),  &
                        xp(maxn), &
                        y(maxn),  &
                        yp(maxn), &
                        p(maxn),  &
                        s(maxn)

    double precision :: emitx,    &
                        gammax,   &
                        emity,    &
                        gammay,   &
                        xsigmax,  &
                        ysigmay

    integer :: j

    double precision :: pi
    double precision :: iix, iiy, phix, phiy
    !
    !-----------------------------------------------------------------------
    !++  Generate particle distribution
    !
    !
    !++  Generate random distribution, assuming optical parameters at IP1
    !
    !++  Calculate the gammas
    write(outlun,*) '  New routine to add the finite beam size in the'
    write(outlun,*) '  other dimension (SR, 08-06-2005).'

    pi=4d0*atan(1d0)

    gammax = (1d0+alphax**2)/betax
    gammay = (1d0+alphay**2)/betay

    do j=1, np
      if ((nex.gt.0d0).and.(ney.eq.0d0)) then

        emitx = emitx0*(nex+((2d0*dble(rand()-0.5))*dex))**2
        xsigmax = sqrt(betax*emitx)
        x(j)   = xsigmax * sin((2d0*pi)*dble(rand()))

        if (rand().gt.0.5) then
          xp(j) = sqrt(emitx/betax-x(j)**2/betax**2)- &
                      (alphax*x(j))/betax
        else
          xp(j) = -1d0*sqrt(emitx/betax-x(j)**2/betax**2)-&
                      (alphax*x(j))/betax
        endif

        phiy = (2d0*pi)*dble(rand())
        iiy = (-1d0*emity0) * log( dble(rand()) )
        y(j) = sqrt((2d0*iiy)*betay) * cos(phiy)
        yp(j) = (-1d0*sqrt((2d0*iiy)/betay)) * (sin(phiy) +     &
                    alphay * cos(phiy))
      elseif ( nex.eq.0d0.and.ney.gt.0d0 ) then

        emity = emity0*(ney+((2d0*dble(rand()-0.5))*dey))**2
        ysigmay = sqrt(betay*emity)
        y(j)   = ysigmay * sin((2d0*pi)*dble(rand()))

        if (rand().gt.0.5) then
          yp(j) = sqrt(emity/betay-y(j)**2/betay**2)- &
                      (alphay*y(j))/betay
        else
          yp(j) = -1d0*sqrt(emity/betay-y(j)**2/betay**2)-&
                      (alphay*y(j))/betay
        endif

        phix = (2d0*pi)*dble(rand())
        iix = (-1d0* emitx0) * log( dble(rand()) )
        x(j) = sqrt((2d0*iix)*betax) * cos(phix)
        xp(j) = (-1d0*sqrt((2d0*iix)/betax)) * (sin(phix) +     &
                    alphax * cos(phix))
      elseif ( nex.eq.0d0.and.ney.eq.0d0 ) then

        phix = (2d0*pi)*dble(rand())
        iix = (-1d0*emitx0) * log( dble(rand()) )
        x(j) = sqrt((2d0*iix)*betax) * cos(phix)
        xp(j) = (-1d0*sqrt((2d0*iix)/betax)) * (sin(phix) +     &
                    alphax * cos(phix))
        phiy = (2d0*pi)*dble(rand())
        iiy = (-1d0*emity0) * log( dble(rand()) )

        y(j) = sqrt((2d0*iiy)*betay) * cos(phiy)
        yp(j) = (-1d0*sqrt((2d0*iiy)/betay)) * (sin(phiy) +     &
                    alphay * cos(phiy))
      else
        write(*,*) "Error - beam parameters not correctly set!"
      endif

    end do

    p = enom
    s = 0d0

    return
  end subroutine beam_makedis_st

  subroutine beam_makedis_radial(np, maxn, alphax, alphay, betax,      &
    betay, emitx0, emity0, enom, nr, ndr, x, xp, y,      &
    yp, p, s)

    !-----------------------------------------------------------------------
    !++  Generate particle distribution
    !
    !
    !++  Generate random distribution, assuming optical parameters at IP1
    !
    integer :: np, maxn
    double precision :: alphax,   &
                        alphay,   &
                        betax,    &
                        betay,    &
                        emitx0,   &
                        emity0,   &
                        enom,     &
                        nr,       &
                        ndr,      &
                        x(maxn),  &
                        xp(maxn), & 
                        y(maxn),  &
                        yp(maxn), & 
                        p(maxn),  &
                        s(maxn)

    double precision :: pi
    double precision :: emitx,    &
                        gammax,   &
                        emity,    &
                        gammay,   &
                        xsigmax,  &
                        ysigmay,  &
                        nex,      &
                        dex,      &
                        ney,      &
                        dey

    integer :: j

    !++  Calculate the gammas
    pi=4d0*atan(1d0)
    gammax = (1d0+alphax**2)/betax
    gammay = (1d0+alphay**2)/betay

    !++  Number of points and generate distribution
    nex = nr/sqrt(2d0)
    dex = ndr/sqrt(2d0)
    ney = nr/sqrt(2d0)
    dey = ndr/sqrt(2d0)

    write(outlun,*)
    write(outlun,*) 'Generation of particle distribution Version 2'
    write(outlun,*)
    write(outlun,*) 'This routine generates particles in that are fully'
    write(outlun,*) 'correlated between X and Y.'
    write(outlun,*)
    write(outlun,*)
    write(outlun,*) 'INFO>  Number of particles   = ', np
    write(outlun,*) 'INFO>  Av number of x sigmas = ', nex
    write(outlun,*) 'INFO>  +- spread in x sigmas = ', dex
    write(outlun,*) 'INFO>  Av number of y sigmas = ', ney
    write(outlun,*) 'INFO>  +- spread in y sigmas = ', dey
    write(outlun,*) 'INFO>  Nominal beam energy   = ', enom
    write(outlun,*) 'INFO>  Sigma_x0 = ', sqrt(betax*emitx0)
    write(outlun,*) 'INFO>  Sigma_y0 = ', sqrt(betay*emity0)
    write(outlun,*)

    do while (j.lt.np)
      j = j + 1

      emitx = emitx0*(nex + ((2d0*dble(rand()-0.5))*dex) )**2
      xsigmax = sqrt(betax*emitx)
      x(j)   = xsigmax * sin((2d0*pi)*dble(rand()))

      if (rand().gt.0.5) then
        xp(j)  = sqrt(emitx/betax-x(j)**2/betax**2)-        &
                     (alphax*x(j))/betax
      else
        xp(j)  = -1d0*sqrt(emitx/betax-x(j)**2/betax**2)-   &
                     (alphax*x(j))/betax
      endif

      emity = emity0*(ney + ((2d0*dble(rand()-0.5))*dey) )**2
      ysigmay = sqrt(betay*emity)
      y(j)   = ysigmay * sin((2d0*pi)*dble(rand()))

      if (rand().gt.0.5) then
        yp(j)  = sqrt(emity/betay-y(j)**2/betay**2)-        &
                     (alphay*y(j))/betay
      else
        yp(j)  = -1d0*sqrt(emity/betay-y(j)**2/betay**2)-   &
                     (alphay*y(j))/betay
      endif

      p(j)   = enom
      s(j)   = 0d0

    end do

    return
  end subroutine beam_makedis_radial

  subroutine beam_readdis(filename, np, maxn, enom, pnom, c, &
                          x, y, xp, yp, s, pc)
  ! Format for the input file:
  !   id: Unique identifier of the particle (integer)
  !   gen: parent ID
  !   weight: Importance of the particle in the simulation (real: >0.0)
  !   x, y, s: Particle position (m)
  !   xp, yp, zp: Particle direction (tangents)
  !   aa, zz: Particle atomic number and charge
  !   m: Rest mass (GeV/c^2)
  !   pc: Particle momentum (GeV)
  !   dt: Time delay with respect to the reference particle (seconds)

    implicit none

    character(len=15) :: filename
    integer :: np, maxn
    double precision :: enom, pnom, c
    double precision :: x(np),  &
                        y(np),  &
                        xp(np), &
                        yp(np), &
                        s(np),  &
                        pc(np)

    integer :: j
    integer :: ierro

    integer :: id
    integer :: gen
    integer :: aa, zz
    double precision :: weigth
    double precision :: m(np)
    double precision :: z(np), zp(np), dt(np)

    write(*,*) "Reading particles from ", filename

    open(unit=outlun, file=filename)

    do j=1,np
      read(outlun,*,iostat=ierro) id, gen, weigth, &
           x(j), y(j), z(j), xp(j), yp(j), zp(j), aa, zz, m(j), pc(j), dt(j)
      if(ierro.ne.0) then
        write(*, *) 'Warning: Error reading particles'
        exit
      endif
    enddo

    np = j - 1
    write(*,*) "Number of particles read = ", np

    x(1:np) = x(1:np) * 1d3
    y(1:np) = y(1:np) * 1d3

    xp(1:np) = xp(1:np) * 1d3
    yp(1:np) = yp(1:np) * 1d3

    pc(1:np) = pc(1:np) * 1d3
    m(1:np)  =  m(1:np) * 1d3

    s(1:np) = - pnom/enom * dt(1:np)*c * 1d3

    close(outlun)

    return
  end subroutine beam_readdis

  double precision function ran_gauss(cut)
!*********************************************************************
!
! RAN_GAUSS - will generate a normal distribution from a uniform
!   distribution between [0,1].
!   See "Communications of the ACM", V. 15 (1972), p. 873.
!
! cut - double precision - cut for distribution in units of sigma
!                the cut must be greater than 0.5
!
!*********************************************************************
    logical :: flag
    double precision :: x, u1, u2, twopi, r,cut

    twopi=8d0*atan(1d0)
1   if (flag) then
      r = dble(rand())
      r = max(r, 0.5d0**32)
      r = min(r, 1d0-0.5d0**32)
      u1 = sqrt(-2d0*log( r ))
      u2 = dble(rand( ))
      x = u1 * cos(twopi*u2)
    else
      x = u1 * sin(twopi*u2)
    endif

    flag = .not. flag

    !  cut the distribution if cut > 0.5
    if (cut .gt. 0.5d0 .and. abs(x) .gt. cut) goto 1

    ran_gauss = x
    return
  end function ran_gauss
end module
