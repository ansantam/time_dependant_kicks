! NON-LINEAR ELEMENTS (1 TO 10)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  if((kz(i).eq.4.or.kz(i).eq.5).and.abs(el(i)).gt.pieni)            &
 &ed(i)=-1d0*ed(i)     

----------------------------------------------------------------
        elseif (kz1.eq.2.or.kz1.eq.5) then
   40     fokm=el(l)*ed(l)
          if(abs(fokm).le.pieni) goto 20
          if(kz1.eq.2) then
            ih1=1
            ih2=2
          else

 elseif (kz1.eq.4.or.kz1.eq.6) then
!-----------------------------------------------------------------------
!  SEKTORMAGNET
!  HORIZONTAL
!-----------------------------------------------------------------------
   60     fokm=el(l)*ed(l)
          if(abs(fokm).le.pieni) goto 20
          if(kz1.eq.4) then
            ih1=1
            ih2=2
          else


elseif (kz1.eq.3) then
!-----------------------------------------------------------------------
!  QUADRUPOLE
!  FOCUSSING
!-----------------------------------------------------------------------
   80   do 90 j=1,napx
            fok(j)=ekv(j,l)*oidpsv(j)
            aek(j)=abs(fok(j))
            hi(j)=sqrt(aek(j))
            fi(j)=el(l)*hi(j)
            if(fok(j).le.zero) then
+if crlibm
              al(1,1,j,l)=cos_rn(fi(j))
+ei
+if .not.crlibm
              al(1,1,j,l)=cos(fi(j))
+ei
+if crlibm
              hi1(j)=sin_rn(fi(j))
+ei
+if .not.crlibm
              hi1(j)=sin(fi(j))
+ei
              if(abs(hi(j)).le.pieni) then
                al(2,1,j,l)=el(l)
              else

 elseif (kz1.eq.7.or.kz1.eq.8) then
!-----------------------------------------------------------------------
!  COMBINED FUNCTION MAGNET HORIZONTAL
!  FOCUSSING
!-----------------------------------------------------------------------
  100     if(kz1.eq.7) then
            do 110 j=1,napx
              fokqv(j)=ekv(j,l)
  110       continue
            ih1=1
            ih2=2
          else
!  COMBINED FUNCTION MAGNET VERTICAL
            do 120 j=1,napx
              fokqv(j)=-ekv(j,l)
  120       continue
            ih1=2
            ih2=1
          endif
          do 130 j=1,napx
            wf(j)=ed(l)/dpsq(j)
!hr01       fok(j)=fokqv(j)/dpd(j)-wf(j)*wf(j)
            fok(j)=fokqv(j)/dpd(j)-wf(j)**2                              !hr01
            afok(j)=abs(fok(j))
            hi(j)=sqrt(afok(j))
            fi(j)=hi(j)*el(l)
            if(afok(j).le.pieni) then
!hr01         as(6,1,j,l)=-rvv(j)*el(l)/c2e3
              as(6,1,j,l)=((-1d0*rvv(j))*el(l))/c2e3                     !hr01
              as(6,2,j,l)=as(6,1,j,l)
+if rvet
              as(1,1,j,l)=el(l)*rvet(j)
+ei
+if .not.rvet
!hr01         as(1,1,j,l)=el(l)*(one-rvv(j))*c1e3
              as(1,1,j,l)=(el(l)*(one-rvv(j)))*c1e3                      !hr01
+ei
            endif

      elseif (kz1.eq.9) then
!-----------------------------------------------------------------------
!  EDGE FOCUSSING
!-----------------------------------------------------------------------
  140     do 150 j=1,napx
            rhoi(j)=ed(l)/dpsq(j)
+if crlibm
!hr01       fok(j)=rhoi(j)*tan_rn(el(l)*rhoi(j)*half)
            fok(j)=rhoi(j)*tan_rn((el(l)*rhoi(j))*half)                  !hr01
+ei
+if .not.crlibm
!hr01       fok(j)=rhoi(j)*tan(el(l)*rhoi(j)*half)
            fok(j)=rhoi(j)*tan((el(l)*rhoi(j))*half)                     !hr01
+ei
            al(3,1,j,l)=fok(j)
            al(3,2,j,l)=-fok(j)
  150     continue
          goto 160
        else
!Eric
! Is really an error but old code went to 160
          goto 160
        endif

! MULTIPOLES (11)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


if(kz(i).eq.11.and.abs(el(i)+one).le.pieni) then
  dki(i,1) = ed(i)
  dki(i,3) = ek(i)
  ed(i) = one
  ek(i) = one
  el(i) = zero
else if(kz(i).eq.11.and.abs(el(i)+two).le.pieni) then
  dki(i,2) = ed(i)
  dki(i,3) = ek(i)
  ed(i) = one
  ek(i) = one
  el(i) = zero
endif

if(kz(ix).eq.11) izu=izu+2*mmul
endif

 if(kz(ix).eq.11) izu=izu+2*mmul
          endif

! CAVITIES (12)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


if(abs(kz(i)).eq.12) then
    if(abs(ed(i)).gt.pieni.and.abs(ek(i)).gt.pieni) then
      ncy2=ncy2+1
      itionc(i)=kz(i)/abs(kz(i))
      kp(i)=6
    endif
    phasc(i)=el(i)
    el(i)=zero
endif

  if(abs(kz(k)).eq.12) then
+if cr
        write(lout,10070) k,bez(k),kz(k),ed(k),ek(k),phasc(k),xpl(k),   &
+ei
+if .not.cr
        write(*,10070) k,bez(k),kz(k),ed(k),ek(k),phasc(k),xpl(k),      &
+ei
     &xrms(k),zpl(k),zrms(k)
        kz(k)=abs(kz(k))
        phasc(k)=phasc(k)*rad
      else
+if cr


  if(abs(kz(ix)).eq.12) ncy=ncy+1
          endif
  725   continue
        do 730 j=1,il
          if(abs(kz(j)).eq.12) then
!hr05       hsyc(j)=two*pi*ek(j)/tlen
            hsyc(j)=((two*pi)*ek(j))/tlen                                !hr05
            if(nvar.eq.5) then
              ition=1
              ed(j)=zero
            endif
          endif
  730   continue
      endif
      goto 110

    if(kz(jj).eq.12) hsyc(jj)=c1m3*hsyc(jj)*itionc(jj)
          if(kz(jj).eq.12) hsyc(jj)=(c1m3*hsyc(jj))*dble(itionc(jj))     !hr01

            if(kz(ix).eq.12) then
+if crlibm
              ejv(j)=ejv(j)+ed(ix)*sin_rn(hsyc(ix)*sigmv(j)+
+ei
+if .not.crlibm
              ejv(j)=ejv(j)+ed(ix)*sin(hsyc(ix)*sigmv(j)+               &
+ei
     &phasc(ix))
            else
+if crlibm
              ejv(j)=ejv(j)+hsy(1)*sin_rn(hsy(3)*sigmv(j))
+ei
+if .not.crlibm
              ejv(j)=ejv(j)+hsy(1)*sin(hsy(3)*sigmv(j))
+ei
            endif

if(kz(ix).eq.12) then
+if crlibm
!hr01         ejv(j)=ejv(j)+ed(ix)*sin_rn(hsyc(ix)*sigmv(j)+phas+
!hr01&phasc(ix))
              ejv(j)=ejv(j)+ed(ix)*sin_rn((hsyc(ix)*sigmv(j)+phas)+     &!hr01
     &phasc(ix))                                                         !hr01
+ei
+if .not.crlibm
!hr01         ejv(j)=ejv(j)+ed(ix)*sin(hsyc(ix)*sigmv(j)+phas+
!hr01&phasc(ix))
              ejv(j)=ejv(j)+ed(ix)*sin((hsyc(ix)*sigmv(j)+phas)+        &
     &phasc(ix))                                                         !hr01
+ei
            else
+if crlibm
              ejv(j)=ejv(j)+hsy(1)*sin_rn(hsy(3)*sigmv(j)+phas)
+ei
+if .not.crlibm
              ejv(j)=ejv(j)+hsy(1)*sin(hsy(3)*sigmv(j)+phas)
+ei
            endif

! AC DIPOLES (16)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


  if(abs(kz(i)).eq.16) then
    if(abs(ed(i)).le.pieni) then
       kz(i)=0
!hr05      ed(i)=0
       ed(i)=0d0                                                     !hr05
!hr05      ek(i)=0
       ek(i)=0d0                                                     !hr05
!hr05      el(i)=0
       el(i)=0d0                                                     !hr05
    else
       acdipph(i)=el(i)
!hr05      el(i)=0
       el(i)=0d0                                                     !hr05
    endif
  endif

!----Insertion for AC dipole
        if(abs(kz(j)).eq.16) then
          nturn1(j)=int(xpl0)
          nturn2(j)=int(xrms0)
          nturn3(j)=int(zpl0)
          nturn4(j)=int(zrms0)
          xpl(j)=0d0
          xrms(j)=0d0
          zpl(j)=0d0
          zrms(j)=0d0
!hr05     if(xrms0.eq.0.and.zpl0.eq.0.and.zrms0.eq.0) then
          if(xrms0.eq.0d0.and.zpl0.eq.0d0.and.zrms0.eq.0d0) then         !hr05
+if cr
            write(lout,*) "ac dipole disregarded (0 length)"
+ei
+if .not.cr
            write(*,*) "ac dipole disregarded (0 length)"
+ei
            kz(j)=0
!hr05       ed(j)=0
            ed(j)=0d0                                                    !hr05
!hr05       ek(j)=0
            ek(j)=0d0                                                    !hr05
          endif
        endif

! BEAM-BEAM (20)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


  if(kz(i).eq.20) then
    ptnfac(i)=el(i)
    el(i)=zero
  endif
  if(abs(el(i)).gt.pieni.and.kz(i).ne.0) ithick=1
  if(i.gt.nele-1) call prror(16)
  if(abs(kz(i)).ne.12) kp(i)=0
  bez(i)=idat
  bez0(i)=idat
  if(ncy2.eq.0) then
    i=i+1
    il=i
    bez(i)=cavi
    bez0(i)=cavi
    kp(i)=6
  else
    il=i
    i=i+1
  endif
  goto 130

! CRABS (23, 26, 27 28)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!--CRABCAVITY
      if(abs(kz(i)).eq.23) then
        if(abs(ed(i)).le.pieni) then
           kz(i)=0
!hr05      ed(i)=0
           ed(i)=0d0                                                     !hr05
!hr05      ek(i)=0
           ek(i)=0d0                                                     !hr05
!hr05      el(i)=0
           el(i)=0d0                                                     !hr05
        else
           crabph(i)=el(i)
!hr05      el(i)=0
           el(i)=0d0                                                     !hr05
        endif
      endif

! JBG RF CC Multipoles
!--CC Mult kick order 2
      if(abs(kz(i)).eq.26) then
        if(abs(ed(i)).le.pieni) then
           kz(i)=0
!hr05      ed(i)=0
           ed(i)=0d0                                                     !hr05
!hr05      ek(i)=0
           ek(i)=0d0                                                     !hr05
!hr05      el(i)=0
           el(i)=0d0                                                     !hr05
        else
           crabph2(i)=el(i)
!hr05      el(i)=0
           el(i)=0d0                                                     !hr05
        endif
      endif

!--CC Mult kick order 3
      if(abs(kz(i)).eq.27) then
        if(abs(ed(i)).le.pieni) then
           kz(i)=0
!hr05      ed(i)=0
           ed(i)=0d0                                                     !hr05
!hr05      ek(i)=0
           ek(i)=0d0                                                     !hr05
!hr05      el(i)=0
           el(i)=0d0                                                     !hr05
        else
           crabph3(i)=el(i)
!hr05      el(i)=0
           el(i)=0d0                                                     !hr05
        endif
      endif
!--CC Mult kick order 4
      if(abs(kz(i)).eq.28) then
        if(abs(ed(i)).le.pieni) then
           kz(i)=0
!hr05      ed(i)=0
           ed(i)=0d0                                                     !hr05
!hr05      ek(i)=0
           ek(i)=0d0                                                     !hr05
!hr05      el(i)=0
           el(i)=0d0                                                     !hr05
        else
           crabph4(i)=el(i)
!hr05      el(i)=0
           el(i)=0d0                                                     !hr05
        endif
      endif

! type 25 not implemented yet
  if(kz(i).eq.25) then
    ed(i)=ed(i)/two
    ek(i)=ek(i)/two
  endif



 if(kz(j).ne.8) elbe(i)=elbe(i)+el(j)
  300 continue
  310 k0=l-1
      goto 220


 if(ilm0(j1).eq.bez(j2)) then
              if(el(j2).ne.zero.or.kz(j2).gt.10) call prror(67)
              ipar(j1)=j2
              goto 540
            endif
  530     continue
          call prror(66)
  540   continue










