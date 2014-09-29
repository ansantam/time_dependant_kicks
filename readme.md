## New function for Crab Cavities with [DYNK](https://twiki.cern.ch/twiki/bin/view/FlukaTeam/CouplingSVNRepositories) ##

### About the `DYNK` block:

```
!     A.Mereghetti, for the FLUKA Team
!     last modified: 03-09-2014
!     COMMON for dynamic kicks
!     always in main code
 
!     in case the DYNK input block is issued, the kick of selected SINGLE
!       ELEMENTs (and all their entries in the accelerator sequence)
!       is modulated turn by turn, according user's specifications
 
!     the user defines a set of basic functions, with their parameters
!     then, for each SINGLE ELEMENT, the user declares how these functions
!       should be combined in order to get the actual profile, and their
!       sequence, including turn numbers
```
 
* `smiv` in the `DYNK` block:

```
      subroutine saveorigsmiv
!
!-----------------------------------------------------------------------
!     A.Mereghetti, for the FLUKA Team
!     last modified: 03-09-2014
!     save original values of smiv variable
!     always in main code
!
!     the idea is to loop over all those SINGLE ELEMENTs flagged for
!         dynamic kicks, and save the original value of smiv;
!     at the same time, check that all the concerned entries in
!         the accelerator structure have the same smiv value;
!
!     NB: originally, the smiv variable is used only in case of non-linear
!         SINGLE ELEMENTs: thus, dynamic kicks can be applied only to
!         non-linear SINGLE ELEMENTs
```

### Error message

When using a lattice with crab cavities installed, and trying to use the DYNK module on a crab element we get the following error:

```
 CALL TO SAVEORIGSMIV 
     ...with checks of possible incosistencies 
 
 
  inconsistency in accelerator structure:
  for SINGLE ELEMENT CRAB5           
  smiv at first occurrence:    541.17124490936453     
  smiv of          168  entry:    451.41929755283127     
 
 
  at least one incosistency in flagging elements
     for dynamic kicks: please check carefully...
 



         ++++++++++++++++++++++++
         +++++ERROR DETECTED+++++
         ++++++++++++++++++++++++
         RUN TERMINATED ABNORMALLY !!!
```

* `smiv(i)`: magnetic kick, composition of the _nominal kick_ `ek(ix)` and the random error.

The values of `smiv` are are assigned in the `maincr` program, in the `sixve.f`module, slightly after the call to the `daten` subroutine. 

* `strack(i)`: BLOC physical length [m] (most frequently referred to as `stracki` in `do` loops.

Values of `strack(i)` are equivalent to `smiv(1,i)`, apart from a power of 1000 depending on the type of non-linear `SINGLE ELEMENT`.


#### Why we get that error message
```
     loop over all the SINGLE ELEMENTs flagged for dynamic kicks:
      do kk=1,NacqDynkSEs
 
        ii=iSEDynks(kk)
 
!       loop over the entries in the accelerator structure:
!          get the value of smiv of the first entry in the structure,
!          and then check that all other entries have the same value
        lfirst = .true.
        do jj=1,iu
          if ( ktrack(jj).ne.1 ) then
!           a SINGLE ELEMENT (skip the blocs)
            if ( ic(jj)-nblo.eq.ii ) then
!             current entry is another instance of the selected
!               SINGLE ELEMENT
              if ( lfirst ) then
                oriSmivSEDynks(kk) = smiv(1,jj)
                lfirst = .false.
              else
                if ( smiv(1,jj).ne.oriSmivSEDynks(kk) ) then
                  write(*,*)''
                  write(*,*)' inconsistency in accelerator structure:'
                  write(*,*)' for SINGLE ELEMENT ',bez(ii)
                  write(*,*)' smiv at first occurrence: ',
     &                                                oriSmivSEDynks(kk)
                  write(*,*)' smiv of ',jj,' entry: ',smiv(1,jj)
                  write(*,*)''
                  lerr= .true.
                endif
              endif
            endif
          endif

(...)
 
!     go to next flagged SINGLE ELEMENT
      enddo
 
!     at least one inconsistency that should be solved by the user:
      if ( lerr ) then
        write(*,*) ''
        write(*,*) ' at least one incosistency in flagging elements'
        write(*,*) '    for dynamic kicks: please check carefully...'
        write(*,*) ''
        call prror(-1)
      endif
```

* `smiv(1,jj)` is the _native_ variable which stores the kick of thin lens.

* `oriSmivSEDynks(kk)` is a variable introduced by Alessio, which simply verifies that every instance of a SINGLE ELEMENT in the accelerator lattice is applied the same kick.

In our case we have `bez(ii)` = CRAB5 , `oriSmivSEDynks(kk)` =  541.17124490936453 , `jj`= 168, `smiv(1,jj)` =  451.41929755283127.

We have the error because 451.41929755283127 not equal to 541.17124490936453

The original code doesn't take into account a kick applied to a nominal kick _with_ multipole terms or random fluctuations. 

What the code does is applying a kick with the same parameters to the elements of the same family, but because on each element of each family it should at least be different, he just propagates one of the kicks to the elements of the same family.

What we need is that each element has it's own kick, applying the new kick to the original values of the variables which store the multipole terms and the random fluctuations.

#### Random Fluctuations

Looking more closely at the assignment of `smiv` in `sixve.f`:

```
          smiv(m,i) = sm(ix) + smizf(i)
          smi(i)    = smiv(m,i)
```

Where `sm(i)` is the first additional datum in element declaration for a _non-linear_ element, interpreted as average multipole strength. Nominal kick.

`smizf(i)` stores the random fluctuation:

```
smizf(i) = zfz(izu)*ek(ix)
```

Where `ek(i)` is the second additional datum in the element declaration, and  `zfz(izu)` contains a randomly sampled number (comes from `rveco = rvec(i)`).

__We have to store `smiv` for each element__.

### Multipole Coefficients

The multipole coefficients are in block 740.

* `im` : counter
* `benkc` : bending strength
* `irm(x)` : family
* `il`: number of single element
* `i` : multipole order
* `m`: always equal to 1
* `nmu(ix)` : maximum multipole order assigned to a certain element (only number 11 elements)
* `aka, bka`: central multipole coefficients. Linopt routine reassigns `aka`.
* `ak0, bk0`: RMS multipole coefficients
* `aaiv, bbiv`: final value of the multipole for a given family, used in the tracking block. Used in _track.f_, _trauthin_ subroutine (`aai`).

__We have to store `aaiv, bbiv` for each element__.

__No propagating__.

### Notes on SixTrack compilation flow

* Astuce files extract a FORTRAN version from the source code (sixtrack.s ---> sixven.f).
* The FORTRAN compiler translates the source code into object code (sixtrack.sh --> track.o, sixve.o).

SixTrack is compiled by `make_six`. This modifies the astuce files by turning on or off certain flags and decks, and calls astuce which generates the fortran files. This is then compiled and linked by the selected compiler.

### daten deck
* __Flag used__: bnlelens, crlibm, hhp, collimat, vvector, iibm, time, debug, cr, llyap, fio.

* __Used in Ast Files__: sixve, sixsc.

```
!  USED DISKS:
!
!  GEOMETRY AND STRENGTH OF THE ACCELERATOR : UNIT  2
!  TRACKING PARAMETER                       : UNIT  3
!  NORMAL PRINTOUT                          : UNIT  6
!  TRACKING DATA                            : UNIT  8
!  DATA FOR SUMMARY OF THE POSTPROCESSING   : UNIT 10
!  AUXILIARY FILE FOR THE INPUT             : UNIT 11
!  ASCII FILE WITH THE HORIZONTAL FFT DATA  : UNIT 14
!  ASCII FILE WITH THE VERTICAL FFT DATA    : UNIT 15
!  METAFILE FOR PLOTTING WITH GKS           : UNIT 20
!
!  CHECKPOINT/RESTART FILES                 : UNIT 95,96
!  OPTIONAL DUMP.DEBUG FILE                 : UNIT 99
!  PROGRESS FILE                            : UNIT 91
!  INTERMEDIATE OUTPUT FILE (LOUT)          : UNIT 92
!  CHECKPOINT/RESTART LOGFILE               : UNIT 93
!  TEMPORARY SCRATCH FILE for C/R           : UNIT 94
```


#### Single Elements Block

```
120 i=1
```
`i` identifies a given single element in the `SING` declaration part.



```
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
```

* `kz(i)` second datum in element declaration, interpreted as __type__. (23 = Crab Cavity)
* `ed(i)` first additional datum in element declaration.
* `ek(i)` second additional datum in element declaration.
* `el(i)` third additional datum in element declaration.
* `pieni`= 1d-38


if(ix.eq.is(1).or.iratioe(ix).eq.is(1)) then
              smiv(m,ncrr)=smi(ncrr)
            else if(ix.eq.is(2).or.iratioe(ix).eq.is(2)) then
              smiv(m,ncrr)=smi(ncrr)

subroutine applydynks(n)
subroutine saveorigsmiv

!       copy it to all occurrences of the current SINGLE ELEMENT
        do jj=1,iu
           if ( ktrack(jj).ne.1 ) then
!             SINGLE ELEMENT, not a BLOCK
              if ( ic(jj)-nblo.eq.ii ) then
!                current entry is another instance of the selected
!                   SINGLE ELEMENT
                 smiv(1,jj) = tmpsmiv
!                same code as in ripple:
                 strack(jj)=smiv(1,jj) *powerscale( ktrack(jj) )
                 strackc(jj)=strack(jj)*tiltc(jj)
                 stracks(jj)=strack(jj)*tilts(jj)
              endif
           endif
!       go to next entry in accelerator structure
        enddo

!     go to next flagged SINGLE ELEMENT
      enddo



      !     loop over all the SINGLE ELEMENTs flagged for dynamic kicks:
      do kk=1,NacqDynkSEs
        ii=iSEDynks(kk)



      i  = identifies a given entry in the sequence (a BLOC of linear elements or single non linear elements)
      j  = identifies a given particle

      Entry i is identified as the single element jj

      jj = current single ELEMENT
      iu = number of entries in the accelerator
      il = number of single ELEMENT
      ktrack = index in GOTO statements, based on the type of current entry
      ic = numerical identifier of the current entry (ix in do loops)
      nblo = number of BLOC max allowed

aaiv(mmul,nmac,nblz)

7203: !hr05         aaiv(k,m,i)=ed(ix)*(ak0(im,k)+zfz(izu)*aka(im,k))/r0a
 7204:               aaiv(k,m,i)=(ed(ix)*(ak0(im,k)+zfz(izu)*aka(im,k)))/r0a    !hr05
 7205:               aai(i,k)=aaiv(k,m,i)
 7206                izu=izu+1
 7207  !hr05         bbiv(k,m,i)=ed(ix)*(bk0(im,k)+zfz(izu)*bka(im,k))/r0a

 do 750 j=1,il
      if(imn.eq.bez(j)) then
        irm(j)=im
        goto 760
      endif