    subroutine interp_profile(nup, ndown, ne_act, t_act, profile, alpha)
      

      implicit none
      
      integer pmne, pmp, pmt, mline
      parameter (pmne = 17)       ! max # of n_e
      parameter (pmt = 7)         ! max # of T
      parameter (pmp = 65)        ! max # of profile points
      parameter (mline = 21)      ! max # of H lines
      
      
      real svcs( pmt, pmne, 0:pmp, mline ) ! VCS profiles
      
      real log_alpha0(mline),   ! log Delta_alpha_min
     +     log_ne0(mline),      ! log n_e_min
     +     log_t0(mline)        ! log T_min
      
      real log_alpha_inc(mline), ! Delta log Delta_alpha
     +     log_t_inc(mline),    ! Delta log n_e
     +     log_ne_inc(mline)    ! Delta log T
      
      integer nl(mline), nu(mline), mp(mline), mne(mline), mt(mline)
      character*1 null(pmt)
      
      integer i, j, k, line, nline, nnl, nnu
      integer ihigh,jhigh,ilow,jlow
      integer npts
      integer ndown,nup
      real ne_act, t_act, ne_high, ne_low, t_high, t_low
      real profile_low,profile_high
      real profile(pmp)
      real ne, f0
      real ireal,jreal
      real xmin,xmax,ymin,ymax,alpha(pmp)
      logical ok
      
      real exp10, x
      exp10( x ) = exp( 2.30258509299405E0*x )
      
C     
C     OPEN DATA FILES
C     
 299  continue
      if (ndown.eq.1) then
         open(unit=15,
     +        file='/Users/sl/code/fortran/cv_code/opacities/lymannd'
     +        ,status='old')
      else if (ndown.eq.2) then
         open(unit=15,
     +        file='/Users/sl/code/fortran/cv_code/opacities/balmernd'
     +        ,status='old')
      else if (ndown.eq.3) then
         open(unit=15,
     +        file='/Users/sl/code/fortran/cv_code/opacities/paschnd'
     +        ,status='old')
      else if (ndown.eq.4) then
         open(unit=15,
     +        file='/Users/sl/code/fortran/cv_code/opacities/bracknd'
     +        ,status='old')
      else 
         write(*,'(a,$)') 
     +        'Lower level invalid - please re-enter: '         
         read(*,*) ndown
         goto 299
      end if
c      open(unit=99,file='output.dat',status='unknown')
c      open(unit=16,file='errors.dat',status='unknown')

c     
C     READ IN VCS ARRAYS
c     
        read( 15, * ) nline
        if( nline .gt. mline ) then
c           write( 99, * ) 'Table too big.  Not more than ', mline,
c     +             ' lines.'
           call exit(1)
        end if
        do i = 1, nline
           read( 15, * ) nl(i), nu(i),
     +                  log_alpha0(i), log_ne0(i), log_t0(i),
     +                  log_alpha_inc(i), log_ne_inc(i), log_t_inc(i),
     +                  mp(i), mne(i), mt(i)
           if( mp(i)  .gt. pmp .or.
     +         mne(i) .gt. pmne .or.
     +         mt(i)  .gt. pmt ) then
c              write( 99, * ) 'Table too big in one of these:'
c              write( 99, * ) 'mp:', mp(i), ' >', pmp
c              write( 99, * ) 'mne:', mne(i), ' >', pmne
c              write( 99, * ) 'mt:', mt(i), ' >', pmt
              call exit(1)
           end if
        end do
        do line = 1, nline
           read( 15, '(3x,i2,4x,i2)' ) nnl, nnu
           if( nnl .ne. nl(line) .or. nnu .ne. nu(line) ) then
c              write( 16, * ) 'Inconsistency in table for', nl(line),
c     +                      ' ->', nu(line)
              call exit(1)
           end if
           read( 15, * ) (( (svcs(i,j,k,line), k = 0, mp(line)),
     +                        i = 1, mt(line)),
     +                          j = 1, mne(line))
        end do
        close(15)
C
C     Select bounding profiles
C
 199    continue
        do line=1,nline
        if (nl(line).eq.ndown .and.nu(line).eq.nup) then
           ok = .true.
           ireal = 1.+(log10(t_act) - log_t0(line))/log_t_inc(line)
           jreal = 1.+(log10(ne_act) - log_ne0(line))/log_ne_inc(line)
           ilow = int(ireal)
           ihigh = ilow+1
           jlow = int(jreal)
           jhigh = jlow+1
           ne_low  = exp10(log_ne0(line) + log_ne_inc(line)*(jlow-1))
           ne_high = exp10(log_ne0(line) + log_ne_inc(line)*(jhigh-1))
           t_low   = exp10(log_t0(line) + log_t_inc(line)*(ilow-1))
           t_high  = exp10(log_t0(line) + log_t_inc(line)*(ihigh-1))
c
c     Check input parameters
c

           if (t_high.gt.exp10(log_t0(line)+log_t_inc(line)*
     +          (mt(line)-1))) then
              write(*,*) 'WARNING:'
              write(*,'(a31,i6)') 'Temp greater than max value of ',
     +             NINT(exp10(log_t0(line)+log_t_inc(line)
     +             *(mt(line)-1)))
              write(*,'(A,$)') 'Please re-enter temp: '
              read(*,*) t_act
              goto 199
           end if

           if (t_low.lt.exp10(log_t0(line))) then
              write(*,*) 'WARNING:'
              write(*,'(a28,i6)') 'Temp less than min value of ',
     +             NINT(exp10(log_t0(line)))
              write(*,'(A,$)') 'Please re-enter temp: '
              read(*,*) t_act
              goto 199
           end if

           if (ne_high.gt.exp10(log_ne0(line)+log_ne_inc(line)*
     +          (mne(line)-1))) then
              write(*,*) 'WARNING:'
              write(*,'(a30,e11.3)') 'Ne greater than max value of ',
     +             (exp10(log_ne0(line)+log_ne_inc(line)
     +             *(mne(line)-1)))
              write(*,'(A,$)') 'Please re-enter electron density: '
              read(*,*) ne_act
              goto 199
           end if

           if (ne_low.lt.exp10(log_ne0(line))) then
              write(*,*) 'WARNING:'
              write(*,'(a27,e11.3)') 'Ne less than min value of ',
     +             (exp10(log_ne0(line)))
              write(*,'(A,$)') 'Please re-enter elctron density: '
              read(*,*) ne_act
              goto 199
           end if

        end if
        end do

        if (.not.ok) then
           write(*,*) 'Upper level is invalid.'
           write(*,'(a,$)') 'Re-enter upper level: '
           read(*,*) nup
           goto 199
        end if
        
        do line=1,nline
           if(nl(line).eq.ndown .and. nu(line).eq.nup) then
           do k=1,mp(line)
C
C     interpolate in electron density (j)
C        
              profile_low=(svcs(ilow,jhigh,k,line) - 
     +             svcs(ilow,jlow,k,line))
     +             /(ne_high-ne_low)*(ne_act-ne_low) +
     +             svcs(ilow,jlow,k,line)
              profile_high=(svcs(ihigh,jhigh,k,line) - 
     +             svcs(ihigh,jlow,k,line))
     +             /(ne_high-ne_low)*(ne_act-ne_low) + 
     +             svcs(ihigh,jlow,k,line)
C     
C     interpolate in temperature (i)
C        
              profile(k) = (profile_high-profile_low)*(t_act-t_low)
     +             /(t_high-t_low) + profile_low
           end do
           end if
        end do

        xmin = 1.e35
        xmax = -1.e35
        ymin = 1.e35
        ymax = -1.e35
        do line = 1, nline 
           if (nl(line).eq.ndown .and. nu(line).eq.nup) then
              npts = mp(line)
              do k=1,mp(line)
                 x = log_alpha0(line) + log_alpha_inc(line)
     +                *(k - 1)
                 f0 = 1.25e-9*ne**(2./3.)
                 alpha(k) = x
                 if(x.lt.xmin) xmin = x
                 if(x.gt.xmax) xmax = x
                 if(profile(k).lt.ymin) ymin = profile(k)
                 if(profile(k).gt.ymax) ymax = profile(k)
              end do
           end if
        end do


    end