module Number_generators

!*********************************************************************
!User accessible routines:
!Subroutines:       
!   RandSeed   sets the random-number seed from the wall clock
!   FixSeed    sets the random-number seed to -99
!              both routines write the seed to standard output
!
!Functions:
!   flat(range)      returns a real from a flat distribution (0.0, .. range]
!                    range (real) sets the max value
!
!   FlatInt(range)   returns an integer from a flat distribution(0, 1, .. range]
!                    range (integer) sets the upper value
!
!   gaussian(mu, sd) returns a real from a Normal(mu, sd)
!                    mu and sd (real) set the mean & SD of the Normal
!
!   goemetric(p)     returns an integer (1..maxint] from a geometric  distribution
!                    p (real) is the probability of success on each trial
!
!   binomial(p)      returns +1 / 0 with p = p(+1)
!
!   exponential(lambda) returns a real from an exponential distribution
!                    lambda is the mean of the distribution
!
!***********************************************************************
! Internal routines:  
!  the random-number seed, idum, is private to the module
!  ran3 (idum):  returns a random number from a uniform distribution (0..1]
!  rnorm(idum):  returns a random number from a Unit Normal distribution
!  expdev(idum): returns a random number from an exponential with mean = 1
!***********************************************************************

 implicit none
 integer, private             :: idum, inext, inextp, inext1, inextp1
 real, private, dimension(55) :: ma, ma2
 logical                      :: switch

 CONTAINS


!***********************************************************************
!  ran3:  generates random values 0..1] 
!  seed (idum) is set to a negative integer on entry
!  From Press et al.,(1992) Numerical recipes in FORTRAN (2nd ed.) CUP 
!
!  This routine runs through the ma vector pre-established when the
! generator is initialized
!***********************************************************************
function ran3(idum)
 implicit none
 real                  :: ran3, mj, mk
 real, parameter       :: MBIG =4000000.0, MSEED =1618033.0, MZ =0.0, FAC = 1.0/MBIG
 integer               :: idum
  
 inext  = inext + 1
 if(inext == 56) inext = 1
 inextp = inextp + 1
 if(inextp == 56)inextp= 1
 
 mj = ma(inext) - ma(inextp)
   
 if(mj < MZ) mj = mj + MBIG
 ma(inext) = mj
 ran3 = mj*FAC
 return
end function ran3




!****************************************************************
! fx_geometric:  returns p(1st success) at trial first_success
!    parameter:  p_success = probability of a success on all trials
!****************************************************************
function fx_geometric(p_success, first_success)
 real, intent (in)      :: p_success
 real                   :: fx_geometric
 integer, intent(in)    :: first_success

 if (first_success < 1) then  ! Can't have a negative number--return
    fx_geometric = 0.0        !impossible result, i.e., Prob = zero
 else
    fx_geometric = ((1.0 - p_success)**(first_success-1)) * p_success
 endif
 return  
end function fx_geometric


!*****************************************************************
!  expdev:  returns a real from an exponential distribution lambda=1
!*****************************************************************
function expdev(idum)
 real                   :: expdev, dum
 integer, intent(inout) :: idum
 
10  dum = ran3(idum) 
      if(dum == 0.0) goto 10
    expdev = -log(dum)
    return
end function expdev          



!******************************************************************
! function exponential(lambda)  returns a deviate from an 
!    exponential distribution with mean = lambda
!******************************************************************
function exponential (lambda)
real, intent(in)     :: lambda
real                 :: exponential

 exponential = expdev(idum)*lambda
 return
end function exponential

 


!******************************************************************
!  geometric: returns an integer 1 ... inf. The values are 
!                 distributed according to a geometric distribution
!  parameter: p_success = probability of a success on each trial 
!******************************************************************    
function geometric(p_success)
 integer           :: j, geometric
 real, intent(in)  :: p_success
 real              :: prob, xs, rn

 j  = 0                              
 xs = 0.0                            ! start cummulative at zero
 rn = ran3(idum)                     ! get random probability value

10 j = j + 1                         ! Search loop: searching the cummulative
   xs = xs + fx_geometric(p_success, j) ! probability of a success on trial j             
 if (xs .le. rn) goto 10                ! Search until the cumulative exceeds
                                        ! rn.  Return the number of failures
 geometric  = j                         ! before the 1st success, i.e.,
end function geometric                  ! where rn falls in the cumulative




!*********************************************************************
! rnorm:  Unit Normal distribution
!*********************************************************************
function rnorm (idum)
 logical, save           :: switch
 data switch /.true./
 real                    :: fac, r, v1, v2, rnorm
 integer, intent(inout)  :: idum
 real,    save           :: rnorm2
 
 if (switch) then
10  v1 = 2.0 * ran3(idum) - 1.0
    v2 = 2.0 * ran3(idum) - 1.0
    r = v1**2 + v2**2
    if ((r .ge. 1.0).or.(r .eq.0)) goto 10

    fac = sqrt(-2.0 * log(r)/r)
    rnorm2 = v1 * fac
    switch = .false.
    rnorm  = v2 * fac
 else
    switch = .true.
    rnorm  = rnorm2
 endif
 return
end function rnorm



!***********************************************************************
! gaussian:  returns a Normal deviate from N(mu, sd)
!***********************************************************************
function gaussian(mu, sd)
 real, intent(in)        :: mu, sd
 real                    :: gaussian
 
 gaussian = (rnorm(idum) * sd) + mu   ! Calculate Normal(mu, sigma)
 return   
end function gaussian



!***********************************************************************
! binomial:  returns a +1 / -1 from a binomial p = (probability of +1)
!***********************************************************************
function binomial(p) 
 real, intent(in)      :: p
 real                  :: binomial
 if (ran3(idum) .lt. p) then
    binomial = +1.0
 else
    binomial = -1.0
 endif      
end function binomial



!**********************************************************************
! flat :  returns a real from a flat probability distribution
! parameter:  range
!**********************************************************************
function flat(range)
 real             :: flat
 real, intent(in) :: range
 
 flat = ran3(idum)*range
 return
end function flat


!**********************************************************************
! FlatInt :  returns an integer from a flat probability distribution
! parameter:  range
!**********************************************************************
function FlatInt(range)
 integer              :: FlatInt
 integer, intent(in)  :: range
 
 FlatInt = int(ran3(idum) * range) + 1
 return
end function FlatInt 


!**********************************************************************
! FixSeed:  Sets a fixed seed (-99) for the random-number routine
! and sets up the ma vector for the random number generator
!**********************************************************************
subroutine FixSeed
 real                  :: mj, mk
 real, parameter       :: MBIG =4000000.0, MSEED =1618033.0, MZ =0.0 
 integer               :: i, ii, k
 idum = -99
 write(*,'(A, I0)') ' Seed for random-number generator = ', idum
 
      mj= MSEED - iabs(idum)
      mj = mod(mj, MBIG)
      ma(55) = mj
      mk=1
      do i=1,54
         ii = mod(21*i, 55)
         ma(ii) = mk
         mk = mj - mk
         if(mk < MZ) mk = mk + MBIG
         mj = ma(ii)
      enddo
      do k = 1, 4
         do i = 1, 55
            ma(i) = ma(i) - ma(1 + mod(i+30, 55) )
            if(ma(i) .lt. MZ) ma(i) = ma(i) + MBIG
         enddo
      enddo
      inext = 0
      inextp = 31
 return
end subroutine FixSeed



!**********************************************************************
!  RandSeed:  Sets a random seed for the random-number routines
!  and sets up the ma vector for the random number generator
!**********************************************************************
subroutine RandSeed 
 implicit none
 real                  :: mj, mk
 real, parameter       :: MBIG =4000000.0, MSEED =1618033.0, MZ =0.0 
 integer               :: i, ii, k
 call system_clock(idum)
 
10 if (idum >  1000000) idum = idum / 10 
   if (idum > 1000000) go to 10
   
   if (idum > 0) idum = idum *(-1)
   
   write(*,'(A, I0)') ' Seed for random-number generator = ', idum

        mj= MSEED - iabs(idum)
        mj = mod(mj, MBIG)
        ma(55) = mj
        mk=1
        do i=1,54
          ii = mod(21*i, 55)
          ma(ii) = mk
          mk = mj - mk
          if(mk < MZ) mk = mk + MBIG
          mj = ma(ii)
        enddo
        do k = 1, 4
          do i = 1, 55
             ma(i) = ma(i) - ma(1 + mod(i+30, 55) )
             if(ma(i) .lt. MZ) ma(i) = ma(i) + MBIG
          enddo
        enddo
        inext = 0
        inextp = 31
   
   return
end subroutine RandSeed


!**********************************************************************
! get_ran_seed:  gets current state of the random-number generator
!**********************************************************************
subroutine get_ran_seed(dummy, MAA)
  integer, intent(out), dimension(2) :: dummy
  real, dimension(55)                :: MAA
    dummy(1) = inext
    dummy(2) = inextp 
    MAA = ma
  return
end subroutine get_ran_seed


!**********************************************************************
! assign_seed:  restores the state of the random-number generator
!**********************************************************************
subroutine assign_seed(iseed, MAA)
  integer, dimension(2), intent(in)  :: iseed
  integer                            :: get_ran_seed
  real, dimension(55)                :: MAA
    inext  = iseed(1)
    inextp = iseed(2) 
    ma     = MAA
    switch = .true.
  return 
end subroutine assign_seed



END MODULE Number_Generators



