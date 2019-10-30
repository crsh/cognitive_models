!********************************************
! Tools for the MinervaAL project:
!
!	subroutine Get_echo (C, P, M, n, r)
!	subroutine Get_stimuli (oMat, n)
!	subroutine Randomize_order (oVec, n)
!	function Similarity (v1, v2, n)
!	function Cosine (v1, v2, n)
!	function Mean (iVec, n)
!	function SEM (iVec, n)
!
!********************************************
module MinervaAL_tools
use Number_generators
implicit none


CONTAINS

!--------------------------
subroutine Get_echo (C, P, M, n, r)
implicit none
  integer		:: j, n 
  real			:: C(:), P(:), M(:,:), r
	
	C = 0.0
	
	do j = 1, 200
		C(j) = flat(r) * binomial(0.5)
	enddo
		
	do j = 1, n
		C(:) = C(:) + Cosine(P(1:180), M(1:180,j), 180)**3 * M(:,j)		
	enddo

	C = C/maxval(abs(C))
		
 return
end subroutine Get_echo

!--------------------------
subroutine Get_stimuli (oMat, n)
implicit none
  integer	:: n, j, v,k
  real		:: oMat(:,:)
  
	oMat = 0.0
	
	v = 0
		do j = 1, n
			do k = 1+v, 20+v
				oMat(k,j) = 1.0
			enddo
		  v = v + 20
		enddo
  
 return
end subroutine Get_stimuli

!-----------------------------
subroutine Randomize_order (oVec, n)
implicit none
 integer		:: i, j, oVec(:), n
 logical		:: b
  oVec = 0
	do i = 1, n 
		do
		  oVec(i) = FlatInt(n) 
		    b = .TRUE.
			  do j = 1, n 
				if (oVec(i) == oVec(j) .and. i /= j) b = .FALSE.
			  enddo
		 if (b) exit
		enddo
	enddo
 return
end subroutine Randomize_order

!------------------------------
function Cosine (v1, v2, n)
 implicit none
	integer		:: n
	real		:: v1(n), v2(n), x(n), y(n), Cosine
	
	Cosine = 0.0
	
	if (sum(abs(v1)) /= 0.0 .and. sum(abs(v2)) /= 0.0) then

	x = v1/sqrt(dot_product(v1,v1))
	y = v2/sqrt(dot_product(v2,v2))
	
	Cosine = dot_product(x,y)
	
	endif
	
 return
end function Cosine

!------------------------------
function Mean (iVec, n)
implicit none
	integer		:: n
	real		:: iVec(n), Mean
	
	Mean = Sum(iVec)/n
	  
  return 
 end function Mean
  
!------------------------------  
function SEM (iVec, n)
 implicit none
 	integer		:: i, n
	real		:: iVec(n), SEM, M, Summ
  
  	M = Mean(iVec, n)
	
	Summ = 0.0
		do i = 1, n
			SUMM = SUMM + (iVec(i) - M)**2
		enddo
	
	SEM = sqrt(Summ/(n-1)) / sqrt(real(n))
  
 return
end function SEM

END MODULE MinervaAL_tools
