SUBROUTINE Gauss (a,inv_a,n)       ! Invert matrix by Gauss method
! --------------------------------------------------------------------
IMPLICIT NONE

INTEGER :: n
REAL :: a(n,n),inv_a(n,n)


! - - - Local Variables - - -
REAL :: b(n,n), c, d, temp(n)
INTEGER :: i, j, k, m, imax(1), ipvt(n)
! - - - - - - - - - - - - - -

b = a
ipvt = (/ (i, i = 1, n) /)

DO k = 1,n
   imax = MAXLOC(ABS(b(k:n,k)))
   m = k-1+imax(1)

   IF (m /= k) THEN
      ipvt( (/m,k/) ) = ipvt( (/k,m/) )
      b((/m,k/),:) = b((/k,m/),:)
   END IF
   d = 1/b(k,k)

   temp = b(:,k)
   DO j = 1, n
      c = b(k,j)*d
      b(:,j) = b(:,j)-temp*c
      b(k,j) = c
   END DO
   b(:,k) = temp*(-d)
   b(k,k) = d
END DO

inv_a(:,ipvt) = b

END SUBROUTINE Gauss