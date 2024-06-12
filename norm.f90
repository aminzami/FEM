subroutine norm(A,n,out)
real :: A(n),out,sum1
integer ::ii,n

sum1=0
do ii=1,n
sum1=sum1+A(ii)**2
end do

out=sqrt(sum1)

end subroutine