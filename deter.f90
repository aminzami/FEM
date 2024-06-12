subroutine det(A,out)

real::A(2,2),out

out=A(2,2)*A(1,1)-A(1,2)*A(2,1)

end subroutine