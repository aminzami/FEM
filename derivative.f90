subroutine derive(kk,x,y,d_var1_var2)

 implicit none
 real,dimension(4) :: x,y,Xquad_p,Yquad_p,zeta_i,eta_i
 real ::d_var1_var2(2,2)
 integer::ii,kk


! Xquad_p=0.577*(/1,1,-1,-1/)
! Yquad_p=0.577*(/-1,1,1,-1/)
! zeta_i=(/1,1,-1,-1/)
! eta_i=(/-1,1,1,-1/)

!   Xquad_p=0.577*(/1,-1,-1,1/)
! Yquad_p=0.577*(/1,1,-1,-1/)
! zeta_i=(/1,-1,-1,1/)
! eta_i=(/1,1,-1,-1/)

  Xquad_p=sqrt(1./3)*(/-1,1,1,-1/)
  Yquad_p=sqrt(1./3)*(/-1,-1,1,1/)
  zeta_i=(/-1,1,1,-1/)
  eta_i=(/-1,-1,1,1/)

!! ye bar x koochak miad ye bar X bozorg avali dar makoose dovomi zarb mishe
d_var1_var2=0


 do ii=1,4
 ! for each quadrature
 
 ! weight function is 1


  d_var1_var2(1,1)=d_var1_var2(1,1)+(1./4)*x(ii) *zeta_i(ii)*(1+ eta_i(ii) *Yquad_p(kk))
!   d_var1_var2(1,2)= d_var1_var2(1,2)+(1./4)*yi(ii) *zeta_i(ii)*(1+ eta_i(ii) *Yquad_p(kk))
!     d_var1_var2(2,1)=d_var1_var(2,1)+(1./4)*xi(ii) *eta_i(ii)*(1+ zeta_i(ii) *Xquad_p(kk))
   d_var1_var2(1,2)= d_var1_var2(1,2)+(1./4)*x(ii) *eta_i(ii)*(1+ zeta_i(ii) *Xquad_p(kk))
     d_var1_var2(2,1)=d_var1_var2(2,1)+(1./4)*y(ii) *zeta_i(ii)*(1+ eta_i(ii) *Yquad_p(kk))
       d_var1_var2(2,2)=d_var1_var2(2,2)+(1./4)*y(ii) *eta_i(ii)*(1+ zeta_i(ii) *Xquad_p(kk))



 
  end do



end subroutine derive