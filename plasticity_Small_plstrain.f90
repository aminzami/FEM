   subroutine plastic_computation(sigma,ep_n,delta_ep,ep_n33_conv,ep_n33,beta_n,beta_33,a_n,norm_ep,C_elasplas,a_nAll,delta_gamma,gamma,strain_all)

    implicit none

    ! Variables


real ::temp,temp_norm,delta_gamma,Tol,Dg,d_gamma2,gamma,norm_ep ,&
             ep_33np1,ep_n33,ep_n33_conv,beta_33,B, temp_gamma,p_strain33
real ,dimension(3) ::sigma,s_np1_trial,beta_n,beta_np1,&
                     zeta_np1_trial,ep_n,strain_all,delta_ep,p_strain,D_strain,S_np1_trial2
real ::n_np1(1,3),n_np1_2dar(1,3)
integer ii,jj,method,nn
integer io2(1,3)
integer,dimension(3,3):: io4, idyad
real :: mu=7.167e10,lambda=9.69e10,young=2e11,nu=0.3,sigma_y=2.5e8, &
       G,k, & !shear and bulk modulus 
       k_n,a_n,a_nAll , H_prim =6e8, K_prim=6e8,zeta33,sigma33
 real :: f_np1_trial,g1,g2 , theta_np1,bar_theta_np1,n3_np1,S33 ,Hyd   !H_np1,H_n=0
 real,dimension(3,3) :: C_elasplas,C_elasplas2,C_elas,ndyad


        C_elas(1,1)=1-nu
        C_elas(2,2)=1-nu
        C_elas(3,3)=(1-2*nu)/2
        C_elas(1,2)=nu
        C_elas(2,1)=nu
        C_elasplas2=(K_prim/(K_prim+young))*(young/((1+nu)*(1-2*nu)) )*C_elas
        C_elas=(young/((1+nu)*(1-2*nu)) )*C_elas

io2(1,1:3)=(/1,1,0/)
k=lambda+(2./3)*mu;

io4(1,1)=1
io4(2,2)=1
io4(3,3)=0.5   !!!! avaz shode az 0.5



 p_strain=ep_n !+delta_ep
 p_strain33=ep_n33_conv +ep_n33

 a_n=a_nAll
        sigma33=(nu)*(sigma(1)+sigma(2)) -2*mu*0*(1+nu)*p_strain33  !!! avaz shode 0 
        Hyd=(1./3)*(sigma(1)+sigma(2)+sigma33)
        S_np1_trial(1:2)=sigma(1:2)-Hyd*(/1,1/)
        S_np1_trial(3)=sigma(3)
        S33=-(S_np1_trial(1)+S_np1_trial(2));  ! *********<<  NOTICE >>> *******

       zeta_np1_trial=S_np1_trial -beta_n;
       zeta33=S33-beta_33
  
       K_n=K_prim*a_n;  !!! inja a_n bayad ba termaye gradient jam she..
       temp_norm=0;

      temp_norm=sqrt(zeta_np1_trial(1)**2 +zeta_np1_trial(2)**2 + 2*(zeta_np1_trial(3)**2) + zeta33**2 ) ! *********  sqrt (3./2) ehtemalan mikhad <<  NOTICE >>> *******
      f_np1_trial=temp_norm -sqrt(2./3)*(sigma_y +K_n); !!! sathe taslim bayad non-local beshe..

      if (f_np1_trial >0) then
      n_np1(1,:)=zeta_np1_trial(:) /temp_norm;

      !!!*********************
     n_np1(1,3)=n_np1(1,3)/2  !!! comment dar amade
      !!!*********************

      n_np1_2dar=n_np1
      n_np1_2dar(1,3)=2*n_np1_2dar(1,3)
      n3_np1=zeta33/temp_norm

    delta_gamma=0;     !!!! ino avaz mikonam
       d_gamma2=1e5  !! test value for enter the loop
       Tol=1e-8;
       g1=1;
       H_prim=0   !!! kinematic Hardening  =0



      do while (abs(delta_gamma-d_gamma2) >tol)

        d_gamma2=delta_gamma;
        a_n=a_n +sqrt(2./3)*delta_gamma;    !!! inja sqrt (2/3) dar a_n hast dar sathe taslim ham hast
        K_n=K_prim *a_n;   !!! sathe taslim bayad non-local beshe..
        g1=-sqrt(2./3)*(K_n+sigma_y) + temp_norm -(2*mu*delta_gamma + sqrt(2./3)*(H_prim*delta_gamma));
        Dg=-(2*mu + sqrt(2./3)*(H_prim) +(2./3)*K_prim)   !!! moshtaghe sathe taslim avaz mishe dar non-local ...
        delta_gamma=delta_gamma -(g1/Dg);
       end do
     
       temp_gamma=gamma+delta_gamma;
       

!!       %Back stress

!    beta_n=beta_n + sqrt(2./3)*(H_prim*delta_gamma)*n_np1(1,:);
!    beta_33=beta_33 + sqrt(2./3)*(H_prim*delta_gamma)*n3_np1;

    delta_ep= delta_gamma*n_np1(1,:);   !! ep_n dar vaghe sefr dar nazar gerefte shode
    ep_n33=  delta_gamma*n3_np1
    norm_ep=norm_ep + sqrt( delta_ep(1)**2 + delta_ep(2)**2 + 2*(delta_ep(3)**2) + ep_n33**2 )
!    a_nAll=a_nAll+a_n

    sigma=sigma-2*mu*delta_gamma *n_np1(1,:);


!    sigma=sigma-(lambda)*ep_n33*(/1,1,0/)-2*mu*delta_gamma *n_np1(1,:);





!     p_strain33=ep_n33_conv +ep_n33

!!!******************************************************  comment shode
!
!        sigma33=(nu)*(sigma(1)+sigma(2))  -2*mu*(1+nu)*p_strain33   !!!! avaz shode  0 
!        Hyd=(1./3)*(sigma(1)+sigma(2)+sigma33)
!        S_np1_trial(1:2)=sigma(1:2)-Hyd*(/1,1/)
!        S_np1_trial(3)=sigma(3)
!
!        
!       S33=-(S_np1_trial(1)+S_np1_trial(2));  ! *********<<  NOTICE >>> *******
!
!       zeta_np1_trial=S_np1_trial -beta_n;
!       zeta33=S33-beta_33
!  
!       K_n=K_prim*a_n;
!       temp_norm=0;
!
!      temp_norm=sqrt(zeta_np1_trial(1)**2 +zeta_np1_trial(2)**2 + 2*(zeta_np1_trial(3)**2) + zeta33**2 ) ! *********  sqrt (3./2) ehtemalan mikhad <<  NOTICE >>> *******
!     
!      n_np1(1,:)=zeta_np1_trial(:) /temp_norm;
!
!      n3_np1=zeta33/temp_norm

!!!******************************************************  comment shode








!    sigma33=sigma33-2*mu*delta_gamma *n3_np1

   
!    
!!    Consistent tangent moduli
   theta_np1=4.0*mu*mu* (delta_gamma/temp_norm);

    bar_theta_np1=2.0*mu*(1.0/(1.0+(K_prim+H_prim)/(3.0*mu))) 


    idyad=matmul(transpose(io2),io2)


     ndyad=matmul(transpose(n_np1),n_np1)
!    C_elasplas=k* idyad + 2*mu*theta_np1*(IO4-(1./3)*idyad) -2*mu*bar_theta_np1*ndyad;

!     C_elasplas=C_elas-(2.0*mu*ndyad/(1.0+(K_prim+H_prim)/(3.0*mu)))  !!!! javab dade behtar bude , max_iteration =12 ,Tol 1e-4

    C_elasplas=C_elas -((bar_theta_np1-theta_np1)*ndyad ) -(4.0*mu*mu*delta_gamma/temp_norm)*Io4 + (4.0*mu*mu*delta_gamma/(3.0*temp_norm))*idyad






      endif



!!!!!!!!!!*********************************************************************

!write (2,'(16f14.8)') delta_gamma
!write (3,'(16f14.8)') gamma


     !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



    end subroutine plastic_computation

