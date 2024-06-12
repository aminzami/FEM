!  FEM2_p1.f90 
!
!  FUNCTIONS:
!  FEM2_p1 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: FEM2_p1
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program FEM2_updated



      

    implicit none
    external sgetrf,sgetri
        ! Variables
  
!    real ::d0,sigma0,d_new,det,deter
    integer::ii,jj,n,kk,counter,ee,e,uu,mm,nn,temp,Count_Xmin=0,count_Xmax=0,ind_quad_S
    integer,ALLOCATABLE::element(:,:),numb(:),ipiv(:),bc_node(:),node_contribution(:,:),vicinity_ele(:,:)

!     integer ::  node_contributon(n,4)=-1;
!  integer  vicinity_ele(e,8)=-1
  integer:: vv,ff
    

 real,dimension(4) :: xi,yi,xit,yit,Xquad_p,Yquad_p,zeta_i,eta_i,ones,zeros
 real, dimension(2,2)::d_x_zeta,d_Xb_zeta,d_xs_zeta,inv_d_xb_zeta=0,inv_d_xs_zeta=0
 real::Ni_zeta(4,2),info,mu=7.167e10,division

 real :: temp_B(3,2)
 real::Bi(3,8)
 real::beta(4,2)

   real::f_int(8),u(8),F(2,2),f_ext(8)=0
   real::K_mat(8,8),k_geo(8,8),S_v(3),cauchy_t(2,2),temp_k_geo(4,4),F_inv(2,2),cauchy_v_temp(3)
   real,allocatable,dimension(:,:) ::inv_Ck_ass,Ck_ass,k_ass,CAUCHY,cauchy_v_all,WORK,cauchy_v_converge,vare_np1,&
                                     ep_n,beta_n,delta_ep_np1,strain_all,Strain_all_conv
   real,allocatable,dimension(:)::Xnode,Ynode,stress_vector,a_n,a_nAll,norm_ep,ep_n33,beta_33,delta_gamma,gamma_all, ep_n33_conv
   real ::Xmin,Xmax,Ymax,Ymin,Tol_inc
   real ::jacub,TAU(2,2),dete,dete1,norm_R_past


    real ::delta_w12,delta_D(3),delta_W(2,2),deltaW_dot_sigma(3)
     
   

   real ,dimension(:),allocatable::load,f_load,delta_u,R_ass,CR_ass,f_new,f_Reac,load_index !2*n
   real ::inc,dt, norm_delta_u_past !!!!!!!!!!!!!!!!
   real::norm_delta_u,norm_u_new,norm_R
   real ::U_local(8),unew_local(8)
   real ::stress_mises,S_dev(2,2),identity(2,2),S_dev3
   real,allocatable,dimension(:)::u_new,u_np1,bc,u_bc !,u_new_past
   real ::nu=0.3,young=2e11,C_elas(3,3)=0,lambda=9.69e10
   real,allocatable,dimension(:,:,:) :: C_elasplas


   
open(5,file=".\displacement_Fortran.plt",action='write',position ="new",Status="unknown")
open(6,file=".\stress_mises_Fortran.plt",action='write',position ="new",Status="unknown")


   !new quad 1bahman
!     n=533
!  e=480

   !!!!!!!!!!!!!!!!!!!harr
   n=341
   e=300
   !!!!!!!!!!!!!!!!!!! number of Nodes and elements

   allocate(C_elasplas(4*e,3,3))


   !tire nazok
! n=136
!   e=99

   allocate(load(2*n),f_load(2*n),delta_u(2*n),R_ass(2*n),CR_ass(2*n),f_new(2*n))
!   allocate(f_Reac(2*n))
   
   load=0
   f_load=0
   delta_u=0
   R_ass=0
   CR_ass=0
   f_new=0

  allocate( node_contribution(n,4),vicinity_ele(e,8))
   node_contribution=-1
   vicinity_ele=-1

   allocate (inv_Ck_ass(2*n,2*n),Ck_ass(2*n,2*n),k_ass(2*n,2*n))
   
!   allocate (invk_a(2*n,2*n))
   allocate (u_new(2*n))

    allocate (u_np1(2*n))
    allocate (element(e,4))
    allocate (Xnode(n))
    allocate (Ynode(n))
    allocate(stress_vector(n))
    allocate(numb(n))
    
   allocate(CAUCHY(e*2,4*2))
   allocate(cauchy_v_all(3,4*e))
   allocate(cauchy_v_converge(3,4*e))
   allocate(strain_all(3,4*e))          !!!
    allocate(strain_all_conv(3,4*e)) 
   allocate(vare_np1(3,4*e),delta_ep_np1(3,4*e),ep_n(3,4*e),beta_n(3,4*e),a_n(4*e),a_nAll(4*e),delta_gamma(4*e),gamma_all(4*e))
   allocate(ep_n33(4*e),beta_33(4*e))
   allocate( ep_n33_conv(4*e))



      allocate(WORK(2*n,2*n),ipiv(2*n))
      work=0
      ipiv=0
 
    cauchy_v_all=0
    cauchy_v_converge=0
    vare_np1=0
    ep_n=0
    beta_n=0
    a_n=0
    a_nAll=0
    ep_n33=0
    beta_33=0
    delta_gamma=0
    gamma_all=0
     ep_n33_conv=0
     strain_all=0
       strain_all_conv=0

        delta_ep_np1=0  
       ep_n33=0


    identity=reshape((/1,0,0,1/),(/2,2/))
    stress_vector=0

    CAUCHY=0
        C_elas(1,1)=1-nu
        C_elas(2,2)=1-nu
        C_elas(3,3)=(1-2*nu)/2
        C_elas(1,2)=nu
        C_elas(2,1)=nu
        C_elas=(young/((1+nu)*(1-2*nu)) )*C_elas

   u_new=0
   u_np1=0

   !!! doroste
!  Xquad_p=0.57735*(/1,-1,-1,1/)
!  Yquad_p=0.57735*(/1,1,-1,-1/)
!  zeta_i=(/1,-1,-1,1/)
!  eta_i=(/1,1,-1,-1/)

  Xquad_p=sqrt(1./3)*(/-1,1,1,-1/)
  Yquad_p=sqrt(1./3)*(/-1,-1,1,1/)
  zeta_i=(/-1,1,1,-1/)
  eta_i=(/-1,-1,1,1/)

  


 u=0

!  OPEN(UNIT=1,FILE='.\Job-longbeam.txt',STATUS='UNKNOWN')

  print *, "Xmin  ",  "Xmax ","        ","Ymin", "Ymax"
   print *, Xmin,Xmax,"        ",Ymin,Ymax
   print *, "norm_u_new"  , "   ", "norm_delta_u"

OPEN(UNIT=1,FILE='.\test-haar-quad.inp',STATUS='UNKNOWN')

!OPEN(UNIT=1,FILE='.\new-quad-1bahman.inp',STATUS='UNKNOWN')

   read (1,*) temp,Xnode(1),Ynode(1)



   Xmin=Xnode(1)
   Xmax=Xnode(1)  
   Ymax=Ynode(1)
   Ymin=Ynode(1)

   do ii=2,n
     read (1,*) temp,Xnode(ii),Ynode(ii)
     if (Xmin > Xnode(ii)) then
      Xmin=Xnode(ii)
      elseif (Xmax < Xnode(ii))  then
      Xmax=Xnode(ii)
      end if

       if (Ymin > Ynode(ii)) then
      Ymin=Ynode(ii)
      elseif (Ymax < Ynode(ii))  then
      Ymax=Ynode(ii)

     end if
    end do

!    Xnode=Xnode*1e3
!    Ynode=Ynode*1e3
!    Xmin=Xmin*1e3
!    Xmax= Xmax*1e3
!    Ymin=Ymin*1e3
!    Ymax=Ymax*1e3

 
   do ii=1,e
   read (1,*) temp,element(ii,1),element(ii,2),element(ii,3),element(ii,4)
   end do

   !*********** Vicinity_element *****************************************


!
!  
!   do nn=1,n
!    ff=1
!    do ee=1,e
!     do ii=1,4
!      if (nn .eq. element(ee,ii) ) then
!        node_contribution(nn,ff)=ee
!        ff=ff+1
!      end if
!     end do
!    end do
!   end do
!   !*************
! 
!
!   
!   do ee=1,e
!   vv=0
!    do jj=1,4
!     do ii=1,4
!     temp=node_contribution(element(ee,jj),ii)
!
!     if (temp .NE. ee  .And. temp .NE. -1) then
!       do kk=1,vv
!        if (temp .eq. vicinity_ele(ee,kk) ) then
!        goto 15
!        end if
!       end do
!       vv=vv+1
!      vicinity_ele(ee,vv)=temp
!      
!   15  end if
!     
!     end do
!    end do
!   end do
!
!
!    open(9,file=".\vicinity_ele.txt",action='write',Status="unknown")  !! for har_quad
!
!   do ee=1,e
!!   write (*,*) ,vicinity_ele(ee,:)
!    write(9,"(8I6.0)") ,vicinity_ele(ee,:)
!   end do
!
!
!   


   !*********** Vicinity_element *****************************************


   load=0
  !**********  Number of Xmin and Xmax
   do jj=1,n
    if (Xnode(jj) .eq. Xmax) then
       count_Xmax=count_Xmax+1
!        load_index=2*jj
!        load(load_index)=1 
    elseif (Xnode(jj) .eq. Xmin) then
        Count_Xmin=Count_Xmin+1
        endif
   end do 
   !***********


   allocate(load_index(count_Xmax))   !!!!!!!!!!!!!!!!!!!!!!!!
!allocate(load_index(1))              !!!!!!!!!!!!!!!!!!!!!!!!
   load_index=0
   mm=1

   !********* Apply load
      do jj=1,n
       if ((Xnode(jj) .eq. Xmax) .and. (Ynode(jj) .eq. Ymax)) then    !!!!!!!!!!!!!!!!!!!!!!!!
!((Xnode(jj) .eq. Xmax) .and. (Ynode(jj) .eq. Ymax))
       load_index(mm)=2*jj !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      load(load_index(mm))=1 
      mm=mm+1
       exit              !!!!! inja avaz shode

       end if
   end do 



   !********** Apply BC

   allocate(bc(2*Count_Xmin),u_bc(2*Count_Xmin),bc_node(Count_Xmin))
   u_bc=0  ! dar in masale tek ye gah

   counter=1
   do jj=1,n
     if (Xnode(jj) .eq. Xmin) then
      bc(2*counter-1:2*counter)=(/2*jj-1,2*jj/)
      bc_node(counter)=Xnode(jj)
      counter=counter+1
     end if

   end do




   !*********************************


!   *******************************************************************
   Ni_zeta=0

 !division = number of Load increment
 ! do ii=1,count_Xmax
! load(load_index(ii))=+(5.5e7/count_Xmax)*load(load_index(II))   !! tension
!
! end do

   load(load_index(1))=-(1e7)*load(load_index(1))   !!! bending


   Tol_inc=1e-4

 do ii=1,4*e
    C_elasplas(ii,:,:)=C_elas
 end do

 f_new=0
 delta_u=0
 f_load=0
 norm_u_new=1
 
  division=1000
  allocate(norm_ep(1000))
  norm_ep=0
  dt=0
  inc=1
  norm_R=0
  dete1=1

do while (abs(f_load(load_index(1))) < abs(load(load_index(1)) ))

2    dt=dt+inc

    R_ass=0
    f_load=dt*(1./division)*load  !CR_ass ! ! increment for load
    R_ass=R_ass-f_load
    CR_ass=R_ass

   counter=0
   norm_delta_u=1.
   norm_u_new=1.
   norm_delta_u_past=1.
   f_new=0

   ! elastic predict plastic element in inc onset =0
    delta_ep_np1=0  ! transform to ep_np1  !!!!! <<<<<<<<<<<<<<   >>>>>>>>>>>>>>>>>
!    delta_gamma=0
    ep_n33=0

!    gamma_all=0
!    beta_n=0
!    beta_33=0
!     a_n=0      !!!! ino avaz kardam baraye sakht shavandeghi hazve asaresh


!  delta_u=0

!  if (dt .eq. 0.48*division) then
!   tol_inc=1e-4
!  end if

    
! u_np1=0


!    do ii=1,4*e
!    C_elasplas(ii,:,:)=C_elas
!    end do



 do while (norm_delta_u/norm_u_new >(Tol_inc) )    

  !Loop for iterate

 counter =counter+1

 k_ass=0
 k_mat=0
 k_geo=0
 Ck_ass=0
 inv_Ck_ass=0
 f_new=0
!  ep_n33= ep_n33_conv    !!! az converge shode estefade mishe



 do ee=1,e
 !Loop for Assemble


do mm=1,4
Xi(mm)=Xnode(element(ee,mm))
Yi(mm)=Ynode(element(ee,mm))

xit(mm)=Xi(mm)+u_new(2*element(ee,mm)-1)
yit(mm)=Yi(mm)+u_new(2*element(ee,mm))

unew_local(2*mm-1)= u_new(2*element(ee,mm)-1)
unew_local(2*mm)= u_new(2*element(ee,mm))

u_local(2*mm-1)=u_new(2*element(ee,mm)-1)-u_np1(2*element(ee,mm)-1)
u_local(2*mm)=u_new(2*element(ee,mm))-u_np1(2*element(ee,mm))
end do 


 do kk=1,4   !Loop fir quarature

    ind_quad_S=4*(ee-1)+kk


 
      do ii=1,4
      Ni_zeta(ii,1:2)= (1./4)*(/  zeta_i(ii)*(1+ eta_i(ii) *Yquad_p(kk)) ,eta_i(ii)*(1+ zeta_i(ii) *Xquad_p(kk)) /)
     end do

 d_xs_zeta=0
 d_xb_zeta=0
inv_d_xb_zeta=0
inv_d_xs_zeta=0

 call derive(kk,xit,yit,d_xs_zeta)

 call derive(kk,Xi,yi,d_Xb_zeta)
 dete=0
 call det(d_Xs_zeta,dete)

!   call Gauss (d_xs_zeta,inv_d_xs_zeta,2)  
inv_d_xs_zeta= d_xs_zeta
      call SGETRF( 2, 2, inv_d_xs_zeta, 2, IPIV, INFO )
     call SGETRI( 2, inv_d_xs_zeta, 2, IPIV, WORK, 2, INFO )

  !xb =X big xs=x small
   beta(1:4,1:2)= matmul(Ni_zeta,inv_d_xs_zeta)       !   beta =d_Ni_x
 

 delta_w12=0

 do ii=1,4
    
   temp_B(1,1)=beta(ii,1)
   temp_B(2,2)=beta(ii,2)
   temp_B(3,1:2)=(/beta(ii,2) ,beta(ii,1) /)   ! dar har ghesmat in B ke local ast ba tavajoh be dada haye daroone an avaz mishavad

   Bi(1:3,(2*ii)-1:2*ii) = temp_B    ! bayad dar halghe jay girad ta 3*8 ijad shavad  => Bi(3*8)^T *C(3*3) *Bi(3*8)*u(8*1)
   delta_w12=delta_w12+ (0.5*(beta(ii,2)*u_local(2*ii-1) - beta(ii,1)*u_local(2*ii) ));
 end do


  delta_D=matmul(Bi,u_local);   ! increment strain
  strain_all(1:3,ind_quad_S)=Strain_all_conv(1:3,ind_quad_S)+delta_D
                                ! increment plastic strain
      
            cauchy_t(1,1)= cauchy_v_all(1,ind_quad_S);
            cauchy_t(2,2)=  cauchy_v_all(2,ind_quad_S);
            cauchy_t(1,2)=  cauchy_v_all(3,ind_quad_S);
            cauchy_t(2,1)=  cauchy_v_all(3,ind_quad_S) ;
        
             deltaW_dot_sigma(1)=2*delta_w12* cauchy_t(1,2);
             deltaW_dot_sigma(2)=-2*delta_w12* cauchy_t(1,2);
             deltaW_dot_sigma(3)=delta_w12* (cauchy_t(2,2)-cauchy_t(1,1));
                
              
!              matmul(C_elasplas( ind_quad_S,1:3,1:3)
!               delta_D-delta_ep_np1(1:3,ind_quad_S)
             
      cauchy_v_all(1:3,ind_quad_S)=matmul(C_elas ,delta_D)+&
                                    deltaW_dot_sigma + cauchy_v_converge(1:3,ind_quad_S)  

!   cauchy_v_all(1:2,ind_quad_S)=cauchy_v_all(1:2,ind_quad_S)-(ep_n33_conv(ind_quad_s)*(lambda)*(/1,1/) )


   


      call plastic_computation(cauchy_v_all(1:3,ind_quad_s),ep_n(1:3,ind_quad_s) ,delta_ep_np1(1:3,ind_quad_s),&
                                ep_n33_conv(ind_quad_s),ep_n33(ind_quad_s),&
                                beta_n(1:3,ind_quad_s),beta_33(ind_quad_s),a_n(ind_quad_s),norm_ep(dt),&
                                c_elasplas( ind_quad_s,:,:),a_nAll(ind_quad_S),delta_gamma(ind_quad_S),gamma_all(ind_quad_S),&
                                strain_all(1:3,ind_quad_S))


7   cauchy_t(1,1)=  cauchy_v_all(1,ind_quad_S)
   cauchy_t(2,2)=  cauchy_v_all(2,ind_quad_S)
   cauchy_t(1,2)=  cauchy_v_all(3,ind_quad_S)
   cauchy_t(2,1)=  cauchy_v_all(3,ind_quad_S) 
   
     k_mat=dete*matmul(matmul(transpose(Bi), c_elasplas( ind_quad_s,:,:)),Bi)
     f_int= dete* matmul(transpose(Bi),cauchy_v_all(1:3,ind_quad_S) )
    
    if (dt .eq. division) then
       CAUCHY(2*ee-1:2*ee,2*kk-1:2*kk)=cauchy_t
    end if

  temp_k_geo=dete*matmul(matmul(beta,cauchy_t),transpose(beta))  ! chun man beta ro 4*2 tarif karde budam
  
  do ii=1,4
   do jj=1,4
     k_geo((2*ii)-1,(2*jj)-1)=temp_k_geo(ii,jj)
     k_geo(2*ii,2*jj)=temp_k_geo(ii,jj)
    end do
   end do
 

   !!!!!Assembly Section

  do nn=1,4

    do mm=1,4
      k_Ass(2*element(ee,nn)-1:2*element(ee,nn),2*element(ee,mm)-1:2*element(ee,mm))=&
           k_Ass(2*element(ee,nn)-1:2*element(ee,nn),2*element(ee,mm)-1:2*element(ee,mm))+&
              K_geo(2*nn-1:2*nn,2*mm-1:2*mm) +K_mat(2*nn-1:2*nn,2*mm-1:2*mm)   ! k_ASS =A
      end do

   
    f_new(2*element(ee,nn)-1:2*element(ee,nn))= f_new(2*element(ee,nn)-1:2*element(ee,nn))+&
                                                 f_int(2*nn-1:2*nn)   ! force vector
   end do 
  

 


 end do !kk

 end do
     !Ready for BC



 R_ass=f_new-f_load


!        R_ass(load_index)=R_ass(load_index)+f_new(load_index)
        
!       do ii=1,size(bc)   !size ,dim,array..should be used
!!        temp=bc(ii,2)
!  !           if bc u(i)
!         do jj=1,2*n
!             R_ass(jj)=R_ass(jj)-k_Ass(jj,ii)*u_bc(ii);
!               
!         end do
!
!       end do

R_ass(bc)=u_bc;
Ck_ass=k_ass;  ! Copy the K_Assemble
   CK_ass(bc,:)=0;
    ck_ass(:,bc)=0;


    do ii=1,size(bc)
     Ck_ass(bc(ii),bc(ii))=1;
    end do

     inv_Ck_ass=Ck_ass
    

      call SGETRF( 2*n, 2*n, inv_Ck_ass, 2*n, IPIV, INFO )
     call SGETRI( 2*n, inv_Ck_ass, 2*n, IPIV, WORK, 2*n, INFO )
   !    call  Gauss (Ck_ass,inv_Ck_ass,2*n)   

     delta_u=-matmul(inv_Ck_ass,R_ass)


     u_new=u_new+delta_u

    call norm(u_new,2*n,norm_u_new)
    call norm(delta_u,2*n,norm_delta_u)
    call norm(R_ass,2*n,norm_R)




   print *, norm_u_new,norm_delta_u,norm_R,counter,dt

   
    !!*****************************************************




   if (counter >100) then 

   do ee=1,e
do kk=1,4
  ind_quad_S=4*(ee-1)+kk       
             cauchy_t(1,1)= cauchy_v_all(1,ind_quad_S);
            cauchy_t(2,2)=  cauchy_v_all(2,ind_quad_S);
            cauchy_t(1,2)=  cauchy_v_all(3,ind_quad_S);
            cauchy_t(2,1)=  cauchy_v_all(3,ind_quad_S) ;

   CAUCHY(2*ee-1:2*ee,2*kk-1:2*kk)=cauchy_t
   end do
   end do
  goto 4
 end if


  end do     !end of loop for while
   
   !in this section convergence is reached


  u_np1=u_new
  ep_n=ep_n+delta_ep_np1
  ep_n33_conv=ep_n33+ep_n33_conv  !!!! avaz shode

  cauchy_v_converge= cauchy_v_all

  Strain_all_conv= strain_all

  a_nAll=a_n
  gamma_all=gamma_all+delta_gamma;


  !!!! anmation  ***********************************************************************

  if (mod(dt,10.0) .eq. 0 ) then


  !!! compute stress

    do ee=1,e
do kk=1,4
  ind_quad_S=4*(ee-1)+kk       
             cauchy_t(1,1)= cauchy_v_all(1,ind_quad_S);
            cauchy_t(2,2)=  cauchy_v_all(2,ind_quad_S);
            cauchy_t(1,2)=  cauchy_v_all(3,ind_quad_S);
            cauchy_t(2,1)=  cauchy_v_all(3,ind_quad_S) ;

   CAUCHY(2*ee-1:2*ee,2*kk-1:2*kk)=cauchy_t
   end do
   end do

    !!! compute stress


4    numb=0

stress_vector=0   !!! dALILE ziad shodane taneshha

   do ee=1,e
   do mm=1,4


     ind_quad_S=4*(ee-1)+mm


     S_dev=CAUCHY(2*ee-1:2*ee,2*mm-1:2*mm)-(1./3.)* (  (1+nu)*(CAUCHY(2*ee-1,2*mm-1)+CAUCHY(2*ee,2*mm)) -2*mu*0*(1+nu)*(ep_n33_conv(ind_quad_S)))*identity   !!! 0 va 1+nu gozashtam 
     S_dev3=-(S_dev(1,1)+S_dev(2,2))
     stress_mises=sqrt((3./2)*(S_dev(1,1)**2+S_dev(2,2)**2+S_dev3**2+ 2*S_dev(1,2)**2))
      stress_vector (element(ee,mm))=stress_vector (element(ee,mm))+ stress_mises   !(1./4)*stress_mises
       numb(element(ee,mm))=numb(element(ee,mm))+1

!       if (ee .eq. e) then
!       ii=1
!
!       end if

   

     end do
   end do


   open(5,file=".\displacement_Fortran.plt",action='write',position ="append",Status="unknown")
write(5,*) ' TITLE = Nonlinear FEM2.tec created by Amin zaami'


   open(6,file=".\stress_mises_Fortran.plt",action='write',position ="append",Status="unknown")
write(6,*) ' TITLE = Nonlinear FEM2.tec created by Amin zaami'

!   open(7,file=".\norm_ep.txt",action='write',position ="append",Status="unknown")
!write(7,*) ' TITLE = Nonlinear FEM2.tec created by Amin zaami- norm_ep'


write(5,*) 'variables=  "X"  "Y"  "Usum" '
write(5,*) 'zone  N=   ' ,n,'  E=  '  ,e, '  F=FEPOINT  ET=QUADRILATERAL '
write(5,*)  'STRANDID=1,', ' SOLUTIONTIME= ' , dt/1000.0  


write(6,*) 'variables=  "X"  "Y"  "Stress_Mises" '
write(6,*) 'zone  N=   ' ,n,'  E=  '  ,e, '  F=FEPOINT  ET=QUADRILATERAL '
write(6,*)  'STRANDID=1,', ' SOLUTIONTIME= ' , dt/1000.0 




do ii=1,n

!if (numb(ii) .eq. 1) then
!numb(ii)=1.25   !! 75% Smoothing
!end if

write(5,*) Xnode(ii)+u_new(2*ii-1) ,Ynode(ii)+u_new(2*ii),sqrt(u_new(2*ii)**2 +u_new(2*ii-1)**2)
write(6,*) Xnode(ii)+u_new(2*ii-1) ,Ynode(ii)+u_new(2*ii),stress_vector(ii)/numb(ii)


end do







do jj=1,e
write(5,*) element (jj,:)
write(6,*) element (jj,:)
end do

!do ii=1,division
!  write(7,*)   norm_ep(ii)    !!!! for one element last
!end do

 

  !!! animation ***********************************************************************
  end if

 end do



!
!
!4    numb=0
!
!
!   do ee=1,e
!   do mm=1,4
!
!
!     ind_quad_S=4*(ee-1)+mm
!
!
!     S_dev=CAUCHY(2*ee-1:2*ee,2*mm-1:2*mm)-(1./3.)* (  (1+nu)*(CAUCHY(2*ee-1,2*mm-1)+CAUCHY(2*ee,2*mm)) -2*(1+nu)*0*(ep_n33_conv(ind_quad_S)))*identity   !!! 0 va 1+nu gozashtam 
!     S_dev3=-(S_dev(1,1)+S_dev(2,2))
!     stress_mises=sqrt((3./2)*(S_dev(1,1)**2+S_dev(2,2)**2+S_dev3**2+ 2*S_dev(1,2)**2))
!      stress_vector (element(ee,mm))=stress_vector (element(ee,mm))+ stress_mises   !(1./4)*stress_mises
!       numb(element(ee,mm))=numb(element(ee,mm))+1
!
!!       if (ee .eq. e) then
!!       ii=1
!!
!!       end if
!
!   
!
!     end do
!   end do
!
!
!   open(5,file=".\displacement.txt",action='write',position ="append",Status="unknown")
!write(5,*) ' TITLE = Nonlinear FEM2.tec created by Amin zaami'
!
!   open(6,file=".\stress_mises.txt",action='write',position ="append",Status="unknown")
!write(6,*) ' TITLE = Nonlinear FEM2.tec created by Amin zaami'
!
!   open(7,file=".\norm_ep.txt",action='write',position ="append",Status="unknown")
!write(7,*) ' TITLE = Nonlinear FEM2.tec created by Amin zaami- norm_ep'
!
!
!write(5,*) 'variables=  "X"  "Y"  "Usum" '
!write(5,*) 'zone  N=   ' ,n,'  E=  '  ,e, '  F=FEPOINT  ET=QUADRILATERAL '
!
!write(6,*) 'variables=  "X"  "Y"  "Stress_Mises" '
!write(6,*) 'zone  N=   ' ,n,'  E=  '  ,e, '  F=FEPOINT  ET=QUADRILATERAL '
!
!
!do ii=1,n
!
!!if (numb(ii) .eq. 1) then
!!numb(ii)=1.25   !! 75% Smoothing
!!end if
!
!write(5,*) Xnode(ii)+u_new(2*ii-1) ,Ynode(ii)+u_new(2*ii),sqrt(u_new(2*ii)**2 +u_new(2*ii-1)**2)
!write(6,*) Xnode(ii)+u_new(2*ii-1) ,Ynode(ii)+u_new(2*ii),stress_vector(ii)/numb(ii)
!
!
!end do
!
!
!
!
!
!
!
!do jj=1,e
!write(5,*) element (jj,:)
!write(6,*) element (jj,:)
!end do
!
!do ii=1,division
!  write(7,*)   norm_ep(ii)    !!!! for one element last
!end do
!
! 


    print *, 'Hello World'

    end program FEM2_p1

