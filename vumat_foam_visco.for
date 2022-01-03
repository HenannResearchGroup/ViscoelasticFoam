!***********************************************************************
!
! User material subroutine (VUMAT) for the large-deformation, 
!  viscoelastic behavior of compressible, elastomeric foams. 
!  This VUMAT is for use with the Dynamic/Explicit step in 
!  Abaqus/Explicit. This VUMAT is not for use in plane stress or in 
!  any other situation in which there are more strain terms than 
!  stress terms.
!
! Xiuqi Li, Dec 2020
! David Henann, July 2021
!
!***********************************************************************
! Usage:
!***********************************************************************
!
!     Material Properties Vector (*user material, constants = nprops)
!       nprops = 15+nbranch*6, where nbranch = props(15) is the 
!       number of non-equilibrium branches
!     --------------------------------------------------------------
!     Equilibrium parameters
!       
!       G0      = props(1)  ! Shear modulus
!       B       = props(2)  ! Bulk modulus
!       Jmin    = props(3)  ! Locking value of J in f-fuction
!       C1      = props(4)  ! C_1 in L-function
!       K10     = props(5)  ! K1_0 in X-function
!       delta_K = props(6)  ! delta_K in X-function
!       X1      = props(7)  ! X1' in X-function
!       X2      = props(8)  ! X2' in X-function
!       C0      = props(9)  ! C_0 in L-function
!       p       = props(10) ! power p in L-function
!       q       = props(11) ! power q in L-function
!       C2      = props(12) ! C_2 in f-function
!       C3      = props(13) ! C_3 in f-function
!       r       = props(14) ! power r in f-function
!
!     Non-equilibrium parameters
!       
!       nbranch = props(15) ! number of non-equilibrium branches
!       
!       do i = 1,nbranch
!         
!         branchFlag = props(16+(i-1)*6)
!         
!         if (branchFlag.eq.1) then
!           non-Newtonian branch
!           Gneq       = props(17+(i-1)*6) ! non-equilibrium shear modulus
!           m          = props(18+(i-1)*6) ! strain-rate sensitivity exponent
!           gamma0_dot = props(19+(i-1)*6) ! reference shear strain rate
!           C1_neq     = props(20+(i-1)*6) ! C_1 in L_neq-function
!           p_neq      = props(21+(i-1)*6) ! power p in L_neq-function
!         elseif (branchFlag.eq.2) then
!           Newtonian branch
!           Gneq       = props(17+(i-1)*6) ! non-equilibrium shear modulus
!           t0         = props(18+(i-1)*6) ! relaxation time
!           *** props(19+(i-1)*6) not used for a Newtonian branch ***
!           C1_neq     = props(20+(i-1)*6) ! C_1 in L_neq-function
!           p_neq      = props(21+(i-1)*6) ! power p in L_neq-function
!         end if
!           
!       end do
!
!     State Variables (*depvar nstatev) 
!       nstatev = 3+nbranch*15, where nbranch = props(15) is the 
!       number of non-equilibrium branches
!     --------------------------------------------------------------
!     statev(1) = K1 -- ln(J), range(-inf,inf)
!     statev(2) = K2 -- amount of distortional deformation, range[0,inf)
!     statev(3) = K3 -- deformation mode, range[-1,1], 
!                           -1 is compression, 0 is shear, 1 is tension
!     
!     do i = 1,nbranch
!       statev(4+(i-1)*15)  = Fv(1,1) -- viscous deformation gradient(1,1)
!       statev(5+(i-1)*15)  = Fv(1,2) -- viscous deformation gradient(1,2)
!       statev(6+(i-1)*15)  = Fv(1,3) -- viscous deformation gradient(1,3)
!       statev(7+(i-1)*15)  = Fv(2,1) -- viscous deformation gradient(2,1)
!       statev(8+(i-1)*15)  = Fv(2,2) -- viscous deformation gradient(2,2)
!       statev(9+(i-1)*15)  = Fv(2,3) -- viscous deformation gradient(2,3)
!       statev(10+(i-1)*15) = Fv(3,1) -- viscous deformation gradient(3,1)
!       statev(11+(i-1)*15) = Fv(3,2) -- viscous deformation gradient(3,2)
!       statev(12+(i-1)*15) = Fv(3,3) -- viscous deformation gradient(3,3)
!       statev(13+(i-1)*15) = Me(1,1) -- Mandel stress(1,1)
!       statev(14+(i-1)*15) = Me(2,2) -- Mandel stress(2,2)
!       statev(15+(i-1)*15) = Me(3,3) -- Mandel stress(3,3)
!       statev(16+(i-1)*15) = Me(2,3) -- Mandel stress(2,3)
!       statev(17+(i-1)*15) = Me(1,3) -- Mandel stress(1,3)
!       statev(18+(i-1)*15) = Me(1,2) -- Mandel stress(1,2)
!     end do
!
!***********************************************************************
      subroutine vumat(
      ! Read only (unmodifiable)variables -
     +  nblock, ndir, nshr, nstatev, nfieldv, nprops, jInfoArray,
     +  stepTime, totalTime, dtArray, cmname, coordMp, charLength,
     +  props, density, strainInc, relSpinInc,
     +  tempOld, stretchOld, defgradOld, fieldOld,
     +  stressOld, stateOld, enerInternOld, enerInelasOld,
     +  tempNew, stretchNew, defgradNew, fieldNew,
      ! Write only (modifiable) variables -
     +  stressNew, stateNew, enerInternNew, enerInelasNew )
      !
      include 'vaba_param.inc'
      !
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     +  charLength(nblock), dtArray(2*(nblock)+1), 
     +  strainInc(nblock,ndir+nshr), relSpinInc(nblock,nshr), 
     +  tempOld(nblock), stretchOld(nblock,ndir+nshr),
     +  defgradOld(nblock,ndir+nshr+nshr),
     +  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     +  stateOld(nblock,nstatev), enerInternOld(nblock),
     +  enerInelasOld(nblock), tempNew(nblock),
     +  stretchNew(nblock,ndir+nshr),
     +  defgradNew(nblock,ndir+nshr+nshr),
     +  fieldNew(nblock,nfieldv),
     +  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     +  enerInternNew(nblock), enerInelasNew(nblock), jInfoArray(*)
      !
      character*80 cmname
      !
      ! Variables defined and used in the VUMAT
      !
      integer i,km,istat,nbranch,dummyFlag
      !
      real*8 properties(nprops),dtime,F_t(3,3),F_tau(3,3),U_t(3,3),
     +  U_tau(3,3),T_tau(3,3),K1,K2,K3,psi,Fv_t(3,3),Me_t(3,3),
     +  T_neq(3,3),Fv_tau(3,3),Me_tau(3,3),psi_neq,vwrkincBranch,
     +  vwrkinc,U_tau_inv(3,3),det_U,R_tau(3,3),stress(3,3)
      !
      real*8 zero,one,two,three,fourth,third,half,six,four
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0,
     +  six=6.d0,third=1.d0/3.d0,half=1.d0/2.d0,fourth=1.d0/4.d0)
      
      
      ! Initialization
      !
      properties = props
      dtime = dtArray(1)
      stress = zero
      
      
      ! Begin the loop over integration points in nblock
      !
      do km = 1,nblock
        !
        ! Copy the old and new deformation gradients into F_t and 
        !   F_tau, respectively. Copy the old and new stretches into 
        !   U_t and U_tau, respectively.
        !
        F_t(1,1) = defgradOld(km,1)
        F_t(2,2) = defgradOld(km,2)
        F_t(3,3) = defgradOld(km,3)
        F_t(1,2) = defgradOld(km,4)
        !
        F_tau(1,1) = defgradNew(km,1)
        F_tau(2,2) = defgradNew(km,2)
        F_tau(3,3) = defgradNew(km,3)
        F_tau(1,2) = defgradNew(km,4)
        !
        U_t(1,1) = stretchOld(km,1)
        U_t(2,2) = stretchOld(km,2)
        U_t(3,3) = stretchOld(km,3)
        U_t(1,2) = stretchOld(km,4)
        !
        U_tau(1,1) = stretchNew(km,1)
        U_tau(2,2) = stretchNew(km,2)
        U_tau(3,3) = stretchNew(km,3)
        U_tau(1,2) = stretchNew(km,4)
        !
        if (nshr.eq.1) then
          !
          F_t(2,3) = zero
          F_t(3,1) = zero
          F_t(2,1) = defgradOld(km,5)
          F_t(3,2) = zero
          F_t(1,3) = zero
          !
          F_tau(2,3) = zero
          F_tau(3,1) = zero
          F_tau(2,1) = defgradNew(km,5)
          F_tau(3,2) = zero
          F_tau(1,3) = zero
          !
          U_t(2,3) = zero
          U_t(3,1) = zero
          U_t(2,1) = U_t(1,2)
          U_t(3,2) = zero
          U_t(1,3) = zero
          !
          U_tau(2,3) = zero
          U_tau(3,1) = zero
          U_tau(2,1) = U_tau(1,2)
          U_tau(3,2) = zero
          U_tau(1,3) = zero
          !
        else
          !
          F_t(2,3) = defgradOld(km,5)
          F_t(3,1) = defgradOld(km,6)
          F_t(2,1) = defgradOld(km,7)
          F_t(3,2) = defgradOld(km,8)
          F_t(1,3) = defgradOld(km,9)
          !
          F_tau(2,3) = defgradNew(km,5)
          F_tau(3,1) = defgradNew(km,6)
          F_tau(2,1) = defgradNew(km,7)
          F_tau(3,2) = defgradNew(km,8)
          F_tau(1,3) = defgradNew(km,9)
          !
          U_t(2,3) = stretchOld(km,5)
          U_t(3,1) = stretchOld(km,6)
          U_t(2,1) = U_t(1,2)
          U_t(3,2) = U_t(2,3)
          U_t(1,3) = U_t(3,1)
          !
          U_tau(2,3) = stretchNew(km,5)
          U_tau(3,1) = stretchNew(km,6)
          U_tau(2,1) = U_tau(1,2)
          U_tau(3,2) = U_tau(2,3)
          U_tau(1,3) = U_tau(3,1)
          !
        end if
        
        
        ! Determine the number of non-equilibrium branches
        !
        nbranch = int(props(15))
        
        
        ! At the start of an Abaqus calculation, the state variables
        !   are passed into the VUMAT subroutine with zero values.
        !   Initialize the state variables. At this point, stepTime 
        !   and totalTime both have a value of zero.
        !
        dummyFlag = 0
        if ((stepTime.eq.zero).and.(totalTime.eq.zero)) then
          !
          stateOld(km,1) = zero
          stateOld(km,2) = zero
          stateOld(km,3) = zero
          !
          do i = 1,nbranch
            stateOld(km,4+(i-1)*15)  = one
            stateOld(km,5+(i-1)*15)  = zero
            stateOld(km,6+(i-1)*15)  = zero
            stateOld(km,7+(i-1)*15)  = zero
            stateOld(km,8+(i-1)*15)  = one
            stateOld(km,9+(i-1)*15)  = zero
            stateOld(km,10+(i-1)*15) = zero
            stateOld(km,11+(i-1)*15) = zero
            stateOld(km,12+(i-1)*15) = one
            stateOld(km,13+(i-1)*15) = zero
            stateOld(km,14+(i-1)*15) = zero
            stateOld(km,15+(i-1)*15) = zero
            stateOld(km,16+(i-1)*15) = zero
            stateOld(km,17+(i-1)*15) = zero
            stateOld(km,18+(i-1)*15) = zero
          end do
          !
          dummyFlag = 1 ! This is a dummy step
          !
        end if
        
        
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        !
        ! Perform the constitutive update
        !
        ! Calculate the stress for the equilibrium branch
        !
        call eqBranch(properties,nprops,dtime,F_tau,
     +                T_tau,K1,K2,K3,psi,istat)
        !
        ! Update the state variables for the equilibrium branch
        !
        stateNew(km,1) = K1
        stateNew(km,2) = K2
        stateNew(km,3) = K3
        
        
        ! Loop over the non-equilibrium branches
        !
        vwrkinc = zero
        do i = 1,nbranch
          ! 
          ! Get the state variables for this non-equilibrium branch
          !
          Fv_t(1,1) = stateOld(km,4+(i-1)*15)
          Fv_t(1,2) = stateOld(km,5+(i-1)*15)
          Fv_t(1,3) = stateOld(km,6+(i-1)*15)
          Fv_t(2,1) = stateOld(km,7+(i-1)*15)
          Fv_t(2,2) = stateOld(km,8+(i-1)*15)
          Fv_t(2,3) = stateOld(km,9+(i-1)*15)
          Fv_t(3,1) = stateOld(km,10+(i-1)*15)
          Fv_t(3,2) = stateOld(km,11+(i-1)*15)
          Fv_t(3,3) = stateOld(km,12+(i-1)*15)
          !
          Me_t(1,1) = stateOld(km,13+(i-1)*15)
          Me_t(2,2) = stateOld(km,14+(i-1)*15)
          Me_t(3,3) = stateOld(km,15+(i-1)*15)
          Me_t(2,3) = stateOld(km,16+(i-1)*15)
          Me_t(3,2) = Me_t(2,3)
          Me_t(1,3) = stateOld(km,17+(i-1)*15)
          Me_t(3,1) = Me_t(1,3)
          Me_t(1,2) = stateOld(km,18+(i-1)*15)
          Me_t(2,1) = Me_t(1,2)
          !
          ! Calculate the stress for this non-equilibrium branch
          !
          call neqBranch(i,properties,nprops,dtime,F_tau,
     +                   Fv_t,Me_t,dummyFlag,
     +                   T_neq,Fv_tau,Me_tau,
     +                   psi_neq,vwrkincBranch,istat)
          !
          ! Update the state variables for this non-equilibrium branch
          !
          stateNew(km,4+(i-1)*15)  = Fv_tau(1,1)
          stateNew(km,5+(i-1)*15)  = Fv_tau(1,2)
          stateNew(km,6+(i-1)*15)  = Fv_tau(1,3)
          stateNew(km,7+(i-1)*15)  = Fv_tau(2,1)
          stateNew(km,8+(i-1)*15)  = Fv_tau(2,2)
          stateNew(km,9+(i-1)*15)  = Fv_tau(2,3)
          stateNew(km,10+(i-1)*15) = Fv_tau(3,1)
          stateNew(km,11+(i-1)*15) = Fv_tau(3,2)
          stateNew(km,12+(i-1)*15) = Fv_tau(3,3)
          !
          stateNew(km,13+(i-1)*15) = Me_tau(1,1)
          stateNew(km,14+(i-1)*15) = Me_tau(2,2)
          stateNew(km,15+(i-1)*15) = Me_tau(3,3)
          stateNew(km,16+(i-1)*15) = Me_tau(2,3)
          stateNew(km,17+(i-1)*15) = Me_tau(1,3)
          stateNew(km,18+(i-1)*15) = Me_tau(1,2)
          !
          ! Add the contributions to the Cauchy stress, free energy
          !  per unit deformed volume, and viscous work increment 
          !  per unit deformed volume due to this branch
          !
          T_tau = T_tau + T_neq
          psi = psi + psi_neq
          vwrkinc = vwrkinc + vwrkincBranch
          !
        end do
        !
        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        
        
        ! Update the stress measure, (R^T) T R, used by Abaqus/Explicit
        !
        call matInv3D(U_tau,U_tau_inv,det_U,istat)
        R_tau = matmul(F_tau,U_tau_inv)
        stress = matmul(matmul(transpose(R_tau),T_tau),R_tau)
        !
        do i = 1,ndir
          stressNew(km,i) = stress(i,i)
        end do
        !
        if (nshr.ne.0) then
          stressNew(km,ndir+1) = stress(1,2)
          if (nshr.ne.1) then
            stressNew(km,ndir+2) = stress(2,3)
            if (nshr.ne.2) then
              stressNew(km,ndir+3) = stress(1,3)
            end if
          end if
        end if
        
        
        ! Update the energy outputs
        !
        enerInternNew(km) = psi/density(km)
        enerInelasNew(km) = enerInelasOld(km) + vwrkinc/density(km)
      
      
      end do
      
      
      return
      end subroutine vumat
      
!***********************************************************************
!     Material subroutines
!***********************************************************************
      
      subroutine eqBranch(
     +        ! Inputs
     +        props,nprops,dtime,F_tau,
     +        ! Outputs
     +        T_tau,K1,K2,K3,psi,istat)
      
      implicit none
      !
      integer nprops,istat
      !
      real*8 props(nprops),dtime,F_tau(3,3),T_tau(3,3),K1,K2,K3,
     +  psi,Iden(3,3),G0,B,Jmin,C1,K10,delta_K,X1,X2,
     +  C0,p,q,C2,C3,r,detF,Rot(3,3),U_tau(3,3),E_tau(3,3),
     +  dev_E(3,3),N(3,3),Y(3,3),f,dfdK1,dfdK1_2,X,dXdK1,dXdK1_2,Lf,
     +  dLdK2,dLdK3,dLdK2dK2,dLdK2dK3,dpsidK1,dpsidK2,dpsidK3
      !
      real*8 zero,one,two,three,fourth,third,half,six,four
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0,
     +     six=6.d0,third=1.d0/3.d0,half=1.d0/2.d0,fourth=1.d0/4.d0)
      
      
      ! Identity matrix
      !
      call onem(Iden)
      
      
      ! Obtain material properties
      !
      G0 = props(1) ! Shear modulus
      B = props(2) ! Bulk modulus
      Jmin = props(3) ! Locking value of J in f-fuction
      C1 = props(4) ! C_1 in L-function
      K10 = props(5) ! K1_0 in X-function
      delta_K = props(6) ! delta_K in X-function
      X1 = props(7) ! X1' in X-function
      X2 = props(8) ! X2' in X-function
      C0 = props(9) ! C_0 in L-function
      p = props(10) ! power p in L-function
      q = props(11) ! power q in L-function
      C2 = props(12) ! C_2 in f-function
      C3 = props(13)! C_3 in f-function
      r = props(14) ! power r in f-function
      
      
      ! Calulate kinematic quantities
      !
      call mdet(F_tau,detF)
      call skinem(F_tau,Rot,U_tau,E_tau,istat)
      E_tau = matmul(matmul(Rot,E_tau),transpose(Rot))
      
      
      ! Calulate the invariants K1, K2, K3 
      !
      K1 = E_tau(1,1) + E_tau(2,2) + E_tau(3,3)	  
      dev_E = E_tau - third*K1*Iden
      K2 = dsqrt(sum(dev_E*dev_E))
      if (K2.eq.zero) then
        ! if K2 == 0, assign an arbitrary value to N, K3, and Y
        N = zero
        K3 = zero
        Y = zero
      else
        ! else compute N, K3, and Y
        N = dev_E/K2
        call mdet(N,K3)
        K3 = K3*three*dsqrt(six)  
        Y  = three*dsqrt(six)*matmul(N,N) 
     +              - dsqrt(six)*Iden - three*K3*N
      end if 
      
      
      ! Compute the free energy function and its derivatives
      !
      call f_fun(f,dfdK1,dfdK1_2,K1,Jmin,C2,C3,r)
      call X_fun(X,dXdK1,dXdK1_2,K1,X1,X2,K10,delta_K)
      call L_fun(Lf,dLdK2,dLdK3,dLdK2dK2,dLdK2dK3,K2,K3,C0,C1,p,q)
      dpsidK1 = G0*dXdK1*K2**two + B*dfdK1
      dpsidK2 = G0*(two*X*K2 + dLdK2)
      dpsidK3 = G0*dLdK3
      psi = G0*(X*(K2**two) + Lf) + B*f
      
      
      ! Compute the Cauchy stress
      !
      if (K2.eq.zero) then
        !
        T_tau = (one/detF)*(dpsidK1*Iden) 
        !
      else
        !
        T_tau = (one/detF)*(dpsidK1*Iden + dpsidK2*N + dpsidK3/K2*Y)
        !
      end if
      
      
      ! Calculate the equilibrium contribution to the free energy
      !   PER UNIT DEFORMED VOLUME.
      !
      psi = psi/detF
      
      
      return
      end subroutine eqBranch 
      
!***********************************************************************
      
      subroutine neqBranch(
     +        ! Inputs
     +        branchNum,props,nprops,dtime,F_tau,Fv_t,Me_t,dummyFlag,
     +        ! Outputs
     +        T_neq,Fv_tau,Me_tau,psi_neq,vwrkinc,istat)
      
      implicit none
      !
      integer branchNum,nprops,dummyFlag,istat,branchFlag
      !
      real*8 props(nprops),dtime,F_tau(3,3),Fv_t(3,3),Me_t(3,3),
     +  T_neq(3,3),Fv_tau(3,3),Me_tau(3,3),psi_neq,vwrkinc,Iden(3,3),
     +  G0,B,Gneq,Bneq,Jmin,C1,K10,delta_K,X1,X2,C0,p,q,C2,C3,r,
     +  trMe_t,Me0_t(3,3),tauBar_t,Nv_t(3,3),m,gamma0_dot,eps0_dot,tau0,
     +  sigma0,gamma_dot_v,eps_dot_v,Dv_t(3,3),t0,eta,kappa,
     +  Fv_tau_inv(3,3),det_Fv_tau,Fe_tau(3,3),det_Fe_tau,Rot(3,3),
     +  Ue_tau(3,3),Ee_tau(3,3),K1,K2,K3,dev_E(3,3),N(3,3),Y(3,3),
     +  f,dfdK1,dfdK1_2,X,dXdK1,dXdK1_2,Lf,dLdK2,dLdK3,dLdK2dK2,
     +  dLdK2dK3,dpsidK1,dpsidK2,dpsidK3
      !
      real*8 zero,one,two,three,fourth,third,half,six,four
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0,
     +     six=6.d0,third=1.d0/3.d0,half=1.d0/2.d0,fourth=1.d0/4.d0)
      
      
      ! Identity matrix
      !
      call onem(Iden)
      
      
      ! Obtain material properties
      !
      G0 = props(1) ! Equilibrium shear modulus
      B = props(2) ! Equilibrium bulk modulus
      Gneq = props(17+(branchNum-1)*6) ! Non-equilibrium shear modulus
      Bneq =  (B/G0)*Gneq ! Non-equilibrium bulk modulus
      Jmin = props(3) ! Locking value of J in f-fuction
      C1 = props(20+(branchNum-1)*6) ! C_1 in L-function
      K10 = props(5) ! K1_0 in X-function
      delta_K = props(6) ! delta_K in X-function
      X1 = props(7) ! X1' in X-function
      X2 = props(8) ! X2' in X-function
      C0 = props(9) ! C_0 in L-function
      p = props(10) ! power p in L-function
      q = props(21+(branchNum-1)*6) ! power q in L-function
      C2 = props(12) ! C_2 in f-function
      C3 = props(13)! C_3 in f-function
      r = props(14) ! power r in f-function
      
      
      ! Calculate quantities related to the Mandel stress
      !
      trMe_t = Me_t(1,1) + Me_t(2,2) + Me_t(3,3) ! Trace of Mandel stress
      Me0_t = Me_t - third*trMe_t*Iden ! Mandel stress deviator
      tauBar_t = dsqrt(sum(Me0_t*Me0_t)/2.d0) ! Equiv. shear stress
      if(tauBar_t.le.zero) then
         Nv_t = zero
      else
         Nv_t = dsqrt(half)*(Me0_t/tauBar_t) ! Deviatoric viscous flow direction
      endif
      
      
      ! Get the branch flag
      !
      branchFlag = int(props(16+(branchNum-1)*6)) ! non-Newtonian = 1, Newtonian = 2


      ! Calculate Dv_t
      !
      if (branchFlag.eq.1) then
        !
        ! non-Newtonian branch
        ! 
        m = props(18+(branchNum-1)*6) ! strain-rate sensitivity exponent
        gamma0_dot = props(19+(branchNum-1)*6) ! reference shear strain rate
        eps0_dot = gamma0_dot ! reference volumetric strain rate
        tau0 = 1.0d5 ! reference shear stress (Pa)
        sigma0 = tau0 ! reference mean stress
        gamma_dot_v = gamma0_dot*abs(tauBar_t/tau0)**(one/m)
        eps_dot_v = eps0_dot*abs(trMe_t/(three*sigma0))**(one/m)
        if (trMe_t>zero) then
          Dv_t = dsqrt(half)*gamma_dot_v*Nv_t + eps_dot_v*Iden/three
        else
          Dv_t = dsqrt(half)*gamma_dot_v*Nv_t - eps_dot_v*Iden/three
        end if
        !
      elseif (branchFlag.eq.2) then
        !
        ! Newtonian branch
        !
        t0 = props(18+(branchNum-1)*6) ! relaxation time
        eta = Gneq*t0 ! shear viscosity
        kappa = eta ! bulk viscosity
        Dv_t = Me0_t/(two*eta) + trMe_t*Iden/(9.d0*kappa)
        !
      else
        !
        write(*,*) 'Invalid branchFlag: branchFlag.ne.1 or 2'
        call xplb_exit
        !
      end if
      
      
      ! For the initial dummy step, freeze viscous deformation
      !
      if (dummyFlag.eq.1) then
        Dv_t = zero
      end if
      
      
      ! Update Fv and Fe
      !
      Fv_tau = matmul((Iden + dtime*Dv_t),Fv_t)
      call matInv3D(Fv_tau,Fv_tau_inv,det_Fv_tau,istat)
      Fe_tau = matmul(F_tau,Fv_tau_inv)
      
      
      ! Calulate kinematic quantities
      !
      call mdet(Fe_tau,det_Fe_tau)
      call skinem(Fe_tau,Rot,Ue_tau,Ee_tau,istat)
      Ee_tau = matmul(matmul(Rot,Ee_tau),transpose(Rot))
      
      
      ! Calulate the invariants K1, K2, K3 
      !
      K1 = Ee_tau(1,1) + Ee_tau(2,2) + Ee_tau(3,3)
      dev_E = Ee_tau - third*K1*Iden
      K2 = dsqrt(sum(dev_E*dev_E))
      if (K2.eq.zero) then
        ! if K2 == 0, assign an arbitrary value to N, K3, and Y
        N = zero
        K3 = zero
        Y = zero
      else
        ! else compute N, K3, and Y
        N = dev_E/K2
        call mdet(N,K3)
        K3 = K3*three*dsqrt(six)  
        Y  = three*dsqrt(six)*matmul(N,N) 
     +              - dsqrt(six)*Iden - three*K3*N
      end if 
      
      
      ! Compute the free energy function and its derivatives
      !
      call f_fun(f,dfdK1,dfdK1_2,K1,Jmin,C2,C3,r)
      call X_fun(X,dXdK1,dXdK1_2,K1,X1,X2,K10,delta_K)
      call L_fun(Lf,dLdK2,dLdK3,dLdK2dK2,dLdK2dK3,K2,K3,C0,C1,p,q)
      dpsidK1 = Gneq*dXdK1*K2**two + Bneq*dfdK1
      dpsidK2 = Gneq*(two*X*K2 + dLdK2)
      dpsidK3 = Gneq*dLdK3
      psi_neq = Gneq*(X*(K2**two) + Lf) + Bneq*f
      
      
      ! Compute the Cauchy stress
      !
      if (K2.eq.zero) then
        !
        T_neq = (one/det_Fe_tau)*(dpsidK1*Iden) 
        !
      else
        !
        T_neq = (one/det_Fe_tau)*(dpsidK1*Iden+dpsidK2*N+dpsidK3/K2*Y)
        !
      end if
      
      
      ! Update the Mandel stress
      !
      Me_tau = matmul(transpose(Rot),matmul(det_Fe_tau*T_neq,Rot))
      
      
      ! Calculate the non-equilibrium contribution to the free energy
      !   PER UNIT DEFORMED VOLUME and the viscous work increment,
      !   also per unit deformed volume.
      !
      psi_neq = psi_neq/det_Fe_tau
      vwrkinc = sum(Me_tau*Dv_t)/det_Fe_tau
      
      
      ! For the initial dummy step, do not update the internal variables
      !
      if (dummyFlag.eq.1) then
        Fv_tau = Fv_t
        Me_tau = Me_t
      end if
      
      
      return
      end subroutine neqBranch 
      
!***********************************************************************
      
      subroutine f_fun(f,dfdK1,dfdK1_2,K1,Jmin,C2,C3,r)
      !
      ! This subroutine calculates the f-function and its derivatives
      !
      implicit none
      !
      real*8 f,dfdK1,dfdK1_2,K1,Jmin,B,C2,C3,r,J,delta_J,
     +    zeroth, first, second, thirdorder, fourthorder
      !
      real*8 zero,one,two,three,fourth,third,sixth,half,six,four,
     + twentyfourth
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,four=4.d0,
     +     six=6.d0,third=1.d0/3.d0,half=1.d0/2.d0,fourth=1.d0/4.d0,
     +     sixth = 1.d0/6.d0, twentyfourth=1.d0/24.d0)


      J = dexp(K1)


      ! For J<(Jmin+delta_J) extrapolate the f-function using a
      !  fourth order polynomial. Compute the coefficients in 
      !  the Taylor Expansion
      !
      delta_J = 0.0023d0
      !
      zeroth = (C3/(r-one))*(Jmin + ((one -
     + Jmin)**r)*((delta_J)**(one - r)) - (Jmin + delta_J)**(one - r)) +
     +  (one/C2**two)*(-one + ((Jmin + delta_J)**C2)
     +  - C2*dlog(Jmin + delta_J))
      !
      first = (((-one + (Jmin + delta_J)**C2)/C2)
     +    + C3*(Jmin + delta_J)
     +  *( -((one - Jmin)**r)*(delta_J**(-r))
     +  + (Jmin + delta_J)**(-r)))
      !
      second = ((Jmin + delta_J)**C2) +
     + C3*(Jmin + delta_J)*((one - r)
     +   *((Jmin + delta_J)**(-r)) + ((one - Jmin)**(r))
     + *(delta_J**(-one - r))*(Jmin*r + (r - one)*delta_J))
      !
      thirdorder = C2*(Jmin + delta_J)**C2 + C3*(Jmin + delta_J)*
     + (((-one + r)**two)*((Jmin + delta_J)**(-r)) - ((one - Jmin)**r)*
     + (delta_J**(-two - r))*((Jmin**two)*r*(one + r) + Jmin*r*(-one +
     + two*r)*delta_J + ((-one + r)**two)*(delta_J**two)))
      !
      fourthorder = (C2**two)*((Jmin + delta_J)**C2) +
     + C3*(Jmin + delta_J)*
     + (-((-one + r)**three)*((Jmin + delta_J)**(-r)) +
     + ((one - Jmin)**r)*(delta_J**(-three - r))*((Jmin**three)*r*
     + (one + r)*(two + r) + three*(Jmin**two)*(r**two)*(one + r)*
     + delta_J + Jmin*r*(one + three*(-one + r)*r)*(delta_J**two) +
     + ((-one + r)**three)*(delta_J**three)))


      ! Compute f
      !
      if ((K1 - dlog(Jmin + delta_J)).ge.zero) then
        f = (dexp(C2*(K1))-C2*K1 - one)/(C2**two) + (C3/(r-one))*
     +      (-J**(one-r) + (J-Jmin)**(one-r)*(one-Jmin)**r + Jmin)
      else
        f = zeroth + first*(K1 - dlog(Jmin + delta_J)) + half*second*
     +     ((K1 - dlog(Jmin + delta_J))**two) + sixth*thirdorder* 
     +     ((K1 - dlog(Jmin + delta_J))**three) +
     +     twentyfourth*fourthorder*((K1 - dlog(Jmin + delta_J))**four)
      end if


      ! Compute first derivative of f wrt K1
      !
      if ((K1 - dlog(Jmin + delta_J)).ge.zero) then
        dfdK1 = (dexp(C2*(K1))-one)/C2 +
     +  C3*J*(J**(-r)-(one-Jmin)**r/(J-Jmin)**r)
      else
        dfdK1 = first + second*(K1 -
     +     dlog(Jmin + delta_J))+ half*thirdorder*
     +     ((K1 - dlog(Jmin + delta_J))**two) +
     +     sixth*fourthorder*((K1 - dlog(Jmin + delta_J))**three)
      end if


      ! Compute second derivative of f wrt K1
      !
      if ((K1 - dlog(Jmin + delta_J)).ge.zero) then
        dfdK1_2 = dexp(C2*K1) + C3*(J**(-r)-(one-Jmin)**(r)
     +     /(J-Jmin)**r)*dexp(K1)
     +     + C3*(-r*J**(-r-one) + r*(one-Jmin)**(r)/(J-Jmin)**(r+one))
     +     *dexp(K1)*dexp(K1)
      else
        dfdK1_2 = second + thirdorder*(K1 - dlog(Jmin + delta_J)) +
     +           half*fourthorder*((K1 - dlog(Jmin + delta_J))**two)
      end if


      return
      end subroutine f_fun
      
!***********************************************************************
      
      subroutine X_fun(X,dXdK1,dXdK1_2,K1,X1,X2,K10,delta_K)
      !
      ! This subroutine calculates the X-function and its derivatives
      !
      implicit none
      !
      real*8 K1,X,dXdK1,X1,X2,delta_K,K10,dXdK1_2
      !
      real*8 zero,one,two,three,fourth,third,half,pi
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,fourth=1.d0/4.d0,
     +     third=1.d0/3.d0,half=1.d0/2.d0, pi=3.141592653d0)


      ! Compute X
      X = (X1+X2)/two*K1 + (X1-X2)/two*delta_K*
     +      dlog(dcosh((K1-K10)/delta_K)/dcosh(K10/delta_K)) + one


      ! Compute first derivative of X wrt K1
      !
      dXdK1 = (X1+X2)/two + (X1-X2)/two*dtanh((K1 - K10)/delta_K)
      
      
      ! Compute second derivative of X wrt K1
      !
      dXdK1_2 = (X1-X2)/two/dcosh((K10-K1)/delta_K)**two/delta_K
      
      
      return
      end subroutine X_fun
      
!***********************************************************************
      
      subroutine L_fun(Lf,dLdK2,dLdK3,dLdK2dK2,dLdK2dK3,K2,K3,C0,C1,p,q)
      !
      ! This subroutine calculates the L-function and its derivatives
      !
      implicit none
      !
      real*8 Lf,dLdK2,dLdK3,dLdK2dK2,dLdK2dK3,K2,K3,C1,C0,p,q
      !
      real*8 zero,one,two,three,fourth,third,half,pi
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,fourth=1.d0/4.d0,
     +     third=1.d0/3.d0,half=1.d0/2.d0, pi=3.141592653d0)
      
      
      ! Compute L
      !
      Lf = C0*K2**p + C1*(one+K3)*K2**q
      
      
      ! Compute first derivative of L wrt K2
      !
      dLdK2 = C0*p*K2**(p-one) + C1*(one+K3)*q*K2**(q-one)
      
      
      ! Compute first derivative of L wrt K3
      !
      dLdK3 = C1*K2**q 
      
      
      ! Compute second derivative of L wrt K2
      !
      dLdK2dK2 = C0*p*(p-one)*K2**(p-two) + 
     +           C1*(one+K3)*q*(q-one)*K2**(q-two)
      
      
      ! Compute mixed derivative of L wrt K2 and K3
      !
      dLdK2dK3 = C1*q*K2**(q-one)
      
      
      return
      end subroutine L_fun
      
!***********************************************************************
!     The next subroutine calculates various kinematical quantities
!      associated with the deformation gradient
!***********************************************************************
      
      subroutine skinem(F,R,U,E,istat)
      !
      ! This subroutine performs the right polar decomposition
      !  F = RU of the deformation gradient F into a rotation
      !  R and the right stretch tensor U.  The logarithmic 
      !  strain E = ln(U) is also returned.
      !
      !	F(3,3):       the deformation gradient; input
      !	detF:         the determinant of F; detF > 0
      !	R(3,3):       the rotation matrix; output
      !	U(3,3):       the right stretch tensor; output
      !	Uinv(3,3):    the inverse of U
      !	C(3,3):       the right Cauchy-Green tensor
      !	omega(3):     the squares of the principal stretches
      ! Ueigval(3):   the principal stretches
      !	eigvec(3,3):  matrix of eigenvectors of U
      !	E(3,3):       the logarithmic strain tensor; output
      ! istat:        success flag, istat=0 for a failed attempt; output
      !
      implicit none
      !
      integer istat
      !
      real*8 F(3,3),C(3,3),omega(3),Ueigval(3),eigvec(3,3),
     +  U(3,3),E(3,3),Uinv(3,3),R(3,3),detF
      
      
      !	Store the identity matrix in R, U, and Uinv
      !
      call onem(R)
      call onem(U)
      call onem(Uinv)
      
      
      ! Store the zero matrix in E
      !
      E = 0.d0
      
      
      ! Check if the determinant of F is greater than zero.
      !  If not, then print a diagnostic and cut back the 
      !  time increment.
      !
      call mdet(F,detF)
      if (detF.le.0.d0) then
        write(*,'(/5X,A/)') '--problem in kinematics-- the',
     +       ' determinant of F is not greater than 0'
        istat = 0
        return
      end if
      
      
      ! Calculate the right Cauchy-Green tensor C
      !
      C = matmul(transpose(F),F)
      
      
      ! Calculate the eigenvalues and eigenvectors of C
      !
      call spectral(C,omega,eigvec,istat)
      
      
      ! Calculate the principal values of U and E
      !
      Ueigval(1) = dsqrt(omega(1))
      Ueigval(2) = dsqrt(omega(2))
      Ueigval(3) = dsqrt(omega(3))
      !
      U(1,1) = Ueigval(1)
      U(2,2) = Ueigval(2)
      U(3,3) = Ueigval(3)
      !
      E(1,1) = dlog(Ueigval(1))
      E(2,2) = dlog(Ueigval(2))
      E(3,3) = dlog(Ueigval(3))
      
      
      ! Calculate the complete tensors U and E
      !
      U = matmul(matmul(eigvec,U),transpose(eigvec))
      E = matmul(matmul(eigvec,E),transpose(eigvec))
      
      
      ! Calculate Uinv
      !
      call matInv3D(U,Uinv,detF,istat)
      
      
      ! calculate R
      !
      R = matmul(F,Uinv)
      
      
      return
      end subroutine skinem
      
!***********************************************************************
!     The following subroutines calculate the spectral
!      decomposition of a symmetric 3 by 3 matrix
!***********************************************************************
      
      subroutine spectral(A,D,V,istat)
      !
      ! This subroutine calculates the eigenvalues and eigenvectors of
      !  a symmetric 3 by 3 matrix A.
      !
      ! The output consists of a vector D containing the three
      !  eigenvalues in ascending order, and a matrix V whose
      !  columns contain the corresponding eigenvectors.
      !
      implicit none
      !
      integer np,nrot,i,j,istat
      parameter(np=3)
      !
      real*8 D(3),V(3,3),A(3,3),E(3,3)
      
      
      E = A
      !
      call jacobi(E,3,np,D,V,nrot,istat)
      call eigsrt(D,V,3,np)
      
      
      return
      end subroutine spectral
      
!***********************************************************************
      
      subroutine jacobi(A,n,np,D,V,nrot,istat)
      !
      ! Computes all eigenvalues and eigenvectors of a real symmetric
      !  matrix A, which is of size n by n, stored in a physical
      !  np by np array.  On output, elements of A above the diagonal
      !  are destroyed, but the diagonal and sub-diagonal are unchanged
      !  and give full information about the original symmetric matrix.
      !  Vector D returns the eigenvalues of A in its first n elements.
      !  V is a matrix with the same logical and physical dimensions as
      !  A whose columns contain, upon output, the normalized
      !  eigenvectors of A.  nrot returns the number of Jacobi rotation
      !  which were required.
      !
      ! This subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer ip,iq,n,nmax,np,nrot,i,j,istat
      parameter (nmax=100)
      !
      real*8 A(np,np),D(np),V(np,np),B(nmax),Z(nmax),
     +  sm,tresh,G,T,H,theta,S,C,tau
      
      
      ! Initialize V to the identity matrix
      !
      call onem(V)
      
      
      ! Initialize B and D to the diagonal of A, and Z to zero.
      !  The vector Z will accumulate terms of the form T*A_PQ as
      !  in equation (11.1.14)
      !
      do ip = 1,n
        B(ip) = A(ip,ip)
        D(ip) = B(ip)
        Z(ip) = 0.d0
      end do
      
      
      ! Begin iteration
      !
      nrot = 0
      do i=1,50
        !
        ! Sum off-diagonal elements
        !
        sm = 0.d0
        do ip=1,n-1
          do iq=ip+1,n
            sm = sm + dabs(A(ip,iq))
          end do
        end do
        !
        ! If sm = 0., then return.  This is the normal return,
        !  which relies on quadratic convergence to machine
        !  underflow.
        !
        if (sm.eq.0.d0) return
        !
        ! In the first three sweeps carry out the PQ rotation only if
        !  |A_PQ| > tresh, where tresh is some threshold value,
        !  see equation (11.1.25).  Thereafter tresh = 0.
        !
        if (i.lt.4) then
          tresh = 0.2d0*sm/n**2
        else
          tresh = 0.d0
        end if
        !
        do ip=1,n-1
          do iq=ip+1,n
            G = 100.d0*dabs(A(ip,iq))
            !
            ! After four sweeps, skip the rotation if the 
            !  off-diagonal element is small.
            !
            if ((i.gt.4).and.(dabs(D(ip))+G.eq.dabs(D(ip)))
     +          .and.(dabs(D(iq))+G.eq.dabs(D(iq)))) then
              A(ip,iq) = 0.d0
            else if (dabs(A(ip,iq)).gt.tresh) then
              H = D(iq) - D(ip)
              if (dabs(H)+G.eq.dabs(H)) then
                !
                ! T = 1./(2.*theta), equation (11.1.10)
                !
                T =A(ip,iq)/H
              else
                theta = 0.5d0*H/A(ip,iq)
                T =1.d0/(dabs(theta)+dsqrt(1.d0+theta**2.d0))
                if (theta.lt.0.d0) T = -T
              end if
              C = 1.d0/dsqrt(1.d0 + T**2.d0)
              S = T*C
              tau = S/(1.d0 + C)
              H = T*A(ip,iq)
              Z(ip) = Z(ip) - H
              Z(iq) = Z(iq) + H
              D(ip) = D(ip) - H
              D(iq) = D(iq) + H
              A(ip,iq) = 0.d0
              !
              ! Case of rotations 1 <= J < P
              !
              do j=1,ip-1
                G = A(j,ip)
                H = A(j,iq)
                A(j,ip) = G - S*(H + G*tau)
                A(j,iq) = H + S*(G - H*tau)
              end do
              !
              ! Case of rotations P < J < Q
              !
              do j=ip+1,iq-1
                G = A(ip,j)
                H = A(j,iq)
                A(ip,j) = G - S*(H + G*tau)
                A(j,iq) = H + S*(G - H*tau)
              end do
              !
              ! Case of rotations Q < J <= N
              !
              do j=iq+1,n
                G = A(ip,j)
                H = A(iq,j)
                A(ip,j) = G - S*(H + G*tau)
                A(iq,j) = H + S*(G - H*tau)
              end do
              do j = 1,n
                G = V(j,ip)
                H = V(j,iq)
                V(j,ip) = G - S*(H + G*tau)
                V(j,iq) = H + S*(G - H*tau)
              end do
              nrot = nrot + 1
            end if
          end do
        end do
        !
        ! Update D with the sum of T*A_PQ, and reinitialize Z
        !
        do ip=1,n
          B(ip) = B(ip) + Z(ip)
          D(ip) = B(ip)
          Z(ip) = 0.d0
        end do
      end do
      
      
      ! If the algorithm has reached this stage, then there
      !  are too many sweeps.  Print a diagnostic and cut the 
      !  time increment.
      !
      write (*,'(/1X,A/)') '50 iterations in jacobi should never happen'
      istat = 0
      
      
      return
      end subroutine jacobi
      
!***********************************************************************
      
      subroutine eigsrt(D,V,n,np)
      !
      ! Given the eigenvalues D and eigenvectors V as output from
      !  jacobi, this subroutine sorts the eigenvales into ascending
      !  order and rearranges the colmns of V accordingly.
      !
      ! The subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer n,np,i,j,k
      !
      real*8 D(np),V(np,np),P
      
      
      do i=1,n-1
        k = i
        P = D(i)
        do j=i+1,n
          if (D(j).ge.P) then
            k = j
            P = D(j)
          end if
        end do
        if (k.ne.i) then
          D(k) = D(i)
          D(i) = P
          do j=1,n
            P = V(j,i)
            V(j,i) = V(j,k)
            V(j,k) = P
          end do
        end if
      end do
      
      
      return
      end subroutine eigsrt
      
!***********************************************************************
!     Utility subroutines
!***********************************************************************
      
      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv
      
      
      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: SUBROUTINE matInv3:'
        write(*,*) 'WARNING: DET of MAT=',DET_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      
      
      return
      end subroutine matInv3D
      
!***********************************************************************
      
      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det
      
      
      det = A(1,1)*A(2,2)*A(3,3) 
     +  + A(1,2)*A(2,3)*A(3,1)
     +  + A(1,3)*A(2,1)*A(3,2)
     +  - A(3,1)*A(2,2)*A(1,3)
     +  - A(3,2)*A(2,3)*A(1,1)
     +  - A(3,3)*A(2,1)*A(1,2)
      
      
      return
      end subroutine mdet
      
!***********************************************************************
      
      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)
      
      
      do i=1,3
        do J=1,3
          if (i .eq. j) then
            A(i,j) = 1.0
          else
            A(i,j) = 0.0
          end if
        end do
      end do
      
      
      return
      end subroutine onem
      
!***********************************************************************
