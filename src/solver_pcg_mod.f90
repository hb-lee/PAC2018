module solver_pcg_mod

   use kinds_mod
   use netcdf_mod
   use blocks
   use distribution
   use domain
   use domain_size
   use constants
   use boundary
   use global_reductions
   use gather_scatter
   use broadcast
   use grid
   use io
   use time_management
   use exit_mod
   use communicate, only: my_task, master_task
   implicit none
   save
   public :: pcg
   real (r8), dimension (nx_block,ny_block,max_blocks_clinic) :: &
      AC,                &! time-independent part of center 9pt weight
      A0_CLINIC           ! time-dependent center weight of 9pt operator
                          !   in baroclinic block distribution

   integer (int_kind) :: &
      solv_sum_iters      ! accumulated no of iterations (diagnostic)

   real (r8) ::  &
      rms_residual        ! residual (also a diagnostic)

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  other operator and preconditioner weights for barotropic operator
!
!-----------------------------------------------------------------------

   real (r8), dimension (nx_block,ny_block,max_blocks_tropic) :: &
      A0,AN,AE,ANE,         &! barotropic (9pt) operator coefficients
      RCALCT_B               ! land mask in barotropic distribution

   real (r8), dimension (:,:,:), allocatable :: &
      PCC,PCN,PCS,PCE,PCW,  &! preconditioner coefficients
      PCNE,PCSE,PCNW,PCSW

!-----------------------------------------------------------------------
!
!  scalar convergence-related variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      lprecond            ! true if computed preconditioner to be used

   real (r8) ::          &
      solv_convrg,       &! convergence error criterion
      sor,               &! for jacobi solver
      resid_norm          ! residual normalization

   integer (int_kind), parameter :: &
      solv_pcg = 1,      &! predefined solver types
      solv_cgr = 2,      &
      solv_jac = 3

   integer (int_kind) :: &
      solv_itype,        &! integer solver method (1=pcg, 2=cgr, 3=jac)
      solv_max_iters,    &! max number of iterations for solver
      solv_ncheck         ! check convergence every ncheck iterations

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      B                         ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      X,X0                  ! on input,  an initial guess for the solution
                         ! on output, solution of the linear system
   
contains

!***********************************************************************
!BOP
! !IROUTINE: pcg
! !INTERFACE:

 subroutine pcg(X,B)

! !DESCRIPTION:
! This routine uses a preconditioned conjugate-gradient solver to
! solve the equation $Ax=b$. Both the operator $A$ and preconditioner
! are nine-point stencils. If no preconditioner has been supplied,
! a diagonal preconditioner is applied. Convergence is checked
! every {\em ncheck} steps.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      B ! right hand side of linear system

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(inout) :: &
      X ! on input, an initial guess for the solution
                         ! on output, solution of the linear system

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      m, &! local iteration counter
      iblock ! local block counter

   real (r8) :: &
      eta0,eta1,rr,rr_old ! scalar inner product results

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
      R, &! residual (b-Ax)
      S, &! conjugate direction vector
      Q,WORK0,WORK1 ! various cg intermediate results

   character (char_len) :: &
      noconvrg ! error message for no convergence

   type (block) :: &
      this_block ! block information for current block

!-----------------------------------------------------------------------
!
! compute initial residual and initialize S
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,this_block)

   do iblock=1,nblocks_tropic
      this_block = get_block(blocks_tropic(iblock),iblock)

      call btrop_operator_opt(S,X,this_block,iblock)
      R(:,:,iblock) = B(:,:,iblock) - S(:,:,iblock)
          WORK0(:,:,iblock) = R(:,:,iblock)*R(:,:,iblock)
      S(:,:,iblock) = c0
          Q(:,:,iblock) = c0
   end do ! block loop

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! initialize fields and scalars
!
!-----------------------------------------------------------------------

   call update_ghost_cells(R, bndy_tropic, field_loc_center, &
                                           field_type_scalar)
   rr = global_sum(WORK0, distrb_tropic, field_loc_center, RCALCT_B)
   rr_old = rr
   eta0 =c1
   solv_sum_iters = solv_max_iters

!-----------------------------------------------------------------------
!
! iterate
!
!-----------------------------------------------------------------------

   iter_loop: do m = 1, solv_max_iters

!-----------------------------------------------------------------------
!
! calculate (PC)r
! diagonal preconditioner if preconditioner not specified
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)
         if (lprecond) then
            call preconditioner(WORK1,R,this_block,iblock)
         else
            call btrop_operator_opt(WORK1, R, this_block, iblock)
         endif
            WORK0(:,:,iblock) = R(:,:,iblock)*WORK1(:,:,iblock)
      end do ! block loop
      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! update conjugate direction vector s
!
!-----------------------------------------------------------------------

      call update_ghost_cells(WORK1,bndy_tropic, field_loc_center,&
                                                    field_type_scalar)
      !*** (r,(PC)r)
      eta1 = global_sum(WORK0, distrb_tropic, field_loc_center, RCALCT_B)

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)
		 S(:,:,iblock) =  R(:,:,iblock) + S(:,:,iblock)*(eta1/eta0)
!-----------------------------------------------------------------------
!
! compute As
!
!-----------------------------------------------------------------------

         Q(:,:,iblock) = WORK1(:,:,iblock) + Q(:,:,iblock)*(eta1/eta0)
         WORK0(:,:,iblock) = Q(:,:,iblock)*Q(:,:,iblock)
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! compute next solution and residual
!
!-----------------------------------------------------------------------
      if (lprecond) &
         call update_ghost_cells(Q, bndy_tropic, field_loc_center, &
                                              field_type_scalar)

      eta0 = eta1
      eta1 = eta0/global_sum(WORK0, distrb_tropic, &
                             field_loc_center, RCALCT_B)

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)

         X(:,:,iblock) = X(:,:,iblock) + eta1*S(:,:,iblock)

         if (mod(m,solv_ncheck) == 0) then

            call btrop_operator_opt(R,X,this_block,iblock)
            R(:,:,iblock) = B(:,:,iblock) - R(:,:,iblock)
            WORK0(:,:,iblock) = R(:,:,iblock)*R(:,:,iblock)
         else
            R(:,:,iblock) = R(:,:,iblock) - eta1*Q(:,:,iblock)
         endif
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! test for convergence
!
!-----------------------------------------------------------------------

        if (mod(m,solv_ncheck) == 0) then

         call update_ghost_cells(R, bndy_tropic, field_loc_center,&
                                                 field_type_scalar)

         rr = global_sum(WORK0, distrb_tropic, &
                         field_loc_center, RCALCT_B) ! (r,r)
		 if (rr > rr_old) then
			solv_sum_iters = m
			exit iter_loop
		 endif
		  rr_old = rr
		 else
		  rr = rr - eta1**2*eta0
		 endif

         if (rr < solv_convrg) then
            ! ljm tuning
            if (my_task == master_task) &
               write(6,*)'pcg_iter_loop:iter#=',m
            solv_sum_iters = m
            exit iter_loop
         endif

   enddo iter_loop

   rms_residual = sqrt(rr*resid_norm)

   if (solv_sum_iters == solv_max_iters) then
      if (solv_convrg /= c0) then
         write(noconvrg,'(a45,i11)') &
           'Barotropic solver not converged at time step ', nsteps_total
         call exit_POP(sigAbort,noconvrg)
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine pcg

!***********************************************************************
!BOP
! !IROUTINE: preconditioner
! !INTERFACE:

 subroutine preconditioner(PX,X,this_block,bid)

! !DESCRIPTION:
! This function applies a precomputed preconditioner as a nine-point
! stencil operator.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      X ! array to be operated on

   type (block), intent(in) :: &
      this_block ! block info for this block

   integer (int_kind), intent(in) :: &
      bid ! local block address for this block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(out) :: &
      PX ! nine point operator result

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j ! dummy counters

!-----------------------------------------------------------------------

   PX(:,:,bid) = c0

   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      PX(i,j,bid) = PCNE(i,j,bid)*X(i+1,j+1,bid) + &
                    PCNW(i,j,bid)*X(i-1,j+1,bid) + &
                    PCSE(i,j,bid)*X(i+1,j-1,bid) + &
                    PCSW(i,j,bid)*X(i-1,j-1,bid) + &
                    PCN (i,j,bid)*X(i ,j+1,bid) + &
                    PCS (i,j,bid)*X(i ,j-1,bid) + &
                    PCE (i,j,bid)*X(i+1,j ,bid) + &
                    PCW (i,j,bid)*X(i-1,j ,bid) + &
                    PCC (i,j,bid)*X(i ,j ,bid)
   end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine preconditioner

!***********************************************************************
!BOP
! !IROUTINE: btrop_operator
! !INTERFACE:

 subroutine btrop_operator(AX,X,this_block,bid)

! !DESCRIPTION:
! This routine applies the nine-point stencil operator for the
! barotropic solver. It takes advantage of some 9pt weights being
! shifted versions of others.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      X ! array to be operated on

   type (block), intent(in) :: &
      this_block ! block info for this block

   integer (int_kind), intent(in) :: &
      bid ! local block address for this block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(out) :: &
      AX ! nine point operator result (Ax)

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j ! dummy counters

!-----------------------------------------------------------------------

   AX(:,:,bid) = c0
   !$OMP PARALLEL DO PRIVATE(i,j)
   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      AX(i,j,bid) = A0 (i ,j ,bid)*X(i ,j ,bid) + &
                    AN (i ,j ,bid)*X(i ,j+1,bid) + &
                    AN (i ,j-1,bid)*X(i ,j-1,bid) + &
                    AE (i ,j ,bid)*X(i+1,j ,bid) + &
                    AE (i-1,j ,bid)*X(i-1,j ,bid) + &
                    ANE(i ,j ,bid)*X(i+1,j+1,bid) + &
                    ANE(i ,j-1,bid)*X(i+1,j-1,bid) + &
                    ANE(i-1,j ,bid)*X(i-1,j+1,bid) + &
                    ANE(i-1,j-1,bid)*X(i-1,j-1,bid)
   end do
   end do
   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine btrop_operator
 
 
 subroutine btrop_operator_opt(AX,X,this_block,bid)

! !DESCRIPTION:
! This routine applies the nine-point stencil operator for the
! barotropic solver. It takes advantage of some 9pt weights being
! shifted versions of others.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(in) :: &
      X ! array to be operated on

   type (block), intent(in) :: &
      this_block ! block info for this block

   integer (int_kind), intent(in) :: &
      bid ! local block address for this block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
      intent(out) :: &
      AX ! nine point operator result (Ax)

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i,j,limit ! dummy counters

!-----------------------------------------------------------------------

    AX(:,:,bid) = c0
   limit = this_block%ie - 3
   do j=this_block%jb,this_block%je
    do i=this_block%ib,limit,4
      AX(i,j,bid) = A0 (i ,j ,bid)*X(i ,j ,bid) + & !center  btropWgtCenter
                    AN (i ,j ,bid)*X(i ,j+1,bid) + & !north
                    AN (i ,j-1,bid)*X(i ,j-1,bid) + & !northeast
                    AE (i ,j ,bid)*X(i+1,j ,bid) + & !east
                    AE (i-1,j ,bid)*X(i-1,j ,bid) + & !east
                    ANE(i ,j ,bid)*X(i+1,j+1,bid) + & !NE
                    ANE(i ,j-1,bid)*X(i+1,j-1,bid) + & !NE
                    ANE(i-1,j ,bid)*X(i-1,j+1,bid) + & !NE
                    ANE(i-1,j-1,bid)*X(i-1,j-1,bid) !NE
      AX(i+1,j,bid) = A0 (i+1 ,j ,bid)*X(i+1 ,j ,bid) + & !center  btropWgtCenter
                    AN (i+1 ,j ,bid)*X(i+1 ,j+1,bid) + & !north
                    AN (i+1 ,j-1,bid)*X(i+1 ,j-1,bid) + & !northeast
                    AE (i+1 ,j ,bid)*X(i+2,j ,bid) + & !east
                    AE (i,j ,bid)*X(i,j ,bid) + & !east
                    ANE(i+1 ,j ,bid)*X(i+2,j+1,bid) + & !NE
                    ANE(i+1 ,j-1,bid)*X(i+2,j-1,bid) + & !NE
                    ANE(i,j ,bid)*X(i,j+1,bid) + & !NE
                    ANE(i,j-1,bid)*X(i,j-1,bid) !NE    
      AX(i+2,j,bid) = A0 (i+2 ,j ,bid)*X(i+2 ,j ,bid) + & !center  btropWgtCenter
                    AN (i+2 ,j ,bid)*X(i+2 ,j+1,bid) + & !north
                    AN (i+2 ,j-1,bid)*X(i+2 ,j-1,bid) + & !northeast
                    AE (i+2 ,j ,bid)*X(i+3,j ,bid) + & !east
                    AE (i+1,j ,bid)*X(i+1,j ,bid) + & !east
                    ANE(i+2 ,j ,bid)*X(i+3,j+1,bid) + & !NE
                    ANE(i+2 ,j-1,bid)*X(i+3,j-1,bid) + & !NE
                    ANE(i+1,j ,bid)*X(i+1,j+1,bid) + & !NE
                    ANE(i+1,j-1,bid)*X(i+1,j-1,bid) !NE
      AX(i+3,j,bid) = A0 (i+3 ,j ,bid)*X(i+3 ,j ,bid) + & !center  btropWgtCenter
                    AN (i+3 ,j ,bid)*X(i+3 ,j+1,bid) + & !north
                    AN (i+3 ,j-1,bid)*X(i+3 ,j-1,bid) + & !northeast
                    AE (i+3 ,j ,bid)*X(i+4,j ,bid) + & !east
                    AE (i+2,j ,bid)*X(i+2,j ,bid) + & !east
                    ANE(i+3 ,j ,bid)*X(i+4,j+1,bid) + & !NE
                    ANE(i+3 ,j-1,bid)*X(i+4,j-1,bid) + & !NE
                    ANE(i+2,j ,bid)*X(i+2,j+1,bid) + & !NE
                    ANE(i+2,j-1,bid)*X(i+2,j-1,bid) !NE
    end do       
    do i=i,this_block%ie
      AX(i,j,bid) = A0 (i ,j ,bid)*X(i ,j ,bid) + & !center  btropWgtCenter
                    AN (i ,j ,bid)*X(i ,j+1,bid) + & !north
                    AN (i ,j-1,bid)*X(i ,j-1,bid) + & !northeast
                    AE (i ,j ,bid)*X(i+1,j ,bid) + & !east
                    AE (i-1,j ,bid)*X(i-1,j ,bid) + & !east
                    ANE(i ,j ,bid)*X(i+1,j+1,bid) + & !NE
                    ANE(i ,j-1,bid)*X(i+1,j-1,bid) + & !NE
                    ANE(i-1,j ,bid)*X(i-1,j+1,bid) + & !NE
                    ANE(i-1,j-1,bid)*X(i-1,j-1,bid) !NE
    end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine btrop_operator_opt

!***********************************************************************

end module solver_pcg_mod
