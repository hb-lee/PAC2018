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
      eta0,eta1,rr ! scalar inner product results

   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: & ! 184 * 104 * 9
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

   do iblock=1,nblocks_tropic     ! nblocks_tropic = 8   local_block_num 
      this_block = get_block(blocks_tropic(iblock),iblock)    ! all_blocks(block_id)    , size(blocks_tropic) = 8/7  0~27:8, 28~55:7  1 <= blocks_tropic(iblock) <= 56

      call btrop_operator(S,X,this_block,iblock)
      R(:,:,iblock) = B(:,:,iblock) - S(:,:,iblock)  !R = B - AX ; block info ref as 3rd demision
      S(:,:,iblock) = c0    ! c0 = 0_r8
   end do ! block loop

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! initialize fields and scalars
!
!-----------------------------------------------------------------------
   ! bndy_tropic maxblocks_ew_snd = 3
   call update_ghost_cells(R, bndy_tropic, field_loc_center, &     !  constant 1
                                           field_type_scalar)   ! call boundary_2d_real  (boundary.f90)
   eta0 =c1
   solv_sum_iters = solv_max_iters    ! solv_max_iters = 1000

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

         if (lprecond) then     ! lprecond = .false.
            call preconditioner(WORK1,R,this_block,iblock)        ! preconditioner(PX,X,this_block,bid)
         else
            where (A0(:,:,iblock) /= c0)
               WORK1(:,:,iblock) = R(:,:,iblock)/A0(:,:,iblock)
            elsewhere
               WORK1(:,:,iblock) = c0
            endwhere
         endif

         WORK0(:,:,iblock) = R(:,:,iblock)*WORK1(:,:,iblock)     ! Work0 = (R * R) / A0;  Work1 = R/A0;  R = B - AX
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! update conjugate direction vector s
!
!-----------------------------------------------------------------------

      if (lprecond) &                ! lprecond = .false.
         call update_ghost_cells(WORK1,bndy_tropic, field_loc_center,&
                                                    field_type_scalar)             
      ! distrb_tropic%local_block_num = 7 / 8
      !*** (r,(PC)r)
      eta1 = global_sum(WORK0, distrb_tropic, field_loc_center, RCALCT_B)    ! global_sum_real (global_reductions.f90)   RCALCT_B 乘法mask  eta1 = global_sum_real;  local_sum + WORK0(i,j,bid)*MASK(i,j,bid) ==> ALLREDUCE

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)

         S(:,:,iblock) = WORK1(:,:,iblock) + S(:,:,iblock)*(eta1/eta0)           ! S = c0 for first iter; Wrok1 = R/A0 ==> S; eta1/eta0 = 1 / global sum;

!-----------------------------------------------------------------------
!
! compute As
!
!-----------------------------------------------------------------------

         call btrop_operator(Q,S,this_block,iblock)      ! btrop_operator(AX,X,this_block,bid)
         WORK0(:,:,iblock) = Q(:,:,iblock)*S(:,:,iblock)     ! A_ * S * S; S = R/A0 + S * (eta/eta0); S = c0 for first iter

      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! compute next solution and residual
!
!-----------------------------------------------------------------------

! Work 0 used to calc global sum, AX used to update ghost cells;

      call update_ghost_cells(Q, bndy_tropic, field_loc_center, &
                                              field_type_scalar)

      eta0 = eta1   ! next iter
      eta1 = eta0/global_sum(WORK0, distrb_tropic, &
                             field_loc_center, RCALCT_B)     ! eta1 = 1 / (global sum) ^ iter time;

      !$OMP PARALLEL DO PRIVATE(iblock,this_block)

      do iblock=1,nblocks_tropic
         this_block = get_block(blocks_tropic(iblock),iblock)

         X(:,:,iblock) = X(:,:,iblock) + eta1*S(:,:,iblock)        ! update the X , init = X0;
         R(:,:,iblock) = R(:,:,iblock) - eta1*Q(:,:,iblock)        ! update the R , init = B - AX;

         if (mod(m,solv_ncheck) == 0) then

            call btrop_operator(R,X,this_block,iblock)
            R(:,:,iblock) = B(:,:,iblock) - R(:,:,iblock)
            WORK0(:,:,iblock) = R(:,:,iblock)*R(:,:,iblock)
         endif
      end do ! block loop

      !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! test for convergence
!
!-----------------------------------------------------------------------

      if (mod(m,solv_ncheck) == 0) then    ! check per 10 times, solv_ncheck    = 10

         call update_ghost_cells(R, bndy_tropic, field_loc_center,&
                                                 field_type_scalar)

         rr = global_sum(WORK0, distrb_tropic, &
                         field_loc_center, RCALCT_B) ! (r,r)

         if (rr < solv_convrg) then    ! solv_convrg    = 3121084.94389626
            ! ljm tuning
            if (my_task == master_task) &
               write(6,*)'pcg_iter_loop:iter#=',m
            solv_sum_iters = m
            exit iter_loop
         endif

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
      PX(i,j,bid) = PCNE(i,j,bid)*X(i+1,j+1,bid) + & !NE  precondNE
                    PCNW(i,j,bid)*X(i-1,j+1,bid) + & !NW
                    PCSE(i,j,bid)*X(i+1,j-1,bid) + & !SE
                    PCSW(i,j,bid)*X(i-1,j-1,bid) + & !SW
                    PCN (i,j,bid)*X(i ,j+1,bid) + & !N
                    PCS (i,j,bid)*X(i ,j-1,bid) + & !S
                    PCE (i,j,bid)*X(i+1,j ,bid) + & !E
                    PCW (i,j,bid)*X(i-1,j ,bid) + & !W
                    PCC (i,j,bid)*X(i ,j ,bid) !Center
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

   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie
      AX(i,j,bid) = A0 (i ,j ,bid)*X(i ,j ,bid) + & !center  btropWgtCenter
                    AN (i ,j ,bid)*X(i ,j+1,bid) + & !north
                    AN (i ,j-1,bid)*X(i ,j-1,bid) + & !north
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

 end subroutine btrop_operator

!***********************************************************************

end module solver_pcg_mod
