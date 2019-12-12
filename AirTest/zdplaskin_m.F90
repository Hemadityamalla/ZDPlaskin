!
! ZDPLASKIN version 2.0a
! (c) 2008, Sergey Pancheshnyi (pancheshnyi@gmail.com)
!
! BOLSIG+
! (c) 2005, Gerjan Hagelaar (gerjan.hagelaar@laplace.univ-tlse.fr)
!
! http://www.zdplaskin.laplace.univ-tlse.fr/
! This software is provided "as is" without warranty and non-commercial use is freely
! granted provided proper reference is made in publications resulting from its use.
! Use of ZDPlasKin in commerical software requires a license.
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Tue Dec 10 15:41:28 2019
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ELEMENTS: E N O 
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE ZDPlasKin
!
!-----------------------------------------------------------------------------------------------------------------------------------
module ZDPlasKin
  use dvode_f90_m, only : vode_opts
  implicit none
  public
!
! config
!
  integer, parameter :: species_max = 44, species_electrons = 44, species_length = 9, reactions_max = 425, reactions_length = 40
  double precision                          :: density(species_max)
  integer                                   :: species_charge(species_max)
  character(species_length)                 :: species_name(species_max)
  character(reactions_length)               :: reaction_sign(reactions_max)
  logical                                   :: lreaction_block(reactions_max)
!
! internal config
!
  double precision, parameter, private      :: vode_atol = 1.00D-10, vode_rtol = 1.00D-05, cfg_rtol = 1.00D-01
  integer, parameter, private               :: vode_neq = species_max + 1
  type (vode_opts), private                 :: vode_options
  integer, private                          :: vode_itask, vode_istate, ifile_unit = 5
  double precision, private                 :: stat_dens(species_max), stat_src(species_max), stat_rrt(reactions_max), stat_time, &
                                               dens_loc(vode_neq,0:3), rrt_loc(reactions_max), tsav = -huge(tsav), &
                                               mach_accur, mach_tiny
  double precision                          :: rrt(reactions_max), mrtm(species_max, reactions_max), ZDPlasKin_cfg(14)
  logical, private                          :: lZDPlasKin_init = .false., lprint, lstat_accum, ldensity_constant, &
                                               density_constant(species_max), lgas_heating
!
! qtplaskin config
!
  logical, private                          :: lqtplaskin, lqtplaskin_first = .true.
  double precision, parameter, private      :: qtplaskin_atol = 1.00D+00, qtplaskin_rtol = 1.00D-02
  character(32), allocatable                :: qtplaskin_user_names(:)
  double precision, allocatable             :: qtplaskin_user_data(:)
!
! physical constants
!
  double precision, parameter, private      :: eV_to_K = 1.16045052d4, q_elem = 1.60217662d-19, k_B = 1.38064852d-23
!
! bolsig+ config
!
  double precision, parameter, private      :: bolsig_rtol = 1.00D-03, bolsig_rtol_half = 3.16D-02, &
                                               bolsig_field_min = 1.00D-01, bolsig_field_max = 1.00D+03, &
                                               bolsig_eecol_frac_def = 1.00D-05
  double precision, private                 :: bolsig_eecol_frac
  integer, parameter, private               :: bolsig_species_max = 20, bolsig_species_length = 6, bolsig_rates_max = 57, &
                                               bolsig_addsect_max = 12 
  character(*), parameter, private          :: bolsigfile = "bolsigdb.dat"
  integer                                   :: bolsig_pointer(bolsig_rates_max) = -1
  integer, private                          :: bolsig_species_index(bolsig_species_max) = -1, bolsig_collisions_max = 0, &
                                               bolsig_addsect(2,bolsig_addsect_max) 
  logical, private                          :: lbolsig_ignore_gas_temp, lbolsig_Maxwell_EEDF
  double precision, allocatable             :: bolsig_rates(:)
  character(bolsig_species_length), private :: bolsig_species(bolsig_species_max)
  interface
    subroutine ZDPlasKin_bolsig_Init(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_Init
    subroutine ZDPlasKin_bolsig_ReadCollisions(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_ReadCollisions
    subroutine ZDPlasKin_bolsig_GetCollisions(i,j)
      integer, intent(out) :: i, j
    end subroutine ZDPlasKin_bolsig_GetCollisions
    subroutine ZDPlasKin_bolsig_GetSpeciesName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetSpeciesName
    subroutine ZDPlasKin_bolsig_GetReactionName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetReactionName
    subroutine ZDPlasKin_bolsig_SolveBoltzmann(i,a,j,b)
      integer, intent(in) :: i, j
      double precision, intent(in)  :: a(i)
      double precision, intent(out) :: b(j)
    end subroutine ZDPlasKin_bolsig_SolveBoltzmann
    subroutine ZDPlasKin_bolsig_GetEEDF(i,a,b)
      integer, intent(in) :: i
      double precision, intent(out) :: a,b
    end subroutine ZDPlasKin_bolsig_GetEEDF
  end interface
!
! data section
!
  data species_charge(1:species_max) &
  / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1, 0, 1,-1,&
    1,-1/
  data species_name(1:species_max) &
  /"N2       ","N2(V1)   ","N2(V2)   ","N2(V3)   ","N2(V4)   ","N2(V5)   ","N2(V6)   ","N2(V7)   ","N2(V8)   ","N2(A3)   ",&
   "N2(B3)   ","N2(A`1)  ","N2(C3)   ","N        ","N(2D)    ","N(2P)    ","N^+      ","N2^+     ","N3^+     ","N4^+     ",&
   "O2       ","O2(V1)   ","O2(V2)   ","O2(V3)   ","O2(V4)   ","O2(A1)   ","O2(B1)   ","O2(4.5EV)","O        ","O(1D)    ",&
   "O(1S)    ","O3       ","O^+      ","O2^+     ","O4^+     ","O^-      ","O2^-     ","O3^-     ","O4^-     ","NO       ",&
   "NO^+     ","NO^-     ","O2^+N2   ","E        "/
  data reaction_sign(1:36) &
  /"bolsig:N2->N2(V1RES)                    ","bolsig:N2->N2(V1)                       ",&
   "bolsig:N2->N2(V2)                       ","bolsig:N2->N2(V3)                       ",&
   "bolsig:N2->N2(V4)                       ","bolsig:N2->N2(V5)                       ",&
   "bolsig:N2->N2(V6)                       ","bolsig:N2->N2(V7)                       ",&
   "bolsig:N2->N2(V8)                       ","bolsig:N2(V1)->N2                       ",&
   "bolsig:N2(V2)->N2                       ","bolsig:N2(V3)->N2                       ",&
   "bolsig:N2(V4)->N2                       ","bolsig:N2(V5)->N2                       ",&
   "bolsig:N2(V6)->N2                       ","bolsig:N2(V7)->N2                       ",&
   "bolsig:N2(V8)->N2                       ","bolsig:O2->O2(V1RES)                    ",&
   "bolsig:O2->O2(V1)                       ","bolsig:O2->O2(V2RES)                    ",&
   "bolsig:O2->O2(V2)                       ","bolsig:O2->O2(V3)                       ",&
   "bolsig:O2->O2(V4)                       ","bolsig:O2(V1)->O2                       ",&
   "bolsig:O2(V2)->O2                       ","bolsig:O2(V3)->O2                       ",&
   "bolsig:O2(V4)->O2                       ","N2(V1)+N2=>N2+N2                        ",&
   "N2(V2)+N2=>N2(V1)+N2                    ","N2(V3)+N2=>N2(V2)+N2                    ",&
   "N2(V4)+N2=>N2(V3)+N2                    ","N2(V5)+N2=>N2(V4)+N2                    ",&
   "N2(V6)+N2=>N2(V5)+N2                    ","N2(V7)+N2=>N2(V6)+N2                    ",&
   "N2(V8)+N2=>N2(V7)+N2                    ","N2+N2=>N2(V1)+N2                        "/
  data reaction_sign(37:72) &
  /"N2(V1)+N2=>N2(V2)+N2                    ","N2(V2)+N2=>N2(V3)+N2                    ",&
   "N2(V3)+N2=>N2(V4)+N2                    ","N2(V4)+N2=>N2(V5)+N2                    ",&
   "N2(V5)+N2=>N2(V6)+N2                    ","N2(V6)+N2=>N2(V7)+N2                    ",&
   "N2(V7)+N2=>N2(V8)+N2                    ","N2(V1)+N=>N2+N                          ",&
   "N2(V2)+N=>N2(V1)+N                      ","N2(V3)+N=>N2(V2)+N                      ",&
   "N2(V4)+N=>N2(V3)+N                      ","N2(V5)+N=>N2(V4)+N                      ",&
   "N2(V6)+N=>N2(V5)+N                      ","N2(V7)+N=>N2(V6)+N                      ",&
   "N2(V8)+N=>N2(V7)+N                      ","N2+N=>N2(V1)+N                          ",&
   "N2(V1)+N=>N2(V2)+N                      ","N2(V2)+N=>N2(V3)+N                      ",&
   "N2(V3)+N=>N2(V4)+N                      ","N2(V4)+N=>N2(V5)+N                      ",&
   "N2(V5)+N=>N2(V6)+N                      ","N2(V6)+N=>N2(V7)+N                      ",&
   "N2(V7)+N=>N2(V8)+N                      ","N2(V1)+O=>N2+O                          ",&
   "N2(V2)+O=>N2(V1)+O                      ","N2(V3)+O=>N2(V2)+O                      ",&
   "N2(V4)+O=>N2(V3)+O                      ","N2(V5)+O=>N2(V4)+O                      ",&
   "N2(V6)+O=>N2(V5)+O                      ","N2(V7)+O=>N2(V6)+O                      ",&
   "N2(V8)+O=>N2(V7)+O                      ","N2+O=>N2(V1)+O                          ",&
   "N2(V1)+O=>N2(V2)+O                      ","N2(V2)+O=>N2(V3)+O                      ",&
   "N2(V3)+O=>N2(V4)+O                      ","N2(V4)+O=>N2(V5)+O                      "/
  data reaction_sign(73:108) &
  /"N2(V5)+O=>N2(V6)+O                      ","N2(V6)+O=>N2(V7)+O                      ",&
   "N2(V7)+O=>N2(V8)+O                      ","O2(V1)+O2=>O2+O2                        ",&
   "O2(V2)+O2=>O2(V1)+O2                    ","O2(V3)+O2=>O2(V2)+O2                    ",&
   "O2(V4)+O2=>O2(V3)+O2                    ","O2+O2=>O2(V1)+O2                        ",&
   "O2(V1)+O2=>O2(V2)+O2                    ","O2(V2)+O2=>O2(V3)+O2                    ",&
   "O2(V3)+O2=>O2(V4)+O2                    ","O2(V1)+O=>O2+O                          ",&
   "O2(V2)+O=>O2(V1)+O                      ","O2(V3)+O=>O2(V2)+O                      ",&
   "O2(V4)+O=>O2(V3)+O                      ","O2+O=>O2(V1)+O                          ",&
   "O2(V1)+O=>O2(V2)+O                      ","O2(V2)+O=>O2(V3)+O                      ",&
   "O2(V3)+O=>O2(V4)+O                      ","bolsig:N2->N2(A3,V5-9)                  ",&
   "bolsig:N2->N2(A3,V10-)                  ","bolsig:N2->N2(B3)                       ",&
   "bolsig:N2->N2(W3)                       ","bolsig:N2->N2(B3)                      ",&
   "bolsig:N2->N2(A1)                      ","bolsig:N2->N2(A1)                       ",&
   "bolsig:N2->N2(W1)                       ","bolsig:N2->N2(C3)                       ",&
   "bolsig:N2->N2(E3)                       ","bolsig:N2->N2(A1)                     ",&
   "bolsig:N2->N2(SUM)                      ","bolsig:N2->N2^+                         ",&
   "bolsig:N2(A3)->N2^+                     ","bolsig:N->N^+                           ",&
   "bolsig:O2->O2(A1)                       ","bolsig:O2->O2(B1)                       "/
  data reaction_sign(109:144) &
  /"bolsig:O2->O2(4.5EV)                    ","bolsig:O2->O2(6.0EV)                    ",&
   "bolsig:O2->O2(8.4EV)                    ","bolsig:O2->O2(9.97EV)                   ",&
   "bolsig:O2->O2^+                         ","bolsig:O2->O^-+O                        ",&
   "bolsig:O2(A1)->O2                       ","bolsig:O2(A1)->O+O                      ",&
   "bolsig:O2(A1)->O2^+                     ","bolsig:O->O(1D)                         ",&
   "bolsig:O->O(1S)                         ","bolsig:O->O^+                           ",&
   "bolsig:NO->NO^+                         ","E+N2^+=>N+N                             ",&
   "E+N2^+=>N+N(2D)                         ","E+N2^+=>N+N(2P)                         ",&
   "E+N3^+=>N2+N                            ","E+N4^+=>N2+N2                           ",&
   "E+O2^+=>O+O                             ","E+O2^+=>O+O(1D)                         ",&
   "E+O2^+=>O+O(1S)                         ","E+O4^+=>O2+O2                           ",&
   "E+NO^+=>O+N                             ","E+NO^+=>O+N(2D)                         ",&
   "E+O2^+N2=>O2+N2                         ","E+N^++E=>N+E                            ",&
   "E+O^++E=>O+E                            ","E+N^++ANY_NEUTRAL=>N+ANY_NEUTRAL        ",&
   "E+O^++ANY_NEUTRAL=>O+ANY_NEUTRAL        ","E+O3=>O2^-+O                            ",&
   "E+O3=>O^-+O2                            ","E+O+O2=>O^-+O2                          ",&
   "E+O+O2=>O2^-+O                          ","E+O3+O2=>O3^-+O2                        ",&
   "E+O2+N2=>O2^-+N2                        ","E+NO+ANY_NEUTRAL=>NO^-+ANY_NEUTRAL      "/
  data reaction_sign(145:180) &
  /"O^-+O=>O2+E                             ","O^-+N=>NO+E                             ",&
   "O^-+O2=>O3+E                            ","O^-+O2(A1)=>O3+E                        ",&
   "O^-+O2(B1)=>O+O2+E                      ","O^-+N2(A3)=>O+N2+E                      ",&
   "O^-+N2(B3)=>O+N2+E                      ","O^-+O3=>O2+O2+E                         ",&
   "O2^-+O=>O3+E                            ","O2^-+O2=>O2+O2+E                        ",&
   "O2^-+O2(A1)=>O2+O2+E                    ","O2^-+O2(B1)=>O2+O2+E                    ",&
   "O2^-+N2=>O2+N2+E                        ","O2^-+N2(A3)=>O2+N2+E                    ",&
   "O2^-+N2(B3)=>O2+N2+E                    ","O3^-+O=>O2+O2+E                         ",&
   "N2(A3)=>N2                              ","N2(B3)=>N2(A3)                          ",&
   "N2(A`1)=>N2                             ","N2(C3)=>N2(B3)                          ",&
   "O2(A1)=>O2                              ","O2(B1)=>O2(A1)                          ",&
   "O2(B1)=>O2                              ","O2(4.5EV)=>O2                           ",&
   "N2(A3)+O=>NO+N(2D)                      ","N2(A3)+O=>N2+O(1S)                      ",&
   "N2(A3)+N=>N2+N                          ","N2(A3)+N=>N2+N(2P)                      ",&
   "N2(A3)+O2=>N2+O+O(1D)                   ","N2(A3)+O2=>N2+O2(A1)                    ",&
   "N2(A3)+O2=>N2+O2(B1)                    ","N2(A3)+N2=>N2+N2                        ",&
   "N2(A3)+NO=>N2+NO                        ","N2(A3)+N2(A3)=>N2+N2(B3)                ",&
   "N2(A3)+N2(A3)=>N2+N2(C3)                ","N2(B3)+N2=>N2(A3)+N2                    "/
  data reaction_sign(181:216) &
  /"N2(B3)+N2=>N2+N2                        ","N2(B3)+O2=>N2+O+O                       ",&
   "N2(B3)+NO=>N2(A3)+NO                    ","N2(C3)+N2=>N2(A`1)+N2                   ",&
   "N2(C3)+O2=>N2+O+O(1S)                   ","N2(A`1)+N2=>N2(B3)+N2                   ",&
   "N2(A`1)+O2=>N2+O+O(1D)                  ","N2(A`1)+NO=>N2+N+O                      ",&
   "N2(A`1)+N2(A3)=>N4^++E                  ","N2(A`1)+N2(A`1)=>N4^++E                 ",&
   "N+N+ANY_NEUTRAL=>N2(A3)+ANY_NEUTRAL     ","N+N+N=>N2(A3)+N                         ",&
   "N+N+O=>N2(A3)+O                         ","N+N+ANY_NEUTRAL=>N2(B3)+ANY_NEUTRAL     ",&
   "N+N+N=>N2(B3)+N                         ","N+N+O=>N2(B3)+O                         ",&
   "N(2D)+O=>N+O(1D)                        ","N(2D)+O2=>NO+O                          ",&
   "N(2D)+NO=>N2+O                          ","N(2D)+N2=>N+N2                          ",&
   "N(2P)+N=>N+N                            ","N(2P)+O=>N+O                            ",&
   "N(2P)+N=>N(2D)+N                        ","N(2P)+N2=>N+N2                          ",&
   "N(2P)+N(2D)=>N2^++E                     ","N(2P)+O2=>NO+O                          ",&
   "N(2P)+NO=>N2(A3)+O                      ","O2(A1)+O=>O2+O                          ",&
   "O2(A1)+N=>NO+O                          ","O2(A1)+O2=>O2+O2                        ",&
   "O2(A1)+N2=>O2+N2                        ","O2(A1)+NO=>O2+NO                        ",&
   "O2(A1)+O3=>O2+O2+O(1D)                  ","O2(A1)+O2(A1)=>O2+O2(B1)                ",&
   "O+O3=>O2+O2(A1)                         ","O2(B1)+O=>O2(A1)+O                      "/
  data reaction_sign(217:252) &
  /"O2(B1)+O=>O2+O(1D)                      ","O2(B1)+O2=>O2(A1)+O2                    ",&
   "O2(B1)+N2=>O2(A1)+N2                    ","O2(B1)+NO=>O2(A1)+NO                    ",&
   "O2(B1)+O3=>O2+O2+O                      ","O2(4.5EV)+O=>O2+O(1S)                   ",&
   "O2(4.5EV)+O2=>O2(B1)+O2(B1)             ","O2(4.5EV)+N2=>O2(B1)+N2                 ",&
   "O+O+ANY_NEUTRAL=>O2(A1)+ANY_NEUTRAL     ","O+O+O2=>O2(A1)+O2                       ",&
   "O+O+O=>O2(A1)+O                         ","O+O+ANY_NEUTRAL=>O2(B1)+ANY_NEUTRAL     ",&
   "O2(A1)+O2(A1)+O2=>O3+O3                 ","O(1D)+O=>O+O                            ",&
   "O(1D)+O2=>O+O2                          ","O(1D)+O2=>O+O2(A1)                      ",&
   "O(1D)+O2=>O+O2(B1)                      ","O(1D)+N2=>O+N2                          ",&
   "O(1D)+O3=>O2+O+O                        ","O(1D)+O3=>O2+O2                         ",&
   "O(1D)+NO=>O2+N                          ","O(1S)+O=>O(1D)+O                        ",&
   "O(1S)+N=>O+N                            ","O(1S)+O2=>O(1D)+O2                      ",&
   "O(1S)+O2=>O+O+O                         ","O(1S)+N2=>O+N2                          ",&
   "O(1S)+O2(A1)=>O+O+O                     ","O(1S)+O2(A1)=>O(1D)+O2(B1)              ",&
   "O(1S)+O2(A1)=>O+O+O                     ","O(1S)+NO=>O+NO                          ",&
   "O(1S)+NO=>O(1D)+NO                      ","O(1S)+O3=>O2+O2                         ",&
   "O(1S)+O3=>O2+O+O(1D)                    ","N+NO=>O+N2                              ",&
   "N+O2=>O+NO                              ","N+O3=>NO+O2                             "/
  data reaction_sign(253:288) &
  /"O+N2=>N+NO                              ","O+NO=>N+O2                              ",&
   "O+O3=>O2+O2(A1)                         ","NO+NO=>N2+O2                            ",&
   "O2+O2=>O+O3                             ","N+N=>N2^++E                             ",&
   "N+O=>NO^++E                             ","N2+ANY_NEUTRAL=>N+N+ANY_NEUTRAL         ",&
   "N2+N=>N+N+N                             ","N2+O=>N+N+O                             ",&
   "O2+ANY_NEUTRAL=>O+O+ANY_NEUTRAL         ","O2+O2=>O+O+O2                           ",&
   "O2+O=>O+O+O                             ","NO+ANY_NEUTRAL=>N+O+ANY_NEUTRAL         ",&
   "NO+N=>N+O+N                             ","NO+O=>N+O+O                             ",&
   "O3+N2=>O2+O+N2                          ","O3+O2=>O2+O+O2                          ",&
   "O3+N=>O2+O+N                            ","O3+O=>O2+O+O                            ",&
   "N+N+N2=>N2+N2                           ","N+N+O2=>N2+O2                           ",&
   "N+N+NO=>N2+NO                           ","N+N+N=>N2+N                             ",&
   "N+N+O=>N2+O                             ","O+O+N2=>O2+N2                           ",&
   "O+O+O2=>O2+O2                           ","O+O+NO=>O2+NO                           ",&
   "O+O+N=>O2+N                             ","O+O+O=>O2+O                             ",&
   "N+O+N2=>NO+N2                           ","N+O+O2=>NO+O2                           ",&
   "N+O+NO=>NO+NO                           ","N+O+N=>NO+N                             ",&
   "N+O+O=>NO+O                             ","O+O2+N2=>O3+N2                          "/
  data reaction_sign(289:324) &
  /"O+O2+O2=>O3+O2                          ","O+O2+NO=>O3+NO                          ",&
   "O+O2+N=>O3+N                            ","O+O2+O=>O3+O                            ",&
   "N^++O=>N+O^+                            ","N^++O2=>O2^++N                          ",&
   "N^++O2=>NO^++O                          ","N^++O2=>O^++NO                          ",&
   "N^++O3=>NO^++O2                         ","N^++NO=>NO^++N                          ",&
   "N^++NO=>N2^++O                          ","N^++NO=>O^++N2                          ",&
   "O^++N2=>NO^++N                          ","O^++O2=>O2^++O                          ",&
   "O^++O3=>O2^++O2                         ","O^++NO=>NO^++O                          ",&
   "O^++NO=>O2^++N                          ","O^++N(2D)=>N^++O                        ",&
   "N2^++O2=>O2^++N2                        ","N2^++O=>NO^++N                          ",&
   "N2^++O=>O^++N2                          ","N2^++O3=>O2^++O+N2                      ",&
   "N2^++N=>N^++N2                          ","N2^++NO=>NO^++N2                        ",&
   "O2^++N2=>NO^++NO                        ","O2^++N=>NO^++O                          ",&
   "O2^++NO=>NO^++O2                        ","N3^++O2=>O2^++N+N2                      ",&
   "N3^++N=>N2^++N2                         ","N3^++NO=>NO^++N+N2                      ",&
   "N4^++N2=>N2^++N2+N2                     ","N4^++O2=>O2^++N2+N2                     ",&
   "N4^++O=>O^++N2+N2                       ","N4^++N=>N^++N2+N2                       ",&
   "N4^++NO=>NO^++N2+N2                     ","O4^++N2=>O2^+N2+O2                      "/
  data reaction_sign(325:360) &
  /"O4^++O2=>O2^++O2+O2                     ","O4^++O2(A1)=>O2^++O2+O2                 ",&
   "O4^++O2(B1)=>O2^++O2+O2                 ","O4^++O=>O2^++O3                         ",&
   "O4^++NO=>NO^++O2+O2                     ","O2^+N2+N2=>O2^++N2+N2                   ",&
   "O2^+N2+O2=>O4^++N2                      ","N^++N2+N2=>N3^++N2                      ",&
   "N^++O+ANY_NEUTRAL=>NO^++ANY_NEUTRAL     ","N^++N+ANY_NEUTRAL=>N2^++ANY_NEUTRAL     ",&
   "O^++N2+ANY_NEUTRAL=>NO^++N+ANY_NEUTRAL  ","O^++O+ANY_NEUTRAL=>O2^++ANY_NEUTRAL     ",&
   "O^++N+ANY_NEUTRAL=>NO^++ANY_NEUTRAL     ","N2^++N2+N2=>N4^++N2                     ",&
   "N2^++N+N2=>N3^++N2                      ","O2^++O2+O2=>O4^++O2                     ",&
   "O2^++N2+N2=>O2^+N2+N2                   ","O^-+O2(A1)=>O2^-+O                      ",&
   "O^-+O3=>O3^-+O                          ","O2^-+O=>O^-+O2                          ",&
   "O2^-+O3=>O3^-+O2                        ","O3^-+O=>O2^-+O2                         ",&
   "NO^-+O2=>O2^-+NO                        ","O4^-+ANY_NEUTRAL=>O2^-+O2+ANY_NEUTRAL   ",&
   "O4^-+O=>O3^-+O2                         ","O4^-+O=>O^-+O2+O2                       ",&
   "O4^-+O2(A1)=>O2^-+O2+O2                 ","O4^-+O2(B1)=>O2^-+O2+O2                 ",&
   "O^-+O2+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL    ","O2^-+O2+ANY_NEUTRAL=>O4^-+ANY_NEUTRAL   ",&
   "O^-+N^+=>O+N                            ","O^-+N2^+=>O+N2                          ",&
   "O^-+O^+=>O+O                            ","O^-+O2^+=>O+O2                          ",&
   "O^-+NO^+=>O+NO                          ","O2^-+N^+=>O2+N                          "/
  data reaction_sign(361:396) &
  /"O2^-+N2^+=>O2+N2                        ","O2^-+O^+=>O2+O                          ",&
   "O2^-+O2^+=>O2+O2                        ","O2^-+NO^+=>O2+NO                        ",&
   "O3^-+N^+=>O3+N                          ","O3^-+N2^+=>O3+N2                        ",&
   "O3^-+O^+=>O3+O                          ","O3^-+O2^+=>O3+O2                        ",&
   "O3^-+NO^+=>O3+NO                        ","NO^-+N^+=>NO+N                          ",&
   "NO^-+N2^+=>NO+N2                        ","NO^-+O^+=>NO+O                          ",&
   "NO^-+O2^+=>NO+O2                        ","NO^-+NO^+=>NO+NO                        ",&
   "O^-+N2^+=>O+N+N                         ","O^-+N3^+=>O+N+N2                        ",&
   "O^-+N4^+=>O+N2+N2                       ","O^-+O2^+=>O+O+O                         ",&
   "O^-+O4^+=>O+O2+O2                       ","O^-+NO^+=>O+N+O                         ",&
   "O^-+O2^+N2=>O+O2+N2                     ","O2^-+N2^+=>O2+N+N                       ",&
   "O2^-+N3^+=>O2+N+N2                      ","O2^-+N4^+=>O2+N2+N2                     ",&
   "O2^-+O2^+=>O2+O+O                       ","O2^-+O4^+=>O2+O2+O2                     ",&
   "O2^-+NO^+=>O2+N+O                       ","O2^-+O2^+N2=>O2+O2+N2                   ",&
   "O3^-+N2^+=>O3+N+N                       ","O3^-+N3^+=>O3+N+N2                      ",&
   "O3^-+N4^+=>O3+N2+N2                     ","O3^-+O2^+=>O3+O+O                       ",&
   "O3^-+O4^+=>O3+O2+O2                     ","O3^-+NO^+=>O3+N+O                       ",&
   "O3^-+O2^+N2=>O3+O2+N2                   ","NO^-+N2^+=>NO+N+N                       "/
  data reaction_sign(397:425) &
  /"NO^-+N3^+=>NO+N+N2                      ","NO^-+N4^+=>NO+N2+N2                     ",&
   "NO^-+O2^+=>NO+O+O                       ","NO^-+O4^+=>NO+O2+O2                     ",&
   "NO^-+NO^+=>NO+N+O                       ","NO^-+O2^+N2=>NO+O2+N2                   ",&
   "O4^-+N^+=>O2+O2+N                       ","O4^-+N2^+=>O2+O2+N2                     ",&
   "O4^-+N3^+=>O2+O2+N2+N                   ","O4^-+N4^+=>O2+O2+N2+N2                  ",&
   "O4^-+O^+=>O2+O2+O                       ","O4^-+O2^+=>O2+O2+O2                     ",&
   "O4^-+O4^+=>O2+O2+O2+O2                  ","O4^-+NO^+=>O2+O2+NO                     ",&
   "O4^-+O2^+N2=>O2+O2+O2+N2                ","O^-+N^++ANY_NEUTRAL=>O+N+ANY_NEUTRAL    ",&
   "O^-+N2^++ANY_NEUTRAL=>O+N2+ANY_NEUTRAL  ","O^-+O^++ANY_NEUTRAL=>O+O+ANY_NEUTRAL    ",&
   "O^-+O2^++ANY_NEUTRAL=>O+O2+ANY_NEUTRAL  ","O^-+NO^++ANY_NEUTRAL=>O+NO+ANY_NEUTRAL  ",&
   "O2^-+N^++ANY_NEUTRAL=>O2+N+ANY_NEUTRAL  ","O2^-+N2^++ANY_NEUTRAL=>O2+N2+ANY_NEUTRAL",&
   "O2^-+O^++ANY_NEUTRAL=>O2+O+ANY_NEUTRAL  ","O2^-+O2^++ANY_NEUTRAL=>O2+O2+ANY_NEUTRAL",&
   "O2^-+NO^++ANY_NEUTRAL=>O2+NO+ANY_NEUTRAL","O^-+N^++ANY_NEUTRAL=>NO+ANY_NEUTRAL     ",&
   "O^-+O^++ANY_NEUTRAL=>O2+ANY_NEUTRAL     ","O^-+O2^++ANY_NEUTRAL=>O3+ANY_NEUTRAL    ",&
   "O2^-+O^++ANY_NEUTRAL=>O3+ANY_NEUTRAL    "/
  data bolsig_species(1:bolsig_species_max) &
  /"N2    ","N2(V1)","N2(V2)","N2(V3)","N2(V4)","N2(V5)","N2(V6)","N2(V7)","N2(V8)","N2(A3)","O2    ","O2(V1)","O2(V2)","O2(V3)",&
   "O2(V4)","O2(A1)","N     ","O     ","O3    ","NO    "/
  data bolsig_addsect(1:2,1:bolsig_addsect_max) &
  / 1, 2, 1, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1, 8, 1, 9,11,12,11,13,11,14,11,15/
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!
! initialization
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_init()
  implicit none
  character(256) :: string
  integer :: i, j, k
  write(*,"(/,A)") "ZDPlasKin (version " // "2.0a" // ") INIT:"
  if( lZDPlasKin_init ) call ZDPlasKin_stop("   ERROR: the ZDPlasKin library has been initialized")
  write(string,*) species_max
  write(*,"(2x,A)")  "species        ... " // trim(adjustl(string))
  write(string,*) reactions_max
  write(*,"(2x,A)")  "reactions      ... " // trim(adjustl(string))
  if(species_max<=0 .or. reactions_max<=0) call ZDPlasKin_stop("   ERROR: wrong preconfig data")
  write(*,"(2x,A,$)") "BOLSIG+ loader ... " // trim(adjustl(bolsigfile)) // " : "
  call ZDPlasKin_bolsig_Init(bolsigfile)
  do i = 1, bolsig_species_max
    call ZDPlasKin_bolsig_ReadCollisions(trim(bolsig_species(i)))
    j = bolsig_collisions_max
    call ZDPlasKin_bolsig_GetCollisions(k,bolsig_collisions_max)
    if(bolsig_collisions_max <= j) then
      write(*,*)
      call ZDPlasKin_stop("ERROR: wrong file or missing data " // &
                        "(" // trim(adjustl(bolsigfile)) // ": <" // trim(bolsig_species(i)) // ">).")
    endif
  enddo
  if(bolsig_species_max /= k) then
    write(*,*)
    call ZDPlasKin_stop("ERROR: internal error in BOLSIG+ loader")
  endif
  write(string,*) bolsig_species_max
  write(*,"(A,$)") trim(adjustl(string)) // " species & "
  write(string,*) bolsig_collisions_max
  write(*,"(A)")   trim(adjustl(string)) // " collisions"
  write(*,"(2x,A,$)") "species  link  ... "
  j = 0
  do i = 1, bolsig_species_max
    j = j + 1
    k = 1
    do while(k<=species_max .and. bolsig_species_index(i)<=0)
      call ZDPlasKin_bolsig_GetSpeciesName(string,i)
      if(trim(species_name(k)) == trim(string)) then
        bolsig_species_index(i) = k
      else
        k = k + 1
      endif
    enddo
    if(bolsig_species_index(i) <= 0) call ZDPlasKin_stop("cannot find species link for <" // trim(string) // ">")
  enddo
  write(string,*) j
  write(*,"(A)") trim(adjustl(string))
  write(*,"(2x,A,$)") "process  link  ... "
  i = 1
  j = 1
  do while(i<=reactions_max .and. j<=bolsig_rates_max)
    if(reaction_sign(i)(1:7) == "bolsig:") then
      k = 1
      do while(k<=bolsig_collisions_max .and. bolsig_pointer(j)<=0)
        call ZDPlasKin_bolsig_GetReactionName(string,k)
        if(trim(string) == trim(reaction_sign(i)(8:))) then
          bolsig_pointer(j) = k
        else
          k = k + 1
        endif
      enddo
      if(bolsig_pointer(j) <= 0) call ZDPlasKin_stop("cannot find processes link for <" // trim(reaction_sign(i)) // ">")
      j = j + 1
    endif
    i = i + 1
  enddo
  if(j <= bolsig_rates_max) then
    call ZDPlasKin_stop("internal error")
  else
    write(string,*) bolsig_rates_max
    write(*,"(A)") trim(adjustl(string))
  endif
  i = 0
  do while((1.0d0+10.0d0**(i-1)) /= 1.0d0)
    i = i - 1
  enddo
  mach_accur = 10.0d0**i
  mach_tiny  = sqrt( tiny(mach_tiny) )
  lZDPlasKin_init = .true.
  call ZDPlasKin_reset()
  write(*,"(A,/)") "ZDPlasKin INIT DONE"
  return
end subroutine ZDPlasKin_init
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using implicit solver dvode_f90
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep(time,dtime)
  use dvode_f90_m, only : dvode_f90
  implicit none
  double precision, intent(in)    ::  time
  double precision, intent(inout) :: dtime
  double precision, save :: densav(vode_neq) = 0.0d0, cfgsav(3) = 0.0d0
  double precision :: tout
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(time < tsav) vode_istate = 1
  tsav = time
  if(dtime > 0.0d0) then
    vode_itask = 1
    tout = time + dtime
    if(dtime < mach_accur*abs(tout)) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: dtime parameter is too small (subroutine ZDPlasKin_timestep)")
  else
    vode_itask = 2
    tout = ( 1.0d0 + mach_accur ) * time + mach_tiny
  endif
  dens_loc(1:species_max,0) = density(:)
  dens_loc(1:species_max,1) = 0.5d0 * ( density(:) + abs( density(:) ) )
  if(any(dens_loc(1:species_max,1) /= densav(1:species_max))) vode_istate = 1
  densav(1:species_max) = dens_loc(1:species_max,1)
  if(vode_istate /= 1 .and. any( abs(cfgsav(:)-ZDPlasKin_cfg(1:3)) > cfg_rtol*abs(cfgsav(:)+ZDPlasKin_cfg(1:3)) )) vode_istate = 1
  cfgsav(:) = ZDPlasKin_cfg(1:3)
  if( lgas_heating ) then
    if(ZDPlasKin_cfg(1) /= densav(species_max+1)) vode_istate = 1
    densav(species_max+1) = ZDPlasKin_cfg(1)
  endif
  call dvode_f90(ZDPlasKin_fex,vode_neq,densav,tsav,tout,vode_itask,vode_istate,vode_options,j_fcn=ZDPlasKin_jex)
  if(vode_istate < 0) then
    write(*,"(A,1pd11.4)") "Tgas   =", ZDPlasKin_cfg(1)
    write(*,"(A,1pd11.4)") "    EN =", ZDPlasKin_cfg(3)
    write(*,"(A,1pd11.4)") "    Te =", ZDPlasKin_cfg(4)
    call ZDPlasKin_stop("ZDPlasKin ERROR: DVODE solver issued an error (subroutine ZDPlasKin_timestep)")
  endif
  if( lgas_heating ) ZDPlasKin_cfg(1) = densav(species_max+1)
  density(:) = dens_loc(1:species_max,0) - dens_loc(1:species_max,1) + densav(1:species_max)
  if(dtime <= 0.0d0) dtime = tsav - time
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using explicit Euler method
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep_explicit(time,dtime,rtol_loc,atol_loc,switch_implicit)
  implicit none
  double precision, intent(in) ::  time, rtol_loc, atol_loc
  double precision, intent(inout) :: dtime
  double precision, optional, intent(in) :: switch_implicit
  double precision :: time_loc, time_end, dtime_loc, dtime_max
  logical, save :: lwarn = .true.
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(rtol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: rtol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  if(atol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: atol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  tsav     = time
  time_loc = 0.0d0
  time_end = 0.5d0 * ( dtime + abs(dtime) ) + mach_tiny
  do while(time_loc < time_end)
    dens_loc(1:species_max,0) = density(:)
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    dens_loc(:,1) = 0.5d0 * ( dens_loc(:,0) + abs( dens_loc(:,0) ) )
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,1),dens_loc(:,2))
    where(dens_loc(:,2) >= 0.0d0)
      dens_loc(:,3) = + dens_loc(:,2) /    ( rtol_loc * dens_loc(:,1) + atol_loc ) 
    elsewhere
      dens_loc(:,3) = - dens_loc(:,2) / min( rtol_loc * dens_loc(:,1) + atol_loc , dens_loc(:,1) + mach_tiny )
    endwhere
    dtime_loc = 1.0d0 / ( maxval( dens_loc(:,3) ) + mach_tiny )
    if(dtime > 0.0d0) then
      dtime_max = dtime - time_loc
      dtime_loc = min( dtime_loc , dtime_max )
      if( present(switch_implicit) ) then
        if(dtime_loc*switch_implicit < dtime_max) then
          if(lprint .and. lwarn) then
            write(*,"(A,/,A,1pd9.2,A)") "ZDPlasKin INFO: low efficiency of Euler method (subroutine ZDPlasKin_timestep_explicit)", &
                        "                ZDPlasKin_timestep subroutine will be used in similar conditions (", switch_implicit, ")"
            lwarn = .false.
          endif
          time_loc = tsav
          density(:) = dens_loc(1:species_max,0)
          if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0)
          call ZDPlasKin_timestep(time_loc,dtime_max)
          return
        endif
      endif
    else
      dtime = dtime_loc
    endif
    time_loc = time_loc + dtime_loc
    tsav     = time     +  time_loc
    density(:) = dens_loc(1:species_max,0) + dtime_loc * dens_loc(1:species_max,2)
    if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0) + dtime_loc * dens_loc(species_max+1,2)
  enddo
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep_explicit
!-----------------------------------------------------------------------------------------------------------------------------------
!
! update BOLSIG+ solution and get electron parameters
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_bolsig_rates(lbolsig_force)
  implicit none
  logical, optional, intent(in) :: lbolsig_force
  logical :: lforce
  integer :: i, j, k
  integer, save :: bolsig_points_max
  logical, save :: lfirst = .true., leecol = .true.
  double precision :: error, density_loc, cfg_loc(6+bolsig_species_max)
  double precision, save :: low_density_limit = bolsig_rtol, bolsig_mesh_a, bolsig_mesh_b
  double precision, save, allocatable :: bolsig_cfg(:,:), bolsig_reslt(:,:)
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    bolsig_mesh_a = 1.0d0 / log( 1.0d0 + bolsig_rtol )
    bolsig_mesh_b = bolsig_mesh_a * log( bolsig_field_min ) - 0.5d0
    bolsig_points_max = int( bolsig_mesh_a * log( bolsig_field_max ) - bolsig_mesh_b )
    allocate(bolsig_rates(bolsig_collisions_max), &
             bolsig_cfg(6+bolsig_species_max,0:bolsig_points_max), &
             bolsig_reslt(10+bolsig_collisions_max,0:bolsig_points_max),stat=i)
    if(i /= 0) call ZDPlasKin_stop("ZDPlasKin ERROR: memory allocation error (subroutine ZDPlasKin_bolsig_rates)")
    bolsig_cfg(:,:) = 0.0d0
    lfirst = .false.
  endif
  if( present(lbolsig_force) ) then
    lforce = lbolsig_force
  else
    lforce = .false.
  endif
  if(ZDPlasKin_cfg(1) <= 0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
  if(.not. lbolsig_Maxwell_EEDF ) then
    if(ZDPlasKin_cfg(2) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_bolsig_rates)")
    if(ZDPlasKin_cfg(3) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_bolsig_rates)")
    ZDPlasKin_cfg(4) = 0.0d0
  else
    if(ZDPlasKin_cfg(4) <= 0.0d0) then
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELECTRON_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
    elseif(lprint .and. ZDPlasKin_cfg(4) < ZDPlasKin_cfg(1)) then
      write(*,"(A)") "ZDPlasKin INFO: ELECTRON_TEMPERATURE is below GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)"
    endif
    ZDPlasKin_cfg(2:3) = 0.0d0
  endif
  density_loc = 0.5d0 * ( sum(density(bolsig_species_index(:))) + sum(abs(density(bolsig_species_index(:)))) )
  if(density_loc <= mach_tiny) then
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined densities configured for BOLSIG+ solver " // &
                                                                     "(subroutine ZDPlasKin_bolsig_rates)")
  elseif( lprint ) then
    call ZDPlasKin_get_density_total(ALL_NEUTRAL=error)
    error = abs( 1.0d0 - density_loc / error )
    if(error > low_density_limit) then
      write(*,"(A,1pd9.2)") "ZDPlasKin INFO: the density of species not configured for BOLSIG+ solver exceeds", error
      low_density_limit = sqrt( low_density_limit )
    endif
  endif
  cfg_loc(1:4) = ZDPlasKin_cfg(1:4)
  cfg_loc(2)   = cfg_loc(2) * 1.0d-6
  cfg_loc(6)   = 0.5d0 * ( density(species_electrons) + abs(density(species_electrons)) )
  cfg_loc(5)   = cfg_loc(6) * 1.0d6
  cfg_loc(6)   = cfg_loc(6) / density_loc
  if(cfg_loc(6) < bolsig_eecol_frac) then
    cfg_loc(6) = 0.0d0
  elseif(lprint .and. leecol) then
    write(*,"(A)") "ZDPlasKin INFO: set electron-electron collisions ON ..."
    leecol = .false.
  endif
  cfg_loc(7:) = 0.5d0 * ( density(bolsig_species_index(:)) + abs(density(bolsig_species_index(:))) ) / density_loc
  if(lbolsig_Maxwell_EEDF .or. bolsig_points_max==0) then
    lforce = .true.
    i = 0
  else
    error = min( max( cfg_loc(3) , bolsig_field_min ) , bolsig_field_max )
    i = int( bolsig_mesh_a * log( error ) - bolsig_mesh_b )
    i = max(0,min(i,bolsig_points_max))
  endif
  if( lforce ) then
    error = 2.0d0
  else
    if(lbolsig_ignore_gas_temp .and. bolsig_cfg(1,i)>0.0d0) then
      cfg_loc(1) = bolsig_cfg(1,i)
      error = 0.0d0
    else
      error = abs( ( cfg_loc(1) - bolsig_cfg(1,i) ) / ( 0.5d0 * ( cfg_loc(1) + bolsig_cfg(1,i) ) + mach_tiny) ) / bolsig_rtol_half
    endif
    if(error <= 1.0d0) then
      error = abs( ( cfg_loc(2) - bolsig_cfg(2,i) ) / ( 0.5d0 * ( cfg_loc(2) + bolsig_cfg(2,i) ) + mach_tiny) ) / bolsig_rtol_half
      if(error <= 1.0d0) then
        error = abs( ( cfg_loc(3) - bolsig_cfg(3,i) ) / ( 0.5d0 * ( cfg_loc(3) + bolsig_cfg(3,i) ) + mach_tiny) ) / bolsig_rtol
        if(error <= 1.0d0) then
          error = abs( ( max(cfg_loc(6),bolsig_eecol_frac) - max(bolsig_cfg(6,i),bolsig_eecol_frac) ) &
           / ( 0.5d0 * ( max(cfg_loc(6),bolsig_eecol_frac) + max(bolsig_cfg(6,i),bolsig_eecol_frac) ) + mach_tiny) ) &
           / bolsig_rtol_half
          if(error <= 1.0d0) error = maxval( abs( cfg_loc(7:) - bolsig_cfg(7:,i) ) ) &
                                     / ( 0.5d0 * maxval( cfg_loc(7:) + bolsig_cfg(7:,i) ) ) / bolsig_rtol
        endif
      endif
    endif
  endif
  if(error > 1.0d0) then
    j = 6 + bolsig_species_max
    k = 10 + bolsig_collisions_max
    bolsig_cfg(:,i) = cfg_loc(:)
    call ZDPlasKin_bolsig_SolveBoltzmann(j,bolsig_cfg(1:j,i),k,bolsig_reslt(1:k,i))
    if(.not. lbolsig_Maxwell_EEDF) then
      bolsig_reslt(2,i) = bolsig_reslt(2, i) * eV_to_K / 1.5d0
    else
      bolsig_reslt(2,i) = cfg_loc(4)
    endif
    bolsig_reslt(3, i) = bolsig_reslt(3, i) * 1.0d-2
    bolsig_reslt(4, i) = bolsig_reslt(4, i) * 1.0d-2 / density_loc
    bolsig_reslt(5, i) = bolsig_reslt(5, i) * 1.0d-2
    bolsig_reslt(6, i) = bolsig_reslt(6, i) * 1.0d-2
    bolsig_reslt(7:,i) = bolsig_reslt(7:,i) * 1.0d6
  endif
  ZDPlasKin_cfg(3:12) = bolsig_reslt(:10,i)
  bolsig_rates(:)     = bolsig_reslt(11:,i)
  return
end subroutine ZDPlasKin_bolsig_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get index of species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_species_index(str,i)
  implicit none
  character(*), intent(in ) :: str
  integer,      intent(out) :: i
  character(species_length) :: string
  integer :: j, istr
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  string = trim(adjustl(str))
  istr   = len_trim(string)
  do i = 1, istr
    j = iachar(string(i:i))
    if(j>=97 .and. j<=122) string(i:i) = achar(j-32)
  enddo
  i = 0
  j = 0
  do while(i==0 .and. j<species_max)
    j = j + 1
    if(string(1:istr) == trim(species_name(j))) i = j
  enddo
  if(i <= 0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: cannot identify species <"//trim(str)//"> (subroutine ZDPlasKin_get_species_index)")
  return
end subroutine ZDPlasKin_get_species_index
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get density for species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in) :: string
  logical, optional, intent(in) :: LDENS_CONST
  double precision, optional, intent(in) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present(DENS) ) density(i) = DENS
  if( present(LDENS_CONST) ) then
    density_constant(i) = LDENS_CONST
    ldensity_constant   = any( density_constant(:) )
  endif
  return
end subroutine ZDPlasKin_set_density
subroutine ZDPlasKin_get_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in ) :: string
  logical, optional, intent(out) :: LDENS_CONST
  double precision, optional, intent(out) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present( DENS)       )  DENS       = density(i)
  if( present(LDENS_CONST) ) LDENS_CONST = density_constant(i)
  return
end subroutine ZDPlasKin_get_density
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get total densities
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_density_total(ALL_SPECIES,ALL_NEUTRAL,ALL_ION_POSITIVE,ALL_ION_NEGATIVE,ALL_CHARGE)
  double precision, optional, intent(out) :: ALL_SPECIES, ALL_NEUTRAL, ALL_ION_POSITIVE, ALL_ION_NEGATIVE, ALL_CHARGE
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ALL_SPECIES)      ) ALL_SPECIES      = sum(density(:))
  if( present(ALL_NEUTRAL)      ) ALL_NEUTRAL      = sum(density(:), mask = species_charge(:)==0)
  if( present(ALL_ION_POSITIVE) ) ALL_ION_POSITIVE = sum(density(:), mask = species_charge(:)>0)
  if( present(ALL_ION_NEGATIVE) ) ALL_ION_NEGATIVE = sum(density(:), mask = species_charge(:)<0) - density(species_electrons)
  if( present(ALL_CHARGE)       ) ALL_CHARGE       = sum(density(:) * dble(species_charge(:)))
  return
end subroutine ZDPlasKin_get_density_total
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get species source terms & reaction rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_rates(SOURCE_TERMS,REACTION_RATES,SOURCE_TERMS_MATRIX,MEAN_DENSITY, &
                               MEAN_SOURCE_TERMS,MEAN_REACTION_RATES,MEAN_SOURCE_TERMS_MATRIX)
  double precision, optional, intent(out) :: SOURCE_TERMS(species_max), REACTION_RATES(reactions_max), &
                                             SOURCE_TERMS_MATRIX(species_max,reactions_max), MEAN_DENSITY(species_max), &
                                             MEAN_SOURCE_TERMS(species_max), MEAN_REACTION_RATES(reactions_max), &
                                             MEAN_SOURCE_TERMS_MATRIX(species_max,reactions_max)
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if(present(SOURCE_TERMS) .or. present(REACTION_RATES) .or. present(SOURCE_TERMS_MATRIX)) then
    dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
    if( present(SOURCE_TERMS)                 ) SOURCE_TERMS(:)   = dens_loc(1:species_max,1)
    if( present(REACTION_RATES)               ) REACTION_RATES(:) = rrt(:)
    if( present(SOURCE_TERMS_MATRIX)          ) call ZDPlasKin_reac_source_matrix(rrt(:),SOURCE_TERMS_MATRIX(:,:))
  endif
  if(present(MEAN_DENSITY)        .or. present(MEAN_SOURCE_TERMS) .or. &
     present(MEAN_REACTION_RATES) .or. present(MEAN_SOURCE_TERMS_MATRIX)) then
    if( lstat_accum ) then
      if(stat_time > 0.0d0) then
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = stat_dens(:) / stat_time
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = stat_src(:)  / stat_time
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = stat_rrt(:)  / stat_time
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) then
          call ZDPlasKin_reac_source_matrix(stat_rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
          MEAN_SOURCE_TERMS_MATRIX(:,:)  =  MEAN_SOURCE_TERMS_MATRIX(:,:) / stat_time
        endif
      else
        dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
        if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
        call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = density(:)
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = dens_loc(1:species_max,1)
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = rrt(:)
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) call ZDPlasKin_reac_source_matrix(rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
      endif
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_get_rates)")
    endif
  endif
  return
end subroutine ZDPlasKin_get_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set config
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_config(ATOL,RTOL,SILENCE_MODE,STAT_ACCUM,QTPLASKIN_SAVE,BOLSIG_EE_FRAC,BOLSIG_IGNORE_GAS_TEMPERATURE)
  use dvode_f90_m, only : set_intermediate_opts
  implicit none
  logical, optional, intent(in) :: SILENCE_MODE, STAT_ACCUM, QTPLASKIN_SAVE, BOLSIG_IGNORE_GAS_TEMPERATURE
  double precision, optional, intent(in) :: ATOL, RTOL, BOLSIG_EE_FRAC
  integer :: i
  logical, save :: lfirst = .true.
  integer, save :: bounded_components(vode_neq)
  double precision :: atol_loc, rtol_loc
  double precision, save :: atol_save = -1.0d0, rtol_save = -1.0d0
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    do i = 1, vode_neq
      bounded_components(i) = i
    enddo
    lfirst = .false.
  endif
  if( present(SILENCE_MODE) ) lprint = ( .not. SILENCE_MODE )
  if( present(BOLSIG_EE_FRAC) ) bolsig_eecol_frac = 0.5d0 * ( BOLSIG_EE_FRAC + abs(BOLSIG_EE_FRAC) )
  if( present(BOLSIG_IGNORE_GAS_TEMPERATURE) ) lbolsig_ignore_gas_temp = BOLSIG_IGNORE_GAS_TEMPERATURE
  if( present(STAT_ACCUM) ) then
    if( lprint ) then
      if(lstat_accum .neqv. STAT_ACCUM) then
        if( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition OFF ..."
        endif
      elseif( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: reset statistic acquisition data ..."
      endif
    endif
    stat_dens(:) = 0.0d0
    stat_src(:)  = 0.0d0
    stat_rrt(:)  = 0.0d0
    stat_time    = 0.0d0
    lstat_accum  = STAT_ACCUM
  endif
  if( present(QTPLASKIN_SAVE) ) then
    if( lprint ) then
      if(lqtplaskin .neqv. QTPLASKIN_SAVE) then
        if( QTPLASKIN_SAVE ) then
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format OFF ..."
        endif
      endif
    endif
    lqtplaskin = QTPLASKIN_SAVE
  endif
  if( present(ATOL) ) then
    atol_loc = ATOL
  else
    atol_loc = atol_save
  endif
  if( present(RTOL) ) then
    rtol_loc = RTOL
  else
    rtol_loc = rtol_save
  endif
  if(min(atol_loc,rtol_loc)<0.0d0 .or. max(atol_loc,rtol_loc)<=0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ATOL/RTOL (ZDPlasKin_set_config)")
  if(atol_loc/=atol_save .or. rtol_loc/=rtol_save) then
    atol_save = atol_loc
    rtol_save = rtol_loc
    if( lprint ) write(*,"(2(A,1pd9.2),A)") "ZDPlasKin INFO: set accuracy", atol_save, " (absolute) &", rtol_save, " (relative)"
    dens_loc(:,0) = 0.0d0
    dens_loc(:,1) = huge(dens_loc)
    vode_options  = set_intermediate_opts(abserr=atol_save,relerr=rtol_save, &
                                          dense_j=.true.,user_supplied_jacobian=.true., &
                                          constrained=bounded_components(:),clower=dens_loc(:,0),cupper=dens_loc(:,1))
    if(vode_istate /= 1) vode_istate = 3
  endif
  return
end subroutine ZDPlasKin_set_config
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get conditions
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,GAS_HEATING,SPEC_HEAT_RATIO,HEAT_SOURCE,SOFT_RESET)
  implicit none
  double precision, optional, intent(in) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, ELEC_TEMPERATURE, &
                                            SPEC_HEAT_RATIO, HEAT_SOURCE
  logical,          optional, intent(in) :: GAS_HEATING, SOFT_RESET
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(GAS_TEMPERATURE) ) then
    if(GAS_TEMPERATURE <= 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(1) = GAS_TEMPERATURE
  endif
  if( present(REDUCED_FREQUENCY) ) then
    if(REDUCED_FREQUENCY < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(2) = REDUCED_FREQUENCY
  endif
  if( present(REDUCED_FIELD) ) then
    if(REDUCED_FIELD < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(3) = REDUCED_FIELD
  endif
  if( present(SOFT_RESET) ) then
    if( SOFT_RESET ) vode_istate = 1
  endif
  if( present(ELEC_TEMPERATURE) ) then
    if(ELEC_TEMPERATURE < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELEC_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    if(ELEC_TEMPERATURE > 0.0d0) then
      lbolsig_Maxwell_EEDF = .true.
    else
      lbolsig_Maxwell_EEDF = .false.
    endif
    ZDPlasKin_cfg(4) = ELEC_TEMPERATURE
  endif
  if( present(GAS_HEATING) ) then
    if(lgas_heating .neqv. GAS_HEATING) then
      if( GAS_HEATING ) then
        if(present(SPEC_HEAT_RATIO) .or. ZDPlasKin_cfg(13)>0.0d0) then
          if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating ON ..."
        else
          ZDPlasKin_cfg(13) = 2.0d0/3.0d0
          if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set gas heating ON; specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
        endif
      else
        if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating OFF ..."
      endif
      lgas_heating = GAS_HEATING
    endif
  endif
  if( present(SPEC_HEAT_RATIO) ) then
    if(SPEC_HEAT_RATIO > 1.0d0) then
      ZDPlasKin_cfg(13) = SPEC_HEAT_RATIO - 1.0d0
      if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong value of SPEC_HEAT_RATIO (subroutine ZDPlasKin_set_conditions)")
    endif
  endif
  if( present(HEAT_SOURCE) ) then
    ZDPlasKin_cfg(14) = HEAT_SOURCE
    if( lprint ) write(*,"(A,1pd9.2,A)") "ZDPlasKin INFO: set heat source =", ZDPlasKin_cfg(14), " W/cm3"
  endif
end subroutine ZDPlasKin_set_conditions
subroutine ZDPlasKin_get_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,ELEC_DRIFT_VELOCITY,ELEC_DIFF_COEFF,ELEC_MOBILITY_N, &
                                    ELEC_MU_EPS_N,ELEC_DIFF_EPS_N,ELEC_FREQUENCY_N, &
                                    ELEC_POWER_N,ELEC_POWER_ELASTIC_N,ELEC_POWER_INELASTIC_N,ELEC_EEDF)
  implicit none
  double precision, optional, intent(out) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, &
                                             ELEC_TEMPERATURE, ELEC_DRIFT_VELOCITY, ELEC_DIFF_COEFF, ELEC_MOBILITY_N, &
                                             ELEC_MU_EPS_N, ELEC_DIFF_EPS_N, ELEC_FREQUENCY_N, &
                                             ELEC_POWER_N, ELEC_POWER_ELASTIC_N, ELEC_POWER_INELASTIC_N
  double precision, optional, dimension(:,:), intent(out) :: ELEC_EEDF
  integer :: i
  double precision :: x,y
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ELEC_EEDF) ) then
    call ZDPlasKin_bolsig_rates(lbolsig_force=.true.)
  else
    call ZDPlasKin_bolsig_rates()
  endif
  if( present(GAS_TEMPERATURE)        ) GAS_TEMPERATURE        = ZDPlasKin_cfg(1)
  if( present(REDUCED_FREQUENCY)      ) REDUCED_FREQUENCY      = ZDPlasKin_cfg(2)
  if( present(REDUCED_FIELD)          ) REDUCED_FIELD          = ZDPlasKin_cfg(3)
  if( present(ELEC_TEMPERATURE)       ) ELEC_TEMPERATURE       = ZDPlasKin_cfg(4)
  if( present(ELEC_DRIFT_VELOCITY)    ) ELEC_DRIFT_VELOCITY    = ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
  if( present(ELEC_DIFF_COEFF)        ) ELEC_DIFF_COEFF        = ZDPlasKin_cfg(6)
  if( present(ELEC_MOBILITY_N)        ) ELEC_MOBILITY_N        = ZDPlasKin_cfg(5)
  if( present(ELEC_MU_EPS_N)          ) ELEC_MU_EPS_N          = ZDPlasKin_cfg(7)
  if( present(ELEC_DIFF_EPS_N)        ) ELEC_DIFF_EPS_N        = ZDPlasKin_cfg(8)
  if( present(ELEC_FREQUENCY_N)       ) ELEC_FREQUENCY_N       = ZDPlasKin_cfg(9)
  if( present(ELEC_POWER_N)           ) ELEC_POWER_N           = ZDPlasKin_cfg(10)
  if( present(ELEC_POWER_ELASTIC_N)   ) ELEC_POWER_ELASTIC_N   = ZDPlasKin_cfg(11)
  if( present(ELEC_POWER_INELASTIC_N) ) ELEC_POWER_INELASTIC_N = ZDPlasKin_cfg(12)
  if( present(ELEC_EEDF) ) then
    ELEC_EEDF = 0d0
  	 if( size(ELEC_EEDF,dim=1) < 2 ) then
      if(lprint) write(*,"(A)") &
  	     "ZDPlasKin WARNING: insufficient first dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
  	  else
  		y = 1.0d0
  		do i = 1, size(ELEC_EEDF,dim=2)
  		  call ZDPlasKin_bolsig_GetEEDF(i,x,y)
  		  if( x >= 0d0 .and. y > 0d0) then
  			ELEC_EEDF(1,i) = x
  			ELEC_EEDF(2,i) = y
  		  else
  			exit
  		  endif
  		enddo
  		if(lprint .and. y>0d0) write(*,"(A)") &
  		  "ZDPlasKin WARNING: insufficient second dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
     endif
  endif
end subroutine ZDPlasKin_get_conditions
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reset
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reset()
  implicit none
  vode_istate         =  1
  density(:)          =  0.0d0
  ZDPlasKin_cfg(:)    =  0.0d0
  ldensity_constant   = .false.
  density_constant(:) = .false.
  lreaction_block(:)  = .false.
  lprint              = .true.
  lstat_accum         = .false.
  lqtplaskin          = .false.
  lgas_heating        = .false.
  bolsig_eecol_frac       = bolsig_eecol_frac_def
  lbolsig_ignore_gas_temp = .false.
  lbolsig_Maxwell_EEDF    = .false.
  write(*,"(A)") "ZDPlasKin INFO: reset data and configuration"
  call ZDPlasKin_set_config(ATOL=vode_atol,RTOL=vode_rtol)
  return
end subroutine ZDPlasKin_reset
!-----------------------------------------------------------------------------------------------------------------------------------
!
! stop
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_stop(string)
  implicit none
  character(*), intent(in) :: string
  if(string /= "") write(*,"(A)") trim(string)
  write(*,"(A,$)") "PRESS ENTER TO EXIT ... "
  read(*,*)
  stop
end subroutine ZDPlasKin_stop
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data to file
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_file(FILE_SPECIES,FILE_REACTIONS,FILE_SOURCE_MATRIX,FILE_UNIT)
  implicit none
  character(*), optional, intent(in) :: FILE_SPECIES, FILE_REACTIONS, FILE_SOURCE_MATRIX
  integer, optional, intent(in) :: FILE_UNIT
  logical :: lerror
  integer :: i
  if( present(FILE_UNIT) ) ifile_unit = FILE_UNIT
  if( present(FILE_SPECIES) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SPECIES)),action="write",err=100)
    do i = 1, species_max
      write(ifile_unit,111,err=100) i, species_name(i)
    enddo
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SPECIES)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
111 format(i2,1x,A9)
  endif
  if( present(FILE_REACTIONS) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_REACTIONS)),action="write",err=200)
    do i = 1, reactions_max
      write(ifile_unit,211,err=200) i, reaction_sign(i)
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_REACTIONS)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
211 format(i3,1x,A40)
  endif
  if( present(FILE_SOURCE_MATRIX) ) then
    if( lstat_accum ) then
      call ZDPlasKin_reac_source_matrix(stat_rrt(:),mrtm(:,:))
      if(stat_time > 0.0d0) mrtm(:,:) = mrtm(:,:) / stat_time
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_write_file)")
    endif
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SOURCE_MATRIX)),action="write",err=300)
    write(ifile_unit,311,err=300) ( i, i = 1, species_max )
    write(ifile_unit,312,err=300) "N", "reaction", ( trim(species_name(i)), i = 1, species_max )
    do i = 1, reactions_max
      write(ifile_unit,313,err=300) i, reaction_sign(i), mrtm(:,i)
    enddo
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SOURCE_MATRIX)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
311 format(451x,44(1x,i9))
312 format(A3,1x,A40,1x,44(1x,A9))
313 format(i3,1x,A40,1x,44(1x,1pd9.2))
  endif
  return
end subroutine ZDPlasKin_write_file
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data in qtplaskin format
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_qtplaskin(time,LFORCE_WRITE)
  implicit none
  double precision, intent(in) :: time
  logical, optional, intent(in) :: LFORCE_WRITE
  integer, parameter :: idef_data = 5
  character(24), parameter :: qtplaskin_names(idef_data) = (/ "Reduced field [Td]      ", "Gas temperature [K]     ", &
                                  "Electron temperature [K]", "Current density [A/cm2] ", "Power density [W/cm3]   " /)
  double precision, save :: densav(0:species_max,2) = -huge(densav)
  double precision :: rtol, cond(idef_data)
  logical, save :: lfirst = .true.
  logical :: lerror
  integer, save :: iuser_data = 0
  integer :: i
  if( time < densav(0,1) ) lfirst = .true.
  if( lfirst ) then
    call ZDPlasKin_write_file(FILE_SPECIES="qt_species_list.txt",FILE_REACTIONS="qt_reactions_list.txt")
    if( allocated(qtplaskin_user_data) ) then
      iuser_data = size(qtplaskin_user_data)
      iuser_data = min(iuser_data,90)
      if( iuser_data > 0 ) then
        if( allocated(qtplaskin_user_names) ) then
          if( size(qtplaskin_user_names) /= iuser_data ) deallocate(qtplaskin_user_names)
        endif
        if( .not. allocated(qtplaskin_user_names) ) then
          allocate(qtplaskin_user_names(iuser_data))
          do i = 1, iuser_data
            write(qtplaskin_user_names(i),"(A,i2.2)") "user defined #", i
          enddo
        endif
      endif
    endif
    lerror = .true.
    open(ifile_unit,file="qt_conditions_list.txt",action="write",err=100)
    do i = 1, idef_data
      write(ifile_unit,"(i3,1x,A)",err=100) i, trim(adjustl(qtplaskin_names(i)))
    enddo
    if( iuser_data > 0 ) then
      do i = 1, iuser_data
        write(ifile_unit,"(i3,1x,A)",err=100) (i+idef_data), trim(adjustl(qtplaskin_user_names(i)))
      enddo
    endif
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions_list.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    rrt(:) = 1.0d0
    call ZDPlasKin_reac_source_matrix(rrt(:),mrtm(:,:))
    open(ifile_unit,file="qt_matrix.txt",action="write",err=200)
    do i = 1, species_max
      write(ifile_unit,"(425(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,44(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_densities.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_conditions.txt",action="write",err=400)
    write(ifile_unit,"(1x,A12,$)",err=400) "Time_s"
    do i = 1, idef_data + iuser_data
      write(ifile_unit,"(11x,i2.2,$)",err=400) i
    enddo
    write(ifile_unit,*,err=400)
    lerror = .false.
400 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_rates.txt",action="write",err=500)
    write(ifile_unit,"(1x,A12,425(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
    lerror = .false.
500 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_rates.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
  endif
  if( present(LFORCE_WRITE) ) then
    if( LFORCE_WRITE ) lfirst = .true.
  endif
  rtol = 10.0d0 ** ( floor( log10( abs(densav(0,1)) + tiny(rtol) ) ) - 6 )
  if( ( time - densav(0,1) ) >= rtol .or. lfirst ) then
    densav(0,2) = time
    densav(1:species_max,2) = density(:)
    where( densav(:,2) < 1.0d-99 ) densav(:,2) = 0.0d0
    if( time > 2.0d0 * densav(0,1) ) then
      rtol = huge(rtol)
    else
      rtol = maxval( abs(densav(1:,1)-densav(1:,2)) / ( abs(densav(1:,1)+densav(1:,2))/2.0d0 + qtplaskin_atol ) )
    endif
    if( rtol > qtplaskin_rtol .or. lfirst ) then
      open(ifile_unit,file="qt_densities.txt",access="append")
      write(ifile_unit,"(1pe15.6,44(1pe13.4))") densav(0,2), densav(1:,2)
      close(ifile_unit)
      open(ifile_unit,file="qt_conditions.txt",access="append")
      cond(1) = ZDPlasKin_cfg(3)
      cond(2) = ZDPlasKin_cfg(1)
      cond(3) = ZDPlasKin_cfg(4)
      cond(4) = q_elem * density(species_electrons) * ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
      call ZDPlasKin_get_density_total(ALL_NEUTRAL=cond(5))
      cond(5) = cond(4) * cond(5) * ZDPlasKin_cfg(3) * 1.0d-17
      where( abs(cond(:)) < 1.0d-99 ) cond(:) = 0.0d0
      write(ifile_unit,"(6(1pe13.4),$)") densav(0,2), cond(:)
      if( iuser_data > 0 ) then
        where( abs(qtplaskin_user_data(1:iuser_data)) < 1.0d-99 ) qtplaskin_user_data(1:iuser_data) = 0.0d0
        write(ifile_unit,"(90(1pe13.4))") qtplaskin_user_data(1:iuser_data)
      else
        write(ifile_unit,*)
      endif
      close(ifile_unit)
      call ZDPlasKin_get_rates(REACTION_RATES=rrt_loc)
      where( abs(rrt_loc(:)) < 1.0d-99 ) rrt_loc(:) = 0.0d0
      open(ifile_unit,file="qt_rates.txt",access="append")
      write(ifile_unit,"(426(1pe13.4))") densav(0,2), rrt_loc(:)
      close(ifile_unit)
      densav(:,1) = densav(:,2)
    endif
  endif
  lfirst = .false.
  lqtplaskin_first = .false.
  return
end subroutine ZDPlasKin_write_qtplaskin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction sensitivity acquisition
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_source_matrix(reac_rate_local,reac_source_local)
  implicit none
  double precision, intent(in)  :: reac_rate_local(reactions_max)
  double precision, intent(out) :: reac_source_local(species_max,reactions_max)
  reac_source_local(:,:) = 0.0d0
  reac_source_local(01,001) = - reac_rate_local(001) 
  reac_source_local(02,001) = + reac_rate_local(001) 
  reac_source_local(01,002) = - reac_rate_local(002) 
  reac_source_local(02,002) = + reac_rate_local(002) 
  reac_source_local(01,003) = - reac_rate_local(003) 
  reac_source_local(03,003) = + reac_rate_local(003) 
  reac_source_local(01,004) = - reac_rate_local(004) 
  reac_source_local(04,004) = + reac_rate_local(004) 
  reac_source_local(01,005) = - reac_rate_local(005) 
  reac_source_local(05,005) = + reac_rate_local(005) 
  reac_source_local(01,006) = - reac_rate_local(006) 
  reac_source_local(06,006) = + reac_rate_local(006) 
  reac_source_local(01,007) = - reac_rate_local(007) 
  reac_source_local(07,007) = + reac_rate_local(007) 
  reac_source_local(01,008) = - reac_rate_local(008) 
  reac_source_local(08,008) = + reac_rate_local(008) 
  reac_source_local(01,009) = - reac_rate_local(009) 
  reac_source_local(09,009) = + reac_rate_local(009) 
  reac_source_local(01,010) = + reac_rate_local(010) 
  reac_source_local(02,010) = - reac_rate_local(010) 
  reac_source_local(01,011) = + reac_rate_local(011) 
  reac_source_local(03,011) = - reac_rate_local(011) 
  reac_source_local(01,012) = + reac_rate_local(012) 
  reac_source_local(04,012) = - reac_rate_local(012) 
  reac_source_local(01,013) = + reac_rate_local(013) 
  reac_source_local(05,013) = - reac_rate_local(013) 
  reac_source_local(01,014) = + reac_rate_local(014) 
  reac_source_local(06,014) = - reac_rate_local(014) 
  reac_source_local(01,015) = + reac_rate_local(015) 
  reac_source_local(07,015) = - reac_rate_local(015) 
  reac_source_local(01,016) = + reac_rate_local(016) 
  reac_source_local(08,016) = - reac_rate_local(016) 
  reac_source_local(01,017) = + reac_rate_local(017) 
  reac_source_local(09,017) = - reac_rate_local(017) 
  reac_source_local(21,018) = - reac_rate_local(018) 
  reac_source_local(22,018) = + reac_rate_local(018) 
  reac_source_local(21,019) = - reac_rate_local(019) 
  reac_source_local(22,019) = + reac_rate_local(019) 
  reac_source_local(21,020) = - reac_rate_local(020) 
  reac_source_local(23,020) = + reac_rate_local(020) 
  reac_source_local(21,021) = - reac_rate_local(021) 
  reac_source_local(23,021) = + reac_rate_local(021) 
  reac_source_local(21,022) = - reac_rate_local(022) 
  reac_source_local(24,022) = + reac_rate_local(022) 
  reac_source_local(21,023) = - reac_rate_local(023) 
  reac_source_local(25,023) = + reac_rate_local(023) 
  reac_source_local(21,024) = + reac_rate_local(024) 
  reac_source_local(22,024) = - reac_rate_local(024) 
  reac_source_local(21,025) = + reac_rate_local(025) 
  reac_source_local(23,025) = - reac_rate_local(025) 
  reac_source_local(21,026) = + reac_rate_local(026) 
  reac_source_local(24,026) = - reac_rate_local(026) 
  reac_source_local(21,027) = + reac_rate_local(027) 
  reac_source_local(25,027) = - reac_rate_local(027) 
  reac_source_local(01,028) = + reac_rate_local(028) 
  reac_source_local(02,028) = - reac_rate_local(028) 
  reac_source_local(02,029) = + reac_rate_local(029) 
  reac_source_local(03,029) = - reac_rate_local(029) 
  reac_source_local(03,030) = + reac_rate_local(030) 
  reac_source_local(04,030) = - reac_rate_local(030) 
  reac_source_local(04,031) = + reac_rate_local(031) 
  reac_source_local(05,031) = - reac_rate_local(031) 
  reac_source_local(05,032) = + reac_rate_local(032) 
  reac_source_local(06,032) = - reac_rate_local(032) 
  reac_source_local(06,033) = + reac_rate_local(033) 
  reac_source_local(07,033) = - reac_rate_local(033) 
  reac_source_local(07,034) = + reac_rate_local(034) 
  reac_source_local(08,034) = - reac_rate_local(034) 
  reac_source_local(08,035) = + reac_rate_local(035) 
  reac_source_local(09,035) = - reac_rate_local(035) 
  reac_source_local(01,036) = - reac_rate_local(036) 
  reac_source_local(02,036) = + reac_rate_local(036) 
  reac_source_local(02,037) = - reac_rate_local(037) 
  reac_source_local(03,037) = + reac_rate_local(037) 
  reac_source_local(03,038) = - reac_rate_local(038) 
  reac_source_local(04,038) = + reac_rate_local(038) 
  reac_source_local(04,039) = - reac_rate_local(039) 
  reac_source_local(05,039) = + reac_rate_local(039) 
  reac_source_local(05,040) = - reac_rate_local(040) 
  reac_source_local(06,040) = + reac_rate_local(040) 
  reac_source_local(06,041) = - reac_rate_local(041) 
  reac_source_local(07,041) = + reac_rate_local(041) 
  reac_source_local(07,042) = - reac_rate_local(042) 
  reac_source_local(08,042) = + reac_rate_local(042) 
  reac_source_local(08,043) = - reac_rate_local(043) 
  reac_source_local(09,043) = + reac_rate_local(043) 
  reac_source_local(01,044) = + reac_rate_local(044) 
  reac_source_local(02,044) = - reac_rate_local(044) 
  reac_source_local(02,045) = + reac_rate_local(045) 
  reac_source_local(03,045) = - reac_rate_local(045) 
  reac_source_local(03,046) = + reac_rate_local(046) 
  reac_source_local(04,046) = - reac_rate_local(046) 
  reac_source_local(04,047) = + reac_rate_local(047) 
  reac_source_local(05,047) = - reac_rate_local(047) 
  reac_source_local(05,048) = + reac_rate_local(048) 
  reac_source_local(06,048) = - reac_rate_local(048) 
  reac_source_local(06,049) = + reac_rate_local(049) 
  reac_source_local(07,049) = - reac_rate_local(049) 
  reac_source_local(07,050) = + reac_rate_local(050) 
  reac_source_local(08,050) = - reac_rate_local(050) 
  reac_source_local(08,051) = + reac_rate_local(051) 
  reac_source_local(09,051) = - reac_rate_local(051) 
  reac_source_local(01,052) = - reac_rate_local(052) 
  reac_source_local(02,052) = + reac_rate_local(052) 
  reac_source_local(02,053) = - reac_rate_local(053) 
  reac_source_local(03,053) = + reac_rate_local(053) 
  reac_source_local(03,054) = - reac_rate_local(054) 
  reac_source_local(04,054) = + reac_rate_local(054) 
  reac_source_local(04,055) = - reac_rate_local(055) 
  reac_source_local(05,055) = + reac_rate_local(055) 
  reac_source_local(05,056) = - reac_rate_local(056) 
  reac_source_local(06,056) = + reac_rate_local(056) 
  reac_source_local(06,057) = - reac_rate_local(057) 
  reac_source_local(07,057) = + reac_rate_local(057) 
  reac_source_local(07,058) = - reac_rate_local(058) 
  reac_source_local(08,058) = + reac_rate_local(058) 
  reac_source_local(08,059) = - reac_rate_local(059) 
  reac_source_local(09,059) = + reac_rate_local(059) 
  reac_source_local(01,060) = + reac_rate_local(060) 
  reac_source_local(02,060) = - reac_rate_local(060) 
  reac_source_local(02,061) = + reac_rate_local(061) 
  reac_source_local(03,061) = - reac_rate_local(061) 
  reac_source_local(03,062) = + reac_rate_local(062) 
  reac_source_local(04,062) = - reac_rate_local(062) 
  reac_source_local(04,063) = + reac_rate_local(063) 
  reac_source_local(05,063) = - reac_rate_local(063) 
  reac_source_local(05,064) = + reac_rate_local(064) 
  reac_source_local(06,064) = - reac_rate_local(064) 
  reac_source_local(06,065) = + reac_rate_local(065) 
  reac_source_local(07,065) = - reac_rate_local(065) 
  reac_source_local(07,066) = + reac_rate_local(066) 
  reac_source_local(08,066) = - reac_rate_local(066) 
  reac_source_local(08,067) = + reac_rate_local(067) 
  reac_source_local(09,067) = - reac_rate_local(067) 
  reac_source_local(01,068) = - reac_rate_local(068) 
  reac_source_local(02,068) = + reac_rate_local(068) 
  reac_source_local(02,069) = - reac_rate_local(069) 
  reac_source_local(03,069) = + reac_rate_local(069) 
  reac_source_local(03,070) = - reac_rate_local(070) 
  reac_source_local(04,070) = + reac_rate_local(070) 
  reac_source_local(04,071) = - reac_rate_local(071) 
  reac_source_local(05,071) = + reac_rate_local(071) 
  reac_source_local(05,072) = - reac_rate_local(072) 
  reac_source_local(06,072) = + reac_rate_local(072) 
  reac_source_local(06,073) = - reac_rate_local(073) 
  reac_source_local(07,073) = + reac_rate_local(073) 
  reac_source_local(07,074) = - reac_rate_local(074) 
  reac_source_local(08,074) = + reac_rate_local(074) 
  reac_source_local(08,075) = - reac_rate_local(075) 
  reac_source_local(09,075) = + reac_rate_local(075) 
  reac_source_local(21,076) = + reac_rate_local(076) 
  reac_source_local(22,076) = - reac_rate_local(076) 
  reac_source_local(22,077) = + reac_rate_local(077) 
  reac_source_local(23,077) = - reac_rate_local(077) 
  reac_source_local(23,078) = + reac_rate_local(078) 
  reac_source_local(24,078) = - reac_rate_local(078) 
  reac_source_local(24,079) = + reac_rate_local(079) 
  reac_source_local(25,079) = - reac_rate_local(079) 
  reac_source_local(21,080) = - reac_rate_local(080) 
  reac_source_local(22,080) = + reac_rate_local(080) 
  reac_source_local(22,081) = - reac_rate_local(081) 
  reac_source_local(23,081) = + reac_rate_local(081) 
  reac_source_local(23,082) = - reac_rate_local(082) 
  reac_source_local(24,082) = + reac_rate_local(082) 
  reac_source_local(24,083) = - reac_rate_local(083) 
  reac_source_local(25,083) = + reac_rate_local(083) 
  reac_source_local(21,084) = + reac_rate_local(084) 
  reac_source_local(22,084) = - reac_rate_local(084) 
  reac_source_local(22,085) = + reac_rate_local(085) 
  reac_source_local(23,085) = - reac_rate_local(085) 
  reac_source_local(23,086) = + reac_rate_local(086) 
  reac_source_local(24,086) = - reac_rate_local(086) 
  reac_source_local(24,087) = + reac_rate_local(087) 
  reac_source_local(25,087) = - reac_rate_local(087) 
  reac_source_local(21,088) = - reac_rate_local(088) 
  reac_source_local(22,088) = + reac_rate_local(088) 
  reac_source_local(22,089) = - reac_rate_local(089) 
  reac_source_local(23,089) = + reac_rate_local(089) 
  reac_source_local(23,090) = - reac_rate_local(090) 
  reac_source_local(24,090) = + reac_rate_local(090) 
  reac_source_local(24,091) = - reac_rate_local(091) 
  reac_source_local(25,091) = + reac_rate_local(091) 
  reac_source_local(01,092) = - reac_rate_local(092) 
  reac_source_local(10,092) = + reac_rate_local(092) 
  reac_source_local(01,093) = - reac_rate_local(093) 
  reac_source_local(10,093) = + reac_rate_local(093) 
  reac_source_local(01,094) = - reac_rate_local(094) 
  reac_source_local(11,094) = + reac_rate_local(094) 
  reac_source_local(01,095) = - reac_rate_local(095) 
  reac_source_local(11,095) = + reac_rate_local(095) 
  reac_source_local(01,096) = - reac_rate_local(096) 
  reac_source_local(11,096) = + reac_rate_local(096) 
  reac_source_local(01,097) = - reac_rate_local(097) 
  reac_source_local(12,097) = + reac_rate_local(097) 
  reac_source_local(01,098) = - reac_rate_local(098) 
  reac_source_local(12,098) = + reac_rate_local(098) 
  reac_source_local(01,099) = - reac_rate_local(099) 
  reac_source_local(12,099) = + reac_rate_local(099) 
  reac_source_local(01,100) = - reac_rate_local(100) 
  reac_source_local(13,100) = + reac_rate_local(100) 
  reac_source_local(01,101) = - reac_rate_local(101) 
  reac_source_local(13,101) = + reac_rate_local(101) 
  reac_source_local(01,102) = - reac_rate_local(102) 
  reac_source_local(13,102) = + reac_rate_local(102) 
  reac_source_local(01,103) = - reac_rate_local(103) 
  reac_source_local(14,103) = + reac_rate_local(103) 
  reac_source_local(15,103) = + reac_rate_local(103) 
  reac_source_local(01,104) = - reac_rate_local(104) 
  reac_source_local(18,104) = + reac_rate_local(104) 
  reac_source_local(44,104) = + reac_rate_local(104) 
  reac_source_local(10,105) = - reac_rate_local(105) 
  reac_source_local(18,105) = + reac_rate_local(105) 
  reac_source_local(44,105) = + reac_rate_local(105) 
  reac_source_local(14,106) = - reac_rate_local(106) 
  reac_source_local(17,106) = + reac_rate_local(106) 
  reac_source_local(44,106) = + reac_rate_local(106) 
  reac_source_local(21,107) = - reac_rate_local(107) 
  reac_source_local(26,107) = + reac_rate_local(107) 
  reac_source_local(21,108) = - reac_rate_local(108) 
  reac_source_local(27,108) = + reac_rate_local(108) 
  reac_source_local(21,109) = - reac_rate_local(109) 
  reac_source_local(28,109) = + reac_rate_local(109) 
  reac_source_local(21,110) = - reac_rate_local(110) 
  reac_source_local(29,110) = + reac_rate_local(110) * 2.d0
  reac_source_local(21,111) = - reac_rate_local(111) 
  reac_source_local(29,111) = + reac_rate_local(111) 
  reac_source_local(30,111) = + reac_rate_local(111) 
  reac_source_local(21,112) = - reac_rate_local(112) 
  reac_source_local(29,112) = + reac_rate_local(112) 
  reac_source_local(31,112) = + reac_rate_local(112) 
  reac_source_local(21,113) = - reac_rate_local(113) 
  reac_source_local(34,113) = + reac_rate_local(113) 
  reac_source_local(44,113) = + reac_rate_local(113) 
  reac_source_local(21,114) = - reac_rate_local(114) 
  reac_source_local(29,114) = + reac_rate_local(114) 
  reac_source_local(36,114) = + reac_rate_local(114) 
  reac_source_local(44,114) = - reac_rate_local(114) 
  reac_source_local(21,115) = + reac_rate_local(115) 
  reac_source_local(26,115) = - reac_rate_local(115) 
  reac_source_local(26,116) = - reac_rate_local(116) 
  reac_source_local(29,116) = + reac_rate_local(116) * 2.d0
  reac_source_local(26,117) = - reac_rate_local(117) 
  reac_source_local(34,117) = + reac_rate_local(117) 
  reac_source_local(44,117) = + reac_rate_local(117) 
  reac_source_local(29,118) = - reac_rate_local(118) 
  reac_source_local(30,118) = + reac_rate_local(118) 
  reac_source_local(29,119) = - reac_rate_local(119) 
  reac_source_local(31,119) = + reac_rate_local(119) 
  reac_source_local(29,120) = - reac_rate_local(120) 
  reac_source_local(33,120) = + reac_rate_local(120) 
  reac_source_local(44,120) = + reac_rate_local(120) 
  reac_source_local(40,121) = - reac_rate_local(121) 
  reac_source_local(41,121) = + reac_rate_local(121) 
  reac_source_local(44,121) = + reac_rate_local(121) 
  reac_source_local(14,122) = + reac_rate_local(122) * 2.d0
  reac_source_local(18,122) = - reac_rate_local(122) 
  reac_source_local(44,122) = - reac_rate_local(122) 
  reac_source_local(14,123) = + reac_rate_local(123) 
  reac_source_local(15,123) = + reac_rate_local(123) 
  reac_source_local(18,123) = - reac_rate_local(123) 
  reac_source_local(44,123) = - reac_rate_local(123) 
  reac_source_local(14,124) = + reac_rate_local(124) 
  reac_source_local(16,124) = + reac_rate_local(124) 
  reac_source_local(18,124) = - reac_rate_local(124) 
  reac_source_local(44,124) = - reac_rate_local(124) 
  reac_source_local(01,125) = + reac_rate_local(125) 
  reac_source_local(14,125) = + reac_rate_local(125) 
  reac_source_local(19,125) = - reac_rate_local(125) 
  reac_source_local(44,125) = - reac_rate_local(125) 
  reac_source_local(01,126) = + reac_rate_local(126) * 2.d0
  reac_source_local(20,126) = - reac_rate_local(126) 
  reac_source_local(44,126) = - reac_rate_local(126) 
  reac_source_local(29,127) = + reac_rate_local(127) * 2.d0
  reac_source_local(34,127) = - reac_rate_local(127) 
  reac_source_local(44,127) = - reac_rate_local(127) 
  reac_source_local(29,128) = + reac_rate_local(128) 
  reac_source_local(30,128) = + reac_rate_local(128) 
  reac_source_local(34,128) = - reac_rate_local(128) 
  reac_source_local(44,128) = - reac_rate_local(128) 
  reac_source_local(29,129) = + reac_rate_local(129) 
  reac_source_local(31,129) = + reac_rate_local(129) 
  reac_source_local(34,129) = - reac_rate_local(129) 
  reac_source_local(44,129) = - reac_rate_local(129) 
  reac_source_local(21,130) = + reac_rate_local(130) * 2.d0
  reac_source_local(35,130) = - reac_rate_local(130) 
  reac_source_local(44,130) = - reac_rate_local(130) 
  reac_source_local(14,131) = + reac_rate_local(131) 
  reac_source_local(29,131) = + reac_rate_local(131) 
  reac_source_local(41,131) = - reac_rate_local(131) 
  reac_source_local(44,131) = - reac_rate_local(131) 
  reac_source_local(15,132) = + reac_rate_local(132) 
  reac_source_local(29,132) = + reac_rate_local(132) 
  reac_source_local(41,132) = - reac_rate_local(132) 
  reac_source_local(44,132) = - reac_rate_local(132) 
  reac_source_local(01,133) = + reac_rate_local(133) 
  reac_source_local(21,133) = + reac_rate_local(133) 
  reac_source_local(43,133) = - reac_rate_local(133) 
  reac_source_local(44,133) = - reac_rate_local(133) 
  reac_source_local(14,134) = + reac_rate_local(134) 
  reac_source_local(17,134) = - reac_rate_local(134) 
  reac_source_local(44,134) = - reac_rate_local(134) 
  reac_source_local(29,135) = + reac_rate_local(135) 
  reac_source_local(33,135) = - reac_rate_local(135) 
  reac_source_local(44,135) = - reac_rate_local(135) 
  reac_source_local(14,136) = + reac_rate_local(136) 
  reac_source_local(17,136) = - reac_rate_local(136) 
  reac_source_local(44,136) = - reac_rate_local(136) 
  reac_source_local(29,137) = + reac_rate_local(137) 
  reac_source_local(33,137) = - reac_rate_local(137) 
  reac_source_local(44,137) = - reac_rate_local(137) 
  reac_source_local(29,138) = + reac_rate_local(138) 
  reac_source_local(32,138) = - reac_rate_local(138) 
  reac_source_local(37,138) = + reac_rate_local(138) 
  reac_source_local(44,138) = - reac_rate_local(138) 
  reac_source_local(21,139) = + reac_rate_local(139) 
  reac_source_local(32,139) = - reac_rate_local(139) 
  reac_source_local(36,139) = + reac_rate_local(139) 
  reac_source_local(44,139) = - reac_rate_local(139) 
  reac_source_local(29,140) = - reac_rate_local(140) 
  reac_source_local(36,140) = + reac_rate_local(140) 
  reac_source_local(44,140) = - reac_rate_local(140) 
  reac_source_local(21,141) = - reac_rate_local(141) 
  reac_source_local(37,141) = + reac_rate_local(141) 
  reac_source_local(44,141) = - reac_rate_local(141) 
  reac_source_local(32,142) = - reac_rate_local(142) 
  reac_source_local(38,142) = + reac_rate_local(142) 
  reac_source_local(44,142) = - reac_rate_local(142) 
  reac_source_local(21,143) = - reac_rate_local(143) 
  reac_source_local(37,143) = + reac_rate_local(143) 
  reac_source_local(44,143) = - reac_rate_local(143) 
  reac_source_local(40,144) = - reac_rate_local(144) 
  reac_source_local(42,144) = + reac_rate_local(144) 
  reac_source_local(44,144) = - reac_rate_local(144) 
  reac_source_local(21,145) = + reac_rate_local(145) 
  reac_source_local(29,145) = - reac_rate_local(145) 
  reac_source_local(36,145) = - reac_rate_local(145) 
  reac_source_local(44,145) = + reac_rate_local(145) 
  reac_source_local(14,146) = - reac_rate_local(146) 
  reac_source_local(36,146) = - reac_rate_local(146) 
  reac_source_local(40,146) = + reac_rate_local(146) 
  reac_source_local(44,146) = + reac_rate_local(146) 
  reac_source_local(21,147) = - reac_rate_local(147) 
  reac_source_local(32,147) = + reac_rate_local(147) 
  reac_source_local(36,147) = - reac_rate_local(147) 
  reac_source_local(44,147) = + reac_rate_local(147) 
  reac_source_local(26,148) = - reac_rate_local(148) 
  reac_source_local(32,148) = + reac_rate_local(148) 
  reac_source_local(36,148) = - reac_rate_local(148) 
  reac_source_local(44,148) = + reac_rate_local(148) 
  reac_source_local(21,149) = + reac_rate_local(149) 
  reac_source_local(27,149) = - reac_rate_local(149) 
  reac_source_local(29,149) = + reac_rate_local(149) 
  reac_source_local(36,149) = - reac_rate_local(149) 
  reac_source_local(44,149) = + reac_rate_local(149) 
  reac_source_local(01,150) = + reac_rate_local(150) 
  reac_source_local(10,150) = - reac_rate_local(150) 
  reac_source_local(29,150) = + reac_rate_local(150) 
  reac_source_local(36,150) = - reac_rate_local(150) 
  reac_source_local(44,150) = + reac_rate_local(150) 
  reac_source_local(01,151) = + reac_rate_local(151) 
  reac_source_local(11,151) = - reac_rate_local(151) 
  reac_source_local(29,151) = + reac_rate_local(151) 
  reac_source_local(36,151) = - reac_rate_local(151) 
  reac_source_local(44,151) = + reac_rate_local(151) 
  reac_source_local(21,152) = + reac_rate_local(152) * 2.d0
  reac_source_local(32,152) = - reac_rate_local(152) 
  reac_source_local(36,152) = - reac_rate_local(152) 
  reac_source_local(44,152) = + reac_rate_local(152) 
  reac_source_local(29,153) = - reac_rate_local(153) 
  reac_source_local(32,153) = + reac_rate_local(153) 
  reac_source_local(37,153) = - reac_rate_local(153) 
  reac_source_local(44,153) = + reac_rate_local(153) 
  reac_source_local(21,154) = + reac_rate_local(154) 
  reac_source_local(37,154) = - reac_rate_local(154) 
  reac_source_local(44,154) = + reac_rate_local(154) 
  reac_source_local(21,155) = + reac_rate_local(155) * 2.d0
  reac_source_local(26,155) = - reac_rate_local(155) 
  reac_source_local(37,155) = - reac_rate_local(155) 
  reac_source_local(44,155) = + reac_rate_local(155) 
  reac_source_local(21,156) = + reac_rate_local(156) * 2.d0
  reac_source_local(27,156) = - reac_rate_local(156) 
  reac_source_local(37,156) = - reac_rate_local(156) 
  reac_source_local(44,156) = + reac_rate_local(156) 
  reac_source_local(21,157) = + reac_rate_local(157) 
  reac_source_local(37,157) = - reac_rate_local(157) 
  reac_source_local(44,157) = + reac_rate_local(157) 
  reac_source_local(01,158) = + reac_rate_local(158) 
  reac_source_local(10,158) = - reac_rate_local(158) 
  reac_source_local(21,158) = + reac_rate_local(158) 
  reac_source_local(37,158) = - reac_rate_local(158) 
  reac_source_local(44,158) = + reac_rate_local(158) 
  reac_source_local(01,159) = + reac_rate_local(159) 
  reac_source_local(11,159) = - reac_rate_local(159) 
  reac_source_local(21,159) = + reac_rate_local(159) 
  reac_source_local(37,159) = - reac_rate_local(159) 
  reac_source_local(44,159) = + reac_rate_local(159) 
  reac_source_local(21,160) = + reac_rate_local(160) * 2.d0
  reac_source_local(29,160) = - reac_rate_local(160) 
  reac_source_local(38,160) = - reac_rate_local(160) 
  reac_source_local(44,160) = + reac_rate_local(160) 
  reac_source_local(01,161) = + reac_rate_local(161) 
  reac_source_local(10,161) = - reac_rate_local(161) 
  reac_source_local(10,162) = + reac_rate_local(162) 
  reac_source_local(11,162) = - reac_rate_local(162) 
  reac_source_local(01,163) = + reac_rate_local(163) 
  reac_source_local(12,163) = - reac_rate_local(163) 
  reac_source_local(11,164) = + reac_rate_local(164) 
  reac_source_local(13,164) = - reac_rate_local(164) 
  reac_source_local(21,165) = + reac_rate_local(165) 
  reac_source_local(26,165) = - reac_rate_local(165) 
  reac_source_local(26,166) = + reac_rate_local(166) 
  reac_source_local(27,166) = - reac_rate_local(166) 
  reac_source_local(21,167) = + reac_rate_local(167) 
  reac_source_local(27,167) = - reac_rate_local(167) 
  reac_source_local(21,168) = + reac_rate_local(168) 
  reac_source_local(28,168) = - reac_rate_local(168) 
  reac_source_local(10,169) = - reac_rate_local(169) 
  reac_source_local(15,169) = + reac_rate_local(169) 
  reac_source_local(29,169) = - reac_rate_local(169) 
  reac_source_local(40,169) = + reac_rate_local(169) 
  reac_source_local(01,170) = + reac_rate_local(170) 
  reac_source_local(10,170) = - reac_rate_local(170) 
  reac_source_local(29,170) = - reac_rate_local(170) 
  reac_source_local(31,170) = + reac_rate_local(170) 
  reac_source_local(01,171) = + reac_rate_local(171) 
  reac_source_local(10,171) = - reac_rate_local(171) 
  reac_source_local(01,172) = + reac_rate_local(172) 
  reac_source_local(10,172) = - reac_rate_local(172) 
  reac_source_local(14,172) = - reac_rate_local(172) 
  reac_source_local(16,172) = + reac_rate_local(172) 
  reac_source_local(01,173) = + reac_rate_local(173) 
  reac_source_local(10,173) = - reac_rate_local(173) 
  reac_source_local(21,173) = - reac_rate_local(173) 
  reac_source_local(29,173) = + reac_rate_local(173) 
  reac_source_local(30,173) = + reac_rate_local(173) 
  reac_source_local(01,174) = + reac_rate_local(174) 
  reac_source_local(10,174) = - reac_rate_local(174) 
  reac_source_local(21,174) = - reac_rate_local(174) 
  reac_source_local(26,174) = + reac_rate_local(174) 
  reac_source_local(01,175) = + reac_rate_local(175) 
  reac_source_local(10,175) = - reac_rate_local(175) 
  reac_source_local(21,175) = - reac_rate_local(175) 
  reac_source_local(27,175) = + reac_rate_local(175) 
  reac_source_local(01,176) = + reac_rate_local(176) 
  reac_source_local(10,176) = - reac_rate_local(176) 
  reac_source_local(01,177) = + reac_rate_local(177) 
  reac_source_local(10,177) = - reac_rate_local(177) 
  reac_source_local(01,178) = + reac_rate_local(178) 
  reac_source_local(10,178) = - reac_rate_local(178) * 2.d0
  reac_source_local(11,178) = + reac_rate_local(178) 
  reac_source_local(01,179) = + reac_rate_local(179) 
  reac_source_local(10,179) = - reac_rate_local(179) * 2.d0
  reac_source_local(13,179) = + reac_rate_local(179) 
  reac_source_local(10,180) = + reac_rate_local(180) 
  reac_source_local(11,180) = - reac_rate_local(180) 
  reac_source_local(01,181) = + reac_rate_local(181) 
  reac_source_local(11,181) = - reac_rate_local(181) 
  reac_source_local(01,182) = + reac_rate_local(182) 
  reac_source_local(11,182) = - reac_rate_local(182) 
  reac_source_local(21,182) = - reac_rate_local(182) 
  reac_source_local(29,182) = + reac_rate_local(182) * 2.d0
  reac_source_local(10,183) = + reac_rate_local(183) 
  reac_source_local(11,183) = - reac_rate_local(183) 
  reac_source_local(12,184) = + reac_rate_local(184) 
  reac_source_local(13,184) = - reac_rate_local(184) 
  reac_source_local(01,185) = + reac_rate_local(185) 
  reac_source_local(13,185) = - reac_rate_local(185) 
  reac_source_local(21,185) = - reac_rate_local(185) 
  reac_source_local(29,185) = + reac_rate_local(185) 
  reac_source_local(31,185) = + reac_rate_local(185) 
  reac_source_local(11,186) = + reac_rate_local(186) 
  reac_source_local(12,186) = - reac_rate_local(186) 
  reac_source_local(01,187) = + reac_rate_local(187) 
  reac_source_local(12,187) = - reac_rate_local(187) 
  reac_source_local(21,187) = - reac_rate_local(187) 
  reac_source_local(29,187) = + reac_rate_local(187) 
  reac_source_local(30,187) = + reac_rate_local(187) 
  reac_source_local(01,188) = + reac_rate_local(188) 
  reac_source_local(12,188) = - reac_rate_local(188) 
  reac_source_local(14,188) = + reac_rate_local(188) 
  reac_source_local(29,188) = + reac_rate_local(188) 
  reac_source_local(40,188) = - reac_rate_local(188) 
  reac_source_local(10,189) = - reac_rate_local(189) 
  reac_source_local(12,189) = - reac_rate_local(189) 
  reac_source_local(20,189) = + reac_rate_local(189) 
  reac_source_local(44,189) = + reac_rate_local(189) 
  reac_source_local(12,190) = - reac_rate_local(190) * 2.d0
  reac_source_local(20,190) = + reac_rate_local(190) 
  reac_source_local(44,190) = + reac_rate_local(190) 
  reac_source_local(10,191) = + reac_rate_local(191) 
  reac_source_local(14,191) = - reac_rate_local(191) * 2.d0
  reac_source_local(10,192) = + reac_rate_local(192) 
  reac_source_local(14,192) = - reac_rate_local(192) * 2.d0
  reac_source_local(10,193) = + reac_rate_local(193) 
  reac_source_local(14,193) = - reac_rate_local(193) * 2.d0
  reac_source_local(11,194) = + reac_rate_local(194) 
  reac_source_local(14,194) = - reac_rate_local(194) * 2.d0
  reac_source_local(11,195) = + reac_rate_local(195) 
  reac_source_local(14,195) = - reac_rate_local(195) * 2.d0
  reac_source_local(11,196) = + reac_rate_local(196) 
  reac_source_local(14,196) = - reac_rate_local(196) * 2.d0
  reac_source_local(14,197) = + reac_rate_local(197) 
  reac_source_local(15,197) = - reac_rate_local(197) 
  reac_source_local(29,197) = - reac_rate_local(197) 
  reac_source_local(30,197) = + reac_rate_local(197) 
  reac_source_local(15,198) = - reac_rate_local(198) 
  reac_source_local(21,198) = - reac_rate_local(198) 
  reac_source_local(29,198) = + reac_rate_local(198) 
  reac_source_local(40,198) = + reac_rate_local(198) 
  reac_source_local(01,199) = + reac_rate_local(199) 
  reac_source_local(15,199) = - reac_rate_local(199) 
  reac_source_local(29,199) = + reac_rate_local(199) 
  reac_source_local(40,199) = - reac_rate_local(199) 
  reac_source_local(14,200) = + reac_rate_local(200) 
  reac_source_local(15,200) = - reac_rate_local(200) 
  reac_source_local(14,201) = + reac_rate_local(201) 
  reac_source_local(16,201) = - reac_rate_local(201) 
  reac_source_local(14,202) = + reac_rate_local(202) 
  reac_source_local(16,202) = - reac_rate_local(202) 
  reac_source_local(15,203) = + reac_rate_local(203) 
  reac_source_local(16,203) = - reac_rate_local(203) 
  reac_source_local(14,204) = + reac_rate_local(204) 
  reac_source_local(16,204) = - reac_rate_local(204) 
  reac_source_local(15,205) = - reac_rate_local(205) 
  reac_source_local(16,205) = - reac_rate_local(205) 
  reac_source_local(18,205) = + reac_rate_local(205) 
  reac_source_local(44,205) = + reac_rate_local(205) 
  reac_source_local(16,206) = - reac_rate_local(206) 
  reac_source_local(21,206) = - reac_rate_local(206) 
  reac_source_local(29,206) = + reac_rate_local(206) 
  reac_source_local(40,206) = + reac_rate_local(206) 
  reac_source_local(10,207) = + reac_rate_local(207) 
  reac_source_local(16,207) = - reac_rate_local(207) 
  reac_source_local(29,207) = + reac_rate_local(207) 
  reac_source_local(40,207) = - reac_rate_local(207) 
  reac_source_local(21,208) = + reac_rate_local(208) 
  reac_source_local(26,208) = - reac_rate_local(208) 
  reac_source_local(14,209) = - reac_rate_local(209) 
  reac_source_local(26,209) = - reac_rate_local(209) 
  reac_source_local(29,209) = + reac_rate_local(209) 
  reac_source_local(40,209) = + reac_rate_local(209) 
  reac_source_local(21,210) = + reac_rate_local(210) 
  reac_source_local(26,210) = - reac_rate_local(210) 
  reac_source_local(21,211) = + reac_rate_local(211) 
  reac_source_local(26,211) = - reac_rate_local(211) 
  reac_source_local(21,212) = + reac_rate_local(212) 
  reac_source_local(26,212) = - reac_rate_local(212) 
  reac_source_local(21,213) = + reac_rate_local(213) * 2.d0
  reac_source_local(26,213) = - reac_rate_local(213) 
  reac_source_local(30,213) = + reac_rate_local(213) 
  reac_source_local(32,213) = - reac_rate_local(213) 
  reac_source_local(21,214) = + reac_rate_local(214) 
  reac_source_local(26,214) = - reac_rate_local(214) * 2.d0
  reac_source_local(27,214) = + reac_rate_local(214) 
  reac_source_local(21,215) = + reac_rate_local(215) 
  reac_source_local(26,215) = + reac_rate_local(215) 
  reac_source_local(29,215) = - reac_rate_local(215) 
  reac_source_local(32,215) = - reac_rate_local(215) 
  reac_source_local(26,216) = + reac_rate_local(216) 
  reac_source_local(27,216) = - reac_rate_local(216) 
  reac_source_local(21,217) = + reac_rate_local(217) 
  reac_source_local(27,217) = - reac_rate_local(217) 
  reac_source_local(29,217) = - reac_rate_local(217) 
  reac_source_local(30,217) = + reac_rate_local(217) 
  reac_source_local(26,218) = + reac_rate_local(218) 
  reac_source_local(27,218) = - reac_rate_local(218) 
  reac_source_local(26,219) = + reac_rate_local(219) 
  reac_source_local(27,219) = - reac_rate_local(219) 
  reac_source_local(26,220) = + reac_rate_local(220) 
  reac_source_local(27,220) = - reac_rate_local(220) 
  reac_source_local(21,221) = + reac_rate_local(221) * 2.d0
  reac_source_local(27,221) = - reac_rate_local(221) 
  reac_source_local(29,221) = + reac_rate_local(221) 
  reac_source_local(32,221) = - reac_rate_local(221) 
  reac_source_local(21,222) = + reac_rate_local(222) 
  reac_source_local(28,222) = - reac_rate_local(222) 
  reac_source_local(29,222) = - reac_rate_local(222) 
  reac_source_local(31,222) = + reac_rate_local(222) 
  reac_source_local(21,223) = - reac_rate_local(223) 
  reac_source_local(27,223) = + reac_rate_local(223) * 2.d0
  reac_source_local(28,223) = - reac_rate_local(223) 
  reac_source_local(27,224) = + reac_rate_local(224) 
  reac_source_local(28,224) = - reac_rate_local(224) 
  reac_source_local(26,225) = + reac_rate_local(225) 
  reac_source_local(29,225) = - reac_rate_local(225) * 2.d0
  reac_source_local(26,226) = + reac_rate_local(226) 
  reac_source_local(29,226) = - reac_rate_local(226) * 2.d0
  reac_source_local(26,227) = + reac_rate_local(227) 
  reac_source_local(29,227) = - reac_rate_local(227) * 2.d0
  reac_source_local(27,228) = + reac_rate_local(228) 
  reac_source_local(29,228) = - reac_rate_local(228) * 2.d0
  reac_source_local(21,229) = - reac_rate_local(229) 
  reac_source_local(26,229) = - reac_rate_local(229) * 2.d0
  reac_source_local(32,229) = + reac_rate_local(229) * 2.d0
  reac_source_local(29,230) = + reac_rate_local(230) 
  reac_source_local(30,230) = - reac_rate_local(230) 
  reac_source_local(29,231) = + reac_rate_local(231) 
  reac_source_local(30,231) = - reac_rate_local(231) 
  reac_source_local(21,232) = - reac_rate_local(232) 
  reac_source_local(26,232) = + reac_rate_local(232) 
  reac_source_local(29,232) = + reac_rate_local(232) 
  reac_source_local(30,232) = - reac_rate_local(232) 
  reac_source_local(21,233) = - reac_rate_local(233) 
  reac_source_local(27,233) = + reac_rate_local(233) 
  reac_source_local(29,233) = + reac_rate_local(233) 
  reac_source_local(30,233) = - reac_rate_local(233) 
  reac_source_local(29,234) = + reac_rate_local(234) 
  reac_source_local(30,234) = - reac_rate_local(234) 
  reac_source_local(21,235) = + reac_rate_local(235) 
  reac_source_local(29,235) = + reac_rate_local(235) * 2.d0
  reac_source_local(30,235) = - reac_rate_local(235) 
  reac_source_local(32,235) = - reac_rate_local(235) 
  reac_source_local(21,236) = + reac_rate_local(236) * 2.d0
  reac_source_local(30,236) = - reac_rate_local(236) 
  reac_source_local(32,236) = - reac_rate_local(236) 
  reac_source_local(14,237) = + reac_rate_local(237) 
  reac_source_local(21,237) = + reac_rate_local(237) 
  reac_source_local(30,237) = - reac_rate_local(237) 
  reac_source_local(40,237) = - reac_rate_local(237) 
  reac_source_local(30,238) = + reac_rate_local(238) 
  reac_source_local(31,238) = - reac_rate_local(238) 
  reac_source_local(29,239) = + reac_rate_local(239) 
  reac_source_local(31,239) = - reac_rate_local(239) 
  reac_source_local(30,240) = + reac_rate_local(240) 
  reac_source_local(31,240) = - reac_rate_local(240) 
  reac_source_local(21,241) = - reac_rate_local(241) 
  reac_source_local(29,241) = + reac_rate_local(241) * 3.d0
  reac_source_local(31,241) = - reac_rate_local(241) 
  reac_source_local(29,242) = + reac_rate_local(242) 
  reac_source_local(31,242) = - reac_rate_local(242) 
  reac_source_local(26,243) = - reac_rate_local(243) 
  reac_source_local(29,243) = + reac_rate_local(243) * 3.d0
  reac_source_local(31,243) = - reac_rate_local(243) 
  reac_source_local(26,244) = - reac_rate_local(244) 
  reac_source_local(27,244) = + reac_rate_local(244) 
  reac_source_local(30,244) = + reac_rate_local(244) 
  reac_source_local(31,244) = - reac_rate_local(244) 
  reac_source_local(26,245) = - reac_rate_local(245) 
  reac_source_local(29,245) = + reac_rate_local(245) * 3.d0
  reac_source_local(31,245) = - reac_rate_local(245) 
  reac_source_local(29,246) = + reac_rate_local(246) 
  reac_source_local(31,246) = - reac_rate_local(246) 
  reac_source_local(30,247) = + reac_rate_local(247) 
  reac_source_local(31,247) = - reac_rate_local(247) 
  reac_source_local(21,248) = + reac_rate_local(248) * 2.d0
  reac_source_local(31,248) = - reac_rate_local(248) 
  reac_source_local(32,248) = - reac_rate_local(248) 
  reac_source_local(21,249) = + reac_rate_local(249) 
  reac_source_local(29,249) = + reac_rate_local(249) 
  reac_source_local(30,249) = + reac_rate_local(249) 
  reac_source_local(31,249) = - reac_rate_local(249) 
  reac_source_local(32,249) = - reac_rate_local(249) 
  reac_source_local(01,250) = + reac_rate_local(250) 
  reac_source_local(14,250) = - reac_rate_local(250) 
  reac_source_local(29,250) = + reac_rate_local(250) 
  reac_source_local(40,250) = - reac_rate_local(250) 
  reac_source_local(14,251) = - reac_rate_local(251) 
  reac_source_local(21,251) = - reac_rate_local(251) 
  reac_source_local(29,251) = + reac_rate_local(251) 
  reac_source_local(40,251) = + reac_rate_local(251) 
  reac_source_local(14,252) = - reac_rate_local(252) 
  reac_source_local(21,252) = + reac_rate_local(252) 
  reac_source_local(32,252) = - reac_rate_local(252) 
  reac_source_local(40,252) = + reac_rate_local(252) 
  reac_source_local(01,253) = - reac_rate_local(253) 
  reac_source_local(14,253) = + reac_rate_local(253) 
  reac_source_local(29,253) = - reac_rate_local(253) 
  reac_source_local(40,253) = + reac_rate_local(253) 
  reac_source_local(14,254) = + reac_rate_local(254) 
  reac_source_local(21,254) = + reac_rate_local(254) 
  reac_source_local(29,254) = - reac_rate_local(254) 
  reac_source_local(40,254) = - reac_rate_local(254) 
  reac_source_local(21,255) = + reac_rate_local(255) 
  reac_source_local(26,255) = + reac_rate_local(255) 
  reac_source_local(29,255) = - reac_rate_local(255) 
  reac_source_local(32,255) = - reac_rate_local(255) 
  reac_source_local(01,256) = + reac_rate_local(256) 
  reac_source_local(21,256) = + reac_rate_local(256) 
  reac_source_local(40,256) = - reac_rate_local(256) * 2.d0
  reac_source_local(21,257) = - reac_rate_local(257) * 2.d0
  reac_source_local(29,257) = + reac_rate_local(257) 
  reac_source_local(32,257) = + reac_rate_local(257) 
  reac_source_local(14,258) = - reac_rate_local(258) * 2.d0
  reac_source_local(18,258) = + reac_rate_local(258) 
  reac_source_local(44,258) = + reac_rate_local(258) 
  reac_source_local(14,259) = - reac_rate_local(259) 
  reac_source_local(29,259) = - reac_rate_local(259) 
  reac_source_local(41,259) = + reac_rate_local(259) 
  reac_source_local(44,259) = + reac_rate_local(259) 
  reac_source_local(01,260) = - reac_rate_local(260) 
  reac_source_local(14,260) = + reac_rate_local(260) * 2.d0
  reac_source_local(01,261) = - reac_rate_local(261) 
  reac_source_local(14,261) = + reac_rate_local(261) * 2.d0
  reac_source_local(01,262) = - reac_rate_local(262) 
  reac_source_local(14,262) = + reac_rate_local(262) * 2.d0
  reac_source_local(21,263) = - reac_rate_local(263) 
  reac_source_local(29,263) = + reac_rate_local(263) * 2.d0
  reac_source_local(21,264) = - reac_rate_local(264) 
  reac_source_local(29,264) = + reac_rate_local(264) * 2.d0
  reac_source_local(21,265) = - reac_rate_local(265) 
  reac_source_local(29,265) = + reac_rate_local(265) * 2.d0
  reac_source_local(14,266) = + reac_rate_local(266) 
  reac_source_local(29,266) = + reac_rate_local(266) 
  reac_source_local(40,266) = - reac_rate_local(266) 
  reac_source_local(14,267) = + reac_rate_local(267) 
  reac_source_local(29,267) = + reac_rate_local(267) 
  reac_source_local(40,267) = - reac_rate_local(267) 
  reac_source_local(14,268) = + reac_rate_local(268) 
  reac_source_local(29,268) = + reac_rate_local(268) 
  reac_source_local(40,268) = - reac_rate_local(268) 
  reac_source_local(21,269) = + reac_rate_local(269) 
  reac_source_local(29,269) = + reac_rate_local(269) 
  reac_source_local(32,269) = - reac_rate_local(269) 
  reac_source_local(21,270) = + reac_rate_local(270) 
  reac_source_local(29,270) = + reac_rate_local(270) 
  reac_source_local(32,270) = - reac_rate_local(270) 
  reac_source_local(21,271) = + reac_rate_local(271) 
  reac_source_local(29,271) = + reac_rate_local(271) 
  reac_source_local(32,271) = - reac_rate_local(271) 
  reac_source_local(21,272) = + reac_rate_local(272) 
  reac_source_local(29,272) = + reac_rate_local(272) 
  reac_source_local(32,272) = - reac_rate_local(272) 
  reac_source_local(01,273) = + reac_rate_local(273) 
  reac_source_local(14,273) = - reac_rate_local(273) * 2.d0
  reac_source_local(01,274) = + reac_rate_local(274) 
  reac_source_local(14,274) = - reac_rate_local(274) * 2.d0
  reac_source_local(01,275) = + reac_rate_local(275) 
  reac_source_local(14,275) = - reac_rate_local(275) * 2.d0
  reac_source_local(01,276) = + reac_rate_local(276) 
  reac_source_local(14,276) = - reac_rate_local(276) * 2.d0
  reac_source_local(01,277) = + reac_rate_local(277) 
  reac_source_local(14,277) = - reac_rate_local(277) * 2.d0
  reac_source_local(21,278) = + reac_rate_local(278) 
  reac_source_local(29,278) = - reac_rate_local(278) * 2.d0
  reac_source_local(21,279) = + reac_rate_local(279) 
  reac_source_local(29,279) = - reac_rate_local(279) * 2.d0
  reac_source_local(21,280) = + reac_rate_local(280) 
  reac_source_local(29,280) = - reac_rate_local(280) * 2.d0
  reac_source_local(21,281) = + reac_rate_local(281) 
  reac_source_local(29,281) = - reac_rate_local(281) * 2.d0
  reac_source_local(21,282) = + reac_rate_local(282) 
  reac_source_local(29,282) = - reac_rate_local(282) * 2.d0
  reac_source_local(14,283) = - reac_rate_local(283) 
  reac_source_local(29,283) = - reac_rate_local(283) 
  reac_source_local(40,283) = + reac_rate_local(283) 
  reac_source_local(14,284) = - reac_rate_local(284) 
  reac_source_local(29,284) = - reac_rate_local(284) 
  reac_source_local(40,284) = + reac_rate_local(284) 
  reac_source_local(14,285) = - reac_rate_local(285) 
  reac_source_local(29,285) = - reac_rate_local(285) 
  reac_source_local(40,285) = + reac_rate_local(285) 
  reac_source_local(14,286) = - reac_rate_local(286) 
  reac_source_local(29,286) = - reac_rate_local(286) 
  reac_source_local(40,286) = + reac_rate_local(286) 
  reac_source_local(14,287) = - reac_rate_local(287) 
  reac_source_local(29,287) = - reac_rate_local(287) 
  reac_source_local(40,287) = + reac_rate_local(287) 
  reac_source_local(21,288) = - reac_rate_local(288) 
  reac_source_local(29,288) = - reac_rate_local(288) 
  reac_source_local(32,288) = + reac_rate_local(288) 
  reac_source_local(21,289) = - reac_rate_local(289) 
  reac_source_local(29,289) = - reac_rate_local(289) 
  reac_source_local(32,289) = + reac_rate_local(289) 
  reac_source_local(21,290) = - reac_rate_local(290) 
  reac_source_local(29,290) = - reac_rate_local(290) 
  reac_source_local(32,290) = + reac_rate_local(290) 
  reac_source_local(21,291) = - reac_rate_local(291) 
  reac_source_local(29,291) = - reac_rate_local(291) 
  reac_source_local(32,291) = + reac_rate_local(291) 
  reac_source_local(21,292) = - reac_rate_local(292) 
  reac_source_local(29,292) = - reac_rate_local(292) 
  reac_source_local(32,292) = + reac_rate_local(292) 
  reac_source_local(14,293) = + reac_rate_local(293) 
  reac_source_local(17,293) = - reac_rate_local(293) 
  reac_source_local(29,293) = - reac_rate_local(293) 
  reac_source_local(33,293) = + reac_rate_local(293) 
  reac_source_local(14,294) = + reac_rate_local(294) 
  reac_source_local(17,294) = - reac_rate_local(294) 
  reac_source_local(21,294) = - reac_rate_local(294) 
  reac_source_local(34,294) = + reac_rate_local(294) 
  reac_source_local(17,295) = - reac_rate_local(295) 
  reac_source_local(21,295) = - reac_rate_local(295) 
  reac_source_local(29,295) = + reac_rate_local(295) 
  reac_source_local(41,295) = + reac_rate_local(295) 
  reac_source_local(17,296) = - reac_rate_local(296) 
  reac_source_local(21,296) = - reac_rate_local(296) 
  reac_source_local(33,296) = + reac_rate_local(296) 
  reac_source_local(40,296) = + reac_rate_local(296) 
  reac_source_local(17,297) = - reac_rate_local(297) 
  reac_source_local(21,297) = + reac_rate_local(297) 
  reac_source_local(32,297) = - reac_rate_local(297) 
  reac_source_local(41,297) = + reac_rate_local(297) 
  reac_source_local(14,298) = + reac_rate_local(298) 
  reac_source_local(17,298) = - reac_rate_local(298) 
  reac_source_local(40,298) = - reac_rate_local(298) 
  reac_source_local(41,298) = + reac_rate_local(298) 
  reac_source_local(17,299) = - reac_rate_local(299) 
  reac_source_local(18,299) = + reac_rate_local(299) 
  reac_source_local(29,299) = + reac_rate_local(299) 
  reac_source_local(40,299) = - reac_rate_local(299) 
  reac_source_local(01,300) = + reac_rate_local(300) 
  reac_source_local(17,300) = - reac_rate_local(300) 
  reac_source_local(33,300) = + reac_rate_local(300) 
  reac_source_local(40,300) = - reac_rate_local(300) 
  reac_source_local(01,301) = - reac_rate_local(301) 
  reac_source_local(14,301) = + reac_rate_local(301) 
  reac_source_local(33,301) = - reac_rate_local(301) 
  reac_source_local(41,301) = + reac_rate_local(301) 
  reac_source_local(21,302) = - reac_rate_local(302) 
  reac_source_local(29,302) = + reac_rate_local(302) 
  reac_source_local(33,302) = - reac_rate_local(302) 
  reac_source_local(34,302) = + reac_rate_local(302) 
  reac_source_local(21,303) = + reac_rate_local(303) 
  reac_source_local(32,303) = - reac_rate_local(303) 
  reac_source_local(33,303) = - reac_rate_local(303) 
  reac_source_local(34,303) = + reac_rate_local(303) 
  reac_source_local(29,304) = + reac_rate_local(304) 
  reac_source_local(33,304) = - reac_rate_local(304) 
  reac_source_local(40,304) = - reac_rate_local(304) 
  reac_source_local(41,304) = + reac_rate_local(304) 
  reac_source_local(14,305) = + reac_rate_local(305) 
  reac_source_local(33,305) = - reac_rate_local(305) 
  reac_source_local(34,305) = + reac_rate_local(305) 
  reac_source_local(40,305) = - reac_rate_local(305) 
  reac_source_local(15,306) = - reac_rate_local(306) 
  reac_source_local(17,306) = + reac_rate_local(306) 
  reac_source_local(29,306) = + reac_rate_local(306) 
  reac_source_local(33,306) = - reac_rate_local(306) 
  reac_source_local(01,307) = + reac_rate_local(307) 
  reac_source_local(18,307) = - reac_rate_local(307) 
  reac_source_local(21,307) = - reac_rate_local(307) 
  reac_source_local(34,307) = + reac_rate_local(307) 
  reac_source_local(14,308) = + reac_rate_local(308) 
  reac_source_local(18,308) = - reac_rate_local(308) 
  reac_source_local(29,308) = - reac_rate_local(308) 
  reac_source_local(41,308) = + reac_rate_local(308) 
  reac_source_local(01,309) = + reac_rate_local(309) 
  reac_source_local(18,309) = - reac_rate_local(309) 
  reac_source_local(29,309) = - reac_rate_local(309) 
  reac_source_local(33,309) = + reac_rate_local(309) 
  reac_source_local(01,310) = + reac_rate_local(310) 
  reac_source_local(18,310) = - reac_rate_local(310) 
  reac_source_local(29,310) = + reac_rate_local(310) 
  reac_source_local(32,310) = - reac_rate_local(310) 
  reac_source_local(34,310) = + reac_rate_local(310) 
  reac_source_local(01,311) = + reac_rate_local(311) 
  reac_source_local(14,311) = - reac_rate_local(311) 
  reac_source_local(17,311) = + reac_rate_local(311) 
  reac_source_local(18,311) = - reac_rate_local(311) 
  reac_source_local(01,312) = + reac_rate_local(312) 
  reac_source_local(18,312) = - reac_rate_local(312) 
  reac_source_local(40,312) = - reac_rate_local(312) 
  reac_source_local(41,312) = + reac_rate_local(312) 
  reac_source_local(01,313) = - reac_rate_local(313) 
  reac_source_local(34,313) = - reac_rate_local(313) 
  reac_source_local(40,313) = + reac_rate_local(313) 
  reac_source_local(41,313) = + reac_rate_local(313) 
  reac_source_local(14,314) = - reac_rate_local(314) 
  reac_source_local(29,314) = + reac_rate_local(314) 
  reac_source_local(34,314) = - reac_rate_local(314) 
  reac_source_local(41,314) = + reac_rate_local(314) 
  reac_source_local(21,315) = + reac_rate_local(315) 
  reac_source_local(34,315) = - reac_rate_local(315) 
  reac_source_local(40,315) = - reac_rate_local(315) 
  reac_source_local(41,315) = + reac_rate_local(315) 
  reac_source_local(01,316) = + reac_rate_local(316) 
  reac_source_local(14,316) = + reac_rate_local(316) 
  reac_source_local(19,316) = - reac_rate_local(316) 
  reac_source_local(21,316) = - reac_rate_local(316) 
  reac_source_local(34,316) = + reac_rate_local(316) 
  reac_source_local(01,317) = + reac_rate_local(317) 
  reac_source_local(14,317) = - reac_rate_local(317) 
  reac_source_local(18,317) = + reac_rate_local(317) 
  reac_source_local(19,317) = - reac_rate_local(317) 
  reac_source_local(01,318) = + reac_rate_local(318) 
  reac_source_local(14,318) = + reac_rate_local(318) 
  reac_source_local(19,318) = - reac_rate_local(318) 
  reac_source_local(40,318) = - reac_rate_local(318) 
  reac_source_local(41,318) = + reac_rate_local(318) 
  reac_source_local(01,319) = + reac_rate_local(319) 
  reac_source_local(18,319) = + reac_rate_local(319) 
  reac_source_local(20,319) = - reac_rate_local(319) 
  reac_source_local(01,320) = + reac_rate_local(320) * 2.d0
  reac_source_local(20,320) = - reac_rate_local(320) 
  reac_source_local(21,320) = - reac_rate_local(320) 
  reac_source_local(34,320) = + reac_rate_local(320) 
  reac_source_local(01,321) = + reac_rate_local(321) * 2.d0
  reac_source_local(20,321) = - reac_rate_local(321) 
  reac_source_local(29,321) = - reac_rate_local(321) 
  reac_source_local(33,321) = + reac_rate_local(321) 
  reac_source_local(01,322) = + reac_rate_local(322) * 2.d0
  reac_source_local(14,322) = - reac_rate_local(322) 
  reac_source_local(17,322) = + reac_rate_local(322) 
  reac_source_local(20,322) = - reac_rate_local(322) 
  reac_source_local(01,323) = + reac_rate_local(323) * 2.d0
  reac_source_local(20,323) = - reac_rate_local(323) 
  reac_source_local(40,323) = - reac_rate_local(323) 
  reac_source_local(41,323) = + reac_rate_local(323) 
  reac_source_local(01,324) = - reac_rate_local(324) 
  reac_source_local(21,324) = + reac_rate_local(324) 
  reac_source_local(35,324) = - reac_rate_local(324) 
  reac_source_local(43,324) = + reac_rate_local(324) 
  reac_source_local(21,325) = + reac_rate_local(325) 
  reac_source_local(34,325) = + reac_rate_local(325) 
  reac_source_local(35,325) = - reac_rate_local(325) 
  reac_source_local(21,326) = + reac_rate_local(326) * 2.d0
  reac_source_local(26,326) = - reac_rate_local(326) 
  reac_source_local(34,326) = + reac_rate_local(326) 
  reac_source_local(35,326) = - reac_rate_local(326) 
  reac_source_local(21,327) = + reac_rate_local(327) * 2.d0
  reac_source_local(27,327) = - reac_rate_local(327) 
  reac_source_local(34,327) = + reac_rate_local(327) 
  reac_source_local(35,327) = - reac_rate_local(327) 
  reac_source_local(29,328) = - reac_rate_local(328) 
  reac_source_local(32,328) = + reac_rate_local(328) 
  reac_source_local(34,328) = + reac_rate_local(328) 
  reac_source_local(35,328) = - reac_rate_local(328) 
  reac_source_local(21,329) = + reac_rate_local(329) * 2.d0
  reac_source_local(35,329) = - reac_rate_local(329) 
  reac_source_local(40,329) = - reac_rate_local(329) 
  reac_source_local(41,329) = + reac_rate_local(329) 
  reac_source_local(01,330) = + reac_rate_local(330) 
  reac_source_local(34,330) = + reac_rate_local(330) 
  reac_source_local(43,330) = - reac_rate_local(330) 
  reac_source_local(01,331) = + reac_rate_local(331) 
  reac_source_local(21,331) = - reac_rate_local(331) 
  reac_source_local(35,331) = + reac_rate_local(331) 
  reac_source_local(43,331) = - reac_rate_local(331) 
  reac_source_local(01,332) = - reac_rate_local(332) 
  reac_source_local(17,332) = - reac_rate_local(332) 
  reac_source_local(19,332) = + reac_rate_local(332) 
  reac_source_local(17,333) = - reac_rate_local(333) 
  reac_source_local(29,333) = - reac_rate_local(333) 
  reac_source_local(41,333) = + reac_rate_local(333) 
  reac_source_local(14,334) = - reac_rate_local(334) 
  reac_source_local(17,334) = - reac_rate_local(334) 
  reac_source_local(18,334) = + reac_rate_local(334) 
  reac_source_local(01,335) = - reac_rate_local(335) 
  reac_source_local(14,335) = + reac_rate_local(335) 
  reac_source_local(33,335) = - reac_rate_local(335) 
  reac_source_local(41,335) = + reac_rate_local(335) 
  reac_source_local(29,336) = - reac_rate_local(336) 
  reac_source_local(33,336) = - reac_rate_local(336) 
  reac_source_local(34,336) = + reac_rate_local(336) 
  reac_source_local(14,337) = - reac_rate_local(337) 
  reac_source_local(33,337) = - reac_rate_local(337) 
  reac_source_local(41,337) = + reac_rate_local(337) 
  reac_source_local(01,338) = - reac_rate_local(338) 
  reac_source_local(18,338) = - reac_rate_local(338) 
  reac_source_local(20,338) = + reac_rate_local(338) 
  reac_source_local(14,339) = - reac_rate_local(339) 
  reac_source_local(18,339) = - reac_rate_local(339) 
  reac_source_local(19,339) = + reac_rate_local(339) 
  reac_source_local(21,340) = - reac_rate_local(340) 
  reac_source_local(34,340) = - reac_rate_local(340) 
  reac_source_local(35,340) = + reac_rate_local(340) 
  reac_source_local(01,341) = - reac_rate_local(341) 
  reac_source_local(34,341) = - reac_rate_local(341) 
  reac_source_local(43,341) = + reac_rate_local(341) 
  reac_source_local(26,342) = - reac_rate_local(342) 
  reac_source_local(29,342) = + reac_rate_local(342) 
  reac_source_local(36,342) = - reac_rate_local(342) 
  reac_source_local(37,342) = + reac_rate_local(342) 
  reac_source_local(29,343) = + reac_rate_local(343) 
  reac_source_local(32,343) = - reac_rate_local(343) 
  reac_source_local(36,343) = - reac_rate_local(343) 
  reac_source_local(38,343) = + reac_rate_local(343) 
  reac_source_local(21,344) = + reac_rate_local(344) 
  reac_source_local(29,344) = - reac_rate_local(344) 
  reac_source_local(36,344) = + reac_rate_local(344) 
  reac_source_local(37,344) = - reac_rate_local(344) 
  reac_source_local(21,345) = + reac_rate_local(345) 
  reac_source_local(32,345) = - reac_rate_local(345) 
  reac_source_local(37,345) = - reac_rate_local(345) 
  reac_source_local(38,345) = + reac_rate_local(345) 
  reac_source_local(21,346) = + reac_rate_local(346) 
  reac_source_local(29,346) = - reac_rate_local(346) 
  reac_source_local(37,346) = + reac_rate_local(346) 
  reac_source_local(38,346) = - reac_rate_local(346) 
  reac_source_local(21,347) = - reac_rate_local(347) 
  reac_source_local(37,347) = + reac_rate_local(347) 
  reac_source_local(40,347) = + reac_rate_local(347) 
  reac_source_local(42,347) = - reac_rate_local(347) 
  reac_source_local(21,348) = + reac_rate_local(348) 
  reac_source_local(37,348) = + reac_rate_local(348) 
  reac_source_local(39,348) = - reac_rate_local(348) 
  reac_source_local(21,349) = + reac_rate_local(349) 
  reac_source_local(29,349) = - reac_rate_local(349) 
  reac_source_local(38,349) = + reac_rate_local(349) 
  reac_source_local(39,349) = - reac_rate_local(349) 
  reac_source_local(21,350) = + reac_rate_local(350) * 2.d0
  reac_source_local(29,350) = - reac_rate_local(350) 
  reac_source_local(36,350) = + reac_rate_local(350) 
  reac_source_local(39,350) = - reac_rate_local(350) 
  reac_source_local(21,351) = + reac_rate_local(351) * 2.d0
  reac_source_local(26,351) = - reac_rate_local(351) 
  reac_source_local(37,351) = + reac_rate_local(351) 
  reac_source_local(39,351) = - reac_rate_local(351) 
  reac_source_local(21,352) = + reac_rate_local(352) * 2.d0
  reac_source_local(27,352) = - reac_rate_local(352) 
  reac_source_local(37,352) = + reac_rate_local(352) 
  reac_source_local(39,352) = - reac_rate_local(352) 
  reac_source_local(21,353) = - reac_rate_local(353) 
  reac_source_local(36,353) = - reac_rate_local(353) 
  reac_source_local(38,353) = + reac_rate_local(353) 
  reac_source_local(21,354) = - reac_rate_local(354) 
  reac_source_local(37,354) = - reac_rate_local(354) 
  reac_source_local(39,354) = + reac_rate_local(354) 
  reac_source_local(14,355) = + reac_rate_local(355) 
  reac_source_local(17,355) = - reac_rate_local(355) 
  reac_source_local(29,355) = + reac_rate_local(355) 
  reac_source_local(36,355) = - reac_rate_local(355) 
  reac_source_local(01,356) = + reac_rate_local(356) 
  reac_source_local(18,356) = - reac_rate_local(356) 
  reac_source_local(29,356) = + reac_rate_local(356) 
  reac_source_local(36,356) = - reac_rate_local(356) 
  reac_source_local(29,357) = + reac_rate_local(357) * 2.d0
  reac_source_local(33,357) = - reac_rate_local(357) 
  reac_source_local(36,357) = - reac_rate_local(357) 
  reac_source_local(21,358) = + reac_rate_local(358) 
  reac_source_local(29,358) = + reac_rate_local(358) 
  reac_source_local(34,358) = - reac_rate_local(358) 
  reac_source_local(36,358) = - reac_rate_local(358) 
  reac_source_local(29,359) = + reac_rate_local(359) 
  reac_source_local(36,359) = - reac_rate_local(359) 
  reac_source_local(40,359) = + reac_rate_local(359) 
  reac_source_local(41,359) = - reac_rate_local(359) 
  reac_source_local(14,360) = + reac_rate_local(360) 
  reac_source_local(17,360) = - reac_rate_local(360) 
  reac_source_local(21,360) = + reac_rate_local(360) 
  reac_source_local(37,360) = - reac_rate_local(360) 
  reac_source_local(01,361) = + reac_rate_local(361) 
  reac_source_local(18,361) = - reac_rate_local(361) 
  reac_source_local(21,361) = + reac_rate_local(361) 
  reac_source_local(37,361) = - reac_rate_local(361) 
  reac_source_local(21,362) = + reac_rate_local(362) 
  reac_source_local(29,362) = + reac_rate_local(362) 
  reac_source_local(33,362) = - reac_rate_local(362) 
  reac_source_local(37,362) = - reac_rate_local(362) 
  reac_source_local(21,363) = + reac_rate_local(363) * 2.d0
  reac_source_local(34,363) = - reac_rate_local(363) 
  reac_source_local(37,363) = - reac_rate_local(363) 
  reac_source_local(21,364) = + reac_rate_local(364) 
  reac_source_local(37,364) = - reac_rate_local(364) 
  reac_source_local(40,364) = + reac_rate_local(364) 
  reac_source_local(41,364) = - reac_rate_local(364) 
  reac_source_local(14,365) = + reac_rate_local(365) 
  reac_source_local(17,365) = - reac_rate_local(365) 
  reac_source_local(32,365) = + reac_rate_local(365) 
  reac_source_local(38,365) = - reac_rate_local(365) 
  reac_source_local(01,366) = + reac_rate_local(366) 
  reac_source_local(18,366) = - reac_rate_local(366) 
  reac_source_local(32,366) = + reac_rate_local(366) 
  reac_source_local(38,366) = - reac_rate_local(366) 
  reac_source_local(29,367) = + reac_rate_local(367) 
  reac_source_local(32,367) = + reac_rate_local(367) 
  reac_source_local(33,367) = - reac_rate_local(367) 
  reac_source_local(38,367) = - reac_rate_local(367) 
  reac_source_local(21,368) = + reac_rate_local(368) 
  reac_source_local(32,368) = + reac_rate_local(368) 
  reac_source_local(34,368) = - reac_rate_local(368) 
  reac_source_local(38,368) = - reac_rate_local(368) 
  reac_source_local(32,369) = + reac_rate_local(369) 
  reac_source_local(38,369) = - reac_rate_local(369) 
  reac_source_local(40,369) = + reac_rate_local(369) 
  reac_source_local(41,369) = - reac_rate_local(369) 
  reac_source_local(14,370) = + reac_rate_local(370) 
  reac_source_local(17,370) = - reac_rate_local(370) 
  reac_source_local(40,370) = + reac_rate_local(370) 
  reac_source_local(42,370) = - reac_rate_local(370) 
  reac_source_local(01,371) = + reac_rate_local(371) 
  reac_source_local(18,371) = - reac_rate_local(371) 
  reac_source_local(40,371) = + reac_rate_local(371) 
  reac_source_local(42,371) = - reac_rate_local(371) 
  reac_source_local(29,372) = + reac_rate_local(372) 
  reac_source_local(33,372) = - reac_rate_local(372) 
  reac_source_local(40,372) = + reac_rate_local(372) 
  reac_source_local(42,372) = - reac_rate_local(372) 
  reac_source_local(21,373) = + reac_rate_local(373) 
  reac_source_local(34,373) = - reac_rate_local(373) 
  reac_source_local(40,373) = + reac_rate_local(373) 
  reac_source_local(42,373) = - reac_rate_local(373) 
  reac_source_local(40,374) = + reac_rate_local(374) * 2.d0
  reac_source_local(41,374) = - reac_rate_local(374) 
  reac_source_local(42,374) = - reac_rate_local(374) 
  reac_source_local(14,375) = + reac_rate_local(375) * 2.d0
  reac_source_local(18,375) = - reac_rate_local(375) 
  reac_source_local(29,375) = + reac_rate_local(375) 
  reac_source_local(36,375) = - reac_rate_local(375) 
  reac_source_local(01,376) = + reac_rate_local(376) 
  reac_source_local(14,376) = + reac_rate_local(376) 
  reac_source_local(19,376) = - reac_rate_local(376) 
  reac_source_local(29,376) = + reac_rate_local(376) 
  reac_source_local(36,376) = - reac_rate_local(376) 
  reac_source_local(01,377) = + reac_rate_local(377) * 2.d0
  reac_source_local(20,377) = - reac_rate_local(377) 
  reac_source_local(29,377) = + reac_rate_local(377) 
  reac_source_local(36,377) = - reac_rate_local(377) 
  reac_source_local(29,378) = + reac_rate_local(378) * 3.d0
  reac_source_local(34,378) = - reac_rate_local(378) 
  reac_source_local(36,378) = - reac_rate_local(378) 
  reac_source_local(21,379) = + reac_rate_local(379) * 2.d0
  reac_source_local(29,379) = + reac_rate_local(379) 
  reac_source_local(35,379) = - reac_rate_local(379) 
  reac_source_local(36,379) = - reac_rate_local(379) 
  reac_source_local(14,380) = + reac_rate_local(380) 
  reac_source_local(29,380) = + reac_rate_local(380) * 2.d0
  reac_source_local(36,380) = - reac_rate_local(380) 
  reac_source_local(41,380) = - reac_rate_local(380) 
  reac_source_local(01,381) = + reac_rate_local(381) 
  reac_source_local(21,381) = + reac_rate_local(381) 
  reac_source_local(29,381) = + reac_rate_local(381) 
  reac_source_local(36,381) = - reac_rate_local(381) 
  reac_source_local(43,381) = - reac_rate_local(381) 
  reac_source_local(14,382) = + reac_rate_local(382) * 2.d0
  reac_source_local(18,382) = - reac_rate_local(382) 
  reac_source_local(21,382) = + reac_rate_local(382) 
  reac_source_local(37,382) = - reac_rate_local(382) 
  reac_source_local(01,383) = + reac_rate_local(383) 
  reac_source_local(14,383) = + reac_rate_local(383) 
  reac_source_local(19,383) = - reac_rate_local(383) 
  reac_source_local(21,383) = + reac_rate_local(383) 
  reac_source_local(37,383) = - reac_rate_local(383) 
  reac_source_local(01,384) = + reac_rate_local(384) * 2.d0
  reac_source_local(20,384) = - reac_rate_local(384) 
  reac_source_local(21,384) = + reac_rate_local(384) 
  reac_source_local(37,384) = - reac_rate_local(384) 
  reac_source_local(21,385) = + reac_rate_local(385) 
  reac_source_local(29,385) = + reac_rate_local(385) * 2.d0
  reac_source_local(34,385) = - reac_rate_local(385) 
  reac_source_local(37,385) = - reac_rate_local(385) 
  reac_source_local(21,386) = + reac_rate_local(386) * 3.d0
  reac_source_local(35,386) = - reac_rate_local(386) 
  reac_source_local(37,386) = - reac_rate_local(386) 
  reac_source_local(14,387) = + reac_rate_local(387) 
  reac_source_local(21,387) = + reac_rate_local(387) 
  reac_source_local(29,387) = + reac_rate_local(387) 
  reac_source_local(37,387) = - reac_rate_local(387) 
  reac_source_local(41,387) = - reac_rate_local(387) 
  reac_source_local(01,388) = + reac_rate_local(388) 
  reac_source_local(21,388) = + reac_rate_local(388) * 2.d0
  reac_source_local(37,388) = - reac_rate_local(388) 
  reac_source_local(43,388) = - reac_rate_local(388) 
  reac_source_local(14,389) = + reac_rate_local(389) * 2.d0
  reac_source_local(18,389) = - reac_rate_local(389) 
  reac_source_local(32,389) = + reac_rate_local(389) 
  reac_source_local(38,389) = - reac_rate_local(389) 
  reac_source_local(01,390) = + reac_rate_local(390) 
  reac_source_local(14,390) = + reac_rate_local(390) 
  reac_source_local(19,390) = - reac_rate_local(390) 
  reac_source_local(32,390) = + reac_rate_local(390) 
  reac_source_local(38,390) = - reac_rate_local(390) 
  reac_source_local(01,391) = + reac_rate_local(391) * 2.d0
  reac_source_local(20,391) = - reac_rate_local(391) 
  reac_source_local(32,391) = + reac_rate_local(391) 
  reac_source_local(38,391) = - reac_rate_local(391) 
  reac_source_local(29,392) = + reac_rate_local(392) * 2.d0
  reac_source_local(32,392) = + reac_rate_local(392) 
  reac_source_local(34,392) = - reac_rate_local(392) 
  reac_source_local(38,392) = - reac_rate_local(392) 
  reac_source_local(21,393) = + reac_rate_local(393) * 2.d0
  reac_source_local(32,393) = + reac_rate_local(393) 
  reac_source_local(35,393) = - reac_rate_local(393) 
  reac_source_local(38,393) = - reac_rate_local(393) 
  reac_source_local(14,394) = + reac_rate_local(394) 
  reac_source_local(29,394) = + reac_rate_local(394) 
  reac_source_local(32,394) = + reac_rate_local(394) 
  reac_source_local(38,394) = - reac_rate_local(394) 
  reac_source_local(41,394) = - reac_rate_local(394) 
  reac_source_local(01,395) = + reac_rate_local(395) 
  reac_source_local(21,395) = + reac_rate_local(395) 
  reac_source_local(32,395) = + reac_rate_local(395) 
  reac_source_local(38,395) = - reac_rate_local(395) 
  reac_source_local(43,395) = - reac_rate_local(395) 
  reac_source_local(14,396) = + reac_rate_local(396) * 2.d0
  reac_source_local(18,396) = - reac_rate_local(396) 
  reac_source_local(40,396) = + reac_rate_local(396) 
  reac_source_local(42,396) = - reac_rate_local(396) 
  reac_source_local(01,397) = + reac_rate_local(397) 
  reac_source_local(14,397) = + reac_rate_local(397) 
  reac_source_local(19,397) = - reac_rate_local(397) 
  reac_source_local(40,397) = + reac_rate_local(397) 
  reac_source_local(42,397) = - reac_rate_local(397) 
  reac_source_local(01,398) = + reac_rate_local(398) * 2.d0
  reac_source_local(20,398) = - reac_rate_local(398) 
  reac_source_local(40,398) = + reac_rate_local(398) 
  reac_source_local(42,398) = - reac_rate_local(398) 
  reac_source_local(29,399) = + reac_rate_local(399) * 2.d0
  reac_source_local(34,399) = - reac_rate_local(399) 
  reac_source_local(40,399) = + reac_rate_local(399) 
  reac_source_local(42,399) = - reac_rate_local(399) 
  reac_source_local(21,400) = + reac_rate_local(400) * 2.d0
  reac_source_local(35,400) = - reac_rate_local(400) 
  reac_source_local(40,400) = + reac_rate_local(400) 
  reac_source_local(42,400) = - reac_rate_local(400) 
  reac_source_local(14,401) = + reac_rate_local(401) 
  reac_source_local(29,401) = + reac_rate_local(401) 
  reac_source_local(40,401) = + reac_rate_local(401) 
  reac_source_local(41,401) = - reac_rate_local(401) 
  reac_source_local(42,401) = - reac_rate_local(401) 
  reac_source_local(01,402) = + reac_rate_local(402) 
  reac_source_local(21,402) = + reac_rate_local(402) 
  reac_source_local(40,402) = + reac_rate_local(402) 
  reac_source_local(42,402) = - reac_rate_local(402) 
  reac_source_local(43,402) = - reac_rate_local(402) 
  reac_source_local(14,403) = + reac_rate_local(403) 
  reac_source_local(17,403) = - reac_rate_local(403) 
  reac_source_local(21,403) = + reac_rate_local(403) * 2.d0
  reac_source_local(39,403) = - reac_rate_local(403) 
  reac_source_local(01,404) = + reac_rate_local(404) 
  reac_source_local(18,404) = - reac_rate_local(404) 
  reac_source_local(21,404) = + reac_rate_local(404) * 2.d0
  reac_source_local(39,404) = - reac_rate_local(404) 
  reac_source_local(01,405) = + reac_rate_local(405) 
  reac_source_local(14,405) = + reac_rate_local(405) 
  reac_source_local(19,405) = - reac_rate_local(405) 
  reac_source_local(21,405) = + reac_rate_local(405) * 2.d0
  reac_source_local(39,405) = - reac_rate_local(405) 
  reac_source_local(01,406) = + reac_rate_local(406) * 2.d0
  reac_source_local(20,406) = - reac_rate_local(406) 
  reac_source_local(21,406) = + reac_rate_local(406) * 2.d0
  reac_source_local(39,406) = - reac_rate_local(406) 
  reac_source_local(21,407) = + reac_rate_local(407) * 2.d0
  reac_source_local(29,407) = + reac_rate_local(407) 
  reac_source_local(33,407) = - reac_rate_local(407) 
  reac_source_local(39,407) = - reac_rate_local(407) 
  reac_source_local(21,408) = + reac_rate_local(408) * 3.d0
  reac_source_local(34,408) = - reac_rate_local(408) 
  reac_source_local(39,408) = - reac_rate_local(408) 
  reac_source_local(21,409) = + reac_rate_local(409) * 4.d0
  reac_source_local(35,409) = - reac_rate_local(409) 
  reac_source_local(39,409) = - reac_rate_local(409) 
  reac_source_local(21,410) = + reac_rate_local(410) * 2.d0
  reac_source_local(39,410) = - reac_rate_local(410) 
  reac_source_local(40,410) = + reac_rate_local(410) 
  reac_source_local(41,410) = - reac_rate_local(410) 
  reac_source_local(01,411) = + reac_rate_local(411) 
  reac_source_local(21,411) = + reac_rate_local(411) * 3.d0
  reac_source_local(39,411) = - reac_rate_local(411) 
  reac_source_local(43,411) = - reac_rate_local(411) 
  reac_source_local(14,412) = + reac_rate_local(412) 
  reac_source_local(17,412) = - reac_rate_local(412) 
  reac_source_local(29,412) = + reac_rate_local(412) 
  reac_source_local(36,412) = - reac_rate_local(412) 
  reac_source_local(01,413) = + reac_rate_local(413) 
  reac_source_local(18,413) = - reac_rate_local(413) 
  reac_source_local(29,413) = + reac_rate_local(413) 
  reac_source_local(36,413) = - reac_rate_local(413) 
  reac_source_local(29,414) = + reac_rate_local(414) * 2.d0
  reac_source_local(33,414) = - reac_rate_local(414) 
  reac_source_local(36,414) = - reac_rate_local(414) 
  reac_source_local(21,415) = + reac_rate_local(415) 
  reac_source_local(29,415) = + reac_rate_local(415) 
  reac_source_local(34,415) = - reac_rate_local(415) 
  reac_source_local(36,415) = - reac_rate_local(415) 
  reac_source_local(29,416) = + reac_rate_local(416) 
  reac_source_local(36,416) = - reac_rate_local(416) 
  reac_source_local(40,416) = + reac_rate_local(416) 
  reac_source_local(41,416) = - reac_rate_local(416) 
  reac_source_local(14,417) = + reac_rate_local(417) 
  reac_source_local(17,417) = - reac_rate_local(417) 
  reac_source_local(21,417) = + reac_rate_local(417) 
  reac_source_local(37,417) = - reac_rate_local(417) 
  reac_source_local(01,418) = + reac_rate_local(418) 
  reac_source_local(18,418) = - reac_rate_local(418) 
  reac_source_local(21,418) = + reac_rate_local(418) 
  reac_source_local(37,418) = - reac_rate_local(418) 
  reac_source_local(21,419) = + reac_rate_local(419) 
  reac_source_local(29,419) = + reac_rate_local(419) 
  reac_source_local(33,419) = - reac_rate_local(419) 
  reac_source_local(37,419) = - reac_rate_local(419) 
  reac_source_local(21,420) = + reac_rate_local(420) * 2.d0
  reac_source_local(34,420) = - reac_rate_local(420) 
  reac_source_local(37,420) = - reac_rate_local(420) 
  reac_source_local(21,421) = + reac_rate_local(421) 
  reac_source_local(37,421) = - reac_rate_local(421) 
  reac_source_local(40,421) = + reac_rate_local(421) 
  reac_source_local(41,421) = - reac_rate_local(421) 
  reac_source_local(17,422) = - reac_rate_local(422) 
  reac_source_local(36,422) = - reac_rate_local(422) 
  reac_source_local(40,422) = + reac_rate_local(422) 
  reac_source_local(21,423) = + reac_rate_local(423) 
  reac_source_local(33,423) = - reac_rate_local(423) 
  reac_source_local(36,423) = - reac_rate_local(423) 
  reac_source_local(32,424) = + reac_rate_local(424) 
  reac_source_local(34,424) = - reac_rate_local(424) 
  reac_source_local(36,424) = - reac_rate_local(424) 
  reac_source_local(32,425) = + reac_rate_local(425) 
  reac_source_local(33,425) = - reac_rate_local(425) 
  reac_source_local(37,425) = - reac_rate_local(425) 
  return
end subroutine ZDPlasKin_reac_source_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction source terms
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_fex(neq,t,y,ydot)
  implicit none
  integer,          intent(in)  :: neq
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: ydot(neq)
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(45)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(01) * density(44) 
  rrt(002) = rrt(002) * density(01) * density(44) 
  rrt(003) = rrt(003) * density(01) * density(44) 
  rrt(004) = rrt(004) * density(01) * density(44) 
  rrt(005) = rrt(005) * density(01) * density(44) 
  rrt(006) = rrt(006) * density(01) * density(44) 
  rrt(007) = rrt(007) * density(01) * density(44) 
  rrt(008) = rrt(008) * density(01) * density(44) 
  rrt(009) = rrt(009) * density(01) * density(44) 
  rrt(010) = rrt(010) * density(02) * density(44) 
  rrt(011) = rrt(011) * density(03) * density(44) 
  rrt(012) = rrt(012) * density(04) * density(44) 
  rrt(013) = rrt(013) * density(05) * density(44) 
  rrt(014) = rrt(014) * density(06) * density(44) 
  rrt(015) = rrt(015) * density(07) * density(44) 
  rrt(016) = rrt(016) * density(08) * density(44) 
  rrt(017) = rrt(017) * density(09) * density(44) 
  rrt(018) = rrt(018) * density(21) * density(44) 
  rrt(019) = rrt(019) * density(21) * density(44) 
  rrt(020) = rrt(020) * density(21) * density(44) 
  rrt(021) = rrt(021) * density(21) * density(44) 
  rrt(022) = rrt(022) * density(21) * density(44) 
  rrt(023) = rrt(023) * density(21) * density(44) 
  rrt(024) = rrt(024) * density(22) * density(44) 
  rrt(025) = rrt(025) * density(23) * density(44) 
  rrt(026) = rrt(026) * density(24) * density(44) 
  rrt(027) = rrt(027) * density(25) * density(44) 
  rrt(028) = rrt(028) * density(01) * density(02) 
  rrt(029) = rrt(029) * density(01) * density(03) 
  rrt(030) = rrt(030) * density(01) * density(04) 
  rrt(031) = rrt(031) * density(01) * density(05) 
  rrt(032) = rrt(032) * density(01) * density(06) 
  rrt(033) = rrt(033) * density(01) * density(07) 
  rrt(034) = rrt(034) * density(01) * density(08) 
  rrt(035) = rrt(035) * density(01) * density(09) 
  rrt(036) = rrt(036) * density(01)**2 
  rrt(037) = rrt(037) * density(01) * density(02) 
  rrt(038) = rrt(038) * density(01) * density(03) 
  rrt(039) = rrt(039) * density(01) * density(04) 
  rrt(040) = rrt(040) * density(01) * density(05) 
  rrt(041) = rrt(041) * density(01) * density(06) 
  rrt(042) = rrt(042) * density(01) * density(07) 
  rrt(043) = rrt(043) * density(01) * density(08) 
  rrt(044) = rrt(044) * density(02) * density(14) 
  rrt(045) = rrt(045) * density(03) * density(14) 
  rrt(046) = rrt(046) * density(04) * density(14) 
  rrt(047) = rrt(047) * density(05) * density(14) 
  rrt(048) = rrt(048) * density(06) * density(14) 
  rrt(049) = rrt(049) * density(07) * density(14) 
  rrt(050) = rrt(050) * density(08) * density(14) 
  rrt(051) = rrt(051) * density(09) * density(14) 
  rrt(052) = rrt(052) * density(01) * density(14) 
  rrt(053) = rrt(053) * density(02) * density(14) 
  rrt(054) = rrt(054) * density(03) * density(14) 
  rrt(055) = rrt(055) * density(04) * density(14) 
  rrt(056) = rrt(056) * density(05) * density(14) 
  rrt(057) = rrt(057) * density(06) * density(14) 
  rrt(058) = rrt(058) * density(07) * density(14) 
  rrt(059) = rrt(059) * density(08) * density(14) 
  rrt(060) = rrt(060) * density(02) * density(29) 
  rrt(061) = rrt(061) * density(03) * density(29) 
  rrt(062) = rrt(062) * density(04) * density(29) 
  rrt(063) = rrt(063) * density(05) * density(29) 
  rrt(064) = rrt(064) * density(06) * density(29) 
  rrt(065) = rrt(065) * density(07) * density(29) 
  rrt(066) = rrt(066) * density(08) * density(29) 
  rrt(067) = rrt(067) * density(09) * density(29) 
  rrt(068) = rrt(068) * density(01) * density(29) 
  rrt(069) = rrt(069) * density(02) * density(29) 
  rrt(070) = rrt(070) * density(03) * density(29) 
  rrt(071) = rrt(071) * density(04) * density(29) 
  rrt(072) = rrt(072) * density(05) * density(29) 
  rrt(073) = rrt(073) * density(06) * density(29) 
  rrt(074) = rrt(074) * density(07) * density(29) 
  rrt(075) = rrt(075) * density(08) * density(29) 
  rrt(076) = rrt(076) * density(21) * density(22) 
  rrt(077) = rrt(077) * density(21) * density(23) 
  rrt(078) = rrt(078) * density(21) * density(24) 
  rrt(079) = rrt(079) * density(21) * density(25) 
  rrt(080) = rrt(080) * density(21)**2 
  rrt(081) = rrt(081) * density(21) * density(22) 
  rrt(082) = rrt(082) * density(21) * density(23) 
  rrt(083) = rrt(083) * density(21) * density(24) 
  rrt(084) = rrt(084) * density(22) * density(29) 
  rrt(085) = rrt(085) * density(23) * density(29) 
  rrt(086) = rrt(086) * density(24) * density(29) 
  rrt(087) = rrt(087) * density(25) * density(29) 
  rrt(088) = rrt(088) * density(21) * density(29) 
  rrt(089) = rrt(089) * density(22) * density(29) 
  rrt(090) = rrt(090) * density(23) * density(29) 
  rrt(091) = rrt(091) * density(24) * density(29) 
  rrt(092) = rrt(092) * density(01) * density(44) 
  rrt(093) = rrt(093) * density(01) * density(44) 
  rrt(094) = rrt(094) * density(01) * density(44) 
  rrt(095) = rrt(095) * density(01) * density(44) 
  rrt(096) = rrt(096) * density(01) * density(44) 
  rrt(097) = rrt(097) * density(01) * density(44) 
  rrt(098) = rrt(098) * density(01) * density(44) 
  rrt(099) = rrt(099) * density(01) * density(44) 
  rrt(100) = rrt(100) * density(01) * density(44) 
  rrt(101) = rrt(101) * density(01) * density(44) 
  rrt(102) = rrt(102) * density(01) * density(44) 
  rrt(103) = rrt(103) * density(01) * density(44) 
  rrt(104) = rrt(104) * density(01) * density(44) 
  rrt(105) = rrt(105) * density(10) * density(44) 
  rrt(106) = rrt(106) * density(14) * density(44) 
  rrt(107) = rrt(107) * density(21) * density(44) 
  rrt(108) = rrt(108) * density(21) * density(44) 
  rrt(109) = rrt(109) * density(21) * density(44) 
  rrt(110) = rrt(110) * density(21) * density(44) 
  rrt(111) = rrt(111) * density(21) * density(44) 
  rrt(112) = rrt(112) * density(21) * density(44) 
  rrt(113) = rrt(113) * density(21) * density(44) 
  rrt(114) = rrt(114) * density(21) * density(44) 
  rrt(115) = rrt(115) * density(26) * density(44) 
  rrt(116) = rrt(116) * density(26) * density(44) 
  rrt(117) = rrt(117) * density(26) * density(44) 
  rrt(118) = rrt(118) * density(29) * density(44) 
  rrt(119) = rrt(119) * density(29) * density(44) 
  rrt(120) = rrt(120) * density(29) * density(44) 
  rrt(121) = rrt(121) * density(40) * density(44) 
  rrt(122) = rrt(122) * density(18) * density(44) 
  rrt(123) = rrt(123) * density(18) * density(44) 
  rrt(124) = rrt(124) * density(18) * density(44) 
  rrt(125) = rrt(125) * density(19) * density(44) 
  rrt(126) = rrt(126) * density(20) * density(44) 
  rrt(127) = rrt(127) * density(34) * density(44) 
  rrt(128) = rrt(128) * density(34) * density(44) 
  rrt(129) = rrt(129) * density(34) * density(44) 
  rrt(130) = rrt(130) * density(35) * density(44) 
  rrt(131) = rrt(131) * density(41) * density(44) 
  rrt(132) = rrt(132) * density(41) * density(44) 
  rrt(133) = rrt(133) * density(43) * density(44) 
  rrt(134) = rrt(134) * density(17) * density(44)**2 
  rrt(135) = rrt(135) * density(33) * density(44)**2 
  rrt(136) = rrt(136) * density(17) * density(44) 
  rrt(137) = rrt(137) * density(33) * density(44) 
  rrt(138) = rrt(138) * density(32) * density(44) 
  rrt(139) = rrt(139) * density(32) * density(44) 
  rrt(140) = rrt(140) * density(21) * density(29) * density(44) 
  rrt(141) = rrt(141) * density(21) * density(29) * density(44) 
  rrt(142) = rrt(142) * density(21) * density(32) * density(44) 
  rrt(143) = rrt(143) * density(01) * density(21) * density(44) 
  rrt(144) = rrt(144) * density(40) * density(44) 
  rrt(145) = rrt(145) * density(29) * density(36) 
  rrt(146) = rrt(146) * density(14) * density(36) 
  rrt(147) = rrt(147) * density(21) * density(36) 
  rrt(148) = rrt(148) * density(26) * density(36) 
  rrt(149) = rrt(149) * density(27) * density(36) 
  rrt(150) = rrt(150) * density(10) * density(36) 
  rrt(151) = rrt(151) * density(11) * density(36) 
  rrt(152) = rrt(152) * density(32) * density(36) 
  rrt(153) = rrt(153) * density(29) * density(37) 
  rrt(154) = rrt(154) * density(21) * density(37) 
  rrt(155) = rrt(155) * density(26) * density(37) 
  rrt(156) = rrt(156) * density(27) * density(37) 
  rrt(157) = rrt(157) * density(01) * density(37) 
  rrt(158) = rrt(158) * density(10) * density(37) 
  rrt(159) = rrt(159) * density(11) * density(37) 
  rrt(160) = rrt(160) * density(29) * density(38) 
  rrt(161) = rrt(161) * density(10) 
  rrt(162) = rrt(162) * density(11) 
  rrt(163) = rrt(163) * density(12) 
  rrt(164) = rrt(164) * density(13) 
  rrt(165) = rrt(165) * density(26) 
  rrt(166) = rrt(166) * density(27) 
  rrt(167) = rrt(167) * density(27) 
  rrt(168) = rrt(168) * density(28) 
  rrt(169) = rrt(169) * density(10) * density(29) 
  rrt(170) = rrt(170) * density(10) * density(29) 
  rrt(171) = rrt(171) * density(10) * density(14) 
  rrt(172) = rrt(172) * density(10) * density(14) 
  rrt(173) = rrt(173) * density(10) * density(21) 
  rrt(174) = rrt(174) * density(10) * density(21) 
  rrt(175) = rrt(175) * density(10) * density(21) 
  rrt(176) = rrt(176) * density(01) * density(10) 
  rrt(177) = rrt(177) * density(10) * density(40) 
  rrt(178) = rrt(178) * density(10)**2 
  rrt(179) = rrt(179) * density(10)**2 
  rrt(180) = rrt(180) * density(01) * density(11) 
  rrt(181) = rrt(181) * density(01) * density(11) 
  rrt(182) = rrt(182) * density(11) * density(21) 
  rrt(183) = rrt(183) * density(11) * density(40) 
  rrt(184) = rrt(184) * density(01) * density(13) 
  rrt(185) = rrt(185) * density(13) * density(21) 
  rrt(186) = rrt(186) * density(01) * density(12) 
  rrt(187) = rrt(187) * density(12) * density(21) 
  rrt(188) = rrt(188) * density(12) * density(40) 
  rrt(189) = rrt(189) * density(10) * density(12) 
  rrt(190) = rrt(190) * density(12)**2 
  rrt(191) = rrt(191) * density(14)**2 
  rrt(192) = rrt(192) * density(14)**3 
  rrt(193) = rrt(193) * density(14)**2 * density(29) 
  rrt(194) = rrt(194) * density(14)**2 
  rrt(195) = rrt(195) * density(14)**3 
  rrt(196) = rrt(196) * density(14)**2 * density(29) 
  rrt(197) = rrt(197) * density(15) * density(29) 
  rrt(198) = rrt(198) * density(15) * density(21) 
  rrt(199) = rrt(199) * density(15) * density(40) 
  rrt(200) = rrt(200) * density(01) * density(15) 
  rrt(201) = rrt(201) * density(14) * density(16) 
  rrt(202) = rrt(202) * density(16) * density(29) 
  rrt(203) = rrt(203) * density(14) * density(16) 
  rrt(204) = rrt(204) * density(01) * density(16) 
  rrt(205) = rrt(205) * density(15) * density(16) 
  rrt(206) = rrt(206) * density(16) * density(21) 
  rrt(207) = rrt(207) * density(16) * density(40) 
  rrt(208) = rrt(208) * density(26) * density(29) 
  rrt(209) = rrt(209) * density(14) * density(26) 
  rrt(210) = rrt(210) * density(21) * density(26) 
  rrt(211) = rrt(211) * density(01) * density(26) 
  rrt(212) = rrt(212) * density(26) * density(40) 
  rrt(213) = rrt(213) * density(26) * density(32) 
  rrt(214) = rrt(214) * density(26)**2 
  rrt(215) = rrt(215) * density(29) * density(32) 
  rrt(216) = rrt(216) * density(27) * density(29) 
  rrt(217) = rrt(217) * density(27) * density(29) 
  rrt(218) = rrt(218) * density(21) * density(27) 
  rrt(219) = rrt(219) * density(01) * density(27) 
  rrt(220) = rrt(220) * density(27) * density(40) 
  rrt(221) = rrt(221) * density(27) * density(32) 
  rrt(222) = rrt(222) * density(28) * density(29) 
  rrt(223) = rrt(223) * density(21) * density(28) 
  rrt(224) = rrt(224) * density(01) * density(28) 
  rrt(225) = rrt(225) * density(29)**2 
  rrt(226) = rrt(226) * density(21) * density(29)**2 
  rrt(227) = rrt(227) * density(29)**3 
  rrt(228) = rrt(228) * density(29)**2 
  rrt(229) = rrt(229) * density(21) * density(26)**2 
  rrt(230) = rrt(230) * density(29) * density(30) 
  rrt(231) = rrt(231) * density(21) * density(30) 
  rrt(232) = rrt(232) * density(21) * density(30) 
  rrt(233) = rrt(233) * density(21) * density(30) 
  rrt(234) = rrt(234) * density(01) * density(30) 
  rrt(235) = rrt(235) * density(30) * density(32) 
  rrt(236) = rrt(236) * density(30) * density(32) 
  rrt(237) = rrt(237) * density(30) * density(40) 
  rrt(238) = rrt(238) * density(29) * density(31) 
  rrt(239) = rrt(239) * density(14) * density(31) 
  rrt(240) = rrt(240) * density(21) * density(31) 
  rrt(241) = rrt(241) * density(21) * density(31) 
  rrt(242) = rrt(242) * density(01) * density(31) 
  rrt(243) = rrt(243) * density(26) * density(31) 
  rrt(244) = rrt(244) * density(26) * density(31) 
  rrt(245) = rrt(245) * density(26) * density(31) 
  rrt(246) = rrt(246) * density(31) * density(40) 
  rrt(247) = rrt(247) * density(31) * density(40) 
  rrt(248) = rrt(248) * density(31) * density(32) 
  rrt(249) = rrt(249) * density(31) * density(32) 
  rrt(250) = rrt(250) * density(14) * density(40) 
  rrt(251) = rrt(251) * density(14) * density(21) 
  rrt(252) = rrt(252) * density(14) * density(32) 
  rrt(253) = rrt(253) * density(01) * density(29) 
  rrt(254) = rrt(254) * density(29) * density(40) 
  rrt(255) = rrt(255) * density(29) * density(32) 
  rrt(256) = rrt(256) * density(40)**2 
  rrt(257) = rrt(257) * density(21)**2 
  rrt(258) = rrt(258) * density(14)**2 
  rrt(259) = rrt(259) * density(14) * density(29) 
  rrt(260) = rrt(260) * density(01) 
  rrt(261) = rrt(261) * density(01) * density(14) 
  rrt(262) = rrt(262) * density(01) * density(29) 
  rrt(263) = rrt(263) * density(21) 
  rrt(264) = rrt(264) * density(21)**2 
  rrt(265) = rrt(265) * density(21) * density(29) 
  rrt(266) = rrt(266) * density(40) 
  rrt(267) = rrt(267) * density(14) * density(40) 
  rrt(268) = rrt(268) * density(29) * density(40) 
  rrt(269) = rrt(269) * density(01) * density(32) 
  rrt(270) = rrt(270) * density(21) * density(32) 
  rrt(271) = rrt(271) * density(14) * density(32) 
  rrt(272) = rrt(272) * density(29) * density(32) 
  rrt(273) = rrt(273) * density(01) * density(14)**2 
  rrt(274) = rrt(274) * density(14)**2 * density(21) 
  rrt(275) = rrt(275) * density(14)**2 * density(40) 
  rrt(276) = rrt(276) * density(14)**3 
  rrt(277) = rrt(277) * density(14)**2 * density(29) 
  rrt(278) = rrt(278) * density(01) * density(29)**2 
  rrt(279) = rrt(279) * density(21) * density(29)**2 
  rrt(280) = rrt(280) * density(29)**2 * density(40) 
  rrt(281) = rrt(281) * density(14) * density(29)**2 
  rrt(282) = rrt(282) * density(29)**3 
  rrt(283) = rrt(283) * density(01) * density(14) * density(29) 
  rrt(284) = rrt(284) * density(14) * density(21) * density(29) 
  rrt(285) = rrt(285) * density(14) * density(29) * density(40) 
  rrt(286) = rrt(286) * density(14)**2 * density(29) 
  rrt(287) = rrt(287) * density(14) * density(29)**2 
  rrt(288) = rrt(288) * density(01) * density(21) * density(29) 
  rrt(289) = rrt(289) * density(21)**2 * density(29) 
  rrt(290) = rrt(290) * density(21) * density(29) * density(40) 
  rrt(291) = rrt(291) * density(14) * density(21) * density(29) 
  rrt(292) = rrt(292) * density(21) * density(29)**2 
  rrt(293) = rrt(293) * density(17) * density(29) 
  rrt(294) = rrt(294) * density(17) * density(21) 
  rrt(295) = rrt(295) * density(17) * density(21) 
  rrt(296) = rrt(296) * density(17) * density(21) 
  rrt(297) = rrt(297) * density(17) * density(32) 
  rrt(298) = rrt(298) * density(17) * density(40) 
  rrt(299) = rrt(299) * density(17) * density(40) 
  rrt(300) = rrt(300) * density(17) * density(40) 
  rrt(301) = rrt(301) * density(01) * density(33) 
  rrt(302) = rrt(302) * density(21) * density(33) 
  rrt(303) = rrt(303) * density(32) * density(33) 
  rrt(304) = rrt(304) * density(33) * density(40) 
  rrt(305) = rrt(305) * density(33) * density(40) 
  rrt(306) = rrt(306) * density(15) * density(33) 
  rrt(307) = rrt(307) * density(18) * density(21) 
  rrt(308) = rrt(308) * density(18) * density(29) 
  rrt(309) = rrt(309) * density(18) * density(29) 
  rrt(310) = rrt(310) * density(18) * density(32) 
  rrt(311) = rrt(311) * density(14) * density(18) 
  rrt(312) = rrt(312) * density(18) * density(40) 
  rrt(313) = rrt(313) * density(01) * density(34) 
  rrt(314) = rrt(314) * density(14) * density(34) 
  rrt(315) = rrt(315) * density(34) * density(40) 
  rrt(316) = rrt(316) * density(19) * density(21) 
  rrt(317) = rrt(317) * density(14) * density(19) 
  rrt(318) = rrt(318) * density(19) * density(40) 
  rrt(319) = rrt(319) * density(01) * density(20) 
  rrt(320) = rrt(320) * density(20) * density(21) 
  rrt(321) = rrt(321) * density(20) * density(29) 
  rrt(322) = rrt(322) * density(14) * density(20) 
  rrt(323) = rrt(323) * density(20) * density(40) 
  rrt(324) = rrt(324) * density(01) * density(35) 
  rrt(325) = rrt(325) * density(21) * density(35) 
  rrt(326) = rrt(326) * density(26) * density(35) 
  rrt(327) = rrt(327) * density(27) * density(35) 
  rrt(328) = rrt(328) * density(29) * density(35) 
  rrt(329) = rrt(329) * density(35) * density(40) 
  rrt(330) = rrt(330) * density(01) * density(43) 
  rrt(331) = rrt(331) * density(21) * density(43) 
  rrt(332) = rrt(332) * density(01)**2 * density(17) 
  rrt(333) = rrt(333) * density(17) * density(29) 
  rrt(334) = rrt(334) * density(14) * density(17) 
  rrt(335) = rrt(335) * density(01) * density(33) 
  rrt(336) = rrt(336) * density(29) * density(33) 
  rrt(337) = rrt(337) * density(14) * density(33) 
  rrt(338) = rrt(338) * density(01)**2 * density(18) 
  rrt(339) = rrt(339) * density(01) * density(14) * density(18) 
  rrt(340) = rrt(340) * density(21)**2 * density(34) 
  rrt(341) = rrt(341) * density(01)**2 * density(34) 
  rrt(342) = rrt(342) * density(26) * density(36) 
  rrt(343) = rrt(343) * density(32) * density(36) 
  rrt(344) = rrt(344) * density(29) * density(37) 
  rrt(345) = rrt(345) * density(32) * density(37) 
  rrt(346) = rrt(346) * density(29) * density(38) 
  rrt(347) = rrt(347) * density(21) * density(42) 
  rrt(348) = rrt(348) * density(39) 
  rrt(349) = rrt(349) * density(29) * density(39) 
  rrt(350) = rrt(350) * density(29) * density(39) 
  rrt(351) = rrt(351) * density(26) * density(39) 
  rrt(352) = rrt(352) * density(27) * density(39) 
  rrt(353) = rrt(353) * density(21) * density(36) 
  rrt(354) = rrt(354) * density(21) * density(37) 
  rrt(355) = rrt(355) * density(17) * density(36) 
  rrt(356) = rrt(356) * density(18) * density(36) 
  rrt(357) = rrt(357) * density(33) * density(36) 
  rrt(358) = rrt(358) * density(34) * density(36) 
  rrt(359) = rrt(359) * density(36) * density(41) 
  rrt(360) = rrt(360) * density(17) * density(37) 
  rrt(361) = rrt(361) * density(18) * density(37) 
  rrt(362) = rrt(362) * density(33) * density(37) 
  rrt(363) = rrt(363) * density(34) * density(37) 
  rrt(364) = rrt(364) * density(37) * density(41) 
  rrt(365) = rrt(365) * density(17) * density(38) 
  rrt(366) = rrt(366) * density(18) * density(38) 
  rrt(367) = rrt(367) * density(33) * density(38) 
  rrt(368) = rrt(368) * density(34) * density(38) 
  rrt(369) = rrt(369) * density(38) * density(41) 
  rrt(370) = rrt(370) * density(17) * density(42) 
  rrt(371) = rrt(371) * density(18) * density(42) 
  rrt(372) = rrt(372) * density(33) * density(42) 
  rrt(373) = rrt(373) * density(34) * density(42) 
  rrt(374) = rrt(374) * density(41) * density(42) 
  rrt(375) = rrt(375) * density(18) * density(36) 
  rrt(376) = rrt(376) * density(19) * density(36) 
  rrt(377) = rrt(377) * density(20) * density(36) 
  rrt(378) = rrt(378) * density(34) * density(36) 
  rrt(379) = rrt(379) * density(35) * density(36) 
  rrt(380) = rrt(380) * density(36) * density(41) 
  rrt(381) = rrt(381) * density(36) * density(43) 
  rrt(382) = rrt(382) * density(18) * density(37) 
  rrt(383) = rrt(383) * density(19) * density(37) 
  rrt(384) = rrt(384) * density(20) * density(37) 
  rrt(385) = rrt(385) * density(34) * density(37) 
  rrt(386) = rrt(386) * density(35) * density(37) 
  rrt(387) = rrt(387) * density(37) * density(41) 
  rrt(388) = rrt(388) * density(37) * density(43) 
  rrt(389) = rrt(389) * density(18) * density(38) 
  rrt(390) = rrt(390) * density(19) * density(38) 
  rrt(391) = rrt(391) * density(20) * density(38) 
  rrt(392) = rrt(392) * density(34) * density(38) 
  rrt(393) = rrt(393) * density(35) * density(38) 
  rrt(394) = rrt(394) * density(38) * density(41) 
  rrt(395) = rrt(395) * density(38) * density(43) 
  rrt(396) = rrt(396) * density(18) * density(42) 
  rrt(397) = rrt(397) * density(19) * density(42) 
  rrt(398) = rrt(398) * density(20) * density(42) 
  rrt(399) = rrt(399) * density(34) * density(42) 
  rrt(400) = rrt(400) * density(35) * density(42) 
  rrt(401) = rrt(401) * density(41) * density(42) 
  rrt(402) = rrt(402) * density(42) * density(43) 
  rrt(403) = rrt(403) * density(17) * density(39) 
  rrt(404) = rrt(404) * density(18) * density(39) 
  rrt(405) = rrt(405) * density(19) * density(39) 
  rrt(406) = rrt(406) * density(20) * density(39) 
  rrt(407) = rrt(407) * density(33) * density(39) 
  rrt(408) = rrt(408) * density(34) * density(39) 
  rrt(409) = rrt(409) * density(35) * density(39) 
  rrt(410) = rrt(410) * density(39) * density(41) 
  rrt(411) = rrt(411) * density(39) * density(43) 
  rrt(412) = rrt(412) * density(17) * density(36) 
  rrt(413) = rrt(413) * density(18) * density(36) 
  rrt(414) = rrt(414) * density(33) * density(36) 
  rrt(415) = rrt(415) * density(34) * density(36) 
  rrt(416) = rrt(416) * density(36) * density(41) 
  rrt(417) = rrt(417) * density(17) * density(37) 
  rrt(418) = rrt(418) * density(18) * density(37) 
  rrt(419) = rrt(419) * density(33) * density(37) 
  rrt(420) = rrt(420) * density(34) * density(37) 
  rrt(421) = rrt(421) * density(37) * density(41) 
  rrt(422) = rrt(422) * density(17) * density(36) 
  rrt(423) = rrt(423) * density(33) * density(36) 
  rrt(424) = rrt(424) * density(34) * density(36) 
  rrt(425) = rrt(425) * density(33) * density(37) 
  ydot(01) = -rrt(001)-rrt(002)-rrt(003)-rrt(004)-rrt(005)-rrt(006)-rrt(007)-rrt(008)-rrt(009)+rrt(010)+rrt(011)+rrt(012)+rrt(013)&
             +rrt(014)+rrt(015)+rrt(016)+rrt(017)+rrt(028)-rrt(036)+rrt(044)-rrt(052)+rrt(060)-rrt(068)-rrt(092)-rrt(093)-rrt(094)&
             -rrt(095)-rrt(096)-rrt(097)-rrt(098)-rrt(099)-rrt(100)-rrt(101)-rrt(102)-rrt(103)-rrt(104)+rrt(125)+  2.d0 * rrt(126)&
             +rrt(133)+rrt(150)+rrt(151)+rrt(158)+rrt(159)+rrt(161)+rrt(163)+rrt(170)+rrt(171)+rrt(172)+rrt(173)+rrt(174)+rrt(175)&
             +rrt(176)+rrt(177)+rrt(178)+rrt(179)+rrt(181)+rrt(182)+rrt(185)+rrt(187)+rrt(188)+rrt(199)+rrt(250)-rrt(253)+rrt(256)&
             -rrt(260)-rrt(261)-rrt(262)+rrt(273)+rrt(274)+rrt(275)+rrt(276)+rrt(277)+rrt(300)-rrt(301)+rrt(307)+rrt(309)+rrt(310)&
             +rrt(311)+rrt(312)-rrt(313)+rrt(316)+rrt(317)+rrt(318)+rrt(319)+  2.d0 * rrt(320)+  2.d0 * rrt(321)+  2.d0 * rrt(322)&
             +  2.d0 * rrt(323)-rrt(324)+rrt(330)+rrt(331)-rrt(332)-rrt(335)-rrt(338)-rrt(341)+rrt(356)+rrt(361)+rrt(366)+rrt(371)&
             +rrt(376)+  2.d0 * rrt(377)+rrt(381)+rrt(383)+  2.d0 * rrt(384)+rrt(388)+rrt(390)+  2.d0 * rrt(391)+rrt(395)+rrt(397)&
             +  2.d0 * rrt(398)+rrt(402)+rrt(404)+rrt(405)+  2.d0 * rrt(406)+rrt(411)+rrt(413)+rrt(418) 
  ydot(02) = +rrt(001)+rrt(002)-rrt(010)-rrt(028)+rrt(029)+rrt(036)-rrt(037)-rrt(044)+rrt(045)+rrt(052)-rrt(053)-rrt(060)+rrt(061)&
             +rrt(068)-rrt(069) 
  ydot(03) = +rrt(003)-rrt(011)-rrt(029)+rrt(030)+rrt(037)-rrt(038)-rrt(045)+rrt(046)+rrt(053)-rrt(054)-rrt(061)+rrt(062)+rrt(069)&
             -rrt(070) 
  ydot(04) = +rrt(004)-rrt(012)-rrt(030)+rrt(031)+rrt(038)-rrt(039)-rrt(046)+rrt(047)+rrt(054)-rrt(055)-rrt(062)+rrt(063)+rrt(070)&
             -rrt(071) 
  ydot(05) = +rrt(005)-rrt(013)-rrt(031)+rrt(032)+rrt(039)-rrt(040)-rrt(047)+rrt(048)+rrt(055)-rrt(056)-rrt(063)+rrt(064)+rrt(071)&
             -rrt(072) 
  ydot(06) = +rrt(006)-rrt(014)-rrt(032)+rrt(033)+rrt(040)-rrt(041)-rrt(048)+rrt(049)+rrt(056)-rrt(057)-rrt(064)+rrt(065)+rrt(072)&
             -rrt(073) 
  ydot(07) = +rrt(007)-rrt(015)-rrt(033)+rrt(034)+rrt(041)-rrt(042)-rrt(049)+rrt(050)+rrt(057)-rrt(058)-rrt(065)+rrt(066)+rrt(073)&
             -rrt(074) 
  ydot(08) = +rrt(008)-rrt(016)-rrt(034)+rrt(035)+rrt(042)-rrt(043)-rrt(050)+rrt(051)+rrt(058)-rrt(059)-rrt(066)+rrt(067)+rrt(074)&
             -rrt(075) 
  ydot(09) = +rrt(009)-rrt(017)-rrt(035)+rrt(043)-rrt(051)+rrt(059)-rrt(067)+rrt(075) 
  ydot(10) = +rrt(092)+rrt(093)-rrt(105)-rrt(150)-rrt(158)-rrt(161)+rrt(162)-rrt(169)-rrt(170)-rrt(171)-rrt(172)-rrt(173)-rrt(174)&
             -rrt(175)-rrt(176)-rrt(177)-  2.d0 * rrt(178)-  2.d0 * rrt(179)+rrt(180)+rrt(183)-rrt(189)+rrt(191)+rrt(192)+rrt(193)&
             +rrt(207) 
  ydot(11) = +rrt(094)+rrt(095)+rrt(096)-rrt(151)-rrt(159)-rrt(162)+rrt(164)+rrt(178)-rrt(180)-rrt(181)-rrt(182)-rrt(183)+rrt(186)&
             +rrt(194)+rrt(195)+rrt(196) 
  ydot(12) = +rrt(097)+rrt(098)+rrt(099)-rrt(163)+rrt(184)-rrt(186)-rrt(187)-rrt(188)-rrt(189)-  2.d0 * rrt(190) 
  ydot(13) = +rrt(100)+rrt(101)+rrt(102)-rrt(164)+rrt(179)-rrt(184)-rrt(185) 
  ydot(14) = +rrt(103)-rrt(106)+  2.d0 * rrt(122)+rrt(123)+rrt(124)+rrt(125)+rrt(131)+rrt(134)+rrt(136)-rrt(146)-rrt(172)+rrt(188)&
             -  2.d0 * rrt(191)-  2.d0 * rrt(192)-  2.d0 * rrt(193)-  2.d0 * rrt(194)-  2.d0 * rrt(195)-  2.d0 * rrt(196)+rrt(197)&
             +rrt(200)+rrt(201)+rrt(202)+rrt(204)-rrt(209)+rrt(237)-rrt(250)-rrt(251)-rrt(252)+rrt(253)+rrt(254)-  2.d0 * rrt(258)&
             -rrt(259)+  2.d0 * rrt(260)+  2.d0 * rrt(261)+  2.d0 * rrt(262)+rrt(266)+rrt(267)+rrt(268)-  2.d0 * rrt(273)&
             -  2.d0 * rrt(274)-  2.d0 * rrt(275)-  2.d0 * rrt(276)-  2.d0 * rrt(277)-rrt(283)-rrt(284)-rrt(285)-rrt(286)-rrt(287)&
             +rrt(293)+rrt(294)+rrt(298)+rrt(301)+rrt(305)+rrt(308)-rrt(311)-rrt(314)+rrt(316)-rrt(317)+rrt(318)-rrt(322)-rrt(334)&
             +rrt(335)-rrt(337)-rrt(339)+rrt(355)+rrt(360)+rrt(365)+rrt(370)+  2.d0 * rrt(375)+rrt(376)+rrt(380)+  2.d0 * rrt(382)&
             +rrt(383)+rrt(387)+  2.d0 * rrt(389)+rrt(390)+rrt(394)+  2.d0 * rrt(396)+rrt(397)+rrt(401)+rrt(403)+rrt(405)+rrt(412)&
             +rrt(417) 
  ydot(15) = +rrt(103)+rrt(123)+rrt(132)+rrt(169)-rrt(197)-rrt(198)-rrt(199)-rrt(200)+rrt(203)-rrt(205)-rrt(306) 
  ydot(16) = +rrt(124)+rrt(172)-rrt(201)-rrt(202)-rrt(203)-rrt(204)-rrt(205)-rrt(206)-rrt(207) 
  ydot(17) = +rrt(106)-rrt(134)-rrt(136)-rrt(293)-rrt(294)-rrt(295)-rrt(296)-rrt(297)-rrt(298)-rrt(299)-rrt(300)+rrt(306)+rrt(311)&
             +rrt(322)-rrt(332)-rrt(333)-rrt(334)-rrt(355)-rrt(360)-rrt(365)-rrt(370)-rrt(403)-rrt(412)-rrt(417)-rrt(422) 
  ydot(18) = +rrt(104)+rrt(105)-rrt(122)-rrt(123)-rrt(124)+rrt(205)+rrt(258)+rrt(299)-rrt(307)-rrt(308)-rrt(309)-rrt(310)-rrt(311)&
             -rrt(312)+rrt(317)+rrt(319)+rrt(334)-rrt(338)-rrt(339)-rrt(356)-rrt(361)-rrt(366)-rrt(371)-rrt(375)-rrt(382)-rrt(389)&
             -rrt(396)-rrt(404)-rrt(413)-rrt(418) 
  ydot(19) = -rrt(125)-rrt(316)-rrt(317)-rrt(318)+rrt(332)+rrt(339)-rrt(376)-rrt(383)-rrt(390)-rrt(397)-rrt(405) 
  ydot(20) = -rrt(126)+rrt(189)+rrt(190)-rrt(319)-rrt(320)-rrt(321)-rrt(322)-rrt(323)+rrt(338)-rrt(377)-rrt(384)-rrt(391)-rrt(398)&
             -rrt(406) 
  ydot(21) = -rrt(018)-rrt(019)-rrt(020)-rrt(021)-rrt(022)-rrt(023)+rrt(024)+rrt(025)+rrt(026)+rrt(027)+rrt(076)-rrt(080)+rrt(084)&
             -rrt(088)-rrt(107)-rrt(108)-rrt(109)-rrt(110)-rrt(111)-rrt(112)-rrt(113)-rrt(114)+rrt(115)+  2.d0 * rrt(130)+rrt(133)&
             +rrt(139)-rrt(141)-rrt(143)+rrt(145)-rrt(147)+rrt(149)+  2.d0 * rrt(152)+rrt(154)+  2.d0 * rrt(155)+  2.d0 * rrt(156)&
             +rrt(157)+rrt(158)+rrt(159)+  2.d0 * rrt(160)+rrt(165)+rrt(167)+rrt(168)-rrt(173)-rrt(174)-rrt(175)-rrt(182)-rrt(185)&
             -rrt(187)-rrt(198)-rrt(206)+rrt(208)+rrt(210)+rrt(211)+rrt(212)+  2.d0 * rrt(213)+rrt(214)+rrt(215)+rrt(217)&
             +  2.d0 * rrt(221)+rrt(222)-rrt(223)-rrt(229)-rrt(232)-rrt(233)+rrt(235)+  2.d0 * rrt(236)+rrt(237)-rrt(241)&
             +  2.d0 * rrt(248)+rrt(249)-rrt(251)+rrt(252)+rrt(254)+rrt(255)+rrt(256)-  2.d0 * rrt(257)-rrt(263)-rrt(264)-rrt(265)&
             +rrt(269)+rrt(270)+rrt(271)+rrt(272)+rrt(278)+rrt(279)+rrt(280)+rrt(281)+rrt(282)-rrt(288)-rrt(289)-rrt(290)-rrt(291)&
             -rrt(292)-rrt(294)-rrt(295)-rrt(296)+rrt(297)-rrt(302)+rrt(303)-rrt(307)+rrt(315)-rrt(316)-rrt(320)+rrt(324)+rrt(325)&
             +  2.d0 * rrt(326)+  2.d0 * rrt(327)+  2.d0 * rrt(329)-rrt(331)-rrt(340)+rrt(344)+rrt(345)+rrt(346)-rrt(347)+rrt(348)&
             +rrt(349)+  2.d0 * rrt(350)+  2.d0 * rrt(351)+  2.d0 * rrt(352)-rrt(353)-rrt(354)+rrt(358)+rrt(360)+rrt(361)+rrt(362)&
             +  2.d0 * rrt(363)+rrt(364)+rrt(368)+rrt(373)+  2.d0 * rrt(379)+rrt(381)+rrt(382)+rrt(383)+rrt(384)+rrt(385)&
             +  3.d0 * rrt(386)+rrt(387)+  2.d0 * rrt(388)+  2.d0 * rrt(393)+rrt(395)+  2.d0 * rrt(400)+rrt(402)+  2.d0 * rrt(403)&
             +  2.d0 * rrt(404)+  2.d0 * rrt(405)+  2.d0 * rrt(406)+  2.d0 * rrt(407)+  3.d0 * rrt(408)+  4.d0 * rrt(409)&
             +  2.d0 * rrt(410)+  3.d0 * rrt(411)+rrt(415)+rrt(417)+rrt(418)+rrt(419)+  2.d0 * rrt(420)+rrt(421)+rrt(423) 
  ydot(22) = +rrt(018)+rrt(019)-rrt(024)-rrt(076)+rrt(077)+rrt(080)-rrt(081)-rrt(084)+rrt(085)+rrt(088)-rrt(089) 
  ydot(23) = +rrt(020)+rrt(021)-rrt(025)-rrt(077)+rrt(078)+rrt(081)-rrt(082)-rrt(085)+rrt(086)+rrt(089)-rrt(090) 
  ydot(24) = +rrt(022)-rrt(026)-rrt(078)+rrt(079)+rrt(082)-rrt(083)-rrt(086)+rrt(087)+rrt(090)-rrt(091) 
  ydot(25) = +rrt(023)-rrt(027)-rrt(079)+rrt(083)-rrt(087)+rrt(091) 
  ydot(26) = +rrt(107)-rrt(115)-rrt(116)-rrt(117)-rrt(148)-rrt(155)-rrt(165)+rrt(166)+rrt(174)-rrt(208)-rrt(209)-rrt(210)-rrt(211)&
             -rrt(212)-rrt(213)-  2.d0 * rrt(214)+rrt(215)+rrt(216)+rrt(218)+rrt(219)+rrt(220)+rrt(225)+rrt(226)+rrt(227)&
             -  2.d0 * rrt(229)+rrt(232)-rrt(243)-rrt(244)-rrt(245)+rrt(255)-rrt(326)-rrt(342)-rrt(351) 
  ydot(27) = +rrt(108)-rrt(149)-rrt(156)-rrt(166)-rrt(167)+rrt(175)+rrt(214)-rrt(216)-rrt(217)-rrt(218)-rrt(219)-rrt(220)-rrt(221)&
             +  2.d0 * rrt(223)+rrt(224)+rrt(228)+rrt(233)+rrt(244)-rrt(327)-rrt(352) 
  ydot(28) = +rrt(109)-rrt(168)-rrt(222)-rrt(223)-rrt(224) 
  ydot(29) = +  2.d0 * rrt(110)+rrt(111)+rrt(112)+rrt(114)+  2.d0 * rrt(116)-rrt(118)-rrt(119)-rrt(120)+  2.d0 * rrt(127)+rrt(128)&
             +rrt(129)+rrt(131)+rrt(132)+rrt(135)+rrt(137)+rrt(138)-rrt(140)-rrt(145)+rrt(149)+rrt(150)+rrt(151)-rrt(153)-rrt(160)&
             -rrt(169)-rrt(170)+rrt(173)+  2.d0 * rrt(182)+rrt(185)+rrt(187)+rrt(188)-rrt(197)+rrt(198)+rrt(199)+rrt(206)+rrt(207)&
             +rrt(209)-rrt(215)-rrt(217)+rrt(221)-rrt(222)-  2.d0 * rrt(225)-  2.d0 * rrt(226)-  2.d0 * rrt(227)-  2.d0 * rrt(228)&
             +rrt(230)+rrt(231)+rrt(232)+rrt(233)+rrt(234)+  2.d0 * rrt(235)+rrt(239)+  3.d0 * rrt(241)+rrt(242)+  3.d0 * rrt(243)&
             +  3.d0 * rrt(245)+rrt(246)+rrt(249)+rrt(250)+rrt(251)-rrt(253)-rrt(254)-rrt(255)+rrt(257)-rrt(259)+  2.d0 * rrt(263)&
             +  2.d0 * rrt(264)+  2.d0 * rrt(265)+rrt(266)+rrt(267)+rrt(268)+rrt(269)+rrt(270)+rrt(271)+rrt(272)-  2.d0 * rrt(278)&
             -  2.d0 * rrt(279)-  2.d0 * rrt(280)-  2.d0 * rrt(281)-  2.d0 * rrt(282)-rrt(283)-rrt(284)-rrt(285)-rrt(286)-rrt(287)&
             -rrt(288)-rrt(289)-rrt(290)-rrt(291)-rrt(292)-rrt(293)+rrt(295)+rrt(299)+rrt(302)+rrt(304)+rrt(306)-rrt(308)-rrt(309)&
             +rrt(310)+rrt(314)-rrt(321)-rrt(328)-rrt(333)-rrt(336)+rrt(342)+rrt(343)-rrt(344)-rrt(346)-rrt(349)-rrt(350)+rrt(355)&
             +rrt(356)+  2.d0 * rrt(357)+rrt(358)+rrt(359)+rrt(362)+rrt(367)+rrt(372)+rrt(375)+rrt(376)+rrt(377)+  3.d0 * rrt(378)&
             +rrt(379)+  2.d0 * rrt(380)+rrt(381)+  2.d0 * rrt(385)+rrt(387)+  2.d0 * rrt(392)+rrt(394)+  2.d0 * rrt(399)+rrt(401)&
             +rrt(407)+rrt(412)+rrt(413)+  2.d0 * rrt(414)+rrt(415)+rrt(416)+rrt(419) 
  ydot(30) = +rrt(111)+rrt(118)+rrt(128)+rrt(173)+rrt(187)+rrt(197)+rrt(213)+rrt(217)-rrt(230)-rrt(231)-rrt(232)-rrt(233)-rrt(234)&
             -rrt(235)-rrt(236)-rrt(237)+rrt(238)+rrt(240)+rrt(244)+rrt(247)+rrt(249) 
  ydot(31) = +rrt(112)+rrt(119)+rrt(129)+rrt(170)+rrt(185)+rrt(222)-rrt(238)-rrt(239)-rrt(240)-rrt(241)-rrt(242)-rrt(243)-rrt(244)&
             -rrt(245)-rrt(246)-rrt(247)-rrt(248)-rrt(249) 
  ydot(32) = -rrt(138)-rrt(139)-rrt(142)+rrt(147)+rrt(148)-rrt(152)+rrt(153)-rrt(213)-rrt(215)-rrt(221)+  2.d0 * rrt(229)-rrt(235)&
             -rrt(236)-rrt(248)-rrt(249)-rrt(252)-rrt(255)+rrt(257)-rrt(269)-rrt(270)-rrt(271)-rrt(272)+rrt(288)+rrt(289)+rrt(290)&
             +rrt(291)+rrt(292)-rrt(297)-rrt(303)-rrt(310)+rrt(328)-rrt(343)-rrt(345)+rrt(365)+rrt(366)+rrt(367)+rrt(368)+rrt(369)&
             +rrt(389)+rrt(390)+rrt(391)+rrt(392)+rrt(393)+rrt(394)+rrt(395)+rrt(424)+rrt(425) 
  ydot(33) = +rrt(120)-rrt(135)-rrt(137)+rrt(293)+rrt(296)+rrt(300)-rrt(301)-rrt(302)-rrt(303)-rrt(304)-rrt(305)-rrt(306)+rrt(309)&
             +rrt(321)-rrt(335)-rrt(336)-rrt(337)-rrt(357)-rrt(362)-rrt(367)-rrt(372)-rrt(407)-rrt(414)-rrt(419)-rrt(423)-rrt(425) 
  ydot(34) = +rrt(113)+rrt(117)-rrt(127)-rrt(128)-rrt(129)+rrt(294)+rrt(302)+rrt(303)+rrt(305)+rrt(307)+rrt(310)-rrt(313)-rrt(314)&
             -rrt(315)+rrt(316)+rrt(320)+rrt(325)+rrt(326)+rrt(327)+rrt(328)+rrt(330)+rrt(336)-rrt(340)-rrt(341)-rrt(358)-rrt(363)&
             -rrt(368)-rrt(373)-rrt(378)-rrt(385)-rrt(392)-rrt(399)-rrt(408)-rrt(415)-rrt(420)-rrt(424) 
  ydot(35) = -rrt(130)-rrt(324)-rrt(325)-rrt(326)-rrt(327)-rrt(328)-rrt(329)+rrt(331)+rrt(340)-rrt(379)-rrt(386)-rrt(393)-rrt(400)&
             -rrt(409) 
  ydot(36) = +rrt(114)+rrt(139)+rrt(140)-rrt(145)-rrt(146)-rrt(147)-rrt(148)-rrt(149)-rrt(150)-rrt(151)-rrt(152)-rrt(342)-rrt(343)&
             +rrt(344)+rrt(350)-rrt(353)-rrt(355)-rrt(356)-rrt(357)-rrt(358)-rrt(359)-rrt(375)-rrt(376)-rrt(377)-rrt(378)-rrt(379)&
             -rrt(380)-rrt(381)-rrt(412)-rrt(413)-rrt(414)-rrt(415)-rrt(416)-rrt(422)-rrt(423)-rrt(424) 
  ydot(37) = +rrt(138)+rrt(141)+rrt(143)-rrt(153)-rrt(154)-rrt(155)-rrt(156)-rrt(157)-rrt(158)-rrt(159)+rrt(342)-rrt(344)-rrt(345)&
             +rrt(346)+rrt(347)+rrt(348)+rrt(351)+rrt(352)-rrt(354)-rrt(360)-rrt(361)-rrt(362)-rrt(363)-rrt(364)-rrt(382)-rrt(383)&
             -rrt(384)-rrt(385)-rrt(386)-rrt(387)-rrt(388)-rrt(417)-rrt(418)-rrt(419)-rrt(420)-rrt(421)-rrt(425) 
  ydot(38) = +rrt(142)-rrt(160)+rrt(343)+rrt(345)-rrt(346)+rrt(349)+rrt(353)-rrt(365)-rrt(366)-rrt(367)-rrt(368)-rrt(369)-rrt(389)&
             -rrt(390)-rrt(391)-rrt(392)-rrt(393)-rrt(394)-rrt(395) 
  ydot(39) = -rrt(348)-rrt(349)-rrt(350)-rrt(351)-rrt(352)+rrt(354)-rrt(403)-rrt(404)-rrt(405)-rrt(406)-rrt(407)-rrt(408)-rrt(409)&
             -rrt(410)-rrt(411) 
  ydot(40) = -rrt(121)-rrt(144)+rrt(146)+rrt(169)-rrt(188)+rrt(198)-rrt(199)+rrt(206)-rrt(207)+rrt(209)-rrt(237)-rrt(250)+rrt(251)&
             +rrt(252)+rrt(253)-rrt(254)-  2.d0 * rrt(256)-rrt(266)-rrt(267)-rrt(268)+rrt(283)+rrt(284)+rrt(285)+rrt(286)+rrt(287)&
             +rrt(296)-rrt(298)-rrt(299)-rrt(300)-rrt(304)-rrt(305)-rrt(312)+rrt(313)-rrt(315)-rrt(318)-rrt(323)-rrt(329)+rrt(347)&
             +rrt(359)+rrt(364)+rrt(369)+rrt(370)+rrt(371)+rrt(372)+rrt(373)+  2.d0 * rrt(374)+rrt(396)+rrt(397)+rrt(398)+rrt(399)&
             +rrt(400)+rrt(401)+rrt(402)+rrt(410)+rrt(416)+rrt(421)+rrt(422) 
  ydot(41) = +rrt(121)-rrt(131)-rrt(132)+rrt(259)+rrt(295)+rrt(297)+rrt(298)+rrt(301)+rrt(304)+rrt(308)+rrt(312)+rrt(313)+rrt(314)&
             +rrt(315)+rrt(318)+rrt(323)+rrt(329)+rrt(333)+rrt(335)+rrt(337)-rrt(359)-rrt(364)-rrt(369)-rrt(374)-rrt(380)-rrt(387)&
             -rrt(394)-rrt(401)-rrt(410)-rrt(416)-rrt(421) 
  ydot(42) = +rrt(144)-rrt(347)-rrt(370)-rrt(371)-rrt(372)-rrt(373)-rrt(374)-rrt(396)-rrt(397)-rrt(398)-rrt(399)-rrt(400)-rrt(401)&
             -rrt(402) 
  ydot(43) = -rrt(133)+rrt(324)-rrt(330)-rrt(331)+rrt(341)-rrt(381)-rrt(388)-rrt(395)-rrt(402)-rrt(411) 
  ydot(44) = +rrt(104)+rrt(105)+rrt(106)+rrt(113)-rrt(114)+rrt(117)+rrt(120)+rrt(121)-rrt(122)-rrt(123)-rrt(124)-rrt(125)-rrt(126)&
             -rrt(127)-rrt(128)-rrt(129)-rrt(130)-rrt(131)-rrt(132)-rrt(133)-rrt(134)-rrt(135)-rrt(136)-rrt(137)-rrt(138)-rrt(139)&
             -rrt(140)-rrt(141)-rrt(142)-rrt(143)-rrt(144)+rrt(145)+rrt(146)+rrt(147)+rrt(148)+rrt(149)+rrt(150)+rrt(151)+rrt(152)&
             +rrt(153)+rrt(154)+rrt(155)+rrt(156)+rrt(157)+rrt(158)+rrt(159)+rrt(160)+rrt(189)+rrt(190)+rrt(205)+rrt(258)+rrt(259) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(45) = 0.0d0
  if( lgas_heating ) then
    ydot(45) = ( ZDPlasKin_cfg(14)/k_B + ydot(45) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(45) = ydot(45) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_fex
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction jacobian
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_jex(neq,t,y,ml,mu,pd,nrpd)
  implicit none
  integer,          intent(in)  :: neq, ml, mu, nrpd
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: pd(nrpd,neq)
  integer                       :: i
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(45)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(01,01) = pd(01,01) - rrt(001) * density(44) 
  pd(01,44) = pd(01,44) - rrt(001) * density(01) 
  pd(02,01) = pd(02,01) + rrt(001) * density(44) 
  pd(02,44) = pd(02,44) + rrt(001) * density(01) 
  pd(01,01) = pd(01,01) - rrt(002) * density(44) 
  pd(01,44) = pd(01,44) - rrt(002) * density(01) 
  pd(02,01) = pd(02,01) + rrt(002) * density(44) 
  pd(02,44) = pd(02,44) + rrt(002) * density(01) 
  pd(01,01) = pd(01,01) - rrt(003) * density(44) 
  pd(01,44) = pd(01,44) - rrt(003) * density(01) 
  pd(03,01) = pd(03,01) + rrt(003) * density(44) 
  pd(03,44) = pd(03,44) + rrt(003) * density(01) 
  pd(01,01) = pd(01,01) - rrt(004) * density(44) 
  pd(01,44) = pd(01,44) - rrt(004) * density(01) 
  pd(04,01) = pd(04,01) + rrt(004) * density(44) 
  pd(04,44) = pd(04,44) + rrt(004) * density(01) 
  pd(01,01) = pd(01,01) - rrt(005) * density(44) 
  pd(01,44) = pd(01,44) - rrt(005) * density(01) 
  pd(05,01) = pd(05,01) + rrt(005) * density(44) 
  pd(05,44) = pd(05,44) + rrt(005) * density(01) 
  pd(01,01) = pd(01,01) - rrt(006) * density(44) 
  pd(01,44) = pd(01,44) - rrt(006) * density(01) 
  pd(06,01) = pd(06,01) + rrt(006) * density(44) 
  pd(06,44) = pd(06,44) + rrt(006) * density(01) 
  pd(01,01) = pd(01,01) - rrt(007) * density(44) 
  pd(01,44) = pd(01,44) - rrt(007) * density(01) 
  pd(07,01) = pd(07,01) + rrt(007) * density(44) 
  pd(07,44) = pd(07,44) + rrt(007) * density(01) 
  pd(01,01) = pd(01,01) - rrt(008) * density(44) 
  pd(01,44) = pd(01,44) - rrt(008) * density(01) 
  pd(08,01) = pd(08,01) + rrt(008) * density(44) 
  pd(08,44) = pd(08,44) + rrt(008) * density(01) 
  pd(01,01) = pd(01,01) - rrt(009) * density(44) 
  pd(01,44) = pd(01,44) - rrt(009) * density(01) 
  pd(09,01) = pd(09,01) + rrt(009) * density(44) 
  pd(09,44) = pd(09,44) + rrt(009) * density(01) 
  pd(01,02) = pd(01,02) + rrt(010) * density(44) 
  pd(01,44) = pd(01,44) + rrt(010) * density(02) 
  pd(02,02) = pd(02,02) - rrt(010) * density(44) 
  pd(02,44) = pd(02,44) - rrt(010) * density(02) 
  pd(01,03) = pd(01,03) + rrt(011) * density(44) 
  pd(01,44) = pd(01,44) + rrt(011) * density(03) 
  pd(03,03) = pd(03,03) - rrt(011) * density(44) 
  pd(03,44) = pd(03,44) - rrt(011) * density(03) 
  pd(01,04) = pd(01,04) + rrt(012) * density(44) 
  pd(01,44) = pd(01,44) + rrt(012) * density(04) 
  pd(04,04) = pd(04,04) - rrt(012) * density(44) 
  pd(04,44) = pd(04,44) - rrt(012) * density(04) 
  pd(01,05) = pd(01,05) + rrt(013) * density(44) 
  pd(01,44) = pd(01,44) + rrt(013) * density(05) 
  pd(05,05) = pd(05,05) - rrt(013) * density(44) 
  pd(05,44) = pd(05,44) - rrt(013) * density(05) 
  pd(01,06) = pd(01,06) + rrt(014) * density(44) 
  pd(01,44) = pd(01,44) + rrt(014) * density(06) 
  pd(06,06) = pd(06,06) - rrt(014) * density(44) 
  pd(06,44) = pd(06,44) - rrt(014) * density(06) 
  pd(01,07) = pd(01,07) + rrt(015) * density(44) 
  pd(01,44) = pd(01,44) + rrt(015) * density(07) 
  pd(07,07) = pd(07,07) - rrt(015) * density(44) 
  pd(07,44) = pd(07,44) - rrt(015) * density(07) 
  pd(01,08) = pd(01,08) + rrt(016) * density(44) 
  pd(01,44) = pd(01,44) + rrt(016) * density(08) 
  pd(08,08) = pd(08,08) - rrt(016) * density(44) 
  pd(08,44) = pd(08,44) - rrt(016) * density(08) 
  pd(01,09) = pd(01,09) + rrt(017) * density(44) 
  pd(01,44) = pd(01,44) + rrt(017) * density(09) 
  pd(09,09) = pd(09,09) - rrt(017) * density(44) 
  pd(09,44) = pd(09,44) - rrt(017) * density(09) 
  pd(21,21) = pd(21,21) - rrt(018) * density(44) 
  pd(21,44) = pd(21,44) - rrt(018) * density(21) 
  pd(22,21) = pd(22,21) + rrt(018) * density(44) 
  pd(22,44) = pd(22,44) + rrt(018) * density(21) 
  pd(21,21) = pd(21,21) - rrt(019) * density(44) 
  pd(21,44) = pd(21,44) - rrt(019) * density(21) 
  pd(22,21) = pd(22,21) + rrt(019) * density(44) 
  pd(22,44) = pd(22,44) + rrt(019) * density(21) 
  pd(21,21) = pd(21,21) - rrt(020) * density(44) 
  pd(21,44) = pd(21,44) - rrt(020) * density(21) 
  pd(23,21) = pd(23,21) + rrt(020) * density(44) 
  pd(23,44) = pd(23,44) + rrt(020) * density(21) 
  pd(21,21) = pd(21,21) - rrt(021) * density(44) 
  pd(21,44) = pd(21,44) - rrt(021) * density(21) 
  pd(23,21) = pd(23,21) + rrt(021) * density(44) 
  pd(23,44) = pd(23,44) + rrt(021) * density(21) 
  pd(21,21) = pd(21,21) - rrt(022) * density(44) 
  pd(21,44) = pd(21,44) - rrt(022) * density(21) 
  pd(24,21) = pd(24,21) + rrt(022) * density(44) 
  pd(24,44) = pd(24,44) + rrt(022) * density(21) 
  pd(21,21) = pd(21,21) - rrt(023) * density(44) 
  pd(21,44) = pd(21,44) - rrt(023) * density(21) 
  pd(25,21) = pd(25,21) + rrt(023) * density(44) 
  pd(25,44) = pd(25,44) + rrt(023) * density(21) 
  pd(21,22) = pd(21,22) + rrt(024) * density(44) 
  pd(21,44) = pd(21,44) + rrt(024) * density(22) 
  pd(22,22) = pd(22,22) - rrt(024) * density(44) 
  pd(22,44) = pd(22,44) - rrt(024) * density(22) 
  pd(21,23) = pd(21,23) + rrt(025) * density(44) 
  pd(21,44) = pd(21,44) + rrt(025) * density(23) 
  pd(23,23) = pd(23,23) - rrt(025) * density(44) 
  pd(23,44) = pd(23,44) - rrt(025) * density(23) 
  pd(21,24) = pd(21,24) + rrt(026) * density(44) 
  pd(21,44) = pd(21,44) + rrt(026) * density(24) 
  pd(24,24) = pd(24,24) - rrt(026) * density(44) 
  pd(24,44) = pd(24,44) - rrt(026) * density(24) 
  pd(21,25) = pd(21,25) + rrt(027) * density(44) 
  pd(21,44) = pd(21,44) + rrt(027) * density(25) 
  pd(25,25) = pd(25,25) - rrt(027) * density(44) 
  pd(25,44) = pd(25,44) - rrt(027) * density(25) 
  pd(01,01) = pd(01,01) + rrt(028) * density(02) 
  pd(01,02) = pd(01,02) + rrt(028) * density(01) 
  pd(02,01) = pd(02,01) - rrt(028) * density(02) 
  pd(02,02) = pd(02,02) - rrt(028) * density(01) 
  pd(02,01) = pd(02,01) + rrt(029) * density(03) 
  pd(02,03) = pd(02,03) + rrt(029) * density(01) 
  pd(03,01) = pd(03,01) - rrt(029) * density(03) 
  pd(03,03) = pd(03,03) - rrt(029) * density(01) 
  pd(03,01) = pd(03,01) + rrt(030) * density(04) 
  pd(03,04) = pd(03,04) + rrt(030) * density(01) 
  pd(04,01) = pd(04,01) - rrt(030) * density(04) 
  pd(04,04) = pd(04,04) - rrt(030) * density(01) 
  pd(04,01) = pd(04,01) + rrt(031) * density(05) 
  pd(04,05) = pd(04,05) + rrt(031) * density(01) 
  pd(05,01) = pd(05,01) - rrt(031) * density(05) 
  pd(05,05) = pd(05,05) - rrt(031) * density(01) 
  pd(05,01) = pd(05,01) + rrt(032) * density(06) 
  pd(05,06) = pd(05,06) + rrt(032) * density(01) 
  pd(06,01) = pd(06,01) - rrt(032) * density(06) 
  pd(06,06) = pd(06,06) - rrt(032) * density(01) 
  pd(06,01) = pd(06,01) + rrt(033) * density(07) 
  pd(06,07) = pd(06,07) + rrt(033) * density(01) 
  pd(07,01) = pd(07,01) - rrt(033) * density(07) 
  pd(07,07) = pd(07,07) - rrt(033) * density(01) 
  pd(07,01) = pd(07,01) + rrt(034) * density(08) 
  pd(07,08) = pd(07,08) + rrt(034) * density(01) 
  pd(08,01) = pd(08,01) - rrt(034) * density(08) 
  pd(08,08) = pd(08,08) - rrt(034) * density(01) 
  pd(08,01) = pd(08,01) + rrt(035) * density(09) 
  pd(08,09) = pd(08,09) + rrt(035) * density(01) 
  pd(09,01) = pd(09,01) - rrt(035) * density(09) 
  pd(09,09) = pd(09,09) - rrt(035) * density(01) 
  pd(01,01) = pd(01,01) - rrt(036) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) + rrt(036) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) - rrt(037) * density(02) 
  pd(02,02) = pd(02,02) - rrt(037) * density(01) 
  pd(03,01) = pd(03,01) + rrt(037) * density(02) 
  pd(03,02) = pd(03,02) + rrt(037) * density(01) 
  pd(03,01) = pd(03,01) - rrt(038) * density(03) 
  pd(03,03) = pd(03,03) - rrt(038) * density(01) 
  pd(04,01) = pd(04,01) + rrt(038) * density(03) 
  pd(04,03) = pd(04,03) + rrt(038) * density(01) 
  pd(04,01) = pd(04,01) - rrt(039) * density(04) 
  pd(04,04) = pd(04,04) - rrt(039) * density(01) 
  pd(05,01) = pd(05,01) + rrt(039) * density(04) 
  pd(05,04) = pd(05,04) + rrt(039) * density(01) 
  pd(05,01) = pd(05,01) - rrt(040) * density(05) 
  pd(05,05) = pd(05,05) - rrt(040) * density(01) 
  pd(06,01) = pd(06,01) + rrt(040) * density(05) 
  pd(06,05) = pd(06,05) + rrt(040) * density(01) 
  pd(06,01) = pd(06,01) - rrt(041) * density(06) 
  pd(06,06) = pd(06,06) - rrt(041) * density(01) 
  pd(07,01) = pd(07,01) + rrt(041) * density(06) 
  pd(07,06) = pd(07,06) + rrt(041) * density(01) 
  pd(07,01) = pd(07,01) - rrt(042) * density(07) 
  pd(07,07) = pd(07,07) - rrt(042) * density(01) 
  pd(08,01) = pd(08,01) + rrt(042) * density(07) 
  pd(08,07) = pd(08,07) + rrt(042) * density(01) 
  pd(08,01) = pd(08,01) - rrt(043) * density(08) 
  pd(08,08) = pd(08,08) - rrt(043) * density(01) 
  pd(09,01) = pd(09,01) + rrt(043) * density(08) 
  pd(09,08) = pd(09,08) + rrt(043) * density(01) 
  pd(01,02) = pd(01,02) + rrt(044) * density(14) 
  pd(01,14) = pd(01,14) + rrt(044) * density(02) 
  pd(02,02) = pd(02,02) - rrt(044) * density(14) 
  pd(02,14) = pd(02,14) - rrt(044) * density(02) 
  pd(02,03) = pd(02,03) + rrt(045) * density(14) 
  pd(02,14) = pd(02,14) + rrt(045) * density(03) 
  pd(03,03) = pd(03,03) - rrt(045) * density(14) 
  pd(03,14) = pd(03,14) - rrt(045) * density(03) 
  pd(03,04) = pd(03,04) + rrt(046) * density(14) 
  pd(03,14) = pd(03,14) + rrt(046) * density(04) 
  pd(04,04) = pd(04,04) - rrt(046) * density(14) 
  pd(04,14) = pd(04,14) - rrt(046) * density(04) 
  pd(04,05) = pd(04,05) + rrt(047) * density(14) 
  pd(04,14) = pd(04,14) + rrt(047) * density(05) 
  pd(05,05) = pd(05,05) - rrt(047) * density(14) 
  pd(05,14) = pd(05,14) - rrt(047) * density(05) 
  pd(05,06) = pd(05,06) + rrt(048) * density(14) 
  pd(05,14) = pd(05,14) + rrt(048) * density(06) 
  pd(06,06) = pd(06,06) - rrt(048) * density(14) 
  pd(06,14) = pd(06,14) - rrt(048) * density(06) 
  pd(06,07) = pd(06,07) + rrt(049) * density(14) 
  pd(06,14) = pd(06,14) + rrt(049) * density(07) 
  pd(07,07) = pd(07,07) - rrt(049) * density(14) 
  pd(07,14) = pd(07,14) - rrt(049) * density(07) 
  pd(07,08) = pd(07,08) + rrt(050) * density(14) 
  pd(07,14) = pd(07,14) + rrt(050) * density(08) 
  pd(08,08) = pd(08,08) - rrt(050) * density(14) 
  pd(08,14) = pd(08,14) - rrt(050) * density(08) 
  pd(08,09) = pd(08,09) + rrt(051) * density(14) 
  pd(08,14) = pd(08,14) + rrt(051) * density(09) 
  pd(09,09) = pd(09,09) - rrt(051) * density(14) 
  pd(09,14) = pd(09,14) - rrt(051) * density(09) 
  pd(01,01) = pd(01,01) - rrt(052) * density(14) 
  pd(01,14) = pd(01,14) - rrt(052) * density(01) 
  pd(02,01) = pd(02,01) + rrt(052) * density(14) 
  pd(02,14) = pd(02,14) + rrt(052) * density(01) 
  pd(02,02) = pd(02,02) - rrt(053) * density(14) 
  pd(02,14) = pd(02,14) - rrt(053) * density(02) 
  pd(03,02) = pd(03,02) + rrt(053) * density(14) 
  pd(03,14) = pd(03,14) + rrt(053) * density(02) 
  pd(03,03) = pd(03,03) - rrt(054) * density(14) 
  pd(03,14) = pd(03,14) - rrt(054) * density(03) 
  pd(04,03) = pd(04,03) + rrt(054) * density(14) 
  pd(04,14) = pd(04,14) + rrt(054) * density(03) 
  pd(04,04) = pd(04,04) - rrt(055) * density(14) 
  pd(04,14) = pd(04,14) - rrt(055) * density(04) 
  pd(05,04) = pd(05,04) + rrt(055) * density(14) 
  pd(05,14) = pd(05,14) + rrt(055) * density(04) 
  pd(05,05) = pd(05,05) - rrt(056) * density(14) 
  pd(05,14) = pd(05,14) - rrt(056) * density(05) 
  pd(06,05) = pd(06,05) + rrt(056) * density(14) 
  pd(06,14) = pd(06,14) + rrt(056) * density(05) 
  pd(06,06) = pd(06,06) - rrt(057) * density(14) 
  pd(06,14) = pd(06,14) - rrt(057) * density(06) 
  pd(07,06) = pd(07,06) + rrt(057) * density(14) 
  pd(07,14) = pd(07,14) + rrt(057) * density(06) 
  pd(07,07) = pd(07,07) - rrt(058) * density(14) 
  pd(07,14) = pd(07,14) - rrt(058) * density(07) 
  pd(08,07) = pd(08,07) + rrt(058) * density(14) 
  pd(08,14) = pd(08,14) + rrt(058) * density(07) 
  pd(08,08) = pd(08,08) - rrt(059) * density(14) 
  pd(08,14) = pd(08,14) - rrt(059) * density(08) 
  pd(09,08) = pd(09,08) + rrt(059) * density(14) 
  pd(09,14) = pd(09,14) + rrt(059) * density(08) 
  pd(01,02) = pd(01,02) + rrt(060) * density(29) 
  pd(01,29) = pd(01,29) + rrt(060) * density(02) 
  pd(02,02) = pd(02,02) - rrt(060) * density(29) 
  pd(02,29) = pd(02,29) - rrt(060) * density(02) 
  pd(02,03) = pd(02,03) + rrt(061) * density(29) 
  pd(02,29) = pd(02,29) + rrt(061) * density(03) 
  pd(03,03) = pd(03,03) - rrt(061) * density(29) 
  pd(03,29) = pd(03,29) - rrt(061) * density(03) 
  pd(03,04) = pd(03,04) + rrt(062) * density(29) 
  pd(03,29) = pd(03,29) + rrt(062) * density(04) 
  pd(04,04) = pd(04,04) - rrt(062) * density(29) 
  pd(04,29) = pd(04,29) - rrt(062) * density(04) 
  pd(04,05) = pd(04,05) + rrt(063) * density(29) 
  pd(04,29) = pd(04,29) + rrt(063) * density(05) 
  pd(05,05) = pd(05,05) - rrt(063) * density(29) 
  pd(05,29) = pd(05,29) - rrt(063) * density(05) 
  pd(05,06) = pd(05,06) + rrt(064) * density(29) 
  pd(05,29) = pd(05,29) + rrt(064) * density(06) 
  pd(06,06) = pd(06,06) - rrt(064) * density(29) 
  pd(06,29) = pd(06,29) - rrt(064) * density(06) 
  pd(06,07) = pd(06,07) + rrt(065) * density(29) 
  pd(06,29) = pd(06,29) + rrt(065) * density(07) 
  pd(07,07) = pd(07,07) - rrt(065) * density(29) 
  pd(07,29) = pd(07,29) - rrt(065) * density(07) 
  pd(07,08) = pd(07,08) + rrt(066) * density(29) 
  pd(07,29) = pd(07,29) + rrt(066) * density(08) 
  pd(08,08) = pd(08,08) - rrt(066) * density(29) 
  pd(08,29) = pd(08,29) - rrt(066) * density(08) 
  pd(08,09) = pd(08,09) + rrt(067) * density(29) 
  pd(08,29) = pd(08,29) + rrt(067) * density(09) 
  pd(09,09) = pd(09,09) - rrt(067) * density(29) 
  pd(09,29) = pd(09,29) - rrt(067) * density(09) 
  pd(01,01) = pd(01,01) - rrt(068) * density(29) 
  pd(01,29) = pd(01,29) - rrt(068) * density(01) 
  pd(02,01) = pd(02,01) + rrt(068) * density(29) 
  pd(02,29) = pd(02,29) + rrt(068) * density(01) 
  pd(02,02) = pd(02,02) - rrt(069) * density(29) 
  pd(02,29) = pd(02,29) - rrt(069) * density(02) 
  pd(03,02) = pd(03,02) + rrt(069) * density(29) 
  pd(03,29) = pd(03,29) + rrt(069) * density(02) 
  pd(03,03) = pd(03,03) - rrt(070) * density(29) 
  pd(03,29) = pd(03,29) - rrt(070) * density(03) 
  pd(04,03) = pd(04,03) + rrt(070) * density(29) 
  pd(04,29) = pd(04,29) + rrt(070) * density(03) 
  pd(04,04) = pd(04,04) - rrt(071) * density(29) 
  pd(04,29) = pd(04,29) - rrt(071) * density(04) 
  pd(05,04) = pd(05,04) + rrt(071) * density(29) 
  pd(05,29) = pd(05,29) + rrt(071) * density(04) 
  pd(05,05) = pd(05,05) - rrt(072) * density(29) 
  pd(05,29) = pd(05,29) - rrt(072) * density(05) 
  pd(06,05) = pd(06,05) + rrt(072) * density(29) 
  pd(06,29) = pd(06,29) + rrt(072) * density(05) 
  pd(06,06) = pd(06,06) - rrt(073) * density(29) 
  pd(06,29) = pd(06,29) - rrt(073) * density(06) 
  pd(07,06) = pd(07,06) + rrt(073) * density(29) 
  pd(07,29) = pd(07,29) + rrt(073) * density(06) 
  pd(07,07) = pd(07,07) - rrt(074) * density(29) 
  pd(07,29) = pd(07,29) - rrt(074) * density(07) 
  pd(08,07) = pd(08,07) + rrt(074) * density(29) 
  pd(08,29) = pd(08,29) + rrt(074) * density(07) 
  pd(08,08) = pd(08,08) - rrt(075) * density(29) 
  pd(08,29) = pd(08,29) - rrt(075) * density(08) 
  pd(09,08) = pd(09,08) + rrt(075) * density(29) 
  pd(09,29) = pd(09,29) + rrt(075) * density(08) 
  pd(21,21) = pd(21,21) + rrt(076) * density(22) 
  pd(21,22) = pd(21,22) + rrt(076) * density(21) 
  pd(22,21) = pd(22,21) - rrt(076) * density(22) 
  pd(22,22) = pd(22,22) - rrt(076) * density(21) 
  pd(22,21) = pd(22,21) + rrt(077) * density(23) 
  pd(22,23) = pd(22,23) + rrt(077) * density(21) 
  pd(23,21) = pd(23,21) - rrt(077) * density(23) 
  pd(23,23) = pd(23,23) - rrt(077) * density(21) 
  pd(23,21) = pd(23,21) + rrt(078) * density(24) 
  pd(23,24) = pd(23,24) + rrt(078) * density(21) 
  pd(24,21) = pd(24,21) - rrt(078) * density(24) 
  pd(24,24) = pd(24,24) - rrt(078) * density(21) 
  pd(24,21) = pd(24,21) + rrt(079) * density(25) 
  pd(24,25) = pd(24,25) + rrt(079) * density(21) 
  pd(25,21) = pd(25,21) - rrt(079) * density(25) 
  pd(25,25) = pd(25,25) - rrt(079) * density(21) 
  pd(21,21) = pd(21,21) - rrt(080) * density(21) * 2.0d0
  pd(22,21) = pd(22,21) + rrt(080) * density(21) * 2.0d0
  pd(22,21) = pd(22,21) - rrt(081) * density(22) 
  pd(22,22) = pd(22,22) - rrt(081) * density(21) 
  pd(23,21) = pd(23,21) + rrt(081) * density(22) 
  pd(23,22) = pd(23,22) + rrt(081) * density(21) 
  pd(23,21) = pd(23,21) - rrt(082) * density(23) 
  pd(23,23) = pd(23,23) - rrt(082) * density(21) 
  pd(24,21) = pd(24,21) + rrt(082) * density(23) 
  pd(24,23) = pd(24,23) + rrt(082) * density(21) 
  pd(24,21) = pd(24,21) - rrt(083) * density(24) 
  pd(24,24) = pd(24,24) - rrt(083) * density(21) 
  pd(25,21) = pd(25,21) + rrt(083) * density(24) 
  pd(25,24) = pd(25,24) + rrt(083) * density(21) 
  pd(21,22) = pd(21,22) + rrt(084) * density(29) 
  pd(21,29) = pd(21,29) + rrt(084) * density(22) 
  pd(22,22) = pd(22,22) - rrt(084) * density(29) 
  pd(22,29) = pd(22,29) - rrt(084) * density(22) 
  pd(22,23) = pd(22,23) + rrt(085) * density(29) 
  pd(22,29) = pd(22,29) + rrt(085) * density(23) 
  pd(23,23) = pd(23,23) - rrt(085) * density(29) 
  pd(23,29) = pd(23,29) - rrt(085) * density(23) 
  pd(23,24) = pd(23,24) + rrt(086) * density(29) 
  pd(23,29) = pd(23,29) + rrt(086) * density(24) 
  pd(24,24) = pd(24,24) - rrt(086) * density(29) 
  pd(24,29) = pd(24,29) - rrt(086) * density(24) 
  pd(24,25) = pd(24,25) + rrt(087) * density(29) 
  pd(24,29) = pd(24,29) + rrt(087) * density(25) 
  pd(25,25) = pd(25,25) - rrt(087) * density(29) 
  pd(25,29) = pd(25,29) - rrt(087) * density(25) 
  pd(21,21) = pd(21,21) - rrt(088) * density(29) 
  pd(21,29) = pd(21,29) - rrt(088) * density(21) 
  pd(22,21) = pd(22,21) + rrt(088) * density(29) 
  pd(22,29) = pd(22,29) + rrt(088) * density(21) 
  pd(22,22) = pd(22,22) - rrt(089) * density(29) 
  pd(22,29) = pd(22,29) - rrt(089) * density(22) 
  pd(23,22) = pd(23,22) + rrt(089) * density(29) 
  pd(23,29) = pd(23,29) + rrt(089) * density(22) 
  pd(23,23) = pd(23,23) - rrt(090) * density(29) 
  pd(23,29) = pd(23,29) - rrt(090) * density(23) 
  pd(24,23) = pd(24,23) + rrt(090) * density(29) 
  pd(24,29) = pd(24,29) + rrt(090) * density(23) 
  pd(24,24) = pd(24,24) - rrt(091) * density(29) 
  pd(24,29) = pd(24,29) - rrt(091) * density(24) 
  pd(25,24) = pd(25,24) + rrt(091) * density(29) 
  pd(25,29) = pd(25,29) + rrt(091) * density(24) 
  pd(01,01) = pd(01,01) - rrt(092) * density(44) 
  pd(01,44) = pd(01,44) - rrt(092) * density(01) 
  pd(10,01) = pd(10,01) + rrt(092) * density(44) 
  pd(10,44) = pd(10,44) + rrt(092) * density(01) 
  pd(01,01) = pd(01,01) - rrt(093) * density(44) 
  pd(01,44) = pd(01,44) - rrt(093) * density(01) 
  pd(10,01) = pd(10,01) + rrt(093) * density(44) 
  pd(10,44) = pd(10,44) + rrt(093) * density(01) 
  pd(01,01) = pd(01,01) - rrt(094) * density(44) 
  pd(01,44) = pd(01,44) - rrt(094) * density(01) 
  pd(11,01) = pd(11,01) + rrt(094) * density(44) 
  pd(11,44) = pd(11,44) + rrt(094) * density(01) 
  pd(01,01) = pd(01,01) - rrt(095) * density(44) 
  pd(01,44) = pd(01,44) - rrt(095) * density(01) 
  pd(11,01) = pd(11,01) + rrt(095) * density(44) 
  pd(11,44) = pd(11,44) + rrt(095) * density(01) 
  pd(01,01) = pd(01,01) - rrt(096) * density(44) 
  pd(01,44) = pd(01,44) - rrt(096) * density(01) 
  pd(11,01) = pd(11,01) + rrt(096) * density(44) 
  pd(11,44) = pd(11,44) + rrt(096) * density(01) 
  pd(01,01) = pd(01,01) - rrt(097) * density(44) 
  pd(01,44) = pd(01,44) - rrt(097) * density(01) 
  pd(12,01) = pd(12,01) + rrt(097) * density(44) 
  pd(12,44) = pd(12,44) + rrt(097) * density(01) 
  pd(01,01) = pd(01,01) - rrt(098) * density(44) 
  pd(01,44) = pd(01,44) - rrt(098) * density(01) 
  pd(12,01) = pd(12,01) + rrt(098) * density(44) 
  pd(12,44) = pd(12,44) + rrt(098) * density(01) 
  pd(01,01) = pd(01,01) - rrt(099) * density(44) 
  pd(01,44) = pd(01,44) - rrt(099) * density(01) 
  pd(12,01) = pd(12,01) + rrt(099) * density(44) 
  pd(12,44) = pd(12,44) + rrt(099) * density(01) 
  pd(01,01) = pd(01,01) - rrt(100) * density(44) 
  pd(01,44) = pd(01,44) - rrt(100) * density(01) 
  pd(13,01) = pd(13,01) + rrt(100) * density(44) 
  pd(13,44) = pd(13,44) + rrt(100) * density(01) 
  pd(01,01) = pd(01,01) - rrt(101) * density(44) 
  pd(01,44) = pd(01,44) - rrt(101) * density(01) 
  pd(13,01) = pd(13,01) + rrt(101) * density(44) 
  pd(13,44) = pd(13,44) + rrt(101) * density(01) 
  pd(01,01) = pd(01,01) - rrt(102) * density(44) 
  pd(01,44) = pd(01,44) - rrt(102) * density(01) 
  pd(13,01) = pd(13,01) + rrt(102) * density(44) 
  pd(13,44) = pd(13,44) + rrt(102) * density(01) 
  pd(01,01) = pd(01,01) - rrt(103) * density(44) 
  pd(01,44) = pd(01,44) - rrt(103) * density(01) 
  pd(14,01) = pd(14,01) + rrt(103) * density(44) 
  pd(14,44) = pd(14,44) + rrt(103) * density(01) 
  pd(15,01) = pd(15,01) + rrt(103) * density(44) 
  pd(15,44) = pd(15,44) + rrt(103) * density(01) 
  pd(01,01) = pd(01,01) - rrt(104) * density(44) 
  pd(01,44) = pd(01,44) - rrt(104) * density(01) 
  pd(18,01) = pd(18,01) + rrt(104) * density(44) 
  pd(18,44) = pd(18,44) + rrt(104) * density(01) 
  pd(44,01) = pd(44,01) + rrt(104) * density(44) 
  pd(44,44) = pd(44,44) + rrt(104) * density(01) 
  pd(10,10) = pd(10,10) - rrt(105) * density(44) 
  pd(10,44) = pd(10,44) - rrt(105) * density(10) 
  pd(18,10) = pd(18,10) + rrt(105) * density(44) 
  pd(18,44) = pd(18,44) + rrt(105) * density(10) 
  pd(44,10) = pd(44,10) + rrt(105) * density(44) 
  pd(44,44) = pd(44,44) + rrt(105) * density(10) 
  pd(14,14) = pd(14,14) - rrt(106) * density(44) 
  pd(14,44) = pd(14,44) - rrt(106) * density(14) 
  pd(17,14) = pd(17,14) + rrt(106) * density(44) 
  pd(17,44) = pd(17,44) + rrt(106) * density(14) 
  pd(44,14) = pd(44,14) + rrt(106) * density(44) 
  pd(44,44) = pd(44,44) + rrt(106) * density(14) 
  pd(21,21) = pd(21,21) - rrt(107) * density(44) 
  pd(21,44) = pd(21,44) - rrt(107) * density(21) 
  pd(26,21) = pd(26,21) + rrt(107) * density(44) 
  pd(26,44) = pd(26,44) + rrt(107) * density(21) 
  pd(21,21) = pd(21,21) - rrt(108) * density(44) 
  pd(21,44) = pd(21,44) - rrt(108) * density(21) 
  pd(27,21) = pd(27,21) + rrt(108) * density(44) 
  pd(27,44) = pd(27,44) + rrt(108) * density(21) 
  pd(21,21) = pd(21,21) - rrt(109) * density(44) 
  pd(21,44) = pd(21,44) - rrt(109) * density(21) 
  pd(28,21) = pd(28,21) + rrt(109) * density(44) 
  pd(28,44) = pd(28,44) + rrt(109) * density(21) 
  pd(21,21) = pd(21,21) - rrt(110) * density(44) 
  pd(21,44) = pd(21,44) - rrt(110) * density(21) 
  pd(29,21) = pd(29,21) + rrt(110) * density(44) * 2.0d0
  pd(29,44) = pd(29,44) + rrt(110) * density(21) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(111) * density(44) 
  pd(21,44) = pd(21,44) - rrt(111) * density(21) 
  pd(29,21) = pd(29,21) + rrt(111) * density(44) 
  pd(29,44) = pd(29,44) + rrt(111) * density(21) 
  pd(30,21) = pd(30,21) + rrt(111) * density(44) 
  pd(30,44) = pd(30,44) + rrt(111) * density(21) 
  pd(21,21) = pd(21,21) - rrt(112) * density(44) 
  pd(21,44) = pd(21,44) - rrt(112) * density(21) 
  pd(29,21) = pd(29,21) + rrt(112) * density(44) 
  pd(29,44) = pd(29,44) + rrt(112) * density(21) 
  pd(31,21) = pd(31,21) + rrt(112) * density(44) 
  pd(31,44) = pd(31,44) + rrt(112) * density(21) 
  pd(21,21) = pd(21,21) - rrt(113) * density(44) 
  pd(21,44) = pd(21,44) - rrt(113) * density(21) 
  pd(34,21) = pd(34,21) + rrt(113) * density(44) 
  pd(34,44) = pd(34,44) + rrt(113) * density(21) 
  pd(44,21) = pd(44,21) + rrt(113) * density(44) 
  pd(44,44) = pd(44,44) + rrt(113) * density(21) 
  pd(21,21) = pd(21,21) - rrt(114) * density(44) 
  pd(21,44) = pd(21,44) - rrt(114) * density(21) 
  pd(29,21) = pd(29,21) + rrt(114) * density(44) 
  pd(29,44) = pd(29,44) + rrt(114) * density(21) 
  pd(36,21) = pd(36,21) + rrt(114) * density(44) 
  pd(36,44) = pd(36,44) + rrt(114) * density(21) 
  pd(44,21) = pd(44,21) - rrt(114) * density(44) 
  pd(44,44) = pd(44,44) - rrt(114) * density(21) 
  pd(21,26) = pd(21,26) + rrt(115) * density(44) 
  pd(21,44) = pd(21,44) + rrt(115) * density(26) 
  pd(26,26) = pd(26,26) - rrt(115) * density(44) 
  pd(26,44) = pd(26,44) - rrt(115) * density(26) 
  pd(26,26) = pd(26,26) - rrt(116) * density(44) 
  pd(26,44) = pd(26,44) - rrt(116) * density(26) 
  pd(29,26) = pd(29,26) + rrt(116) * density(44) * 2.0d0
  pd(29,44) = pd(29,44) + rrt(116) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(117) * density(44) 
  pd(26,44) = pd(26,44) - rrt(117) * density(26) 
  pd(34,26) = pd(34,26) + rrt(117) * density(44) 
  pd(34,44) = pd(34,44) + rrt(117) * density(26) 
  pd(44,26) = pd(44,26) + rrt(117) * density(44) 
  pd(44,44) = pd(44,44) + rrt(117) * density(26) 
  pd(29,29) = pd(29,29) - rrt(118) * density(44) 
  pd(29,44) = pd(29,44) - rrt(118) * density(29) 
  pd(30,29) = pd(30,29) + rrt(118) * density(44) 
  pd(30,44) = pd(30,44) + rrt(118) * density(29) 
  pd(29,29) = pd(29,29) - rrt(119) * density(44) 
  pd(29,44) = pd(29,44) - rrt(119) * density(29) 
  pd(31,29) = pd(31,29) + rrt(119) * density(44) 
  pd(31,44) = pd(31,44) + rrt(119) * density(29) 
  pd(29,29) = pd(29,29) - rrt(120) * density(44) 
  pd(29,44) = pd(29,44) - rrt(120) * density(29) 
  pd(33,29) = pd(33,29) + rrt(120) * density(44) 
  pd(33,44) = pd(33,44) + rrt(120) * density(29) 
  pd(44,29) = pd(44,29) + rrt(120) * density(44) 
  pd(44,44) = pd(44,44) + rrt(120) * density(29) 
  pd(40,40) = pd(40,40) - rrt(121) * density(44) 
  pd(40,44) = pd(40,44) - rrt(121) * density(40) 
  pd(41,40) = pd(41,40) + rrt(121) * density(44) 
  pd(41,44) = pd(41,44) + rrt(121) * density(40) 
  pd(44,40) = pd(44,40) + rrt(121) * density(44) 
  pd(44,44) = pd(44,44) + rrt(121) * density(40) 
  pd(14,18) = pd(14,18) + rrt(122) * density(44) * 2.0d0
  pd(14,44) = pd(14,44) + rrt(122) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(122) * density(44) 
  pd(18,44) = pd(18,44) - rrt(122) * density(18) 
  pd(44,18) = pd(44,18) - rrt(122) * density(44) 
  pd(44,44) = pd(44,44) - rrt(122) * density(18) 
  pd(14,18) = pd(14,18) + rrt(123) * density(44) 
  pd(14,44) = pd(14,44) + rrt(123) * density(18) 
  pd(15,18) = pd(15,18) + rrt(123) * density(44) 
  pd(15,44) = pd(15,44) + rrt(123) * density(18) 
  pd(18,18) = pd(18,18) - rrt(123) * density(44) 
  pd(18,44) = pd(18,44) - rrt(123) * density(18) 
  pd(44,18) = pd(44,18) - rrt(123) * density(44) 
  pd(44,44) = pd(44,44) - rrt(123) * density(18) 
  pd(14,18) = pd(14,18) + rrt(124) * density(44) 
  pd(14,44) = pd(14,44) + rrt(124) * density(18) 
  pd(16,18) = pd(16,18) + rrt(124) * density(44) 
  pd(16,44) = pd(16,44) + rrt(124) * density(18) 
  pd(18,18) = pd(18,18) - rrt(124) * density(44) 
  pd(18,44) = pd(18,44) - rrt(124) * density(18) 
  pd(44,18) = pd(44,18) - rrt(124) * density(44) 
  pd(44,44) = pd(44,44) - rrt(124) * density(18) 
  pd(01,19) = pd(01,19) + rrt(125) * density(44) 
  pd(01,44) = pd(01,44) + rrt(125) * density(19) 
  pd(14,19) = pd(14,19) + rrt(125) * density(44) 
  pd(14,44) = pd(14,44) + rrt(125) * density(19) 
  pd(19,19) = pd(19,19) - rrt(125) * density(44) 
  pd(19,44) = pd(19,44) - rrt(125) * density(19) 
  pd(44,19) = pd(44,19) - rrt(125) * density(44) 
  pd(44,44) = pd(44,44) - rrt(125) * density(19) 
  pd(01,20) = pd(01,20) + rrt(126) * density(44) * 2.0d0
  pd(01,44) = pd(01,44) + rrt(126) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(126) * density(44) 
  pd(20,44) = pd(20,44) - rrt(126) * density(20) 
  pd(44,20) = pd(44,20) - rrt(126) * density(44) 
  pd(44,44) = pd(44,44) - rrt(126) * density(20) 
  pd(29,34) = pd(29,34) + rrt(127) * density(44) * 2.0d0
  pd(29,44) = pd(29,44) + rrt(127) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(127) * density(44) 
  pd(34,44) = pd(34,44) - rrt(127) * density(34) 
  pd(44,34) = pd(44,34) - rrt(127) * density(44) 
  pd(44,44) = pd(44,44) - rrt(127) * density(34) 
  pd(29,34) = pd(29,34) + rrt(128) * density(44) 
  pd(29,44) = pd(29,44) + rrt(128) * density(34) 
  pd(30,34) = pd(30,34) + rrt(128) * density(44) 
  pd(30,44) = pd(30,44) + rrt(128) * density(34) 
  pd(34,34) = pd(34,34) - rrt(128) * density(44) 
  pd(34,44) = pd(34,44) - rrt(128) * density(34) 
  pd(44,34) = pd(44,34) - rrt(128) * density(44) 
  pd(44,44) = pd(44,44) - rrt(128) * density(34) 
  pd(29,34) = pd(29,34) + rrt(129) * density(44) 
  pd(29,44) = pd(29,44) + rrt(129) * density(34) 
  pd(31,34) = pd(31,34) + rrt(129) * density(44) 
  pd(31,44) = pd(31,44) + rrt(129) * density(34) 
  pd(34,34) = pd(34,34) - rrt(129) * density(44) 
  pd(34,44) = pd(34,44) - rrt(129) * density(34) 
  pd(44,34) = pd(44,34) - rrt(129) * density(44) 
  pd(44,44) = pd(44,44) - rrt(129) * density(34) 
  pd(21,35) = pd(21,35) + rrt(130) * density(44) * 2.0d0
  pd(21,44) = pd(21,44) + rrt(130) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(130) * density(44) 
  pd(35,44) = pd(35,44) - rrt(130) * density(35) 
  pd(44,35) = pd(44,35) - rrt(130) * density(44) 
  pd(44,44) = pd(44,44) - rrt(130) * density(35) 
  pd(14,41) = pd(14,41) + rrt(131) * density(44) 
  pd(14,44) = pd(14,44) + rrt(131) * density(41) 
  pd(29,41) = pd(29,41) + rrt(131) * density(44) 
  pd(29,44) = pd(29,44) + rrt(131) * density(41) 
  pd(41,41) = pd(41,41) - rrt(131) * density(44) 
  pd(41,44) = pd(41,44) - rrt(131) * density(41) 
  pd(44,41) = pd(44,41) - rrt(131) * density(44) 
  pd(44,44) = pd(44,44) - rrt(131) * density(41) 
  pd(15,41) = pd(15,41) + rrt(132) * density(44) 
  pd(15,44) = pd(15,44) + rrt(132) * density(41) 
  pd(29,41) = pd(29,41) + rrt(132) * density(44) 
  pd(29,44) = pd(29,44) + rrt(132) * density(41) 
  pd(41,41) = pd(41,41) - rrt(132) * density(44) 
  pd(41,44) = pd(41,44) - rrt(132) * density(41) 
  pd(44,41) = pd(44,41) - rrt(132) * density(44) 
  pd(44,44) = pd(44,44) - rrt(132) * density(41) 
  pd(01,43) = pd(01,43) + rrt(133) * density(44) 
  pd(01,44) = pd(01,44) + rrt(133) * density(43) 
  pd(21,43) = pd(21,43) + rrt(133) * density(44) 
  pd(21,44) = pd(21,44) + rrt(133) * density(43) 
  pd(43,43) = pd(43,43) - rrt(133) * density(44) 
  pd(43,44) = pd(43,44) - rrt(133) * density(43) 
  pd(44,43) = pd(44,43) - rrt(133) * density(44) 
  pd(44,44) = pd(44,44) - rrt(133) * density(43) 
  pd(14,17) = pd(14,17) + rrt(134) * density(44)**2 
  pd(14,44) = pd(14,44) + rrt(134) * density(17) * density(44) * 2.0d0
  pd(17,17) = pd(17,17) - rrt(134) * density(44)**2 
  pd(17,44) = pd(17,44) - rrt(134) * density(17) * density(44) * 2.0d0
  pd(44,17) = pd(44,17) - rrt(134) * density(44)**2 
  pd(44,44) = pd(44,44) - rrt(134) * density(17) * density(44) * 2.0d0
  pd(29,33) = pd(29,33) + rrt(135) * density(44)**2 
  pd(29,44) = pd(29,44) + rrt(135) * density(33) * density(44) * 2.0d0
  pd(33,33) = pd(33,33) - rrt(135) * density(44)**2 
  pd(33,44) = pd(33,44) - rrt(135) * density(33) * density(44) * 2.0d0
  pd(44,33) = pd(44,33) - rrt(135) * density(44)**2 
  pd(44,44) = pd(44,44) - rrt(135) * density(33) * density(44) * 2.0d0
  pd(14,17) = pd(14,17) + rrt(136) * density(44) 
  pd(14,44) = pd(14,44) + rrt(136) * density(17) 
  pd(17,17) = pd(17,17) - rrt(136) * density(44) 
  pd(17,44) = pd(17,44) - rrt(136) * density(17) 
  pd(44,17) = pd(44,17) - rrt(136) * density(44) 
  pd(44,44) = pd(44,44) - rrt(136) * density(17) 
  pd(29,33) = pd(29,33) + rrt(137) * density(44) 
  pd(29,44) = pd(29,44) + rrt(137) * density(33) 
  pd(33,33) = pd(33,33) - rrt(137) * density(44) 
  pd(33,44) = pd(33,44) - rrt(137) * density(33) 
  pd(44,33) = pd(44,33) - rrt(137) * density(44) 
  pd(44,44) = pd(44,44) - rrt(137) * density(33) 
  pd(29,32) = pd(29,32) + rrt(138) * density(44) 
  pd(29,44) = pd(29,44) + rrt(138) * density(32) 
  pd(32,32) = pd(32,32) - rrt(138) * density(44) 
  pd(32,44) = pd(32,44) - rrt(138) * density(32) 
  pd(37,32) = pd(37,32) + rrt(138) * density(44) 
  pd(37,44) = pd(37,44) + rrt(138) * density(32) 
  pd(44,32) = pd(44,32) - rrt(138) * density(44) 
  pd(44,44) = pd(44,44) - rrt(138) * density(32) 
  pd(21,32) = pd(21,32) + rrt(139) * density(44) 
  pd(21,44) = pd(21,44) + rrt(139) * density(32) 
  pd(32,32) = pd(32,32) - rrt(139) * density(44) 
  pd(32,44) = pd(32,44) - rrt(139) * density(32) 
  pd(36,32) = pd(36,32) + rrt(139) * density(44) 
  pd(36,44) = pd(36,44) + rrt(139) * density(32) 
  pd(44,32) = pd(44,32) - rrt(139) * density(44) 
  pd(44,44) = pd(44,44) - rrt(139) * density(32) 
  pd(29,21) = pd(29,21) - rrt(140) * density(29) * density(44) 
  pd(29,29) = pd(29,29) - rrt(140) * density(21) * density(44) 
  pd(29,44) = pd(29,44) - rrt(140) * density(21) * density(29) 
  pd(36,21) = pd(36,21) + rrt(140) * density(29) * density(44) 
  pd(36,29) = pd(36,29) + rrt(140) * density(21) * density(44) 
  pd(36,44) = pd(36,44) + rrt(140) * density(21) * density(29) 
  pd(44,21) = pd(44,21) - rrt(140) * density(29) * density(44) 
  pd(44,29) = pd(44,29) - rrt(140) * density(21) * density(44) 
  pd(44,44) = pd(44,44) - rrt(140) * density(21) * density(29) 
  pd(21,21) = pd(21,21) - rrt(141) * density(29) * density(44) 
  pd(21,29) = pd(21,29) - rrt(141) * density(21) * density(44) 
  pd(21,44) = pd(21,44) - rrt(141) * density(21) * density(29) 
  pd(37,21) = pd(37,21) + rrt(141) * density(29) * density(44) 
  pd(37,29) = pd(37,29) + rrt(141) * density(21) * density(44) 
  pd(37,44) = pd(37,44) + rrt(141) * density(21) * density(29) 
  pd(44,21) = pd(44,21) - rrt(141) * density(29) * density(44) 
  pd(44,29) = pd(44,29) - rrt(141) * density(21) * density(44) 
  pd(44,44) = pd(44,44) - rrt(141) * density(21) * density(29) 
  pd(32,21) = pd(32,21) - rrt(142) * density(32) * density(44) 
  pd(32,32) = pd(32,32) - rrt(142) * density(21) * density(44) 
  pd(32,44) = pd(32,44) - rrt(142) * density(21) * density(32) 
  pd(38,21) = pd(38,21) + rrt(142) * density(32) * density(44) 
  pd(38,32) = pd(38,32) + rrt(142) * density(21) * density(44) 
  pd(38,44) = pd(38,44) + rrt(142) * density(21) * density(32) 
  pd(44,21) = pd(44,21) - rrt(142) * density(32) * density(44) 
  pd(44,32) = pd(44,32) - rrt(142) * density(21) * density(44) 
  pd(44,44) = pd(44,44) - rrt(142) * density(21) * density(32) 
  pd(21,01) = pd(21,01) - rrt(143) * density(21) * density(44) 
  pd(21,21) = pd(21,21) - rrt(143) * density(01) * density(44) 
  pd(21,44) = pd(21,44) - rrt(143) * density(01) * density(21) 
  pd(37,01) = pd(37,01) + rrt(143) * density(21) * density(44) 
  pd(37,21) = pd(37,21) + rrt(143) * density(01) * density(44) 
  pd(37,44) = pd(37,44) + rrt(143) * density(01) * density(21) 
  pd(44,01) = pd(44,01) - rrt(143) * density(21) * density(44) 
  pd(44,21) = pd(44,21) - rrt(143) * density(01) * density(44) 
  pd(44,44) = pd(44,44) - rrt(143) * density(01) * density(21) 
  pd(40,40) = pd(40,40) - rrt(144) * density(44) 
  pd(40,44) = pd(40,44) - rrt(144) * density(40) 
  pd(42,40) = pd(42,40) + rrt(144) * density(44) 
  pd(42,44) = pd(42,44) + rrt(144) * density(40) 
  pd(44,40) = pd(44,40) - rrt(144) * density(44) 
  pd(44,44) = pd(44,44) - rrt(144) * density(40) 
  pd(21,29) = pd(21,29) + rrt(145) * density(36) 
  pd(21,36) = pd(21,36) + rrt(145) * density(29) 
  pd(29,29) = pd(29,29) - rrt(145) * density(36) 
  pd(29,36) = pd(29,36) - rrt(145) * density(29) 
  pd(36,29) = pd(36,29) - rrt(145) * density(36) 
  pd(36,36) = pd(36,36) - rrt(145) * density(29) 
  pd(44,29) = pd(44,29) + rrt(145) * density(36) 
  pd(44,36) = pd(44,36) + rrt(145) * density(29) 
  pd(14,14) = pd(14,14) - rrt(146) * density(36) 
  pd(14,36) = pd(14,36) - rrt(146) * density(14) 
  pd(36,14) = pd(36,14) - rrt(146) * density(36) 
  pd(36,36) = pd(36,36) - rrt(146) * density(14) 
  pd(40,14) = pd(40,14) + rrt(146) * density(36) 
  pd(40,36) = pd(40,36) + rrt(146) * density(14) 
  pd(44,14) = pd(44,14) + rrt(146) * density(36) 
  pd(44,36) = pd(44,36) + rrt(146) * density(14) 
  pd(21,21) = pd(21,21) - rrt(147) * density(36) 
  pd(21,36) = pd(21,36) - rrt(147) * density(21) 
  pd(32,21) = pd(32,21) + rrt(147) * density(36) 
  pd(32,36) = pd(32,36) + rrt(147) * density(21) 
  pd(36,21) = pd(36,21) - rrt(147) * density(36) 
  pd(36,36) = pd(36,36) - rrt(147) * density(21) 
  pd(44,21) = pd(44,21) + rrt(147) * density(36) 
  pd(44,36) = pd(44,36) + rrt(147) * density(21) 
  pd(26,26) = pd(26,26) - rrt(148) * density(36) 
  pd(26,36) = pd(26,36) - rrt(148) * density(26) 
  pd(32,26) = pd(32,26) + rrt(148) * density(36) 
  pd(32,36) = pd(32,36) + rrt(148) * density(26) 
  pd(36,26) = pd(36,26) - rrt(148) * density(36) 
  pd(36,36) = pd(36,36) - rrt(148) * density(26) 
  pd(44,26) = pd(44,26) + rrt(148) * density(36) 
  pd(44,36) = pd(44,36) + rrt(148) * density(26) 
  pd(21,27) = pd(21,27) + rrt(149) * density(36) 
  pd(21,36) = pd(21,36) + rrt(149) * density(27) 
  pd(27,27) = pd(27,27) - rrt(149) * density(36) 
  pd(27,36) = pd(27,36) - rrt(149) * density(27) 
  pd(29,27) = pd(29,27) + rrt(149) * density(36) 
  pd(29,36) = pd(29,36) + rrt(149) * density(27) 
  pd(36,27) = pd(36,27) - rrt(149) * density(36) 
  pd(36,36) = pd(36,36) - rrt(149) * density(27) 
  pd(44,27) = pd(44,27) + rrt(149) * density(36) 
  pd(44,36) = pd(44,36) + rrt(149) * density(27) 
  pd(01,10) = pd(01,10) + rrt(150) * density(36) 
  pd(01,36) = pd(01,36) + rrt(150) * density(10) 
  pd(10,10) = pd(10,10) - rrt(150) * density(36) 
  pd(10,36) = pd(10,36) - rrt(150) * density(10) 
  pd(29,10) = pd(29,10) + rrt(150) * density(36) 
  pd(29,36) = pd(29,36) + rrt(150) * density(10) 
  pd(36,10) = pd(36,10) - rrt(150) * density(36) 
  pd(36,36) = pd(36,36) - rrt(150) * density(10) 
  pd(44,10) = pd(44,10) + rrt(150) * density(36) 
  pd(44,36) = pd(44,36) + rrt(150) * density(10) 
  pd(01,11) = pd(01,11) + rrt(151) * density(36) 
  pd(01,36) = pd(01,36) + rrt(151) * density(11) 
  pd(11,11) = pd(11,11) - rrt(151) * density(36) 
  pd(11,36) = pd(11,36) - rrt(151) * density(11) 
  pd(29,11) = pd(29,11) + rrt(151) * density(36) 
  pd(29,36) = pd(29,36) + rrt(151) * density(11) 
  pd(36,11) = pd(36,11) - rrt(151) * density(36) 
  pd(36,36) = pd(36,36) - rrt(151) * density(11) 
  pd(44,11) = pd(44,11) + rrt(151) * density(36) 
  pd(44,36) = pd(44,36) + rrt(151) * density(11) 
  pd(21,32) = pd(21,32) + rrt(152) * density(36) * 2.0d0
  pd(21,36) = pd(21,36) + rrt(152) * density(32) * 2.0d0
  pd(32,32) = pd(32,32) - rrt(152) * density(36) 
  pd(32,36) = pd(32,36) - rrt(152) * density(32) 
  pd(36,32) = pd(36,32) - rrt(152) * density(36) 
  pd(36,36) = pd(36,36) - rrt(152) * density(32) 
  pd(44,32) = pd(44,32) + rrt(152) * density(36) 
  pd(44,36) = pd(44,36) + rrt(152) * density(32) 
  pd(29,29) = pd(29,29) - rrt(153) * density(37) 
  pd(29,37) = pd(29,37) - rrt(153) * density(29) 
  pd(32,29) = pd(32,29) + rrt(153) * density(37) 
  pd(32,37) = pd(32,37) + rrt(153) * density(29) 
  pd(37,29) = pd(37,29) - rrt(153) * density(37) 
  pd(37,37) = pd(37,37) - rrt(153) * density(29) 
  pd(44,29) = pd(44,29) + rrt(153) * density(37) 
  pd(44,37) = pd(44,37) + rrt(153) * density(29) 
  pd(21,21) = pd(21,21) + rrt(154) * density(37) 
  pd(21,37) = pd(21,37) + rrt(154) * density(21) 
  pd(37,21) = pd(37,21) - rrt(154) * density(37) 
  pd(37,37) = pd(37,37) - rrt(154) * density(21) 
  pd(44,21) = pd(44,21) + rrt(154) * density(37) 
  pd(44,37) = pd(44,37) + rrt(154) * density(21) 
  pd(21,26) = pd(21,26) + rrt(155) * density(37) * 2.0d0
  pd(21,37) = pd(21,37) + rrt(155) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(155) * density(37) 
  pd(26,37) = pd(26,37) - rrt(155) * density(26) 
  pd(37,26) = pd(37,26) - rrt(155) * density(37) 
  pd(37,37) = pd(37,37) - rrt(155) * density(26) 
  pd(44,26) = pd(44,26) + rrt(155) * density(37) 
  pd(44,37) = pd(44,37) + rrt(155) * density(26) 
  pd(21,27) = pd(21,27) + rrt(156) * density(37) * 2.0d0
  pd(21,37) = pd(21,37) + rrt(156) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(156) * density(37) 
  pd(27,37) = pd(27,37) - rrt(156) * density(27) 
  pd(37,27) = pd(37,27) - rrt(156) * density(37) 
  pd(37,37) = pd(37,37) - rrt(156) * density(27) 
  pd(44,27) = pd(44,27) + rrt(156) * density(37) 
  pd(44,37) = pd(44,37) + rrt(156) * density(27) 
  pd(21,01) = pd(21,01) + rrt(157) * density(37) 
  pd(21,37) = pd(21,37) + rrt(157) * density(01) 
  pd(37,01) = pd(37,01) - rrt(157) * density(37) 
  pd(37,37) = pd(37,37) - rrt(157) * density(01) 
  pd(44,01) = pd(44,01) + rrt(157) * density(37) 
  pd(44,37) = pd(44,37) + rrt(157) * density(01) 
  pd(01,10) = pd(01,10) + rrt(158) * density(37) 
  pd(01,37) = pd(01,37) + rrt(158) * density(10) 
  pd(10,10) = pd(10,10) - rrt(158) * density(37) 
  pd(10,37) = pd(10,37) - rrt(158) * density(10) 
  pd(21,10) = pd(21,10) + rrt(158) * density(37) 
  pd(21,37) = pd(21,37) + rrt(158) * density(10) 
  pd(37,10) = pd(37,10) - rrt(158) * density(37) 
  pd(37,37) = pd(37,37) - rrt(158) * density(10) 
  pd(44,10) = pd(44,10) + rrt(158) * density(37) 
  pd(44,37) = pd(44,37) + rrt(158) * density(10) 
  pd(01,11) = pd(01,11) + rrt(159) * density(37) 
  pd(01,37) = pd(01,37) + rrt(159) * density(11) 
  pd(11,11) = pd(11,11) - rrt(159) * density(37) 
  pd(11,37) = pd(11,37) - rrt(159) * density(11) 
  pd(21,11) = pd(21,11) + rrt(159) * density(37) 
  pd(21,37) = pd(21,37) + rrt(159) * density(11) 
  pd(37,11) = pd(37,11) - rrt(159) * density(37) 
  pd(37,37) = pd(37,37) - rrt(159) * density(11) 
  pd(44,11) = pd(44,11) + rrt(159) * density(37) 
  pd(44,37) = pd(44,37) + rrt(159) * density(11) 
  pd(21,29) = pd(21,29) + rrt(160) * density(38) * 2.0d0
  pd(21,38) = pd(21,38) + rrt(160) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(160) * density(38) 
  pd(29,38) = pd(29,38) - rrt(160) * density(29) 
  pd(38,29) = pd(38,29) - rrt(160) * density(38) 
  pd(38,38) = pd(38,38) - rrt(160) * density(29) 
  pd(44,29) = pd(44,29) + rrt(160) * density(38) 
  pd(44,38) = pd(44,38) + rrt(160) * density(29) 
  pd(01,10) = pd(01,10) + rrt(161) 
  pd(10,10) = pd(10,10) - rrt(161) 
  pd(10,11) = pd(10,11) + rrt(162) 
  pd(11,11) = pd(11,11) - rrt(162) 
  pd(01,12) = pd(01,12) + rrt(163) 
  pd(12,12) = pd(12,12) - rrt(163) 
  pd(11,13) = pd(11,13) + rrt(164) 
  pd(13,13) = pd(13,13) - rrt(164) 
  pd(21,26) = pd(21,26) + rrt(165) 
  pd(26,26) = pd(26,26) - rrt(165) 
  pd(26,27) = pd(26,27) + rrt(166) 
  pd(27,27) = pd(27,27) - rrt(166) 
  pd(21,27) = pd(21,27) + rrt(167) 
  pd(27,27) = pd(27,27) - rrt(167) 
  pd(21,28) = pd(21,28) + rrt(168) 
  pd(28,28) = pd(28,28) - rrt(168) 
  pd(10,10) = pd(10,10) - rrt(169) * density(29) 
  pd(10,29) = pd(10,29) - rrt(169) * density(10) 
  pd(15,10) = pd(15,10) + rrt(169) * density(29) 
  pd(15,29) = pd(15,29) + rrt(169) * density(10) 
  pd(29,10) = pd(29,10) - rrt(169) * density(29) 
  pd(29,29) = pd(29,29) - rrt(169) * density(10) 
  pd(40,10) = pd(40,10) + rrt(169) * density(29) 
  pd(40,29) = pd(40,29) + rrt(169) * density(10) 
  pd(01,10) = pd(01,10) + rrt(170) * density(29) 
  pd(01,29) = pd(01,29) + rrt(170) * density(10) 
  pd(10,10) = pd(10,10) - rrt(170) * density(29) 
  pd(10,29) = pd(10,29) - rrt(170) * density(10) 
  pd(29,10) = pd(29,10) - rrt(170) * density(29) 
  pd(29,29) = pd(29,29) - rrt(170) * density(10) 
  pd(31,10) = pd(31,10) + rrt(170) * density(29) 
  pd(31,29) = pd(31,29) + rrt(170) * density(10) 
  pd(01,10) = pd(01,10) + rrt(171) * density(14) 
  pd(01,14) = pd(01,14) + rrt(171) * density(10) 
  pd(10,10) = pd(10,10) - rrt(171) * density(14) 
  pd(10,14) = pd(10,14) - rrt(171) * density(10) 
  pd(01,10) = pd(01,10) + rrt(172) * density(14) 
  pd(01,14) = pd(01,14) + rrt(172) * density(10) 
  pd(10,10) = pd(10,10) - rrt(172) * density(14) 
  pd(10,14) = pd(10,14) - rrt(172) * density(10) 
  pd(14,10) = pd(14,10) - rrt(172) * density(14) 
  pd(14,14) = pd(14,14) - rrt(172) * density(10) 
  pd(16,10) = pd(16,10) + rrt(172) * density(14) 
  pd(16,14) = pd(16,14) + rrt(172) * density(10) 
  pd(01,10) = pd(01,10) + rrt(173) * density(21) 
  pd(01,21) = pd(01,21) + rrt(173) * density(10) 
  pd(10,10) = pd(10,10) - rrt(173) * density(21) 
  pd(10,21) = pd(10,21) - rrt(173) * density(10) 
  pd(21,10) = pd(21,10) - rrt(173) * density(21) 
  pd(21,21) = pd(21,21) - rrt(173) * density(10) 
  pd(29,10) = pd(29,10) + rrt(173) * density(21) 
  pd(29,21) = pd(29,21) + rrt(173) * density(10) 
  pd(30,10) = pd(30,10) + rrt(173) * density(21) 
  pd(30,21) = pd(30,21) + rrt(173) * density(10) 
  pd(01,10) = pd(01,10) + rrt(174) * density(21) 
  pd(01,21) = pd(01,21) + rrt(174) * density(10) 
  pd(10,10) = pd(10,10) - rrt(174) * density(21) 
  pd(10,21) = pd(10,21) - rrt(174) * density(10) 
  pd(21,10) = pd(21,10) - rrt(174) * density(21) 
  pd(21,21) = pd(21,21) - rrt(174) * density(10) 
  pd(26,10) = pd(26,10) + rrt(174) * density(21) 
  pd(26,21) = pd(26,21) + rrt(174) * density(10) 
  pd(01,10) = pd(01,10) + rrt(175) * density(21) 
  pd(01,21) = pd(01,21) + rrt(175) * density(10) 
  pd(10,10) = pd(10,10) - rrt(175) * density(21) 
  pd(10,21) = pd(10,21) - rrt(175) * density(10) 
  pd(21,10) = pd(21,10) - rrt(175) * density(21) 
  pd(21,21) = pd(21,21) - rrt(175) * density(10) 
  pd(27,10) = pd(27,10) + rrt(175) * density(21) 
  pd(27,21) = pd(27,21) + rrt(175) * density(10) 
  pd(01,01) = pd(01,01) + rrt(176) * density(10) 
  pd(01,10) = pd(01,10) + rrt(176) * density(01) 
  pd(10,01) = pd(10,01) - rrt(176) * density(10) 
  pd(10,10) = pd(10,10) - rrt(176) * density(01) 
  pd(01,10) = pd(01,10) + rrt(177) * density(40) 
  pd(01,40) = pd(01,40) + rrt(177) * density(10) 
  pd(10,10) = pd(10,10) - rrt(177) * density(40) 
  pd(10,40) = pd(10,40) - rrt(177) * density(10) 
  pd(01,10) = pd(01,10) + rrt(178) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(178) * density(10) * 4.0d0
  pd(11,10) = pd(11,10) + rrt(178) * density(10) * 2.0d0
  pd(01,10) = pd(01,10) + rrt(179) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(179) * density(10) * 4.0d0
  pd(13,10) = pd(13,10) + rrt(179) * density(10) * 2.0d0
  pd(10,01) = pd(10,01) + rrt(180) * density(11) 
  pd(10,11) = pd(10,11) + rrt(180) * density(01) 
  pd(11,01) = pd(11,01) - rrt(180) * density(11) 
  pd(11,11) = pd(11,11) - rrt(180) * density(01) 
  pd(01,01) = pd(01,01) + rrt(181) * density(11) 
  pd(01,11) = pd(01,11) + rrt(181) * density(01) 
  pd(11,01) = pd(11,01) - rrt(181) * density(11) 
  pd(11,11) = pd(11,11) - rrt(181) * density(01) 
  pd(01,11) = pd(01,11) + rrt(182) * density(21) 
  pd(01,21) = pd(01,21) + rrt(182) * density(11) 
  pd(11,11) = pd(11,11) - rrt(182) * density(21) 
  pd(11,21) = pd(11,21) - rrt(182) * density(11) 
  pd(21,11) = pd(21,11) - rrt(182) * density(21) 
  pd(21,21) = pd(21,21) - rrt(182) * density(11) 
  pd(29,11) = pd(29,11) + rrt(182) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(182) * density(11) * 2.0d0
  pd(10,11) = pd(10,11) + rrt(183) * density(40) 
  pd(10,40) = pd(10,40) + rrt(183) * density(11) 
  pd(11,11) = pd(11,11) - rrt(183) * density(40) 
  pd(11,40) = pd(11,40) - rrt(183) * density(11) 
  pd(12,01) = pd(12,01) + rrt(184) * density(13) 
  pd(12,13) = pd(12,13) + rrt(184) * density(01) 
  pd(13,01) = pd(13,01) - rrt(184) * density(13) 
  pd(13,13) = pd(13,13) - rrt(184) * density(01) 
  pd(01,13) = pd(01,13) + rrt(185) * density(21) 
  pd(01,21) = pd(01,21) + rrt(185) * density(13) 
  pd(13,13) = pd(13,13) - rrt(185) * density(21) 
  pd(13,21) = pd(13,21) - rrt(185) * density(13) 
  pd(21,13) = pd(21,13) - rrt(185) * density(21) 
  pd(21,21) = pd(21,21) - rrt(185) * density(13) 
  pd(29,13) = pd(29,13) + rrt(185) * density(21) 
  pd(29,21) = pd(29,21) + rrt(185) * density(13) 
  pd(31,13) = pd(31,13) + rrt(185) * density(21) 
  pd(31,21) = pd(31,21) + rrt(185) * density(13) 
  pd(11,01) = pd(11,01) + rrt(186) * density(12) 
  pd(11,12) = pd(11,12) + rrt(186) * density(01) 
  pd(12,01) = pd(12,01) - rrt(186) * density(12) 
  pd(12,12) = pd(12,12) - rrt(186) * density(01) 
  pd(01,12) = pd(01,12) + rrt(187) * density(21) 
  pd(01,21) = pd(01,21) + rrt(187) * density(12) 
  pd(12,12) = pd(12,12) - rrt(187) * density(21) 
  pd(12,21) = pd(12,21) - rrt(187) * density(12) 
  pd(21,12) = pd(21,12) - rrt(187) * density(21) 
  pd(21,21) = pd(21,21) - rrt(187) * density(12) 
  pd(29,12) = pd(29,12) + rrt(187) * density(21) 
  pd(29,21) = pd(29,21) + rrt(187) * density(12) 
  pd(30,12) = pd(30,12) + rrt(187) * density(21) 
  pd(30,21) = pd(30,21) + rrt(187) * density(12) 
  pd(01,12) = pd(01,12) + rrt(188) * density(40) 
  pd(01,40) = pd(01,40) + rrt(188) * density(12) 
  pd(12,12) = pd(12,12) - rrt(188) * density(40) 
  pd(12,40) = pd(12,40) - rrt(188) * density(12) 
  pd(14,12) = pd(14,12) + rrt(188) * density(40) 
  pd(14,40) = pd(14,40) + rrt(188) * density(12) 
  pd(29,12) = pd(29,12) + rrt(188) * density(40) 
  pd(29,40) = pd(29,40) + rrt(188) * density(12) 
  pd(40,12) = pd(40,12) - rrt(188) * density(40) 
  pd(40,40) = pd(40,40) - rrt(188) * density(12) 
  pd(10,10) = pd(10,10) - rrt(189) * density(12) 
  pd(10,12) = pd(10,12) - rrt(189) * density(10) 
  pd(12,10) = pd(12,10) - rrt(189) * density(12) 
  pd(12,12) = pd(12,12) - rrt(189) * density(10) 
  pd(20,10) = pd(20,10) + rrt(189) * density(12) 
  pd(20,12) = pd(20,12) + rrt(189) * density(10) 
  pd(44,10) = pd(44,10) + rrt(189) * density(12) 
  pd(44,12) = pd(44,12) + rrt(189) * density(10) 
  pd(12,12) = pd(12,12) - rrt(190) * density(12) * 4.0d0
  pd(20,12) = pd(20,12) + rrt(190) * density(12) * 2.0d0
  pd(44,12) = pd(44,12) + rrt(190) * density(12) * 2.0d0
  pd(10,14) = pd(10,14) + rrt(191) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(191) * density(14) * 4.0d0
  pd(10,14) = pd(10,14) + rrt(192) * density(14)**2 * 3.0d0
  pd(14,14) = pd(14,14) - rrt(192) * density(14)**2 * 6.0d0
  pd(10,14) = pd(10,14) + rrt(193) * density(14) * density(29) * 2.0d0
  pd(10,29) = pd(10,29) + rrt(193) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(193) * density(14) * density(29) * 4.0d0
  pd(14,29) = pd(14,29) - rrt(193) * density(14)**2 * 2.0d0
  pd(11,14) = pd(11,14) + rrt(194) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(194) * density(14) * 4.0d0
  pd(11,14) = pd(11,14) + rrt(195) * density(14)**2 * 3.0d0
  pd(14,14) = pd(14,14) - rrt(195) * density(14)**2 * 6.0d0
  pd(11,14) = pd(11,14) + rrt(196) * density(14) * density(29) * 2.0d0
  pd(11,29) = pd(11,29) + rrt(196) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(196) * density(14) * density(29) * 4.0d0
  pd(14,29) = pd(14,29) - rrt(196) * density(14)**2 * 2.0d0
  pd(14,15) = pd(14,15) + rrt(197) * density(29) 
  pd(14,29) = pd(14,29) + rrt(197) * density(15) 
  pd(15,15) = pd(15,15) - rrt(197) * density(29) 
  pd(15,29) = pd(15,29) - rrt(197) * density(15) 
  pd(29,15) = pd(29,15) - rrt(197) * density(29) 
  pd(29,29) = pd(29,29) - rrt(197) * density(15) 
  pd(30,15) = pd(30,15) + rrt(197) * density(29) 
  pd(30,29) = pd(30,29) + rrt(197) * density(15) 
  pd(15,15) = pd(15,15) - rrt(198) * density(21) 
  pd(15,21) = pd(15,21) - rrt(198) * density(15) 
  pd(21,15) = pd(21,15) - rrt(198) * density(21) 
  pd(21,21) = pd(21,21) - rrt(198) * density(15) 
  pd(29,15) = pd(29,15) + rrt(198) * density(21) 
  pd(29,21) = pd(29,21) + rrt(198) * density(15) 
  pd(40,15) = pd(40,15) + rrt(198) * density(21) 
  pd(40,21) = pd(40,21) + rrt(198) * density(15) 
  pd(01,15) = pd(01,15) + rrt(199) * density(40) 
  pd(01,40) = pd(01,40) + rrt(199) * density(15) 
  pd(15,15) = pd(15,15) - rrt(199) * density(40) 
  pd(15,40) = pd(15,40) - rrt(199) * density(15) 
  pd(29,15) = pd(29,15) + rrt(199) * density(40) 
  pd(29,40) = pd(29,40) + rrt(199) * density(15) 
  pd(40,15) = pd(40,15) - rrt(199) * density(40) 
  pd(40,40) = pd(40,40) - rrt(199) * density(15) 
  pd(14,01) = pd(14,01) + rrt(200) * density(15) 
  pd(14,15) = pd(14,15) + rrt(200) * density(01) 
  pd(15,01) = pd(15,01) - rrt(200) * density(15) 
  pd(15,15) = pd(15,15) - rrt(200) * density(01) 
  pd(14,14) = pd(14,14) + rrt(201) * density(16) 
  pd(14,16) = pd(14,16) + rrt(201) * density(14) 
  pd(16,14) = pd(16,14) - rrt(201) * density(16) 
  pd(16,16) = pd(16,16) - rrt(201) * density(14) 
  pd(14,16) = pd(14,16) + rrt(202) * density(29) 
  pd(14,29) = pd(14,29) + rrt(202) * density(16) 
  pd(16,16) = pd(16,16) - rrt(202) * density(29) 
  pd(16,29) = pd(16,29) - rrt(202) * density(16) 
  pd(15,14) = pd(15,14) + rrt(203) * density(16) 
  pd(15,16) = pd(15,16) + rrt(203) * density(14) 
  pd(16,14) = pd(16,14) - rrt(203) * density(16) 
  pd(16,16) = pd(16,16) - rrt(203) * density(14) 
  pd(14,01) = pd(14,01) + rrt(204) * density(16) 
  pd(14,16) = pd(14,16) + rrt(204) * density(01) 
  pd(16,01) = pd(16,01) - rrt(204) * density(16) 
  pd(16,16) = pd(16,16) - rrt(204) * density(01) 
  pd(15,15) = pd(15,15) - rrt(205) * density(16) 
  pd(15,16) = pd(15,16) - rrt(205) * density(15) 
  pd(16,15) = pd(16,15) - rrt(205) * density(16) 
  pd(16,16) = pd(16,16) - rrt(205) * density(15) 
  pd(18,15) = pd(18,15) + rrt(205) * density(16) 
  pd(18,16) = pd(18,16) + rrt(205) * density(15) 
  pd(44,15) = pd(44,15) + rrt(205) * density(16) 
  pd(44,16) = pd(44,16) + rrt(205) * density(15) 
  pd(16,16) = pd(16,16) - rrt(206) * density(21) 
  pd(16,21) = pd(16,21) - rrt(206) * density(16) 
  pd(21,16) = pd(21,16) - rrt(206) * density(21) 
  pd(21,21) = pd(21,21) - rrt(206) * density(16) 
  pd(29,16) = pd(29,16) + rrt(206) * density(21) 
  pd(29,21) = pd(29,21) + rrt(206) * density(16) 
  pd(40,16) = pd(40,16) + rrt(206) * density(21) 
  pd(40,21) = pd(40,21) + rrt(206) * density(16) 
  pd(10,16) = pd(10,16) + rrt(207) * density(40) 
  pd(10,40) = pd(10,40) + rrt(207) * density(16) 
  pd(16,16) = pd(16,16) - rrt(207) * density(40) 
  pd(16,40) = pd(16,40) - rrt(207) * density(16) 
  pd(29,16) = pd(29,16) + rrt(207) * density(40) 
  pd(29,40) = pd(29,40) + rrt(207) * density(16) 
  pd(40,16) = pd(40,16) - rrt(207) * density(40) 
  pd(40,40) = pd(40,40) - rrt(207) * density(16) 
  pd(21,26) = pd(21,26) + rrt(208) * density(29) 
  pd(21,29) = pd(21,29) + rrt(208) * density(26) 
  pd(26,26) = pd(26,26) - rrt(208) * density(29) 
  pd(26,29) = pd(26,29) - rrt(208) * density(26) 
  pd(14,14) = pd(14,14) - rrt(209) * density(26) 
  pd(14,26) = pd(14,26) - rrt(209) * density(14) 
  pd(26,14) = pd(26,14) - rrt(209) * density(26) 
  pd(26,26) = pd(26,26) - rrt(209) * density(14) 
  pd(29,14) = pd(29,14) + rrt(209) * density(26) 
  pd(29,26) = pd(29,26) + rrt(209) * density(14) 
  pd(40,14) = pd(40,14) + rrt(209) * density(26) 
  pd(40,26) = pd(40,26) + rrt(209) * density(14) 
  pd(21,21) = pd(21,21) + rrt(210) * density(26) 
  pd(21,26) = pd(21,26) + rrt(210) * density(21) 
  pd(26,21) = pd(26,21) - rrt(210) * density(26) 
  pd(26,26) = pd(26,26) - rrt(210) * density(21) 
  pd(21,01) = pd(21,01) + rrt(211) * density(26) 
  pd(21,26) = pd(21,26) + rrt(211) * density(01) 
  pd(26,01) = pd(26,01) - rrt(211) * density(26) 
  pd(26,26) = pd(26,26) - rrt(211) * density(01) 
  pd(21,26) = pd(21,26) + rrt(212) * density(40) 
  pd(21,40) = pd(21,40) + rrt(212) * density(26) 
  pd(26,26) = pd(26,26) - rrt(212) * density(40) 
  pd(26,40) = pd(26,40) - rrt(212) * density(26) 
  pd(21,26) = pd(21,26) + rrt(213) * density(32) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(213) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(213) * density(32) 
  pd(26,32) = pd(26,32) - rrt(213) * density(26) 
  pd(30,26) = pd(30,26) + rrt(213) * density(32) 
  pd(30,32) = pd(30,32) + rrt(213) * density(26) 
  pd(32,26) = pd(32,26) - rrt(213) * density(32) 
  pd(32,32) = pd(32,32) - rrt(213) * density(26) 
  pd(21,26) = pd(21,26) + rrt(214) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(214) * density(26) * 4.0d0
  pd(27,26) = pd(27,26) + rrt(214) * density(26) * 2.0d0
  pd(21,29) = pd(21,29) + rrt(215) * density(32) 
  pd(21,32) = pd(21,32) + rrt(215) * density(29) 
  pd(26,29) = pd(26,29) + rrt(215) * density(32) 
  pd(26,32) = pd(26,32) + rrt(215) * density(29) 
  pd(29,29) = pd(29,29) - rrt(215) * density(32) 
  pd(29,32) = pd(29,32) - rrt(215) * density(29) 
  pd(32,29) = pd(32,29) - rrt(215) * density(32) 
  pd(32,32) = pd(32,32) - rrt(215) * density(29) 
  pd(26,27) = pd(26,27) + rrt(216) * density(29) 
  pd(26,29) = pd(26,29) + rrt(216) * density(27) 
  pd(27,27) = pd(27,27) - rrt(216) * density(29) 
  pd(27,29) = pd(27,29) - rrt(216) * density(27) 
  pd(21,27) = pd(21,27) + rrt(217) * density(29) 
  pd(21,29) = pd(21,29) + rrt(217) * density(27) 
  pd(27,27) = pd(27,27) - rrt(217) * density(29) 
  pd(27,29) = pd(27,29) - rrt(217) * density(27) 
  pd(29,27) = pd(29,27) - rrt(217) * density(29) 
  pd(29,29) = pd(29,29) - rrt(217) * density(27) 
  pd(30,27) = pd(30,27) + rrt(217) * density(29) 
  pd(30,29) = pd(30,29) + rrt(217) * density(27) 
  pd(26,21) = pd(26,21) + rrt(218) * density(27) 
  pd(26,27) = pd(26,27) + rrt(218) * density(21) 
  pd(27,21) = pd(27,21) - rrt(218) * density(27) 
  pd(27,27) = pd(27,27) - rrt(218) * density(21) 
  pd(26,01) = pd(26,01) + rrt(219) * density(27) 
  pd(26,27) = pd(26,27) + rrt(219) * density(01) 
  pd(27,01) = pd(27,01) - rrt(219) * density(27) 
  pd(27,27) = pd(27,27) - rrt(219) * density(01) 
  pd(26,27) = pd(26,27) + rrt(220) * density(40) 
  pd(26,40) = pd(26,40) + rrt(220) * density(27) 
  pd(27,27) = pd(27,27) - rrt(220) * density(40) 
  pd(27,40) = pd(27,40) - rrt(220) * density(27) 
  pd(21,27) = pd(21,27) + rrt(221) * density(32) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(221) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(221) * density(32) 
  pd(27,32) = pd(27,32) - rrt(221) * density(27) 
  pd(29,27) = pd(29,27) + rrt(221) * density(32) 
  pd(29,32) = pd(29,32) + rrt(221) * density(27) 
  pd(32,27) = pd(32,27) - rrt(221) * density(32) 
  pd(32,32) = pd(32,32) - rrt(221) * density(27) 
  pd(21,28) = pd(21,28) + rrt(222) * density(29) 
  pd(21,29) = pd(21,29) + rrt(222) * density(28) 
  pd(28,28) = pd(28,28) - rrt(222) * density(29) 
  pd(28,29) = pd(28,29) - rrt(222) * density(28) 
  pd(29,28) = pd(29,28) - rrt(222) * density(29) 
  pd(29,29) = pd(29,29) - rrt(222) * density(28) 
  pd(31,28) = pd(31,28) + rrt(222) * density(29) 
  pd(31,29) = pd(31,29) + rrt(222) * density(28) 
  pd(21,21) = pd(21,21) - rrt(223) * density(28) 
  pd(21,28) = pd(21,28) - rrt(223) * density(21) 
  pd(27,21) = pd(27,21) + rrt(223) * density(28) * 2.0d0
  pd(27,28) = pd(27,28) + rrt(223) * density(21) * 2.0d0
  pd(28,21) = pd(28,21) - rrt(223) * density(28) 
  pd(28,28) = pd(28,28) - rrt(223) * density(21) 
  pd(27,01) = pd(27,01) + rrt(224) * density(28) 
  pd(27,28) = pd(27,28) + rrt(224) * density(01) 
  pd(28,01) = pd(28,01) - rrt(224) * density(28) 
  pd(28,28) = pd(28,28) - rrt(224) * density(01) 
  pd(26,29) = pd(26,29) + rrt(225) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(225) * density(29) * 4.0d0
  pd(26,21) = pd(26,21) + rrt(226) * density(29)**2 
  pd(26,29) = pd(26,29) + rrt(226) * density(21) * density(29) * 2.0d0
  pd(29,21) = pd(29,21) - rrt(226) * density(29)**2 * 2.0d0
  pd(29,29) = pd(29,29) - rrt(226) * density(21) * density(29) * 4.0d0
  pd(26,29) = pd(26,29) + rrt(227) * density(29)**2 * 3.0d0
  pd(29,29) = pd(29,29) - rrt(227) * density(29)**2 * 6.0d0
  pd(27,29) = pd(27,29) + rrt(228) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(228) * density(29) * 4.0d0
  pd(21,21) = pd(21,21) - rrt(229) * density(26)**2 
  pd(21,26) = pd(21,26) - rrt(229) * density(21) * density(26) * 2.0d0
  pd(26,21) = pd(26,21) - rrt(229) * density(26)**2 * 2.0d0
  pd(26,26) = pd(26,26) - rrt(229) * density(21) * density(26) * 4.0d0
  pd(32,21) = pd(32,21) + rrt(229) * density(26)**2 * 2.0d0
  pd(32,26) = pd(32,26) + rrt(229) * density(21) * density(26) * 4.0d0
  pd(29,29) = pd(29,29) + rrt(230) * density(30) 
  pd(29,30) = pd(29,30) + rrt(230) * density(29) 
  pd(30,29) = pd(30,29) - rrt(230) * density(30) 
  pd(30,30) = pd(30,30) - rrt(230) * density(29) 
  pd(29,21) = pd(29,21) + rrt(231) * density(30) 
  pd(29,30) = pd(29,30) + rrt(231) * density(21) 
  pd(30,21) = pd(30,21) - rrt(231) * density(30) 
  pd(30,30) = pd(30,30) - rrt(231) * density(21) 
  pd(21,21) = pd(21,21) - rrt(232) * density(30) 
  pd(21,30) = pd(21,30) - rrt(232) * density(21) 
  pd(26,21) = pd(26,21) + rrt(232) * density(30) 
  pd(26,30) = pd(26,30) + rrt(232) * density(21) 
  pd(29,21) = pd(29,21) + rrt(232) * density(30) 
  pd(29,30) = pd(29,30) + rrt(232) * density(21) 
  pd(30,21) = pd(30,21) - rrt(232) * density(30) 
  pd(30,30) = pd(30,30) - rrt(232) * density(21) 
  pd(21,21) = pd(21,21) - rrt(233) * density(30) 
  pd(21,30) = pd(21,30) - rrt(233) * density(21) 
  pd(27,21) = pd(27,21) + rrt(233) * density(30) 
  pd(27,30) = pd(27,30) + rrt(233) * density(21) 
  pd(29,21) = pd(29,21) + rrt(233) * density(30) 
  pd(29,30) = pd(29,30) + rrt(233) * density(21) 
  pd(30,21) = pd(30,21) - rrt(233) * density(30) 
  pd(30,30) = pd(30,30) - rrt(233) * density(21) 
  pd(29,01) = pd(29,01) + rrt(234) * density(30) 
  pd(29,30) = pd(29,30) + rrt(234) * density(01) 
  pd(30,01) = pd(30,01) - rrt(234) * density(30) 
  pd(30,30) = pd(30,30) - rrt(234) * density(01) 
  pd(21,30) = pd(21,30) + rrt(235) * density(32) 
  pd(21,32) = pd(21,32) + rrt(235) * density(30) 
  pd(29,30) = pd(29,30) + rrt(235) * density(32) * 2.0d0
  pd(29,32) = pd(29,32) + rrt(235) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(235) * density(32) 
  pd(30,32) = pd(30,32) - rrt(235) * density(30) 
  pd(32,30) = pd(32,30) - rrt(235) * density(32) 
  pd(32,32) = pd(32,32) - rrt(235) * density(30) 
  pd(21,30) = pd(21,30) + rrt(236) * density(32) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(236) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(236) * density(32) 
  pd(30,32) = pd(30,32) - rrt(236) * density(30) 
  pd(32,30) = pd(32,30) - rrt(236) * density(32) 
  pd(32,32) = pd(32,32) - rrt(236) * density(30) 
  pd(14,30) = pd(14,30) + rrt(237) * density(40) 
  pd(14,40) = pd(14,40) + rrt(237) * density(30) 
  pd(21,30) = pd(21,30) + rrt(237) * density(40) 
  pd(21,40) = pd(21,40) + rrt(237) * density(30) 
  pd(30,30) = pd(30,30) - rrt(237) * density(40) 
  pd(30,40) = pd(30,40) - rrt(237) * density(30) 
  pd(40,30) = pd(40,30) - rrt(237) * density(40) 
  pd(40,40) = pd(40,40) - rrt(237) * density(30) 
  pd(30,29) = pd(30,29) + rrt(238) * density(31) 
  pd(30,31) = pd(30,31) + rrt(238) * density(29) 
  pd(31,29) = pd(31,29) - rrt(238) * density(31) 
  pd(31,31) = pd(31,31) - rrt(238) * density(29) 
  pd(29,14) = pd(29,14) + rrt(239) * density(31) 
  pd(29,31) = pd(29,31) + rrt(239) * density(14) 
  pd(31,14) = pd(31,14) - rrt(239) * density(31) 
  pd(31,31) = pd(31,31) - rrt(239) * density(14) 
  pd(30,21) = pd(30,21) + rrt(240) * density(31) 
  pd(30,31) = pd(30,31) + rrt(240) * density(21) 
  pd(31,21) = pd(31,21) - rrt(240) * density(31) 
  pd(31,31) = pd(31,31) - rrt(240) * density(21) 
  pd(21,21) = pd(21,21) - rrt(241) * density(31) 
  pd(21,31) = pd(21,31) - rrt(241) * density(21) 
  pd(29,21) = pd(29,21) + rrt(241) * density(31) * 3.0d0
  pd(29,31) = pd(29,31) + rrt(241) * density(21) * 3.0d0
  pd(31,21) = pd(31,21) - rrt(241) * density(31) 
  pd(31,31) = pd(31,31) - rrt(241) * density(21) 
  pd(29,01) = pd(29,01) + rrt(242) * density(31) 
  pd(29,31) = pd(29,31) + rrt(242) * density(01) 
  pd(31,01) = pd(31,01) - rrt(242) * density(31) 
  pd(31,31) = pd(31,31) - rrt(242) * density(01) 
  pd(26,26) = pd(26,26) - rrt(243) * density(31) 
  pd(26,31) = pd(26,31) - rrt(243) * density(26) 
  pd(29,26) = pd(29,26) + rrt(243) * density(31) * 3.0d0
  pd(29,31) = pd(29,31) + rrt(243) * density(26) * 3.0d0
  pd(31,26) = pd(31,26) - rrt(243) * density(31) 
  pd(31,31) = pd(31,31) - rrt(243) * density(26) 
  pd(26,26) = pd(26,26) - rrt(244) * density(31) 
  pd(26,31) = pd(26,31) - rrt(244) * density(26) 
  pd(27,26) = pd(27,26) + rrt(244) * density(31) 
  pd(27,31) = pd(27,31) + rrt(244) * density(26) 
  pd(30,26) = pd(30,26) + rrt(244) * density(31) 
  pd(30,31) = pd(30,31) + rrt(244) * density(26) 
  pd(31,26) = pd(31,26) - rrt(244) * density(31) 
  pd(31,31) = pd(31,31) - rrt(244) * density(26) 
  pd(26,26) = pd(26,26) - rrt(245) * density(31) 
  pd(26,31) = pd(26,31) - rrt(245) * density(26) 
  pd(29,26) = pd(29,26) + rrt(245) * density(31) * 3.0d0
  pd(29,31) = pd(29,31) + rrt(245) * density(26) * 3.0d0
  pd(31,26) = pd(31,26) - rrt(245) * density(31) 
  pd(31,31) = pd(31,31) - rrt(245) * density(26) 
  pd(29,31) = pd(29,31) + rrt(246) * density(40) 
  pd(29,40) = pd(29,40) + rrt(246) * density(31) 
  pd(31,31) = pd(31,31) - rrt(246) * density(40) 
  pd(31,40) = pd(31,40) - rrt(246) * density(31) 
  pd(30,31) = pd(30,31) + rrt(247) * density(40) 
  pd(30,40) = pd(30,40) + rrt(247) * density(31) 
  pd(31,31) = pd(31,31) - rrt(247) * density(40) 
  pd(31,40) = pd(31,40) - rrt(247) * density(31) 
  pd(21,31) = pd(21,31) + rrt(248) * density(32) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(248) * density(31) * 2.0d0
  pd(31,31) = pd(31,31) - rrt(248) * density(32) 
  pd(31,32) = pd(31,32) - rrt(248) * density(31) 
  pd(32,31) = pd(32,31) - rrt(248) * density(32) 
  pd(32,32) = pd(32,32) - rrt(248) * density(31) 
  pd(21,31) = pd(21,31) + rrt(249) * density(32) 
  pd(21,32) = pd(21,32) + rrt(249) * density(31) 
  pd(29,31) = pd(29,31) + rrt(249) * density(32) 
  pd(29,32) = pd(29,32) + rrt(249) * density(31) 
  pd(30,31) = pd(30,31) + rrt(249) * density(32) 
  pd(30,32) = pd(30,32) + rrt(249) * density(31) 
  pd(31,31) = pd(31,31) - rrt(249) * density(32) 
  pd(31,32) = pd(31,32) - rrt(249) * density(31) 
  pd(32,31) = pd(32,31) - rrt(249) * density(32) 
  pd(32,32) = pd(32,32) - rrt(249) * density(31) 
  pd(01,14) = pd(01,14) + rrt(250) * density(40) 
  pd(01,40) = pd(01,40) + rrt(250) * density(14) 
  pd(14,14) = pd(14,14) - rrt(250) * density(40) 
  pd(14,40) = pd(14,40) - rrt(250) * density(14) 
  pd(29,14) = pd(29,14) + rrt(250) * density(40) 
  pd(29,40) = pd(29,40) + rrt(250) * density(14) 
  pd(40,14) = pd(40,14) - rrt(250) * density(40) 
  pd(40,40) = pd(40,40) - rrt(250) * density(14) 
  pd(14,14) = pd(14,14) - rrt(251) * density(21) 
  pd(14,21) = pd(14,21) - rrt(251) * density(14) 
  pd(21,14) = pd(21,14) - rrt(251) * density(21) 
  pd(21,21) = pd(21,21) - rrt(251) * density(14) 
  pd(29,14) = pd(29,14) + rrt(251) * density(21) 
  pd(29,21) = pd(29,21) + rrt(251) * density(14) 
  pd(40,14) = pd(40,14) + rrt(251) * density(21) 
  pd(40,21) = pd(40,21) + rrt(251) * density(14) 
  pd(14,14) = pd(14,14) - rrt(252) * density(32) 
  pd(14,32) = pd(14,32) - rrt(252) * density(14) 
  pd(21,14) = pd(21,14) + rrt(252) * density(32) 
  pd(21,32) = pd(21,32) + rrt(252) * density(14) 
  pd(32,14) = pd(32,14) - rrt(252) * density(32) 
  pd(32,32) = pd(32,32) - rrt(252) * density(14) 
  pd(40,14) = pd(40,14) + rrt(252) * density(32) 
  pd(40,32) = pd(40,32) + rrt(252) * density(14) 
  pd(01,01) = pd(01,01) - rrt(253) * density(29) 
  pd(01,29) = pd(01,29) - rrt(253) * density(01) 
  pd(14,01) = pd(14,01) + rrt(253) * density(29) 
  pd(14,29) = pd(14,29) + rrt(253) * density(01) 
  pd(29,01) = pd(29,01) - rrt(253) * density(29) 
  pd(29,29) = pd(29,29) - rrt(253) * density(01) 
  pd(40,01) = pd(40,01) + rrt(253) * density(29) 
  pd(40,29) = pd(40,29) + rrt(253) * density(01) 
  pd(14,29) = pd(14,29) + rrt(254) * density(40) 
  pd(14,40) = pd(14,40) + rrt(254) * density(29) 
  pd(21,29) = pd(21,29) + rrt(254) * density(40) 
  pd(21,40) = pd(21,40) + rrt(254) * density(29) 
  pd(29,29) = pd(29,29) - rrt(254) * density(40) 
  pd(29,40) = pd(29,40) - rrt(254) * density(29) 
  pd(40,29) = pd(40,29) - rrt(254) * density(40) 
  pd(40,40) = pd(40,40) - rrt(254) * density(29) 
  pd(21,29) = pd(21,29) + rrt(255) * density(32) 
  pd(21,32) = pd(21,32) + rrt(255) * density(29) 
  pd(26,29) = pd(26,29) + rrt(255) * density(32) 
  pd(26,32) = pd(26,32) + rrt(255) * density(29) 
  pd(29,29) = pd(29,29) - rrt(255) * density(32) 
  pd(29,32) = pd(29,32) - rrt(255) * density(29) 
  pd(32,29) = pd(32,29) - rrt(255) * density(32) 
  pd(32,32) = pd(32,32) - rrt(255) * density(29) 
  pd(01,40) = pd(01,40) + rrt(256) * density(40) * 2.0d0
  pd(21,40) = pd(21,40) + rrt(256) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(256) * density(40) * 4.0d0
  pd(21,21) = pd(21,21) - rrt(257) * density(21) * 4.0d0
  pd(29,21) = pd(29,21) + rrt(257) * density(21) * 2.0d0
  pd(32,21) = pd(32,21) + rrt(257) * density(21) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(258) * density(14) * 4.0d0
  pd(18,14) = pd(18,14) + rrt(258) * density(14) * 2.0d0
  pd(44,14) = pd(44,14) + rrt(258) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(259) * density(29) 
  pd(14,29) = pd(14,29) - rrt(259) * density(14) 
  pd(29,14) = pd(29,14) - rrt(259) * density(29) 
  pd(29,29) = pd(29,29) - rrt(259) * density(14) 
  pd(41,14) = pd(41,14) + rrt(259) * density(29) 
  pd(41,29) = pd(41,29) + rrt(259) * density(14) 
  pd(44,14) = pd(44,14) + rrt(259) * density(29) 
  pd(44,29) = pd(44,29) + rrt(259) * density(14) 
  pd(01,01) = pd(01,01) - rrt(260) 
  pd(14,01) = pd(14,01) + rrt(260) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(261) * density(14) 
  pd(01,14) = pd(01,14) - rrt(261) * density(01) 
  pd(14,01) = pd(14,01) + rrt(261) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) + rrt(261) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(262) * density(29) 
  pd(01,29) = pd(01,29) - rrt(262) * density(01) 
  pd(14,01) = pd(14,01) + rrt(262) * density(29) * 2.0d0
  pd(14,29) = pd(14,29) + rrt(262) * density(01) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(263) 
  pd(29,21) = pd(29,21) + rrt(263) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(264) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(264) * density(21) * 4.0d0
  pd(21,21) = pd(21,21) - rrt(265) * density(29) 
  pd(21,29) = pd(21,29) - rrt(265) * density(21) 
  pd(29,21) = pd(29,21) + rrt(265) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) + rrt(265) * density(21) * 2.0d0
  pd(14,40) = pd(14,40) + rrt(266) 
  pd(29,40) = pd(29,40) + rrt(266) 
  pd(40,40) = pd(40,40) - rrt(266) 
  pd(14,14) = pd(14,14) + rrt(267) * density(40) 
  pd(14,40) = pd(14,40) + rrt(267) * density(14) 
  pd(29,14) = pd(29,14) + rrt(267) * density(40) 
  pd(29,40) = pd(29,40) + rrt(267) * density(14) 
  pd(40,14) = pd(40,14) - rrt(267) * density(40) 
  pd(40,40) = pd(40,40) - rrt(267) * density(14) 
  pd(14,29) = pd(14,29) + rrt(268) * density(40) 
  pd(14,40) = pd(14,40) + rrt(268) * density(29) 
  pd(29,29) = pd(29,29) + rrt(268) * density(40) 
  pd(29,40) = pd(29,40) + rrt(268) * density(29) 
  pd(40,29) = pd(40,29) - rrt(268) * density(40) 
  pd(40,40) = pd(40,40) - rrt(268) * density(29) 
  pd(21,01) = pd(21,01) + rrt(269) * density(32) 
  pd(21,32) = pd(21,32) + rrt(269) * density(01) 
  pd(29,01) = pd(29,01) + rrt(269) * density(32) 
  pd(29,32) = pd(29,32) + rrt(269) * density(01) 
  pd(32,01) = pd(32,01) - rrt(269) * density(32) 
  pd(32,32) = pd(32,32) - rrt(269) * density(01) 
  pd(21,21) = pd(21,21) + rrt(270) * density(32) 
  pd(21,32) = pd(21,32) + rrt(270) * density(21) 
  pd(29,21) = pd(29,21) + rrt(270) * density(32) 
  pd(29,32) = pd(29,32) + rrt(270) * density(21) 
  pd(32,21) = pd(32,21) - rrt(270) * density(32) 
  pd(32,32) = pd(32,32) - rrt(270) * density(21) 
  pd(21,14) = pd(21,14) + rrt(271) * density(32) 
  pd(21,32) = pd(21,32) + rrt(271) * density(14) 
  pd(29,14) = pd(29,14) + rrt(271) * density(32) 
  pd(29,32) = pd(29,32) + rrt(271) * density(14) 
  pd(32,14) = pd(32,14) - rrt(271) * density(32) 
  pd(32,32) = pd(32,32) - rrt(271) * density(14) 
  pd(21,29) = pd(21,29) + rrt(272) * density(32) 
  pd(21,32) = pd(21,32) + rrt(272) * density(29) 
  pd(29,29) = pd(29,29) + rrt(272) * density(32) 
  pd(29,32) = pd(29,32) + rrt(272) * density(29) 
  pd(32,29) = pd(32,29) - rrt(272) * density(32) 
  pd(32,32) = pd(32,32) - rrt(272) * density(29) 
  pd(01,01) = pd(01,01) + rrt(273) * density(14)**2 
  pd(01,14) = pd(01,14) + rrt(273) * density(01) * density(14) * 2.0d0
  pd(14,01) = pd(14,01) - rrt(273) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(273) * density(01) * density(14) * 4.0d0
  pd(01,14) = pd(01,14) + rrt(274) * density(14) * density(21) * 2.0d0
  pd(01,21) = pd(01,21) + rrt(274) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(274) * density(14) * density(21) * 4.0d0
  pd(14,21) = pd(14,21) - rrt(274) * density(14)**2 * 2.0d0
  pd(01,14) = pd(01,14) + rrt(275) * density(14) * density(40) * 2.0d0
  pd(01,40) = pd(01,40) + rrt(275) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(275) * density(14) * density(40) * 4.0d0
  pd(14,40) = pd(14,40) - rrt(275) * density(14)**2 * 2.0d0
  pd(01,14) = pd(01,14) + rrt(276) * density(14)**2 * 3.0d0
  pd(14,14) = pd(14,14) - rrt(276) * density(14)**2 * 6.0d0
  pd(01,14) = pd(01,14) + rrt(277) * density(14) * density(29) * 2.0d0
  pd(01,29) = pd(01,29) + rrt(277) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(277) * density(14) * density(29) * 4.0d0
  pd(14,29) = pd(14,29) - rrt(277) * density(14)**2 * 2.0d0
  pd(21,01) = pd(21,01) + rrt(278) * density(29)**2 
  pd(21,29) = pd(21,29) + rrt(278) * density(01) * density(29) * 2.0d0
  pd(29,01) = pd(29,01) - rrt(278) * density(29)**2 * 2.0d0
  pd(29,29) = pd(29,29) - rrt(278) * density(01) * density(29) * 4.0d0
  pd(21,21) = pd(21,21) + rrt(279) * density(29)**2 
  pd(21,29) = pd(21,29) + rrt(279) * density(21) * density(29) * 2.0d0
  pd(29,21) = pd(29,21) - rrt(279) * density(29)**2 * 2.0d0
  pd(29,29) = pd(29,29) - rrt(279) * density(21) * density(29) * 4.0d0
  pd(21,29) = pd(21,29) + rrt(280) * density(29) * density(40) * 2.0d0
  pd(21,40) = pd(21,40) + rrt(280) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(280) * density(29) * density(40) * 4.0d0
  pd(29,40) = pd(29,40) - rrt(280) * density(29)**2 * 2.0d0
  pd(21,14) = pd(21,14) + rrt(281) * density(29)**2 
  pd(21,29) = pd(21,29) + rrt(281) * density(14) * density(29) * 2.0d0
  pd(29,14) = pd(29,14) - rrt(281) * density(29)**2 * 2.0d0
  pd(29,29) = pd(29,29) - rrt(281) * density(14) * density(29) * 4.0d0
  pd(21,29) = pd(21,29) + rrt(282) * density(29)**2 * 3.0d0
  pd(29,29) = pd(29,29) - rrt(282) * density(29)**2 * 6.0d0
  pd(14,01) = pd(14,01) - rrt(283) * density(14) * density(29) 
  pd(14,14) = pd(14,14) - rrt(283) * density(01) * density(29) 
  pd(14,29) = pd(14,29) - rrt(283) * density(01) * density(14) 
  pd(29,01) = pd(29,01) - rrt(283) * density(14) * density(29) 
  pd(29,14) = pd(29,14) - rrt(283) * density(01) * density(29) 
  pd(29,29) = pd(29,29) - rrt(283) * density(01) * density(14) 
  pd(40,01) = pd(40,01) + rrt(283) * density(14) * density(29) 
  pd(40,14) = pd(40,14) + rrt(283) * density(01) * density(29) 
  pd(40,29) = pd(40,29) + rrt(283) * density(01) * density(14) 
  pd(14,14) = pd(14,14) - rrt(284) * density(21) * density(29) 
  pd(14,21) = pd(14,21) - rrt(284) * density(14) * density(29) 
  pd(14,29) = pd(14,29) - rrt(284) * density(14) * density(21) 
  pd(29,14) = pd(29,14) - rrt(284) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(284) * density(14) * density(29) 
  pd(29,29) = pd(29,29) - rrt(284) * density(14) * density(21) 
  pd(40,14) = pd(40,14) + rrt(284) * density(21) * density(29) 
  pd(40,21) = pd(40,21) + rrt(284) * density(14) * density(29) 
  pd(40,29) = pd(40,29) + rrt(284) * density(14) * density(21) 
  pd(14,14) = pd(14,14) - rrt(285) * density(29) * density(40) 
  pd(14,29) = pd(14,29) - rrt(285) * density(14) * density(40) 
  pd(14,40) = pd(14,40) - rrt(285) * density(14) * density(29) 
  pd(29,14) = pd(29,14) - rrt(285) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(285) * density(14) * density(40) 
  pd(29,40) = pd(29,40) - rrt(285) * density(14) * density(29) 
  pd(40,14) = pd(40,14) + rrt(285) * density(29) * density(40) 
  pd(40,29) = pd(40,29) + rrt(285) * density(14) * density(40) 
  pd(40,40) = pd(40,40) + rrt(285) * density(14) * density(29) 
  pd(14,14) = pd(14,14) - rrt(286) * density(14) * density(29) * 2.0d0
  pd(14,29) = pd(14,29) - rrt(286) * density(14)**2 
  pd(29,14) = pd(29,14) - rrt(286) * density(14) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(286) * density(14)**2 
  pd(40,14) = pd(40,14) + rrt(286) * density(14) * density(29) * 2.0d0
  pd(40,29) = pd(40,29) + rrt(286) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(287) * density(29)**2 
  pd(14,29) = pd(14,29) - rrt(287) * density(14) * density(29) * 2.0d0
  pd(29,14) = pd(29,14) - rrt(287) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(287) * density(14) * density(29) * 2.0d0
  pd(40,14) = pd(40,14) + rrt(287) * density(29)**2 
  pd(40,29) = pd(40,29) + rrt(287) * density(14) * density(29) * 2.0d0
  pd(21,01) = pd(21,01) - rrt(288) * density(21) * density(29) 
  pd(21,21) = pd(21,21) - rrt(288) * density(01) * density(29) 
  pd(21,29) = pd(21,29) - rrt(288) * density(01) * density(21) 
  pd(29,01) = pd(29,01) - rrt(288) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(288) * density(01) * density(29) 
  pd(29,29) = pd(29,29) - rrt(288) * density(01) * density(21) 
  pd(32,01) = pd(32,01) + rrt(288) * density(21) * density(29) 
  pd(32,21) = pd(32,21) + rrt(288) * density(01) * density(29) 
  pd(32,29) = pd(32,29) + rrt(288) * density(01) * density(21) 
  pd(21,21) = pd(21,21) - rrt(289) * density(21) * density(29) * 2.0d0
  pd(21,29) = pd(21,29) - rrt(289) * density(21)**2 
  pd(29,21) = pd(29,21) - rrt(289) * density(21) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(289) * density(21)**2 
  pd(32,21) = pd(32,21) + rrt(289) * density(21) * density(29) * 2.0d0
  pd(32,29) = pd(32,29) + rrt(289) * density(21)**2 
  pd(21,21) = pd(21,21) - rrt(290) * density(29) * density(40) 
  pd(21,29) = pd(21,29) - rrt(290) * density(21) * density(40) 
  pd(21,40) = pd(21,40) - rrt(290) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(290) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(290) * density(21) * density(40) 
  pd(29,40) = pd(29,40) - rrt(290) * density(21) * density(29) 
  pd(32,21) = pd(32,21) + rrt(290) * density(29) * density(40) 
  pd(32,29) = pd(32,29) + rrt(290) * density(21) * density(40) 
  pd(32,40) = pd(32,40) + rrt(290) * density(21) * density(29) 
  pd(21,14) = pd(21,14) - rrt(291) * density(21) * density(29) 
  pd(21,21) = pd(21,21) - rrt(291) * density(14) * density(29) 
  pd(21,29) = pd(21,29) - rrt(291) * density(14) * density(21) 
  pd(29,14) = pd(29,14) - rrt(291) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(291) * density(14) * density(29) 
  pd(29,29) = pd(29,29) - rrt(291) * density(14) * density(21) 
  pd(32,14) = pd(32,14) + rrt(291) * density(21) * density(29) 
  pd(32,21) = pd(32,21) + rrt(291) * density(14) * density(29) 
  pd(32,29) = pd(32,29) + rrt(291) * density(14) * density(21) 
  pd(21,21) = pd(21,21) - rrt(292) * density(29)**2 
  pd(21,29) = pd(21,29) - rrt(292) * density(21) * density(29) * 2.0d0
  pd(29,21) = pd(29,21) - rrt(292) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(292) * density(21) * density(29) * 2.0d0
  pd(32,21) = pd(32,21) + rrt(292) * density(29)**2 
  pd(32,29) = pd(32,29) + rrt(292) * density(21) * density(29) * 2.0d0
  pd(14,17) = pd(14,17) + rrt(293) * density(29) 
  pd(14,29) = pd(14,29) + rrt(293) * density(17) 
  pd(17,17) = pd(17,17) - rrt(293) * density(29) 
  pd(17,29) = pd(17,29) - rrt(293) * density(17) 
  pd(29,17) = pd(29,17) - rrt(293) * density(29) 
  pd(29,29) = pd(29,29) - rrt(293) * density(17) 
  pd(33,17) = pd(33,17) + rrt(293) * density(29) 
  pd(33,29) = pd(33,29) + rrt(293) * density(17) 
  pd(14,17) = pd(14,17) + rrt(294) * density(21) 
  pd(14,21) = pd(14,21) + rrt(294) * density(17) 
  pd(17,17) = pd(17,17) - rrt(294) * density(21) 
  pd(17,21) = pd(17,21) - rrt(294) * density(17) 
  pd(21,17) = pd(21,17) - rrt(294) * density(21) 
  pd(21,21) = pd(21,21) - rrt(294) * density(17) 
  pd(34,17) = pd(34,17) + rrt(294) * density(21) 
  pd(34,21) = pd(34,21) + rrt(294) * density(17) 
  pd(17,17) = pd(17,17) - rrt(295) * density(21) 
  pd(17,21) = pd(17,21) - rrt(295) * density(17) 
  pd(21,17) = pd(21,17) - rrt(295) * density(21) 
  pd(21,21) = pd(21,21) - rrt(295) * density(17) 
  pd(29,17) = pd(29,17) + rrt(295) * density(21) 
  pd(29,21) = pd(29,21) + rrt(295) * density(17) 
  pd(41,17) = pd(41,17) + rrt(295) * density(21) 
  pd(41,21) = pd(41,21) + rrt(295) * density(17) 
  pd(17,17) = pd(17,17) - rrt(296) * density(21) 
  pd(17,21) = pd(17,21) - rrt(296) * density(17) 
  pd(21,17) = pd(21,17) - rrt(296) * density(21) 
  pd(21,21) = pd(21,21) - rrt(296) * density(17) 
  pd(33,17) = pd(33,17) + rrt(296) * density(21) 
  pd(33,21) = pd(33,21) + rrt(296) * density(17) 
  pd(40,17) = pd(40,17) + rrt(296) * density(21) 
  pd(40,21) = pd(40,21) + rrt(296) * density(17) 
  pd(17,17) = pd(17,17) - rrt(297) * density(32) 
  pd(17,32) = pd(17,32) - rrt(297) * density(17) 
  pd(21,17) = pd(21,17) + rrt(297) * density(32) 
  pd(21,32) = pd(21,32) + rrt(297) * density(17) 
  pd(32,17) = pd(32,17) - rrt(297) * density(32) 
  pd(32,32) = pd(32,32) - rrt(297) * density(17) 
  pd(41,17) = pd(41,17) + rrt(297) * density(32) 
  pd(41,32) = pd(41,32) + rrt(297) * density(17) 
  pd(14,17) = pd(14,17) + rrt(298) * density(40) 
  pd(14,40) = pd(14,40) + rrt(298) * density(17) 
  pd(17,17) = pd(17,17) - rrt(298) * density(40) 
  pd(17,40) = pd(17,40) - rrt(298) * density(17) 
  pd(40,17) = pd(40,17) - rrt(298) * density(40) 
  pd(40,40) = pd(40,40) - rrt(298) * density(17) 
  pd(41,17) = pd(41,17) + rrt(298) * density(40) 
  pd(41,40) = pd(41,40) + rrt(298) * density(17) 
  pd(17,17) = pd(17,17) - rrt(299) * density(40) 
  pd(17,40) = pd(17,40) - rrt(299) * density(17) 
  pd(18,17) = pd(18,17) + rrt(299) * density(40) 
  pd(18,40) = pd(18,40) + rrt(299) * density(17) 
  pd(29,17) = pd(29,17) + rrt(299) * density(40) 
  pd(29,40) = pd(29,40) + rrt(299) * density(17) 
  pd(40,17) = pd(40,17) - rrt(299) * density(40) 
  pd(40,40) = pd(40,40) - rrt(299) * density(17) 
  pd(01,17) = pd(01,17) + rrt(300) * density(40) 
  pd(01,40) = pd(01,40) + rrt(300) * density(17) 
  pd(17,17) = pd(17,17) - rrt(300) * density(40) 
  pd(17,40) = pd(17,40) - rrt(300) * density(17) 
  pd(33,17) = pd(33,17) + rrt(300) * density(40) 
  pd(33,40) = pd(33,40) + rrt(300) * density(17) 
  pd(40,17) = pd(40,17) - rrt(300) * density(40) 
  pd(40,40) = pd(40,40) - rrt(300) * density(17) 
  pd(01,01) = pd(01,01) - rrt(301) * density(33) 
  pd(01,33) = pd(01,33) - rrt(301) * density(01) 
  pd(14,01) = pd(14,01) + rrt(301) * density(33) 
  pd(14,33) = pd(14,33) + rrt(301) * density(01) 
  pd(33,01) = pd(33,01) - rrt(301) * density(33) 
  pd(33,33) = pd(33,33) - rrt(301) * density(01) 
  pd(41,01) = pd(41,01) + rrt(301) * density(33) 
  pd(41,33) = pd(41,33) + rrt(301) * density(01) 
  pd(21,21) = pd(21,21) - rrt(302) * density(33) 
  pd(21,33) = pd(21,33) - rrt(302) * density(21) 
  pd(29,21) = pd(29,21) + rrt(302) * density(33) 
  pd(29,33) = pd(29,33) + rrt(302) * density(21) 
  pd(33,21) = pd(33,21) - rrt(302) * density(33) 
  pd(33,33) = pd(33,33) - rrt(302) * density(21) 
  pd(34,21) = pd(34,21) + rrt(302) * density(33) 
  pd(34,33) = pd(34,33) + rrt(302) * density(21) 
  pd(21,32) = pd(21,32) + rrt(303) * density(33) 
  pd(21,33) = pd(21,33) + rrt(303) * density(32) 
  pd(32,32) = pd(32,32) - rrt(303) * density(33) 
  pd(32,33) = pd(32,33) - rrt(303) * density(32) 
  pd(33,32) = pd(33,32) - rrt(303) * density(33) 
  pd(33,33) = pd(33,33) - rrt(303) * density(32) 
  pd(34,32) = pd(34,32) + rrt(303) * density(33) 
  pd(34,33) = pd(34,33) + rrt(303) * density(32) 
  pd(29,33) = pd(29,33) + rrt(304) * density(40) 
  pd(29,40) = pd(29,40) + rrt(304) * density(33) 
  pd(33,33) = pd(33,33) - rrt(304) * density(40) 
  pd(33,40) = pd(33,40) - rrt(304) * density(33) 
  pd(40,33) = pd(40,33) - rrt(304) * density(40) 
  pd(40,40) = pd(40,40) - rrt(304) * density(33) 
  pd(41,33) = pd(41,33) + rrt(304) * density(40) 
  pd(41,40) = pd(41,40) + rrt(304) * density(33) 
  pd(14,33) = pd(14,33) + rrt(305) * density(40) 
  pd(14,40) = pd(14,40) + rrt(305) * density(33) 
  pd(33,33) = pd(33,33) - rrt(305) * density(40) 
  pd(33,40) = pd(33,40) - rrt(305) * density(33) 
  pd(34,33) = pd(34,33) + rrt(305) * density(40) 
  pd(34,40) = pd(34,40) + rrt(305) * density(33) 
  pd(40,33) = pd(40,33) - rrt(305) * density(40) 
  pd(40,40) = pd(40,40) - rrt(305) * density(33) 
  pd(15,15) = pd(15,15) - rrt(306) * density(33) 
  pd(15,33) = pd(15,33) - rrt(306) * density(15) 
  pd(17,15) = pd(17,15) + rrt(306) * density(33) 
  pd(17,33) = pd(17,33) + rrt(306) * density(15) 
  pd(29,15) = pd(29,15) + rrt(306) * density(33) 
  pd(29,33) = pd(29,33) + rrt(306) * density(15) 
  pd(33,15) = pd(33,15) - rrt(306) * density(33) 
  pd(33,33) = pd(33,33) - rrt(306) * density(15) 
  pd(01,18) = pd(01,18) + rrt(307) * density(21) 
  pd(01,21) = pd(01,21) + rrt(307) * density(18) 
  pd(18,18) = pd(18,18) - rrt(307) * density(21) 
  pd(18,21) = pd(18,21) - rrt(307) * density(18) 
  pd(21,18) = pd(21,18) - rrt(307) * density(21) 
  pd(21,21) = pd(21,21) - rrt(307) * density(18) 
  pd(34,18) = pd(34,18) + rrt(307) * density(21) 
  pd(34,21) = pd(34,21) + rrt(307) * density(18) 
  pd(14,18) = pd(14,18) + rrt(308) * density(29) 
  pd(14,29) = pd(14,29) + rrt(308) * density(18) 
  pd(18,18) = pd(18,18) - rrt(308) * density(29) 
  pd(18,29) = pd(18,29) - rrt(308) * density(18) 
  pd(29,18) = pd(29,18) - rrt(308) * density(29) 
  pd(29,29) = pd(29,29) - rrt(308) * density(18) 
  pd(41,18) = pd(41,18) + rrt(308) * density(29) 
  pd(41,29) = pd(41,29) + rrt(308) * density(18) 
  pd(01,18) = pd(01,18) + rrt(309) * density(29) 
  pd(01,29) = pd(01,29) + rrt(309) * density(18) 
  pd(18,18) = pd(18,18) - rrt(309) * density(29) 
  pd(18,29) = pd(18,29) - rrt(309) * density(18) 
  pd(29,18) = pd(29,18) - rrt(309) * density(29) 
  pd(29,29) = pd(29,29) - rrt(309) * density(18) 
  pd(33,18) = pd(33,18) + rrt(309) * density(29) 
  pd(33,29) = pd(33,29) + rrt(309) * density(18) 
  pd(01,18) = pd(01,18) + rrt(310) * density(32) 
  pd(01,32) = pd(01,32) + rrt(310) * density(18) 
  pd(18,18) = pd(18,18) - rrt(310) * density(32) 
  pd(18,32) = pd(18,32) - rrt(310) * density(18) 
  pd(29,18) = pd(29,18) + rrt(310) * density(32) 
  pd(29,32) = pd(29,32) + rrt(310) * density(18) 
  pd(32,18) = pd(32,18) - rrt(310) * density(32) 
  pd(32,32) = pd(32,32) - rrt(310) * density(18) 
  pd(34,18) = pd(34,18) + rrt(310) * density(32) 
  pd(34,32) = pd(34,32) + rrt(310) * density(18) 
  pd(01,14) = pd(01,14) + rrt(311) * density(18) 
  pd(01,18) = pd(01,18) + rrt(311) * density(14) 
  pd(14,14) = pd(14,14) - rrt(311) * density(18) 
  pd(14,18) = pd(14,18) - rrt(311) * density(14) 
  pd(17,14) = pd(17,14) + rrt(311) * density(18) 
  pd(17,18) = pd(17,18) + rrt(311) * density(14) 
  pd(18,14) = pd(18,14) - rrt(311) * density(18) 
  pd(18,18) = pd(18,18) - rrt(311) * density(14) 
  pd(01,18) = pd(01,18) + rrt(312) * density(40) 
  pd(01,40) = pd(01,40) + rrt(312) * density(18) 
  pd(18,18) = pd(18,18) - rrt(312) * density(40) 
  pd(18,40) = pd(18,40) - rrt(312) * density(18) 
  pd(40,18) = pd(40,18) - rrt(312) * density(40) 
  pd(40,40) = pd(40,40) - rrt(312) * density(18) 
  pd(41,18) = pd(41,18) + rrt(312) * density(40) 
  pd(41,40) = pd(41,40) + rrt(312) * density(18) 
  pd(01,01) = pd(01,01) - rrt(313) * density(34) 
  pd(01,34) = pd(01,34) - rrt(313) * density(01) 
  pd(34,01) = pd(34,01) - rrt(313) * density(34) 
  pd(34,34) = pd(34,34) - rrt(313) * density(01) 
  pd(40,01) = pd(40,01) + rrt(313) * density(34) 
  pd(40,34) = pd(40,34) + rrt(313) * density(01) 
  pd(41,01) = pd(41,01) + rrt(313) * density(34) 
  pd(41,34) = pd(41,34) + rrt(313) * density(01) 
  pd(14,14) = pd(14,14) - rrt(314) * density(34) 
  pd(14,34) = pd(14,34) - rrt(314) * density(14) 
  pd(29,14) = pd(29,14) + rrt(314) * density(34) 
  pd(29,34) = pd(29,34) + rrt(314) * density(14) 
  pd(34,14) = pd(34,14) - rrt(314) * density(34) 
  pd(34,34) = pd(34,34) - rrt(314) * density(14) 
  pd(41,14) = pd(41,14) + rrt(314) * density(34) 
  pd(41,34) = pd(41,34) + rrt(314) * density(14) 
  pd(21,34) = pd(21,34) + rrt(315) * density(40) 
  pd(21,40) = pd(21,40) + rrt(315) * density(34) 
  pd(34,34) = pd(34,34) - rrt(315) * density(40) 
  pd(34,40) = pd(34,40) - rrt(315) * density(34) 
  pd(40,34) = pd(40,34) - rrt(315) * density(40) 
  pd(40,40) = pd(40,40) - rrt(315) * density(34) 
  pd(41,34) = pd(41,34) + rrt(315) * density(40) 
  pd(41,40) = pd(41,40) + rrt(315) * density(34) 
  pd(01,19) = pd(01,19) + rrt(316) * density(21) 
  pd(01,21) = pd(01,21) + rrt(316) * density(19) 
  pd(14,19) = pd(14,19) + rrt(316) * density(21) 
  pd(14,21) = pd(14,21) + rrt(316) * density(19) 
  pd(19,19) = pd(19,19) - rrt(316) * density(21) 
  pd(19,21) = pd(19,21) - rrt(316) * density(19) 
  pd(21,19) = pd(21,19) - rrt(316) * density(21) 
  pd(21,21) = pd(21,21) - rrt(316) * density(19) 
  pd(34,19) = pd(34,19) + rrt(316) * density(21) 
  pd(34,21) = pd(34,21) + rrt(316) * density(19) 
  pd(01,14) = pd(01,14) + rrt(317) * density(19) 
  pd(01,19) = pd(01,19) + rrt(317) * density(14) 
  pd(14,14) = pd(14,14) - rrt(317) * density(19) 
  pd(14,19) = pd(14,19) - rrt(317) * density(14) 
  pd(18,14) = pd(18,14) + rrt(317) * density(19) 
  pd(18,19) = pd(18,19) + rrt(317) * density(14) 
  pd(19,14) = pd(19,14) - rrt(317) * density(19) 
  pd(19,19) = pd(19,19) - rrt(317) * density(14) 
  pd(01,19) = pd(01,19) + rrt(318) * density(40) 
  pd(01,40) = pd(01,40) + rrt(318) * density(19) 
  pd(14,19) = pd(14,19) + rrt(318) * density(40) 
  pd(14,40) = pd(14,40) + rrt(318) * density(19) 
  pd(19,19) = pd(19,19) - rrt(318) * density(40) 
  pd(19,40) = pd(19,40) - rrt(318) * density(19) 
  pd(40,19) = pd(40,19) - rrt(318) * density(40) 
  pd(40,40) = pd(40,40) - rrt(318) * density(19) 
  pd(41,19) = pd(41,19) + rrt(318) * density(40) 
  pd(41,40) = pd(41,40) + rrt(318) * density(19) 
  pd(01,01) = pd(01,01) + rrt(319) * density(20) 
  pd(01,20) = pd(01,20) + rrt(319) * density(01) 
  pd(18,01) = pd(18,01) + rrt(319) * density(20) 
  pd(18,20) = pd(18,20) + rrt(319) * density(01) 
  pd(20,01) = pd(20,01) - rrt(319) * density(20) 
  pd(20,20) = pd(20,20) - rrt(319) * density(01) 
  pd(01,20) = pd(01,20) + rrt(320) * density(21) * 2.0d0
  pd(01,21) = pd(01,21) + rrt(320) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(320) * density(21) 
  pd(20,21) = pd(20,21) - rrt(320) * density(20) 
  pd(21,20) = pd(21,20) - rrt(320) * density(21) 
  pd(21,21) = pd(21,21) - rrt(320) * density(20) 
  pd(34,20) = pd(34,20) + rrt(320) * density(21) 
  pd(34,21) = pd(34,21) + rrt(320) * density(20) 
  pd(01,20) = pd(01,20) + rrt(321) * density(29) * 2.0d0
  pd(01,29) = pd(01,29) + rrt(321) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(321) * density(29) 
  pd(20,29) = pd(20,29) - rrt(321) * density(20) 
  pd(29,20) = pd(29,20) - rrt(321) * density(29) 
  pd(29,29) = pd(29,29) - rrt(321) * density(20) 
  pd(33,20) = pd(33,20) + rrt(321) * density(29) 
  pd(33,29) = pd(33,29) + rrt(321) * density(20) 
  pd(01,14) = pd(01,14) + rrt(322) * density(20) * 2.0d0
  pd(01,20) = pd(01,20) + rrt(322) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(322) * density(20) 
  pd(14,20) = pd(14,20) - rrt(322) * density(14) 
  pd(17,14) = pd(17,14) + rrt(322) * density(20) 
  pd(17,20) = pd(17,20) + rrt(322) * density(14) 
  pd(20,14) = pd(20,14) - rrt(322) * density(20) 
  pd(20,20) = pd(20,20) - rrt(322) * density(14) 
  pd(01,20) = pd(01,20) + rrt(323) * density(40) * 2.0d0
  pd(01,40) = pd(01,40) + rrt(323) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(323) * density(40) 
  pd(20,40) = pd(20,40) - rrt(323) * density(20) 
  pd(40,20) = pd(40,20) - rrt(323) * density(40) 
  pd(40,40) = pd(40,40) - rrt(323) * density(20) 
  pd(41,20) = pd(41,20) + rrt(323) * density(40) 
  pd(41,40) = pd(41,40) + rrt(323) * density(20) 
  pd(01,01) = pd(01,01) - rrt(324) * density(35) 
  pd(01,35) = pd(01,35) - rrt(324) * density(01) 
  pd(21,01) = pd(21,01) + rrt(324) * density(35) 
  pd(21,35) = pd(21,35) + rrt(324) * density(01) 
  pd(35,01) = pd(35,01) - rrt(324) * density(35) 
  pd(35,35) = pd(35,35) - rrt(324) * density(01) 
  pd(43,01) = pd(43,01) + rrt(324) * density(35) 
  pd(43,35) = pd(43,35) + rrt(324) * density(01) 
  pd(21,21) = pd(21,21) + rrt(325) * density(35) 
  pd(21,35) = pd(21,35) + rrt(325) * density(21) 
  pd(34,21) = pd(34,21) + rrt(325) * density(35) 
  pd(34,35) = pd(34,35) + rrt(325) * density(21) 
  pd(35,21) = pd(35,21) - rrt(325) * density(35) 
  pd(35,35) = pd(35,35) - rrt(325) * density(21) 
  pd(21,26) = pd(21,26) + rrt(326) * density(35) * 2.0d0
  pd(21,35) = pd(21,35) + rrt(326) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(326) * density(35) 
  pd(26,35) = pd(26,35) - rrt(326) * density(26) 
  pd(34,26) = pd(34,26) + rrt(326) * density(35) 
  pd(34,35) = pd(34,35) + rrt(326) * density(26) 
  pd(35,26) = pd(35,26) - rrt(326) * density(35) 
  pd(35,35) = pd(35,35) - rrt(326) * density(26) 
  pd(21,27) = pd(21,27) + rrt(327) * density(35) * 2.0d0
  pd(21,35) = pd(21,35) + rrt(327) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(327) * density(35) 
  pd(27,35) = pd(27,35) - rrt(327) * density(27) 
  pd(34,27) = pd(34,27) + rrt(327) * density(35) 
  pd(34,35) = pd(34,35) + rrt(327) * density(27) 
  pd(35,27) = pd(35,27) - rrt(327) * density(35) 
  pd(35,35) = pd(35,35) - rrt(327) * density(27) 
  pd(29,29) = pd(29,29) - rrt(328) * density(35) 
  pd(29,35) = pd(29,35) - rrt(328) * density(29) 
  pd(32,29) = pd(32,29) + rrt(328) * density(35) 
  pd(32,35) = pd(32,35) + rrt(328) * density(29) 
  pd(34,29) = pd(34,29) + rrt(328) * density(35) 
  pd(34,35) = pd(34,35) + rrt(328) * density(29) 
  pd(35,29) = pd(35,29) - rrt(328) * density(35) 
  pd(35,35) = pd(35,35) - rrt(328) * density(29) 
  pd(21,35) = pd(21,35) + rrt(329) * density(40) * 2.0d0
  pd(21,40) = pd(21,40) + rrt(329) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(329) * density(40) 
  pd(35,40) = pd(35,40) - rrt(329) * density(35) 
  pd(40,35) = pd(40,35) - rrt(329) * density(40) 
  pd(40,40) = pd(40,40) - rrt(329) * density(35) 
  pd(41,35) = pd(41,35) + rrt(329) * density(40) 
  pd(41,40) = pd(41,40) + rrt(329) * density(35) 
  pd(01,01) = pd(01,01) + rrt(330) * density(43) 
  pd(01,43) = pd(01,43) + rrt(330) * density(01) 
  pd(34,01) = pd(34,01) + rrt(330) * density(43) 
  pd(34,43) = pd(34,43) + rrt(330) * density(01) 
  pd(43,01) = pd(43,01) - rrt(330) * density(43) 
  pd(43,43) = pd(43,43) - rrt(330) * density(01) 
  pd(01,21) = pd(01,21) + rrt(331) * density(43) 
  pd(01,43) = pd(01,43) + rrt(331) * density(21) 
  pd(21,21) = pd(21,21) - rrt(331) * density(43) 
  pd(21,43) = pd(21,43) - rrt(331) * density(21) 
  pd(35,21) = pd(35,21) + rrt(331) * density(43) 
  pd(35,43) = pd(35,43) + rrt(331) * density(21) 
  pd(43,21) = pd(43,21) - rrt(331) * density(43) 
  pd(43,43) = pd(43,43) - rrt(331) * density(21) 
  pd(01,01) = pd(01,01) - rrt(332) * density(01) * density(17) * 2.0d0
  pd(01,17) = pd(01,17) - rrt(332) * density(01)**2 
  pd(17,01) = pd(17,01) - rrt(332) * density(01) * density(17) * 2.0d0
  pd(17,17) = pd(17,17) - rrt(332) * density(01)**2 
  pd(19,01) = pd(19,01) + rrt(332) * density(01) * density(17) * 2.0d0
  pd(19,17) = pd(19,17) + rrt(332) * density(01)**2 
  pd(17,17) = pd(17,17) - rrt(333) * density(29) 
  pd(17,29) = pd(17,29) - rrt(333) * density(17) 
  pd(29,17) = pd(29,17) - rrt(333) * density(29) 
  pd(29,29) = pd(29,29) - rrt(333) * density(17) 
  pd(41,17) = pd(41,17) + rrt(333) * density(29) 
  pd(41,29) = pd(41,29) + rrt(333) * density(17) 
  pd(14,14) = pd(14,14) - rrt(334) * density(17) 
  pd(14,17) = pd(14,17) - rrt(334) * density(14) 
  pd(17,14) = pd(17,14) - rrt(334) * density(17) 
  pd(17,17) = pd(17,17) - rrt(334) * density(14) 
  pd(18,14) = pd(18,14) + rrt(334) * density(17) 
  pd(18,17) = pd(18,17) + rrt(334) * density(14) 
  pd(01,01) = pd(01,01) - rrt(335) * density(33) 
  pd(01,33) = pd(01,33) - rrt(335) * density(01) 
  pd(14,01) = pd(14,01) + rrt(335) * density(33) 
  pd(14,33) = pd(14,33) + rrt(335) * density(01) 
  pd(33,01) = pd(33,01) - rrt(335) * density(33) 
  pd(33,33) = pd(33,33) - rrt(335) * density(01) 
  pd(41,01) = pd(41,01) + rrt(335) * density(33) 
  pd(41,33) = pd(41,33) + rrt(335) * density(01) 
  pd(29,29) = pd(29,29) - rrt(336) * density(33) 
  pd(29,33) = pd(29,33) - rrt(336) * density(29) 
  pd(33,29) = pd(33,29) - rrt(336) * density(33) 
  pd(33,33) = pd(33,33) - rrt(336) * density(29) 
  pd(34,29) = pd(34,29) + rrt(336) * density(33) 
  pd(34,33) = pd(34,33) + rrt(336) * density(29) 
  pd(14,14) = pd(14,14) - rrt(337) * density(33) 
  pd(14,33) = pd(14,33) - rrt(337) * density(14) 
  pd(33,14) = pd(33,14) - rrt(337) * density(33) 
  pd(33,33) = pd(33,33) - rrt(337) * density(14) 
  pd(41,14) = pd(41,14) + rrt(337) * density(33) 
  pd(41,33) = pd(41,33) + rrt(337) * density(14) 
  pd(01,01) = pd(01,01) - rrt(338) * density(01) * density(18) * 2.0d0
  pd(01,18) = pd(01,18) - rrt(338) * density(01)**2 
  pd(18,01) = pd(18,01) - rrt(338) * density(01) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(338) * density(01)**2 
  pd(20,01) = pd(20,01) + rrt(338) * density(01) * density(18) * 2.0d0
  pd(20,18) = pd(20,18) + rrt(338) * density(01)**2 
  pd(14,01) = pd(14,01) - rrt(339) * density(14) * density(18) 
  pd(14,14) = pd(14,14) - rrt(339) * density(01) * density(18) 
  pd(14,18) = pd(14,18) - rrt(339) * density(01) * density(14) 
  pd(18,01) = pd(18,01) - rrt(339) * density(14) * density(18) 
  pd(18,14) = pd(18,14) - rrt(339) * density(01) * density(18) 
  pd(18,18) = pd(18,18) - rrt(339) * density(01) * density(14) 
  pd(19,01) = pd(19,01) + rrt(339) * density(14) * density(18) 
  pd(19,14) = pd(19,14) + rrt(339) * density(01) * density(18) 
  pd(19,18) = pd(19,18) + rrt(339) * density(01) * density(14) 
  pd(21,21) = pd(21,21) - rrt(340) * density(21) * density(34) * 2.0d0
  pd(21,34) = pd(21,34) - rrt(340) * density(21)**2 
  pd(34,21) = pd(34,21) - rrt(340) * density(21) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(340) * density(21)**2 
  pd(35,21) = pd(35,21) + rrt(340) * density(21) * density(34) * 2.0d0
  pd(35,34) = pd(35,34) + rrt(340) * density(21)**2 
  pd(01,01) = pd(01,01) - rrt(341) * density(01) * density(34) * 2.0d0
  pd(01,34) = pd(01,34) - rrt(341) * density(01)**2 
  pd(34,01) = pd(34,01) - rrt(341) * density(01) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(341) * density(01)**2 
  pd(43,01) = pd(43,01) + rrt(341) * density(01) * density(34) * 2.0d0
  pd(43,34) = pd(43,34) + rrt(341) * density(01)**2 
  pd(26,26) = pd(26,26) - rrt(342) * density(36) 
  pd(26,36) = pd(26,36) - rrt(342) * density(26) 
  pd(29,26) = pd(29,26) + rrt(342) * density(36) 
  pd(29,36) = pd(29,36) + rrt(342) * density(26) 
  pd(36,26) = pd(36,26) - rrt(342) * density(36) 
  pd(36,36) = pd(36,36) - rrt(342) * density(26) 
  pd(37,26) = pd(37,26) + rrt(342) * density(36) 
  pd(37,36) = pd(37,36) + rrt(342) * density(26) 
  pd(29,32) = pd(29,32) + rrt(343) * density(36) 
  pd(29,36) = pd(29,36) + rrt(343) * density(32) 
  pd(32,32) = pd(32,32) - rrt(343) * density(36) 
  pd(32,36) = pd(32,36) - rrt(343) * density(32) 
  pd(36,32) = pd(36,32) - rrt(343) * density(36) 
  pd(36,36) = pd(36,36) - rrt(343) * density(32) 
  pd(38,32) = pd(38,32) + rrt(343) * density(36) 
  pd(38,36) = pd(38,36) + rrt(343) * density(32) 
  pd(21,29) = pd(21,29) + rrt(344) * density(37) 
  pd(21,37) = pd(21,37) + rrt(344) * density(29) 
  pd(29,29) = pd(29,29) - rrt(344) * density(37) 
  pd(29,37) = pd(29,37) - rrt(344) * density(29) 
  pd(36,29) = pd(36,29) + rrt(344) * density(37) 
  pd(36,37) = pd(36,37) + rrt(344) * density(29) 
  pd(37,29) = pd(37,29) - rrt(344) * density(37) 
  pd(37,37) = pd(37,37) - rrt(344) * density(29) 
  pd(21,32) = pd(21,32) + rrt(345) * density(37) 
  pd(21,37) = pd(21,37) + rrt(345) * density(32) 
  pd(32,32) = pd(32,32) - rrt(345) * density(37) 
  pd(32,37) = pd(32,37) - rrt(345) * density(32) 
  pd(37,32) = pd(37,32) - rrt(345) * density(37) 
  pd(37,37) = pd(37,37) - rrt(345) * density(32) 
  pd(38,32) = pd(38,32) + rrt(345) * density(37) 
  pd(38,37) = pd(38,37) + rrt(345) * density(32) 
  pd(21,29) = pd(21,29) + rrt(346) * density(38) 
  pd(21,38) = pd(21,38) + rrt(346) * density(29) 
  pd(29,29) = pd(29,29) - rrt(346) * density(38) 
  pd(29,38) = pd(29,38) - rrt(346) * density(29) 
  pd(37,29) = pd(37,29) + rrt(346) * density(38) 
  pd(37,38) = pd(37,38) + rrt(346) * density(29) 
  pd(38,29) = pd(38,29) - rrt(346) * density(38) 
  pd(38,38) = pd(38,38) - rrt(346) * density(29) 
  pd(21,21) = pd(21,21) - rrt(347) * density(42) 
  pd(21,42) = pd(21,42) - rrt(347) * density(21) 
  pd(37,21) = pd(37,21) + rrt(347) * density(42) 
  pd(37,42) = pd(37,42) + rrt(347) * density(21) 
  pd(40,21) = pd(40,21) + rrt(347) * density(42) 
  pd(40,42) = pd(40,42) + rrt(347) * density(21) 
  pd(42,21) = pd(42,21) - rrt(347) * density(42) 
  pd(42,42) = pd(42,42) - rrt(347) * density(21) 
  pd(21,39) = pd(21,39) + rrt(348) 
  pd(37,39) = pd(37,39) + rrt(348) 
  pd(39,39) = pd(39,39) - rrt(348) 
  pd(21,29) = pd(21,29) + rrt(349) * density(39) 
  pd(21,39) = pd(21,39) + rrt(349) * density(29) 
  pd(29,29) = pd(29,29) - rrt(349) * density(39) 
  pd(29,39) = pd(29,39) - rrt(349) * density(29) 
  pd(38,29) = pd(38,29) + rrt(349) * density(39) 
  pd(38,39) = pd(38,39) + rrt(349) * density(29) 
  pd(39,29) = pd(39,29) - rrt(349) * density(39) 
  pd(39,39) = pd(39,39) - rrt(349) * density(29) 
  pd(21,29) = pd(21,29) + rrt(350) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(350) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(350) * density(39) 
  pd(29,39) = pd(29,39) - rrt(350) * density(29) 
  pd(36,29) = pd(36,29) + rrt(350) * density(39) 
  pd(36,39) = pd(36,39) + rrt(350) * density(29) 
  pd(39,29) = pd(39,29) - rrt(350) * density(39) 
  pd(39,39) = pd(39,39) - rrt(350) * density(29) 
  pd(21,26) = pd(21,26) + rrt(351) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(351) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(351) * density(39) 
  pd(26,39) = pd(26,39) - rrt(351) * density(26) 
  pd(37,26) = pd(37,26) + rrt(351) * density(39) 
  pd(37,39) = pd(37,39) + rrt(351) * density(26) 
  pd(39,26) = pd(39,26) - rrt(351) * density(39) 
  pd(39,39) = pd(39,39) - rrt(351) * density(26) 
  pd(21,27) = pd(21,27) + rrt(352) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(352) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(352) * density(39) 
  pd(27,39) = pd(27,39) - rrt(352) * density(27) 
  pd(37,27) = pd(37,27) + rrt(352) * density(39) 
  pd(37,39) = pd(37,39) + rrt(352) * density(27) 
  pd(39,27) = pd(39,27) - rrt(352) * density(39) 
  pd(39,39) = pd(39,39) - rrt(352) * density(27) 
  pd(21,21) = pd(21,21) - rrt(353) * density(36) 
  pd(21,36) = pd(21,36) - rrt(353) * density(21) 
  pd(36,21) = pd(36,21) - rrt(353) * density(36) 
  pd(36,36) = pd(36,36) - rrt(353) * density(21) 
  pd(38,21) = pd(38,21) + rrt(353) * density(36) 
  pd(38,36) = pd(38,36) + rrt(353) * density(21) 
  pd(21,21) = pd(21,21) - rrt(354) * density(37) 
  pd(21,37) = pd(21,37) - rrt(354) * density(21) 
  pd(37,21) = pd(37,21) - rrt(354) * density(37) 
  pd(37,37) = pd(37,37) - rrt(354) * density(21) 
  pd(39,21) = pd(39,21) + rrt(354) * density(37) 
  pd(39,37) = pd(39,37) + rrt(354) * density(21) 
  pd(14,17) = pd(14,17) + rrt(355) * density(36) 
  pd(14,36) = pd(14,36) + rrt(355) * density(17) 
  pd(17,17) = pd(17,17) - rrt(355) * density(36) 
  pd(17,36) = pd(17,36) - rrt(355) * density(17) 
  pd(29,17) = pd(29,17) + rrt(355) * density(36) 
  pd(29,36) = pd(29,36) + rrt(355) * density(17) 
  pd(36,17) = pd(36,17) - rrt(355) * density(36) 
  pd(36,36) = pd(36,36) - rrt(355) * density(17) 
  pd(01,18) = pd(01,18) + rrt(356) * density(36) 
  pd(01,36) = pd(01,36) + rrt(356) * density(18) 
  pd(18,18) = pd(18,18) - rrt(356) * density(36) 
  pd(18,36) = pd(18,36) - rrt(356) * density(18) 
  pd(29,18) = pd(29,18) + rrt(356) * density(36) 
  pd(29,36) = pd(29,36) + rrt(356) * density(18) 
  pd(36,18) = pd(36,18) - rrt(356) * density(36) 
  pd(36,36) = pd(36,36) - rrt(356) * density(18) 
  pd(29,33) = pd(29,33) + rrt(357) * density(36) * 2.0d0
  pd(29,36) = pd(29,36) + rrt(357) * density(33) * 2.0d0
  pd(33,33) = pd(33,33) - rrt(357) * density(36) 
  pd(33,36) = pd(33,36) - rrt(357) * density(33) 
  pd(36,33) = pd(36,33) - rrt(357) * density(36) 
  pd(36,36) = pd(36,36) - rrt(357) * density(33) 
  pd(21,34) = pd(21,34) + rrt(358) * density(36) 
  pd(21,36) = pd(21,36) + rrt(358) * density(34) 
  pd(29,34) = pd(29,34) + rrt(358) * density(36) 
  pd(29,36) = pd(29,36) + rrt(358) * density(34) 
  pd(34,34) = pd(34,34) - rrt(358) * density(36) 
  pd(34,36) = pd(34,36) - rrt(358) * density(34) 
  pd(36,34) = pd(36,34) - rrt(358) * density(36) 
  pd(36,36) = pd(36,36) - rrt(358) * density(34) 
  pd(29,36) = pd(29,36) + rrt(359) * density(41) 
  pd(29,41) = pd(29,41) + rrt(359) * density(36) 
  pd(36,36) = pd(36,36) - rrt(359) * density(41) 
  pd(36,41) = pd(36,41) - rrt(359) * density(36) 
  pd(40,36) = pd(40,36) + rrt(359) * density(41) 
  pd(40,41) = pd(40,41) + rrt(359) * density(36) 
  pd(41,36) = pd(41,36) - rrt(359) * density(41) 
  pd(41,41) = pd(41,41) - rrt(359) * density(36) 
  pd(14,17) = pd(14,17) + rrt(360) * density(37) 
  pd(14,37) = pd(14,37) + rrt(360) * density(17) 
  pd(17,17) = pd(17,17) - rrt(360) * density(37) 
  pd(17,37) = pd(17,37) - rrt(360) * density(17) 
  pd(21,17) = pd(21,17) + rrt(360) * density(37) 
  pd(21,37) = pd(21,37) + rrt(360) * density(17) 
  pd(37,17) = pd(37,17) - rrt(360) * density(37) 
  pd(37,37) = pd(37,37) - rrt(360) * density(17) 
  pd(01,18) = pd(01,18) + rrt(361) * density(37) 
  pd(01,37) = pd(01,37) + rrt(361) * density(18) 
  pd(18,18) = pd(18,18) - rrt(361) * density(37) 
  pd(18,37) = pd(18,37) - rrt(361) * density(18) 
  pd(21,18) = pd(21,18) + rrt(361) * density(37) 
  pd(21,37) = pd(21,37) + rrt(361) * density(18) 
  pd(37,18) = pd(37,18) - rrt(361) * density(37) 
  pd(37,37) = pd(37,37) - rrt(361) * density(18) 
  pd(21,33) = pd(21,33) + rrt(362) * density(37) 
  pd(21,37) = pd(21,37) + rrt(362) * density(33) 
  pd(29,33) = pd(29,33) + rrt(362) * density(37) 
  pd(29,37) = pd(29,37) + rrt(362) * density(33) 
  pd(33,33) = pd(33,33) - rrt(362) * density(37) 
  pd(33,37) = pd(33,37) - rrt(362) * density(33) 
  pd(37,33) = pd(37,33) - rrt(362) * density(37) 
  pd(37,37) = pd(37,37) - rrt(362) * density(33) 
  pd(21,34) = pd(21,34) + rrt(363) * density(37) * 2.0d0
  pd(21,37) = pd(21,37) + rrt(363) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(363) * density(37) 
  pd(34,37) = pd(34,37) - rrt(363) * density(34) 
  pd(37,34) = pd(37,34) - rrt(363) * density(37) 
  pd(37,37) = pd(37,37) - rrt(363) * density(34) 
  pd(21,37) = pd(21,37) + rrt(364) * density(41) 
  pd(21,41) = pd(21,41) + rrt(364) * density(37) 
  pd(37,37) = pd(37,37) - rrt(364) * density(41) 
  pd(37,41) = pd(37,41) - rrt(364) * density(37) 
  pd(40,37) = pd(40,37) + rrt(364) * density(41) 
  pd(40,41) = pd(40,41) + rrt(364) * density(37) 
  pd(41,37) = pd(41,37) - rrt(364) * density(41) 
  pd(41,41) = pd(41,41) - rrt(364) * density(37) 
  pd(14,17) = pd(14,17) + rrt(365) * density(38) 
  pd(14,38) = pd(14,38) + rrt(365) * density(17) 
  pd(17,17) = pd(17,17) - rrt(365) * density(38) 
  pd(17,38) = pd(17,38) - rrt(365) * density(17) 
  pd(32,17) = pd(32,17) + rrt(365) * density(38) 
  pd(32,38) = pd(32,38) + rrt(365) * density(17) 
  pd(38,17) = pd(38,17) - rrt(365) * density(38) 
  pd(38,38) = pd(38,38) - rrt(365) * density(17) 
  pd(01,18) = pd(01,18) + rrt(366) * density(38) 
  pd(01,38) = pd(01,38) + rrt(366) * density(18) 
  pd(18,18) = pd(18,18) - rrt(366) * density(38) 
  pd(18,38) = pd(18,38) - rrt(366) * density(18) 
  pd(32,18) = pd(32,18) + rrt(366) * density(38) 
  pd(32,38) = pd(32,38) + rrt(366) * density(18) 
  pd(38,18) = pd(38,18) - rrt(366) * density(38) 
  pd(38,38) = pd(38,38) - rrt(366) * density(18) 
  pd(29,33) = pd(29,33) + rrt(367) * density(38) 
  pd(29,38) = pd(29,38) + rrt(367) * density(33) 
  pd(32,33) = pd(32,33) + rrt(367) * density(38) 
  pd(32,38) = pd(32,38) + rrt(367) * density(33) 
  pd(33,33) = pd(33,33) - rrt(367) * density(38) 
  pd(33,38) = pd(33,38) - rrt(367) * density(33) 
  pd(38,33) = pd(38,33) - rrt(367) * density(38) 
  pd(38,38) = pd(38,38) - rrt(367) * density(33) 
  pd(21,34) = pd(21,34) + rrt(368) * density(38) 
  pd(21,38) = pd(21,38) + rrt(368) * density(34) 
  pd(32,34) = pd(32,34) + rrt(368) * density(38) 
  pd(32,38) = pd(32,38) + rrt(368) * density(34) 
  pd(34,34) = pd(34,34) - rrt(368) * density(38) 
  pd(34,38) = pd(34,38) - rrt(368) * density(34) 
  pd(38,34) = pd(38,34) - rrt(368) * density(38) 
  pd(38,38) = pd(38,38) - rrt(368) * density(34) 
  pd(32,38) = pd(32,38) + rrt(369) * density(41) 
  pd(32,41) = pd(32,41) + rrt(369) * density(38) 
  pd(38,38) = pd(38,38) - rrt(369) * density(41) 
  pd(38,41) = pd(38,41) - rrt(369) * density(38) 
  pd(40,38) = pd(40,38) + rrt(369) * density(41) 
  pd(40,41) = pd(40,41) + rrt(369) * density(38) 
  pd(41,38) = pd(41,38) - rrt(369) * density(41) 
  pd(41,41) = pd(41,41) - rrt(369) * density(38) 
  pd(14,17) = pd(14,17) + rrt(370) * density(42) 
  pd(14,42) = pd(14,42) + rrt(370) * density(17) 
  pd(17,17) = pd(17,17) - rrt(370) * density(42) 
  pd(17,42) = pd(17,42) - rrt(370) * density(17) 
  pd(40,17) = pd(40,17) + rrt(370) * density(42) 
  pd(40,42) = pd(40,42) + rrt(370) * density(17) 
  pd(42,17) = pd(42,17) - rrt(370) * density(42) 
  pd(42,42) = pd(42,42) - rrt(370) * density(17) 
  pd(01,18) = pd(01,18) + rrt(371) * density(42) 
  pd(01,42) = pd(01,42) + rrt(371) * density(18) 
  pd(18,18) = pd(18,18) - rrt(371) * density(42) 
  pd(18,42) = pd(18,42) - rrt(371) * density(18) 
  pd(40,18) = pd(40,18) + rrt(371) * density(42) 
  pd(40,42) = pd(40,42) + rrt(371) * density(18) 
  pd(42,18) = pd(42,18) - rrt(371) * density(42) 
  pd(42,42) = pd(42,42) - rrt(371) * density(18) 
  pd(29,33) = pd(29,33) + rrt(372) * density(42) 
  pd(29,42) = pd(29,42) + rrt(372) * density(33) 
  pd(33,33) = pd(33,33) - rrt(372) * density(42) 
  pd(33,42) = pd(33,42) - rrt(372) * density(33) 
  pd(40,33) = pd(40,33) + rrt(372) * density(42) 
  pd(40,42) = pd(40,42) + rrt(372) * density(33) 
  pd(42,33) = pd(42,33) - rrt(372) * density(42) 
  pd(42,42) = pd(42,42) - rrt(372) * density(33) 
  pd(21,34) = pd(21,34) + rrt(373) * density(42) 
  pd(21,42) = pd(21,42) + rrt(373) * density(34) 
  pd(34,34) = pd(34,34) - rrt(373) * density(42) 
  pd(34,42) = pd(34,42) - rrt(373) * density(34) 
  pd(40,34) = pd(40,34) + rrt(373) * density(42) 
  pd(40,42) = pd(40,42) + rrt(373) * density(34) 
  pd(42,34) = pd(42,34) - rrt(373) * density(42) 
  pd(42,42) = pd(42,42) - rrt(373) * density(34) 
  pd(40,41) = pd(40,41) + rrt(374) * density(42) * 2.0d0
  pd(40,42) = pd(40,42) + rrt(374) * density(41) * 2.0d0
  pd(41,41) = pd(41,41) - rrt(374) * density(42) 
  pd(41,42) = pd(41,42) - rrt(374) * density(41) 
  pd(42,41) = pd(42,41) - rrt(374) * density(42) 
  pd(42,42) = pd(42,42) - rrt(374) * density(41) 
  pd(14,18) = pd(14,18) + rrt(375) * density(36) * 2.0d0
  pd(14,36) = pd(14,36) + rrt(375) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(375) * density(36) 
  pd(18,36) = pd(18,36) - rrt(375) * density(18) 
  pd(29,18) = pd(29,18) + rrt(375) * density(36) 
  pd(29,36) = pd(29,36) + rrt(375) * density(18) 
  pd(36,18) = pd(36,18) - rrt(375) * density(36) 
  pd(36,36) = pd(36,36) - rrt(375) * density(18) 
  pd(01,19) = pd(01,19) + rrt(376) * density(36) 
  pd(01,36) = pd(01,36) + rrt(376) * density(19) 
  pd(14,19) = pd(14,19) + rrt(376) * density(36) 
  pd(14,36) = pd(14,36) + rrt(376) * density(19) 
  pd(19,19) = pd(19,19) - rrt(376) * density(36) 
  pd(19,36) = pd(19,36) - rrt(376) * density(19) 
  pd(29,19) = pd(29,19) + rrt(376) * density(36) 
  pd(29,36) = pd(29,36) + rrt(376) * density(19) 
  pd(36,19) = pd(36,19) - rrt(376) * density(36) 
  pd(36,36) = pd(36,36) - rrt(376) * density(19) 
  pd(01,20) = pd(01,20) + rrt(377) * density(36) * 2.0d0
  pd(01,36) = pd(01,36) + rrt(377) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(377) * density(36) 
  pd(20,36) = pd(20,36) - rrt(377) * density(20) 
  pd(29,20) = pd(29,20) + rrt(377) * density(36) 
  pd(29,36) = pd(29,36) + rrt(377) * density(20) 
  pd(36,20) = pd(36,20) - rrt(377) * density(36) 
  pd(36,36) = pd(36,36) - rrt(377) * density(20) 
  pd(29,34) = pd(29,34) + rrt(378) * density(36) * 3.0d0
  pd(29,36) = pd(29,36) + rrt(378) * density(34) * 3.0d0
  pd(34,34) = pd(34,34) - rrt(378) * density(36) 
  pd(34,36) = pd(34,36) - rrt(378) * density(34) 
  pd(36,34) = pd(36,34) - rrt(378) * density(36) 
  pd(36,36) = pd(36,36) - rrt(378) * density(34) 
  pd(21,35) = pd(21,35) + rrt(379) * density(36) * 2.0d0
  pd(21,36) = pd(21,36) + rrt(379) * density(35) * 2.0d0
  pd(29,35) = pd(29,35) + rrt(379) * density(36) 
  pd(29,36) = pd(29,36) + rrt(379) * density(35) 
  pd(35,35) = pd(35,35) - rrt(379) * density(36) 
  pd(35,36) = pd(35,36) - rrt(379) * density(35) 
  pd(36,35) = pd(36,35) - rrt(379) * density(36) 
  pd(36,36) = pd(36,36) - rrt(379) * density(35) 
  pd(14,36) = pd(14,36) + rrt(380) * density(41) 
  pd(14,41) = pd(14,41) + rrt(380) * density(36) 
  pd(29,36) = pd(29,36) + rrt(380) * density(41) * 2.0d0
  pd(29,41) = pd(29,41) + rrt(380) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(380) * density(41) 
  pd(36,41) = pd(36,41) - rrt(380) * density(36) 
  pd(41,36) = pd(41,36) - rrt(380) * density(41) 
  pd(41,41) = pd(41,41) - rrt(380) * density(36) 
  pd(01,36) = pd(01,36) + rrt(381) * density(43) 
  pd(01,43) = pd(01,43) + rrt(381) * density(36) 
  pd(21,36) = pd(21,36) + rrt(381) * density(43) 
  pd(21,43) = pd(21,43) + rrt(381) * density(36) 
  pd(29,36) = pd(29,36) + rrt(381) * density(43) 
  pd(29,43) = pd(29,43) + rrt(381) * density(36) 
  pd(36,36) = pd(36,36) - rrt(381) * density(43) 
  pd(36,43) = pd(36,43) - rrt(381) * density(36) 
  pd(43,36) = pd(43,36) - rrt(381) * density(43) 
  pd(43,43) = pd(43,43) - rrt(381) * density(36) 
  pd(14,18) = pd(14,18) + rrt(382) * density(37) * 2.0d0
  pd(14,37) = pd(14,37) + rrt(382) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(382) * density(37) 
  pd(18,37) = pd(18,37) - rrt(382) * density(18) 
  pd(21,18) = pd(21,18) + rrt(382) * density(37) 
  pd(21,37) = pd(21,37) + rrt(382) * density(18) 
  pd(37,18) = pd(37,18) - rrt(382) * density(37) 
  pd(37,37) = pd(37,37) - rrt(382) * density(18) 
  pd(01,19) = pd(01,19) + rrt(383) * density(37) 
  pd(01,37) = pd(01,37) + rrt(383) * density(19) 
  pd(14,19) = pd(14,19) + rrt(383) * density(37) 
  pd(14,37) = pd(14,37) + rrt(383) * density(19) 
  pd(19,19) = pd(19,19) - rrt(383) * density(37) 
  pd(19,37) = pd(19,37) - rrt(383) * density(19) 
  pd(21,19) = pd(21,19) + rrt(383) * density(37) 
  pd(21,37) = pd(21,37) + rrt(383) * density(19) 
  pd(37,19) = pd(37,19) - rrt(383) * density(37) 
  pd(37,37) = pd(37,37) - rrt(383) * density(19) 
  pd(01,20) = pd(01,20) + rrt(384) * density(37) * 2.0d0
  pd(01,37) = pd(01,37) + rrt(384) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(384) * density(37) 
  pd(20,37) = pd(20,37) - rrt(384) * density(20) 
  pd(21,20) = pd(21,20) + rrt(384) * density(37) 
  pd(21,37) = pd(21,37) + rrt(384) * density(20) 
  pd(37,20) = pd(37,20) - rrt(384) * density(37) 
  pd(37,37) = pd(37,37) - rrt(384) * density(20) 
  pd(21,34) = pd(21,34) + rrt(385) * density(37) 
  pd(21,37) = pd(21,37) + rrt(385) * density(34) 
  pd(29,34) = pd(29,34) + rrt(385) * density(37) * 2.0d0
  pd(29,37) = pd(29,37) + rrt(385) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(385) * density(37) 
  pd(34,37) = pd(34,37) - rrt(385) * density(34) 
  pd(37,34) = pd(37,34) - rrt(385) * density(37) 
  pd(37,37) = pd(37,37) - rrt(385) * density(34) 
  pd(21,35) = pd(21,35) + rrt(386) * density(37) * 3.0d0
  pd(21,37) = pd(21,37) + rrt(386) * density(35) * 3.0d0
  pd(35,35) = pd(35,35) - rrt(386) * density(37) 
  pd(35,37) = pd(35,37) - rrt(386) * density(35) 
  pd(37,35) = pd(37,35) - rrt(386) * density(37) 
  pd(37,37) = pd(37,37) - rrt(386) * density(35) 
  pd(14,37) = pd(14,37) + rrt(387) * density(41) 
  pd(14,41) = pd(14,41) + rrt(387) * density(37) 
  pd(21,37) = pd(21,37) + rrt(387) * density(41) 
  pd(21,41) = pd(21,41) + rrt(387) * density(37) 
  pd(29,37) = pd(29,37) + rrt(387) * density(41) 
  pd(29,41) = pd(29,41) + rrt(387) * density(37) 
  pd(37,37) = pd(37,37) - rrt(387) * density(41) 
  pd(37,41) = pd(37,41) - rrt(387) * density(37) 
  pd(41,37) = pd(41,37) - rrt(387) * density(41) 
  pd(41,41) = pd(41,41) - rrt(387) * density(37) 
  pd(01,37) = pd(01,37) + rrt(388) * density(43) 
  pd(01,43) = pd(01,43) + rrt(388) * density(37) 
  pd(21,37) = pd(21,37) + rrt(388) * density(43) * 2.0d0
  pd(21,43) = pd(21,43) + rrt(388) * density(37) * 2.0d0
  pd(37,37) = pd(37,37) - rrt(388) * density(43) 
  pd(37,43) = pd(37,43) - rrt(388) * density(37) 
  pd(43,37) = pd(43,37) - rrt(388) * density(43) 
  pd(43,43) = pd(43,43) - rrt(388) * density(37) 
  pd(14,18) = pd(14,18) + rrt(389) * density(38) * 2.0d0
  pd(14,38) = pd(14,38) + rrt(389) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(389) * density(38) 
  pd(18,38) = pd(18,38) - rrt(389) * density(18) 
  pd(32,18) = pd(32,18) + rrt(389) * density(38) 
  pd(32,38) = pd(32,38) + rrt(389) * density(18) 
  pd(38,18) = pd(38,18) - rrt(389) * density(38) 
  pd(38,38) = pd(38,38) - rrt(389) * density(18) 
  pd(01,19) = pd(01,19) + rrt(390) * density(38) 
  pd(01,38) = pd(01,38) + rrt(390) * density(19) 
  pd(14,19) = pd(14,19) + rrt(390) * density(38) 
  pd(14,38) = pd(14,38) + rrt(390) * density(19) 
  pd(19,19) = pd(19,19) - rrt(390) * density(38) 
  pd(19,38) = pd(19,38) - rrt(390) * density(19) 
  pd(32,19) = pd(32,19) + rrt(390) * density(38) 
  pd(32,38) = pd(32,38) + rrt(390) * density(19) 
  pd(38,19) = pd(38,19) - rrt(390) * density(38) 
  pd(38,38) = pd(38,38) - rrt(390) * density(19) 
  pd(01,20) = pd(01,20) + rrt(391) * density(38) * 2.0d0
  pd(01,38) = pd(01,38) + rrt(391) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(391) * density(38) 
  pd(20,38) = pd(20,38) - rrt(391) * density(20) 
  pd(32,20) = pd(32,20) + rrt(391) * density(38) 
  pd(32,38) = pd(32,38) + rrt(391) * density(20) 
  pd(38,20) = pd(38,20) - rrt(391) * density(38) 
  pd(38,38) = pd(38,38) - rrt(391) * density(20) 
  pd(29,34) = pd(29,34) + rrt(392) * density(38) * 2.0d0
  pd(29,38) = pd(29,38) + rrt(392) * density(34) * 2.0d0
  pd(32,34) = pd(32,34) + rrt(392) * density(38) 
  pd(32,38) = pd(32,38) + rrt(392) * density(34) 
  pd(34,34) = pd(34,34) - rrt(392) * density(38) 
  pd(34,38) = pd(34,38) - rrt(392) * density(34) 
  pd(38,34) = pd(38,34) - rrt(392) * density(38) 
  pd(38,38) = pd(38,38) - rrt(392) * density(34) 
  pd(21,35) = pd(21,35) + rrt(393) * density(38) * 2.0d0
  pd(21,38) = pd(21,38) + rrt(393) * density(35) * 2.0d0
  pd(32,35) = pd(32,35) + rrt(393) * density(38) 
  pd(32,38) = pd(32,38) + rrt(393) * density(35) 
  pd(35,35) = pd(35,35) - rrt(393) * density(38) 
  pd(35,38) = pd(35,38) - rrt(393) * density(35) 
  pd(38,35) = pd(38,35) - rrt(393) * density(38) 
  pd(38,38) = pd(38,38) - rrt(393) * density(35) 
  pd(14,38) = pd(14,38) + rrt(394) * density(41) 
  pd(14,41) = pd(14,41) + rrt(394) * density(38) 
  pd(29,38) = pd(29,38) + rrt(394) * density(41) 
  pd(29,41) = pd(29,41) + rrt(394) * density(38) 
  pd(32,38) = pd(32,38) + rrt(394) * density(41) 
  pd(32,41) = pd(32,41) + rrt(394) * density(38) 
  pd(38,38) = pd(38,38) - rrt(394) * density(41) 
  pd(38,41) = pd(38,41) - rrt(394) * density(38) 
  pd(41,38) = pd(41,38) - rrt(394) * density(41) 
  pd(41,41) = pd(41,41) - rrt(394) * density(38) 
  pd(01,38) = pd(01,38) + rrt(395) * density(43) 
  pd(01,43) = pd(01,43) + rrt(395) * density(38) 
  pd(21,38) = pd(21,38) + rrt(395) * density(43) 
  pd(21,43) = pd(21,43) + rrt(395) * density(38) 
  pd(32,38) = pd(32,38) + rrt(395) * density(43) 
  pd(32,43) = pd(32,43) + rrt(395) * density(38) 
  pd(38,38) = pd(38,38) - rrt(395) * density(43) 
  pd(38,43) = pd(38,43) - rrt(395) * density(38) 
  pd(43,38) = pd(43,38) - rrt(395) * density(43) 
  pd(43,43) = pd(43,43) - rrt(395) * density(38) 
  pd(14,18) = pd(14,18) + rrt(396) * density(42) * 2.0d0
  pd(14,42) = pd(14,42) + rrt(396) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(396) * density(42) 
  pd(18,42) = pd(18,42) - rrt(396) * density(18) 
  pd(40,18) = pd(40,18) + rrt(396) * density(42) 
  pd(40,42) = pd(40,42) + rrt(396) * density(18) 
  pd(42,18) = pd(42,18) - rrt(396) * density(42) 
  pd(42,42) = pd(42,42) - rrt(396) * density(18) 
  pd(01,19) = pd(01,19) + rrt(397) * density(42) 
  pd(01,42) = pd(01,42) + rrt(397) * density(19) 
  pd(14,19) = pd(14,19) + rrt(397) * density(42) 
  pd(14,42) = pd(14,42) + rrt(397) * density(19) 
  pd(19,19) = pd(19,19) - rrt(397) * density(42) 
  pd(19,42) = pd(19,42) - rrt(397) * density(19) 
  pd(40,19) = pd(40,19) + rrt(397) * density(42) 
  pd(40,42) = pd(40,42) + rrt(397) * density(19) 
  pd(42,19) = pd(42,19) - rrt(397) * density(42) 
  pd(42,42) = pd(42,42) - rrt(397) * density(19) 
  pd(01,20) = pd(01,20) + rrt(398) * density(42) * 2.0d0
  pd(01,42) = pd(01,42) + rrt(398) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(398) * density(42) 
  pd(20,42) = pd(20,42) - rrt(398) * density(20) 
  pd(40,20) = pd(40,20) + rrt(398) * density(42) 
  pd(40,42) = pd(40,42) + rrt(398) * density(20) 
  pd(42,20) = pd(42,20) - rrt(398) * density(42) 
  pd(42,42) = pd(42,42) - rrt(398) * density(20) 
  pd(29,34) = pd(29,34) + rrt(399) * density(42) * 2.0d0
  pd(29,42) = pd(29,42) + rrt(399) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(399) * density(42) 
  pd(34,42) = pd(34,42) - rrt(399) * density(34) 
  pd(40,34) = pd(40,34) + rrt(399) * density(42) 
  pd(40,42) = pd(40,42) + rrt(399) * density(34) 
  pd(42,34) = pd(42,34) - rrt(399) * density(42) 
  pd(42,42) = pd(42,42) - rrt(399) * density(34) 
  pd(21,35) = pd(21,35) + rrt(400) * density(42) * 2.0d0
  pd(21,42) = pd(21,42) + rrt(400) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(400) * density(42) 
  pd(35,42) = pd(35,42) - rrt(400) * density(35) 
  pd(40,35) = pd(40,35) + rrt(400) * density(42) 
  pd(40,42) = pd(40,42) + rrt(400) * density(35) 
  pd(42,35) = pd(42,35) - rrt(400) * density(42) 
  pd(42,42) = pd(42,42) - rrt(400) * density(35) 
  pd(14,41) = pd(14,41) + rrt(401) * density(42) 
  pd(14,42) = pd(14,42) + rrt(401) * density(41) 
  pd(29,41) = pd(29,41) + rrt(401) * density(42) 
  pd(29,42) = pd(29,42) + rrt(401) * density(41) 
  pd(40,41) = pd(40,41) + rrt(401) * density(42) 
  pd(40,42) = pd(40,42) + rrt(401) * density(41) 
  pd(41,41) = pd(41,41) - rrt(401) * density(42) 
  pd(41,42) = pd(41,42) - rrt(401) * density(41) 
  pd(42,41) = pd(42,41) - rrt(401) * density(42) 
  pd(42,42) = pd(42,42) - rrt(401) * density(41) 
  pd(01,42) = pd(01,42) + rrt(402) * density(43) 
  pd(01,43) = pd(01,43) + rrt(402) * density(42) 
  pd(21,42) = pd(21,42) + rrt(402) * density(43) 
  pd(21,43) = pd(21,43) + rrt(402) * density(42) 
  pd(40,42) = pd(40,42) + rrt(402) * density(43) 
  pd(40,43) = pd(40,43) + rrt(402) * density(42) 
  pd(42,42) = pd(42,42) - rrt(402) * density(43) 
  pd(42,43) = pd(42,43) - rrt(402) * density(42) 
  pd(43,42) = pd(43,42) - rrt(402) * density(43) 
  pd(43,43) = pd(43,43) - rrt(402) * density(42) 
  pd(14,17) = pd(14,17) + rrt(403) * density(39) 
  pd(14,39) = pd(14,39) + rrt(403) * density(17) 
  pd(17,17) = pd(17,17) - rrt(403) * density(39) 
  pd(17,39) = pd(17,39) - rrt(403) * density(17) 
  pd(21,17) = pd(21,17) + rrt(403) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(403) * density(17) * 2.0d0
  pd(39,17) = pd(39,17) - rrt(403) * density(39) 
  pd(39,39) = pd(39,39) - rrt(403) * density(17) 
  pd(01,18) = pd(01,18) + rrt(404) * density(39) 
  pd(01,39) = pd(01,39) + rrt(404) * density(18) 
  pd(18,18) = pd(18,18) - rrt(404) * density(39) 
  pd(18,39) = pd(18,39) - rrt(404) * density(18) 
  pd(21,18) = pd(21,18) + rrt(404) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(404) * density(18) * 2.0d0
  pd(39,18) = pd(39,18) - rrt(404) * density(39) 
  pd(39,39) = pd(39,39) - rrt(404) * density(18) 
  pd(01,19) = pd(01,19) + rrt(405) * density(39) 
  pd(01,39) = pd(01,39) + rrt(405) * density(19) 
  pd(14,19) = pd(14,19) + rrt(405) * density(39) 
  pd(14,39) = pd(14,39) + rrt(405) * density(19) 
  pd(19,19) = pd(19,19) - rrt(405) * density(39) 
  pd(19,39) = pd(19,39) - rrt(405) * density(19) 
  pd(21,19) = pd(21,19) + rrt(405) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(405) * density(19) * 2.0d0
  pd(39,19) = pd(39,19) - rrt(405) * density(39) 
  pd(39,39) = pd(39,39) - rrt(405) * density(19) 
  pd(01,20) = pd(01,20) + rrt(406) * density(39) * 2.0d0
  pd(01,39) = pd(01,39) + rrt(406) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(406) * density(39) 
  pd(20,39) = pd(20,39) - rrt(406) * density(20) 
  pd(21,20) = pd(21,20) + rrt(406) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(406) * density(20) * 2.0d0
  pd(39,20) = pd(39,20) - rrt(406) * density(39) 
  pd(39,39) = pd(39,39) - rrt(406) * density(20) 
  pd(21,33) = pd(21,33) + rrt(407) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(407) * density(33) * 2.0d0
  pd(29,33) = pd(29,33) + rrt(407) * density(39) 
  pd(29,39) = pd(29,39) + rrt(407) * density(33) 
  pd(33,33) = pd(33,33) - rrt(407) * density(39) 
  pd(33,39) = pd(33,39) - rrt(407) * density(33) 
  pd(39,33) = pd(39,33) - rrt(407) * density(39) 
  pd(39,39) = pd(39,39) - rrt(407) * density(33) 
  pd(21,34) = pd(21,34) + rrt(408) * density(39) * 3.0d0
  pd(21,39) = pd(21,39) + rrt(408) * density(34) * 3.0d0
  pd(34,34) = pd(34,34) - rrt(408) * density(39) 
  pd(34,39) = pd(34,39) - rrt(408) * density(34) 
  pd(39,34) = pd(39,34) - rrt(408) * density(39) 
  pd(39,39) = pd(39,39) - rrt(408) * density(34) 
  pd(21,35) = pd(21,35) + rrt(409) * density(39) * 4.0d0
  pd(21,39) = pd(21,39) + rrt(409) * density(35) * 4.0d0
  pd(35,35) = pd(35,35) - rrt(409) * density(39) 
  pd(35,39) = pd(35,39) - rrt(409) * density(35) 
  pd(39,35) = pd(39,35) - rrt(409) * density(39) 
  pd(39,39) = pd(39,39) - rrt(409) * density(35) 
  pd(21,39) = pd(21,39) + rrt(410) * density(41) * 2.0d0
  pd(21,41) = pd(21,41) + rrt(410) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(410) * density(41) 
  pd(39,41) = pd(39,41) - rrt(410) * density(39) 
  pd(40,39) = pd(40,39) + rrt(410) * density(41) 
  pd(40,41) = pd(40,41) + rrt(410) * density(39) 
  pd(41,39) = pd(41,39) - rrt(410) * density(41) 
  pd(41,41) = pd(41,41) - rrt(410) * density(39) 
  pd(01,39) = pd(01,39) + rrt(411) * density(43) 
  pd(01,43) = pd(01,43) + rrt(411) * density(39) 
  pd(21,39) = pd(21,39) + rrt(411) * density(43) * 3.0d0
  pd(21,43) = pd(21,43) + rrt(411) * density(39) * 3.0d0
  pd(39,39) = pd(39,39) - rrt(411) * density(43) 
  pd(39,43) = pd(39,43) - rrt(411) * density(39) 
  pd(43,39) = pd(43,39) - rrt(411) * density(43) 
  pd(43,43) = pd(43,43) - rrt(411) * density(39) 
  pd(14,17) = pd(14,17) + rrt(412) * density(36) 
  pd(14,36) = pd(14,36) + rrt(412) * density(17) 
  pd(17,17) = pd(17,17) - rrt(412) * density(36) 
  pd(17,36) = pd(17,36) - rrt(412) * density(17) 
  pd(29,17) = pd(29,17) + rrt(412) * density(36) 
  pd(29,36) = pd(29,36) + rrt(412) * density(17) 
  pd(36,17) = pd(36,17) - rrt(412) * density(36) 
  pd(36,36) = pd(36,36) - rrt(412) * density(17) 
  pd(01,18) = pd(01,18) + rrt(413) * density(36) 
  pd(01,36) = pd(01,36) + rrt(413) * density(18) 
  pd(18,18) = pd(18,18) - rrt(413) * density(36) 
  pd(18,36) = pd(18,36) - rrt(413) * density(18) 
  pd(29,18) = pd(29,18) + rrt(413) * density(36) 
  pd(29,36) = pd(29,36) + rrt(413) * density(18) 
  pd(36,18) = pd(36,18) - rrt(413) * density(36) 
  pd(36,36) = pd(36,36) - rrt(413) * density(18) 
  pd(29,33) = pd(29,33) + rrt(414) * density(36) * 2.0d0
  pd(29,36) = pd(29,36) + rrt(414) * density(33) * 2.0d0
  pd(33,33) = pd(33,33) - rrt(414) * density(36) 
  pd(33,36) = pd(33,36) - rrt(414) * density(33) 
  pd(36,33) = pd(36,33) - rrt(414) * density(36) 
  pd(36,36) = pd(36,36) - rrt(414) * density(33) 
  pd(21,34) = pd(21,34) + rrt(415) * density(36) 
  pd(21,36) = pd(21,36) + rrt(415) * density(34) 
  pd(29,34) = pd(29,34) + rrt(415) * density(36) 
  pd(29,36) = pd(29,36) + rrt(415) * density(34) 
  pd(34,34) = pd(34,34) - rrt(415) * density(36) 
  pd(34,36) = pd(34,36) - rrt(415) * density(34) 
  pd(36,34) = pd(36,34) - rrt(415) * density(36) 
  pd(36,36) = pd(36,36) - rrt(415) * density(34) 
  pd(29,36) = pd(29,36) + rrt(416) * density(41) 
  pd(29,41) = pd(29,41) + rrt(416) * density(36) 
  pd(36,36) = pd(36,36) - rrt(416) * density(41) 
  pd(36,41) = pd(36,41) - rrt(416) * density(36) 
  pd(40,36) = pd(40,36) + rrt(416) * density(41) 
  pd(40,41) = pd(40,41) + rrt(416) * density(36) 
  pd(41,36) = pd(41,36) - rrt(416) * density(41) 
  pd(41,41) = pd(41,41) - rrt(416) * density(36) 
  pd(14,17) = pd(14,17) + rrt(417) * density(37) 
  pd(14,37) = pd(14,37) + rrt(417) * density(17) 
  pd(17,17) = pd(17,17) - rrt(417) * density(37) 
  pd(17,37) = pd(17,37) - rrt(417) * density(17) 
  pd(21,17) = pd(21,17) + rrt(417) * density(37) 
  pd(21,37) = pd(21,37) + rrt(417) * density(17) 
  pd(37,17) = pd(37,17) - rrt(417) * density(37) 
  pd(37,37) = pd(37,37) - rrt(417) * density(17) 
  pd(01,18) = pd(01,18) + rrt(418) * density(37) 
  pd(01,37) = pd(01,37) + rrt(418) * density(18) 
  pd(18,18) = pd(18,18) - rrt(418) * density(37) 
  pd(18,37) = pd(18,37) - rrt(418) * density(18) 
  pd(21,18) = pd(21,18) + rrt(418) * density(37) 
  pd(21,37) = pd(21,37) + rrt(418) * density(18) 
  pd(37,18) = pd(37,18) - rrt(418) * density(37) 
  pd(37,37) = pd(37,37) - rrt(418) * density(18) 
  pd(21,33) = pd(21,33) + rrt(419) * density(37) 
  pd(21,37) = pd(21,37) + rrt(419) * density(33) 
  pd(29,33) = pd(29,33) + rrt(419) * density(37) 
  pd(29,37) = pd(29,37) + rrt(419) * density(33) 
  pd(33,33) = pd(33,33) - rrt(419) * density(37) 
  pd(33,37) = pd(33,37) - rrt(419) * density(33) 
  pd(37,33) = pd(37,33) - rrt(419) * density(37) 
  pd(37,37) = pd(37,37) - rrt(419) * density(33) 
  pd(21,34) = pd(21,34) + rrt(420) * density(37) * 2.0d0
  pd(21,37) = pd(21,37) + rrt(420) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(420) * density(37) 
  pd(34,37) = pd(34,37) - rrt(420) * density(34) 
  pd(37,34) = pd(37,34) - rrt(420) * density(37) 
  pd(37,37) = pd(37,37) - rrt(420) * density(34) 
  pd(21,37) = pd(21,37) + rrt(421) * density(41) 
  pd(21,41) = pd(21,41) + rrt(421) * density(37) 
  pd(37,37) = pd(37,37) - rrt(421) * density(41) 
  pd(37,41) = pd(37,41) - rrt(421) * density(37) 
  pd(40,37) = pd(40,37) + rrt(421) * density(41) 
  pd(40,41) = pd(40,41) + rrt(421) * density(37) 
  pd(41,37) = pd(41,37) - rrt(421) * density(41) 
  pd(41,41) = pd(41,41) - rrt(421) * density(37) 
  pd(17,17) = pd(17,17) - rrt(422) * density(36) 
  pd(17,36) = pd(17,36) - rrt(422) * density(17) 
  pd(36,17) = pd(36,17) - rrt(422) * density(36) 
  pd(36,36) = pd(36,36) - rrt(422) * density(17) 
  pd(40,17) = pd(40,17) + rrt(422) * density(36) 
  pd(40,36) = pd(40,36) + rrt(422) * density(17) 
  pd(21,33) = pd(21,33) + rrt(423) * density(36) 
  pd(21,36) = pd(21,36) + rrt(423) * density(33) 
  pd(33,33) = pd(33,33) - rrt(423) * density(36) 
  pd(33,36) = pd(33,36) - rrt(423) * density(33) 
  pd(36,33) = pd(36,33) - rrt(423) * density(36) 
  pd(36,36) = pd(36,36) - rrt(423) * density(33) 
  pd(32,34) = pd(32,34) + rrt(424) * density(36) 
  pd(32,36) = pd(32,36) + rrt(424) * density(34) 
  pd(34,34) = pd(34,34) - rrt(424) * density(36) 
  pd(34,36) = pd(34,36) - rrt(424) * density(34) 
  pd(36,34) = pd(36,34) - rrt(424) * density(36) 
  pd(36,36) = pd(36,36) - rrt(424) * density(34) 
  pd(32,33) = pd(32,33) + rrt(425) * density(37) 
  pd(32,37) = pd(32,37) + rrt(425) * density(33) 
  pd(33,33) = pd(33,33) - rrt(425) * density(37) 
  pd(33,37) = pd(33,37) - rrt(425) * density(33) 
  pd(37,33) = pd(37,33) - rrt(425) * density(37) 
  pd(37,37) = pd(37,37) - rrt(425) * density(33) 
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(45,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(45,:) = pd(45,:) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_jex
end module ZDPlasKin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction constant rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_rates(Time)
  use ZDPlasKin, only : ZDPlasKin_bolsig_rates, bolsig_rates, bolsig_pointer, ZDPlasKin_cfg, ZDPlasKin_get_density_total, &
                        lreaction_block, rrt
  implicit none
  double precision, intent(in) :: Time
  double precision :: Tgas
  double precision :: EN
  double precision :: Te
  double precision :: ANY_NEUTRAL
  DOUBLE PRECISION :: DTION, TIONN, TIONN2, TIONN3, TIONN4, TEFFN, TEFFN2, TEFFN3, TEFFN4 ! K
  DOUBLE PRECISION, PARAMETER :: ENERGY_VIBN2 = 0.290D0*11605.0D0, ENERGY_VIBO2 = 0.190D0*11605.0D0 ! K
  DOUBLE PRECISION :: QVIBN2, KVT10_N2N2, KVT01_N2N2, KVT10_N2N, KVT01_N2N, KVT10_N2O, KVT01_N2O ! CM3.S-1
  DOUBLE PRECISION :: QVIBO2, KVT10_O2O2, KVT01_O2O2, KVT10_O2O, KVT01_O2O ! CM3.S-1
  call ZDPlasKin_bolsig_rates()
  Tgas = ZDPlasKin_cfg(1)
  EN  = ZDPlasKin_cfg(3)
  Te  = ZDPlasKin_cfg(4)
  call ZDPlasKin_get_density_total(ALL_NEUTRAL=ANY_NEUTRAL)
  DTION = 2.0D0 / ( 3.0D0 * 1.3807D-16 ) * 1.6605D-24 * ( 1.0D-17 * EN )**2
  TIONN = TGAS + DTION * 14.0D0 * 8.0D19**2
  TIONN2 = TGAS + DTION * 28.0D0 * 4.1D19**2
  TIONN3 = TGAS + DTION * 42.0D0 * 6.1D19**2
  TIONN4 = TGAS + DTION * 56.0D0 * 7.0D19**2
  TEFFN = ( TIONN + 0.5D0 * TGAS ) / ( 1.0D0 + 0.5D0 )
  TEFFN2 = ( TIONN2 + 1.0D0 * TGAS ) / ( 1.0D0 + 1.0D0 )
  TEFFN3 = ( TIONN3 + 1.5D0 * TGAS ) / ( 1.0D0 + 1.5D0 )
  TEFFN4 = ( TIONN4 + 2.0D0 * TGAS ) / ( 1.0D0 + 2.0D0 )
  QVIBN2 = EXP( - ENERGY_VIBN2 / TGAS )
  KVT10_N2N2 = 7.80D-12 * TGAS * EXP( - 218.0 / TGAS**(1.0/3.0) + 690.0 / TGAS ) / ( 1.0 - QVIBN2 )
  KVT10_N2N = 4.00D-16 * ( TGAS / 300.0D0 )**0.5
  KVT10_N2O = 1.20D-13 * EXP( - 27.6 / TGAS**(1.0/3.0) )
  KVT01_N2N2 = KVT10_N2N2 * QVIBN2
  KVT01_N2N = KVT10_N2N * QVIBN2
  KVT01_N2O = KVT10_N2O * QVIBN2
  QVIBO2 = EXP( - ENERGY_VIBO2 / TGAS )
  KVT10_O2O2 = 1.35D-12 * TGAS * EXP( - 137.9 / TGAS**(1.0/3.0) ) / ( 1.0 - QVIBO2 )
  KVT10_O2O = 4.50D-15 * TGAS
  KVT01_O2O2 = KVT10_O2O2 * QVIBO2
  KVT01_O2O = KVT10_O2O * QVIBO2
  rrt(001) = bolsig_rates(bolsig_pointer(1))
  rrt(002) = bolsig_rates(bolsig_pointer(2))
  rrt(003) = bolsig_rates(bolsig_pointer(3))
  rrt(004) = bolsig_rates(bolsig_pointer(4))
  rrt(005) = bolsig_rates(bolsig_pointer(5))
  rrt(006) = bolsig_rates(bolsig_pointer(6))
  rrt(007) = bolsig_rates(bolsig_pointer(7))
  rrt(008) = bolsig_rates(bolsig_pointer(8))
  rrt(009) = bolsig_rates(bolsig_pointer(9))
  rrt(010) = bolsig_rates(bolsig_pointer(10))
  rrt(011) = bolsig_rates(bolsig_pointer(11))
  rrt(012) = bolsig_rates(bolsig_pointer(12))
  rrt(013) = bolsig_rates(bolsig_pointer(13))
  rrt(014) = bolsig_rates(bolsig_pointer(14))
  rrt(015) = bolsig_rates(bolsig_pointer(15))
  rrt(016) = bolsig_rates(bolsig_pointer(16))
  rrt(017) = bolsig_rates(bolsig_pointer(17))
  rrt(018) = bolsig_rates(bolsig_pointer(18))
  rrt(019) = bolsig_rates(bolsig_pointer(19))
  rrt(020) = bolsig_rates(bolsig_pointer(20))
  rrt(021) = bolsig_rates(bolsig_pointer(21))
  rrt(022) = bolsig_rates(bolsig_pointer(22))
  rrt(023) = bolsig_rates(bolsig_pointer(23))
  rrt(024) = bolsig_rates(bolsig_pointer(24))
  rrt(025) = bolsig_rates(bolsig_pointer(25))
  rrt(026) = bolsig_rates(bolsig_pointer(26))
  rrt(027) = bolsig_rates(bolsig_pointer(27))
  rrt(028) = KVT10_N2N2
  rrt(029) = KVT10_N2N2*2.0D0
  rrt(030) = KVT10_N2N2*3.0D0
  rrt(031) = KVT10_N2N2*4.0D0
  rrt(032) = KVT10_N2N2*5.0D0
  rrt(033) = KVT10_N2N2*6.0D0
  rrt(034) = KVT10_N2N2*7.0D0
  rrt(035) = KVT10_N2N2*8.0D0
  rrt(036) = KVT01_N2N2
  rrt(037) = KVT01_N2N2*2.0D0
  rrt(038) = KVT01_N2N2*3.0D0
  rrt(039) = KVT01_N2N2*4.0D0
  rrt(040) = KVT01_N2N2*5.0D0
  rrt(041) = KVT01_N2N2*6.0D0
  rrt(042) = KVT01_N2N2*7.0D0
  rrt(043) = KVT01_N2N2*8.0D0
  rrt(044) = KVT10_N2N
  rrt(045) = KVT10_N2N*2.0D0
  rrt(046) = KVT10_N2N*3.0D0
  rrt(047) = KVT10_N2N*4.0D0
  rrt(048) = KVT10_N2N*5.0D0
  rrt(049) = KVT10_N2N*6.0D0
  rrt(050) = KVT10_N2N*7.0D0
  rrt(051) = KVT10_N2N*8.0D0
  rrt(052) = KVT01_N2N
  rrt(053) = KVT01_N2N*2.0D0
  rrt(054) = KVT01_N2N*3.0D0
  rrt(055) = KVT01_N2N*4.0D0
  rrt(056) = KVT01_N2N*5.0D0
  rrt(057) = KVT01_N2N*6.0D0
  rrt(058) = KVT01_N2N*7.0D0
  rrt(059) = KVT01_N2N*8.0D0
  rrt(060) = KVT10_N2O
  rrt(061) = KVT10_N2O*2.0D0
  rrt(062) = KVT10_N2O*3.0D0
  rrt(063) = KVT10_N2O*4.0D0
  rrt(064) = KVT10_N2O*5.0D0
  rrt(065) = KVT10_N2O*6.0D0
  rrt(066) = KVT10_N2O*7.0D0
  rrt(067) = KVT10_N2O*8.0D0
  rrt(068) = KVT01_N2O
  rrt(069) = KVT01_N2O*2.0D0
  rrt(070) = KVT01_N2O*3.0D0
  rrt(071) = KVT01_N2O*4.0D0
  rrt(072) = KVT01_N2O*5.0D0
  rrt(073) = KVT01_N2O*6.0D0
  rrt(074) = KVT01_N2O*7.0D0
  rrt(075) = KVT01_N2O*8.0D0
  rrt(076) = KVT10_O2O2
  rrt(077) = KVT10_O2O2*2.0D0
  rrt(078) = KVT10_O2O2*3.0D0
  rrt(079) = KVT10_O2O2*4.0D0
  rrt(080) = KVT01_O2O2
  rrt(081) = KVT01_O2O2*2.0D0
  rrt(082) = KVT01_O2O2*3.0D0
  rrt(083) = KVT01_O2O2*4.0D0
  rrt(084) = KVT10_O2O
  rrt(085) = KVT10_O2O*2.0D0
  rrt(086) = KVT10_O2O*3.0D0
  rrt(087) = KVT10_O2O*4.0D0
  rrt(088) = KVT01_O2O
  rrt(089) = KVT01_O2O*2.0D0
  rrt(090) = KVT01_O2O*3.0D0
  rrt(091) = KVT01_O2O*4.0D0
  rrt(092) = bolsig_rates(bolsig_pointer(28))
  rrt(093) = bolsig_rates(bolsig_pointer(29))
  rrt(094) = bolsig_rates(bolsig_pointer(30))
  rrt(095) = bolsig_rates(bolsig_pointer(31))
  rrt(096) = bolsig_rates(bolsig_pointer(32))
  rrt(097) = bolsig_rates(bolsig_pointer(33))
  rrt(098) = bolsig_rates(bolsig_pointer(34))
  rrt(099) = bolsig_rates(bolsig_pointer(35))
  rrt(100) = bolsig_rates(bolsig_pointer(36))
  rrt(101) = bolsig_rates(bolsig_pointer(37))
  rrt(102) = bolsig_rates(bolsig_pointer(38))
  rrt(103) = bolsig_rates(bolsig_pointer(39))
  rrt(104) = bolsig_rates(bolsig_pointer(40))
  rrt(105) = bolsig_rates(bolsig_pointer(41))
  rrt(106) = bolsig_rates(bolsig_pointer(42))
  rrt(107) = bolsig_rates(bolsig_pointer(43))
  rrt(108) = bolsig_rates(bolsig_pointer(44))
  rrt(109) = bolsig_rates(bolsig_pointer(45))
  rrt(110) = bolsig_rates(bolsig_pointer(46))
  rrt(111) = bolsig_rates(bolsig_pointer(47))
  rrt(112) = bolsig_rates(bolsig_pointer(48))
  rrt(113) = bolsig_rates(bolsig_pointer(49))
  rrt(114) = bolsig_rates(bolsig_pointer(50))
  rrt(115) = bolsig_rates(bolsig_pointer(51))
  rrt(116) = bolsig_rates(bolsig_pointer(52))
  rrt(117) = bolsig_rates(bolsig_pointer(53))
  rrt(118) = bolsig_rates(bolsig_pointer(54))
  rrt(119) = bolsig_rates(bolsig_pointer(55))
  rrt(120) = bolsig_rates(bolsig_pointer(56))
  rrt(121) = bolsig_rates(bolsig_pointer(57))
  rrt(122) = 1.8D-7*(300.0D0/TE)**0.39*0.50D0
  rrt(123) = 1.8D-7*(300.0D0/TE)**0.39*0.45D0
  rrt(124) = 1.8D-7*(300.0D0/TE)**0.39*0.05D0
  rrt(125) = 2.0D-7*(300.0D0/TE)**0.5
  rrt(126) = 2.3D-6*(300.0D0/TE)**0.53
  rrt(127) = 2.7D-7*(300.0D0/TE)**0.7*0.55D0
  rrt(128) = 2.7D-7*(300.0D0/TE)**0.7*0.40D0
  rrt(129) = 2.7D-7*(300.0D0/TE)**0.7*0.05D0
  rrt(130) = 1.4D-6*(300.0D0/TE)**0.5
  rrt(131) = 4.2D-7*(300.0D0/TE)**0.85*0.20D0
  rrt(132) = 4.2D-7*(300.0D0/TE)**0.85*0.80D0
  rrt(133) = 1.3D-6*(300.0D0/TE)**0.5
  rrt(134) = 7.0D-20*(300.0D0/TE)**4.5
  rrt(135) = rrt(134)
  rrt(136) = 6.0D-27*(300.0D0/TE)**1.5*ANY_NEUTRAL
  rrt(137) = rrt(136)
  rrt(138) = 1.0D-9
  rrt(139) = 1.0D-11
  rrt(140) = 1.0D-31
  rrt(141) = 1.0D-31
  rrt(142) = 1.0D-31
  rrt(143) = 1.1D-31*(300.0D0/TE)**2*EXP(-70.0D0/TGAS)*EXP(1500.0D0*(TE-TGAS)/(TE*TGAS))
  rrt(144) = 1.0D-30*ANY_NEUTRAL
  rrt(145) = 1.4D-10
  rrt(146) = 2.6D-10
  rrt(147) = 5.0D-15
  rrt(148) = 3.0D-10
  rrt(149) = 6.9D-10
  rrt(150) = 2.2D-9
  rrt(151) = 1.9D-9
  rrt(152) = 3.0D-10
  rrt(153) = 1.5D-10
  rrt(154) = 2.7D-10*(TEFFN2/300.0D0)**0.5*EXP(-5590.0D0/TEFFN2)
  rrt(155) = 2.0D-10
  rrt(156) = 3.6D-10
  rrt(157) = 1.9D-12*(TEFFN2/300.0D0)**0.5*EXP(-4990.0D0/TEFFN2)
  rrt(158) = 2.1D-9
  rrt(159) = 2.5D-9
  rrt(160) = 3.0D-10
  rrt(161) = 0.50D0
  rrt(162) = 1.34D5
  rrt(163) = 1.0D2
  rrt(164) = 2.45D7
  rrt(165) = 2.6D-4
  rrt(166) = 1.5D-3
  rrt(167) = 8.5D-2
  rrt(168) = 11.0D0
  rrt(169) = 7.0D-12
  rrt(170) = 2.1D-11
  rrt(171) = 2.0D-12
  rrt(172) = 4.0D-11*(300.0D0/TGAS)**0.667
  rrt(173) = 2.1D-12*(TGAS/300.0D0)**0.55
  rrt(174) = 2.0D-13*(TGAS/300.0D0)**0.55
  rrt(175) = rrt(174)
  rrt(176) = 3.0D-16
  rrt(177) = 6.9D-11
  rrt(178) = 3.0D-10
  rrt(179) = 1.5D-10
  rrt(180) = 3.0D-11
  rrt(181) = 2.0D-12
  rrt(182) = 3.0D-10
  rrt(183) = 2.4D-10
  rrt(184) = 1.0D-11
  rrt(185) = 3.0D-10
  rrt(186) = 1.9D-13
  rrt(187) = 2.8D-11
  rrt(188) = 3.6D-10
  rrt(189) = 3.2D-12
  rrt(190) = 1.0D-11
  rrt(191) = 1.7D-33*ANY_NEUTRAL
  rrt(192) = 1.0D-32-1.7D-33
  rrt(193) = 1.0D-32-1.7D-33
  rrt(194) = 2.4D-33*ANY_NEUTRAL
  rrt(195) = 1.4D-32-2.4D-33
  rrt(196) = 1.4D-32-2.4D-33
  rrt(197) = 4.0D-13
  rrt(198) = 5.2D-12
  rrt(199) = 1.8D-10
  rrt(200) = 1.0D-13*EXP(-510.0D0/TGAS)
  rrt(201) = 1.8D-12
  rrt(202) = 1.0D-12
  rrt(203) = 6.0D-13
  rrt(204) = 6.0D-14
  rrt(205) = 1.0D-13
  rrt(206) = 2.6D-12
  rrt(207) = 3.0D-11
  rrt(208) = 7.0D-16
  rrt(209) = 2.0D-14*EXP(-600.0D0/TGAS)
  rrt(210) = 3.8D-18*EXP(-205.0D0/TGAS)
  rrt(211) = 3.0D-21
  rrt(212) = 2.5D-11
  rrt(213) = 5.2D-11*EXP(-2840.0D0/TGAS)
  rrt(214) = 7.0D-28*TGAS**3.8*EXP(700.0D0/TGAS)
  rrt(215) = 1.0D-11*EXP(-2300.0D0/TGAS)
  rrt(216) = 8.1D-14
  rrt(217) = 3.4D-11*(300.0D0/TGAS)**0.1*EXP(-4200.0D0/TGAS)
  rrt(218) = 4.3D-22*TGAS**2.4*EXP(-281.0D0/TGAS)
  rrt(219) = 1.7D-15*(TGAS/300.0D0)
  rrt(220) = 6.0D-14
  rrt(221) = 2.2D-11
  rrt(222) = 9.0D-12
  rrt(223) = 3.0D-13
  rrt(224) = 9.0D-15
  rrt(225) = 0.07*6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*ANY_NEUTRAL
  rrt(226) = 0.07*6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*(5.9D0-1.0D0)
  rrt(227) = 0.07*6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*(21.D0-1.0D0)
  rrt(228) = 0.01*6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*ANY_NEUTRAL
  rrt(229) = 1.0D-31
  rrt(230) = 8.0D-12
  rrt(231) = 6.4D-12*EXP(67.0D0/TGAS)
  rrt(232) = 1.0D-12
  rrt(233) = 2.6D-11*EXP(67.0D0/TGAS)
  rrt(234) = 2.3D-11
  rrt(235) = 1.2D-10
  rrt(236) = 1.2D-10
  rrt(237) = 1.7D-10
  rrt(238) = 5.0D-11*EXP(-300.0D0/TGAS)
  rrt(239) = 1.0D-12
  rrt(240) = 1.3D-12*EXP(-850.0D0/TGAS)
  rrt(241) = 3.0D-12*EXP(-850.0D0/TGAS)
  rrt(242) = 1.0D-17
  rrt(243) = 1.1D-10
  rrt(244) = 2.9D-11
  rrt(245) = 3.2D-11
  rrt(246) = 2.9D-10
  rrt(247) = 5.1D-10
  rrt(248) = 2.9D-10
  rrt(249) = 2.9D-10
  rrt(250) = 1.8D-11*(TGAS/300.0)**0.5
  rrt(251) = 3.2D-12*(TGAS/300.0)*EXP(-3150.0D0/TGAS)
  rrt(252) = 2.0D-16
  rrt(253) = 3.0D-10*EXP(-38370.D0/TGAS)
  rrt(254) = 7.5D-12*(TGAS/300.0)*EXP(-19500.0D0/TGAS)
  rrt(255) = 2.0D-11*EXP(-2280.0D0/TGAS)
  rrt(256) = 5.1D-13*EXP(-33660.0D0/TGAS)
  rrt(257) = 2.0D-11*EXP(-49800.0D0/TGAS)
  rrt(258) = 2.7D-11*EXP(-6.74D4/TGAS)
  rrt(259) = 1.6D-12*(TGAS/300.0D0)**0.5*(0.19D0+8.6D0*TGAS)*EXP(-32000.0D0/TGAS)
  rrt(260) = 5.4D-8*(1.0D0-EXP(-3354.0D0/TGAS))*EXP(-113200.0D0/TGAS)*ANY_NEUTRAL
  rrt(261) = 5.4D-8*(1.0D0-EXP(-3354.0D0/TGAS))*EXP(-113200.0D0/TGAS)*(6.6D0-1.0D0)
  rrt(262) = rrt(261)
  rrt(263) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*ANY_NEUTRAL
  rrt(264) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*(5.9D0-1.0D0)
  rrt(265) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*(21.D0-1.0D0)
  rrt(266) = 8.7D-9*EXP(-75994.0D0/TGAS)*ANY_NEUTRAL
  rrt(267) = 8.7D-9*EXP(-75994.0D0/TGAS)*(20.0D0-1.0D0)
  rrt(268) = rrt(267)
  rrt(269) = 6.6D-10*EXP(-11600.0D0/TGAS)
  rrt(270) = 6.6D-10*EXP(-11600.0D0/TGAS)*0.38D0
  rrt(271) = 6.6D-10*EXP(-11600.0D0/TGAS)*6.3D0*EXP(170.0D0/TGAS)
  rrt(272) = rrt(271)
  rrt(273) = MAX(8.3D-34*EXP(500.0D0/TGAS),1.91D-33)
  rrt(274) = 1.8D-33*EXP(435.0D0/TGAS)
  rrt(275) = rrt(274)
  rrt(276) = 1.8D-33*EXP(435.0D0/TGAS)*3.0D0
  rrt(277) = rrt(276)
  rrt(278) = MAX(2.8D-34*EXP(720.0D0/TGAS),1.0D-33*(300.0D0/TGAS)**0.41)
  rrt(279) = 4.0D-33*(300.0D0/TGAS)**0.41
  rrt(280) = 4.0D-33*(300.0D0/TGAS)**0.41*0.17D0
  rrt(281) = 4.0D-33*(300.0D0/TGAS)**0.41*0.8D0
  rrt(282) = 4.0D-33*(300.0D0/TGAS)**0.41*3.6D0
  rrt(283) = 1.0D-32*(300.0D0/TGAS)**0.5
  rrt(284) = rrt(283)
  rrt(285) = 1.8D-31*(300.0D0/TGAS)
  rrt(286) = rrt(285)
  rrt(287) = rrt(285)
  rrt(288) = MAX(5.8D-34*(300.0D0/TGAS)**2.8,5.4D-34*(300.0D0/TGAS)**1.9)
  rrt(289) = 7.6D-34*(300.0D0/TGAS)**1.9
  rrt(290) = rrt(289)
  rrt(291) = MIN(3.9D-33*(300.0D0/TGAS)**1.9,1.1D-34*EXP(1060.0D0/TGAS))
  rrt(292) = rrt(291)
  rrt(293) = 1.0D-12
  rrt(294) = 2.8D-10
  rrt(295) = 2.5D-10
  rrt(296) = 2.8D-11
  rrt(297) = 5.0D-10
  rrt(298) = 8.0D-10
  rrt(299) = 3.0D-12
  rrt(300) = 1.0D-12
  rrt(301) = (1.5D0-2.0D-3*TEFFN+9.6D-7*TEFFN**2)*1.0D-12
  rrt(302) = 2.0D-11*(300.0D0/TEFFN)**0.5
  rrt(303) = 1.0D-10
  rrt(304) = 2.4D-11
  rrt(305) = 3.0D-12
  rrt(306) = 1.3D-10
  rrt(307) = 6.0D-11*(300.0D0/TEFFN2)**0.5
  rrt(308) = 1.3D-10*(300.0D0/TEFFN2)**0.5
  rrt(309) = 1.0D-11*(300.0D0/TEFFN2)**0.2
  rrt(310) = 1.0D-10
  rrt(311) = 7.2D-13*(TEFFN2/300.0D0)
  rrt(312) = 3.3D-10
  rrt(313) = 1.0D-17
  rrt(314) = 1.2D-10
  rrt(315) = 6.3D-10
  rrt(316) = 2.3D-11
  rrt(317) = 6.6D-11
  rrt(318) = 7.0D-11
  rrt(319) = MIN(2.1D-16*EXP(TEFFN4/121.0D0),1.0D-10)
  rrt(320) = 2.5D-10
  rrt(321) = 2.5D-10
  rrt(322) = 1.0D-11
  rrt(323) = 4.0D-10
  rrt(324) = 4.6D-12*(TEFFN4/300.0D0)**2.5*EXP(-2650.0D0/TEFFN4)
  rrt(325) = 3.3D-6*(300.0D0/TEFFN4)**4*EXP(-5030.0D0/TEFFN4)
  rrt(326) = 1.0D-10
  rrt(327) = 1.0D-10
  rrt(328) = 3.0D-10
  rrt(329) = 1.0D-10
  rrt(330) = 1.1D-6*(300.0D0/TEFFN4)**5.3*EXP(-2360.0D0/TEFFN4)
  rrt(331) = 1.0D-9
  rrt(332) = 1.7D-29*(300.0D0/TEFFN)**2.1
  rrt(333) = 1.0D-29*ANY_NEUTRAL
  rrt(334) = rrt(333)
  rrt(335) = 6.0D-29*(300.0D0/TEFFN)**2*ANY_NEUTRAL
  rrt(336) = rrt(333)
  rrt(337) = rrt(333)
  rrt(338) = 5.2D-29*(300.0D0/TEFFN2)**2.2
  rrt(339) = 9.0D-30*EXP(400.0D0/TEFFN2)
  rrt(340) = 2.4D-30*(300.0D0/TEFFN2)**3.2
  rrt(341) = 9.0D-31*(300.0D0/TEFFN2)**2
  rrt(342) = 1.0D-10
  rrt(343) = 8.0D-10
  rrt(344) = 3.3D-10
  rrt(345) = 3.5D-10
  rrt(346) = 3.2D-10
  rrt(347) = 5.0D-10
  rrt(348) = 1.0D-10*EXP(-1044.0D0/TEFFN4)*ANY_NEUTRAL
  rrt(349) = 4.0D-10
  rrt(350) = 3.0D-10
  rrt(351) = 1.0D-10
  rrt(352) = 1.0D-10
  rrt(353) = 1.1D-30*(300.0D0/TEFFN)*ANY_NEUTRAL
  rrt(354) = 3.5D-31*(300.0D0/TEFFN2)*ANY_NEUTRAL
  rrt(355) = 2.0D-7*(300.0D0/TIONN)
  rrt(356) = rrt(355)
  rrt(357) = rrt(355)
  rrt(358) = rrt(355)
  rrt(359) = rrt(355)
  rrt(360) = 2.0D-7*(300.0D0/TIONN2)
  rrt(361) = rrt(360)
  rrt(362) = rrt(360)
  rrt(363) = rrt(360)
  rrt(364) = rrt(360)
  rrt(365) = 2.0D-7*(300.0D0/TIONN3)
  rrt(366) = rrt(365)
  rrt(367) = rrt(365)
  rrt(368) = rrt(365)
  rrt(369) = rrt(365)
  rrt(370) = rrt(360)
  rrt(371) = rrt(360)
  rrt(372) = rrt(360)
  rrt(373) = rrt(360)
  rrt(374) = rrt(360)
  rrt(375) = 1.0D-7
  rrt(376) = 1.0D-7
  rrt(377) = 1.0D-7
  rrt(378) = 1.0D-7
  rrt(379) = 1.0D-7
  rrt(380) = 1.0D-7
  rrt(381) = 1.0D-7
  rrt(382) = 1.0D-7
  rrt(383) = 1.0D-7
  rrt(384) = 1.0D-7
  rrt(385) = 1.0D-7
  rrt(386) = 1.0D-7
  rrt(387) = 1.0D-7
  rrt(388) = 1.0D-7
  rrt(389) = 1.0D-7
  rrt(390) = 1.0D-7
  rrt(391) = 1.0D-7
  rrt(392) = 1.0D-7
  rrt(393) = 1.0D-7
  rrt(394) = 1.0D-7
  rrt(395) = 1.0D-7
  rrt(396) = 1.0D-7
  rrt(397) = 1.0D-7
  rrt(398) = 1.0D-7
  rrt(399) = 1.0D-7
  rrt(400) = 1.0D-7
  rrt(401) = 1.0D-7
  rrt(402) = 1.0D-7
  rrt(403) = 1.0D-7
  rrt(404) = 1.0D-7
  rrt(405) = 1.0D-7
  rrt(406) = 1.0D-7
  rrt(407) = 1.0D-7
  rrt(408) = 1.0D-7
  rrt(409) = 1.0D-7
  rrt(410) = 1.0D-7
  rrt(411) = 1.0D-7
  rrt(412) = 2.0D-25*(300.0D0/TIONN)**2.5*ANY_NEUTRAL
  rrt(413) = rrt(412)
  rrt(414) = rrt(412)
  rrt(415) = rrt(412)
  rrt(416) = rrt(412)
  rrt(417) = 2.0D-25*(300.0D0/TIONN2)**2.5*ANY_NEUTRAL
  rrt(418) = rrt(417)
  rrt(419) = rrt(417)
  rrt(420) = rrt(417)
  rrt(421) = rrt(417)
  rrt(422) = rrt(412)
  rrt(423) = rrt(412)
  rrt(424) = rrt(412)
  rrt(425) = rrt(412)
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
