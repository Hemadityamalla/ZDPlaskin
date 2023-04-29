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
! Sun May  2 18:06:22 2021
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
  integer, parameter :: species_max = 53, species_electrons = 53, species_length = 9, reactions_max = 625, reactions_length = 44
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
  integer, parameter, private               :: bolsig_species_max = 2, bolsig_species_length = 2, bolsig_rates_max = 33 
  character(*), parameter, private          :: bolsigfile = "bolsigdb.dat"
  integer                                   :: bolsig_pointer(bolsig_rates_max) = -1
  integer, private                          :: bolsig_species_index(bolsig_species_max) = -1, bolsig_collisions_max = 0 
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
  / 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0,&
    0, 0, 1, 1, 1,-1,-1,-1,-1, 1,-1/
  data species_name(1:species_max) &
  /"N2       ","N2(V1)   ","N2(V2)   ","N2(V3)   ","N2(V4)   ","N2(V5)   ","N2(V6)   ","N2(V7)   ","N2(V8)   ","N2(A3)   ",&
   "N2(B3)   ","N2(A`1)  ","N2(C3)   ","N        ","N(2D)    ","N(2P)    ","N^+      ","N2^+     ","N3^+     ","N4^+     ",&
   "O2       ","O2(V1)   ","O2(V2)   ","O2(V3)   ","O2(V4)   ","O2(A1)   ","O2(B1)   ","O2(4.5EV)","O        ","O(1D)    ",&
   "O(1S)    ","O3       ","O^+      ","O2^+     ","O4^+     ","O^-      ","O2^-     ","O3^-     ","O4^-     ","NO       ",&
   "N2O      ","NO2      ","NO3      ","N2O5     ","NO^+     ","N2O^+    ","NO2^+    ","NO^-     ","N2O^-    ","NO2^-    ",&
   "NO3^-    ","O2^+N2   ","E        "/
  data reaction_sign(1:36) &
  /"N2(V1)+N2=>N2+N2                            ","N2(V2)+N2=>N2(V1)+N2                        ",&
   "N2(V3)+N2=>N2(V2)+N2                        ","N2(V4)+N2=>N2(V3)+N2                        ",&
   "N2(V5)+N2=>N2(V4)+N2                        ","N2(V6)+N2=>N2(V5)+N2                        ",&
   "N2(V7)+N2=>N2(V6)+N2                        ","N2(V8)+N2=>N2(V7)+N2                        ",&
   "N2+N2=>N2(V1)+N2                            ","N2(V1)+N2=>N2(V2)+N2                        ",&
   "N2(V2)+N2=>N2(V3)+N2                        ","N2(V3)+N2=>N2(V4)+N2                        ",&
   "N2(V4)+N2=>N2(V5)+N2                        ","N2(V5)+N2=>N2(V6)+N2                        ",&
   "N2(V6)+N2=>N2(V7)+N2                        ","N2(V7)+N2=>N2(V8)+N2                        ",&
   "N2(V1)+N=>N2+N                              ","N2(V2)+N=>N2(V1)+N                          ",&
   "N2(V3)+N=>N2(V2)+N                          ","N2(V4)+N=>N2(V3)+N                          ",&
   "N2(V5)+N=>N2(V4)+N                          ","N2(V6)+N=>N2(V5)+N                          ",&
   "N2(V7)+N=>N2(V6)+N                          ","N2(V8)+N=>N2(V7)+N                          ",&
   "N2+N=>N2(V1)+N                              ","N2(V1)+N=>N2(V2)+N                          ",&
   "N2(V2)+N=>N2(V3)+N                          ","N2(V3)+N=>N2(V4)+N                          ",&
   "N2(V4)+N=>N2(V5)+N                          ","N2(V5)+N=>N2(V6)+N                          ",&
   "N2(V6)+N=>N2(V7)+N                          ","N2(V7)+N=>N2(V8)+N                          ",&
   "N2(V1)+O=>N2+O                              ","N2(V2)+O=>N2(V1)+O                          ",&
   "N2(V3)+O=>N2(V2)+O                          ","N2(V4)+O=>N2(V3)+O                          "/
  data reaction_sign(37:72) &
  /"N2(V5)+O=>N2(V4)+O                          ","N2(V6)+O=>N2(V5)+O                          ",&
   "N2(V7)+O=>N2(V6)+O                          ","N2(V8)+O=>N2(V7)+O                          ",&
   "N2+O=>N2(V1)+O                              ","N2(V1)+O=>N2(V2)+O                          ",&
   "N2(V2)+O=>N2(V3)+O                          ","N2(V3)+O=>N2(V4)+O                          ",&
   "N2(V4)+O=>N2(V5)+O                          ","N2(V5)+O=>N2(V6)+O                          ",&
   "N2(V6)+O=>N2(V7)+O                          ","N2(V7)+O=>N2(V8)+O                          ",&
   "O2(V1)+O2=>O2+O2                            ","O2(V2)+O2=>O2(V1)+O2                        ",&
   "O2(V3)+O2=>O2(V2)+O2                        ","O2(V4)+O2=>O2(V3)+O2                        ",&
   "O2+O2=>O2(V1)+O2                            ","O2(V1)+O2=>O2(V2)+O2                        ",&
   "O2(V2)+O2=>O2(V3)+O2                        ","O2(V3)+O2=>O2(V4)+O2                        ",&
   "O2(V1)+O=>O2+O                              ","O2(V2)+O=>O2(V1)+O                          ",&
   "O2(V3)+O=>O2(V2)+O                          ","O2(V4)+O=>O2(V3)+O                          ",&
   "O2+O=>O2(V1)+O                              ","O2(V1)+O=>O2(V2)+O                          ",&
   "O2(V2)+O=>O2(V3)+O                          ","O2(V3)+O=>O2(V4)+O                          ",&
   "bolsig:N2->N2(A3SU)                         ","bolsig:N2->N2(A'1)                          ",&
   "bolsig:N2->N2(A1)                           ","bolsig:N2->N2(W1)                           ",&
   "bolsig:N2->N2(C3)                           ","bolsig:N2->N2(E3)                           ",&
   "bolsig:N2->N2(A''1)                         ","bolsig:N2->N2(SUM)                          "/
  data reaction_sign(73:108) &
  /"bolsig:O2->O2(A1)                           ","bolsig:O2->O2(B1)                           ",&
   "bolsig:O2->O2(4.5EV)                        ","bolsig:O2->O2(6.0EV)                        ",&
   "bolsig:O2->O2(8.4EV)                        ","bolsig:O2->O2(9.97EV)                       ",&
   "bolsig:O2(A1)->O+O                          ","bolsig:O->O(1D)                             ",&
   "bolsig:O->O(1S)                             ","bolsig:N2(A3)->N2                           ",&
   "bolsig:O2(A1)->O2                           ","bolsig:N->N^+                               ",&
   "bolsig:O->O^+                               ","bolsig:N2->N2^+                             ",&
   "bolsig:N2(A3)->N2^+                         ","bolsig:O2->O2^+                             ",&
   "bolsig:O2(A1)->O2^+                         ","bolsig:NO->NO^+                             ",&
   "bolsig:N2O->N2O^+                           ","E+N2^+=>N+N                                 ",&
   "E+N2^+=>N+N(2D)                             ","E+N2^+=>N+N(2P)                             ",&
   "E+O2^+=>O+O                                 ","E+O2^+=>O+O(1D)                             ",&
   "E+O2^+=>O+O(1S)                             ","E+NO^+=>O+N                                 ",&
   "E+NO^+=>O+N(2D)                             ","E+N3^+=>N2+N                                ",&
   "E+N4^+=>N2+N2                               ","E+N2O^+=>N2+O                               ",&
   "E+NO2^+=>NO+O                               ","E+O4^+=>O2+O2                               ",&
   "E+O2^+N2=>O2+N2                             ","E+N^++E=>N+E                                ",&
   "E+O^++E=>O+E                                ","E+N^++ANY_NEUTRAL=>N+ANY_NEUTRAL            "/
  data reaction_sign(109:144) &
  /"E+O^++ANY_NEUTRAL=>O+ANY_NEUTRAL            ","bolsig:O2->O^-+O                            ",&
   "bolsig:NO->O^-+N                            ","bolsig:O3->O^-+O2                           ",&
   "bolsig:O3->O2^-+O                           ","bolsig:N2O->NO^-+N                          ",&
   "bolsig:O2+O2->O2^-+O2                       ","E+NO2=>O^-+NO                               ",&
   "E+O+O2=>O^-+O2                              ","E+O+O2=>O2^-+O                              ",&
   "E+O3+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL          ","E+NO+ANY_NEUTRAL=>NO^-+ANY_NEUTRAL          ",&
   "E+N2O+ANY_NEUTRAL=>N2O^-+ANY_NEUTRAL        ","E+O2+N2=>O2^-+N2                            ",&
   "O^-+O=>O2+E                                 ","O^-+N=>NO+E                                 ",&
   "O^-+NO=>NO2+E                               ","O^-+N2=>N2O+E                               ",&
   "O^-+O2=>O3+E                                ","O^-+O2(A1)=>O3+E                            ",&
   "O^-+O2(B1)=>O+O2+E                          ","O^-+N2(A3)=>O+N2+E                          ",&
   "O^-+N2(B3)=>O+N2+E                          ","O^-+O3=>O2+O2+E                             ",&
   "O2^-+O=>O3+E                                ","O2^-+N=>NO2+E                               ",&
   "O2^-+O2=>O2+O2+E                            ","O2^-+O2(A1)=>O2+O2+E                        ",&
   "O2^-+O2(B1)=>O2+O2+E                        ","O2^-+N2=>O2+N2+E                            ",&
   "O2^-+N2(A3)=>O2+N2+E                        ","O2^-+N2(B3)=>O2+N2+E                        ",&
   "O3^-+O=>O2+O2+E                             ","NO^-+N=>N2O+E                               ",&
   "O3^-+N=>NO+O2+E                             ","N2O^-+N=>NO+N2+E                            "/
  data reaction_sign(145:180) &
  /"NO2^-+N=>NO+NO+E                            ","NO3^-+N=>NO+NO2+E                           ",&
   "NO^-+O=>NO2+E                               ","N2O^-+O=>NO+NO+E                            ",&
   "NO2^-+O=>NO+O2+E                            ","NO3^-+O=>NO+O3+E                            ",&
   "O3^-+N2(A3)=>O3+N2+E                        ","NO^-+N2(A3)=>NO+N2+E                        ",&
   "N2O^-+N2(A3)=>N2O+N2+E                      ","NO2^-+N2(A3)=>NO2+N2+E                      ",&
   "NO3^-+N2(A3)=>NO3+N2+E                      ","O3^-+N2(B3)=>O3+N2+E                        ",&
   "NO^-+N2(B3)=>NO+N2+E                        ","N2O^-+N2(B3)=>N2O+N2+E                      ",&
   "NO2^-+N2(B3)=>NO2+N2+E                      ","NO3^-+N2(B3)=>NO3+N2+E                      ",&
   "N2(A3)=>N2                                  ","N2(B3)=>N2(A3)                              ",&
   "N2(A`1)=>N2                                 ","N2(C3)=>N2(B3)                              ",&
   "O2(A1)=>O2                                  ","O2(B1)=>O2(A1)                              ",&
   "O2(B1)=>O2                                  ","O2(4.5EV)=>O2                               ",&
   "N2(A3)+O=>NO+N(2D)                          ","N2(A3)+O=>N2+O(1S)                          ",&
   "N2(A3)+N=>N2+N                              ","N2(A3)+N=>N2+N(2P)                          ",&
   "N2(A3)+O2=>N2+O+O(1D)                       ","N2(A3)+O2=>N2+O2(A1)                        ",&
   "N2(A3)+O2=>N2+O2(B1)                        ","N2(A3)+O2=>N2O+O                            ",&
   "N2(A3)+N2=>N2+N2                            ","N2(A3)+NO=>N2+NO                            ",&
   "N2(A3)+N2O=>N2+N+NO                         ","N2(A3)+NO2=>N2+O+NO                         "/
  data reaction_sign(181:216) &
  /"N2(A3)+N2(A3)=>N2+N2(B3)                    ","N2(A3)+N2(A3)=>N2+N2(C3)                    ",&
   "N2(B3)+N2=>N2(A3)+N2                        ","N2(B3)+N2=>N2+N2                            ",&
   "N2(B3)+O2=>N2+O+O                           ","N2(B3)+NO=>N2(A3)+NO                        ",&
   "N2(C3)+N2=>N2(A`1)+N2                       ","N2(C3)+O2=>N2+O+O(1S)                       ",&
   "N2(A`1)+N2=>N2(B3)+N2                       ","N2(A`1)+O2=>N2+O+O                          ",&
   "N2(A`1)+NO=>N2+N+O                          ","N2(A`1)+N2(A3)=>N4^++E                      ",&
   "N2(A`1)+N2(A`1)=>N4^++E                     ","N+N+N2=>N2(A3)+N2                           ",&
   "N+N+O2=>N2(A3)+O2                           ","N+N+NO=>N2(A3)+NO                           ",&
   "N+N+N=>N2(A3)+N                             ","N+N+O=>N2(A3)+O                             ",&
   "N+N+N2=>N2(B3)+N2                           ","N+N+O2=>N2(B3)+O2                           ",&
   "N+N+NO=>N2(B3)+NO                           ","N+N+N=>N2(B3)+N                             ",&
   "N+N+O=>N2(B3)+O                             ","N(2D)+O=>N+O(1D)                            ",&
   "N(2D)+O2=>NO+O                              ","N(2D)+NO=>N2+O                              ",&
   "N(2D)+N2O=>NO+N2                            ","N(2D)+N2=>N+N2                              ",&
   "N(2P)+N=>N+N                                ","N(2P)+O=>N+O                                ",&
   "N(2P)+N=>N(2D)+N                            ","N(2P)+N2=>N+N2                              ",&
   "N(2P)+N(2D)=>N2^++E                         ","N(2P)+O2=>NO+O                              ",&
   "N(2P)+NO=>N2(A3)+O                          ","O2(A1)+O=>O2+O                              "/
  data reaction_sign(217:252) &
  /"O2(A1)+N=>NO+O                              ","O2(A1)+O2=>O2+O2                            ",&
   "O2(A1)+N2=>O2+N2                            ","O2(A1)+NO=>O2+NO                            ",&
   "O2(A1)+O3=>O2+O2+O(1D)                      ","O2(A1)+O2(A1)=>O2+O2(B1)                    ",&
   "O+O3=>O2+O2(A1)                             ","O2(B1)+O=>O2(A1)+O                          ",&
   "O2(B1)+O=>O2+O(1D)                          ","O2(B1)+O2=>O2(A1)+O2                        ",&
   "O2(B1)+N2=>O2(A1)+N2                        ","O2(B1)+NO=>O2(A1)+NO                        ",&
   "O2(B1)+O3=>O2+O2+O                          ","O2(4.5EV)+O=>O2+O(1S)                       ",&
   "O2(4.5EV)+O2=>O2(B1)+O2(B1)                 ","O2(4.5EV)+N2=>O2(B1)+N2                     ",&
   "O(1D)+O=>O+O                                ","O(1D)+O2=>O+O2                              ",&
   "O(1D)+O2=>O+O2(A1)                          ","O(1D)+O2=>O+O2(B1)                          ",&
   "O(1D)+N2=>O+N2                              ","O(1D)+O3=>O2+O+O                            ",&
   "O(1D)+O3=>O2+O2                             ","O(1D)+NO=>O2+N                              ",&
   "O(1D)+N2O=>NO+NO                            ","O(1D)+N2O=>O2+N2                            ",&
   "O(1S)+O=>O(1D)+O                            ","O(1S)+N=>O+N                                ",&
   "O(1S)+O2=>O(1D)+O2                          ","O(1S)+O2=>O+O+O                             ",&
   "O(1S)+N2=>O+N2                              ","O(1S)+O2(A1)=>O+O2(4.5EV)                   ",&
   "O(1S)+O2(A1)=>O(1D)+O2(B1)                  ","O(1S)+O2(A1)=>O+O+O                         ",&
   "O(1S)+NO=>O+NO                              ","O(1S)+NO=>O(1D)+NO                          "/
  data reaction_sign(253:288) &
  /"O(1S)+O3=>O2+O2                             ","O(1S)+O3=>O2+O+O(1D)                        ",&
   "O(1S)+N2O=>O+N2O                            ","O(1S)+N2O=>O(1D)+N2O                        ",&
   "N+NO=>O+N2                                  ","N+O2=>O+NO                                  ",&
   "N+NO2=>O+O+N2                               ","N+NO2=>O+N2O                                ",&
   "N+NO2=>N2+O2                                ","N+NO2=>NO+NO                                ",&
   "O+N2=>N+NO                                  ","O+NO=>N+O2                                  ",&
   "O+NO=>NO2                                   ","O+N2O=>N2+O2                                ",&
   "O+N2O=>NO+NO                                ","O+NO2=>NO+O2                                ",&
   "O+NO3=>O2+NO2                               ","N2+O2=>O+N2O                                ",&
   "NO+NO=>N+NO2                                ","NO+NO=>O+N2O                                ",&
   "NO+NO=>N2+O2                                ","NO+O2=>O+NO2                                ",&
   "NO+O3=>O2+NO2                               ","NO+N2O=>N2+NO2                              ",&
   "NO+NO3=>NO2+NO2                             ","O2+O2=>O+O3                                 ",&
   "O2+NO2=>NO+O3                               ","NO2+NO2=>NO+NO+O2                           ",&
   "NO2+NO2=>NO+NO3                             ","NO2+O3=>O2+NO3                              ",&
   "NO2+NO3=>NO+NO2+O2                          ","NO3+O2=>NO2+O3                              ",&
   "NO3+NO3=>O2+NO2+NO2                         ","N+N=>N2^++E                                 ",&
   "N+O=>NO^++E                                 ","N2+N2=>N+N+N2                               "/
  data reaction_sign(289:324) &
  /"N2+O2=>N+N+O2                               ","N2+NO=>N+N+NO                               ",&
   "N2+O=>N+N+O                                 ","N2+N=>N+N+N                                 ",&
   "O2+N2=>O+O+N2                               ","O2+O2=>O+O+O2                               ",&
   "O2+O=>O+O+O                                 ","O2+N=>O+O+N                                 ",&
   "O2+NO=>O+O+NO                               ","NO+N2=>N+O+N2                               ",&
   "NO+O2=>N+O+O2                               ","NO+O=>N+O+O                                 ",&
   "NO+N=>N+O+N                                 ","NO+NO=>N+O+NO                               ",&
   "O3+N2=>O2+O+N2                              ","O3+O2=>O2+O+O2                              ",&
   "O3+N=>O2+O+N                                ","O3+O=>O2+O+O                                ",&
   "N2O+N2=>N2+O+N2                             ","N2O+O2=>N2+O+O2                             ",&
   "N2O+NO=>N2+O+NO                             ","N2O+N2O=>N2+O+N2O                           ",&
   "NO2+N2=>NO+O+N2                             ","NO2+O2=>NO+O+O2                             ",&
   "NO2+NO=>NO+O+NO                             ","NO2+NO2=>NO+O+NO2                           ",&
   "NO3+N2=>NO2+O+N2                            ","NO3+O2=>NO2+O+O2                            ",&
   "NO3+NO=>NO2+O+NO                            ","NO3+N=>NO2+O+N                              ",&
   "NO3+O=>NO2+O+O                              ","NO3+N2=>NO+O2+N2                            ",&
   "NO3+O2=>NO+O2+O2                            ","NO3+NO=>NO+O2+NO                            ",&
   "NO3+N=>NO+O2+N                              ","NO3+O=>NO+O2+O                              "/
  data reaction_sign(325:360) &
  /"N2O5+ANY_NEUTRAL=>NO2+NO3+ANY_NEUTRAL       ","N+N+N2=>N2+N2                               ",&
   "N+N+O2=>N2+O2                               ","N+N+NO=>N2+NO                               ",&
   "N+N+N=>N2+N                                 ","N+N+O=>N2+O                                 ",&
   "O+O+N2=>O2+N2                               ","O+O+O2=>O2+O2                               ",&
   "O+O+N=>O2+N                                 ","O+O+O=>O2+O                                 ",&
   "O+O+NO=>O2+NO                               ","N+O+N2=>NO+N2                               ",&
   "N+O+O2=>NO+O2                               ","N+O+N=>NO+N                                 ",&
   "N+O+O=>NO+O                                 ","N+O+NO=>NO+NO                               ",&
   "O+O2+N2=>O3+N2                              ","O+O2+O2=>O3+O2                              ",&
   "O+O2+NO=>O3+NO                              ","O+O2+N=>O3+N                                ",&
   "O+O2+O=>O3+O                                ","O+N2+ANY_NEUTRAL=>N2O+ANY_NEUTRAL           ",&
   "O+NO+N2=>NO2+N2                             ","O+NO+O2=>NO2+O2                             ",&
   "O+NO+NO=>NO2+NO                             ","O+NO2+N2=>NO3+N2                            ",&
   "O+NO2+O2=>NO3+O2                            ","O+NO2+N=>NO3+N                              ",&
   "O+NO2+O=>NO3+O                              ","O+NO2+NO=>NO3+NO                            ",&
   "NO2+NO3+ANY_NEUTRAL=>N2O5+ANY_NEUTRAL       ","N^++O=>N+O^+                                ",&
   "N^++O2=>O2^++N                              ","N^++O2=>NO^++O                              ",&
   "N^++O2=>O^++NO                              ","N^++O3=>NO^++O2                             "/
  data reaction_sign(361:396) &
  /"N^++NO=>NO^++N                              ","N^++NO=>N2^++O                              ",&
   "N^++NO=>O^++N2                              ","N^++N2O=>NO^++N2                            ",&
   "O^++N2=>NO^++N                              ","O^++O2=>O2^++O                              ",&
   "O^++O3=>O2^++O2                             ","O^++NO=>NO^++O                              ",&
   "O^++NO=>O2^++N                              ","O^++N(2D)=>N^++O                            ",&
   "O^++N2O=>NO^++NO                            ","O^++N2O=>N2O^++O                            ",&
   "O^++N2O=>O2^++N2                            ","O^++NO2=>NO2^++O                            ",&
   "N2^++O2=>O2^++N2                            ","N2^++O=>NO^++N                              ",&
   "N2^++O3=>O2^++O+N2                          ","N2^++N=>N^++N2                              ",&
   "N2^++NO=>NO^++N2                            ","N2^++N2O=>N2O^++N2                          ",&
   "N2^++N2O=>NO^++N+N2                         ","O2^++N2=>NO^++NO                            ",&
   "O2^++N=>NO^++O                              ","O2^++NO=>NO^++O2                            ",&
   "O2^++NO2=>NO^++O3                           ","O2^++NO2=>NO2^++O2                          ",&
   "N3^++O2=>O2^++N+N2                          ","N3^++O2=>NO2^++N2                           ",&
   "N3^++N=>N2^++N2                             ","N3^++NO=>NO^++N+N2                          ",&
   "N3^++NO=>N2O^++N2                           ","NO2^++NO=>NO^++NO2                          ",&
   "N2O^++NO=>NO^++N2O                          ","N4^++N2=>N2^++N2+N2                         ",&
   "N4^++O2=>O2^++N2+N2                         ","N4^++O=>O^++N2+N2                           "/
  data reaction_sign(397:432) &
  /"N4^++N=>N^++N2+N2                           ","N4^++NO=>NO^++N2+N2                         ",&
   "O4^++N2=>O2^+N2+O2                          ","O4^++O2=>O2^++O2+O2                         ",&
   "O4^++O2(A1)=>O2^++O2+O2                     ","O4^++O2(B1)=>O2^++O2+O2                     ",&
   "O4^++O=>O2^++O3                             ","O4^++NO=>NO^++O2+O2                         ",&
   "O2^+N2+N2=>O2^++N2+N2                       ","O2^+N2+O2=>O4^++N2                          ",&
   "N^++N2+N2=>N3^++N2                          ","N^++O+ANY_NEUTRAL=>NO^++ANY_NEUTRAL         ",&
   "N^++N+ANY_NEUTRAL=>N2^++ANY_NEUTRAL         ","O^++N2+ANY_NEUTRAL=>NO^++N+ANY_NEUTRAL      ",&
   "O^++O+ANY_NEUTRAL=>O2^++ANY_NEUTRAL         ","O^++N+ANY_NEUTRAL=>NO^++ANY_NEUTRAL         ",&
   "N2^++N2+N2=>N4^++N2                         ","N2^++N+N2=>N3^++N2                          ",&
   "O2^++O2+O2=>O4^++O2                         ","O2^++N2+N2=>O2^+N2+N2                       ",&
   "O^-+O2(A1)=>O2^-+O                          ","O^-+O3=>O3^-+O                              ",&
   "O^-+NO2=>NO2^-+O                            ","O^-+N2O=>NO^-+NO                            ",&
   "O^-+N2O=>N2O^-+O                            ","O2^-+O=>O^-+O2                              ",&
   "O2^-+O3=>O3^-+O2                            ","O2^-+NO2=>NO2^-+O2                          ",&
   "O2^-+NO3=>NO3^-+O2                          ","O3^-+O=>O2^-+O2                             ",&
   "O3^-+NO=>NO3^-+O                            ","O3^-+NO=>NO2^-+O2                           ",&
   "O3^-+NO2=>NO2^-+O3                          ","O3^-+NO2=>NO3^-+O2                          ",&
   "O3^-+NO3=>NO3^-+O3                          ","NO^-+O2=>O2^-+NO                            "/
  data reaction_sign(433:468) &
  /"NO^-+NO2=>NO2^-+NO                          ","NO^-+N2O=>NO2^-+N2                          ",&
   "NO2^-+O3=>NO3^-+O2                          ","NO2^-+NO2=>NO3^-+NO                         ",&
   "NO2^-+NO3=>NO3^-+NO2                        ","NO2^-+N2O5=>NO3^-+NO2+NO2                   ",&
   "NO3^-+NO=>NO2^-+NO2                         ","O4^-+N2=>O2^-+O2+N2                         ",&
   "O4^-+O2=>O2^-+O2+O2                         ","O4^-+O=>O3^-+O2                             ",&
   "O4^-+O=>O^-+O2+O2                           ","O4^-+O2(A1)=>O2^-+O2+O2                     ",&
   "O4^-+O2(B1)=>O2^-+O2+O2                     ","O4^-+NO=>NO3^-+O2                           ",&
   "O^-+O2+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL        ","O^-+NO+ANY_NEUTRAL=>NO2^-+ANY_NEUTRAL       ",&
   "O2^-+O2+ANY_NEUTRAL=>O4^-+ANY_NEUTRAL       ","O^-+N^+=>O+N                                ",&
   "O^-+N2^+=>O+N2                              ","O^-+O^+=>O+O                                ",&
   "O^-+O2^+=>O+O2                              ","O^-+NO^+=>O+NO                              ",&
   "O^-+N2O^+=>O+N2O                            ","O^-+NO2^+=>O+NO2                            ",&
   "O2^-+N^+=>O2+N                              ","O2^-+N2^+=>O2+N2                            ",&
   "O2^-+O^+=>O2+O                              ","O2^-+O2^+=>O2+O2                            ",&
   "O2^-+NO^+=>O2+NO                            ","O2^-+N2O^+=>O2+N2O                          ",&
   "O2^-+NO2^+=>O2+NO2                          ","O3^-+N^+=>O3+N                              ",&
   "O3^-+N2^+=>O3+N2                            ","O3^-+O^+=>O3+O                              ",&
   "O3^-+O2^+=>O3+O2                            ","O3^-+NO^+=>O3+NO                            "/
  data reaction_sign(469:504) &
  /"O3^-+N2O^+=>O3+N2O                          ","O3^-+NO2^+=>O3+NO2                          ",&
   "NO^-+N^+=>NO+N                              ","NO^-+N2^+=>NO+N2                            ",&
   "NO^-+O^+=>NO+O                              ","NO^-+O2^+=>NO+O2                            ",&
   "NO^-+NO^+=>NO+NO                            ","NO^-+N2O^+=>NO+N2O                          ",&
   "NO^-+NO2^+=>NO+NO2                          ","N2O^-+N^+=>N2O+N                            ",&
   "N2O^-+N2^+=>N2O+N2                          ","N2O^-+O^+=>N2O+O                            ",&
   "N2O^-+O2^+=>N2O+O2                          ","N2O^-+NO^+=>N2O+NO                          ",&
   "N2O^-+N2O^+=>N2O+N2O                        ","N2O^-+NO2^+=>N2O+NO2                        ",&
   "NO2^-+N^+=>NO2+N                            ","NO2^-+N2^+=>NO2+N2                          ",&
   "NO2^-+O^+=>NO2+O                            ","NO2^-+O2^+=>NO2+O2                          ",&
   "NO2^-+NO^+=>NO2+NO                          ","NO2^-+N2O^+=>NO2+N2O                        ",&
   "NO2^-+NO2^+=>NO2+NO2                        ","NO3^-+N^+=>NO3+N                            ",&
   "NO3^-+N2^+=>NO3+N2                          ","NO3^-+O^+=>NO3+O                            ",&
   "NO3^-+O2^+=>NO3+O2                          ","NO3^-+NO^+=>NO3+NO                          ",&
   "NO3^-+N2O^+=>NO3+N2O                        ","NO3^-+NO2^+=>NO3+NO2                        ",&
   "O^-+N2^+=>O+N+N                             ","O^-+N3^+=>O+N+N2                            ",&
   "O^-+N4^+=>O+N2+N2                           ","O^-+O2^+=>O+O+O                             ",&
   "O^-+O4^+=>O+O2+O2                           ","O^-+NO^+=>O+N+O                             "/
  data reaction_sign(505:540) &
  /"O^-+N2O^+=>O+N2+O                           ","O^-+NO2^+=>O+N+O2                           ",&
   "O^-+O2^+N2=>O+O2+N2                         ","O2^-+N2^+=>O2+N+N                           ",&
   "O2^-+N3^+=>O2+N+N2                          ","O2^-+N4^+=>O2+N2+N2                         ",&
   "O2^-+O2^+=>O2+O+O                           ","O2^-+O4^+=>O2+O2+O2                         ",&
   "O2^-+NO^+=>O2+N+O                           ","O2^-+N2O^+=>O2+N2+O                         ",&
   "O2^-+NO2^+=>O2+N+O2                         ","O2^-+O2^+N2=>O2+O2+N2                       ",&
   "O3^-+N2^+=>O3+N+N                           ","O3^-+N3^+=>O3+N+N2                          ",&
   "O3^-+N4^+=>O3+N2+N2                         ","O3^-+O2^+=>O3+O+O                           ",&
   "O3^-+O4^+=>O3+O2+O2                         ","O3^-+NO^+=>O3+N+O                           ",&
   "O3^-+N2O^+=>O3+N2+O                         ","O3^-+NO2^+=>O3+N+O2                         ",&
   "O3^-+O2^+N2=>O3+O2+N2                       ","NO^-+N2^+=>NO+N+N                           ",&
   "NO^-+N3^+=>NO+N+N2                          ","NO^-+N4^+=>NO+N2+N2                         ",&
   "NO^-+O2^+=>NO+O+O                           ","NO^-+O4^+=>NO+O2+O2                         ",&
   "NO^-+NO^+=>NO+N+O                           ","NO^-+N2O^+=>NO+N2+O                         ",&
   "NO^-+NO2^+=>NO+N+O2                         ","NO^-+O2^+N2=>NO+O2+N2                       ",&
   "N2O^-+N2^+=>N2O+N+N                         ","N2O^-+N3^+=>N2O+N+N2                        ",&
   "N2O^-+N4^+=>N2O+N2+N2                       ","N2O^-+O2^+=>N2O+O+O                         ",&
   "N2O^-+O4^+=>N2O+O2+O2                       ","N2O^-+NO^+=>N2O+N+O                         "/
  data reaction_sign(541:576) &
  /"N2O^-+N2O^+=>N2O+N2+O                       ","N2O^-+NO2^+=>N2O+N+O2                       ",&
   "N2O^-+O2^+N2=>N2O+O2+N2                     ","NO2^-+N2^+=>NO2+N+N                         ",&
   "NO2^-+N3^+=>NO2+N+N2                        ","NO2^-+N4^+=>NO2+N2+N2                       ",&
   "NO2^-+O2^+=>NO2+O+O                         ","NO2^-+O4^+=>NO2+O2+O2                       ",&
   "NO2^-+NO^+=>NO2+N+O                         ","NO2^-+N2O^+=>NO2+N2+O                       ",&
   "NO2^-+NO2^+=>NO2+N+O2                       ","NO2^-+O2^+N2=>NO2+O2+N2                     ",&
   "NO3^-+N2^+=>NO3+N+N                         ","NO3^-+N3^+=>NO3+N+N2                        ",&
   "NO3^-+N4^+=>NO3+N2+N2                       ","NO3^-+O2^+=>NO3+O+O                         ",&
   "NO3^-+O4^+=>NO3+O2+O2                       ","NO3^-+NO^+=>NO3+N+O                         ",&
   "NO3^-+N2O^+=>NO3+N2+O                       ","NO3^-+NO2^+=>NO3+N+O2                       ",&
   "NO3^-+O2^+N2=>NO3+O2+N2                     ","O4^-+N^+=>O2+O2+N                           ",&
   "O4^-+N2^+=>O2+O2+N2                         ","O4^-+O^+=>O2+O2+O                           ",&
   "O4^-+O2^+=>O2+O2+O2                         ","O4^-+NO^+=>O2+O2+NO                         ",&
   "O4^-+N2O^+=>O2+O2+N2O                       ","O4^-+NO2^+=>O2+O2+NO2                       ",&
   "O4^-+N3^+=>O2+O2+N2+N                       ","O4^-+N4^+=>O2+O2+N2+N2                      ",&
   "O4^-+O4^+=>O2+O2+O2+O2                      ","O4^-+O2^+N2=>O2+O2+O2+N2                    ",&
   "O^-+N^++ANY_NEUTRAL=>O+N+ANY_NEUTRAL        ","O^-+N2^++ANY_NEUTRAL=>O+N2+ANY_NEUTRAL      ",&
   "O^-+O^++ANY_NEUTRAL=>O+O+ANY_NEUTRAL        ","O^-+O2^++ANY_NEUTRAL=>O+O2+ANY_NEUTRAL      "/
  data reaction_sign(577:612) &
  /"O^-+NO^++ANY_NEUTRAL=>O+NO+ANY_NEUTRAL      ","O2^-+N^++ANY_NEUTRAL=>O2+N+ANY_NEUTRAL      ",&
   "O2^-+N2^++ANY_NEUTRAL=>O2+N2+ANY_NEUTRAL    ","O2^-+O^++ANY_NEUTRAL=>O2+O+ANY_NEUTRAL      ",&
   "O2^-+O2^++ANY_NEUTRAL=>O2+O2+ANY_NEUTRAL    ","O2^-+NO^++ANY_NEUTRAL=>O2+NO+ANY_NEUTRAL    ",&
   "O^-+N^++ANY_NEUTRAL=>NO+ANY_NEUTRAL         ","O^-+N2^++ANY_NEUTRAL=>N2O+ANY_NEUTRAL       ",&
   "O^-+O^++ANY_NEUTRAL=>O2+ANY_NEUTRAL         ","O^-+O2^++ANY_NEUTRAL=>O3+ANY_NEUTRAL        ",&
   "O^-+NO^++ANY_NEUTRAL=>NO2+ANY_NEUTRAL       ","O2^-+N^++ANY_NEUTRAL=>NO2+ANY_NEUTRAL       ",&
   "O2^-+O^++ANY_NEUTRAL=>O3+ANY_NEUTRAL        ","O2^-+NO^++ANY_NEUTRAL=>NO3+ANY_NEUTRAL      ",&
   "O3^-+N^++ANY_NEUTRAL=>O3+N+ANY_NEUTRAL      ","O3^-+N2^++ANY_NEUTRAL=>O3+N2+ANY_NEUTRAL    ",&
   "O3^-+O^++ANY_NEUTRAL=>O3+O+ANY_NEUTRAL      ","O3^-+O2^++ANY_NEUTRAL=>O3+O2+ANY_NEUTRAL    ",&
   "O3^-+NO^++ANY_NEUTRAL=>O3+NO+ANY_NEUTRAL    ","O3^-+N2O^++ANY_NEUTRAL=>O3+N2O+ANY_NEUTRAL  ",&
   "O3^-+NO2^++ANY_NEUTRAL=>O3+NO2+ANY_NEUTRAL  ","NO^-+N^++ANY_NEUTRAL=>NO+N+ANY_NEUTRAL      ",&
   "NO^-+N2^++ANY_NEUTRAL=>NO+N2+ANY_NEUTRAL    ","NO^-+O^++ANY_NEUTRAL=>NO+O+ANY_NEUTRAL      ",&
   "NO^-+O2^++ANY_NEUTRAL=>NO+O2+ANY_NEUTRAL    ","NO^-+NO^++ANY_NEUTRAL=>NO+NO+ANY_NEUTRAL    ",&
   "NO^-+N2O^++ANY_NEUTRAL=>NO+N2O+ANY_NEUTRAL  ","NO^-+NO2^++ANY_NEUTRAL=>NO+NO2+ANY_NEUTRAL  ",&
   "N2O^-+N^++ANY_NEUTRAL=>N2O+N+ANY_NEUTRAL    ","N2O^-+N2^++ANY_NEUTRAL=>N2O+N2+ANY_NEUTRAL  ",&
   "N2O^-+O^++ANY_NEUTRAL=>N2O+O+ANY_NEUTRAL    ","N2O^-+O2^++ANY_NEUTRAL=>N2O+O2+ANY_NEUTRAL  ",&
   "N2O^-+NO^++ANY_NEUTRAL=>N2O+NO+ANY_NEUTRAL  ","N2O^-+N2O^++ANY_NEUTRAL=>N2O+N2O+ANY_NEUTRAL",&
   "N2O^-+NO2^++ANY_NEUTRAL=>N2O+NO2+ANY_NEUTRAL","NO2^-+N^++ANY_NEUTRAL=>NO2+N+ANY_NEUTRAL    "/
  data reaction_sign(613:625) &
  /"NO2^-+N2^++ANY_NEUTRAL=>NO2+N2+ANY_NEUTRAL  ","NO2^-+O^++ANY_NEUTRAL=>NO2+O+ANY_NEUTRAL    ",&
   "NO2^-+O2^++ANY_NEUTRAL=>NO2+O2+ANY_NEUTRAL  ","NO2^-+NO^++ANY_NEUTRAL=>NO2+NO+ANY_NEUTRAL  ",&
   "NO2^-+N2O^++ANY_NEUTRAL=>NO2+N2O+ANY_NEUTRAL","NO2^-+NO2^++ANY_NEUTRAL=>NO2+NO2+ANY_NEUTRAL",&
   "NO3^-+N^++ANY_NEUTRAL=>NO3+N+ANY_NEUTRAL    ","NO3^-+N2^++ANY_NEUTRAL=>NO3+N2+ANY_NEUTRAL  ",&
   "NO3^-+O^++ANY_NEUTRAL=>NO3+O+ANY_NEUTRAL    ","NO3^-+O2^++ANY_NEUTRAL=>NO3+O2+ANY_NEUTRAL  ",&
   "NO3^-+NO^++ANY_NEUTRAL=>NO3+NO+ANY_NEUTRAL  ","NO3^-+N2O^++ANY_NEUTRAL=>NO3+N2O+ANY_NEUTRAL",&
   "NO3^-+NO2^++ANY_NEUTRAL=>NO3+NO2+ANY_NEUTRAL"/
  data bolsig_species(1:bolsig_species_max) &
  /"N2","O2"/
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
211 format(i3,1x,A44)
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
311 format(491x,53(1x,i9))
312 format(A3,1x,A44,1x,53(1x,A9))
313 format(i3,1x,A44,1x,53(1x,1pd9.2))
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
      write(ifile_unit,"(625(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,53(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
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
    write(ifile_unit,"(1x,A12,625(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
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
      write(ifile_unit,"(1pe15.6,53(1pe13.4))") densav(0,2), densav(1:,2)
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
      write(ifile_unit,"(626(1pe13.4))") densav(0,2), rrt_loc(:)
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
  reac_source_local(01,001) = + reac_rate_local(001) 
  reac_source_local(02,001) = - reac_rate_local(001) 
  reac_source_local(02,002) = + reac_rate_local(002) 
  reac_source_local(03,002) = - reac_rate_local(002) 
  reac_source_local(03,003) = + reac_rate_local(003) 
  reac_source_local(04,003) = - reac_rate_local(003) 
  reac_source_local(04,004) = + reac_rate_local(004) 
  reac_source_local(05,004) = - reac_rate_local(004) 
  reac_source_local(05,005) = + reac_rate_local(005) 
  reac_source_local(06,005) = - reac_rate_local(005) 
  reac_source_local(06,006) = + reac_rate_local(006) 
  reac_source_local(07,006) = - reac_rate_local(006) 
  reac_source_local(07,007) = + reac_rate_local(007) 
  reac_source_local(08,007) = - reac_rate_local(007) 
  reac_source_local(08,008) = + reac_rate_local(008) 
  reac_source_local(09,008) = - reac_rate_local(008) 
  reac_source_local(01,009) = - reac_rate_local(009) 
  reac_source_local(02,009) = + reac_rate_local(009) 
  reac_source_local(02,010) = - reac_rate_local(010) 
  reac_source_local(03,010) = + reac_rate_local(010) 
  reac_source_local(03,011) = - reac_rate_local(011) 
  reac_source_local(04,011) = + reac_rate_local(011) 
  reac_source_local(04,012) = - reac_rate_local(012) 
  reac_source_local(05,012) = + reac_rate_local(012) 
  reac_source_local(05,013) = - reac_rate_local(013) 
  reac_source_local(06,013) = + reac_rate_local(013) 
  reac_source_local(06,014) = - reac_rate_local(014) 
  reac_source_local(07,014) = + reac_rate_local(014) 
  reac_source_local(07,015) = - reac_rate_local(015) 
  reac_source_local(08,015) = + reac_rate_local(015) 
  reac_source_local(08,016) = - reac_rate_local(016) 
  reac_source_local(09,016) = + reac_rate_local(016) 
  reac_source_local(01,017) = + reac_rate_local(017) 
  reac_source_local(02,017) = - reac_rate_local(017) 
  reac_source_local(02,018) = + reac_rate_local(018) 
  reac_source_local(03,018) = - reac_rate_local(018) 
  reac_source_local(03,019) = + reac_rate_local(019) 
  reac_source_local(04,019) = - reac_rate_local(019) 
  reac_source_local(04,020) = + reac_rate_local(020) 
  reac_source_local(05,020) = - reac_rate_local(020) 
  reac_source_local(05,021) = + reac_rate_local(021) 
  reac_source_local(06,021) = - reac_rate_local(021) 
  reac_source_local(06,022) = + reac_rate_local(022) 
  reac_source_local(07,022) = - reac_rate_local(022) 
  reac_source_local(07,023) = + reac_rate_local(023) 
  reac_source_local(08,023) = - reac_rate_local(023) 
  reac_source_local(08,024) = + reac_rate_local(024) 
  reac_source_local(09,024) = - reac_rate_local(024) 
  reac_source_local(01,025) = - reac_rate_local(025) 
  reac_source_local(02,025) = + reac_rate_local(025) 
  reac_source_local(02,026) = - reac_rate_local(026) 
  reac_source_local(03,026) = + reac_rate_local(026) 
  reac_source_local(03,027) = - reac_rate_local(027) 
  reac_source_local(04,027) = + reac_rate_local(027) 
  reac_source_local(04,028) = - reac_rate_local(028) 
  reac_source_local(05,028) = + reac_rate_local(028) 
  reac_source_local(05,029) = - reac_rate_local(029) 
  reac_source_local(06,029) = + reac_rate_local(029) 
  reac_source_local(06,030) = - reac_rate_local(030) 
  reac_source_local(07,030) = + reac_rate_local(030) 
  reac_source_local(07,031) = - reac_rate_local(031) 
  reac_source_local(08,031) = + reac_rate_local(031) 
  reac_source_local(08,032) = - reac_rate_local(032) 
  reac_source_local(09,032) = + reac_rate_local(032) 
  reac_source_local(01,033) = + reac_rate_local(033) 
  reac_source_local(02,033) = - reac_rate_local(033) 
  reac_source_local(02,034) = + reac_rate_local(034) 
  reac_source_local(03,034) = - reac_rate_local(034) 
  reac_source_local(03,035) = + reac_rate_local(035) 
  reac_source_local(04,035) = - reac_rate_local(035) 
  reac_source_local(04,036) = + reac_rate_local(036) 
  reac_source_local(05,036) = - reac_rate_local(036) 
  reac_source_local(05,037) = + reac_rate_local(037) 
  reac_source_local(06,037) = - reac_rate_local(037) 
  reac_source_local(06,038) = + reac_rate_local(038) 
  reac_source_local(07,038) = - reac_rate_local(038) 
  reac_source_local(07,039) = + reac_rate_local(039) 
  reac_source_local(08,039) = - reac_rate_local(039) 
  reac_source_local(08,040) = + reac_rate_local(040) 
  reac_source_local(09,040) = - reac_rate_local(040) 
  reac_source_local(01,041) = - reac_rate_local(041) 
  reac_source_local(02,041) = + reac_rate_local(041) 
  reac_source_local(02,042) = - reac_rate_local(042) 
  reac_source_local(03,042) = + reac_rate_local(042) 
  reac_source_local(03,043) = - reac_rate_local(043) 
  reac_source_local(04,043) = + reac_rate_local(043) 
  reac_source_local(04,044) = - reac_rate_local(044) 
  reac_source_local(05,044) = + reac_rate_local(044) 
  reac_source_local(05,045) = - reac_rate_local(045) 
  reac_source_local(06,045) = + reac_rate_local(045) 
  reac_source_local(06,046) = - reac_rate_local(046) 
  reac_source_local(07,046) = + reac_rate_local(046) 
  reac_source_local(07,047) = - reac_rate_local(047) 
  reac_source_local(08,047) = + reac_rate_local(047) 
  reac_source_local(08,048) = - reac_rate_local(048) 
  reac_source_local(09,048) = + reac_rate_local(048) 
  reac_source_local(21,049) = + reac_rate_local(049) 
  reac_source_local(22,049) = - reac_rate_local(049) 
  reac_source_local(22,050) = + reac_rate_local(050) 
  reac_source_local(23,050) = - reac_rate_local(050) 
  reac_source_local(23,051) = + reac_rate_local(051) 
  reac_source_local(24,051) = - reac_rate_local(051) 
  reac_source_local(24,052) = + reac_rate_local(052) 
  reac_source_local(25,052) = - reac_rate_local(052) 
  reac_source_local(21,053) = - reac_rate_local(053) 
  reac_source_local(22,053) = + reac_rate_local(053) 
  reac_source_local(22,054) = - reac_rate_local(054) 
  reac_source_local(23,054) = + reac_rate_local(054) 
  reac_source_local(23,055) = - reac_rate_local(055) 
  reac_source_local(24,055) = + reac_rate_local(055) 
  reac_source_local(24,056) = - reac_rate_local(056) 
  reac_source_local(25,056) = + reac_rate_local(056) 
  reac_source_local(21,057) = + reac_rate_local(057) 
  reac_source_local(22,057) = - reac_rate_local(057) 
  reac_source_local(22,058) = + reac_rate_local(058) 
  reac_source_local(23,058) = - reac_rate_local(058) 
  reac_source_local(23,059) = + reac_rate_local(059) 
  reac_source_local(24,059) = - reac_rate_local(059) 
  reac_source_local(24,060) = + reac_rate_local(060) 
  reac_source_local(25,060) = - reac_rate_local(060) 
  reac_source_local(21,061) = - reac_rate_local(061) 
  reac_source_local(22,061) = + reac_rate_local(061) 
  reac_source_local(22,062) = - reac_rate_local(062) 
  reac_source_local(23,062) = + reac_rate_local(062) 
  reac_source_local(23,063) = - reac_rate_local(063) 
  reac_source_local(24,063) = + reac_rate_local(063) 
  reac_source_local(24,064) = - reac_rate_local(064) 
  reac_source_local(25,064) = + reac_rate_local(064) 
  reac_source_local(01,065) = - reac_rate_local(065) 
  reac_source_local(10,065) = + reac_rate_local(065) 
  reac_source_local(01,066) = - reac_rate_local(066) 
  reac_source_local(12,066) = + reac_rate_local(066) 
  reac_source_local(01,067) = - reac_rate_local(067) 
  reac_source_local(12,067) = + reac_rate_local(067) 
  reac_source_local(01,068) = - reac_rate_local(068) 
  reac_source_local(12,068) = + reac_rate_local(068) 
  reac_source_local(01,069) = - reac_rate_local(069) 
  reac_source_local(13,069) = + reac_rate_local(069) 
  reac_source_local(01,070) = - reac_rate_local(070) 
  reac_source_local(13,070) = + reac_rate_local(070) 
  reac_source_local(01,071) = - reac_rate_local(071) 
  reac_source_local(13,071) = + reac_rate_local(071) 
  reac_source_local(01,072) = - reac_rate_local(072) 
  reac_source_local(14,072) = + reac_rate_local(072) 
  reac_source_local(15,072) = + reac_rate_local(072) 
  reac_source_local(21,073) = - reac_rate_local(073) 
  reac_source_local(26,073) = + reac_rate_local(073) 
  reac_source_local(21,074) = - reac_rate_local(074) 
  reac_source_local(27,074) = + reac_rate_local(074) 
  reac_source_local(21,075) = - reac_rate_local(075) 
  reac_source_local(28,075) = + reac_rate_local(075) 
  reac_source_local(21,076) = - reac_rate_local(076) 
  reac_source_local(29,076) = + reac_rate_local(076) * 2.d0
  reac_source_local(21,077) = - reac_rate_local(077) 
  reac_source_local(29,077) = + reac_rate_local(077) 
  reac_source_local(30,077) = + reac_rate_local(077) 
  reac_source_local(21,078) = - reac_rate_local(078) 
  reac_source_local(29,078) = + reac_rate_local(078) 
  reac_source_local(31,078) = + reac_rate_local(078) 
  reac_source_local(26,079) = - reac_rate_local(079) 
  reac_source_local(29,079) = + reac_rate_local(079) * 2.d0
  reac_source_local(29,080) = - reac_rate_local(080) 
  reac_source_local(30,080) = + reac_rate_local(080) 
  reac_source_local(29,081) = - reac_rate_local(081) 
  reac_source_local(31,081) = + reac_rate_local(081) 
  reac_source_local(01,082) = + reac_rate_local(082) 
  reac_source_local(10,082) = - reac_rate_local(082) 
  reac_source_local(21,083) = + reac_rate_local(083) 
  reac_source_local(26,083) = - reac_rate_local(083) 
  reac_source_local(14,084) = - reac_rate_local(084) 
  reac_source_local(17,084) = + reac_rate_local(084) 
  reac_source_local(53,084) = + reac_rate_local(084) 
  reac_source_local(29,085) = - reac_rate_local(085) 
  reac_source_local(33,085) = + reac_rate_local(085) 
  reac_source_local(53,085) = + reac_rate_local(085) 
  reac_source_local(01,086) = - reac_rate_local(086) 
  reac_source_local(18,086) = + reac_rate_local(086) 
  reac_source_local(53,086) = + reac_rate_local(086) 
  reac_source_local(10,087) = - reac_rate_local(087) 
  reac_source_local(18,087) = + reac_rate_local(087) 
  reac_source_local(53,087) = + reac_rate_local(087) 
  reac_source_local(21,088) = - reac_rate_local(088) 
  reac_source_local(34,088) = + reac_rate_local(088) 
  reac_source_local(53,088) = + reac_rate_local(088) 
  reac_source_local(26,089) = - reac_rate_local(089) 
  reac_source_local(34,089) = + reac_rate_local(089) 
  reac_source_local(53,089) = + reac_rate_local(089) 
  reac_source_local(40,090) = - reac_rate_local(090) 
  reac_source_local(45,090) = + reac_rate_local(090) 
  reac_source_local(53,090) = + reac_rate_local(090) 
  reac_source_local(41,091) = - reac_rate_local(091) 
  reac_source_local(46,091) = + reac_rate_local(091) 
  reac_source_local(53,091) = + reac_rate_local(091) 
  reac_source_local(14,092) = + reac_rate_local(092) * 2.d0
  reac_source_local(18,092) = - reac_rate_local(092) 
  reac_source_local(53,092) = - reac_rate_local(092) 
  reac_source_local(14,093) = + reac_rate_local(093) 
  reac_source_local(15,093) = + reac_rate_local(093) 
  reac_source_local(18,093) = - reac_rate_local(093) 
  reac_source_local(53,093) = - reac_rate_local(093) 
  reac_source_local(14,094) = + reac_rate_local(094) 
  reac_source_local(16,094) = + reac_rate_local(094) 
  reac_source_local(18,094) = - reac_rate_local(094) 
  reac_source_local(53,094) = - reac_rate_local(094) 
  reac_source_local(29,095) = + reac_rate_local(095) * 2.d0
  reac_source_local(34,095) = - reac_rate_local(095) 
  reac_source_local(53,095) = - reac_rate_local(095) 
  reac_source_local(29,096) = + reac_rate_local(096) 
  reac_source_local(30,096) = + reac_rate_local(096) 
  reac_source_local(34,096) = - reac_rate_local(096) 
  reac_source_local(53,096) = - reac_rate_local(096) 
  reac_source_local(29,097) = + reac_rate_local(097) 
  reac_source_local(31,097) = + reac_rate_local(097) 
  reac_source_local(34,097) = - reac_rate_local(097) 
  reac_source_local(53,097) = - reac_rate_local(097) 
  reac_source_local(14,098) = + reac_rate_local(098) 
  reac_source_local(29,098) = + reac_rate_local(098) 
  reac_source_local(45,098) = - reac_rate_local(098) 
  reac_source_local(53,098) = - reac_rate_local(098) 
  reac_source_local(15,099) = + reac_rate_local(099) 
  reac_source_local(29,099) = + reac_rate_local(099) 
  reac_source_local(45,099) = - reac_rate_local(099) 
  reac_source_local(53,099) = - reac_rate_local(099) 
  reac_source_local(01,100) = + reac_rate_local(100) 
  reac_source_local(14,100) = + reac_rate_local(100) 
  reac_source_local(19,100) = - reac_rate_local(100) 
  reac_source_local(53,100) = - reac_rate_local(100) 
  reac_source_local(01,101) = + reac_rate_local(101) * 2.d0
  reac_source_local(20,101) = - reac_rate_local(101) 
  reac_source_local(53,101) = - reac_rate_local(101) 
  reac_source_local(01,102) = + reac_rate_local(102) 
  reac_source_local(29,102) = + reac_rate_local(102) 
  reac_source_local(46,102) = - reac_rate_local(102) 
  reac_source_local(53,102) = - reac_rate_local(102) 
  reac_source_local(29,103) = + reac_rate_local(103) 
  reac_source_local(40,103) = + reac_rate_local(103) 
  reac_source_local(47,103) = - reac_rate_local(103) 
  reac_source_local(53,103) = - reac_rate_local(103) 
  reac_source_local(21,104) = + reac_rate_local(104) * 2.d0
  reac_source_local(35,104) = - reac_rate_local(104) 
  reac_source_local(53,104) = - reac_rate_local(104) 
  reac_source_local(01,105) = + reac_rate_local(105) 
  reac_source_local(21,105) = + reac_rate_local(105) 
  reac_source_local(52,105) = - reac_rate_local(105) 
  reac_source_local(53,105) = - reac_rate_local(105) 
  reac_source_local(14,106) = + reac_rate_local(106) 
  reac_source_local(17,106) = - reac_rate_local(106) 
  reac_source_local(53,106) = - reac_rate_local(106) 
  reac_source_local(29,107) = + reac_rate_local(107) 
  reac_source_local(33,107) = - reac_rate_local(107) 
  reac_source_local(53,107) = - reac_rate_local(107) 
  reac_source_local(14,108) = + reac_rate_local(108) 
  reac_source_local(17,108) = - reac_rate_local(108) 
  reac_source_local(53,108) = - reac_rate_local(108) 
  reac_source_local(29,109) = + reac_rate_local(109) 
  reac_source_local(33,109) = - reac_rate_local(109) 
  reac_source_local(53,109) = - reac_rate_local(109) 
  reac_source_local(21,110) = - reac_rate_local(110) 
  reac_source_local(29,110) = + reac_rate_local(110) 
  reac_source_local(36,110) = + reac_rate_local(110) 
  reac_source_local(53,110) = - reac_rate_local(110) 
  reac_source_local(14,111) = + reac_rate_local(111) 
  reac_source_local(36,111) = + reac_rate_local(111) 
  reac_source_local(40,111) = - reac_rate_local(111) 
  reac_source_local(53,111) = - reac_rate_local(111) 
  reac_source_local(21,112) = + reac_rate_local(112) 
  reac_source_local(32,112) = - reac_rate_local(112) 
  reac_source_local(36,112) = + reac_rate_local(112) 
  reac_source_local(53,112) = - reac_rate_local(112) 
  reac_source_local(29,113) = + reac_rate_local(113) 
  reac_source_local(32,113) = - reac_rate_local(113) 
  reac_source_local(37,113) = + reac_rate_local(113) 
  reac_source_local(53,113) = - reac_rate_local(113) 
  reac_source_local(14,114) = + reac_rate_local(114) 
  reac_source_local(41,114) = - reac_rate_local(114) 
  reac_source_local(48,114) = + reac_rate_local(114) 
  reac_source_local(53,114) = - reac_rate_local(114) 
  reac_source_local(21,115) = - reac_rate_local(115) 
  reac_source_local(37,115) = + reac_rate_local(115) 
  reac_source_local(53,115) = - reac_rate_local(115) 
  reac_source_local(36,116) = + reac_rate_local(116) 
  reac_source_local(40,116) = + reac_rate_local(116) 
  reac_source_local(42,116) = - reac_rate_local(116) 
  reac_source_local(53,116) = - reac_rate_local(116) 
  reac_source_local(29,117) = - reac_rate_local(117) 
  reac_source_local(36,117) = + reac_rate_local(117) 
  reac_source_local(53,117) = - reac_rate_local(117) 
  reac_source_local(21,118) = - reac_rate_local(118) 
  reac_source_local(37,118) = + reac_rate_local(118) 
  reac_source_local(53,118) = - reac_rate_local(118) 
  reac_source_local(32,119) = - reac_rate_local(119) 
  reac_source_local(38,119) = + reac_rate_local(119) 
  reac_source_local(53,119) = - reac_rate_local(119) 
  reac_source_local(40,120) = - reac_rate_local(120) 
  reac_source_local(48,120) = + reac_rate_local(120) 
  reac_source_local(53,120) = - reac_rate_local(120) 
  reac_source_local(41,121) = - reac_rate_local(121) 
  reac_source_local(49,121) = + reac_rate_local(121) 
  reac_source_local(53,121) = - reac_rate_local(121) 
  reac_source_local(21,122) = - reac_rate_local(122) 
  reac_source_local(37,122) = + reac_rate_local(122) 
  reac_source_local(53,122) = - reac_rate_local(122) 
  reac_source_local(21,123) = + reac_rate_local(123) 
  reac_source_local(29,123) = - reac_rate_local(123) 
  reac_source_local(36,123) = - reac_rate_local(123) 
  reac_source_local(53,123) = + reac_rate_local(123) 
  reac_source_local(14,124) = - reac_rate_local(124) 
  reac_source_local(36,124) = - reac_rate_local(124) 
  reac_source_local(40,124) = + reac_rate_local(124) 
  reac_source_local(53,124) = + reac_rate_local(124) 
  reac_source_local(36,125) = - reac_rate_local(125) 
  reac_source_local(40,125) = - reac_rate_local(125) 
  reac_source_local(42,125) = + reac_rate_local(125) 
  reac_source_local(53,125) = + reac_rate_local(125) 
  reac_source_local(01,126) = - reac_rate_local(126) 
  reac_source_local(36,126) = - reac_rate_local(126) 
  reac_source_local(41,126) = + reac_rate_local(126) 
  reac_source_local(53,126) = + reac_rate_local(126) 
  reac_source_local(21,127) = - reac_rate_local(127) 
  reac_source_local(32,127) = + reac_rate_local(127) 
  reac_source_local(36,127) = - reac_rate_local(127) 
  reac_source_local(53,127) = + reac_rate_local(127) 
  reac_source_local(26,128) = - reac_rate_local(128) 
  reac_source_local(32,128) = + reac_rate_local(128) 
  reac_source_local(36,128) = - reac_rate_local(128) 
  reac_source_local(53,128) = + reac_rate_local(128) 
  reac_source_local(21,129) = + reac_rate_local(129) 
  reac_source_local(27,129) = - reac_rate_local(129) 
  reac_source_local(29,129) = + reac_rate_local(129) 
  reac_source_local(36,129) = - reac_rate_local(129) 
  reac_source_local(53,129) = + reac_rate_local(129) 
  reac_source_local(01,130) = + reac_rate_local(130) 
  reac_source_local(10,130) = - reac_rate_local(130) 
  reac_source_local(29,130) = + reac_rate_local(130) 
  reac_source_local(36,130) = - reac_rate_local(130) 
  reac_source_local(53,130) = + reac_rate_local(130) 
  reac_source_local(01,131) = + reac_rate_local(131) 
  reac_source_local(11,131) = - reac_rate_local(131) 
  reac_source_local(29,131) = + reac_rate_local(131) 
  reac_source_local(36,131) = - reac_rate_local(131) 
  reac_source_local(53,131) = + reac_rate_local(131) 
  reac_source_local(21,132) = + reac_rate_local(132) * 2.d0
  reac_source_local(32,132) = - reac_rate_local(132) 
  reac_source_local(36,132) = - reac_rate_local(132) 
  reac_source_local(53,132) = + reac_rate_local(132) 
  reac_source_local(29,133) = - reac_rate_local(133) 
  reac_source_local(32,133) = + reac_rate_local(133) 
  reac_source_local(37,133) = - reac_rate_local(133) 
  reac_source_local(53,133) = + reac_rate_local(133) 
  reac_source_local(14,134) = - reac_rate_local(134) 
  reac_source_local(37,134) = - reac_rate_local(134) 
  reac_source_local(42,134) = + reac_rate_local(134) 
  reac_source_local(53,134) = + reac_rate_local(134) 
  reac_source_local(21,135) = + reac_rate_local(135) 
  reac_source_local(37,135) = - reac_rate_local(135) 
  reac_source_local(53,135) = + reac_rate_local(135) 
  reac_source_local(21,136) = + reac_rate_local(136) * 2.d0
  reac_source_local(26,136) = - reac_rate_local(136) 
  reac_source_local(37,136) = - reac_rate_local(136) 
  reac_source_local(53,136) = + reac_rate_local(136) 
  reac_source_local(21,137) = + reac_rate_local(137) * 2.d0
  reac_source_local(27,137) = - reac_rate_local(137) 
  reac_source_local(37,137) = - reac_rate_local(137) 
  reac_source_local(53,137) = + reac_rate_local(137) 
  reac_source_local(21,138) = + reac_rate_local(138) 
  reac_source_local(37,138) = - reac_rate_local(138) 
  reac_source_local(53,138) = + reac_rate_local(138) 
  reac_source_local(01,139) = + reac_rate_local(139) 
  reac_source_local(10,139) = - reac_rate_local(139) 
  reac_source_local(21,139) = + reac_rate_local(139) 
  reac_source_local(37,139) = - reac_rate_local(139) 
  reac_source_local(53,139) = + reac_rate_local(139) 
  reac_source_local(01,140) = + reac_rate_local(140) 
  reac_source_local(11,140) = - reac_rate_local(140) 
  reac_source_local(21,140) = + reac_rate_local(140) 
  reac_source_local(37,140) = - reac_rate_local(140) 
  reac_source_local(53,140) = + reac_rate_local(140) 
  reac_source_local(21,141) = + reac_rate_local(141) * 2.d0
  reac_source_local(29,141) = - reac_rate_local(141) 
  reac_source_local(38,141) = - reac_rate_local(141) 
  reac_source_local(53,141) = + reac_rate_local(141) 
  reac_source_local(14,142) = - reac_rate_local(142) 
  reac_source_local(41,142) = + reac_rate_local(142) 
  reac_source_local(48,142) = - reac_rate_local(142) 
  reac_source_local(53,142) = + reac_rate_local(142) 
  reac_source_local(14,143) = - reac_rate_local(143) 
  reac_source_local(21,143) = + reac_rate_local(143) 
  reac_source_local(38,143) = - reac_rate_local(143) 
  reac_source_local(40,143) = + reac_rate_local(143) 
  reac_source_local(53,143) = + reac_rate_local(143) 
  reac_source_local(01,144) = + reac_rate_local(144) 
  reac_source_local(14,144) = - reac_rate_local(144) 
  reac_source_local(40,144) = + reac_rate_local(144) 
  reac_source_local(49,144) = - reac_rate_local(144) 
  reac_source_local(53,144) = + reac_rate_local(144) 
  reac_source_local(14,145) = - reac_rate_local(145) 
  reac_source_local(40,145) = + reac_rate_local(145) * 2.d0
  reac_source_local(50,145) = - reac_rate_local(145) 
  reac_source_local(53,145) = + reac_rate_local(145) 
  reac_source_local(14,146) = - reac_rate_local(146) 
  reac_source_local(40,146) = + reac_rate_local(146) 
  reac_source_local(42,146) = + reac_rate_local(146) 
  reac_source_local(51,146) = - reac_rate_local(146) 
  reac_source_local(53,146) = + reac_rate_local(146) 
  reac_source_local(29,147) = - reac_rate_local(147) 
  reac_source_local(42,147) = + reac_rate_local(147) 
  reac_source_local(48,147) = - reac_rate_local(147) 
  reac_source_local(53,147) = + reac_rate_local(147) 
  reac_source_local(29,148) = - reac_rate_local(148) 
  reac_source_local(40,148) = + reac_rate_local(148) * 2.d0
  reac_source_local(49,148) = - reac_rate_local(148) 
  reac_source_local(53,148) = + reac_rate_local(148) 
  reac_source_local(21,149) = + reac_rate_local(149) 
  reac_source_local(29,149) = - reac_rate_local(149) 
  reac_source_local(40,149) = + reac_rate_local(149) 
  reac_source_local(50,149) = - reac_rate_local(149) 
  reac_source_local(53,149) = + reac_rate_local(149) 
  reac_source_local(29,150) = - reac_rate_local(150) 
  reac_source_local(32,150) = + reac_rate_local(150) 
  reac_source_local(40,150) = + reac_rate_local(150) 
  reac_source_local(51,150) = - reac_rate_local(150) 
  reac_source_local(53,150) = + reac_rate_local(150) 
  reac_source_local(01,151) = + reac_rate_local(151) 
  reac_source_local(10,151) = - reac_rate_local(151) 
  reac_source_local(32,151) = + reac_rate_local(151) 
  reac_source_local(38,151) = - reac_rate_local(151) 
  reac_source_local(53,151) = + reac_rate_local(151) 
  reac_source_local(01,152) = + reac_rate_local(152) 
  reac_source_local(10,152) = - reac_rate_local(152) 
  reac_source_local(40,152) = + reac_rate_local(152) 
  reac_source_local(48,152) = - reac_rate_local(152) 
  reac_source_local(53,152) = + reac_rate_local(152) 
  reac_source_local(01,153) = + reac_rate_local(153) 
  reac_source_local(10,153) = - reac_rate_local(153) 
  reac_source_local(41,153) = + reac_rate_local(153) 
  reac_source_local(49,153) = - reac_rate_local(153) 
  reac_source_local(53,153) = + reac_rate_local(153) 
  reac_source_local(01,154) = + reac_rate_local(154) 
  reac_source_local(10,154) = - reac_rate_local(154) 
  reac_source_local(42,154) = + reac_rate_local(154) 
  reac_source_local(50,154) = - reac_rate_local(154) 
  reac_source_local(53,154) = + reac_rate_local(154) 
  reac_source_local(01,155) = + reac_rate_local(155) 
  reac_source_local(10,155) = - reac_rate_local(155) 
  reac_source_local(43,155) = + reac_rate_local(155) 
  reac_source_local(51,155) = - reac_rate_local(155) 
  reac_source_local(53,155) = + reac_rate_local(155) 
  reac_source_local(01,156) = + reac_rate_local(156) 
  reac_source_local(11,156) = - reac_rate_local(156) 
  reac_source_local(32,156) = + reac_rate_local(156) 
  reac_source_local(38,156) = - reac_rate_local(156) 
  reac_source_local(53,156) = + reac_rate_local(156) 
  reac_source_local(01,157) = + reac_rate_local(157) 
  reac_source_local(11,157) = - reac_rate_local(157) 
  reac_source_local(40,157) = + reac_rate_local(157) 
  reac_source_local(48,157) = - reac_rate_local(157) 
  reac_source_local(53,157) = + reac_rate_local(157) 
  reac_source_local(01,158) = + reac_rate_local(158) 
  reac_source_local(11,158) = - reac_rate_local(158) 
  reac_source_local(41,158) = + reac_rate_local(158) 
  reac_source_local(49,158) = - reac_rate_local(158) 
  reac_source_local(53,158) = + reac_rate_local(158) 
  reac_source_local(01,159) = + reac_rate_local(159) 
  reac_source_local(11,159) = - reac_rate_local(159) 
  reac_source_local(42,159) = + reac_rate_local(159) 
  reac_source_local(50,159) = - reac_rate_local(159) 
  reac_source_local(53,159) = + reac_rate_local(159) 
  reac_source_local(01,160) = + reac_rate_local(160) 
  reac_source_local(11,160) = - reac_rate_local(160) 
  reac_source_local(43,160) = + reac_rate_local(160) 
  reac_source_local(51,160) = - reac_rate_local(160) 
  reac_source_local(53,160) = + reac_rate_local(160) 
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
  reac_source_local(10,176) = - reac_rate_local(176) 
  reac_source_local(21,176) = - reac_rate_local(176) 
  reac_source_local(29,176) = + reac_rate_local(176) 
  reac_source_local(41,176) = + reac_rate_local(176) 
  reac_source_local(01,177) = + reac_rate_local(177) 
  reac_source_local(10,177) = - reac_rate_local(177) 
  reac_source_local(01,178) = + reac_rate_local(178) 
  reac_source_local(10,178) = - reac_rate_local(178) 
  reac_source_local(01,179) = + reac_rate_local(179) 
  reac_source_local(10,179) = - reac_rate_local(179) 
  reac_source_local(14,179) = + reac_rate_local(179) 
  reac_source_local(40,179) = + reac_rate_local(179) 
  reac_source_local(41,179) = - reac_rate_local(179) 
  reac_source_local(01,180) = + reac_rate_local(180) 
  reac_source_local(10,180) = - reac_rate_local(180) 
  reac_source_local(29,180) = + reac_rate_local(180) 
  reac_source_local(40,180) = + reac_rate_local(180) 
  reac_source_local(42,180) = - reac_rate_local(180) 
  reac_source_local(01,181) = + reac_rate_local(181) 
  reac_source_local(10,181) = - reac_rate_local(181) * 2.d0
  reac_source_local(11,181) = + reac_rate_local(181) 
  reac_source_local(01,182) = + reac_rate_local(182) 
  reac_source_local(10,182) = - reac_rate_local(182) * 2.d0
  reac_source_local(13,182) = + reac_rate_local(182) 
  reac_source_local(10,183) = + reac_rate_local(183) 
  reac_source_local(11,183) = - reac_rate_local(183) 
  reac_source_local(01,184) = + reac_rate_local(184) 
  reac_source_local(11,184) = - reac_rate_local(184) 
  reac_source_local(01,185) = + reac_rate_local(185) 
  reac_source_local(11,185) = - reac_rate_local(185) 
  reac_source_local(21,185) = - reac_rate_local(185) 
  reac_source_local(29,185) = + reac_rate_local(185) * 2.d0
  reac_source_local(10,186) = + reac_rate_local(186) 
  reac_source_local(11,186) = - reac_rate_local(186) 
  reac_source_local(12,187) = + reac_rate_local(187) 
  reac_source_local(13,187) = - reac_rate_local(187) 
  reac_source_local(01,188) = + reac_rate_local(188) 
  reac_source_local(13,188) = - reac_rate_local(188) 
  reac_source_local(21,188) = - reac_rate_local(188) 
  reac_source_local(29,188) = + reac_rate_local(188) 
  reac_source_local(31,188) = + reac_rate_local(188) 
  reac_source_local(11,189) = + reac_rate_local(189) 
  reac_source_local(12,189) = - reac_rate_local(189) 
  reac_source_local(01,190) = + reac_rate_local(190) 
  reac_source_local(12,190) = - reac_rate_local(190) 
  reac_source_local(21,190) = - reac_rate_local(190) 
  reac_source_local(29,190) = + reac_rate_local(190) * 2.d0
  reac_source_local(01,191) = + reac_rate_local(191) 
  reac_source_local(12,191) = - reac_rate_local(191) 
  reac_source_local(14,191) = + reac_rate_local(191) 
  reac_source_local(29,191) = + reac_rate_local(191) 
  reac_source_local(40,191) = - reac_rate_local(191) 
  reac_source_local(10,192) = - reac_rate_local(192) 
  reac_source_local(12,192) = - reac_rate_local(192) 
  reac_source_local(20,192) = + reac_rate_local(192) 
  reac_source_local(53,192) = + reac_rate_local(192) 
  reac_source_local(12,193) = - reac_rate_local(193) * 2.d0
  reac_source_local(20,193) = + reac_rate_local(193) 
  reac_source_local(53,193) = + reac_rate_local(193) 
  reac_source_local(10,194) = + reac_rate_local(194) 
  reac_source_local(14,194) = - reac_rate_local(194) * 2.d0
  reac_source_local(10,195) = + reac_rate_local(195) 
  reac_source_local(14,195) = - reac_rate_local(195) * 2.d0
  reac_source_local(10,196) = + reac_rate_local(196) 
  reac_source_local(14,196) = - reac_rate_local(196) * 2.d0
  reac_source_local(10,197) = + reac_rate_local(197) 
  reac_source_local(14,197) = - reac_rate_local(197) * 2.d0
  reac_source_local(10,198) = + reac_rate_local(198) 
  reac_source_local(14,198) = - reac_rate_local(198) * 2.d0
  reac_source_local(11,199) = + reac_rate_local(199) 
  reac_source_local(14,199) = - reac_rate_local(199) * 2.d0
  reac_source_local(11,200) = + reac_rate_local(200) 
  reac_source_local(14,200) = - reac_rate_local(200) * 2.d0
  reac_source_local(11,201) = + reac_rate_local(201) 
  reac_source_local(14,201) = - reac_rate_local(201) * 2.d0
  reac_source_local(11,202) = + reac_rate_local(202) 
  reac_source_local(14,202) = - reac_rate_local(202) * 2.d0
  reac_source_local(11,203) = + reac_rate_local(203) 
  reac_source_local(14,203) = - reac_rate_local(203) * 2.d0
  reac_source_local(14,204) = + reac_rate_local(204) 
  reac_source_local(15,204) = - reac_rate_local(204) 
  reac_source_local(29,204) = - reac_rate_local(204) 
  reac_source_local(30,204) = + reac_rate_local(204) 
  reac_source_local(15,205) = - reac_rate_local(205) 
  reac_source_local(21,205) = - reac_rate_local(205) 
  reac_source_local(29,205) = + reac_rate_local(205) 
  reac_source_local(40,205) = + reac_rate_local(205) 
  reac_source_local(01,206) = + reac_rate_local(206) 
  reac_source_local(15,206) = - reac_rate_local(206) 
  reac_source_local(29,206) = + reac_rate_local(206) 
  reac_source_local(40,206) = - reac_rate_local(206) 
  reac_source_local(01,207) = + reac_rate_local(207) 
  reac_source_local(15,207) = - reac_rate_local(207) 
  reac_source_local(40,207) = + reac_rate_local(207) 
  reac_source_local(41,207) = - reac_rate_local(207) 
  reac_source_local(14,208) = + reac_rate_local(208) 
  reac_source_local(15,208) = - reac_rate_local(208) 
  reac_source_local(14,209) = + reac_rate_local(209) 
  reac_source_local(16,209) = - reac_rate_local(209) 
  reac_source_local(14,210) = + reac_rate_local(210) 
  reac_source_local(16,210) = - reac_rate_local(210) 
  reac_source_local(15,211) = + reac_rate_local(211) 
  reac_source_local(16,211) = - reac_rate_local(211) 
  reac_source_local(14,212) = + reac_rate_local(212) 
  reac_source_local(16,212) = - reac_rate_local(212) 
  reac_source_local(15,213) = - reac_rate_local(213) 
  reac_source_local(16,213) = - reac_rate_local(213) 
  reac_source_local(18,213) = + reac_rate_local(213) 
  reac_source_local(53,213) = + reac_rate_local(213) 
  reac_source_local(16,214) = - reac_rate_local(214) 
  reac_source_local(21,214) = - reac_rate_local(214) 
  reac_source_local(29,214) = + reac_rate_local(214) 
  reac_source_local(40,214) = + reac_rate_local(214) 
  reac_source_local(10,215) = + reac_rate_local(215) 
  reac_source_local(16,215) = - reac_rate_local(215) 
  reac_source_local(29,215) = + reac_rate_local(215) 
  reac_source_local(40,215) = - reac_rate_local(215) 
  reac_source_local(21,216) = + reac_rate_local(216) 
  reac_source_local(26,216) = - reac_rate_local(216) 
  reac_source_local(14,217) = - reac_rate_local(217) 
  reac_source_local(26,217) = - reac_rate_local(217) 
  reac_source_local(29,217) = + reac_rate_local(217) 
  reac_source_local(40,217) = + reac_rate_local(217) 
  reac_source_local(21,218) = + reac_rate_local(218) 
  reac_source_local(26,218) = - reac_rate_local(218) 
  reac_source_local(21,219) = + reac_rate_local(219) 
  reac_source_local(26,219) = - reac_rate_local(219) 
  reac_source_local(21,220) = + reac_rate_local(220) 
  reac_source_local(26,220) = - reac_rate_local(220) 
  reac_source_local(21,221) = + reac_rate_local(221) * 2.d0
  reac_source_local(26,221) = - reac_rate_local(221) 
  reac_source_local(30,221) = + reac_rate_local(221) 
  reac_source_local(32,221) = - reac_rate_local(221) 
  reac_source_local(21,222) = + reac_rate_local(222) 
  reac_source_local(26,222) = - reac_rate_local(222) * 2.d0
  reac_source_local(27,222) = + reac_rate_local(222) 
  reac_source_local(21,223) = + reac_rate_local(223) 
  reac_source_local(26,223) = + reac_rate_local(223) 
  reac_source_local(29,223) = - reac_rate_local(223) 
  reac_source_local(32,223) = - reac_rate_local(223) 
  reac_source_local(26,224) = + reac_rate_local(224) 
  reac_source_local(27,224) = - reac_rate_local(224) 
  reac_source_local(21,225) = + reac_rate_local(225) 
  reac_source_local(27,225) = - reac_rate_local(225) 
  reac_source_local(29,225) = - reac_rate_local(225) 
  reac_source_local(30,225) = + reac_rate_local(225) 
  reac_source_local(26,226) = + reac_rate_local(226) 
  reac_source_local(27,226) = - reac_rate_local(226) 
  reac_source_local(26,227) = + reac_rate_local(227) 
  reac_source_local(27,227) = - reac_rate_local(227) 
  reac_source_local(26,228) = + reac_rate_local(228) 
  reac_source_local(27,228) = - reac_rate_local(228) 
  reac_source_local(21,229) = + reac_rate_local(229) * 2.d0
  reac_source_local(27,229) = - reac_rate_local(229) 
  reac_source_local(29,229) = + reac_rate_local(229) 
  reac_source_local(32,229) = - reac_rate_local(229) 
  reac_source_local(21,230) = + reac_rate_local(230) 
  reac_source_local(28,230) = - reac_rate_local(230) 
  reac_source_local(29,230) = - reac_rate_local(230) 
  reac_source_local(31,230) = + reac_rate_local(230) 
  reac_source_local(21,231) = - reac_rate_local(231) 
  reac_source_local(27,231) = + reac_rate_local(231) * 2.d0
  reac_source_local(28,231) = - reac_rate_local(231) 
  reac_source_local(27,232) = + reac_rate_local(232) 
  reac_source_local(28,232) = - reac_rate_local(232) 
  reac_source_local(29,233) = + reac_rate_local(233) 
  reac_source_local(30,233) = - reac_rate_local(233) 
  reac_source_local(29,234) = + reac_rate_local(234) 
  reac_source_local(30,234) = - reac_rate_local(234) 
  reac_source_local(21,235) = - reac_rate_local(235) 
  reac_source_local(26,235) = + reac_rate_local(235) 
  reac_source_local(29,235) = + reac_rate_local(235) 
  reac_source_local(30,235) = - reac_rate_local(235) 
  reac_source_local(21,236) = - reac_rate_local(236) 
  reac_source_local(27,236) = + reac_rate_local(236) 
  reac_source_local(29,236) = + reac_rate_local(236) 
  reac_source_local(30,236) = - reac_rate_local(236) 
  reac_source_local(29,237) = + reac_rate_local(237) 
  reac_source_local(30,237) = - reac_rate_local(237) 
  reac_source_local(21,238) = + reac_rate_local(238) 
  reac_source_local(29,238) = + reac_rate_local(238) * 2.d0
  reac_source_local(30,238) = - reac_rate_local(238) 
  reac_source_local(32,238) = - reac_rate_local(238) 
  reac_source_local(21,239) = + reac_rate_local(239) * 2.d0
  reac_source_local(30,239) = - reac_rate_local(239) 
  reac_source_local(32,239) = - reac_rate_local(239) 
  reac_source_local(14,240) = + reac_rate_local(240) 
  reac_source_local(21,240) = + reac_rate_local(240) 
  reac_source_local(30,240) = - reac_rate_local(240) 
  reac_source_local(40,240) = - reac_rate_local(240) 
  reac_source_local(30,241) = - reac_rate_local(241) 
  reac_source_local(40,241) = + reac_rate_local(241) * 2.d0
  reac_source_local(41,241) = - reac_rate_local(241) 
  reac_source_local(01,242) = + reac_rate_local(242) 
  reac_source_local(21,242) = + reac_rate_local(242) 
  reac_source_local(30,242) = - reac_rate_local(242) 
  reac_source_local(41,242) = - reac_rate_local(242) 
  reac_source_local(30,243) = + reac_rate_local(243) 
  reac_source_local(31,243) = - reac_rate_local(243) 
  reac_source_local(29,244) = + reac_rate_local(244) 
  reac_source_local(31,244) = - reac_rate_local(244) 
  reac_source_local(30,245) = + reac_rate_local(245) 
  reac_source_local(31,245) = - reac_rate_local(245) 
  reac_source_local(21,246) = - reac_rate_local(246) 
  reac_source_local(29,246) = + reac_rate_local(246) * 3.d0
  reac_source_local(31,246) = - reac_rate_local(246) 
  reac_source_local(29,247) = + reac_rate_local(247) 
  reac_source_local(31,247) = - reac_rate_local(247) 
  reac_source_local(26,248) = - reac_rate_local(248) 
  reac_source_local(28,248) = + reac_rate_local(248) 
  reac_source_local(29,248) = + reac_rate_local(248) 
  reac_source_local(31,248) = - reac_rate_local(248) 
  reac_source_local(26,249) = - reac_rate_local(249) 
  reac_source_local(27,249) = + reac_rate_local(249) 
  reac_source_local(30,249) = + reac_rate_local(249) 
  reac_source_local(31,249) = - reac_rate_local(249) 
  reac_source_local(26,250) = - reac_rate_local(250) 
  reac_source_local(29,250) = + reac_rate_local(250) * 3.d0
  reac_source_local(31,250) = - reac_rate_local(250) 
  reac_source_local(29,251) = + reac_rate_local(251) 
  reac_source_local(31,251) = - reac_rate_local(251) 
  reac_source_local(30,252) = + reac_rate_local(252) 
  reac_source_local(31,252) = - reac_rate_local(252) 
  reac_source_local(21,253) = + reac_rate_local(253) * 2.d0
  reac_source_local(31,253) = - reac_rate_local(253) 
  reac_source_local(32,253) = - reac_rate_local(253) 
  reac_source_local(21,254) = + reac_rate_local(254) 
  reac_source_local(29,254) = + reac_rate_local(254) 
  reac_source_local(30,254) = + reac_rate_local(254) 
  reac_source_local(31,254) = - reac_rate_local(254) 
  reac_source_local(32,254) = - reac_rate_local(254) 
  reac_source_local(29,255) = + reac_rate_local(255) 
  reac_source_local(31,255) = - reac_rate_local(255) 
  reac_source_local(30,256) = + reac_rate_local(256) 
  reac_source_local(31,256) = - reac_rate_local(256) 
  reac_source_local(01,257) = + reac_rate_local(257) 
  reac_source_local(14,257) = - reac_rate_local(257) 
  reac_source_local(29,257) = + reac_rate_local(257) 
  reac_source_local(40,257) = - reac_rate_local(257) 
  reac_source_local(14,258) = - reac_rate_local(258) 
  reac_source_local(21,258) = - reac_rate_local(258) 
  reac_source_local(29,258) = + reac_rate_local(258) 
  reac_source_local(40,258) = + reac_rate_local(258) 
  reac_source_local(01,259) = + reac_rate_local(259) 
  reac_source_local(14,259) = - reac_rate_local(259) 
  reac_source_local(29,259) = + reac_rate_local(259) * 2.d0
  reac_source_local(42,259) = - reac_rate_local(259) 
  reac_source_local(14,260) = - reac_rate_local(260) 
  reac_source_local(29,260) = + reac_rate_local(260) 
  reac_source_local(41,260) = + reac_rate_local(260) 
  reac_source_local(42,260) = - reac_rate_local(260) 
  reac_source_local(01,261) = + reac_rate_local(261) 
  reac_source_local(14,261) = - reac_rate_local(261) 
  reac_source_local(21,261) = + reac_rate_local(261) 
  reac_source_local(42,261) = - reac_rate_local(261) 
  reac_source_local(14,262) = - reac_rate_local(262) 
  reac_source_local(40,262) = + reac_rate_local(262) * 2.d0
  reac_source_local(42,262) = - reac_rate_local(262) 
  reac_source_local(01,263) = - reac_rate_local(263) 
  reac_source_local(14,263) = + reac_rate_local(263) 
  reac_source_local(29,263) = - reac_rate_local(263) 
  reac_source_local(40,263) = + reac_rate_local(263) 
  reac_source_local(14,264) = + reac_rate_local(264) 
  reac_source_local(21,264) = + reac_rate_local(264) 
  reac_source_local(29,264) = - reac_rate_local(264) 
  reac_source_local(40,264) = - reac_rate_local(264) 
  reac_source_local(29,265) = - reac_rate_local(265) 
  reac_source_local(40,265) = - reac_rate_local(265) 
  reac_source_local(42,265) = + reac_rate_local(265) 
  reac_source_local(01,266) = + reac_rate_local(266) 
  reac_source_local(21,266) = + reac_rate_local(266) 
  reac_source_local(29,266) = - reac_rate_local(266) 
  reac_source_local(41,266) = - reac_rate_local(266) 
  reac_source_local(29,267) = - reac_rate_local(267) 
  reac_source_local(40,267) = + reac_rate_local(267) * 2.d0
  reac_source_local(41,267) = - reac_rate_local(267) 
  reac_source_local(21,268) = + reac_rate_local(268) 
  reac_source_local(29,268) = - reac_rate_local(268) 
  reac_source_local(40,268) = + reac_rate_local(268) 
  reac_source_local(42,268) = - reac_rate_local(268) 
  reac_source_local(21,269) = + reac_rate_local(269) 
  reac_source_local(29,269) = - reac_rate_local(269) 
  reac_source_local(42,269) = + reac_rate_local(269) 
  reac_source_local(43,269) = - reac_rate_local(269) 
  reac_source_local(01,270) = - reac_rate_local(270) 
  reac_source_local(21,270) = - reac_rate_local(270) 
  reac_source_local(29,270) = + reac_rate_local(270) 
  reac_source_local(41,270) = + reac_rate_local(270) 
  reac_source_local(14,271) = + reac_rate_local(271) 
  reac_source_local(40,271) = - reac_rate_local(271) * 2.d0
  reac_source_local(42,271) = + reac_rate_local(271) 
  reac_source_local(29,272) = + reac_rate_local(272) 
  reac_source_local(40,272) = - reac_rate_local(272) * 2.d0
  reac_source_local(41,272) = + reac_rate_local(272) 
  reac_source_local(01,273) = + reac_rate_local(273) 
  reac_source_local(21,273) = + reac_rate_local(273) 
  reac_source_local(40,273) = - reac_rate_local(273) * 2.d0
  reac_source_local(21,274) = - reac_rate_local(274) 
  reac_source_local(29,274) = + reac_rate_local(274) 
  reac_source_local(40,274) = - reac_rate_local(274) 
  reac_source_local(42,274) = + reac_rate_local(274) 
  reac_source_local(21,275) = + reac_rate_local(275) 
  reac_source_local(32,275) = - reac_rate_local(275) 
  reac_source_local(40,275) = - reac_rate_local(275) 
  reac_source_local(42,275) = + reac_rate_local(275) 
  reac_source_local(01,276) = + reac_rate_local(276) 
  reac_source_local(40,276) = - reac_rate_local(276) 
  reac_source_local(41,276) = - reac_rate_local(276) 
  reac_source_local(42,276) = + reac_rate_local(276) 
  reac_source_local(40,277) = - reac_rate_local(277) 
  reac_source_local(42,277) = + reac_rate_local(277) * 2.d0
  reac_source_local(43,277) = - reac_rate_local(277) 
  reac_source_local(21,278) = - reac_rate_local(278) * 2.d0
  reac_source_local(29,278) = + reac_rate_local(278) 
  reac_source_local(32,278) = + reac_rate_local(278) 
  reac_source_local(21,279) = - reac_rate_local(279) 
  reac_source_local(32,279) = + reac_rate_local(279) 
  reac_source_local(40,279) = + reac_rate_local(279) 
  reac_source_local(42,279) = - reac_rate_local(279) 
  reac_source_local(21,280) = + reac_rate_local(280) 
  reac_source_local(40,280) = + reac_rate_local(280) * 2.d0
  reac_source_local(42,280) = - reac_rate_local(280) * 2.d0
  reac_source_local(40,281) = + reac_rate_local(281) 
  reac_source_local(42,281) = - reac_rate_local(281) * 2.d0
  reac_source_local(43,281) = + reac_rate_local(281) 
  reac_source_local(21,282) = + reac_rate_local(282) 
  reac_source_local(32,282) = - reac_rate_local(282) 
  reac_source_local(42,282) = - reac_rate_local(282) 
  reac_source_local(43,282) = + reac_rate_local(282) 
  reac_source_local(21,283) = + reac_rate_local(283) 
  reac_source_local(40,283) = + reac_rate_local(283) 
  reac_source_local(43,283) = - reac_rate_local(283) 
  reac_source_local(21,284) = - reac_rate_local(284) 
  reac_source_local(32,284) = + reac_rate_local(284) 
  reac_source_local(42,284) = + reac_rate_local(284) 
  reac_source_local(43,284) = - reac_rate_local(284) 
  reac_source_local(21,285) = + reac_rate_local(285) 
  reac_source_local(42,285) = + reac_rate_local(285) * 2.d0
  reac_source_local(43,285) = - reac_rate_local(285) * 2.d0
  reac_source_local(14,286) = - reac_rate_local(286) * 2.d0
  reac_source_local(18,286) = + reac_rate_local(286) 
  reac_source_local(53,286) = + reac_rate_local(286) 
  reac_source_local(14,287) = - reac_rate_local(287) 
  reac_source_local(29,287) = - reac_rate_local(287) 
  reac_source_local(45,287) = + reac_rate_local(287) 
  reac_source_local(53,287) = + reac_rate_local(287) 
  reac_source_local(01,288) = - reac_rate_local(288) 
  reac_source_local(14,288) = + reac_rate_local(288) * 2.d0
  reac_source_local(01,289) = - reac_rate_local(289) 
  reac_source_local(14,289) = + reac_rate_local(289) * 2.d0
  reac_source_local(01,290) = - reac_rate_local(290) 
  reac_source_local(14,290) = + reac_rate_local(290) * 2.d0
  reac_source_local(01,291) = - reac_rate_local(291) 
  reac_source_local(14,291) = + reac_rate_local(291) * 2.d0
  reac_source_local(01,292) = - reac_rate_local(292) 
  reac_source_local(14,292) = + reac_rate_local(292) * 2.d0
  reac_source_local(21,293) = - reac_rate_local(293) 
  reac_source_local(29,293) = + reac_rate_local(293) * 2.d0
  reac_source_local(21,294) = - reac_rate_local(294) 
  reac_source_local(29,294) = + reac_rate_local(294) * 2.d0
  reac_source_local(21,295) = - reac_rate_local(295) 
  reac_source_local(29,295) = + reac_rate_local(295) * 2.d0
  reac_source_local(21,296) = - reac_rate_local(296) 
  reac_source_local(29,296) = + reac_rate_local(296) * 2.d0
  reac_source_local(21,297) = - reac_rate_local(297) 
  reac_source_local(29,297) = + reac_rate_local(297) * 2.d0
  reac_source_local(14,298) = + reac_rate_local(298) 
  reac_source_local(29,298) = + reac_rate_local(298) 
  reac_source_local(40,298) = - reac_rate_local(298) 
  reac_source_local(14,299) = + reac_rate_local(299) 
  reac_source_local(29,299) = + reac_rate_local(299) 
  reac_source_local(40,299) = - reac_rate_local(299) 
  reac_source_local(14,300) = + reac_rate_local(300) 
  reac_source_local(29,300) = + reac_rate_local(300) 
  reac_source_local(40,300) = - reac_rate_local(300) 
  reac_source_local(14,301) = + reac_rate_local(301) 
  reac_source_local(29,301) = + reac_rate_local(301) 
  reac_source_local(40,301) = - reac_rate_local(301) 
  reac_source_local(14,302) = + reac_rate_local(302) 
  reac_source_local(29,302) = + reac_rate_local(302) 
  reac_source_local(40,302) = - reac_rate_local(302) 
  reac_source_local(21,303) = + reac_rate_local(303) 
  reac_source_local(29,303) = + reac_rate_local(303) 
  reac_source_local(32,303) = - reac_rate_local(303) 
  reac_source_local(21,304) = + reac_rate_local(304) 
  reac_source_local(29,304) = + reac_rate_local(304) 
  reac_source_local(32,304) = - reac_rate_local(304) 
  reac_source_local(21,305) = + reac_rate_local(305) 
  reac_source_local(29,305) = + reac_rate_local(305) 
  reac_source_local(32,305) = - reac_rate_local(305) 
  reac_source_local(21,306) = + reac_rate_local(306) 
  reac_source_local(29,306) = + reac_rate_local(306) 
  reac_source_local(32,306) = - reac_rate_local(306) 
  reac_source_local(01,307) = + reac_rate_local(307) 
  reac_source_local(29,307) = + reac_rate_local(307) 
  reac_source_local(41,307) = - reac_rate_local(307) 
  reac_source_local(01,308) = + reac_rate_local(308) 
  reac_source_local(29,308) = + reac_rate_local(308) 
  reac_source_local(41,308) = - reac_rate_local(308) 
  reac_source_local(01,309) = + reac_rate_local(309) 
  reac_source_local(29,309) = + reac_rate_local(309) 
  reac_source_local(41,309) = - reac_rate_local(309) 
  reac_source_local(01,310) = + reac_rate_local(310) 
  reac_source_local(29,310) = + reac_rate_local(310) 
  reac_source_local(41,310) = - reac_rate_local(310) 
  reac_source_local(29,311) = + reac_rate_local(311) 
  reac_source_local(40,311) = + reac_rate_local(311) 
  reac_source_local(42,311) = - reac_rate_local(311) 
  reac_source_local(29,312) = + reac_rate_local(312) 
  reac_source_local(40,312) = + reac_rate_local(312) 
  reac_source_local(42,312) = - reac_rate_local(312) 
  reac_source_local(29,313) = + reac_rate_local(313) 
  reac_source_local(40,313) = + reac_rate_local(313) 
  reac_source_local(42,313) = - reac_rate_local(313) 
  reac_source_local(29,314) = + reac_rate_local(314) 
  reac_source_local(40,314) = + reac_rate_local(314) 
  reac_source_local(42,314) = - reac_rate_local(314) 
  reac_source_local(29,315) = + reac_rate_local(315) 
  reac_source_local(42,315) = + reac_rate_local(315) 
  reac_source_local(43,315) = - reac_rate_local(315) 
  reac_source_local(29,316) = + reac_rate_local(316) 
  reac_source_local(42,316) = + reac_rate_local(316) 
  reac_source_local(43,316) = - reac_rate_local(316) 
  reac_source_local(29,317) = + reac_rate_local(317) 
  reac_source_local(42,317) = + reac_rate_local(317) 
  reac_source_local(43,317) = - reac_rate_local(317) 
  reac_source_local(29,318) = + reac_rate_local(318) 
  reac_source_local(42,318) = + reac_rate_local(318) 
  reac_source_local(43,318) = - reac_rate_local(318) 
  reac_source_local(29,319) = + reac_rate_local(319) 
  reac_source_local(42,319) = + reac_rate_local(319) 
  reac_source_local(43,319) = - reac_rate_local(319) 
  reac_source_local(21,320) = + reac_rate_local(320) 
  reac_source_local(40,320) = + reac_rate_local(320) 
  reac_source_local(43,320) = - reac_rate_local(320) 
  reac_source_local(21,321) = + reac_rate_local(321) 
  reac_source_local(40,321) = + reac_rate_local(321) 
  reac_source_local(43,321) = - reac_rate_local(321) 
  reac_source_local(21,322) = + reac_rate_local(322) 
  reac_source_local(40,322) = + reac_rate_local(322) 
  reac_source_local(43,322) = - reac_rate_local(322) 
  reac_source_local(21,323) = + reac_rate_local(323) 
  reac_source_local(40,323) = + reac_rate_local(323) 
  reac_source_local(43,323) = - reac_rate_local(323) 
  reac_source_local(21,324) = + reac_rate_local(324) 
  reac_source_local(40,324) = + reac_rate_local(324) 
  reac_source_local(43,324) = - reac_rate_local(324) 
  reac_source_local(42,325) = + reac_rate_local(325) 
  reac_source_local(43,325) = + reac_rate_local(325) 
  reac_source_local(44,325) = - reac_rate_local(325) 
  reac_source_local(01,326) = + reac_rate_local(326) 
  reac_source_local(14,326) = - reac_rate_local(326) * 2.d0
  reac_source_local(01,327) = + reac_rate_local(327) 
  reac_source_local(14,327) = - reac_rate_local(327) * 2.d0
  reac_source_local(01,328) = + reac_rate_local(328) 
  reac_source_local(14,328) = - reac_rate_local(328) * 2.d0
  reac_source_local(01,329) = + reac_rate_local(329) 
  reac_source_local(14,329) = - reac_rate_local(329) * 2.d0
  reac_source_local(01,330) = + reac_rate_local(330) 
  reac_source_local(14,330) = - reac_rate_local(330) * 2.d0
  reac_source_local(21,331) = + reac_rate_local(331) 
  reac_source_local(29,331) = - reac_rate_local(331) * 2.d0
  reac_source_local(21,332) = + reac_rate_local(332) 
  reac_source_local(29,332) = - reac_rate_local(332) * 2.d0
  reac_source_local(21,333) = + reac_rate_local(333) 
  reac_source_local(29,333) = - reac_rate_local(333) * 2.d0
  reac_source_local(21,334) = + reac_rate_local(334) 
  reac_source_local(29,334) = - reac_rate_local(334) * 2.d0
  reac_source_local(21,335) = + reac_rate_local(335) 
  reac_source_local(29,335) = - reac_rate_local(335) * 2.d0
  reac_source_local(14,336) = - reac_rate_local(336) 
  reac_source_local(29,336) = - reac_rate_local(336) 
  reac_source_local(40,336) = + reac_rate_local(336) 
  reac_source_local(14,337) = - reac_rate_local(337) 
  reac_source_local(29,337) = - reac_rate_local(337) 
  reac_source_local(40,337) = + reac_rate_local(337) 
  reac_source_local(14,338) = - reac_rate_local(338) 
  reac_source_local(29,338) = - reac_rate_local(338) 
  reac_source_local(40,338) = + reac_rate_local(338) 
  reac_source_local(14,339) = - reac_rate_local(339) 
  reac_source_local(29,339) = - reac_rate_local(339) 
  reac_source_local(40,339) = + reac_rate_local(339) 
  reac_source_local(14,340) = - reac_rate_local(340) 
  reac_source_local(29,340) = - reac_rate_local(340) 
  reac_source_local(40,340) = + reac_rate_local(340) 
  reac_source_local(21,341) = - reac_rate_local(341) 
  reac_source_local(29,341) = - reac_rate_local(341) 
  reac_source_local(32,341) = + reac_rate_local(341) 
  reac_source_local(21,342) = - reac_rate_local(342) 
  reac_source_local(29,342) = - reac_rate_local(342) 
  reac_source_local(32,342) = + reac_rate_local(342) 
  reac_source_local(21,343) = - reac_rate_local(343) 
  reac_source_local(29,343) = - reac_rate_local(343) 
  reac_source_local(32,343) = + reac_rate_local(343) 
  reac_source_local(21,344) = - reac_rate_local(344) 
  reac_source_local(29,344) = - reac_rate_local(344) 
  reac_source_local(32,344) = + reac_rate_local(344) 
  reac_source_local(21,345) = - reac_rate_local(345) 
  reac_source_local(29,345) = - reac_rate_local(345) 
  reac_source_local(32,345) = + reac_rate_local(345) 
  reac_source_local(01,346) = - reac_rate_local(346) 
  reac_source_local(29,346) = - reac_rate_local(346) 
  reac_source_local(41,346) = + reac_rate_local(346) 
  reac_source_local(29,347) = - reac_rate_local(347) 
  reac_source_local(40,347) = - reac_rate_local(347) 
  reac_source_local(42,347) = + reac_rate_local(347) 
  reac_source_local(29,348) = - reac_rate_local(348) 
  reac_source_local(40,348) = - reac_rate_local(348) 
  reac_source_local(42,348) = + reac_rate_local(348) 
  reac_source_local(29,349) = - reac_rate_local(349) 
  reac_source_local(40,349) = - reac_rate_local(349) 
  reac_source_local(42,349) = + reac_rate_local(349) 
  reac_source_local(29,350) = - reac_rate_local(350) 
  reac_source_local(42,350) = - reac_rate_local(350) 
  reac_source_local(43,350) = + reac_rate_local(350) 
  reac_source_local(29,351) = - reac_rate_local(351) 
  reac_source_local(42,351) = - reac_rate_local(351) 
  reac_source_local(43,351) = + reac_rate_local(351) 
  reac_source_local(29,352) = - reac_rate_local(352) 
  reac_source_local(42,352) = - reac_rate_local(352) 
  reac_source_local(43,352) = + reac_rate_local(352) 
  reac_source_local(29,353) = - reac_rate_local(353) 
  reac_source_local(42,353) = - reac_rate_local(353) 
  reac_source_local(43,353) = + reac_rate_local(353) 
  reac_source_local(29,354) = - reac_rate_local(354) 
  reac_source_local(42,354) = - reac_rate_local(354) 
  reac_source_local(43,354) = + reac_rate_local(354) 
  reac_source_local(42,355) = - reac_rate_local(355) 
  reac_source_local(43,355) = - reac_rate_local(355) 
  reac_source_local(44,355) = + reac_rate_local(355) 
  reac_source_local(14,356) = + reac_rate_local(356) 
  reac_source_local(17,356) = - reac_rate_local(356) 
  reac_source_local(29,356) = - reac_rate_local(356) 
  reac_source_local(33,356) = + reac_rate_local(356) 
  reac_source_local(14,357) = + reac_rate_local(357) 
  reac_source_local(17,357) = - reac_rate_local(357) 
  reac_source_local(21,357) = - reac_rate_local(357) 
  reac_source_local(34,357) = + reac_rate_local(357) 
  reac_source_local(17,358) = - reac_rate_local(358) 
  reac_source_local(21,358) = - reac_rate_local(358) 
  reac_source_local(29,358) = + reac_rate_local(358) 
  reac_source_local(45,358) = + reac_rate_local(358) 
  reac_source_local(17,359) = - reac_rate_local(359) 
  reac_source_local(21,359) = - reac_rate_local(359) 
  reac_source_local(33,359) = + reac_rate_local(359) 
  reac_source_local(40,359) = + reac_rate_local(359) 
  reac_source_local(17,360) = - reac_rate_local(360) 
  reac_source_local(21,360) = + reac_rate_local(360) 
  reac_source_local(32,360) = - reac_rate_local(360) 
  reac_source_local(45,360) = + reac_rate_local(360) 
  reac_source_local(14,361) = + reac_rate_local(361) 
  reac_source_local(17,361) = - reac_rate_local(361) 
  reac_source_local(40,361) = - reac_rate_local(361) 
  reac_source_local(45,361) = + reac_rate_local(361) 
  reac_source_local(17,362) = - reac_rate_local(362) 
  reac_source_local(18,362) = + reac_rate_local(362) 
  reac_source_local(29,362) = + reac_rate_local(362) 
  reac_source_local(40,362) = - reac_rate_local(362) 
  reac_source_local(01,363) = + reac_rate_local(363) 
  reac_source_local(17,363) = - reac_rate_local(363) 
  reac_source_local(33,363) = + reac_rate_local(363) 
  reac_source_local(40,363) = - reac_rate_local(363) 
  reac_source_local(01,364) = + reac_rate_local(364) 
  reac_source_local(17,364) = - reac_rate_local(364) 
  reac_source_local(41,364) = - reac_rate_local(364) 
  reac_source_local(45,364) = + reac_rate_local(364) 
  reac_source_local(01,365) = - reac_rate_local(365) 
  reac_source_local(14,365) = + reac_rate_local(365) 
  reac_source_local(33,365) = - reac_rate_local(365) 
  reac_source_local(45,365) = + reac_rate_local(365) 
  reac_source_local(21,366) = - reac_rate_local(366) 
  reac_source_local(29,366) = + reac_rate_local(366) 
  reac_source_local(33,366) = - reac_rate_local(366) 
  reac_source_local(34,366) = + reac_rate_local(366) 
  reac_source_local(21,367) = + reac_rate_local(367) 
  reac_source_local(32,367) = - reac_rate_local(367) 
  reac_source_local(33,367) = - reac_rate_local(367) 
  reac_source_local(34,367) = + reac_rate_local(367) 
  reac_source_local(29,368) = + reac_rate_local(368) 
  reac_source_local(33,368) = - reac_rate_local(368) 
  reac_source_local(40,368) = - reac_rate_local(368) 
  reac_source_local(45,368) = + reac_rate_local(368) 
  reac_source_local(14,369) = + reac_rate_local(369) 
  reac_source_local(33,369) = - reac_rate_local(369) 
  reac_source_local(34,369) = + reac_rate_local(369) 
  reac_source_local(40,369) = - reac_rate_local(369) 
  reac_source_local(15,370) = - reac_rate_local(370) 
  reac_source_local(17,370) = + reac_rate_local(370) 
  reac_source_local(29,370) = + reac_rate_local(370) 
  reac_source_local(33,370) = - reac_rate_local(370) 
  reac_source_local(33,371) = - reac_rate_local(371) 
  reac_source_local(40,371) = + reac_rate_local(371) 
  reac_source_local(41,371) = - reac_rate_local(371) 
  reac_source_local(45,371) = + reac_rate_local(371) 
  reac_source_local(29,372) = + reac_rate_local(372) 
  reac_source_local(33,372) = - reac_rate_local(372) 
  reac_source_local(41,372) = - reac_rate_local(372) 
  reac_source_local(46,372) = + reac_rate_local(372) 
  reac_source_local(01,373) = + reac_rate_local(373) 
  reac_source_local(33,373) = - reac_rate_local(373) 
  reac_source_local(34,373) = + reac_rate_local(373) 
  reac_source_local(41,373) = - reac_rate_local(373) 
  reac_source_local(29,374) = + reac_rate_local(374) 
  reac_source_local(33,374) = - reac_rate_local(374) 
  reac_source_local(42,374) = - reac_rate_local(374) 
  reac_source_local(47,374) = + reac_rate_local(374) 
  reac_source_local(01,375) = + reac_rate_local(375) 
  reac_source_local(18,375) = - reac_rate_local(375) 
  reac_source_local(21,375) = - reac_rate_local(375) 
  reac_source_local(34,375) = + reac_rate_local(375) 
  reac_source_local(14,376) = + reac_rate_local(376) 
  reac_source_local(18,376) = - reac_rate_local(376) 
  reac_source_local(29,376) = - reac_rate_local(376) 
  reac_source_local(45,376) = + reac_rate_local(376) 
  reac_source_local(01,377) = + reac_rate_local(377) 
  reac_source_local(18,377) = - reac_rate_local(377) 
  reac_source_local(29,377) = + reac_rate_local(377) 
  reac_source_local(32,377) = - reac_rate_local(377) 
  reac_source_local(34,377) = + reac_rate_local(377) 
  reac_source_local(01,378) = + reac_rate_local(378) 
  reac_source_local(14,378) = - reac_rate_local(378) 
  reac_source_local(17,378) = + reac_rate_local(378) 
  reac_source_local(18,378) = - reac_rate_local(378) 
  reac_source_local(01,379) = + reac_rate_local(379) 
  reac_source_local(18,379) = - reac_rate_local(379) 
  reac_source_local(40,379) = - reac_rate_local(379) 
  reac_source_local(45,379) = + reac_rate_local(379) 
  reac_source_local(01,380) = + reac_rate_local(380) 
  reac_source_local(18,380) = - reac_rate_local(380) 
  reac_source_local(41,380) = - reac_rate_local(380) 
  reac_source_local(46,380) = + reac_rate_local(380) 
  reac_source_local(01,381) = + reac_rate_local(381) 
  reac_source_local(14,381) = + reac_rate_local(381) 
  reac_source_local(18,381) = - reac_rate_local(381) 
  reac_source_local(41,381) = - reac_rate_local(381) 
  reac_source_local(45,381) = + reac_rate_local(381) 
  reac_source_local(01,382) = - reac_rate_local(382) 
  reac_source_local(34,382) = - reac_rate_local(382) 
  reac_source_local(40,382) = + reac_rate_local(382) 
  reac_source_local(45,382) = + reac_rate_local(382) 
  reac_source_local(14,383) = - reac_rate_local(383) 
  reac_source_local(29,383) = + reac_rate_local(383) 
  reac_source_local(34,383) = - reac_rate_local(383) 
  reac_source_local(45,383) = + reac_rate_local(383) 
  reac_source_local(21,384) = + reac_rate_local(384) 
  reac_source_local(34,384) = - reac_rate_local(384) 
  reac_source_local(40,384) = - reac_rate_local(384) 
  reac_source_local(45,384) = + reac_rate_local(384) 
  reac_source_local(32,385) = + reac_rate_local(385) 
  reac_source_local(34,385) = - reac_rate_local(385) 
  reac_source_local(42,385) = - reac_rate_local(385) 
  reac_source_local(45,385) = + reac_rate_local(385) 
  reac_source_local(21,386) = + reac_rate_local(386) 
  reac_source_local(34,386) = - reac_rate_local(386) 
  reac_source_local(42,386) = - reac_rate_local(386) 
  reac_source_local(47,386) = + reac_rate_local(386) 
  reac_source_local(01,387) = + reac_rate_local(387) 
  reac_source_local(14,387) = + reac_rate_local(387) 
  reac_source_local(19,387) = - reac_rate_local(387) 
  reac_source_local(21,387) = - reac_rate_local(387) 
  reac_source_local(34,387) = + reac_rate_local(387) 
  reac_source_local(01,388) = + reac_rate_local(388) 
  reac_source_local(19,388) = - reac_rate_local(388) 
  reac_source_local(21,388) = - reac_rate_local(388) 
  reac_source_local(47,388) = + reac_rate_local(388) 
  reac_source_local(01,389) = + reac_rate_local(389) 
  reac_source_local(14,389) = - reac_rate_local(389) 
  reac_source_local(18,389) = + reac_rate_local(389) 
  reac_source_local(19,389) = - reac_rate_local(389) 
  reac_source_local(01,390) = + reac_rate_local(390) 
  reac_source_local(14,390) = + reac_rate_local(390) 
  reac_source_local(19,390) = - reac_rate_local(390) 
  reac_source_local(40,390) = - reac_rate_local(390) 
  reac_source_local(45,390) = + reac_rate_local(390) 
  reac_source_local(01,391) = + reac_rate_local(391) 
  reac_source_local(19,391) = - reac_rate_local(391) 
  reac_source_local(40,391) = - reac_rate_local(391) 
  reac_source_local(46,391) = + reac_rate_local(391) 
  reac_source_local(40,392) = - reac_rate_local(392) 
  reac_source_local(42,392) = + reac_rate_local(392) 
  reac_source_local(45,392) = + reac_rate_local(392) 
  reac_source_local(47,392) = - reac_rate_local(392) 
  reac_source_local(40,393) = - reac_rate_local(393) 
  reac_source_local(41,393) = + reac_rate_local(393) 
  reac_source_local(45,393) = + reac_rate_local(393) 
  reac_source_local(46,393) = - reac_rate_local(393) 
  reac_source_local(01,394) = + reac_rate_local(394) 
  reac_source_local(18,394) = + reac_rate_local(394) 
  reac_source_local(20,394) = - reac_rate_local(394) 
  reac_source_local(01,395) = + reac_rate_local(395) * 2.d0
  reac_source_local(20,395) = - reac_rate_local(395) 
  reac_source_local(21,395) = - reac_rate_local(395) 
  reac_source_local(34,395) = + reac_rate_local(395) 
  reac_source_local(01,396) = + reac_rate_local(396) * 2.d0
  reac_source_local(20,396) = - reac_rate_local(396) 
  reac_source_local(29,396) = - reac_rate_local(396) 
  reac_source_local(33,396) = + reac_rate_local(396) 
  reac_source_local(01,397) = + reac_rate_local(397) * 2.d0
  reac_source_local(14,397) = - reac_rate_local(397) 
  reac_source_local(17,397) = + reac_rate_local(397) 
  reac_source_local(20,397) = - reac_rate_local(397) 
  reac_source_local(01,398) = + reac_rate_local(398) * 2.d0
  reac_source_local(20,398) = - reac_rate_local(398) 
  reac_source_local(40,398) = - reac_rate_local(398) 
  reac_source_local(45,398) = + reac_rate_local(398) 
  reac_source_local(01,399) = - reac_rate_local(399) 
  reac_source_local(21,399) = + reac_rate_local(399) 
  reac_source_local(35,399) = - reac_rate_local(399) 
  reac_source_local(52,399) = + reac_rate_local(399) 
  reac_source_local(21,400) = + reac_rate_local(400) 
  reac_source_local(34,400) = + reac_rate_local(400) 
  reac_source_local(35,400) = - reac_rate_local(400) 
  reac_source_local(21,401) = + reac_rate_local(401) * 2.d0
  reac_source_local(26,401) = - reac_rate_local(401) 
  reac_source_local(34,401) = + reac_rate_local(401) 
  reac_source_local(35,401) = - reac_rate_local(401) 
  reac_source_local(21,402) = + reac_rate_local(402) * 2.d0
  reac_source_local(27,402) = - reac_rate_local(402) 
  reac_source_local(34,402) = + reac_rate_local(402) 
  reac_source_local(35,402) = - reac_rate_local(402) 
  reac_source_local(29,403) = - reac_rate_local(403) 
  reac_source_local(32,403) = + reac_rate_local(403) 
  reac_source_local(34,403) = + reac_rate_local(403) 
  reac_source_local(35,403) = - reac_rate_local(403) 
  reac_source_local(21,404) = + reac_rate_local(404) * 2.d0
  reac_source_local(35,404) = - reac_rate_local(404) 
  reac_source_local(40,404) = - reac_rate_local(404) 
  reac_source_local(45,404) = + reac_rate_local(404) 
  reac_source_local(01,405) = + reac_rate_local(405) 
  reac_source_local(34,405) = + reac_rate_local(405) 
  reac_source_local(52,405) = - reac_rate_local(405) 
  reac_source_local(01,406) = + reac_rate_local(406) 
  reac_source_local(21,406) = - reac_rate_local(406) 
  reac_source_local(35,406) = + reac_rate_local(406) 
  reac_source_local(52,406) = - reac_rate_local(406) 
  reac_source_local(01,407) = - reac_rate_local(407) 
  reac_source_local(17,407) = - reac_rate_local(407) 
  reac_source_local(19,407) = + reac_rate_local(407) 
  reac_source_local(17,408) = - reac_rate_local(408) 
  reac_source_local(29,408) = - reac_rate_local(408) 
  reac_source_local(45,408) = + reac_rate_local(408) 
  reac_source_local(14,409) = - reac_rate_local(409) 
  reac_source_local(17,409) = - reac_rate_local(409) 
  reac_source_local(18,409) = + reac_rate_local(409) 
  reac_source_local(01,410) = - reac_rate_local(410) 
  reac_source_local(14,410) = + reac_rate_local(410) 
  reac_source_local(33,410) = - reac_rate_local(410) 
  reac_source_local(45,410) = + reac_rate_local(410) 
  reac_source_local(29,411) = - reac_rate_local(411) 
  reac_source_local(33,411) = - reac_rate_local(411) 
  reac_source_local(34,411) = + reac_rate_local(411) 
  reac_source_local(14,412) = - reac_rate_local(412) 
  reac_source_local(33,412) = - reac_rate_local(412) 
  reac_source_local(45,412) = + reac_rate_local(412) 
  reac_source_local(01,413) = - reac_rate_local(413) 
  reac_source_local(18,413) = - reac_rate_local(413) 
  reac_source_local(20,413) = + reac_rate_local(413) 
  reac_source_local(14,414) = - reac_rate_local(414) 
  reac_source_local(18,414) = - reac_rate_local(414) 
  reac_source_local(19,414) = + reac_rate_local(414) 
  reac_source_local(21,415) = - reac_rate_local(415) 
  reac_source_local(34,415) = - reac_rate_local(415) 
  reac_source_local(35,415) = + reac_rate_local(415) 
  reac_source_local(01,416) = - reac_rate_local(416) 
  reac_source_local(34,416) = - reac_rate_local(416) 
  reac_source_local(52,416) = + reac_rate_local(416) 
  reac_source_local(26,417) = - reac_rate_local(417) 
  reac_source_local(29,417) = + reac_rate_local(417) 
  reac_source_local(36,417) = - reac_rate_local(417) 
  reac_source_local(37,417) = + reac_rate_local(417) 
  reac_source_local(29,418) = + reac_rate_local(418) 
  reac_source_local(32,418) = - reac_rate_local(418) 
  reac_source_local(36,418) = - reac_rate_local(418) 
  reac_source_local(38,418) = + reac_rate_local(418) 
  reac_source_local(29,419) = + reac_rate_local(419) 
  reac_source_local(36,419) = - reac_rate_local(419) 
  reac_source_local(42,419) = - reac_rate_local(419) 
  reac_source_local(50,419) = + reac_rate_local(419) 
  reac_source_local(36,420) = - reac_rate_local(420) 
  reac_source_local(40,420) = + reac_rate_local(420) 
  reac_source_local(41,420) = - reac_rate_local(420) 
  reac_source_local(48,420) = + reac_rate_local(420) 
  reac_source_local(29,421) = + reac_rate_local(421) 
  reac_source_local(36,421) = - reac_rate_local(421) 
  reac_source_local(41,421) = - reac_rate_local(421) 
  reac_source_local(49,421) = + reac_rate_local(421) 
  reac_source_local(21,422) = + reac_rate_local(422) 
  reac_source_local(29,422) = - reac_rate_local(422) 
  reac_source_local(36,422) = + reac_rate_local(422) 
  reac_source_local(37,422) = - reac_rate_local(422) 
  reac_source_local(21,423) = + reac_rate_local(423) 
  reac_source_local(32,423) = - reac_rate_local(423) 
  reac_source_local(37,423) = - reac_rate_local(423) 
  reac_source_local(38,423) = + reac_rate_local(423) 
  reac_source_local(21,424) = + reac_rate_local(424) 
  reac_source_local(37,424) = - reac_rate_local(424) 
  reac_source_local(42,424) = - reac_rate_local(424) 
  reac_source_local(50,424) = + reac_rate_local(424) 
  reac_source_local(21,425) = + reac_rate_local(425) 
  reac_source_local(37,425) = - reac_rate_local(425) 
  reac_source_local(43,425) = - reac_rate_local(425) 
  reac_source_local(51,425) = + reac_rate_local(425) 
  reac_source_local(21,426) = + reac_rate_local(426) 
  reac_source_local(29,426) = - reac_rate_local(426) 
  reac_source_local(37,426) = + reac_rate_local(426) 
  reac_source_local(38,426) = - reac_rate_local(426) 
  reac_source_local(29,427) = + reac_rate_local(427) 
  reac_source_local(38,427) = - reac_rate_local(427) 
  reac_source_local(40,427) = - reac_rate_local(427) 
  reac_source_local(51,427) = + reac_rate_local(427) 
  reac_source_local(21,428) = + reac_rate_local(428) 
  reac_source_local(38,428) = - reac_rate_local(428) 
  reac_source_local(40,428) = - reac_rate_local(428) 
  reac_source_local(50,428) = + reac_rate_local(428) 
  reac_source_local(32,429) = + reac_rate_local(429) 
  reac_source_local(38,429) = - reac_rate_local(429) 
  reac_source_local(42,429) = - reac_rate_local(429) 
  reac_source_local(50,429) = + reac_rate_local(429) 
  reac_source_local(21,430) = + reac_rate_local(430) 
  reac_source_local(38,430) = - reac_rate_local(430) 
  reac_source_local(42,430) = - reac_rate_local(430) 
  reac_source_local(51,430) = + reac_rate_local(430) 
  reac_source_local(32,431) = + reac_rate_local(431) 
  reac_source_local(38,431) = - reac_rate_local(431) 
  reac_source_local(43,431) = - reac_rate_local(431) 
  reac_source_local(51,431) = + reac_rate_local(431) 
  reac_source_local(21,432) = - reac_rate_local(432) 
  reac_source_local(37,432) = + reac_rate_local(432) 
  reac_source_local(40,432) = + reac_rate_local(432) 
  reac_source_local(48,432) = - reac_rate_local(432) 
  reac_source_local(40,433) = + reac_rate_local(433) 
  reac_source_local(42,433) = - reac_rate_local(433) 
  reac_source_local(48,433) = - reac_rate_local(433) 
  reac_source_local(50,433) = + reac_rate_local(433) 
  reac_source_local(01,434) = + reac_rate_local(434) 
  reac_source_local(41,434) = - reac_rate_local(434) 
  reac_source_local(48,434) = - reac_rate_local(434) 
  reac_source_local(50,434) = + reac_rate_local(434) 
  reac_source_local(21,435) = + reac_rate_local(435) 
  reac_source_local(32,435) = - reac_rate_local(435) 
  reac_source_local(50,435) = - reac_rate_local(435) 
  reac_source_local(51,435) = + reac_rate_local(435) 
  reac_source_local(40,436) = + reac_rate_local(436) 
  reac_source_local(42,436) = - reac_rate_local(436) 
  reac_source_local(50,436) = - reac_rate_local(436) 
  reac_source_local(51,436) = + reac_rate_local(436) 
  reac_source_local(42,437) = + reac_rate_local(437) 
  reac_source_local(43,437) = - reac_rate_local(437) 
  reac_source_local(50,437) = - reac_rate_local(437) 
  reac_source_local(51,437) = + reac_rate_local(437) 
  reac_source_local(42,438) = + reac_rate_local(438) * 2.d0
  reac_source_local(44,438) = - reac_rate_local(438) 
  reac_source_local(50,438) = - reac_rate_local(438) 
  reac_source_local(51,438) = + reac_rate_local(438) 
  reac_source_local(40,439) = - reac_rate_local(439) 
  reac_source_local(42,439) = + reac_rate_local(439) 
  reac_source_local(50,439) = + reac_rate_local(439) 
  reac_source_local(51,439) = - reac_rate_local(439) 
  reac_source_local(21,440) = + reac_rate_local(440) 
  reac_source_local(37,440) = + reac_rate_local(440) 
  reac_source_local(39,440) = - reac_rate_local(440) 
  reac_source_local(21,441) = + reac_rate_local(441) 
  reac_source_local(37,441) = + reac_rate_local(441) 
  reac_source_local(39,441) = - reac_rate_local(441) 
  reac_source_local(21,442) = + reac_rate_local(442) 
  reac_source_local(29,442) = - reac_rate_local(442) 
  reac_source_local(38,442) = + reac_rate_local(442) 
  reac_source_local(39,442) = - reac_rate_local(442) 
  reac_source_local(21,443) = + reac_rate_local(443) * 2.d0
  reac_source_local(29,443) = - reac_rate_local(443) 
  reac_source_local(36,443) = + reac_rate_local(443) 
  reac_source_local(39,443) = - reac_rate_local(443) 
  reac_source_local(21,444) = + reac_rate_local(444) * 2.d0
  reac_source_local(26,444) = - reac_rate_local(444) 
  reac_source_local(37,444) = + reac_rate_local(444) 
  reac_source_local(39,444) = - reac_rate_local(444) 
  reac_source_local(21,445) = + reac_rate_local(445) * 2.d0
  reac_source_local(27,445) = - reac_rate_local(445) 
  reac_source_local(37,445) = + reac_rate_local(445) 
  reac_source_local(39,445) = - reac_rate_local(445) 
  reac_source_local(21,446) = + reac_rate_local(446) 
  reac_source_local(39,446) = - reac_rate_local(446) 
  reac_source_local(40,446) = - reac_rate_local(446) 
  reac_source_local(51,446) = + reac_rate_local(446) 
  reac_source_local(21,447) = - reac_rate_local(447) 
  reac_source_local(36,447) = - reac_rate_local(447) 
  reac_source_local(38,447) = + reac_rate_local(447) 
  reac_source_local(36,448) = - reac_rate_local(448) 
  reac_source_local(40,448) = - reac_rate_local(448) 
  reac_source_local(50,448) = + reac_rate_local(448) 
  reac_source_local(21,449) = - reac_rate_local(449) 
  reac_source_local(37,449) = - reac_rate_local(449) 
  reac_source_local(39,449) = + reac_rate_local(449) 
  reac_source_local(14,450) = + reac_rate_local(450) 
  reac_source_local(17,450) = - reac_rate_local(450) 
  reac_source_local(29,450) = + reac_rate_local(450) 
  reac_source_local(36,450) = - reac_rate_local(450) 
  reac_source_local(01,451) = + reac_rate_local(451) 
  reac_source_local(18,451) = - reac_rate_local(451) 
  reac_source_local(29,451) = + reac_rate_local(451) 
  reac_source_local(36,451) = - reac_rate_local(451) 
  reac_source_local(29,452) = + reac_rate_local(452) * 2.d0
  reac_source_local(33,452) = - reac_rate_local(452) 
  reac_source_local(36,452) = - reac_rate_local(452) 
  reac_source_local(21,453) = + reac_rate_local(453) 
  reac_source_local(29,453) = + reac_rate_local(453) 
  reac_source_local(34,453) = - reac_rate_local(453) 
  reac_source_local(36,453) = - reac_rate_local(453) 
  reac_source_local(29,454) = + reac_rate_local(454) 
  reac_source_local(36,454) = - reac_rate_local(454) 
  reac_source_local(40,454) = + reac_rate_local(454) 
  reac_source_local(45,454) = - reac_rate_local(454) 
  reac_source_local(29,455) = + reac_rate_local(455) 
  reac_source_local(36,455) = - reac_rate_local(455) 
  reac_source_local(41,455) = + reac_rate_local(455) 
  reac_source_local(46,455) = - reac_rate_local(455) 
  reac_source_local(29,456) = + reac_rate_local(456) 
  reac_source_local(36,456) = - reac_rate_local(456) 
  reac_source_local(42,456) = + reac_rate_local(456) 
  reac_source_local(47,456) = - reac_rate_local(456) 
  reac_source_local(14,457) = + reac_rate_local(457) 
  reac_source_local(17,457) = - reac_rate_local(457) 
  reac_source_local(21,457) = + reac_rate_local(457) 
  reac_source_local(37,457) = - reac_rate_local(457) 
  reac_source_local(01,458) = + reac_rate_local(458) 
  reac_source_local(18,458) = - reac_rate_local(458) 
  reac_source_local(21,458) = + reac_rate_local(458) 
  reac_source_local(37,458) = - reac_rate_local(458) 
  reac_source_local(21,459) = + reac_rate_local(459) 
  reac_source_local(29,459) = + reac_rate_local(459) 
  reac_source_local(33,459) = - reac_rate_local(459) 
  reac_source_local(37,459) = - reac_rate_local(459) 
  reac_source_local(21,460) = + reac_rate_local(460) * 2.d0
  reac_source_local(34,460) = - reac_rate_local(460) 
  reac_source_local(37,460) = - reac_rate_local(460) 
  reac_source_local(21,461) = + reac_rate_local(461) 
  reac_source_local(37,461) = - reac_rate_local(461) 
  reac_source_local(40,461) = + reac_rate_local(461) 
  reac_source_local(45,461) = - reac_rate_local(461) 
  reac_source_local(21,462) = + reac_rate_local(462) 
  reac_source_local(37,462) = - reac_rate_local(462) 
  reac_source_local(41,462) = + reac_rate_local(462) 
  reac_source_local(46,462) = - reac_rate_local(462) 
  reac_source_local(21,463) = + reac_rate_local(463) 
  reac_source_local(37,463) = - reac_rate_local(463) 
  reac_source_local(42,463) = + reac_rate_local(463) 
  reac_source_local(47,463) = - reac_rate_local(463) 
  reac_source_local(14,464) = + reac_rate_local(464) 
  reac_source_local(17,464) = - reac_rate_local(464) 
  reac_source_local(32,464) = + reac_rate_local(464) 
  reac_source_local(38,464) = - reac_rate_local(464) 
  reac_source_local(01,465) = + reac_rate_local(465) 
  reac_source_local(18,465) = - reac_rate_local(465) 
  reac_source_local(32,465) = + reac_rate_local(465) 
  reac_source_local(38,465) = - reac_rate_local(465) 
  reac_source_local(29,466) = + reac_rate_local(466) 
  reac_source_local(32,466) = + reac_rate_local(466) 
  reac_source_local(33,466) = - reac_rate_local(466) 
  reac_source_local(38,466) = - reac_rate_local(466) 
  reac_source_local(21,467) = + reac_rate_local(467) 
  reac_source_local(32,467) = + reac_rate_local(467) 
  reac_source_local(34,467) = - reac_rate_local(467) 
  reac_source_local(38,467) = - reac_rate_local(467) 
  reac_source_local(32,468) = + reac_rate_local(468) 
  reac_source_local(38,468) = - reac_rate_local(468) 
  reac_source_local(40,468) = + reac_rate_local(468) 
  reac_source_local(45,468) = - reac_rate_local(468) 
  reac_source_local(32,469) = + reac_rate_local(469) 
  reac_source_local(38,469) = - reac_rate_local(469) 
  reac_source_local(41,469) = + reac_rate_local(469) 
  reac_source_local(46,469) = - reac_rate_local(469) 
  reac_source_local(32,470) = + reac_rate_local(470) 
  reac_source_local(38,470) = - reac_rate_local(470) 
  reac_source_local(42,470) = + reac_rate_local(470) 
  reac_source_local(47,470) = - reac_rate_local(470) 
  reac_source_local(14,471) = + reac_rate_local(471) 
  reac_source_local(17,471) = - reac_rate_local(471) 
  reac_source_local(40,471) = + reac_rate_local(471) 
  reac_source_local(48,471) = - reac_rate_local(471) 
  reac_source_local(01,472) = + reac_rate_local(472) 
  reac_source_local(18,472) = - reac_rate_local(472) 
  reac_source_local(40,472) = + reac_rate_local(472) 
  reac_source_local(48,472) = - reac_rate_local(472) 
  reac_source_local(29,473) = + reac_rate_local(473) 
  reac_source_local(33,473) = - reac_rate_local(473) 
  reac_source_local(40,473) = + reac_rate_local(473) 
  reac_source_local(48,473) = - reac_rate_local(473) 
  reac_source_local(21,474) = + reac_rate_local(474) 
  reac_source_local(34,474) = - reac_rate_local(474) 
  reac_source_local(40,474) = + reac_rate_local(474) 
  reac_source_local(48,474) = - reac_rate_local(474) 
  reac_source_local(40,475) = + reac_rate_local(475) * 2.d0
  reac_source_local(45,475) = - reac_rate_local(475) 
  reac_source_local(48,475) = - reac_rate_local(475) 
  reac_source_local(40,476) = + reac_rate_local(476) 
  reac_source_local(41,476) = + reac_rate_local(476) 
  reac_source_local(46,476) = - reac_rate_local(476) 
  reac_source_local(48,476) = - reac_rate_local(476) 
  reac_source_local(40,477) = + reac_rate_local(477) 
  reac_source_local(42,477) = + reac_rate_local(477) 
  reac_source_local(47,477) = - reac_rate_local(477) 
  reac_source_local(48,477) = - reac_rate_local(477) 
  reac_source_local(14,478) = + reac_rate_local(478) 
  reac_source_local(17,478) = - reac_rate_local(478) 
  reac_source_local(41,478) = + reac_rate_local(478) 
  reac_source_local(49,478) = - reac_rate_local(478) 
  reac_source_local(01,479) = + reac_rate_local(479) 
  reac_source_local(18,479) = - reac_rate_local(479) 
  reac_source_local(41,479) = + reac_rate_local(479) 
  reac_source_local(49,479) = - reac_rate_local(479) 
  reac_source_local(29,480) = + reac_rate_local(480) 
  reac_source_local(33,480) = - reac_rate_local(480) 
  reac_source_local(41,480) = + reac_rate_local(480) 
  reac_source_local(49,480) = - reac_rate_local(480) 
  reac_source_local(21,481) = + reac_rate_local(481) 
  reac_source_local(34,481) = - reac_rate_local(481) 
  reac_source_local(41,481) = + reac_rate_local(481) 
  reac_source_local(49,481) = - reac_rate_local(481) 
  reac_source_local(40,482) = + reac_rate_local(482) 
  reac_source_local(41,482) = + reac_rate_local(482) 
  reac_source_local(45,482) = - reac_rate_local(482) 
  reac_source_local(49,482) = - reac_rate_local(482) 
  reac_source_local(41,483) = + reac_rate_local(483) * 2.d0
  reac_source_local(46,483) = - reac_rate_local(483) 
  reac_source_local(49,483) = - reac_rate_local(483) 
  reac_source_local(41,484) = + reac_rate_local(484) 
  reac_source_local(42,484) = + reac_rate_local(484) 
  reac_source_local(47,484) = - reac_rate_local(484) 
  reac_source_local(49,484) = - reac_rate_local(484) 
  reac_source_local(14,485) = + reac_rate_local(485) 
  reac_source_local(17,485) = - reac_rate_local(485) 
  reac_source_local(42,485) = + reac_rate_local(485) 
  reac_source_local(50,485) = - reac_rate_local(485) 
  reac_source_local(01,486) = + reac_rate_local(486) 
  reac_source_local(18,486) = - reac_rate_local(486) 
  reac_source_local(42,486) = + reac_rate_local(486) 
  reac_source_local(50,486) = - reac_rate_local(486) 
  reac_source_local(29,487) = + reac_rate_local(487) 
  reac_source_local(33,487) = - reac_rate_local(487) 
  reac_source_local(42,487) = + reac_rate_local(487) 
  reac_source_local(50,487) = - reac_rate_local(487) 
  reac_source_local(21,488) = + reac_rate_local(488) 
  reac_source_local(34,488) = - reac_rate_local(488) 
  reac_source_local(42,488) = + reac_rate_local(488) 
  reac_source_local(50,488) = - reac_rate_local(488) 
  reac_source_local(40,489) = + reac_rate_local(489) 
  reac_source_local(42,489) = + reac_rate_local(489) 
  reac_source_local(45,489) = - reac_rate_local(489) 
  reac_source_local(50,489) = - reac_rate_local(489) 
  reac_source_local(41,490) = + reac_rate_local(490) 
  reac_source_local(42,490) = + reac_rate_local(490) 
  reac_source_local(46,490) = - reac_rate_local(490) 
  reac_source_local(50,490) = - reac_rate_local(490) 
  reac_source_local(42,491) = + reac_rate_local(491) * 2.d0
  reac_source_local(47,491) = - reac_rate_local(491) 
  reac_source_local(50,491) = - reac_rate_local(491) 
  reac_source_local(14,492) = + reac_rate_local(492) 
  reac_source_local(17,492) = - reac_rate_local(492) 
  reac_source_local(43,492) = + reac_rate_local(492) 
  reac_source_local(51,492) = - reac_rate_local(492) 
  reac_source_local(01,493) = + reac_rate_local(493) 
  reac_source_local(18,493) = - reac_rate_local(493) 
  reac_source_local(43,493) = + reac_rate_local(493) 
  reac_source_local(51,493) = - reac_rate_local(493) 
  reac_source_local(29,494) = + reac_rate_local(494) 
  reac_source_local(33,494) = - reac_rate_local(494) 
  reac_source_local(43,494) = + reac_rate_local(494) 
  reac_source_local(51,494) = - reac_rate_local(494) 
  reac_source_local(21,495) = + reac_rate_local(495) 
  reac_source_local(34,495) = - reac_rate_local(495) 
  reac_source_local(43,495) = + reac_rate_local(495) 
  reac_source_local(51,495) = - reac_rate_local(495) 
  reac_source_local(40,496) = + reac_rate_local(496) 
  reac_source_local(43,496) = + reac_rate_local(496) 
  reac_source_local(45,496) = - reac_rate_local(496) 
  reac_source_local(51,496) = - reac_rate_local(496) 
  reac_source_local(41,497) = + reac_rate_local(497) 
  reac_source_local(43,497) = + reac_rate_local(497) 
  reac_source_local(46,497) = - reac_rate_local(497) 
  reac_source_local(51,497) = - reac_rate_local(497) 
  reac_source_local(42,498) = + reac_rate_local(498) 
  reac_source_local(43,498) = + reac_rate_local(498) 
  reac_source_local(47,498) = - reac_rate_local(498) 
  reac_source_local(51,498) = - reac_rate_local(498) 
  reac_source_local(14,499) = + reac_rate_local(499) * 2.d0
  reac_source_local(18,499) = - reac_rate_local(499) 
  reac_source_local(29,499) = + reac_rate_local(499) 
  reac_source_local(36,499) = - reac_rate_local(499) 
  reac_source_local(01,500) = + reac_rate_local(500) 
  reac_source_local(14,500) = + reac_rate_local(500) 
  reac_source_local(19,500) = - reac_rate_local(500) 
  reac_source_local(29,500) = + reac_rate_local(500) 
  reac_source_local(36,500) = - reac_rate_local(500) 
  reac_source_local(01,501) = + reac_rate_local(501) * 2.d0
  reac_source_local(20,501) = - reac_rate_local(501) 
  reac_source_local(29,501) = + reac_rate_local(501) 
  reac_source_local(36,501) = - reac_rate_local(501) 
  reac_source_local(29,502) = + reac_rate_local(502) * 3.d0
  reac_source_local(34,502) = - reac_rate_local(502) 
  reac_source_local(36,502) = - reac_rate_local(502) 
  reac_source_local(21,503) = + reac_rate_local(503) * 2.d0
  reac_source_local(29,503) = + reac_rate_local(503) 
  reac_source_local(35,503) = - reac_rate_local(503) 
  reac_source_local(36,503) = - reac_rate_local(503) 
  reac_source_local(14,504) = + reac_rate_local(504) 
  reac_source_local(29,504) = + reac_rate_local(504) * 2.d0
  reac_source_local(36,504) = - reac_rate_local(504) 
  reac_source_local(45,504) = - reac_rate_local(504) 
  reac_source_local(01,505) = + reac_rate_local(505) 
  reac_source_local(29,505) = + reac_rate_local(505) * 2.d0
  reac_source_local(36,505) = - reac_rate_local(505) 
  reac_source_local(46,505) = - reac_rate_local(505) 
  reac_source_local(14,506) = + reac_rate_local(506) 
  reac_source_local(21,506) = + reac_rate_local(506) 
  reac_source_local(29,506) = + reac_rate_local(506) 
  reac_source_local(36,506) = - reac_rate_local(506) 
  reac_source_local(47,506) = - reac_rate_local(506) 
  reac_source_local(01,507) = + reac_rate_local(507) 
  reac_source_local(21,507) = + reac_rate_local(507) 
  reac_source_local(29,507) = + reac_rate_local(507) 
  reac_source_local(36,507) = - reac_rate_local(507) 
  reac_source_local(52,507) = - reac_rate_local(507) 
  reac_source_local(14,508) = + reac_rate_local(508) * 2.d0
  reac_source_local(18,508) = - reac_rate_local(508) 
  reac_source_local(21,508) = + reac_rate_local(508) 
  reac_source_local(37,508) = - reac_rate_local(508) 
  reac_source_local(01,509) = + reac_rate_local(509) 
  reac_source_local(14,509) = + reac_rate_local(509) 
  reac_source_local(19,509) = - reac_rate_local(509) 
  reac_source_local(21,509) = + reac_rate_local(509) 
  reac_source_local(37,509) = - reac_rate_local(509) 
  reac_source_local(01,510) = + reac_rate_local(510) * 2.d0
  reac_source_local(20,510) = - reac_rate_local(510) 
  reac_source_local(21,510) = + reac_rate_local(510) 
  reac_source_local(37,510) = - reac_rate_local(510) 
  reac_source_local(21,511) = + reac_rate_local(511) 
  reac_source_local(29,511) = + reac_rate_local(511) * 2.d0
  reac_source_local(34,511) = - reac_rate_local(511) 
  reac_source_local(37,511) = - reac_rate_local(511) 
  reac_source_local(21,512) = + reac_rate_local(512) * 3.d0
  reac_source_local(35,512) = - reac_rate_local(512) 
  reac_source_local(37,512) = - reac_rate_local(512) 
  reac_source_local(14,513) = + reac_rate_local(513) 
  reac_source_local(21,513) = + reac_rate_local(513) 
  reac_source_local(29,513) = + reac_rate_local(513) 
  reac_source_local(37,513) = - reac_rate_local(513) 
  reac_source_local(45,513) = - reac_rate_local(513) 
  reac_source_local(01,514) = + reac_rate_local(514) 
  reac_source_local(21,514) = + reac_rate_local(514) 
  reac_source_local(29,514) = + reac_rate_local(514) 
  reac_source_local(37,514) = - reac_rate_local(514) 
  reac_source_local(46,514) = - reac_rate_local(514) 
  reac_source_local(14,515) = + reac_rate_local(515) 
  reac_source_local(21,515) = + reac_rate_local(515) * 2.d0
  reac_source_local(37,515) = - reac_rate_local(515) 
  reac_source_local(47,515) = - reac_rate_local(515) 
  reac_source_local(01,516) = + reac_rate_local(516) 
  reac_source_local(21,516) = + reac_rate_local(516) * 2.d0
  reac_source_local(37,516) = - reac_rate_local(516) 
  reac_source_local(52,516) = - reac_rate_local(516) 
  reac_source_local(14,517) = + reac_rate_local(517) * 2.d0
  reac_source_local(18,517) = - reac_rate_local(517) 
  reac_source_local(32,517) = + reac_rate_local(517) 
  reac_source_local(38,517) = - reac_rate_local(517) 
  reac_source_local(01,518) = + reac_rate_local(518) 
  reac_source_local(14,518) = + reac_rate_local(518) 
  reac_source_local(19,518) = - reac_rate_local(518) 
  reac_source_local(32,518) = + reac_rate_local(518) 
  reac_source_local(38,518) = - reac_rate_local(518) 
  reac_source_local(01,519) = + reac_rate_local(519) * 2.d0
  reac_source_local(20,519) = - reac_rate_local(519) 
  reac_source_local(32,519) = + reac_rate_local(519) 
  reac_source_local(38,519) = - reac_rate_local(519) 
  reac_source_local(29,520) = + reac_rate_local(520) * 2.d0
  reac_source_local(32,520) = + reac_rate_local(520) 
  reac_source_local(34,520) = - reac_rate_local(520) 
  reac_source_local(38,520) = - reac_rate_local(520) 
  reac_source_local(21,521) = + reac_rate_local(521) * 2.d0
  reac_source_local(32,521) = + reac_rate_local(521) 
  reac_source_local(35,521) = - reac_rate_local(521) 
  reac_source_local(38,521) = - reac_rate_local(521) 
  reac_source_local(14,522) = + reac_rate_local(522) 
  reac_source_local(29,522) = + reac_rate_local(522) 
  reac_source_local(32,522) = + reac_rate_local(522) 
  reac_source_local(38,522) = - reac_rate_local(522) 
  reac_source_local(45,522) = - reac_rate_local(522) 
  reac_source_local(01,523) = + reac_rate_local(523) 
  reac_source_local(29,523) = + reac_rate_local(523) 
  reac_source_local(32,523) = + reac_rate_local(523) 
  reac_source_local(38,523) = - reac_rate_local(523) 
  reac_source_local(46,523) = - reac_rate_local(523) 
  reac_source_local(14,524) = + reac_rate_local(524) 
  reac_source_local(21,524) = + reac_rate_local(524) 
  reac_source_local(32,524) = + reac_rate_local(524) 
  reac_source_local(38,524) = - reac_rate_local(524) 
  reac_source_local(47,524) = - reac_rate_local(524) 
  reac_source_local(01,525) = + reac_rate_local(525) 
  reac_source_local(21,525) = + reac_rate_local(525) 
  reac_source_local(32,525) = + reac_rate_local(525) 
  reac_source_local(38,525) = - reac_rate_local(525) 
  reac_source_local(52,525) = - reac_rate_local(525) 
  reac_source_local(14,526) = + reac_rate_local(526) * 2.d0
  reac_source_local(18,526) = - reac_rate_local(526) 
  reac_source_local(40,526) = + reac_rate_local(526) 
  reac_source_local(48,526) = - reac_rate_local(526) 
  reac_source_local(01,527) = + reac_rate_local(527) 
  reac_source_local(14,527) = + reac_rate_local(527) 
  reac_source_local(19,527) = - reac_rate_local(527) 
  reac_source_local(40,527) = + reac_rate_local(527) 
  reac_source_local(48,527) = - reac_rate_local(527) 
  reac_source_local(01,528) = + reac_rate_local(528) * 2.d0
  reac_source_local(20,528) = - reac_rate_local(528) 
  reac_source_local(40,528) = + reac_rate_local(528) 
  reac_source_local(48,528) = - reac_rate_local(528) 
  reac_source_local(29,529) = + reac_rate_local(529) * 2.d0
  reac_source_local(34,529) = - reac_rate_local(529) 
  reac_source_local(40,529) = + reac_rate_local(529) 
  reac_source_local(48,529) = - reac_rate_local(529) 
  reac_source_local(21,530) = + reac_rate_local(530) * 2.d0
  reac_source_local(35,530) = - reac_rate_local(530) 
  reac_source_local(40,530) = + reac_rate_local(530) 
  reac_source_local(48,530) = - reac_rate_local(530) 
  reac_source_local(14,531) = + reac_rate_local(531) 
  reac_source_local(29,531) = + reac_rate_local(531) 
  reac_source_local(40,531) = + reac_rate_local(531) 
  reac_source_local(45,531) = - reac_rate_local(531) 
  reac_source_local(48,531) = - reac_rate_local(531) 
  reac_source_local(01,532) = + reac_rate_local(532) 
  reac_source_local(29,532) = + reac_rate_local(532) 
  reac_source_local(40,532) = + reac_rate_local(532) 
  reac_source_local(46,532) = - reac_rate_local(532) 
  reac_source_local(48,532) = - reac_rate_local(532) 
  reac_source_local(14,533) = + reac_rate_local(533) 
  reac_source_local(21,533) = + reac_rate_local(533) 
  reac_source_local(40,533) = + reac_rate_local(533) 
  reac_source_local(47,533) = - reac_rate_local(533) 
  reac_source_local(48,533) = - reac_rate_local(533) 
  reac_source_local(01,534) = + reac_rate_local(534) 
  reac_source_local(21,534) = + reac_rate_local(534) 
  reac_source_local(40,534) = + reac_rate_local(534) 
  reac_source_local(48,534) = - reac_rate_local(534) 
  reac_source_local(52,534) = - reac_rate_local(534) 
  reac_source_local(14,535) = + reac_rate_local(535) * 2.d0
  reac_source_local(18,535) = - reac_rate_local(535) 
  reac_source_local(41,535) = + reac_rate_local(535) 
  reac_source_local(49,535) = - reac_rate_local(535) 
  reac_source_local(01,536) = + reac_rate_local(536) 
  reac_source_local(14,536) = + reac_rate_local(536) 
  reac_source_local(19,536) = - reac_rate_local(536) 
  reac_source_local(41,536) = + reac_rate_local(536) 
  reac_source_local(49,536) = - reac_rate_local(536) 
  reac_source_local(01,537) = + reac_rate_local(537) * 2.d0
  reac_source_local(20,537) = - reac_rate_local(537) 
  reac_source_local(41,537) = + reac_rate_local(537) 
  reac_source_local(49,537) = - reac_rate_local(537) 
  reac_source_local(29,538) = + reac_rate_local(538) * 2.d0
  reac_source_local(34,538) = - reac_rate_local(538) 
  reac_source_local(41,538) = + reac_rate_local(538) 
  reac_source_local(49,538) = - reac_rate_local(538) 
  reac_source_local(21,539) = + reac_rate_local(539) * 2.d0
  reac_source_local(35,539) = - reac_rate_local(539) 
  reac_source_local(41,539) = + reac_rate_local(539) 
  reac_source_local(49,539) = - reac_rate_local(539) 
  reac_source_local(14,540) = + reac_rate_local(540) 
  reac_source_local(29,540) = + reac_rate_local(540) 
  reac_source_local(41,540) = + reac_rate_local(540) 
  reac_source_local(45,540) = - reac_rate_local(540) 
  reac_source_local(49,540) = - reac_rate_local(540) 
  reac_source_local(01,541) = + reac_rate_local(541) 
  reac_source_local(29,541) = + reac_rate_local(541) 
  reac_source_local(41,541) = + reac_rate_local(541) 
  reac_source_local(46,541) = - reac_rate_local(541) 
  reac_source_local(49,541) = - reac_rate_local(541) 
  reac_source_local(14,542) = + reac_rate_local(542) 
  reac_source_local(21,542) = + reac_rate_local(542) 
  reac_source_local(41,542) = + reac_rate_local(542) 
  reac_source_local(47,542) = - reac_rate_local(542) 
  reac_source_local(49,542) = - reac_rate_local(542) 
  reac_source_local(01,543) = + reac_rate_local(543) 
  reac_source_local(21,543) = + reac_rate_local(543) 
  reac_source_local(41,543) = + reac_rate_local(543) 
  reac_source_local(49,543) = - reac_rate_local(543) 
  reac_source_local(52,543) = - reac_rate_local(543) 
  reac_source_local(14,544) = + reac_rate_local(544) * 2.d0
  reac_source_local(18,544) = - reac_rate_local(544) 
  reac_source_local(42,544) = + reac_rate_local(544) 
  reac_source_local(50,544) = - reac_rate_local(544) 
  reac_source_local(01,545) = + reac_rate_local(545) 
  reac_source_local(14,545) = + reac_rate_local(545) 
  reac_source_local(19,545) = - reac_rate_local(545) 
  reac_source_local(42,545) = + reac_rate_local(545) 
  reac_source_local(50,545) = - reac_rate_local(545) 
  reac_source_local(01,546) = + reac_rate_local(546) * 2.d0
  reac_source_local(20,546) = - reac_rate_local(546) 
  reac_source_local(42,546) = + reac_rate_local(546) 
  reac_source_local(50,546) = - reac_rate_local(546) 
  reac_source_local(29,547) = + reac_rate_local(547) * 2.d0
  reac_source_local(34,547) = - reac_rate_local(547) 
  reac_source_local(42,547) = + reac_rate_local(547) 
  reac_source_local(50,547) = - reac_rate_local(547) 
  reac_source_local(21,548) = + reac_rate_local(548) * 2.d0
  reac_source_local(35,548) = - reac_rate_local(548) 
  reac_source_local(42,548) = + reac_rate_local(548) 
  reac_source_local(50,548) = - reac_rate_local(548) 
  reac_source_local(14,549) = + reac_rate_local(549) 
  reac_source_local(29,549) = + reac_rate_local(549) 
  reac_source_local(42,549) = + reac_rate_local(549) 
  reac_source_local(45,549) = - reac_rate_local(549) 
  reac_source_local(50,549) = - reac_rate_local(549) 
  reac_source_local(01,550) = + reac_rate_local(550) 
  reac_source_local(29,550) = + reac_rate_local(550) 
  reac_source_local(42,550) = + reac_rate_local(550) 
  reac_source_local(46,550) = - reac_rate_local(550) 
  reac_source_local(50,550) = - reac_rate_local(550) 
  reac_source_local(14,551) = + reac_rate_local(551) 
  reac_source_local(21,551) = + reac_rate_local(551) 
  reac_source_local(42,551) = + reac_rate_local(551) 
  reac_source_local(47,551) = - reac_rate_local(551) 
  reac_source_local(50,551) = - reac_rate_local(551) 
  reac_source_local(01,552) = + reac_rate_local(552) 
  reac_source_local(21,552) = + reac_rate_local(552) 
  reac_source_local(42,552) = + reac_rate_local(552) 
  reac_source_local(50,552) = - reac_rate_local(552) 
  reac_source_local(52,552) = - reac_rate_local(552) 
  reac_source_local(14,553) = + reac_rate_local(553) * 2.d0
  reac_source_local(18,553) = - reac_rate_local(553) 
  reac_source_local(43,553) = + reac_rate_local(553) 
  reac_source_local(51,553) = - reac_rate_local(553) 
  reac_source_local(01,554) = + reac_rate_local(554) 
  reac_source_local(14,554) = + reac_rate_local(554) 
  reac_source_local(19,554) = - reac_rate_local(554) 
  reac_source_local(43,554) = + reac_rate_local(554) 
  reac_source_local(51,554) = - reac_rate_local(554) 
  reac_source_local(01,555) = + reac_rate_local(555) * 2.d0
  reac_source_local(20,555) = - reac_rate_local(555) 
  reac_source_local(43,555) = + reac_rate_local(555) 
  reac_source_local(51,555) = - reac_rate_local(555) 
  reac_source_local(29,556) = + reac_rate_local(556) * 2.d0
  reac_source_local(34,556) = - reac_rate_local(556) 
  reac_source_local(43,556) = + reac_rate_local(556) 
  reac_source_local(51,556) = - reac_rate_local(556) 
  reac_source_local(21,557) = + reac_rate_local(557) * 2.d0
  reac_source_local(35,557) = - reac_rate_local(557) 
  reac_source_local(43,557) = + reac_rate_local(557) 
  reac_source_local(51,557) = - reac_rate_local(557) 
  reac_source_local(14,558) = + reac_rate_local(558) 
  reac_source_local(29,558) = + reac_rate_local(558) 
  reac_source_local(43,558) = + reac_rate_local(558) 
  reac_source_local(45,558) = - reac_rate_local(558) 
  reac_source_local(51,558) = - reac_rate_local(558) 
  reac_source_local(01,559) = + reac_rate_local(559) 
  reac_source_local(29,559) = + reac_rate_local(559) 
  reac_source_local(43,559) = + reac_rate_local(559) 
  reac_source_local(46,559) = - reac_rate_local(559) 
  reac_source_local(51,559) = - reac_rate_local(559) 
  reac_source_local(14,560) = + reac_rate_local(560) 
  reac_source_local(21,560) = + reac_rate_local(560) 
  reac_source_local(43,560) = + reac_rate_local(560) 
  reac_source_local(47,560) = - reac_rate_local(560) 
  reac_source_local(51,560) = - reac_rate_local(560) 
  reac_source_local(01,561) = + reac_rate_local(561) 
  reac_source_local(21,561) = + reac_rate_local(561) 
  reac_source_local(43,561) = + reac_rate_local(561) 
  reac_source_local(51,561) = - reac_rate_local(561) 
  reac_source_local(52,561) = - reac_rate_local(561) 
  reac_source_local(14,562) = + reac_rate_local(562) 
  reac_source_local(17,562) = - reac_rate_local(562) 
  reac_source_local(21,562) = + reac_rate_local(562) * 2.d0
  reac_source_local(39,562) = - reac_rate_local(562) 
  reac_source_local(01,563) = + reac_rate_local(563) 
  reac_source_local(18,563) = - reac_rate_local(563) 
  reac_source_local(21,563) = + reac_rate_local(563) * 2.d0
  reac_source_local(39,563) = - reac_rate_local(563) 
  reac_source_local(21,564) = + reac_rate_local(564) * 2.d0
  reac_source_local(29,564) = + reac_rate_local(564) 
  reac_source_local(33,564) = - reac_rate_local(564) 
  reac_source_local(39,564) = - reac_rate_local(564) 
  reac_source_local(21,565) = + reac_rate_local(565) * 3.d0
  reac_source_local(34,565) = - reac_rate_local(565) 
  reac_source_local(39,565) = - reac_rate_local(565) 
  reac_source_local(21,566) = + reac_rate_local(566) * 2.d0
  reac_source_local(39,566) = - reac_rate_local(566) 
  reac_source_local(40,566) = + reac_rate_local(566) 
  reac_source_local(45,566) = - reac_rate_local(566) 
  reac_source_local(21,567) = + reac_rate_local(567) * 2.d0
  reac_source_local(39,567) = - reac_rate_local(567) 
  reac_source_local(41,567) = + reac_rate_local(567) 
  reac_source_local(46,567) = - reac_rate_local(567) 
  reac_source_local(21,568) = + reac_rate_local(568) * 2.d0
  reac_source_local(39,568) = - reac_rate_local(568) 
  reac_source_local(42,568) = + reac_rate_local(568) 
  reac_source_local(47,568) = - reac_rate_local(568) 
  reac_source_local(01,569) = + reac_rate_local(569) 
  reac_source_local(14,569) = + reac_rate_local(569) 
  reac_source_local(19,569) = - reac_rate_local(569) 
  reac_source_local(21,569) = + reac_rate_local(569) * 2.d0
  reac_source_local(39,569) = - reac_rate_local(569) 
  reac_source_local(01,570) = + reac_rate_local(570) * 2.d0
  reac_source_local(20,570) = - reac_rate_local(570) 
  reac_source_local(21,570) = + reac_rate_local(570) * 2.d0
  reac_source_local(39,570) = - reac_rate_local(570) 
  reac_source_local(21,571) = + reac_rate_local(571) * 4.d0
  reac_source_local(35,571) = - reac_rate_local(571) 
  reac_source_local(39,571) = - reac_rate_local(571) 
  reac_source_local(01,572) = + reac_rate_local(572) 
  reac_source_local(21,572) = + reac_rate_local(572) * 3.d0
  reac_source_local(39,572) = - reac_rate_local(572) 
  reac_source_local(52,572) = - reac_rate_local(572) 
  reac_source_local(14,573) = + reac_rate_local(573) 
  reac_source_local(17,573) = - reac_rate_local(573) 
  reac_source_local(29,573) = + reac_rate_local(573) 
  reac_source_local(36,573) = - reac_rate_local(573) 
  reac_source_local(01,574) = + reac_rate_local(574) 
  reac_source_local(18,574) = - reac_rate_local(574) 
  reac_source_local(29,574) = + reac_rate_local(574) 
  reac_source_local(36,574) = - reac_rate_local(574) 
  reac_source_local(29,575) = + reac_rate_local(575) * 2.d0
  reac_source_local(33,575) = - reac_rate_local(575) 
  reac_source_local(36,575) = - reac_rate_local(575) 
  reac_source_local(21,576) = + reac_rate_local(576) 
  reac_source_local(29,576) = + reac_rate_local(576) 
  reac_source_local(34,576) = - reac_rate_local(576) 
  reac_source_local(36,576) = - reac_rate_local(576) 
  reac_source_local(29,577) = + reac_rate_local(577) 
  reac_source_local(36,577) = - reac_rate_local(577) 
  reac_source_local(40,577) = + reac_rate_local(577) 
  reac_source_local(45,577) = - reac_rate_local(577) 
  reac_source_local(14,578) = + reac_rate_local(578) 
  reac_source_local(17,578) = - reac_rate_local(578) 
  reac_source_local(21,578) = + reac_rate_local(578) 
  reac_source_local(37,578) = - reac_rate_local(578) 
  reac_source_local(01,579) = + reac_rate_local(579) 
  reac_source_local(18,579) = - reac_rate_local(579) 
  reac_source_local(21,579) = + reac_rate_local(579) 
  reac_source_local(37,579) = - reac_rate_local(579) 
  reac_source_local(21,580) = + reac_rate_local(580) 
  reac_source_local(29,580) = + reac_rate_local(580) 
  reac_source_local(33,580) = - reac_rate_local(580) 
  reac_source_local(37,580) = - reac_rate_local(580) 
  reac_source_local(21,581) = + reac_rate_local(581) * 2.d0
  reac_source_local(34,581) = - reac_rate_local(581) 
  reac_source_local(37,581) = - reac_rate_local(581) 
  reac_source_local(21,582) = + reac_rate_local(582) 
  reac_source_local(37,582) = - reac_rate_local(582) 
  reac_source_local(40,582) = + reac_rate_local(582) 
  reac_source_local(45,582) = - reac_rate_local(582) 
  reac_source_local(17,583) = - reac_rate_local(583) 
  reac_source_local(36,583) = - reac_rate_local(583) 
  reac_source_local(40,583) = + reac_rate_local(583) 
  reac_source_local(18,584) = - reac_rate_local(584) 
  reac_source_local(36,584) = - reac_rate_local(584) 
  reac_source_local(41,584) = + reac_rate_local(584) 
  reac_source_local(21,585) = + reac_rate_local(585) 
  reac_source_local(33,585) = - reac_rate_local(585) 
  reac_source_local(36,585) = - reac_rate_local(585) 
  reac_source_local(32,586) = + reac_rate_local(586) 
  reac_source_local(34,586) = - reac_rate_local(586) 
  reac_source_local(36,586) = - reac_rate_local(586) 
  reac_source_local(36,587) = - reac_rate_local(587) 
  reac_source_local(42,587) = + reac_rate_local(587) 
  reac_source_local(45,587) = - reac_rate_local(587) 
  reac_source_local(17,588) = - reac_rate_local(588) 
  reac_source_local(37,588) = - reac_rate_local(588) 
  reac_source_local(42,588) = + reac_rate_local(588) 
  reac_source_local(32,589) = + reac_rate_local(589) 
  reac_source_local(33,589) = - reac_rate_local(589) 
  reac_source_local(37,589) = - reac_rate_local(589) 
  reac_source_local(37,590) = - reac_rate_local(590) 
  reac_source_local(43,590) = + reac_rate_local(590) 
  reac_source_local(45,590) = - reac_rate_local(590) 
  reac_source_local(14,591) = + reac_rate_local(591) 
  reac_source_local(17,591) = - reac_rate_local(591) 
  reac_source_local(32,591) = + reac_rate_local(591) 
  reac_source_local(38,591) = - reac_rate_local(591) 
  reac_source_local(01,592) = + reac_rate_local(592) 
  reac_source_local(18,592) = - reac_rate_local(592) 
  reac_source_local(32,592) = + reac_rate_local(592) 
  reac_source_local(38,592) = - reac_rate_local(592) 
  reac_source_local(29,593) = + reac_rate_local(593) 
  reac_source_local(32,593) = + reac_rate_local(593) 
  reac_source_local(33,593) = - reac_rate_local(593) 
  reac_source_local(38,593) = - reac_rate_local(593) 
  reac_source_local(21,594) = + reac_rate_local(594) 
  reac_source_local(32,594) = + reac_rate_local(594) 
  reac_source_local(34,594) = - reac_rate_local(594) 
  reac_source_local(38,594) = - reac_rate_local(594) 
  reac_source_local(32,595) = + reac_rate_local(595) 
  reac_source_local(38,595) = - reac_rate_local(595) 
  reac_source_local(40,595) = + reac_rate_local(595) 
  reac_source_local(45,595) = - reac_rate_local(595) 
  reac_source_local(32,596) = + reac_rate_local(596) 
  reac_source_local(38,596) = - reac_rate_local(596) 
  reac_source_local(41,596) = + reac_rate_local(596) 
  reac_source_local(46,596) = - reac_rate_local(596) 
  reac_source_local(32,597) = + reac_rate_local(597) 
  reac_source_local(38,597) = - reac_rate_local(597) 
  reac_source_local(42,597) = + reac_rate_local(597) 
  reac_source_local(47,597) = - reac_rate_local(597) 
  reac_source_local(14,598) = + reac_rate_local(598) 
  reac_source_local(17,598) = - reac_rate_local(598) 
  reac_source_local(40,598) = + reac_rate_local(598) 
  reac_source_local(48,598) = - reac_rate_local(598) 
  reac_source_local(01,599) = + reac_rate_local(599) 
  reac_source_local(18,599) = - reac_rate_local(599) 
  reac_source_local(40,599) = + reac_rate_local(599) 
  reac_source_local(48,599) = - reac_rate_local(599) 
  reac_source_local(29,600) = + reac_rate_local(600) 
  reac_source_local(33,600) = - reac_rate_local(600) 
  reac_source_local(40,600) = + reac_rate_local(600) 
  reac_source_local(48,600) = - reac_rate_local(600) 
  reac_source_local(21,601) = + reac_rate_local(601) 
  reac_source_local(34,601) = - reac_rate_local(601) 
  reac_source_local(40,601) = + reac_rate_local(601) 
  reac_source_local(48,601) = - reac_rate_local(601) 
  reac_source_local(40,602) = + reac_rate_local(602) * 2.d0
  reac_source_local(45,602) = - reac_rate_local(602) 
  reac_source_local(48,602) = - reac_rate_local(602) 
  reac_source_local(40,603) = + reac_rate_local(603) 
  reac_source_local(41,603) = + reac_rate_local(603) 
  reac_source_local(46,603) = - reac_rate_local(603) 
  reac_source_local(48,603) = - reac_rate_local(603) 
  reac_source_local(40,604) = + reac_rate_local(604) 
  reac_source_local(42,604) = + reac_rate_local(604) 
  reac_source_local(47,604) = - reac_rate_local(604) 
  reac_source_local(48,604) = - reac_rate_local(604) 
  reac_source_local(14,605) = + reac_rate_local(605) 
  reac_source_local(17,605) = - reac_rate_local(605) 
  reac_source_local(41,605) = + reac_rate_local(605) 
  reac_source_local(49,605) = - reac_rate_local(605) 
  reac_source_local(01,606) = + reac_rate_local(606) 
  reac_source_local(18,606) = - reac_rate_local(606) 
  reac_source_local(41,606) = + reac_rate_local(606) 
  reac_source_local(49,606) = - reac_rate_local(606) 
  reac_source_local(29,607) = + reac_rate_local(607) 
  reac_source_local(33,607) = - reac_rate_local(607) 
  reac_source_local(41,607) = + reac_rate_local(607) 
  reac_source_local(49,607) = - reac_rate_local(607) 
  reac_source_local(21,608) = + reac_rate_local(608) 
  reac_source_local(34,608) = - reac_rate_local(608) 
  reac_source_local(41,608) = + reac_rate_local(608) 
  reac_source_local(49,608) = - reac_rate_local(608) 
  reac_source_local(40,609) = + reac_rate_local(609) 
  reac_source_local(41,609) = + reac_rate_local(609) 
  reac_source_local(45,609) = - reac_rate_local(609) 
  reac_source_local(49,609) = - reac_rate_local(609) 
  reac_source_local(41,610) = + reac_rate_local(610) * 2.d0
  reac_source_local(46,610) = - reac_rate_local(610) 
  reac_source_local(49,610) = - reac_rate_local(610) 
  reac_source_local(41,611) = + reac_rate_local(611) 
  reac_source_local(42,611) = + reac_rate_local(611) 
  reac_source_local(47,611) = - reac_rate_local(611) 
  reac_source_local(49,611) = - reac_rate_local(611) 
  reac_source_local(14,612) = + reac_rate_local(612) 
  reac_source_local(17,612) = - reac_rate_local(612) 
  reac_source_local(42,612) = + reac_rate_local(612) 
  reac_source_local(50,612) = - reac_rate_local(612) 
  reac_source_local(01,613) = + reac_rate_local(613) 
  reac_source_local(18,613) = - reac_rate_local(613) 
  reac_source_local(42,613) = + reac_rate_local(613) 
  reac_source_local(50,613) = - reac_rate_local(613) 
  reac_source_local(29,614) = + reac_rate_local(614) 
  reac_source_local(33,614) = - reac_rate_local(614) 
  reac_source_local(42,614) = + reac_rate_local(614) 
  reac_source_local(50,614) = - reac_rate_local(614) 
  reac_source_local(21,615) = + reac_rate_local(615) 
  reac_source_local(34,615) = - reac_rate_local(615) 
  reac_source_local(42,615) = + reac_rate_local(615) 
  reac_source_local(50,615) = - reac_rate_local(615) 
  reac_source_local(40,616) = + reac_rate_local(616) 
  reac_source_local(42,616) = + reac_rate_local(616) 
  reac_source_local(45,616) = - reac_rate_local(616) 
  reac_source_local(50,616) = - reac_rate_local(616) 
  reac_source_local(41,617) = + reac_rate_local(617) 
  reac_source_local(42,617) = + reac_rate_local(617) 
  reac_source_local(46,617) = - reac_rate_local(617) 
  reac_source_local(50,617) = - reac_rate_local(617) 
  reac_source_local(42,618) = + reac_rate_local(618) * 2.d0
  reac_source_local(47,618) = - reac_rate_local(618) 
  reac_source_local(50,618) = - reac_rate_local(618) 
  reac_source_local(14,619) = + reac_rate_local(619) 
  reac_source_local(17,619) = - reac_rate_local(619) 
  reac_source_local(43,619) = + reac_rate_local(619) 
  reac_source_local(51,619) = - reac_rate_local(619) 
  reac_source_local(01,620) = + reac_rate_local(620) 
  reac_source_local(18,620) = - reac_rate_local(620) 
  reac_source_local(43,620) = + reac_rate_local(620) 
  reac_source_local(51,620) = - reac_rate_local(620) 
  reac_source_local(29,621) = + reac_rate_local(621) 
  reac_source_local(33,621) = - reac_rate_local(621) 
  reac_source_local(43,621) = + reac_rate_local(621) 
  reac_source_local(51,621) = - reac_rate_local(621) 
  reac_source_local(21,622) = + reac_rate_local(622) 
  reac_source_local(34,622) = - reac_rate_local(622) 
  reac_source_local(43,622) = + reac_rate_local(622) 
  reac_source_local(51,622) = - reac_rate_local(622) 
  reac_source_local(40,623) = + reac_rate_local(623) 
  reac_source_local(43,623) = + reac_rate_local(623) 
  reac_source_local(45,623) = - reac_rate_local(623) 
  reac_source_local(51,623) = - reac_rate_local(623) 
  reac_source_local(41,624) = + reac_rate_local(624) 
  reac_source_local(43,624) = + reac_rate_local(624) 
  reac_source_local(46,624) = - reac_rate_local(624) 
  reac_source_local(51,624) = - reac_rate_local(624) 
  reac_source_local(42,625) = + reac_rate_local(625) 
  reac_source_local(43,625) = + reac_rate_local(625) 
  reac_source_local(47,625) = - reac_rate_local(625) 
  reac_source_local(51,625) = - reac_rate_local(625) 
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(54)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(01) * density(02) 
  rrt(002) = rrt(002) * density(01) * density(03) 
  rrt(003) = rrt(003) * density(01) * density(04) 
  rrt(004) = rrt(004) * density(01) * density(05) 
  rrt(005) = rrt(005) * density(01) * density(06) 
  rrt(006) = rrt(006) * density(01) * density(07) 
  rrt(007) = rrt(007) * density(01) * density(08) 
  rrt(008) = rrt(008) * density(01) * density(09) 
  rrt(009) = rrt(009) * density(01)**2 
  rrt(010) = rrt(010) * density(01) * density(02) 
  rrt(011) = rrt(011) * density(01) * density(03) 
  rrt(012) = rrt(012) * density(01) * density(04) 
  rrt(013) = rrt(013) * density(01) * density(05) 
  rrt(014) = rrt(014) * density(01) * density(06) 
  rrt(015) = rrt(015) * density(01) * density(07) 
  rrt(016) = rrt(016) * density(01) * density(08) 
  rrt(017) = rrt(017) * density(02) * density(14) 
  rrt(018) = rrt(018) * density(03) * density(14) 
  rrt(019) = rrt(019) * density(04) * density(14) 
  rrt(020) = rrt(020) * density(05) * density(14) 
  rrt(021) = rrt(021) * density(06) * density(14) 
  rrt(022) = rrt(022) * density(07) * density(14) 
  rrt(023) = rrt(023) * density(08) * density(14) 
  rrt(024) = rrt(024) * density(09) * density(14) 
  rrt(025) = rrt(025) * density(01) * density(14) 
  rrt(026) = rrt(026) * density(02) * density(14) 
  rrt(027) = rrt(027) * density(03) * density(14) 
  rrt(028) = rrt(028) * density(04) * density(14) 
  rrt(029) = rrt(029) * density(05) * density(14) 
  rrt(030) = rrt(030) * density(06) * density(14) 
  rrt(031) = rrt(031) * density(07) * density(14) 
  rrt(032) = rrt(032) * density(08) * density(14) 
  rrt(033) = rrt(033) * density(02) * density(29) 
  rrt(034) = rrt(034) * density(03) * density(29) 
  rrt(035) = rrt(035) * density(04) * density(29) 
  rrt(036) = rrt(036) * density(05) * density(29) 
  rrt(037) = rrt(037) * density(06) * density(29) 
  rrt(038) = rrt(038) * density(07) * density(29) 
  rrt(039) = rrt(039) * density(08) * density(29) 
  rrt(040) = rrt(040) * density(09) * density(29) 
  rrt(041) = rrt(041) * density(01) * density(29) 
  rrt(042) = rrt(042) * density(02) * density(29) 
  rrt(043) = rrt(043) * density(03) * density(29) 
  rrt(044) = rrt(044) * density(04) * density(29) 
  rrt(045) = rrt(045) * density(05) * density(29) 
  rrt(046) = rrt(046) * density(06) * density(29) 
  rrt(047) = rrt(047) * density(07) * density(29) 
  rrt(048) = rrt(048) * density(08) * density(29) 
  rrt(049) = rrt(049) * density(21) * density(22) 
  rrt(050) = rrt(050) * density(21) * density(23) 
  rrt(051) = rrt(051) * density(21) * density(24) 
  rrt(052) = rrt(052) * density(21) * density(25) 
  rrt(053) = rrt(053) * density(21)**2 
  rrt(054) = rrt(054) * density(21) * density(22) 
  rrt(055) = rrt(055) * density(21) * density(23) 
  rrt(056) = rrt(056) * density(21) * density(24) 
  rrt(057) = rrt(057) * density(22) * density(29) 
  rrt(058) = rrt(058) * density(23) * density(29) 
  rrt(059) = rrt(059) * density(24) * density(29) 
  rrt(060) = rrt(060) * density(25) * density(29) 
  rrt(061) = rrt(061) * density(21) * density(29) 
  rrt(062) = rrt(062) * density(22) * density(29) 
  rrt(063) = rrt(063) * density(23) * density(29) 
  rrt(064) = rrt(064) * density(24) * density(29) 
  rrt(065) = rrt(065) * density(01) * density(53) 
  rrt(066) = rrt(066) * density(01) * density(53) 
  rrt(067) = rrt(067) * density(01) * density(53) 
  rrt(068) = rrt(068) * density(01) * density(53) 
  rrt(069) = rrt(069) * density(01) * density(53) 
  rrt(070) = rrt(070) * density(01) * density(53) 
  rrt(071) = rrt(071) * density(01) * density(53) 
  rrt(072) = rrt(072) * density(01) * density(53) 
  rrt(073) = rrt(073) * density(21) * density(53) 
  rrt(074) = rrt(074) * density(21) * density(53) 
  rrt(075) = rrt(075) * density(21) * density(53) 
  rrt(076) = rrt(076) * density(21) * density(53) 
  rrt(077) = rrt(077) * density(21) * density(53) 
  rrt(078) = rrt(078) * density(21) * density(53) 
  rrt(079) = rrt(079) * density(26) * density(53) 
  rrt(080) = rrt(080) * density(29) * density(53) 
  rrt(081) = rrt(081) * density(29) * density(53) 
  rrt(082) = rrt(082) * density(10) * density(53) 
  rrt(083) = rrt(083) * density(26) * density(53) 
  rrt(084) = rrt(084) * density(14) * density(53) 
  rrt(085) = rrt(085) * density(29) * density(53) 
  rrt(086) = rrt(086) * density(01) * density(53) 
  rrt(087) = rrt(087) * density(10) * density(53) 
  rrt(088) = rrt(088) * density(21) * density(53) 
  rrt(089) = rrt(089) * density(26) * density(53) 
  rrt(090) = rrt(090) * density(40) * density(53) 
  rrt(091) = rrt(091) * density(41) * density(53) 
  rrt(092) = rrt(092) * density(18) * density(53) 
  rrt(093) = rrt(093) * density(18) * density(53) 
  rrt(094) = rrt(094) * density(18) * density(53) 
  rrt(095) = rrt(095) * density(34) * density(53) 
  rrt(096) = rrt(096) * density(34) * density(53) 
  rrt(097) = rrt(097) * density(34) * density(53) 
  rrt(098) = rrt(098) * density(45) * density(53) 
  rrt(099) = rrt(099) * density(45) * density(53) 
  rrt(100) = rrt(100) * density(19) * density(53) 
  rrt(101) = rrt(101) * density(20) * density(53) 
  rrt(102) = rrt(102) * density(46) * density(53) 
  rrt(103) = rrt(103) * density(47) * density(53) 
  rrt(104) = rrt(104) * density(35) * density(53) 
  rrt(105) = rrt(105) * density(52) * density(53) 
  rrt(106) = rrt(106) * density(17) * density(53)**2 
  rrt(107) = rrt(107) * density(33) * density(53)**2 
  rrt(108) = rrt(108) * density(17) * density(53) 
  rrt(109) = rrt(109) * density(33) * density(53) 
  rrt(110) = rrt(110) * density(21) * density(53) 
  rrt(111) = rrt(111) * density(40) * density(53) 
  rrt(112) = rrt(112) * density(32) * density(53) 
  rrt(113) = rrt(113) * density(32) * density(53) 
  rrt(114) = rrt(114) * density(41) * density(53) 
  rrt(115) = rrt(115) * density(21)**2 * density(53) 
  rrt(116) = rrt(116) * density(42) * density(53) 
  rrt(117) = rrt(117) * density(21) * density(29) * density(53) 
  rrt(118) = rrt(118) * density(21) * density(29) * density(53) 
  rrt(119) = rrt(119) * density(32) * density(53) 
  rrt(120) = rrt(120) * density(40) * density(53) 
  rrt(121) = rrt(121) * density(41) * density(53) 
  rrt(122) = rrt(122) * density(01) * density(21) * density(53) 
  rrt(123) = rrt(123) * density(29) * density(36) 
  rrt(124) = rrt(124) * density(14) * density(36) 
  rrt(125) = rrt(125) * density(36) * density(40) 
  rrt(126) = rrt(126) * density(01) * density(36) 
  rrt(127) = rrt(127) * density(21) * density(36) 
  rrt(128) = rrt(128) * density(26) * density(36) 
  rrt(129) = rrt(129) * density(27) * density(36) 
  rrt(130) = rrt(130) * density(10) * density(36) 
  rrt(131) = rrt(131) * density(11) * density(36) 
  rrt(132) = rrt(132) * density(32) * density(36) 
  rrt(133) = rrt(133) * density(29) * density(37) 
  rrt(134) = rrt(134) * density(14) * density(37) 
  rrt(135) = rrt(135) * density(21) * density(37) 
  rrt(136) = rrt(136) * density(26) * density(37) 
  rrt(137) = rrt(137) * density(27) * density(37) 
  rrt(138) = rrt(138) * density(01) * density(37) 
  rrt(139) = rrt(139) * density(10) * density(37) 
  rrt(140) = rrt(140) * density(11) * density(37) 
  rrt(141) = rrt(141) * density(29) * density(38) 
  rrt(142) = rrt(142) * density(14) * density(48) 
  rrt(143) = rrt(143) * density(14) * density(38) 
  rrt(144) = rrt(144) * density(14) * density(49) 
  rrt(145) = rrt(145) * density(14) * density(50) 
  rrt(146) = rrt(146) * density(14) * density(51) 
  rrt(147) = rrt(147) * density(29) * density(48) 
  rrt(148) = rrt(148) * density(29) * density(49) 
  rrt(149) = rrt(149) * density(29) * density(50) 
  rrt(150) = rrt(150) * density(29) * density(51) 
  rrt(151) = rrt(151) * density(10) * density(38) 
  rrt(152) = rrt(152) * density(10) * density(48) 
  rrt(153) = rrt(153) * density(10) * density(49) 
  rrt(154) = rrt(154) * density(10) * density(50) 
  rrt(155) = rrt(155) * density(10) * density(51) 
  rrt(156) = rrt(156) * density(11) * density(38) 
  rrt(157) = rrt(157) * density(11) * density(48) 
  rrt(158) = rrt(158) * density(11) * density(49) 
  rrt(159) = rrt(159) * density(11) * density(50) 
  rrt(160) = rrt(160) * density(11) * density(51) 
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
  rrt(176) = rrt(176) * density(10) * density(21) 
  rrt(177) = rrt(177) * density(01) * density(10) 
  rrt(178) = rrt(178) * density(10) * density(40) 
  rrt(179) = rrt(179) * density(10) * density(41) 
  rrt(180) = rrt(180) * density(10) * density(42) 
  rrt(181) = rrt(181) * density(10)**2 
  rrt(182) = rrt(182) * density(10)**2 
  rrt(183) = rrt(183) * density(01) * density(11) 
  rrt(184) = rrt(184) * density(01) * density(11) 
  rrt(185) = rrt(185) * density(11) * density(21) 
  rrt(186) = rrt(186) * density(11) * density(40) 
  rrt(187) = rrt(187) * density(01) * density(13) 
  rrt(188) = rrt(188) * density(13) * density(21) 
  rrt(189) = rrt(189) * density(01) * density(12) 
  rrt(190) = rrt(190) * density(12) * density(21) 
  rrt(191) = rrt(191) * density(12) * density(40) 
  rrt(192) = rrt(192) * density(10) * density(12) 
  rrt(193) = rrt(193) * density(12)**2 
  rrt(194) = rrt(194) * density(01) * density(14)**2 
  rrt(195) = rrt(195) * density(14)**2 * density(21) 
  rrt(196) = rrt(196) * density(14)**2 * density(40) 
  rrt(197) = rrt(197) * density(14)**3 
  rrt(198) = rrt(198) * density(14)**2 * density(29) 
  rrt(199) = rrt(199) * density(01) * density(14)**2 
  rrt(200) = rrt(200) * density(14)**2 * density(21) 
  rrt(201) = rrt(201) * density(14)**2 * density(40) 
  rrt(202) = rrt(202) * density(14)**3 
  rrt(203) = rrt(203) * density(14)**2 * density(29) 
  rrt(204) = rrt(204) * density(15) * density(29) 
  rrt(205) = rrt(205) * density(15) * density(21) 
  rrt(206) = rrt(206) * density(15) * density(40) 
  rrt(207) = rrt(207) * density(15) * density(41) 
  rrt(208) = rrt(208) * density(01) * density(15) 
  rrt(209) = rrt(209) * density(14) * density(16) 
  rrt(210) = rrt(210) * density(16) * density(29) 
  rrt(211) = rrt(211) * density(14) * density(16) 
  rrt(212) = rrt(212) * density(01) * density(16) 
  rrt(213) = rrt(213) * density(15) * density(16) 
  rrt(214) = rrt(214) * density(16) * density(21) 
  rrt(215) = rrt(215) * density(16) * density(40) 
  rrt(216) = rrt(216) * density(26) * density(29) 
  rrt(217) = rrt(217) * density(14) * density(26) 
  rrt(218) = rrt(218) * density(21) * density(26) 
  rrt(219) = rrt(219) * density(01) * density(26) 
  rrt(220) = rrt(220) * density(26) * density(40) 
  rrt(221) = rrt(221) * density(26) * density(32) 
  rrt(222) = rrt(222) * density(26)**2 
  rrt(223) = rrt(223) * density(29) * density(32) 
  rrt(224) = rrt(224) * density(27) * density(29) 
  rrt(225) = rrt(225) * density(27) * density(29) 
  rrt(226) = rrt(226) * density(21) * density(27) 
  rrt(227) = rrt(227) * density(01) * density(27) 
  rrt(228) = rrt(228) * density(27) * density(40) 
  rrt(229) = rrt(229) * density(27) * density(32) 
  rrt(230) = rrt(230) * density(28) * density(29) 
  rrt(231) = rrt(231) * density(21) * density(28) 
  rrt(232) = rrt(232) * density(01) * density(28) 
  rrt(233) = rrt(233) * density(29) * density(30) 
  rrt(234) = rrt(234) * density(21) * density(30) 
  rrt(235) = rrt(235) * density(21) * density(30) 
  rrt(236) = rrt(236) * density(21) * density(30) 
  rrt(237) = rrt(237) * density(01) * density(30) 
  rrt(238) = rrt(238) * density(30) * density(32) 
  rrt(239) = rrt(239) * density(30) * density(32) 
  rrt(240) = rrt(240) * density(30) * density(40) 
  rrt(241) = rrt(241) * density(30) * density(41) 
  rrt(242) = rrt(242) * density(30) * density(41) 
  rrt(243) = rrt(243) * density(29) * density(31) 
  rrt(244) = rrt(244) * density(14) * density(31) 
  rrt(245) = rrt(245) * density(21) * density(31) 
  rrt(246) = rrt(246) * density(21) * density(31) 
  rrt(247) = rrt(247) * density(01) * density(31) 
  rrt(248) = rrt(248) * density(26) * density(31) 
  rrt(249) = rrt(249) * density(26) * density(31) 
  rrt(250) = rrt(250) * density(26) * density(31) 
  rrt(251) = rrt(251) * density(31) * density(40) 
  rrt(252) = rrt(252) * density(31) * density(40) 
  rrt(253) = rrt(253) * density(31) * density(32) 
  rrt(254) = rrt(254) * density(31) * density(32) 
  rrt(255) = rrt(255) * density(31) * density(41) 
  rrt(256) = rrt(256) * density(31) * density(41) 
  rrt(257) = rrt(257) * density(14) * density(40) 
  rrt(258) = rrt(258) * density(14) * density(21) 
  rrt(259) = rrt(259) * density(14) * density(42) 
  rrt(260) = rrt(260) * density(14) * density(42) 
  rrt(261) = rrt(261) * density(14) * density(42) 
  rrt(262) = rrt(262) * density(14) * density(42) 
  rrt(263) = rrt(263) * density(01) * density(29) 
  rrt(264) = rrt(264) * density(29) * density(40) 
  rrt(265) = rrt(265) * density(29) * density(40) 
  rrt(266) = rrt(266) * density(29) * density(41) 
  rrt(267) = rrt(267) * density(29) * density(41) 
  rrt(268) = rrt(268) * density(29) * density(42) 
  rrt(269) = rrt(269) * density(29) * density(43) 
  rrt(270) = rrt(270) * density(01) * density(21) 
  rrt(271) = rrt(271) * density(40)**2 
  rrt(272) = rrt(272) * density(40)**2 
  rrt(273) = rrt(273) * density(40)**2 
  rrt(274) = rrt(274) * density(21) * density(40) 
  rrt(275) = rrt(275) * density(32) * density(40) 
  rrt(276) = rrt(276) * density(40) * density(41) 
  rrt(277) = rrt(277) * density(40) * density(43) 
  rrt(278) = rrt(278) * density(21)**2 
  rrt(279) = rrt(279) * density(21) * density(42) 
  rrt(280) = rrt(280) * density(42)**2 
  rrt(281) = rrt(281) * density(42)**2 
  rrt(282) = rrt(282) * density(32) * density(42) 
  rrt(283) = rrt(283) * density(42) * density(43) 
  rrt(284) = rrt(284) * density(21) * density(43) 
  rrt(285) = rrt(285) * density(43)**2 
  rrt(286) = rrt(286) * density(14)**2 
  rrt(287) = rrt(287) * density(14) * density(29) 
  rrt(288) = rrt(288) * density(01)**2 
  rrt(289) = rrt(289) * density(01) * density(21) 
  rrt(290) = rrt(290) * density(01) * density(40) 
  rrt(291) = rrt(291) * density(01) * density(29) 
  rrt(292) = rrt(292) * density(01) * density(14) 
  rrt(293) = rrt(293) * density(01) * density(21) 
  rrt(294) = rrt(294) * density(21)**2 
  rrt(295) = rrt(295) * density(21) * density(29) 
  rrt(296) = rrt(296) * density(14) * density(21) 
  rrt(297) = rrt(297) * density(21) * density(40) 
  rrt(298) = rrt(298) * density(01) * density(40) 
  rrt(299) = rrt(299) * density(21) * density(40) 
  rrt(300) = rrt(300) * density(29) * density(40) 
  rrt(301) = rrt(301) * density(14) * density(40) 
  rrt(302) = rrt(302) * density(40)**2 
  rrt(303) = rrt(303) * density(01) * density(32) 
  rrt(304) = rrt(304) * density(21) * density(32) 
  rrt(305) = rrt(305) * density(14) * density(32) 
  rrt(306) = rrt(306) * density(29) * density(32) 
  rrt(307) = rrt(307) * density(01) * density(41) 
  rrt(308) = rrt(308) * density(21) * density(41) 
  rrt(309) = rrt(309) * density(40) * density(41) 
  rrt(310) = rrt(310) * density(41)**2 
  rrt(311) = rrt(311) * density(01) * density(42) 
  rrt(312) = rrt(312) * density(21) * density(42) 
  rrt(313) = rrt(313) * density(40) * density(42) 
  rrt(314) = rrt(314) * density(42)**2 
  rrt(315) = rrt(315) * density(01) * density(43) 
  rrt(316) = rrt(316) * density(21) * density(43) 
  rrt(317) = rrt(317) * density(40) * density(43) 
  rrt(318) = rrt(318) * density(14) * density(43) 
  rrt(319) = rrt(319) * density(29) * density(43) 
  rrt(320) = rrt(320) * density(01) * density(43) 
  rrt(321) = rrt(321) * density(21) * density(43) 
  rrt(322) = rrt(322) * density(40) * density(43) 
  rrt(323) = rrt(323) * density(14) * density(43) 
  rrt(324) = rrt(324) * density(29) * density(43) 
  rrt(325) = rrt(325) * density(44) 
  rrt(326) = rrt(326) * density(01) * density(14)**2 
  rrt(327) = rrt(327) * density(14)**2 * density(21) 
  rrt(328) = rrt(328) * density(14)**2 * density(40) 
  rrt(329) = rrt(329) * density(14)**3 
  rrt(330) = rrt(330) * density(14)**2 * density(29) 
  rrt(331) = rrt(331) * density(01) * density(29)**2 
  rrt(332) = rrt(332) * density(21) * density(29)**2 
  rrt(333) = rrt(333) * density(14) * density(29)**2 
  rrt(334) = rrt(334) * density(29)**3 
  rrt(335) = rrt(335) * density(29)**2 * density(40) 
  rrt(336) = rrt(336) * density(01) * density(14) * density(29) 
  rrt(337) = rrt(337) * density(14) * density(21) * density(29) 
  rrt(338) = rrt(338) * density(14)**2 * density(29) 
  rrt(339) = rrt(339) * density(14) * density(29)**2 
  rrt(340) = rrt(340) * density(14) * density(29) * density(40) 
  rrt(341) = rrt(341) * density(01) * density(21) * density(29) 
  rrt(342) = rrt(342) * density(21)**2 * density(29) 
  rrt(343) = rrt(343) * density(21) * density(29) * density(40) 
  rrt(344) = rrt(344) * density(14) * density(21) * density(29) 
  rrt(345) = rrt(345) * density(21) * density(29)**2 
  rrt(346) = rrt(346) * density(01) * density(29) 
  rrt(347) = rrt(347) * density(01) * density(29) * density(40) 
  rrt(348) = rrt(348) * density(21) * density(29) * density(40) 
  rrt(349) = rrt(349) * density(29) * density(40)**2 
  rrt(350) = rrt(350) * density(01) * density(29) * density(42) 
  rrt(351) = rrt(351) * density(21) * density(29) * density(42) 
  rrt(352) = rrt(352) * density(14) * density(29) * density(42) 
  rrt(353) = rrt(353) * density(29)**2 * density(42) 
  rrt(354) = rrt(354) * density(29) * density(40) * density(42) 
  rrt(355) = rrt(355) * density(42) * density(43) 
  rrt(356) = rrt(356) * density(17) * density(29) 
  rrt(357) = rrt(357) * density(17) * density(21) 
  rrt(358) = rrt(358) * density(17) * density(21) 
  rrt(359) = rrt(359) * density(17) * density(21) 
  rrt(360) = rrt(360) * density(17) * density(32) 
  rrt(361) = rrt(361) * density(17) * density(40) 
  rrt(362) = rrt(362) * density(17) * density(40) 
  rrt(363) = rrt(363) * density(17) * density(40) 
  rrt(364) = rrt(364) * density(17) * density(41) 
  rrt(365) = rrt(365) * density(01) * density(33) 
  rrt(366) = rrt(366) * density(21) * density(33) 
  rrt(367) = rrt(367) * density(32) * density(33) 
  rrt(368) = rrt(368) * density(33) * density(40) 
  rrt(369) = rrt(369) * density(33) * density(40) 
  rrt(370) = rrt(370) * density(15) * density(33) 
  rrt(371) = rrt(371) * density(33) * density(41) 
  rrt(372) = rrt(372) * density(33) * density(41) 
  rrt(373) = rrt(373) * density(33) * density(41) 
  rrt(374) = rrt(374) * density(33) * density(42) 
  rrt(375) = rrt(375) * density(18) * density(21) 
  rrt(376) = rrt(376) * density(18) * density(29) 
  rrt(377) = rrt(377) * density(18) * density(32) 
  rrt(378) = rrt(378) * density(14) * density(18) 
  rrt(379) = rrt(379) * density(18) * density(40) 
  rrt(380) = rrt(380) * density(18) * density(41) 
  rrt(381) = rrt(381) * density(18) * density(41) 
  rrt(382) = rrt(382) * density(01) * density(34) 
  rrt(383) = rrt(383) * density(14) * density(34) 
  rrt(384) = rrt(384) * density(34) * density(40) 
  rrt(385) = rrt(385) * density(34) * density(42) 
  rrt(386) = rrt(386) * density(34) * density(42) 
  rrt(387) = rrt(387) * density(19) * density(21) 
  rrt(388) = rrt(388) * density(19) * density(21) 
  rrt(389) = rrt(389) * density(14) * density(19) 
  rrt(390) = rrt(390) * density(19) * density(40) 
  rrt(391) = rrt(391) * density(19) * density(40) 
  rrt(392) = rrt(392) * density(40) * density(47) 
  rrt(393) = rrt(393) * density(40) * density(46) 
  rrt(394) = rrt(394) * density(01) * density(20) 
  rrt(395) = rrt(395) * density(20) * density(21) 
  rrt(396) = rrt(396) * density(20) * density(29) 
  rrt(397) = rrt(397) * density(14) * density(20) 
  rrt(398) = rrt(398) * density(20) * density(40) 
  rrt(399) = rrt(399) * density(01) * density(35) 
  rrt(400) = rrt(400) * density(21) * density(35) 
  rrt(401) = rrt(401) * density(26) * density(35) 
  rrt(402) = rrt(402) * density(27) * density(35) 
  rrt(403) = rrt(403) * density(29) * density(35) 
  rrt(404) = rrt(404) * density(35) * density(40) 
  rrt(405) = rrt(405) * density(01) * density(52) 
  rrt(406) = rrt(406) * density(21) * density(52) 
  rrt(407) = rrt(407) * density(01)**2 * density(17) 
  rrt(408) = rrt(408) * density(17) * density(29) 
  rrt(409) = rrt(409) * density(14) * density(17) 
  rrt(410) = rrt(410) * density(01) * density(33) 
  rrt(411) = rrt(411) * density(29) * density(33) 
  rrt(412) = rrt(412) * density(14) * density(33) 
  rrt(413) = rrt(413) * density(01)**2 * density(18) 
  rrt(414) = rrt(414) * density(01) * density(14) * density(18) 
  rrt(415) = rrt(415) * density(21)**2 * density(34) 
  rrt(416) = rrt(416) * density(01)**2 * density(34) 
  rrt(417) = rrt(417) * density(26) * density(36) 
  rrt(418) = rrt(418) * density(32) * density(36) 
  rrt(419) = rrt(419) * density(36) * density(42) 
  rrt(420) = rrt(420) * density(36) * density(41) 
  rrt(421) = rrt(421) * density(36) * density(41) 
  rrt(422) = rrt(422) * density(29) * density(37) 
  rrt(423) = rrt(423) * density(32) * density(37) 
  rrt(424) = rrt(424) * density(37) * density(42) 
  rrt(425) = rrt(425) * density(37) * density(43) 
  rrt(426) = rrt(426) * density(29) * density(38) 
  rrt(427) = rrt(427) * density(38) * density(40) 
  rrt(428) = rrt(428) * density(38) * density(40) 
  rrt(429) = rrt(429) * density(38) * density(42) 
  rrt(430) = rrt(430) * density(38) * density(42) 
  rrt(431) = rrt(431) * density(38) * density(43) 
  rrt(432) = rrt(432) * density(21) * density(48) 
  rrt(433) = rrt(433) * density(42) * density(48) 
  rrt(434) = rrt(434) * density(41) * density(48) 
  rrt(435) = rrt(435) * density(32) * density(50) 
  rrt(436) = rrt(436) * density(42) * density(50) 
  rrt(437) = rrt(437) * density(43) * density(50) 
  rrt(438) = rrt(438) * density(44) * density(50) 
  rrt(439) = rrt(439) * density(40) * density(51) 
  rrt(440) = rrt(440) * density(01) * density(39) 
  rrt(441) = rrt(441) * density(21) * density(39) 
  rrt(442) = rrt(442) * density(29) * density(39) 
  rrt(443) = rrt(443) * density(29) * density(39) 
  rrt(444) = rrt(444) * density(26) * density(39) 
  rrt(445) = rrt(445) * density(27) * density(39) 
  rrt(446) = rrt(446) * density(39) * density(40) 
  rrt(447) = rrt(447) * density(21) * density(36) 
  rrt(448) = rrt(448) * density(36) * density(40) 
  rrt(449) = rrt(449) * density(21) * density(37) 
  rrt(450) = rrt(450) * density(17) * density(36) 
  rrt(451) = rrt(451) * density(18) * density(36) 
  rrt(452) = rrt(452) * density(33) * density(36) 
  rrt(453) = rrt(453) * density(34) * density(36) 
  rrt(454) = rrt(454) * density(36) * density(45) 
  rrt(455) = rrt(455) * density(36) * density(46) 
  rrt(456) = rrt(456) * density(36) * density(47) 
  rrt(457) = rrt(457) * density(17) * density(37) 
  rrt(458) = rrt(458) * density(18) * density(37) 
  rrt(459) = rrt(459) * density(33) * density(37) 
  rrt(460) = rrt(460) * density(34) * density(37) 
  rrt(461) = rrt(461) * density(37) * density(45) 
  rrt(462) = rrt(462) * density(37) * density(46) 
  rrt(463) = rrt(463) * density(37) * density(47) 
  rrt(464) = rrt(464) * density(17) * density(38) 
  rrt(465) = rrt(465) * density(18) * density(38) 
  rrt(466) = rrt(466) * density(33) * density(38) 
  rrt(467) = rrt(467) * density(34) * density(38) 
  rrt(468) = rrt(468) * density(38) * density(45) 
  rrt(469) = rrt(469) * density(38) * density(46) 
  rrt(470) = rrt(470) * density(38) * density(47) 
  rrt(471) = rrt(471) * density(17) * density(48) 
  rrt(472) = rrt(472) * density(18) * density(48) 
  rrt(473) = rrt(473) * density(33) * density(48) 
  rrt(474) = rrt(474) * density(34) * density(48) 
  rrt(475) = rrt(475) * density(45) * density(48) 
  rrt(476) = rrt(476) * density(46) * density(48) 
  rrt(477) = rrt(477) * density(47) * density(48) 
  rrt(478) = rrt(478) * density(17) * density(49) 
  rrt(479) = rrt(479) * density(18) * density(49) 
  rrt(480) = rrt(480) * density(33) * density(49) 
  rrt(481) = rrt(481) * density(34) * density(49) 
  rrt(482) = rrt(482) * density(45) * density(49) 
  rrt(483) = rrt(483) * density(46) * density(49) 
  rrt(484) = rrt(484) * density(47) * density(49) 
  rrt(485) = rrt(485) * density(17) * density(50) 
  rrt(486) = rrt(486) * density(18) * density(50) 
  rrt(487) = rrt(487) * density(33) * density(50) 
  rrt(488) = rrt(488) * density(34) * density(50) 
  rrt(489) = rrt(489) * density(45) * density(50) 
  rrt(490) = rrt(490) * density(46) * density(50) 
  rrt(491) = rrt(491) * density(47) * density(50) 
  rrt(492) = rrt(492) * density(17) * density(51) 
  rrt(493) = rrt(493) * density(18) * density(51) 
  rrt(494) = rrt(494) * density(33) * density(51) 
  rrt(495) = rrt(495) * density(34) * density(51) 
  rrt(496) = rrt(496) * density(45) * density(51) 
  rrt(497) = rrt(497) * density(46) * density(51) 
  rrt(498) = rrt(498) * density(47) * density(51) 
  rrt(499) = rrt(499) * density(18) * density(36) 
  rrt(500) = rrt(500) * density(19) * density(36) 
  rrt(501) = rrt(501) * density(20) * density(36) 
  rrt(502) = rrt(502) * density(34) * density(36) 
  rrt(503) = rrt(503) * density(35) * density(36) 
  rrt(504) = rrt(504) * density(36) * density(45) 
  rrt(505) = rrt(505) * density(36) * density(46) 
  rrt(506) = rrt(506) * density(36) * density(47) 
  rrt(507) = rrt(507) * density(36) * density(52) 
  rrt(508) = rrt(508) * density(18) * density(37) 
  rrt(509) = rrt(509) * density(19) * density(37) 
  rrt(510) = rrt(510) * density(20) * density(37) 
  rrt(511) = rrt(511) * density(34) * density(37) 
  rrt(512) = rrt(512) * density(35) * density(37) 
  rrt(513) = rrt(513) * density(37) * density(45) 
  rrt(514) = rrt(514) * density(37) * density(46) 
  rrt(515) = rrt(515) * density(37) * density(47) 
  rrt(516) = rrt(516) * density(37) * density(52) 
  rrt(517) = rrt(517) * density(18) * density(38) 
  rrt(518) = rrt(518) * density(19) * density(38) 
  rrt(519) = rrt(519) * density(20) * density(38) 
  rrt(520) = rrt(520) * density(34) * density(38) 
  rrt(521) = rrt(521) * density(35) * density(38) 
  rrt(522) = rrt(522) * density(38) * density(45) 
  rrt(523) = rrt(523) * density(38) * density(46) 
  rrt(524) = rrt(524) * density(38) * density(47) 
  rrt(525) = rrt(525) * density(38) * density(52) 
  rrt(526) = rrt(526) * density(18) * density(48) 
  rrt(527) = rrt(527) * density(19) * density(48) 
  rrt(528) = rrt(528) * density(20) * density(48) 
  rrt(529) = rrt(529) * density(34) * density(48) 
  rrt(530) = rrt(530) * density(35) * density(48) 
  rrt(531) = rrt(531) * density(45) * density(48) 
  rrt(532) = rrt(532) * density(46) * density(48) 
  rrt(533) = rrt(533) * density(47) * density(48) 
  rrt(534) = rrt(534) * density(48) * density(52) 
  rrt(535) = rrt(535) * density(18) * density(49) 
  rrt(536) = rrt(536) * density(19) * density(49) 
  rrt(537) = rrt(537) * density(20) * density(49) 
  rrt(538) = rrt(538) * density(34) * density(49) 
  rrt(539) = rrt(539) * density(35) * density(49) 
  rrt(540) = rrt(540) * density(45) * density(49) 
  rrt(541) = rrt(541) * density(46) * density(49) 
  rrt(542) = rrt(542) * density(47) * density(49) 
  rrt(543) = rrt(543) * density(49) * density(52) 
  rrt(544) = rrt(544) * density(18) * density(50) 
  rrt(545) = rrt(545) * density(19) * density(50) 
  rrt(546) = rrt(546) * density(20) * density(50) 
  rrt(547) = rrt(547) * density(34) * density(50) 
  rrt(548) = rrt(548) * density(35) * density(50) 
  rrt(549) = rrt(549) * density(45) * density(50) 
  rrt(550) = rrt(550) * density(46) * density(50) 
  rrt(551) = rrt(551) * density(47) * density(50) 
  rrt(552) = rrt(552) * density(50) * density(52) 
  rrt(553) = rrt(553) * density(18) * density(51) 
  rrt(554) = rrt(554) * density(19) * density(51) 
  rrt(555) = rrt(555) * density(20) * density(51) 
  rrt(556) = rrt(556) * density(34) * density(51) 
  rrt(557) = rrt(557) * density(35) * density(51) 
  rrt(558) = rrt(558) * density(45) * density(51) 
  rrt(559) = rrt(559) * density(46) * density(51) 
  rrt(560) = rrt(560) * density(47) * density(51) 
  rrt(561) = rrt(561) * density(51) * density(52) 
  rrt(562) = rrt(562) * density(17) * density(39) 
  rrt(563) = rrt(563) * density(18) * density(39) 
  rrt(564) = rrt(564) * density(33) * density(39) 
  rrt(565) = rrt(565) * density(34) * density(39) 
  rrt(566) = rrt(566) * density(39) * density(45) 
  rrt(567) = rrt(567) * density(39) * density(46) 
  rrt(568) = rrt(568) * density(39) * density(47) 
  rrt(569) = rrt(569) * density(19) * density(39) 
  rrt(570) = rrt(570) * density(20) * density(39) 
  rrt(571) = rrt(571) * density(35) * density(39) 
  rrt(572) = rrt(572) * density(39) * density(52) 
  rrt(573) = rrt(573) * density(17) * density(36) 
  rrt(574) = rrt(574) * density(18) * density(36) 
  rrt(575) = rrt(575) * density(33) * density(36) 
  rrt(576) = rrt(576) * density(34) * density(36) 
  rrt(577) = rrt(577) * density(36) * density(45) 
  rrt(578) = rrt(578) * density(17) * density(37) 
  rrt(579) = rrt(579) * density(18) * density(37) 
  rrt(580) = rrt(580) * density(33) * density(37) 
  rrt(581) = rrt(581) * density(34) * density(37) 
  rrt(582) = rrt(582) * density(37) * density(45) 
  rrt(583) = rrt(583) * density(17) * density(36) 
  rrt(584) = rrt(584) * density(18) * density(36) 
  rrt(585) = rrt(585) * density(33) * density(36) 
  rrt(586) = rrt(586) * density(34) * density(36) 
  rrt(587) = rrt(587) * density(36) * density(45) 
  rrt(588) = rrt(588) * density(17) * density(37) 
  rrt(589) = rrt(589) * density(33) * density(37) 
  rrt(590) = rrt(590) * density(37) * density(45) 
  rrt(591) = rrt(591) * density(17) * density(38) 
  rrt(592) = rrt(592) * density(18) * density(38) 
  rrt(593) = rrt(593) * density(33) * density(38) 
  rrt(594) = rrt(594) * density(34) * density(38) 
  rrt(595) = rrt(595) * density(38) * density(45) 
  rrt(596) = rrt(596) * density(38) * density(46) 
  rrt(597) = rrt(597) * density(38) * density(47) 
  rrt(598) = rrt(598) * density(17) * density(48) 
  rrt(599) = rrt(599) * density(18) * density(48) 
  rrt(600) = rrt(600) * density(33) * density(48) 
  rrt(601) = rrt(601) * density(34) * density(48) 
  rrt(602) = rrt(602) * density(45) * density(48) 
  rrt(603) = rrt(603) * density(46) * density(48) 
  rrt(604) = rrt(604) * density(47) * density(48) 
  rrt(605) = rrt(605) * density(17) * density(49) 
  rrt(606) = rrt(606) * density(18) * density(49) 
  rrt(607) = rrt(607) * density(33) * density(49) 
  rrt(608) = rrt(608) * density(34) * density(49) 
  rrt(609) = rrt(609) * density(45) * density(49) 
  rrt(610) = rrt(610) * density(46) * density(49) 
  rrt(611) = rrt(611) * density(47) * density(49) 
  rrt(612) = rrt(612) * density(17) * density(50) 
  rrt(613) = rrt(613) * density(18) * density(50) 
  rrt(614) = rrt(614) * density(33) * density(50) 
  rrt(615) = rrt(615) * density(34) * density(50) 
  rrt(616) = rrt(616) * density(45) * density(50) 
  rrt(617) = rrt(617) * density(46) * density(50) 
  rrt(618) = rrt(618) * density(47) * density(50) 
  rrt(619) = rrt(619) * density(17) * density(51) 
  rrt(620) = rrt(620) * density(18) * density(51) 
  rrt(621) = rrt(621) * density(33) * density(51) 
  rrt(622) = rrt(622) * density(34) * density(51) 
  rrt(623) = rrt(623) * density(45) * density(51) 
  rrt(624) = rrt(624) * density(46) * density(51) 
  rrt(625) = rrt(625) * density(47) * density(51) 
  ydot(01) = +rrt(001)-rrt(009)+rrt(017)-rrt(025)+rrt(033)-rrt(041)-rrt(065)-rrt(066)-rrt(067)-rrt(068)-rrt(069)-rrt(070)-rrt(071)&
             -rrt(072)+rrt(082)-rrt(086)+rrt(100)+  2.d0 * rrt(101)+rrt(102)+rrt(105)-rrt(126)+rrt(130)+rrt(131)+rrt(139)+rrt(140)&
             +rrt(144)+rrt(151)+rrt(152)+rrt(153)+rrt(154)+rrt(155)+rrt(156)+rrt(157)+rrt(158)+rrt(159)+rrt(160)+rrt(161)+rrt(163)&
             +rrt(170)+rrt(171)+rrt(172)+rrt(173)+rrt(174)+rrt(175)+rrt(177)+rrt(178)+rrt(179)+rrt(180)+rrt(181)+rrt(182)+rrt(184)&
             +rrt(185)+rrt(188)+rrt(190)+rrt(191)+rrt(206)+rrt(207)+rrt(242)+rrt(257)+rrt(259)+rrt(261)-rrt(263)+rrt(266)-rrt(270)&
             +rrt(273)+rrt(276)-rrt(288)-rrt(289)-rrt(290)-rrt(291)-rrt(292)+rrt(307)+rrt(308)+rrt(309)+rrt(310)+rrt(326)+rrt(327)&
             +rrt(328)+rrt(329)+rrt(330)-rrt(346)+rrt(363)+rrt(364)-rrt(365)+rrt(373)+rrt(375)+rrt(377)+rrt(378)+rrt(379)+rrt(380)&
             +rrt(381)-rrt(382)+rrt(387)+rrt(388)+rrt(389)+rrt(390)+rrt(391)+rrt(394)+  2.d0 * rrt(395)+  2.d0 * rrt(396)&
             +  2.d0 * rrt(397)+  2.d0 * rrt(398)-rrt(399)+rrt(405)+rrt(406)-rrt(407)-rrt(410)-rrt(413)-rrt(416)+rrt(434)+rrt(451)&
             +rrt(458)+rrt(465)+rrt(472)+rrt(479)+rrt(486)+rrt(493)+rrt(500)+  2.d0 * rrt(501)+rrt(505)+rrt(507)+rrt(509)&
             +  2.d0 * rrt(510)+rrt(514)+rrt(516)+rrt(518)+  2.d0 * rrt(519)+rrt(523)+rrt(525)+rrt(527)+  2.d0 * rrt(528)+rrt(532)&
             +rrt(534)+rrt(536)+  2.d0 * rrt(537)+rrt(541)+rrt(543)+rrt(545)+  2.d0 * rrt(546)+rrt(550)+rrt(552)+rrt(554)&
             +  2.d0 * rrt(555)+rrt(559)+rrt(561)+rrt(563)+rrt(569)+  2.d0 * rrt(570)+rrt(572)+rrt(574)+rrt(579)+rrt(592)+rrt(599)&
             +rrt(606)+rrt(613)+rrt(620) 
  ydot(02) = -rrt(001)+rrt(002)+rrt(009)-rrt(010)-rrt(017)+rrt(018)+rrt(025)-rrt(026)-rrt(033)+rrt(034)+rrt(041)-rrt(042) 
  ydot(03) = -rrt(002)+rrt(003)+rrt(010)-rrt(011)-rrt(018)+rrt(019)+rrt(026)-rrt(027)-rrt(034)+rrt(035)+rrt(042)-rrt(043) 
  ydot(04) = -rrt(003)+rrt(004)+rrt(011)-rrt(012)-rrt(019)+rrt(020)+rrt(027)-rrt(028)-rrt(035)+rrt(036)+rrt(043)-rrt(044) 
  ydot(05) = -rrt(004)+rrt(005)+rrt(012)-rrt(013)-rrt(020)+rrt(021)+rrt(028)-rrt(029)-rrt(036)+rrt(037)+rrt(044)-rrt(045) 
  ydot(06) = -rrt(005)+rrt(006)+rrt(013)-rrt(014)-rrt(021)+rrt(022)+rrt(029)-rrt(030)-rrt(037)+rrt(038)+rrt(045)-rrt(046) 
  ydot(07) = -rrt(006)+rrt(007)+rrt(014)-rrt(015)-rrt(022)+rrt(023)+rrt(030)-rrt(031)-rrt(038)+rrt(039)+rrt(046)-rrt(047) 
  ydot(08) = -rrt(007)+rrt(008)+rrt(015)-rrt(016)-rrt(023)+rrt(024)+rrt(031)-rrt(032)-rrt(039)+rrt(040)+rrt(047)-rrt(048) 
  ydot(09) = -rrt(008)+rrt(016)-rrt(024)+rrt(032)-rrt(040)+rrt(048) 
  ydot(10) = +rrt(065)-rrt(082)-rrt(087)-rrt(130)-rrt(139)-rrt(151)-rrt(152)-rrt(153)-rrt(154)-rrt(155)-rrt(161)+rrt(162)-rrt(169)&
             -rrt(170)-rrt(171)-rrt(172)-rrt(173)-rrt(174)-rrt(175)-rrt(176)-rrt(177)-rrt(178)-rrt(179)-rrt(180)-  2.d0 * rrt(181)&
             -  2.d0 * rrt(182)+rrt(183)+rrt(186)-rrt(192)+rrt(194)+rrt(195)+rrt(196)+rrt(197)+rrt(198)+rrt(215) 
  ydot(11) = -rrt(131)-rrt(140)-rrt(156)-rrt(157)-rrt(158)-rrt(159)-rrt(160)-rrt(162)+rrt(164)+rrt(181)-rrt(183)-rrt(184)-rrt(185)&
             -rrt(186)+rrt(189)+rrt(199)+rrt(200)+rrt(201)+rrt(202)+rrt(203) 
  ydot(12) = +rrt(066)+rrt(067)+rrt(068)-rrt(163)+rrt(187)-rrt(189)-rrt(190)-rrt(191)-rrt(192)-  2.d0 * rrt(193) 
  ydot(13) = +rrt(069)+rrt(070)+rrt(071)-rrt(164)+rrt(182)-rrt(187)-rrt(188) 
  ydot(14) = +rrt(072)-rrt(084)+  2.d0 * rrt(092)+rrt(093)+rrt(094)+rrt(098)+rrt(100)+rrt(106)+rrt(108)+rrt(111)+rrt(114)-rrt(124)&
             -rrt(134)-rrt(142)-rrt(143)-rrt(144)-rrt(145)-rrt(146)-rrt(172)+rrt(179)+rrt(191)-  2.d0 * rrt(194)-  2.d0 * rrt(195)&
             -  2.d0 * rrt(196)-  2.d0 * rrt(197)-  2.d0 * rrt(198)-  2.d0 * rrt(199)-  2.d0 * rrt(200)-  2.d0 * rrt(201)&
             -  2.d0 * rrt(202)-  2.d0 * rrt(203)+rrt(204)+rrt(208)+rrt(209)+rrt(210)+rrt(212)-rrt(217)+rrt(240)-rrt(257)-rrt(258)&
             -rrt(259)-rrt(260)-rrt(261)-rrt(262)+rrt(263)+rrt(264)+rrt(271)-  2.d0 * rrt(286)-rrt(287)+  2.d0 * rrt(288)&
             +  2.d0 * rrt(289)+  2.d0 * rrt(290)+  2.d0 * rrt(291)+  2.d0 * rrt(292)+rrt(298)+rrt(299)+rrt(300)+rrt(301)+rrt(302)&
             -  2.d0 * rrt(326)-  2.d0 * rrt(327)-  2.d0 * rrt(328)-  2.d0 * rrt(329)-  2.d0 * rrt(330)-rrt(336)-rrt(337)-rrt(338)&
             -rrt(339)-rrt(340)+rrt(356)+rrt(357)+rrt(361)+rrt(365)+rrt(369)+rrt(376)-rrt(378)+rrt(381)-rrt(383)+rrt(387)-rrt(389)&
             +rrt(390)-rrt(397)-rrt(409)+rrt(410)-rrt(412)-rrt(414)+rrt(450)+rrt(457)+rrt(464)+rrt(471)+rrt(478)+rrt(485)+rrt(492)&
             +  2.d0 * rrt(499)+rrt(500)+rrt(504)+rrt(506)+  2.d0 * rrt(508)+rrt(509)+rrt(513)+rrt(515)+  2.d0 * rrt(517)+rrt(518)&
             +rrt(522)+rrt(524)+  2.d0 * rrt(526)+rrt(527)+rrt(531)+rrt(533)+  2.d0 * rrt(535)+rrt(536)+rrt(540)+rrt(542)&
             +  2.d0 * rrt(544)+rrt(545)+rrt(549)+rrt(551)+  2.d0 * rrt(553)+rrt(554)+rrt(558)+rrt(560)+rrt(562)+rrt(569)+rrt(573)&
             +rrt(578)+rrt(591)+rrt(598)+rrt(605)+rrt(612)+rrt(619) 
  ydot(15) = +rrt(072)+rrt(093)+rrt(099)+rrt(169)-rrt(204)-rrt(205)-rrt(206)-rrt(207)-rrt(208)+rrt(211)-rrt(213)-rrt(370) 
  ydot(16) = +rrt(094)+rrt(172)-rrt(209)-rrt(210)-rrt(211)-rrt(212)-rrt(213)-rrt(214)-rrt(215) 
  ydot(17) = +rrt(084)-rrt(106)-rrt(108)-rrt(356)-rrt(357)-rrt(358)-rrt(359)-rrt(360)-rrt(361)-rrt(362)-rrt(363)-rrt(364)+rrt(370)&
             +rrt(378)+rrt(397)-rrt(407)-rrt(408)-rrt(409)-rrt(450)-rrt(457)-rrt(464)-rrt(471)-rrt(478)-rrt(485)-rrt(492)-rrt(562)&
             -rrt(573)-rrt(578)-rrt(583)-rrt(588)-rrt(591)-rrt(598)-rrt(605)-rrt(612)-rrt(619) 
  ydot(18) = +rrt(086)+rrt(087)-rrt(092)-rrt(093)-rrt(094)+rrt(213)+rrt(286)+rrt(362)-rrt(375)-rrt(376)-rrt(377)-rrt(378)-rrt(379)&
             -rrt(380)-rrt(381)+rrt(389)+rrt(394)+rrt(409)-rrt(413)-rrt(414)-rrt(451)-rrt(458)-rrt(465)-rrt(472)-rrt(479)-rrt(486)&
             -rrt(493)-rrt(499)-rrt(508)-rrt(517)-rrt(526)-rrt(535)-rrt(544)-rrt(553)-rrt(563)-rrt(574)-rrt(579)-rrt(584)-rrt(592)&
             -rrt(599)-rrt(606)-rrt(613)-rrt(620) 
  ydot(19) = -rrt(100)-rrt(387)-rrt(388)-rrt(389)-rrt(390)-rrt(391)+rrt(407)+rrt(414)-rrt(500)-rrt(509)-rrt(518)-rrt(527)-rrt(536)&
             -rrt(545)-rrt(554)-rrt(569) 
  ydot(20) = -rrt(101)+rrt(192)+rrt(193)-rrt(394)-rrt(395)-rrt(396)-rrt(397)-rrt(398)+rrt(413)-rrt(501)-rrt(510)-rrt(519)-rrt(528)&
             -rrt(537)-rrt(546)-rrt(555)-rrt(570) 
  ydot(21) = +rrt(049)-rrt(053)+rrt(057)-rrt(061)-rrt(073)-rrt(074)-rrt(075)-rrt(076)-rrt(077)-rrt(078)+rrt(083)-rrt(088)&
             +  2.d0 * rrt(104)+rrt(105)-rrt(110)+rrt(112)-rrt(115)-rrt(118)-rrt(122)+rrt(123)-rrt(127)+rrt(129)+  2.d0 * rrt(132)&
             +rrt(135)+  2.d0 * rrt(136)+  2.d0 * rrt(137)+rrt(138)+rrt(139)+rrt(140)+  2.d0 * rrt(141)+rrt(143)+rrt(149)+rrt(165)&
             +rrt(167)+rrt(168)-rrt(173)-rrt(174)-rrt(175)-rrt(176)-rrt(185)-rrt(188)-rrt(190)-rrt(205)-rrt(214)+rrt(216)+rrt(218)&
             +rrt(219)+rrt(220)+  2.d0 * rrt(221)+rrt(222)+rrt(223)+rrt(225)+  2.d0 * rrt(229)+rrt(230)-rrt(231)-rrt(235)-rrt(236)&
             +rrt(238)+  2.d0 * rrt(239)+rrt(240)+rrt(242)-rrt(246)+  2.d0 * rrt(253)+rrt(254)-rrt(258)+rrt(261)+rrt(264)+rrt(266)&
             +rrt(268)+rrt(269)-rrt(270)+rrt(273)-rrt(274)+rrt(275)-  2.d0 * rrt(278)-rrt(279)+rrt(280)+rrt(282)+rrt(283)-rrt(284)&
             +rrt(285)-rrt(293)-rrt(294)-rrt(295)-rrt(296)-rrt(297)+rrt(303)+rrt(304)+rrt(305)+rrt(306)+rrt(320)+rrt(321)+rrt(322)&
             +rrt(323)+rrt(324)+rrt(331)+rrt(332)+rrt(333)+rrt(334)+rrt(335)-rrt(341)-rrt(342)-rrt(343)-rrt(344)-rrt(345)-rrt(357)&
             -rrt(358)-rrt(359)+rrt(360)-rrt(366)+rrt(367)-rrt(375)+rrt(384)+rrt(386)-rrt(387)-rrt(388)-rrt(395)+rrt(399)+rrt(400)&
             +  2.d0 * rrt(401)+  2.d0 * rrt(402)+  2.d0 * rrt(404)-rrt(406)-rrt(415)+rrt(422)+rrt(423)+rrt(424)+rrt(425)+rrt(426)&
             +rrt(428)+rrt(430)-rrt(432)+rrt(435)+rrt(440)+rrt(441)+rrt(442)+  2.d0 * rrt(443)+  2.d0 * rrt(444)+  2.d0 * rrt(445)&
             +rrt(446)-rrt(447)-rrt(449)+rrt(453)+rrt(457)+rrt(458)+rrt(459)+  2.d0 * rrt(460)+rrt(461)+rrt(462)+rrt(463)+rrt(467)&
             +rrt(474)+rrt(481)+rrt(488)+rrt(495)+  2.d0 * rrt(503)+rrt(506)+rrt(507)+rrt(508)+rrt(509)+rrt(510)+rrt(511)&
             +  3.d0 * rrt(512)+rrt(513)+rrt(514)+  2.d0 * rrt(515)+  2.d0 * rrt(516)+  2.d0 * rrt(521)+rrt(524)+rrt(525)&
             +  2.d0 * rrt(530)+rrt(533)+rrt(534)+  2.d0 * rrt(539)+rrt(542)+rrt(543)+  2.d0 * rrt(548)+rrt(551)+rrt(552)&
             +  2.d0 * rrt(557)+rrt(560)+rrt(561)+  2.d0 * rrt(562)+  2.d0 * rrt(563)+  2.d0 * rrt(564)+  3.d0 * rrt(565)&
             +  2.d0 * rrt(566)+  2.d0 * rrt(567)+  2.d0 * rrt(568)+  2.d0 * rrt(569)+  2.d0 * rrt(570)+  4.d0 * rrt(571)&
             +  3.d0 * rrt(572)+rrt(576)+rrt(578)+rrt(579)+rrt(580)+  2.d0 * rrt(581)+rrt(582)+rrt(585)+rrt(594)+rrt(601)+rrt(608)
  ydot(21) = ydot(21) &
             +rrt(615)+rrt(622) 
  ydot(22) = -rrt(049)+rrt(050)+rrt(053)-rrt(054)-rrt(057)+rrt(058)+rrt(061)-rrt(062) 
  ydot(23) = -rrt(050)+rrt(051)+rrt(054)-rrt(055)-rrt(058)+rrt(059)+rrt(062)-rrt(063) 
  ydot(24) = -rrt(051)+rrt(052)+rrt(055)-rrt(056)-rrt(059)+rrt(060)+rrt(063)-rrt(064) 
  ydot(25) = -rrt(052)+rrt(056)-rrt(060)+rrt(064) 
  ydot(26) = +rrt(073)-rrt(079)-rrt(083)-rrt(089)-rrt(128)-rrt(136)-rrt(165)+rrt(166)+rrt(174)-rrt(216)-rrt(217)-rrt(218)-rrt(219)&
             -rrt(220)-rrt(221)-  2.d0 * rrt(222)+rrt(223)+rrt(224)+rrt(226)+rrt(227)+rrt(228)+rrt(235)-rrt(248)-rrt(249)-rrt(250)&
             -rrt(401)-rrt(417)-rrt(444) 
  ydot(27) = +rrt(074)-rrt(129)-rrt(137)-rrt(166)-rrt(167)+rrt(175)+rrt(222)-rrt(224)-rrt(225)-rrt(226)-rrt(227)-rrt(228)-rrt(229)&
             +  2.d0 * rrt(231)+rrt(232)+rrt(236)+rrt(249)-rrt(402)-rrt(445) 
  ydot(28) = +rrt(075)-rrt(168)-rrt(230)-rrt(231)-rrt(232)+rrt(248) 
  ydot(29) = +  2.d0 * rrt(076)+rrt(077)+rrt(078)+  2.d0 * rrt(079)-rrt(080)-rrt(081)-rrt(085)+  2.d0 * rrt(095)+rrt(096)+rrt(097)&
             +rrt(098)+rrt(099)+rrt(102)+rrt(103)+rrt(107)+rrt(109)+rrt(110)+rrt(113)-rrt(117)-rrt(123)+rrt(129)+rrt(130)+rrt(131)&
             -rrt(133)-rrt(141)-rrt(147)-rrt(148)-rrt(149)-rrt(150)-rrt(169)-rrt(170)+rrt(173)+rrt(176)+rrt(180)+  2.d0 * rrt(185)&
             +rrt(188)+  2.d0 * rrt(190)+rrt(191)-rrt(204)+rrt(205)+rrt(206)+rrt(214)+rrt(215)+rrt(217)-rrt(223)-rrt(225)+rrt(229)&
             -rrt(230)+rrt(233)+rrt(234)+rrt(235)+rrt(236)+rrt(237)+  2.d0 * rrt(238)+rrt(244)+  3.d0 * rrt(246)+rrt(247)+rrt(248)&
             +  3.d0 * rrt(250)+rrt(251)+rrt(254)+rrt(255)+rrt(257)+rrt(258)+  2.d0 * rrt(259)+rrt(260)-rrt(263)-rrt(264)-rrt(265)&
             -rrt(266)-rrt(267)-rrt(268)-rrt(269)+rrt(270)+rrt(272)+rrt(274)+rrt(278)-rrt(287)+  2.d0 * rrt(293)+  2.d0 * rrt(294)&
             +  2.d0 * rrt(295)+  2.d0 * rrt(296)+  2.d0 * rrt(297)+rrt(298)+rrt(299)+rrt(300)+rrt(301)+rrt(302)+rrt(303)+rrt(304)&
             +rrt(305)+rrt(306)+rrt(307)+rrt(308)+rrt(309)+rrt(310)+rrt(311)+rrt(312)+rrt(313)+rrt(314)+rrt(315)+rrt(316)+rrt(317)&
             +rrt(318)+rrt(319)-  2.d0 * rrt(331)-  2.d0 * rrt(332)-  2.d0 * rrt(333)-  2.d0 * rrt(334)-  2.d0 * rrt(335)-rrt(336)&
             -rrt(337)-rrt(338)-rrt(339)-rrt(340)-rrt(341)-rrt(342)-rrt(343)-rrt(344)-rrt(345)-rrt(346)-rrt(347)-rrt(348)-rrt(349)&
             -rrt(350)-rrt(351)-rrt(352)-rrt(353)-rrt(354)-rrt(356)+rrt(358)+rrt(362)+rrt(366)+rrt(368)+rrt(370)+rrt(372)+rrt(374)&
             -rrt(376)+rrt(377)+rrt(383)-rrt(396)-rrt(403)-rrt(408)-rrt(411)+rrt(417)+rrt(418)+rrt(419)+rrt(421)-rrt(422)-rrt(426)&
             +rrt(427)-rrt(442)-rrt(443)+rrt(450)+rrt(451)+  2.d0 * rrt(452)+rrt(453)+rrt(454)+rrt(455)+rrt(456)+rrt(459)+rrt(466)&
             +rrt(473)+rrt(480)+rrt(487)+rrt(494)+rrt(499)+rrt(500)+rrt(501)+  3.d0 * rrt(502)+rrt(503)+  2.d0 * rrt(504)&
             +  2.d0 * rrt(505)+rrt(506)+rrt(507)+  2.d0 * rrt(511)+rrt(513)+rrt(514)+  2.d0 * rrt(520)+rrt(522)+rrt(523)&
             +  2.d0 * rrt(529)+rrt(531)+rrt(532)+  2.d0 * rrt(538)+rrt(540)+rrt(541)+  2.d0 * rrt(547)+rrt(549)+rrt(550)&
             +  2.d0 * rrt(556)+rrt(558)+rrt(559)+rrt(564)+rrt(573)+rrt(574)+  2.d0 * rrt(575)+rrt(576)+rrt(577)+rrt(580)+rrt(593)&
             +rrt(600)+rrt(607)+rrt(614)+rrt(621) 
  ydot(30) = +rrt(077)+rrt(080)+rrt(096)+rrt(173)+rrt(204)+rrt(221)+rrt(225)-rrt(233)-rrt(234)-rrt(235)-rrt(236)-rrt(237)-rrt(238)&
             -rrt(239)-rrt(240)-rrt(241)-rrt(242)+rrt(243)+rrt(245)+rrt(249)+rrt(252)+rrt(254)+rrt(256) 
  ydot(31) = +rrt(078)+rrt(081)+rrt(097)+rrt(170)+rrt(188)+rrt(230)-rrt(243)-rrt(244)-rrt(245)-rrt(246)-rrt(247)-rrt(248)-rrt(249)&
             -rrt(250)-rrt(251)-rrt(252)-rrt(253)-rrt(254)-rrt(255)-rrt(256) 
  ydot(32) = -rrt(112)-rrt(113)-rrt(119)+rrt(127)+rrt(128)-rrt(132)+rrt(133)+rrt(150)+rrt(151)+rrt(156)-rrt(221)-rrt(223)-rrt(229)&
             -rrt(238)-rrt(239)-rrt(253)-rrt(254)-rrt(275)+rrt(278)+rrt(279)-rrt(282)+rrt(284)-rrt(303)-rrt(304)-rrt(305)-rrt(306)&
             +rrt(341)+rrt(342)+rrt(343)+rrt(344)+rrt(345)-rrt(360)-rrt(367)-rrt(377)+rrt(385)+rrt(403)-rrt(418)-rrt(423)+rrt(429)&
             +rrt(431)-rrt(435)+rrt(464)+rrt(465)+rrt(466)+rrt(467)+rrt(468)+rrt(469)+rrt(470)+rrt(517)+rrt(518)+rrt(519)+rrt(520)&
             +rrt(521)+rrt(522)+rrt(523)+rrt(524)+rrt(525)+rrt(586)+rrt(589)+rrt(591)+rrt(592)+rrt(593)+rrt(594)+rrt(595)+rrt(596)&
             +rrt(597) 
  ydot(33) = +rrt(085)-rrt(107)-rrt(109)+rrt(356)+rrt(359)+rrt(363)-rrt(365)-rrt(366)-rrt(367)-rrt(368)-rrt(369)-rrt(370)-rrt(371)&
             -rrt(372)-rrt(373)-rrt(374)+rrt(396)-rrt(410)-rrt(411)-rrt(412)-rrt(452)-rrt(459)-rrt(466)-rrt(473)-rrt(480)-rrt(487)&
             -rrt(494)-rrt(564)-rrt(575)-rrt(580)-rrt(585)-rrt(589)-rrt(593)-rrt(600)-rrt(607)-rrt(614)-rrt(621) 
  ydot(34) = +rrt(088)+rrt(089)-rrt(095)-rrt(096)-rrt(097)+rrt(357)+rrt(366)+rrt(367)+rrt(369)+rrt(373)+rrt(375)+rrt(377)-rrt(382)&
             -rrt(383)-rrt(384)-rrt(385)-rrt(386)+rrt(387)+rrt(395)+rrt(400)+rrt(401)+rrt(402)+rrt(403)+rrt(405)+rrt(411)-rrt(415)&
             -rrt(416)-rrt(453)-rrt(460)-rrt(467)-rrt(474)-rrt(481)-rrt(488)-rrt(495)-rrt(502)-rrt(511)-rrt(520)-rrt(529)-rrt(538)&
             -rrt(547)-rrt(556)-rrt(565)-rrt(576)-rrt(581)-rrt(586)-rrt(594)-rrt(601)-rrt(608)-rrt(615)-rrt(622) 
  ydot(35) = -rrt(104)-rrt(399)-rrt(400)-rrt(401)-rrt(402)-rrt(403)-rrt(404)+rrt(406)+rrt(415)-rrt(503)-rrt(512)-rrt(521)-rrt(530)&
             -rrt(539)-rrt(548)-rrt(557)-rrt(571) 
  ydot(36) = +rrt(110)+rrt(111)+rrt(112)+rrt(116)+rrt(117)-rrt(123)-rrt(124)-rrt(125)-rrt(126)-rrt(127)-rrt(128)-rrt(129)-rrt(130)&
             -rrt(131)-rrt(132)-rrt(417)-rrt(418)-rrt(419)-rrt(420)-rrt(421)+rrt(422)+rrt(443)-rrt(447)-rrt(448)-rrt(450)-rrt(451)&
             -rrt(452)-rrt(453)-rrt(454)-rrt(455)-rrt(456)-rrt(499)-rrt(500)-rrt(501)-rrt(502)-rrt(503)-rrt(504)-rrt(505)-rrt(506)&
             -rrt(507)-rrt(573)-rrt(574)-rrt(575)-rrt(576)-rrt(577)-rrt(583)-rrt(584)-rrt(585)-rrt(586)-rrt(587) 
  ydot(37) = +rrt(113)+rrt(115)+rrt(118)+rrt(122)-rrt(133)-rrt(134)-rrt(135)-rrt(136)-rrt(137)-rrt(138)-rrt(139)-rrt(140)+rrt(417)&
             -rrt(422)-rrt(423)-rrt(424)-rrt(425)+rrt(426)+rrt(432)+rrt(440)+rrt(441)+rrt(444)+rrt(445)-rrt(449)-rrt(457)-rrt(458)&
             -rrt(459)-rrt(460)-rrt(461)-rrt(462)-rrt(463)-rrt(508)-rrt(509)-rrt(510)-rrt(511)-rrt(512)-rrt(513)-rrt(514)-rrt(515)&
             -rrt(516)-rrt(578)-rrt(579)-rrt(580)-rrt(581)-rrt(582)-rrt(588)-rrt(589)-rrt(590) 
  ydot(38) = +rrt(119)-rrt(141)-rrt(143)-rrt(151)-rrt(156)+rrt(418)+rrt(423)-rrt(426)-rrt(427)-rrt(428)-rrt(429)-rrt(430)-rrt(431)&
             +rrt(442)+rrt(447)-rrt(464)-rrt(465)-rrt(466)-rrt(467)-rrt(468)-rrt(469)-rrt(470)-rrt(517)-rrt(518)-rrt(519)-rrt(520)&
             -rrt(521)-rrt(522)-rrt(523)-rrt(524)-rrt(525)-rrt(591)-rrt(592)-rrt(593)-rrt(594)-rrt(595)-rrt(596)-rrt(597) 
  ydot(39) = -rrt(440)-rrt(441)-rrt(442)-rrt(443)-rrt(444)-rrt(445)-rrt(446)+rrt(449)-rrt(562)-rrt(563)-rrt(564)-rrt(565)-rrt(566)&
             -rrt(567)-rrt(568)-rrt(569)-rrt(570)-rrt(571)-rrt(572) 
  ydot(40) = -rrt(090)+rrt(103)-rrt(111)+rrt(116)-rrt(120)+rrt(124)-rrt(125)+rrt(143)+rrt(144)+  2.d0 * rrt(145)+rrt(146)&
             +  2.d0 * rrt(148)+rrt(149)+rrt(150)+rrt(152)+rrt(157)+rrt(169)+rrt(179)+rrt(180)-rrt(191)+rrt(205)-rrt(206)+rrt(207)&
             +rrt(214)-rrt(215)+rrt(217)-rrt(240)+  2.d0 * rrt(241)-rrt(257)+rrt(258)+  2.d0 * rrt(262)+rrt(263)-rrt(264)-rrt(265)&
             +  2.d0 * rrt(267)+rrt(268)-  2.d0 * rrt(271)-  2.d0 * rrt(272)-  2.d0 * rrt(273)-rrt(274)-rrt(275)-rrt(276)-rrt(277)&
             +rrt(279)+  2.d0 * rrt(280)+rrt(281)+rrt(283)-rrt(298)-rrt(299)-rrt(300)-rrt(301)-rrt(302)+rrt(311)+rrt(312)+rrt(313)&
             +rrt(314)+rrt(320)+rrt(321)+rrt(322)+rrt(323)+rrt(324)+rrt(336)+rrt(337)+rrt(338)+rrt(339)+rrt(340)-rrt(347)-rrt(348)&
             -rrt(349)+rrt(359)-rrt(361)-rrt(362)-rrt(363)-rrt(368)-rrt(369)+rrt(371)-rrt(379)+rrt(382)-rrt(384)-rrt(390)-rrt(391)&
             -rrt(392)-rrt(393)-rrt(398)-rrt(404)+rrt(420)-rrt(427)-rrt(428)+rrt(432)+rrt(433)+rrt(436)-rrt(439)-rrt(446)-rrt(448)&
             +rrt(454)+rrt(461)+rrt(468)+rrt(471)+rrt(472)+rrt(473)+rrt(474)+  2.d0 * rrt(475)+rrt(476)+rrt(477)+rrt(482)+rrt(489)&
             +rrt(496)+rrt(526)+rrt(527)+rrt(528)+rrt(529)+rrt(530)+rrt(531)+rrt(532)+rrt(533)+rrt(534)+rrt(566)+rrt(577)+rrt(582)&
             +rrt(583)+rrt(595)+rrt(598)+rrt(599)+rrt(600)+rrt(601)+  2.d0 * rrt(602)+rrt(603)+rrt(604)+rrt(609)+rrt(616)+rrt(623) 
  ydot(41) = -rrt(091)-rrt(114)-rrt(121)+rrt(126)+rrt(142)+rrt(153)+rrt(158)+rrt(176)-rrt(179)-rrt(207)-rrt(241)-rrt(242)+rrt(260)&
             -rrt(266)-rrt(267)+rrt(270)+rrt(272)-rrt(276)-rrt(307)-rrt(308)-rrt(309)-rrt(310)+rrt(346)-rrt(364)-rrt(371)-rrt(372)&
             -rrt(373)-rrt(380)-rrt(381)+rrt(393)-rrt(420)-rrt(421)-rrt(434)+rrt(455)+rrt(462)+rrt(469)+rrt(476)+rrt(478)+rrt(479)&
             +rrt(480)+rrt(481)+rrt(482)+  2.d0 * rrt(483)+rrt(484)+rrt(490)+rrt(497)+rrt(535)+rrt(536)+rrt(537)+rrt(538)+rrt(539)&
             +rrt(540)+rrt(541)+rrt(542)+rrt(543)+rrt(567)+rrt(584)+rrt(596)+rrt(603)+rrt(605)+rrt(606)+rrt(607)+rrt(608)+rrt(609)&
             +  2.d0 * rrt(610)+rrt(611)+rrt(617)+rrt(624) 
  ydot(42) = -rrt(116)+rrt(125)+rrt(134)+rrt(146)+rrt(147)+rrt(154)+rrt(159)-rrt(180)-rrt(259)-rrt(260)-rrt(261)-rrt(262)+rrt(265)&
             -rrt(268)+rrt(269)+rrt(271)+rrt(274)+rrt(275)+rrt(276)+  2.d0 * rrt(277)-rrt(279)-  2.d0 * rrt(280)-  2.d0 * rrt(281)&
             -rrt(282)+rrt(284)+  2.d0 * rrt(285)-rrt(311)-rrt(312)-rrt(313)-rrt(314)+rrt(315)+rrt(316)+rrt(317)+rrt(318)+rrt(319)&
             +rrt(325)+rrt(347)+rrt(348)+rrt(349)-rrt(350)-rrt(351)-rrt(352)-rrt(353)-rrt(354)-rrt(355)-rrt(374)-rrt(385)-rrt(386)&
             +rrt(392)-rrt(419)-rrt(424)-rrt(429)-rrt(430)-rrt(433)-rrt(436)+rrt(437)+  2.d0 * rrt(438)+rrt(439)+rrt(456)+rrt(463)&
             +rrt(470)+rrt(477)+rrt(484)+rrt(485)+rrt(486)+rrt(487)+rrt(488)+rrt(489)+rrt(490)+  2.d0 * rrt(491)+rrt(498)+rrt(544)&
             +rrt(545)+rrt(546)+rrt(547)+rrt(548)+rrt(549)+rrt(550)+rrt(551)+rrt(552)+rrt(568)+rrt(587)+rrt(588)+rrt(597)+rrt(604)&
             +rrt(611)+rrt(612)+rrt(613)+rrt(614)+rrt(615)+rrt(616)+rrt(617)+  2.d0 * rrt(618)+rrt(625) 
  ydot(43) = +rrt(155)+rrt(160)-rrt(269)-rrt(277)+rrt(281)+rrt(282)-rrt(283)-rrt(284)-  2.d0 * rrt(285)-rrt(315)-rrt(316)-rrt(317)&
             -rrt(318)-rrt(319)-rrt(320)-rrt(321)-rrt(322)-rrt(323)-rrt(324)+rrt(325)+rrt(350)+rrt(351)+rrt(352)+rrt(353)+rrt(354)&
             -rrt(355)-rrt(425)-rrt(431)-rrt(437)+rrt(492)+rrt(493)+rrt(494)+rrt(495)+rrt(496)+rrt(497)+rrt(498)+rrt(553)+rrt(554)&
             +rrt(555)+rrt(556)+rrt(557)+rrt(558)+rrt(559)+rrt(560)+rrt(561)+rrt(590)+rrt(619)+rrt(620)+rrt(621)+rrt(622)+rrt(623)&
             +rrt(624)+rrt(625) 
  ydot(44) = -rrt(325)+rrt(355)-rrt(438) 
  ydot(45) = +rrt(090)-rrt(098)-rrt(099)+rrt(287)+rrt(358)+rrt(360)+rrt(361)+rrt(364)+rrt(365)+rrt(368)+rrt(371)+rrt(376)+rrt(379)&
             +rrt(381)+rrt(382)+rrt(383)+rrt(384)+rrt(385)+rrt(390)+rrt(392)+rrt(393)+rrt(398)+rrt(404)+rrt(408)+rrt(410)+rrt(412)&
             -rrt(454)-rrt(461)-rrt(468)-rrt(475)-rrt(482)-rrt(489)-rrt(496)-rrt(504)-rrt(513)-rrt(522)-rrt(531)-rrt(540)-rrt(549)&
             -rrt(558)-rrt(566)-rrt(577)-rrt(582)-rrt(587)-rrt(590)-rrt(595)-rrt(602)-rrt(609)-rrt(616)-rrt(623) 
  ydot(46) = +rrt(091)-rrt(102)+rrt(372)+rrt(380)+rrt(391)-rrt(393)-rrt(455)-rrt(462)-rrt(469)-rrt(476)-rrt(483)-rrt(490)-rrt(497)&
             -rrt(505)-rrt(514)-rrt(523)-rrt(532)-rrt(541)-rrt(550)-rrt(559)-rrt(567)-rrt(596)-rrt(603)-rrt(610)-rrt(617)-rrt(624) 
  ydot(47) = -rrt(103)+rrt(374)+rrt(386)+rrt(388)-rrt(392)-rrt(456)-rrt(463)-rrt(470)-rrt(477)-rrt(484)-rrt(491)-rrt(498)-rrt(506)&
             -rrt(515)-rrt(524)-rrt(533)-rrt(542)-rrt(551)-rrt(560)-rrt(568)-rrt(597)-rrt(604)-rrt(611)-rrt(618)-rrt(625) 
  ydot(48) = +rrt(114)+rrt(120)-rrt(142)-rrt(147)-rrt(152)-rrt(157)+rrt(420)-rrt(432)-rrt(433)-rrt(434)-rrt(471)-rrt(472)-rrt(473)&
             -rrt(474)-rrt(475)-rrt(476)-rrt(477)-rrt(526)-rrt(527)-rrt(528)-rrt(529)-rrt(530)-rrt(531)-rrt(532)-rrt(533)-rrt(534)&
             -rrt(598)-rrt(599)-rrt(600)-rrt(601)-rrt(602)-rrt(603)-rrt(604) 
  ydot(49) = +rrt(121)-rrt(144)-rrt(148)-rrt(153)-rrt(158)+rrt(421)-rrt(478)-rrt(479)-rrt(480)-rrt(481)-rrt(482)-rrt(483)-rrt(484)&
             -rrt(535)-rrt(536)-rrt(537)-rrt(538)-rrt(539)-rrt(540)-rrt(541)-rrt(542)-rrt(543)-rrt(605)-rrt(606)-rrt(607)-rrt(608)&
             -rrt(609)-rrt(610)-rrt(611) 
  ydot(50) = -rrt(145)-rrt(149)-rrt(154)-rrt(159)+rrt(419)+rrt(424)+rrt(428)+rrt(429)+rrt(433)+rrt(434)-rrt(435)-rrt(436)-rrt(437)&
             -rrt(438)+rrt(439)+rrt(448)-rrt(485)-rrt(486)-rrt(487)-rrt(488)-rrt(489)-rrt(490)-rrt(491)-rrt(544)-rrt(545)-rrt(546)&
             -rrt(547)-rrt(548)-rrt(549)-rrt(550)-rrt(551)-rrt(552)-rrt(612)-rrt(613)-rrt(614)-rrt(615)-rrt(616)-rrt(617)-rrt(618) 
  ydot(51) = -rrt(146)-rrt(150)-rrt(155)-rrt(160)+rrt(425)+rrt(427)+rrt(430)+rrt(431)+rrt(435)+rrt(436)+rrt(437)+rrt(438)-rrt(439)&
             +rrt(446)-rrt(492)-rrt(493)-rrt(494)-rrt(495)-rrt(496)-rrt(497)-rrt(498)-rrt(553)-rrt(554)-rrt(555)-rrt(556)-rrt(557)&
             -rrt(558)-rrt(559)-rrt(560)-rrt(561)-rrt(619)-rrt(620)-rrt(621)-rrt(622)-rrt(623)-rrt(624)-rrt(625) 
  ydot(52) = -rrt(105)+rrt(399)-rrt(405)-rrt(406)+rrt(416)-rrt(507)-rrt(516)-rrt(525)-rrt(534)-rrt(543)-rrt(552)-rrt(561)-rrt(572) 
  ydot(53) = +rrt(084)+rrt(085)+rrt(086)+rrt(087)+rrt(088)+rrt(089)+rrt(090)+rrt(091)-rrt(092)-rrt(093)-rrt(094)-rrt(095)-rrt(096)&
             -rrt(097)-rrt(098)-rrt(099)-rrt(100)-rrt(101)-rrt(102)-rrt(103)-rrt(104)-rrt(105)-rrt(106)-rrt(107)-rrt(108)-rrt(109)&
             -rrt(110)-rrt(111)-rrt(112)-rrt(113)-rrt(114)-rrt(115)-rrt(116)-rrt(117)-rrt(118)-rrt(119)-rrt(120)-rrt(121)-rrt(122)&
             +rrt(123)+rrt(124)+rrt(125)+rrt(126)+rrt(127)+rrt(128)+rrt(129)+rrt(130)+rrt(131)+rrt(132)+rrt(133)+rrt(134)+rrt(135)&
             +rrt(136)+rrt(137)+rrt(138)+rrt(139)+rrt(140)+rrt(141)+rrt(142)+rrt(143)+rrt(144)+rrt(145)+rrt(146)+rrt(147)+rrt(148)&
             +rrt(149)+rrt(150)+rrt(151)+rrt(152)+rrt(153)+rrt(154)+rrt(155)+rrt(156)+rrt(157)+rrt(158)+rrt(159)+rrt(160)+rrt(192)&
             +rrt(193)+rrt(213)+rrt(286)+rrt(287) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(54) = 0.0d0
  if( lgas_heating ) then
    ydot(54) = ( ZDPlasKin_cfg(14)/k_B + ydot(54) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(54) = ydot(54) * ZDPlasKin_cfg(13)
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
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(54)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(01,01) = pd(01,01) + rrt(001) * density(02) 
  pd(01,02) = pd(01,02) + rrt(001) * density(01) 
  pd(02,01) = pd(02,01) - rrt(001) * density(02) 
  pd(02,02) = pd(02,02) - rrt(001) * density(01) 
  pd(02,01) = pd(02,01) + rrt(002) * density(03) 
  pd(02,03) = pd(02,03) + rrt(002) * density(01) 
  pd(03,01) = pd(03,01) - rrt(002) * density(03) 
  pd(03,03) = pd(03,03) - rrt(002) * density(01) 
  pd(03,01) = pd(03,01) + rrt(003) * density(04) 
  pd(03,04) = pd(03,04) + rrt(003) * density(01) 
  pd(04,01) = pd(04,01) - rrt(003) * density(04) 
  pd(04,04) = pd(04,04) - rrt(003) * density(01) 
  pd(04,01) = pd(04,01) + rrt(004) * density(05) 
  pd(04,05) = pd(04,05) + rrt(004) * density(01) 
  pd(05,01) = pd(05,01) - rrt(004) * density(05) 
  pd(05,05) = pd(05,05) - rrt(004) * density(01) 
  pd(05,01) = pd(05,01) + rrt(005) * density(06) 
  pd(05,06) = pd(05,06) + rrt(005) * density(01) 
  pd(06,01) = pd(06,01) - rrt(005) * density(06) 
  pd(06,06) = pd(06,06) - rrt(005) * density(01) 
  pd(06,01) = pd(06,01) + rrt(006) * density(07) 
  pd(06,07) = pd(06,07) + rrt(006) * density(01) 
  pd(07,01) = pd(07,01) - rrt(006) * density(07) 
  pd(07,07) = pd(07,07) - rrt(006) * density(01) 
  pd(07,01) = pd(07,01) + rrt(007) * density(08) 
  pd(07,08) = pd(07,08) + rrt(007) * density(01) 
  pd(08,01) = pd(08,01) - rrt(007) * density(08) 
  pd(08,08) = pd(08,08) - rrt(007) * density(01) 
  pd(08,01) = pd(08,01) + rrt(008) * density(09) 
  pd(08,09) = pd(08,09) + rrt(008) * density(01) 
  pd(09,01) = pd(09,01) - rrt(008) * density(09) 
  pd(09,09) = pd(09,09) - rrt(008) * density(01) 
  pd(01,01) = pd(01,01) - rrt(009) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) + rrt(009) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) - rrt(010) * density(02) 
  pd(02,02) = pd(02,02) - rrt(010) * density(01) 
  pd(03,01) = pd(03,01) + rrt(010) * density(02) 
  pd(03,02) = pd(03,02) + rrt(010) * density(01) 
  pd(03,01) = pd(03,01) - rrt(011) * density(03) 
  pd(03,03) = pd(03,03) - rrt(011) * density(01) 
  pd(04,01) = pd(04,01) + rrt(011) * density(03) 
  pd(04,03) = pd(04,03) + rrt(011) * density(01) 
  pd(04,01) = pd(04,01) - rrt(012) * density(04) 
  pd(04,04) = pd(04,04) - rrt(012) * density(01) 
  pd(05,01) = pd(05,01) + rrt(012) * density(04) 
  pd(05,04) = pd(05,04) + rrt(012) * density(01) 
  pd(05,01) = pd(05,01) - rrt(013) * density(05) 
  pd(05,05) = pd(05,05) - rrt(013) * density(01) 
  pd(06,01) = pd(06,01) + rrt(013) * density(05) 
  pd(06,05) = pd(06,05) + rrt(013) * density(01) 
  pd(06,01) = pd(06,01) - rrt(014) * density(06) 
  pd(06,06) = pd(06,06) - rrt(014) * density(01) 
  pd(07,01) = pd(07,01) + rrt(014) * density(06) 
  pd(07,06) = pd(07,06) + rrt(014) * density(01) 
  pd(07,01) = pd(07,01) - rrt(015) * density(07) 
  pd(07,07) = pd(07,07) - rrt(015) * density(01) 
  pd(08,01) = pd(08,01) + rrt(015) * density(07) 
  pd(08,07) = pd(08,07) + rrt(015) * density(01) 
  pd(08,01) = pd(08,01) - rrt(016) * density(08) 
  pd(08,08) = pd(08,08) - rrt(016) * density(01) 
  pd(09,01) = pd(09,01) + rrt(016) * density(08) 
  pd(09,08) = pd(09,08) + rrt(016) * density(01) 
  pd(01,02) = pd(01,02) + rrt(017) * density(14) 
  pd(01,14) = pd(01,14) + rrt(017) * density(02) 
  pd(02,02) = pd(02,02) - rrt(017) * density(14) 
  pd(02,14) = pd(02,14) - rrt(017) * density(02) 
  pd(02,03) = pd(02,03) + rrt(018) * density(14) 
  pd(02,14) = pd(02,14) + rrt(018) * density(03) 
  pd(03,03) = pd(03,03) - rrt(018) * density(14) 
  pd(03,14) = pd(03,14) - rrt(018) * density(03) 
  pd(03,04) = pd(03,04) + rrt(019) * density(14) 
  pd(03,14) = pd(03,14) + rrt(019) * density(04) 
  pd(04,04) = pd(04,04) - rrt(019) * density(14) 
  pd(04,14) = pd(04,14) - rrt(019) * density(04) 
  pd(04,05) = pd(04,05) + rrt(020) * density(14) 
  pd(04,14) = pd(04,14) + rrt(020) * density(05) 
  pd(05,05) = pd(05,05) - rrt(020) * density(14) 
  pd(05,14) = pd(05,14) - rrt(020) * density(05) 
  pd(05,06) = pd(05,06) + rrt(021) * density(14) 
  pd(05,14) = pd(05,14) + rrt(021) * density(06) 
  pd(06,06) = pd(06,06) - rrt(021) * density(14) 
  pd(06,14) = pd(06,14) - rrt(021) * density(06) 
  pd(06,07) = pd(06,07) + rrt(022) * density(14) 
  pd(06,14) = pd(06,14) + rrt(022) * density(07) 
  pd(07,07) = pd(07,07) - rrt(022) * density(14) 
  pd(07,14) = pd(07,14) - rrt(022) * density(07) 
  pd(07,08) = pd(07,08) + rrt(023) * density(14) 
  pd(07,14) = pd(07,14) + rrt(023) * density(08) 
  pd(08,08) = pd(08,08) - rrt(023) * density(14) 
  pd(08,14) = pd(08,14) - rrt(023) * density(08) 
  pd(08,09) = pd(08,09) + rrt(024) * density(14) 
  pd(08,14) = pd(08,14) + rrt(024) * density(09) 
  pd(09,09) = pd(09,09) - rrt(024) * density(14) 
  pd(09,14) = pd(09,14) - rrt(024) * density(09) 
  pd(01,01) = pd(01,01) - rrt(025) * density(14) 
  pd(01,14) = pd(01,14) - rrt(025) * density(01) 
  pd(02,01) = pd(02,01) + rrt(025) * density(14) 
  pd(02,14) = pd(02,14) + rrt(025) * density(01) 
  pd(02,02) = pd(02,02) - rrt(026) * density(14) 
  pd(02,14) = pd(02,14) - rrt(026) * density(02) 
  pd(03,02) = pd(03,02) + rrt(026) * density(14) 
  pd(03,14) = pd(03,14) + rrt(026) * density(02) 
  pd(03,03) = pd(03,03) - rrt(027) * density(14) 
  pd(03,14) = pd(03,14) - rrt(027) * density(03) 
  pd(04,03) = pd(04,03) + rrt(027) * density(14) 
  pd(04,14) = pd(04,14) + rrt(027) * density(03) 
  pd(04,04) = pd(04,04) - rrt(028) * density(14) 
  pd(04,14) = pd(04,14) - rrt(028) * density(04) 
  pd(05,04) = pd(05,04) + rrt(028) * density(14) 
  pd(05,14) = pd(05,14) + rrt(028) * density(04) 
  pd(05,05) = pd(05,05) - rrt(029) * density(14) 
  pd(05,14) = pd(05,14) - rrt(029) * density(05) 
  pd(06,05) = pd(06,05) + rrt(029) * density(14) 
  pd(06,14) = pd(06,14) + rrt(029) * density(05) 
  pd(06,06) = pd(06,06) - rrt(030) * density(14) 
  pd(06,14) = pd(06,14) - rrt(030) * density(06) 
  pd(07,06) = pd(07,06) + rrt(030) * density(14) 
  pd(07,14) = pd(07,14) + rrt(030) * density(06) 
  pd(07,07) = pd(07,07) - rrt(031) * density(14) 
  pd(07,14) = pd(07,14) - rrt(031) * density(07) 
  pd(08,07) = pd(08,07) + rrt(031) * density(14) 
  pd(08,14) = pd(08,14) + rrt(031) * density(07) 
  pd(08,08) = pd(08,08) - rrt(032) * density(14) 
  pd(08,14) = pd(08,14) - rrt(032) * density(08) 
  pd(09,08) = pd(09,08) + rrt(032) * density(14) 
  pd(09,14) = pd(09,14) + rrt(032) * density(08) 
  pd(01,02) = pd(01,02) + rrt(033) * density(29) 
  pd(01,29) = pd(01,29) + rrt(033) * density(02) 
  pd(02,02) = pd(02,02) - rrt(033) * density(29) 
  pd(02,29) = pd(02,29) - rrt(033) * density(02) 
  pd(02,03) = pd(02,03) + rrt(034) * density(29) 
  pd(02,29) = pd(02,29) + rrt(034) * density(03) 
  pd(03,03) = pd(03,03) - rrt(034) * density(29) 
  pd(03,29) = pd(03,29) - rrt(034) * density(03) 
  pd(03,04) = pd(03,04) + rrt(035) * density(29) 
  pd(03,29) = pd(03,29) + rrt(035) * density(04) 
  pd(04,04) = pd(04,04) - rrt(035) * density(29) 
  pd(04,29) = pd(04,29) - rrt(035) * density(04) 
  pd(04,05) = pd(04,05) + rrt(036) * density(29) 
  pd(04,29) = pd(04,29) + rrt(036) * density(05) 
  pd(05,05) = pd(05,05) - rrt(036) * density(29) 
  pd(05,29) = pd(05,29) - rrt(036) * density(05) 
  pd(05,06) = pd(05,06) + rrt(037) * density(29) 
  pd(05,29) = pd(05,29) + rrt(037) * density(06) 
  pd(06,06) = pd(06,06) - rrt(037) * density(29) 
  pd(06,29) = pd(06,29) - rrt(037) * density(06) 
  pd(06,07) = pd(06,07) + rrt(038) * density(29) 
  pd(06,29) = pd(06,29) + rrt(038) * density(07) 
  pd(07,07) = pd(07,07) - rrt(038) * density(29) 
  pd(07,29) = pd(07,29) - rrt(038) * density(07) 
  pd(07,08) = pd(07,08) + rrt(039) * density(29) 
  pd(07,29) = pd(07,29) + rrt(039) * density(08) 
  pd(08,08) = pd(08,08) - rrt(039) * density(29) 
  pd(08,29) = pd(08,29) - rrt(039) * density(08) 
  pd(08,09) = pd(08,09) + rrt(040) * density(29) 
  pd(08,29) = pd(08,29) + rrt(040) * density(09) 
  pd(09,09) = pd(09,09) - rrt(040) * density(29) 
  pd(09,29) = pd(09,29) - rrt(040) * density(09) 
  pd(01,01) = pd(01,01) - rrt(041) * density(29) 
  pd(01,29) = pd(01,29) - rrt(041) * density(01) 
  pd(02,01) = pd(02,01) + rrt(041) * density(29) 
  pd(02,29) = pd(02,29) + rrt(041) * density(01) 
  pd(02,02) = pd(02,02) - rrt(042) * density(29) 
  pd(02,29) = pd(02,29) - rrt(042) * density(02) 
  pd(03,02) = pd(03,02) + rrt(042) * density(29) 
  pd(03,29) = pd(03,29) + rrt(042) * density(02) 
  pd(03,03) = pd(03,03) - rrt(043) * density(29) 
  pd(03,29) = pd(03,29) - rrt(043) * density(03) 
  pd(04,03) = pd(04,03) + rrt(043) * density(29) 
  pd(04,29) = pd(04,29) + rrt(043) * density(03) 
  pd(04,04) = pd(04,04) - rrt(044) * density(29) 
  pd(04,29) = pd(04,29) - rrt(044) * density(04) 
  pd(05,04) = pd(05,04) + rrt(044) * density(29) 
  pd(05,29) = pd(05,29) + rrt(044) * density(04) 
  pd(05,05) = pd(05,05) - rrt(045) * density(29) 
  pd(05,29) = pd(05,29) - rrt(045) * density(05) 
  pd(06,05) = pd(06,05) + rrt(045) * density(29) 
  pd(06,29) = pd(06,29) + rrt(045) * density(05) 
  pd(06,06) = pd(06,06) - rrt(046) * density(29) 
  pd(06,29) = pd(06,29) - rrt(046) * density(06) 
  pd(07,06) = pd(07,06) + rrt(046) * density(29) 
  pd(07,29) = pd(07,29) + rrt(046) * density(06) 
  pd(07,07) = pd(07,07) - rrt(047) * density(29) 
  pd(07,29) = pd(07,29) - rrt(047) * density(07) 
  pd(08,07) = pd(08,07) + rrt(047) * density(29) 
  pd(08,29) = pd(08,29) + rrt(047) * density(07) 
  pd(08,08) = pd(08,08) - rrt(048) * density(29) 
  pd(08,29) = pd(08,29) - rrt(048) * density(08) 
  pd(09,08) = pd(09,08) + rrt(048) * density(29) 
  pd(09,29) = pd(09,29) + rrt(048) * density(08) 
  pd(21,21) = pd(21,21) + rrt(049) * density(22) 
  pd(21,22) = pd(21,22) + rrt(049) * density(21) 
  pd(22,21) = pd(22,21) - rrt(049) * density(22) 
  pd(22,22) = pd(22,22) - rrt(049) * density(21) 
  pd(22,21) = pd(22,21) + rrt(050) * density(23) 
  pd(22,23) = pd(22,23) + rrt(050) * density(21) 
  pd(23,21) = pd(23,21) - rrt(050) * density(23) 
  pd(23,23) = pd(23,23) - rrt(050) * density(21) 
  pd(23,21) = pd(23,21) + rrt(051) * density(24) 
  pd(23,24) = pd(23,24) + rrt(051) * density(21) 
  pd(24,21) = pd(24,21) - rrt(051) * density(24) 
  pd(24,24) = pd(24,24) - rrt(051) * density(21) 
  pd(24,21) = pd(24,21) + rrt(052) * density(25) 
  pd(24,25) = pd(24,25) + rrt(052) * density(21) 
  pd(25,21) = pd(25,21) - rrt(052) * density(25) 
  pd(25,25) = pd(25,25) - rrt(052) * density(21) 
  pd(21,21) = pd(21,21) - rrt(053) * density(21) * 2.0d0
  pd(22,21) = pd(22,21) + rrt(053) * density(21) * 2.0d0
  pd(22,21) = pd(22,21) - rrt(054) * density(22) 
  pd(22,22) = pd(22,22) - rrt(054) * density(21) 
  pd(23,21) = pd(23,21) + rrt(054) * density(22) 
  pd(23,22) = pd(23,22) + rrt(054) * density(21) 
  pd(23,21) = pd(23,21) - rrt(055) * density(23) 
  pd(23,23) = pd(23,23) - rrt(055) * density(21) 
  pd(24,21) = pd(24,21) + rrt(055) * density(23) 
  pd(24,23) = pd(24,23) + rrt(055) * density(21) 
  pd(24,21) = pd(24,21) - rrt(056) * density(24) 
  pd(24,24) = pd(24,24) - rrt(056) * density(21) 
  pd(25,21) = pd(25,21) + rrt(056) * density(24) 
  pd(25,24) = pd(25,24) + rrt(056) * density(21) 
  pd(21,22) = pd(21,22) + rrt(057) * density(29) 
  pd(21,29) = pd(21,29) + rrt(057) * density(22) 
  pd(22,22) = pd(22,22) - rrt(057) * density(29) 
  pd(22,29) = pd(22,29) - rrt(057) * density(22) 
  pd(22,23) = pd(22,23) + rrt(058) * density(29) 
  pd(22,29) = pd(22,29) + rrt(058) * density(23) 
  pd(23,23) = pd(23,23) - rrt(058) * density(29) 
  pd(23,29) = pd(23,29) - rrt(058) * density(23) 
  pd(23,24) = pd(23,24) + rrt(059) * density(29) 
  pd(23,29) = pd(23,29) + rrt(059) * density(24) 
  pd(24,24) = pd(24,24) - rrt(059) * density(29) 
  pd(24,29) = pd(24,29) - rrt(059) * density(24) 
  pd(24,25) = pd(24,25) + rrt(060) * density(29) 
  pd(24,29) = pd(24,29) + rrt(060) * density(25) 
  pd(25,25) = pd(25,25) - rrt(060) * density(29) 
  pd(25,29) = pd(25,29) - rrt(060) * density(25) 
  pd(21,21) = pd(21,21) - rrt(061) * density(29) 
  pd(21,29) = pd(21,29) - rrt(061) * density(21) 
  pd(22,21) = pd(22,21) + rrt(061) * density(29) 
  pd(22,29) = pd(22,29) + rrt(061) * density(21) 
  pd(22,22) = pd(22,22) - rrt(062) * density(29) 
  pd(22,29) = pd(22,29) - rrt(062) * density(22) 
  pd(23,22) = pd(23,22) + rrt(062) * density(29) 
  pd(23,29) = pd(23,29) + rrt(062) * density(22) 
  pd(23,23) = pd(23,23) - rrt(063) * density(29) 
  pd(23,29) = pd(23,29) - rrt(063) * density(23) 
  pd(24,23) = pd(24,23) + rrt(063) * density(29) 
  pd(24,29) = pd(24,29) + rrt(063) * density(23) 
  pd(24,24) = pd(24,24) - rrt(064) * density(29) 
  pd(24,29) = pd(24,29) - rrt(064) * density(24) 
  pd(25,24) = pd(25,24) + rrt(064) * density(29) 
  pd(25,29) = pd(25,29) + rrt(064) * density(24) 
  pd(01,01) = pd(01,01) - rrt(065) * density(53) 
  pd(01,53) = pd(01,53) - rrt(065) * density(01) 
  pd(10,01) = pd(10,01) + rrt(065) * density(53) 
  pd(10,53) = pd(10,53) + rrt(065) * density(01) 
  pd(01,01) = pd(01,01) - rrt(066) * density(53) 
  pd(01,53) = pd(01,53) - rrt(066) * density(01) 
  pd(12,01) = pd(12,01) + rrt(066) * density(53) 
  pd(12,53) = pd(12,53) + rrt(066) * density(01) 
  pd(01,01) = pd(01,01) - rrt(067) * density(53) 
  pd(01,53) = pd(01,53) - rrt(067) * density(01) 
  pd(12,01) = pd(12,01) + rrt(067) * density(53) 
  pd(12,53) = pd(12,53) + rrt(067) * density(01) 
  pd(01,01) = pd(01,01) - rrt(068) * density(53) 
  pd(01,53) = pd(01,53) - rrt(068) * density(01) 
  pd(12,01) = pd(12,01) + rrt(068) * density(53) 
  pd(12,53) = pd(12,53) + rrt(068) * density(01) 
  pd(01,01) = pd(01,01) - rrt(069) * density(53) 
  pd(01,53) = pd(01,53) - rrt(069) * density(01) 
  pd(13,01) = pd(13,01) + rrt(069) * density(53) 
  pd(13,53) = pd(13,53) + rrt(069) * density(01) 
  pd(01,01) = pd(01,01) - rrt(070) * density(53) 
  pd(01,53) = pd(01,53) - rrt(070) * density(01) 
  pd(13,01) = pd(13,01) + rrt(070) * density(53) 
  pd(13,53) = pd(13,53) + rrt(070) * density(01) 
  pd(01,01) = pd(01,01) - rrt(071) * density(53) 
  pd(01,53) = pd(01,53) - rrt(071) * density(01) 
  pd(13,01) = pd(13,01) + rrt(071) * density(53) 
  pd(13,53) = pd(13,53) + rrt(071) * density(01) 
  pd(01,01) = pd(01,01) - rrt(072) * density(53) 
  pd(01,53) = pd(01,53) - rrt(072) * density(01) 
  pd(14,01) = pd(14,01) + rrt(072) * density(53) 
  pd(14,53) = pd(14,53) + rrt(072) * density(01) 
  pd(15,01) = pd(15,01) + rrt(072) * density(53) 
  pd(15,53) = pd(15,53) + rrt(072) * density(01) 
  pd(21,21) = pd(21,21) - rrt(073) * density(53) 
  pd(21,53) = pd(21,53) - rrt(073) * density(21) 
  pd(26,21) = pd(26,21) + rrt(073) * density(53) 
  pd(26,53) = pd(26,53) + rrt(073) * density(21) 
  pd(21,21) = pd(21,21) - rrt(074) * density(53) 
  pd(21,53) = pd(21,53) - rrt(074) * density(21) 
  pd(27,21) = pd(27,21) + rrt(074) * density(53) 
  pd(27,53) = pd(27,53) + rrt(074) * density(21) 
  pd(21,21) = pd(21,21) - rrt(075) * density(53) 
  pd(21,53) = pd(21,53) - rrt(075) * density(21) 
  pd(28,21) = pd(28,21) + rrt(075) * density(53) 
  pd(28,53) = pd(28,53) + rrt(075) * density(21) 
  pd(21,21) = pd(21,21) - rrt(076) * density(53) 
  pd(21,53) = pd(21,53) - rrt(076) * density(21) 
  pd(29,21) = pd(29,21) + rrt(076) * density(53) * 2.0d0
  pd(29,53) = pd(29,53) + rrt(076) * density(21) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(077) * density(53) 
  pd(21,53) = pd(21,53) - rrt(077) * density(21) 
  pd(29,21) = pd(29,21) + rrt(077) * density(53) 
  pd(29,53) = pd(29,53) + rrt(077) * density(21) 
  pd(30,21) = pd(30,21) + rrt(077) * density(53) 
  pd(30,53) = pd(30,53) + rrt(077) * density(21) 
  pd(21,21) = pd(21,21) - rrt(078) * density(53) 
  pd(21,53) = pd(21,53) - rrt(078) * density(21) 
  pd(29,21) = pd(29,21) + rrt(078) * density(53) 
  pd(29,53) = pd(29,53) + rrt(078) * density(21) 
  pd(31,21) = pd(31,21) + rrt(078) * density(53) 
  pd(31,53) = pd(31,53) + rrt(078) * density(21) 
  pd(26,26) = pd(26,26) - rrt(079) * density(53) 
  pd(26,53) = pd(26,53) - rrt(079) * density(26) 
  pd(29,26) = pd(29,26) + rrt(079) * density(53) * 2.0d0
  pd(29,53) = pd(29,53) + rrt(079) * density(26) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(080) * density(53) 
  pd(29,53) = pd(29,53) - rrt(080) * density(29) 
  pd(30,29) = pd(30,29) + rrt(080) * density(53) 
  pd(30,53) = pd(30,53) + rrt(080) * density(29) 
  pd(29,29) = pd(29,29) - rrt(081) * density(53) 
  pd(29,53) = pd(29,53) - rrt(081) * density(29) 
  pd(31,29) = pd(31,29) + rrt(081) * density(53) 
  pd(31,53) = pd(31,53) + rrt(081) * density(29) 
  pd(01,10) = pd(01,10) + rrt(082) * density(53) 
  pd(01,53) = pd(01,53) + rrt(082) * density(10) 
  pd(10,10) = pd(10,10) - rrt(082) * density(53) 
  pd(10,53) = pd(10,53) - rrt(082) * density(10) 
  pd(21,26) = pd(21,26) + rrt(083) * density(53) 
  pd(21,53) = pd(21,53) + rrt(083) * density(26) 
  pd(26,26) = pd(26,26) - rrt(083) * density(53) 
  pd(26,53) = pd(26,53) - rrt(083) * density(26) 
  pd(14,14) = pd(14,14) - rrt(084) * density(53) 
  pd(14,53) = pd(14,53) - rrt(084) * density(14) 
  pd(17,14) = pd(17,14) + rrt(084) * density(53) 
  pd(17,53) = pd(17,53) + rrt(084) * density(14) 
  pd(53,14) = pd(53,14) + rrt(084) * density(53) 
  pd(53,53) = pd(53,53) + rrt(084) * density(14) 
  pd(29,29) = pd(29,29) - rrt(085) * density(53) 
  pd(29,53) = pd(29,53) - rrt(085) * density(29) 
  pd(33,29) = pd(33,29) + rrt(085) * density(53) 
  pd(33,53) = pd(33,53) + rrt(085) * density(29) 
  pd(53,29) = pd(53,29) + rrt(085) * density(53) 
  pd(53,53) = pd(53,53) + rrt(085) * density(29) 
  pd(01,01) = pd(01,01) - rrt(086) * density(53) 
  pd(01,53) = pd(01,53) - rrt(086) * density(01) 
  pd(18,01) = pd(18,01) + rrt(086) * density(53) 
  pd(18,53) = pd(18,53) + rrt(086) * density(01) 
  pd(53,01) = pd(53,01) + rrt(086) * density(53) 
  pd(53,53) = pd(53,53) + rrt(086) * density(01) 
  pd(10,10) = pd(10,10) - rrt(087) * density(53) 
  pd(10,53) = pd(10,53) - rrt(087) * density(10) 
  pd(18,10) = pd(18,10) + rrt(087) * density(53) 
  pd(18,53) = pd(18,53) + rrt(087) * density(10) 
  pd(53,10) = pd(53,10) + rrt(087) * density(53) 
  pd(53,53) = pd(53,53) + rrt(087) * density(10) 
  pd(21,21) = pd(21,21) - rrt(088) * density(53) 
  pd(21,53) = pd(21,53) - rrt(088) * density(21) 
  pd(34,21) = pd(34,21) + rrt(088) * density(53) 
  pd(34,53) = pd(34,53) + rrt(088) * density(21) 
  pd(53,21) = pd(53,21) + rrt(088) * density(53) 
  pd(53,53) = pd(53,53) + rrt(088) * density(21) 
  pd(26,26) = pd(26,26) - rrt(089) * density(53) 
  pd(26,53) = pd(26,53) - rrt(089) * density(26) 
  pd(34,26) = pd(34,26) + rrt(089) * density(53) 
  pd(34,53) = pd(34,53) + rrt(089) * density(26) 
  pd(53,26) = pd(53,26) + rrt(089) * density(53) 
  pd(53,53) = pd(53,53) + rrt(089) * density(26) 
  pd(40,40) = pd(40,40) - rrt(090) * density(53) 
  pd(40,53) = pd(40,53) - rrt(090) * density(40) 
  pd(45,40) = pd(45,40) + rrt(090) * density(53) 
  pd(45,53) = pd(45,53) + rrt(090) * density(40) 
  pd(53,40) = pd(53,40) + rrt(090) * density(53) 
  pd(53,53) = pd(53,53) + rrt(090) * density(40) 
  pd(41,41) = pd(41,41) - rrt(091) * density(53) 
  pd(41,53) = pd(41,53) - rrt(091) * density(41) 
  pd(46,41) = pd(46,41) + rrt(091) * density(53) 
  pd(46,53) = pd(46,53) + rrt(091) * density(41) 
  pd(53,41) = pd(53,41) + rrt(091) * density(53) 
  pd(53,53) = pd(53,53) + rrt(091) * density(41) 
  pd(14,18) = pd(14,18) + rrt(092) * density(53) * 2.0d0
  pd(14,53) = pd(14,53) + rrt(092) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(092) * density(53) 
  pd(18,53) = pd(18,53) - rrt(092) * density(18) 
  pd(53,18) = pd(53,18) - rrt(092) * density(53) 
  pd(53,53) = pd(53,53) - rrt(092) * density(18) 
  pd(14,18) = pd(14,18) + rrt(093) * density(53) 
  pd(14,53) = pd(14,53) + rrt(093) * density(18) 
  pd(15,18) = pd(15,18) + rrt(093) * density(53) 
  pd(15,53) = pd(15,53) + rrt(093) * density(18) 
  pd(18,18) = pd(18,18) - rrt(093) * density(53) 
  pd(18,53) = pd(18,53) - rrt(093) * density(18) 
  pd(53,18) = pd(53,18) - rrt(093) * density(53) 
  pd(53,53) = pd(53,53) - rrt(093) * density(18) 
  pd(14,18) = pd(14,18) + rrt(094) * density(53) 
  pd(14,53) = pd(14,53) + rrt(094) * density(18) 
  pd(16,18) = pd(16,18) + rrt(094) * density(53) 
  pd(16,53) = pd(16,53) + rrt(094) * density(18) 
  pd(18,18) = pd(18,18) - rrt(094) * density(53) 
  pd(18,53) = pd(18,53) - rrt(094) * density(18) 
  pd(53,18) = pd(53,18) - rrt(094) * density(53) 
  pd(53,53) = pd(53,53) - rrt(094) * density(18) 
  pd(29,34) = pd(29,34) + rrt(095) * density(53) * 2.0d0
  pd(29,53) = pd(29,53) + rrt(095) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(095) * density(53) 
  pd(34,53) = pd(34,53) - rrt(095) * density(34) 
  pd(53,34) = pd(53,34) - rrt(095) * density(53) 
  pd(53,53) = pd(53,53) - rrt(095) * density(34) 
  pd(29,34) = pd(29,34) + rrt(096) * density(53) 
  pd(29,53) = pd(29,53) + rrt(096) * density(34) 
  pd(30,34) = pd(30,34) + rrt(096) * density(53) 
  pd(30,53) = pd(30,53) + rrt(096) * density(34) 
  pd(34,34) = pd(34,34) - rrt(096) * density(53) 
  pd(34,53) = pd(34,53) - rrt(096) * density(34) 
  pd(53,34) = pd(53,34) - rrt(096) * density(53) 
  pd(53,53) = pd(53,53) - rrt(096) * density(34) 
  pd(29,34) = pd(29,34) + rrt(097) * density(53) 
  pd(29,53) = pd(29,53) + rrt(097) * density(34) 
  pd(31,34) = pd(31,34) + rrt(097) * density(53) 
  pd(31,53) = pd(31,53) + rrt(097) * density(34) 
  pd(34,34) = pd(34,34) - rrt(097) * density(53) 
  pd(34,53) = pd(34,53) - rrt(097) * density(34) 
  pd(53,34) = pd(53,34) - rrt(097) * density(53) 
  pd(53,53) = pd(53,53) - rrt(097) * density(34) 
  pd(14,45) = pd(14,45) + rrt(098) * density(53) 
  pd(14,53) = pd(14,53) + rrt(098) * density(45) 
  pd(29,45) = pd(29,45) + rrt(098) * density(53) 
  pd(29,53) = pd(29,53) + rrt(098) * density(45) 
  pd(45,45) = pd(45,45) - rrt(098) * density(53) 
  pd(45,53) = pd(45,53) - rrt(098) * density(45) 
  pd(53,45) = pd(53,45) - rrt(098) * density(53) 
  pd(53,53) = pd(53,53) - rrt(098) * density(45) 
  pd(15,45) = pd(15,45) + rrt(099) * density(53) 
  pd(15,53) = pd(15,53) + rrt(099) * density(45) 
  pd(29,45) = pd(29,45) + rrt(099) * density(53) 
  pd(29,53) = pd(29,53) + rrt(099) * density(45) 
  pd(45,45) = pd(45,45) - rrt(099) * density(53) 
  pd(45,53) = pd(45,53) - rrt(099) * density(45) 
  pd(53,45) = pd(53,45) - rrt(099) * density(53) 
  pd(53,53) = pd(53,53) - rrt(099) * density(45) 
  pd(01,19) = pd(01,19) + rrt(100) * density(53) 
  pd(01,53) = pd(01,53) + rrt(100) * density(19) 
  pd(14,19) = pd(14,19) + rrt(100) * density(53) 
  pd(14,53) = pd(14,53) + rrt(100) * density(19) 
  pd(19,19) = pd(19,19) - rrt(100) * density(53) 
  pd(19,53) = pd(19,53) - rrt(100) * density(19) 
  pd(53,19) = pd(53,19) - rrt(100) * density(53) 
  pd(53,53) = pd(53,53) - rrt(100) * density(19) 
  pd(01,20) = pd(01,20) + rrt(101) * density(53) * 2.0d0
  pd(01,53) = pd(01,53) + rrt(101) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(101) * density(53) 
  pd(20,53) = pd(20,53) - rrt(101) * density(20) 
  pd(53,20) = pd(53,20) - rrt(101) * density(53) 
  pd(53,53) = pd(53,53) - rrt(101) * density(20) 
  pd(01,46) = pd(01,46) + rrt(102) * density(53) 
  pd(01,53) = pd(01,53) + rrt(102) * density(46) 
  pd(29,46) = pd(29,46) + rrt(102) * density(53) 
  pd(29,53) = pd(29,53) + rrt(102) * density(46) 
  pd(46,46) = pd(46,46) - rrt(102) * density(53) 
  pd(46,53) = pd(46,53) - rrt(102) * density(46) 
  pd(53,46) = pd(53,46) - rrt(102) * density(53) 
  pd(53,53) = pd(53,53) - rrt(102) * density(46) 
  pd(29,47) = pd(29,47) + rrt(103) * density(53) 
  pd(29,53) = pd(29,53) + rrt(103) * density(47) 
  pd(40,47) = pd(40,47) + rrt(103) * density(53) 
  pd(40,53) = pd(40,53) + rrt(103) * density(47) 
  pd(47,47) = pd(47,47) - rrt(103) * density(53) 
  pd(47,53) = pd(47,53) - rrt(103) * density(47) 
  pd(53,47) = pd(53,47) - rrt(103) * density(53) 
  pd(53,53) = pd(53,53) - rrt(103) * density(47) 
  pd(21,35) = pd(21,35) + rrt(104) * density(53) * 2.0d0
  pd(21,53) = pd(21,53) + rrt(104) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(104) * density(53) 
  pd(35,53) = pd(35,53) - rrt(104) * density(35) 
  pd(53,35) = pd(53,35) - rrt(104) * density(53) 
  pd(53,53) = pd(53,53) - rrt(104) * density(35) 
  pd(01,52) = pd(01,52) + rrt(105) * density(53) 
  pd(01,53) = pd(01,53) + rrt(105) * density(52) 
  pd(21,52) = pd(21,52) + rrt(105) * density(53) 
  pd(21,53) = pd(21,53) + rrt(105) * density(52) 
  pd(52,52) = pd(52,52) - rrt(105) * density(53) 
  pd(52,53) = pd(52,53) - rrt(105) * density(52) 
  pd(53,52) = pd(53,52) - rrt(105) * density(53) 
  pd(53,53) = pd(53,53) - rrt(105) * density(52) 
  pd(14,17) = pd(14,17) + rrt(106) * density(53)**2 
  pd(14,53) = pd(14,53) + rrt(106) * density(17) * density(53) * 2.0d0
  pd(17,17) = pd(17,17) - rrt(106) * density(53)**2 
  pd(17,53) = pd(17,53) - rrt(106) * density(17) * density(53) * 2.0d0
  pd(53,17) = pd(53,17) - rrt(106) * density(53)**2 
  pd(53,53) = pd(53,53) - rrt(106) * density(17) * density(53) * 2.0d0
  pd(29,33) = pd(29,33) + rrt(107) * density(53)**2 
  pd(29,53) = pd(29,53) + rrt(107) * density(33) * density(53) * 2.0d0
  pd(33,33) = pd(33,33) - rrt(107) * density(53)**2 
  pd(33,53) = pd(33,53) - rrt(107) * density(33) * density(53) * 2.0d0
  pd(53,33) = pd(53,33) - rrt(107) * density(53)**2 
  pd(53,53) = pd(53,53) - rrt(107) * density(33) * density(53) * 2.0d0
  pd(14,17) = pd(14,17) + rrt(108) * density(53) 
  pd(14,53) = pd(14,53) + rrt(108) * density(17) 
  pd(17,17) = pd(17,17) - rrt(108) * density(53) 
  pd(17,53) = pd(17,53) - rrt(108) * density(17) 
  pd(53,17) = pd(53,17) - rrt(108) * density(53) 
  pd(53,53) = pd(53,53) - rrt(108) * density(17) 
  pd(29,33) = pd(29,33) + rrt(109) * density(53) 
  pd(29,53) = pd(29,53) + rrt(109) * density(33) 
  pd(33,33) = pd(33,33) - rrt(109) * density(53) 
  pd(33,53) = pd(33,53) - rrt(109) * density(33) 
  pd(53,33) = pd(53,33) - rrt(109) * density(53) 
  pd(53,53) = pd(53,53) - rrt(109) * density(33) 
  pd(21,21) = pd(21,21) - rrt(110) * density(53) 
  pd(21,53) = pd(21,53) - rrt(110) * density(21) 
  pd(29,21) = pd(29,21) + rrt(110) * density(53) 
  pd(29,53) = pd(29,53) + rrt(110) * density(21) 
  pd(36,21) = pd(36,21) + rrt(110) * density(53) 
  pd(36,53) = pd(36,53) + rrt(110) * density(21) 
  pd(53,21) = pd(53,21) - rrt(110) * density(53) 
  pd(53,53) = pd(53,53) - rrt(110) * density(21) 
  pd(14,40) = pd(14,40) + rrt(111) * density(53) 
  pd(14,53) = pd(14,53) + rrt(111) * density(40) 
  pd(36,40) = pd(36,40) + rrt(111) * density(53) 
  pd(36,53) = pd(36,53) + rrt(111) * density(40) 
  pd(40,40) = pd(40,40) - rrt(111) * density(53) 
  pd(40,53) = pd(40,53) - rrt(111) * density(40) 
  pd(53,40) = pd(53,40) - rrt(111) * density(53) 
  pd(53,53) = pd(53,53) - rrt(111) * density(40) 
  pd(21,32) = pd(21,32) + rrt(112) * density(53) 
  pd(21,53) = pd(21,53) + rrt(112) * density(32) 
  pd(32,32) = pd(32,32) - rrt(112) * density(53) 
  pd(32,53) = pd(32,53) - rrt(112) * density(32) 
  pd(36,32) = pd(36,32) + rrt(112) * density(53) 
  pd(36,53) = pd(36,53) + rrt(112) * density(32) 
  pd(53,32) = pd(53,32) - rrt(112) * density(53) 
  pd(53,53) = pd(53,53) - rrt(112) * density(32) 
  pd(29,32) = pd(29,32) + rrt(113) * density(53) 
  pd(29,53) = pd(29,53) + rrt(113) * density(32) 
  pd(32,32) = pd(32,32) - rrt(113) * density(53) 
  pd(32,53) = pd(32,53) - rrt(113) * density(32) 
  pd(37,32) = pd(37,32) + rrt(113) * density(53) 
  pd(37,53) = pd(37,53) + rrt(113) * density(32) 
  pd(53,32) = pd(53,32) - rrt(113) * density(53) 
  pd(53,53) = pd(53,53) - rrt(113) * density(32) 
  pd(14,41) = pd(14,41) + rrt(114) * density(53) 
  pd(14,53) = pd(14,53) + rrt(114) * density(41) 
  pd(41,41) = pd(41,41) - rrt(114) * density(53) 
  pd(41,53) = pd(41,53) - rrt(114) * density(41) 
  pd(48,41) = pd(48,41) + rrt(114) * density(53) 
  pd(48,53) = pd(48,53) + rrt(114) * density(41) 
  pd(53,41) = pd(53,41) - rrt(114) * density(53) 
  pd(53,53) = pd(53,53) - rrt(114) * density(41) 
  pd(21,21) = pd(21,21) - rrt(115) * density(21) * density(53) * 2.0d0
  pd(21,53) = pd(21,53) - rrt(115) * density(21)**2 
  pd(37,21) = pd(37,21) + rrt(115) * density(21) * density(53) * 2.0d0
  pd(37,53) = pd(37,53) + rrt(115) * density(21)**2 
  pd(53,21) = pd(53,21) - rrt(115) * density(21) * density(53) * 2.0d0
  pd(53,53) = pd(53,53) - rrt(115) * density(21)**2 
  pd(36,42) = pd(36,42) + rrt(116) * density(53) 
  pd(36,53) = pd(36,53) + rrt(116) * density(42) 
  pd(40,42) = pd(40,42) + rrt(116) * density(53) 
  pd(40,53) = pd(40,53) + rrt(116) * density(42) 
  pd(42,42) = pd(42,42) - rrt(116) * density(53) 
  pd(42,53) = pd(42,53) - rrt(116) * density(42) 
  pd(53,42) = pd(53,42) - rrt(116) * density(53) 
  pd(53,53) = pd(53,53) - rrt(116) * density(42) 
  pd(29,21) = pd(29,21) - rrt(117) * density(29) * density(53) 
  pd(29,29) = pd(29,29) - rrt(117) * density(21) * density(53) 
  pd(29,53) = pd(29,53) - rrt(117) * density(21) * density(29) 
  pd(36,21) = pd(36,21) + rrt(117) * density(29) * density(53) 
  pd(36,29) = pd(36,29) + rrt(117) * density(21) * density(53) 
  pd(36,53) = pd(36,53) + rrt(117) * density(21) * density(29) 
  pd(53,21) = pd(53,21) - rrt(117) * density(29) * density(53) 
  pd(53,29) = pd(53,29) - rrt(117) * density(21) * density(53) 
  pd(53,53) = pd(53,53) - rrt(117) * density(21) * density(29) 
  pd(21,21) = pd(21,21) - rrt(118) * density(29) * density(53) 
  pd(21,29) = pd(21,29) - rrt(118) * density(21) * density(53) 
  pd(21,53) = pd(21,53) - rrt(118) * density(21) * density(29) 
  pd(37,21) = pd(37,21) + rrt(118) * density(29) * density(53) 
  pd(37,29) = pd(37,29) + rrt(118) * density(21) * density(53) 
  pd(37,53) = pd(37,53) + rrt(118) * density(21) * density(29) 
  pd(53,21) = pd(53,21) - rrt(118) * density(29) * density(53) 
  pd(53,29) = pd(53,29) - rrt(118) * density(21) * density(53) 
  pd(53,53) = pd(53,53) - rrt(118) * density(21) * density(29) 
  pd(32,32) = pd(32,32) - rrt(119) * density(53) 
  pd(32,53) = pd(32,53) - rrt(119) * density(32) 
  pd(38,32) = pd(38,32) + rrt(119) * density(53) 
  pd(38,53) = pd(38,53) + rrt(119) * density(32) 
  pd(53,32) = pd(53,32) - rrt(119) * density(53) 
  pd(53,53) = pd(53,53) - rrt(119) * density(32) 
  pd(40,40) = pd(40,40) - rrt(120) * density(53) 
  pd(40,53) = pd(40,53) - rrt(120) * density(40) 
  pd(48,40) = pd(48,40) + rrt(120) * density(53) 
  pd(48,53) = pd(48,53) + rrt(120) * density(40) 
  pd(53,40) = pd(53,40) - rrt(120) * density(53) 
  pd(53,53) = pd(53,53) - rrt(120) * density(40) 
  pd(41,41) = pd(41,41) - rrt(121) * density(53) 
  pd(41,53) = pd(41,53) - rrt(121) * density(41) 
  pd(49,41) = pd(49,41) + rrt(121) * density(53) 
  pd(49,53) = pd(49,53) + rrt(121) * density(41) 
  pd(53,41) = pd(53,41) - rrt(121) * density(53) 
  pd(53,53) = pd(53,53) - rrt(121) * density(41) 
  pd(21,01) = pd(21,01) - rrt(122) * density(21) * density(53) 
  pd(21,21) = pd(21,21) - rrt(122) * density(01) * density(53) 
  pd(21,53) = pd(21,53) - rrt(122) * density(01) * density(21) 
  pd(37,01) = pd(37,01) + rrt(122) * density(21) * density(53) 
  pd(37,21) = pd(37,21) + rrt(122) * density(01) * density(53) 
  pd(37,53) = pd(37,53) + rrt(122) * density(01) * density(21) 
  pd(53,01) = pd(53,01) - rrt(122) * density(21) * density(53) 
  pd(53,21) = pd(53,21) - rrt(122) * density(01) * density(53) 
  pd(53,53) = pd(53,53) - rrt(122) * density(01) * density(21) 
  pd(21,29) = pd(21,29) + rrt(123) * density(36) 
  pd(21,36) = pd(21,36) + rrt(123) * density(29) 
  pd(29,29) = pd(29,29) - rrt(123) * density(36) 
  pd(29,36) = pd(29,36) - rrt(123) * density(29) 
  pd(36,29) = pd(36,29) - rrt(123) * density(36) 
  pd(36,36) = pd(36,36) - rrt(123) * density(29) 
  pd(53,29) = pd(53,29) + rrt(123) * density(36) 
  pd(53,36) = pd(53,36) + rrt(123) * density(29) 
  pd(14,14) = pd(14,14) - rrt(124) * density(36) 
  pd(14,36) = pd(14,36) - rrt(124) * density(14) 
  pd(36,14) = pd(36,14) - rrt(124) * density(36) 
  pd(36,36) = pd(36,36) - rrt(124) * density(14) 
  pd(40,14) = pd(40,14) + rrt(124) * density(36) 
  pd(40,36) = pd(40,36) + rrt(124) * density(14) 
  pd(53,14) = pd(53,14) + rrt(124) * density(36) 
  pd(53,36) = pd(53,36) + rrt(124) * density(14) 
  pd(36,36) = pd(36,36) - rrt(125) * density(40) 
  pd(36,40) = pd(36,40) - rrt(125) * density(36) 
  pd(40,36) = pd(40,36) - rrt(125) * density(40) 
  pd(40,40) = pd(40,40) - rrt(125) * density(36) 
  pd(42,36) = pd(42,36) + rrt(125) * density(40) 
  pd(42,40) = pd(42,40) + rrt(125) * density(36) 
  pd(53,36) = pd(53,36) + rrt(125) * density(40) 
  pd(53,40) = pd(53,40) + rrt(125) * density(36) 
  pd(01,01) = pd(01,01) - rrt(126) * density(36) 
  pd(01,36) = pd(01,36) - rrt(126) * density(01) 
  pd(36,01) = pd(36,01) - rrt(126) * density(36) 
  pd(36,36) = pd(36,36) - rrt(126) * density(01) 
  pd(41,01) = pd(41,01) + rrt(126) * density(36) 
  pd(41,36) = pd(41,36) + rrt(126) * density(01) 
  pd(53,01) = pd(53,01) + rrt(126) * density(36) 
  pd(53,36) = pd(53,36) + rrt(126) * density(01) 
  pd(21,21) = pd(21,21) - rrt(127) * density(36) 
  pd(21,36) = pd(21,36) - rrt(127) * density(21) 
  pd(32,21) = pd(32,21) + rrt(127) * density(36) 
  pd(32,36) = pd(32,36) + rrt(127) * density(21) 
  pd(36,21) = pd(36,21) - rrt(127) * density(36) 
  pd(36,36) = pd(36,36) - rrt(127) * density(21) 
  pd(53,21) = pd(53,21) + rrt(127) * density(36) 
  pd(53,36) = pd(53,36) + rrt(127) * density(21) 
  pd(26,26) = pd(26,26) - rrt(128) * density(36) 
  pd(26,36) = pd(26,36) - rrt(128) * density(26) 
  pd(32,26) = pd(32,26) + rrt(128) * density(36) 
  pd(32,36) = pd(32,36) + rrt(128) * density(26) 
  pd(36,26) = pd(36,26) - rrt(128) * density(36) 
  pd(36,36) = pd(36,36) - rrt(128) * density(26) 
  pd(53,26) = pd(53,26) + rrt(128) * density(36) 
  pd(53,36) = pd(53,36) + rrt(128) * density(26) 
  pd(21,27) = pd(21,27) + rrt(129) * density(36) 
  pd(21,36) = pd(21,36) + rrt(129) * density(27) 
  pd(27,27) = pd(27,27) - rrt(129) * density(36) 
  pd(27,36) = pd(27,36) - rrt(129) * density(27) 
  pd(29,27) = pd(29,27) + rrt(129) * density(36) 
  pd(29,36) = pd(29,36) + rrt(129) * density(27) 
  pd(36,27) = pd(36,27) - rrt(129) * density(36) 
  pd(36,36) = pd(36,36) - rrt(129) * density(27) 
  pd(53,27) = pd(53,27) + rrt(129) * density(36) 
  pd(53,36) = pd(53,36) + rrt(129) * density(27) 
  pd(01,10) = pd(01,10) + rrt(130) * density(36) 
  pd(01,36) = pd(01,36) + rrt(130) * density(10) 
  pd(10,10) = pd(10,10) - rrt(130) * density(36) 
  pd(10,36) = pd(10,36) - rrt(130) * density(10) 
  pd(29,10) = pd(29,10) + rrt(130) * density(36) 
  pd(29,36) = pd(29,36) + rrt(130) * density(10) 
  pd(36,10) = pd(36,10) - rrt(130) * density(36) 
  pd(36,36) = pd(36,36) - rrt(130) * density(10) 
  pd(53,10) = pd(53,10) + rrt(130) * density(36) 
  pd(53,36) = pd(53,36) + rrt(130) * density(10) 
  pd(01,11) = pd(01,11) + rrt(131) * density(36) 
  pd(01,36) = pd(01,36) + rrt(131) * density(11) 
  pd(11,11) = pd(11,11) - rrt(131) * density(36) 
  pd(11,36) = pd(11,36) - rrt(131) * density(11) 
  pd(29,11) = pd(29,11) + rrt(131) * density(36) 
  pd(29,36) = pd(29,36) + rrt(131) * density(11) 
  pd(36,11) = pd(36,11) - rrt(131) * density(36) 
  pd(36,36) = pd(36,36) - rrt(131) * density(11) 
  pd(53,11) = pd(53,11) + rrt(131) * density(36) 
  pd(53,36) = pd(53,36) + rrt(131) * density(11) 
  pd(21,32) = pd(21,32) + rrt(132) * density(36) * 2.0d0
  pd(21,36) = pd(21,36) + rrt(132) * density(32) * 2.0d0
  pd(32,32) = pd(32,32) - rrt(132) * density(36) 
  pd(32,36) = pd(32,36) - rrt(132) * density(32) 
  pd(36,32) = pd(36,32) - rrt(132) * density(36) 
  pd(36,36) = pd(36,36) - rrt(132) * density(32) 
  pd(53,32) = pd(53,32) + rrt(132) * density(36) 
  pd(53,36) = pd(53,36) + rrt(132) * density(32) 
  pd(29,29) = pd(29,29) - rrt(133) * density(37) 
  pd(29,37) = pd(29,37) - rrt(133) * density(29) 
  pd(32,29) = pd(32,29) + rrt(133) * density(37) 
  pd(32,37) = pd(32,37) + rrt(133) * density(29) 
  pd(37,29) = pd(37,29) - rrt(133) * density(37) 
  pd(37,37) = pd(37,37) - rrt(133) * density(29) 
  pd(53,29) = pd(53,29) + rrt(133) * density(37) 
  pd(53,37) = pd(53,37) + rrt(133) * density(29) 
  pd(14,14) = pd(14,14) - rrt(134) * density(37) 
  pd(14,37) = pd(14,37) - rrt(134) * density(14) 
  pd(37,14) = pd(37,14) - rrt(134) * density(37) 
  pd(37,37) = pd(37,37) - rrt(134) * density(14) 
  pd(42,14) = pd(42,14) + rrt(134) * density(37) 
  pd(42,37) = pd(42,37) + rrt(134) * density(14) 
  pd(53,14) = pd(53,14) + rrt(134) * density(37) 
  pd(53,37) = pd(53,37) + rrt(134) * density(14) 
  pd(21,21) = pd(21,21) + rrt(135) * density(37) 
  pd(21,37) = pd(21,37) + rrt(135) * density(21) 
  pd(37,21) = pd(37,21) - rrt(135) * density(37) 
  pd(37,37) = pd(37,37) - rrt(135) * density(21) 
  pd(53,21) = pd(53,21) + rrt(135) * density(37) 
  pd(53,37) = pd(53,37) + rrt(135) * density(21) 
  pd(21,26) = pd(21,26) + rrt(136) * density(37) * 2.0d0
  pd(21,37) = pd(21,37) + rrt(136) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(136) * density(37) 
  pd(26,37) = pd(26,37) - rrt(136) * density(26) 
  pd(37,26) = pd(37,26) - rrt(136) * density(37) 
  pd(37,37) = pd(37,37) - rrt(136) * density(26) 
  pd(53,26) = pd(53,26) + rrt(136) * density(37) 
  pd(53,37) = pd(53,37) + rrt(136) * density(26) 
  pd(21,27) = pd(21,27) + rrt(137) * density(37) * 2.0d0
  pd(21,37) = pd(21,37) + rrt(137) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(137) * density(37) 
  pd(27,37) = pd(27,37) - rrt(137) * density(27) 
  pd(37,27) = pd(37,27) - rrt(137) * density(37) 
  pd(37,37) = pd(37,37) - rrt(137) * density(27) 
  pd(53,27) = pd(53,27) + rrt(137) * density(37) 
  pd(53,37) = pd(53,37) + rrt(137) * density(27) 
  pd(21,01) = pd(21,01) + rrt(138) * density(37) 
  pd(21,37) = pd(21,37) + rrt(138) * density(01) 
  pd(37,01) = pd(37,01) - rrt(138) * density(37) 
  pd(37,37) = pd(37,37) - rrt(138) * density(01) 
  pd(53,01) = pd(53,01) + rrt(138) * density(37) 
  pd(53,37) = pd(53,37) + rrt(138) * density(01) 
  pd(01,10) = pd(01,10) + rrt(139) * density(37) 
  pd(01,37) = pd(01,37) + rrt(139) * density(10) 
  pd(10,10) = pd(10,10) - rrt(139) * density(37) 
  pd(10,37) = pd(10,37) - rrt(139) * density(10) 
  pd(21,10) = pd(21,10) + rrt(139) * density(37) 
  pd(21,37) = pd(21,37) + rrt(139) * density(10) 
  pd(37,10) = pd(37,10) - rrt(139) * density(37) 
  pd(37,37) = pd(37,37) - rrt(139) * density(10) 
  pd(53,10) = pd(53,10) + rrt(139) * density(37) 
  pd(53,37) = pd(53,37) + rrt(139) * density(10) 
  pd(01,11) = pd(01,11) + rrt(140) * density(37) 
  pd(01,37) = pd(01,37) + rrt(140) * density(11) 
  pd(11,11) = pd(11,11) - rrt(140) * density(37) 
  pd(11,37) = pd(11,37) - rrt(140) * density(11) 
  pd(21,11) = pd(21,11) + rrt(140) * density(37) 
  pd(21,37) = pd(21,37) + rrt(140) * density(11) 
  pd(37,11) = pd(37,11) - rrt(140) * density(37) 
  pd(37,37) = pd(37,37) - rrt(140) * density(11) 
  pd(53,11) = pd(53,11) + rrt(140) * density(37) 
  pd(53,37) = pd(53,37) + rrt(140) * density(11) 
  pd(21,29) = pd(21,29) + rrt(141) * density(38) * 2.0d0
  pd(21,38) = pd(21,38) + rrt(141) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(141) * density(38) 
  pd(29,38) = pd(29,38) - rrt(141) * density(29) 
  pd(38,29) = pd(38,29) - rrt(141) * density(38) 
  pd(38,38) = pd(38,38) - rrt(141) * density(29) 
  pd(53,29) = pd(53,29) + rrt(141) * density(38) 
  pd(53,38) = pd(53,38) + rrt(141) * density(29) 
  pd(14,14) = pd(14,14) - rrt(142) * density(48) 
  pd(14,48) = pd(14,48) - rrt(142) * density(14) 
  pd(41,14) = pd(41,14) + rrt(142) * density(48) 
  pd(41,48) = pd(41,48) + rrt(142) * density(14) 
  pd(48,14) = pd(48,14) - rrt(142) * density(48) 
  pd(48,48) = pd(48,48) - rrt(142) * density(14) 
  pd(53,14) = pd(53,14) + rrt(142) * density(48) 
  pd(53,48) = pd(53,48) + rrt(142) * density(14) 
  pd(14,14) = pd(14,14) - rrt(143) * density(38) 
  pd(14,38) = pd(14,38) - rrt(143) * density(14) 
  pd(21,14) = pd(21,14) + rrt(143) * density(38) 
  pd(21,38) = pd(21,38) + rrt(143) * density(14) 
  pd(38,14) = pd(38,14) - rrt(143) * density(38) 
  pd(38,38) = pd(38,38) - rrt(143) * density(14) 
  pd(40,14) = pd(40,14) + rrt(143) * density(38) 
  pd(40,38) = pd(40,38) + rrt(143) * density(14) 
  pd(53,14) = pd(53,14) + rrt(143) * density(38) 
  pd(53,38) = pd(53,38) + rrt(143) * density(14) 
  pd(01,14) = pd(01,14) + rrt(144) * density(49) 
  pd(01,49) = pd(01,49) + rrt(144) * density(14) 
  pd(14,14) = pd(14,14) - rrt(144) * density(49) 
  pd(14,49) = pd(14,49) - rrt(144) * density(14) 
  pd(40,14) = pd(40,14) + rrt(144) * density(49) 
  pd(40,49) = pd(40,49) + rrt(144) * density(14) 
  pd(49,14) = pd(49,14) - rrt(144) * density(49) 
  pd(49,49) = pd(49,49) - rrt(144) * density(14) 
  pd(53,14) = pd(53,14) + rrt(144) * density(49) 
  pd(53,49) = pd(53,49) + rrt(144) * density(14) 
  pd(14,14) = pd(14,14) - rrt(145) * density(50) 
  pd(14,50) = pd(14,50) - rrt(145) * density(14) 
  pd(40,14) = pd(40,14) + rrt(145) * density(50) * 2.0d0
  pd(40,50) = pd(40,50) + rrt(145) * density(14) * 2.0d0
  pd(50,14) = pd(50,14) - rrt(145) * density(50) 
  pd(50,50) = pd(50,50) - rrt(145) * density(14) 
  pd(53,14) = pd(53,14) + rrt(145) * density(50) 
  pd(53,50) = pd(53,50) + rrt(145) * density(14) 
  pd(14,14) = pd(14,14) - rrt(146) * density(51) 
  pd(14,51) = pd(14,51) - rrt(146) * density(14) 
  pd(40,14) = pd(40,14) + rrt(146) * density(51) 
  pd(40,51) = pd(40,51) + rrt(146) * density(14) 
  pd(42,14) = pd(42,14) + rrt(146) * density(51) 
  pd(42,51) = pd(42,51) + rrt(146) * density(14) 
  pd(51,14) = pd(51,14) - rrt(146) * density(51) 
  pd(51,51) = pd(51,51) - rrt(146) * density(14) 
  pd(53,14) = pd(53,14) + rrt(146) * density(51) 
  pd(53,51) = pd(53,51) + rrt(146) * density(14) 
  pd(29,29) = pd(29,29) - rrt(147) * density(48) 
  pd(29,48) = pd(29,48) - rrt(147) * density(29) 
  pd(42,29) = pd(42,29) + rrt(147) * density(48) 
  pd(42,48) = pd(42,48) + rrt(147) * density(29) 
  pd(48,29) = pd(48,29) - rrt(147) * density(48) 
  pd(48,48) = pd(48,48) - rrt(147) * density(29) 
  pd(53,29) = pd(53,29) + rrt(147) * density(48) 
  pd(53,48) = pd(53,48) + rrt(147) * density(29) 
  pd(29,29) = pd(29,29) - rrt(148) * density(49) 
  pd(29,49) = pd(29,49) - rrt(148) * density(29) 
  pd(40,29) = pd(40,29) + rrt(148) * density(49) * 2.0d0
  pd(40,49) = pd(40,49) + rrt(148) * density(29) * 2.0d0
  pd(49,29) = pd(49,29) - rrt(148) * density(49) 
  pd(49,49) = pd(49,49) - rrt(148) * density(29) 
  pd(53,29) = pd(53,29) + rrt(148) * density(49) 
  pd(53,49) = pd(53,49) + rrt(148) * density(29) 
  pd(21,29) = pd(21,29) + rrt(149) * density(50) 
  pd(21,50) = pd(21,50) + rrt(149) * density(29) 
  pd(29,29) = pd(29,29) - rrt(149) * density(50) 
  pd(29,50) = pd(29,50) - rrt(149) * density(29) 
  pd(40,29) = pd(40,29) + rrt(149) * density(50) 
  pd(40,50) = pd(40,50) + rrt(149) * density(29) 
  pd(50,29) = pd(50,29) - rrt(149) * density(50) 
  pd(50,50) = pd(50,50) - rrt(149) * density(29) 
  pd(53,29) = pd(53,29) + rrt(149) * density(50) 
  pd(53,50) = pd(53,50) + rrt(149) * density(29) 
  pd(29,29) = pd(29,29) - rrt(150) * density(51) 
  pd(29,51) = pd(29,51) - rrt(150) * density(29) 
  pd(32,29) = pd(32,29) + rrt(150) * density(51) 
  pd(32,51) = pd(32,51) + rrt(150) * density(29) 
  pd(40,29) = pd(40,29) + rrt(150) * density(51) 
  pd(40,51) = pd(40,51) + rrt(150) * density(29) 
  pd(51,29) = pd(51,29) - rrt(150) * density(51) 
  pd(51,51) = pd(51,51) - rrt(150) * density(29) 
  pd(53,29) = pd(53,29) + rrt(150) * density(51) 
  pd(53,51) = pd(53,51) + rrt(150) * density(29) 
  pd(01,10) = pd(01,10) + rrt(151) * density(38) 
  pd(01,38) = pd(01,38) + rrt(151) * density(10) 
  pd(10,10) = pd(10,10) - rrt(151) * density(38) 
  pd(10,38) = pd(10,38) - rrt(151) * density(10) 
  pd(32,10) = pd(32,10) + rrt(151) * density(38) 
  pd(32,38) = pd(32,38) + rrt(151) * density(10) 
  pd(38,10) = pd(38,10) - rrt(151) * density(38) 
  pd(38,38) = pd(38,38) - rrt(151) * density(10) 
  pd(53,10) = pd(53,10) + rrt(151) * density(38) 
  pd(53,38) = pd(53,38) + rrt(151) * density(10) 
  pd(01,10) = pd(01,10) + rrt(152) * density(48) 
  pd(01,48) = pd(01,48) + rrt(152) * density(10) 
  pd(10,10) = pd(10,10) - rrt(152) * density(48) 
  pd(10,48) = pd(10,48) - rrt(152) * density(10) 
  pd(40,10) = pd(40,10) + rrt(152) * density(48) 
  pd(40,48) = pd(40,48) + rrt(152) * density(10) 
  pd(48,10) = pd(48,10) - rrt(152) * density(48) 
  pd(48,48) = pd(48,48) - rrt(152) * density(10) 
  pd(53,10) = pd(53,10) + rrt(152) * density(48) 
  pd(53,48) = pd(53,48) + rrt(152) * density(10) 
  pd(01,10) = pd(01,10) + rrt(153) * density(49) 
  pd(01,49) = pd(01,49) + rrt(153) * density(10) 
  pd(10,10) = pd(10,10) - rrt(153) * density(49) 
  pd(10,49) = pd(10,49) - rrt(153) * density(10) 
  pd(41,10) = pd(41,10) + rrt(153) * density(49) 
  pd(41,49) = pd(41,49) + rrt(153) * density(10) 
  pd(49,10) = pd(49,10) - rrt(153) * density(49) 
  pd(49,49) = pd(49,49) - rrt(153) * density(10) 
  pd(53,10) = pd(53,10) + rrt(153) * density(49) 
  pd(53,49) = pd(53,49) + rrt(153) * density(10) 
  pd(01,10) = pd(01,10) + rrt(154) * density(50) 
  pd(01,50) = pd(01,50) + rrt(154) * density(10) 
  pd(10,10) = pd(10,10) - rrt(154) * density(50) 
  pd(10,50) = pd(10,50) - rrt(154) * density(10) 
  pd(42,10) = pd(42,10) + rrt(154) * density(50) 
  pd(42,50) = pd(42,50) + rrt(154) * density(10) 
  pd(50,10) = pd(50,10) - rrt(154) * density(50) 
  pd(50,50) = pd(50,50) - rrt(154) * density(10) 
  pd(53,10) = pd(53,10) + rrt(154) * density(50) 
  pd(53,50) = pd(53,50) + rrt(154) * density(10) 
  pd(01,10) = pd(01,10) + rrt(155) * density(51) 
  pd(01,51) = pd(01,51) + rrt(155) * density(10) 
  pd(10,10) = pd(10,10) - rrt(155) * density(51) 
  pd(10,51) = pd(10,51) - rrt(155) * density(10) 
  pd(43,10) = pd(43,10) + rrt(155) * density(51) 
  pd(43,51) = pd(43,51) + rrt(155) * density(10) 
  pd(51,10) = pd(51,10) - rrt(155) * density(51) 
  pd(51,51) = pd(51,51) - rrt(155) * density(10) 
  pd(53,10) = pd(53,10) + rrt(155) * density(51) 
  pd(53,51) = pd(53,51) + rrt(155) * density(10) 
  pd(01,11) = pd(01,11) + rrt(156) * density(38) 
  pd(01,38) = pd(01,38) + rrt(156) * density(11) 
  pd(11,11) = pd(11,11) - rrt(156) * density(38) 
  pd(11,38) = pd(11,38) - rrt(156) * density(11) 
  pd(32,11) = pd(32,11) + rrt(156) * density(38) 
  pd(32,38) = pd(32,38) + rrt(156) * density(11) 
  pd(38,11) = pd(38,11) - rrt(156) * density(38) 
  pd(38,38) = pd(38,38) - rrt(156) * density(11) 
  pd(53,11) = pd(53,11) + rrt(156) * density(38) 
  pd(53,38) = pd(53,38) + rrt(156) * density(11) 
  pd(01,11) = pd(01,11) + rrt(157) * density(48) 
  pd(01,48) = pd(01,48) + rrt(157) * density(11) 
  pd(11,11) = pd(11,11) - rrt(157) * density(48) 
  pd(11,48) = pd(11,48) - rrt(157) * density(11) 
  pd(40,11) = pd(40,11) + rrt(157) * density(48) 
  pd(40,48) = pd(40,48) + rrt(157) * density(11) 
  pd(48,11) = pd(48,11) - rrt(157) * density(48) 
  pd(48,48) = pd(48,48) - rrt(157) * density(11) 
  pd(53,11) = pd(53,11) + rrt(157) * density(48) 
  pd(53,48) = pd(53,48) + rrt(157) * density(11) 
  pd(01,11) = pd(01,11) + rrt(158) * density(49) 
  pd(01,49) = pd(01,49) + rrt(158) * density(11) 
  pd(11,11) = pd(11,11) - rrt(158) * density(49) 
  pd(11,49) = pd(11,49) - rrt(158) * density(11) 
  pd(41,11) = pd(41,11) + rrt(158) * density(49) 
  pd(41,49) = pd(41,49) + rrt(158) * density(11) 
  pd(49,11) = pd(49,11) - rrt(158) * density(49) 
  pd(49,49) = pd(49,49) - rrt(158) * density(11) 
  pd(53,11) = pd(53,11) + rrt(158) * density(49) 
  pd(53,49) = pd(53,49) + rrt(158) * density(11) 
  pd(01,11) = pd(01,11) + rrt(159) * density(50) 
  pd(01,50) = pd(01,50) + rrt(159) * density(11) 
  pd(11,11) = pd(11,11) - rrt(159) * density(50) 
  pd(11,50) = pd(11,50) - rrt(159) * density(11) 
  pd(42,11) = pd(42,11) + rrt(159) * density(50) 
  pd(42,50) = pd(42,50) + rrt(159) * density(11) 
  pd(50,11) = pd(50,11) - rrt(159) * density(50) 
  pd(50,50) = pd(50,50) - rrt(159) * density(11) 
  pd(53,11) = pd(53,11) + rrt(159) * density(50) 
  pd(53,50) = pd(53,50) + rrt(159) * density(11) 
  pd(01,11) = pd(01,11) + rrt(160) * density(51) 
  pd(01,51) = pd(01,51) + rrt(160) * density(11) 
  pd(11,11) = pd(11,11) - rrt(160) * density(51) 
  pd(11,51) = pd(11,51) - rrt(160) * density(11) 
  pd(43,11) = pd(43,11) + rrt(160) * density(51) 
  pd(43,51) = pd(43,51) + rrt(160) * density(11) 
  pd(51,11) = pd(51,11) - rrt(160) * density(51) 
  pd(51,51) = pd(51,51) - rrt(160) * density(11) 
  pd(53,11) = pd(53,11) + rrt(160) * density(51) 
  pd(53,51) = pd(53,51) + rrt(160) * density(11) 
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
  pd(10,10) = pd(10,10) - rrt(176) * density(21) 
  pd(10,21) = pd(10,21) - rrt(176) * density(10) 
  pd(21,10) = pd(21,10) - rrt(176) * density(21) 
  pd(21,21) = pd(21,21) - rrt(176) * density(10) 
  pd(29,10) = pd(29,10) + rrt(176) * density(21) 
  pd(29,21) = pd(29,21) + rrt(176) * density(10) 
  pd(41,10) = pd(41,10) + rrt(176) * density(21) 
  pd(41,21) = pd(41,21) + rrt(176) * density(10) 
  pd(01,01) = pd(01,01) + rrt(177) * density(10) 
  pd(01,10) = pd(01,10) + rrt(177) * density(01) 
  pd(10,01) = pd(10,01) - rrt(177) * density(10) 
  pd(10,10) = pd(10,10) - rrt(177) * density(01) 
  pd(01,10) = pd(01,10) + rrt(178) * density(40) 
  pd(01,40) = pd(01,40) + rrt(178) * density(10) 
  pd(10,10) = pd(10,10) - rrt(178) * density(40) 
  pd(10,40) = pd(10,40) - rrt(178) * density(10) 
  pd(01,10) = pd(01,10) + rrt(179) * density(41) 
  pd(01,41) = pd(01,41) + rrt(179) * density(10) 
  pd(10,10) = pd(10,10) - rrt(179) * density(41) 
  pd(10,41) = pd(10,41) - rrt(179) * density(10) 
  pd(14,10) = pd(14,10) + rrt(179) * density(41) 
  pd(14,41) = pd(14,41) + rrt(179) * density(10) 
  pd(40,10) = pd(40,10) + rrt(179) * density(41) 
  pd(40,41) = pd(40,41) + rrt(179) * density(10) 
  pd(41,10) = pd(41,10) - rrt(179) * density(41) 
  pd(41,41) = pd(41,41) - rrt(179) * density(10) 
  pd(01,10) = pd(01,10) + rrt(180) * density(42) 
  pd(01,42) = pd(01,42) + rrt(180) * density(10) 
  pd(10,10) = pd(10,10) - rrt(180) * density(42) 
  pd(10,42) = pd(10,42) - rrt(180) * density(10) 
  pd(29,10) = pd(29,10) + rrt(180) * density(42) 
  pd(29,42) = pd(29,42) + rrt(180) * density(10) 
  pd(40,10) = pd(40,10) + rrt(180) * density(42) 
  pd(40,42) = pd(40,42) + rrt(180) * density(10) 
  pd(42,10) = pd(42,10) - rrt(180) * density(42) 
  pd(42,42) = pd(42,42) - rrt(180) * density(10) 
  pd(01,10) = pd(01,10) + rrt(181) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(181) * density(10) * 4.0d0
  pd(11,10) = pd(11,10) + rrt(181) * density(10) * 2.0d0
  pd(01,10) = pd(01,10) + rrt(182) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(182) * density(10) * 4.0d0
  pd(13,10) = pd(13,10) + rrt(182) * density(10) * 2.0d0
  pd(10,01) = pd(10,01) + rrt(183) * density(11) 
  pd(10,11) = pd(10,11) + rrt(183) * density(01) 
  pd(11,01) = pd(11,01) - rrt(183) * density(11) 
  pd(11,11) = pd(11,11) - rrt(183) * density(01) 
  pd(01,01) = pd(01,01) + rrt(184) * density(11) 
  pd(01,11) = pd(01,11) + rrt(184) * density(01) 
  pd(11,01) = pd(11,01) - rrt(184) * density(11) 
  pd(11,11) = pd(11,11) - rrt(184) * density(01) 
  pd(01,11) = pd(01,11) + rrt(185) * density(21) 
  pd(01,21) = pd(01,21) + rrt(185) * density(11) 
  pd(11,11) = pd(11,11) - rrt(185) * density(21) 
  pd(11,21) = pd(11,21) - rrt(185) * density(11) 
  pd(21,11) = pd(21,11) - rrt(185) * density(21) 
  pd(21,21) = pd(21,21) - rrt(185) * density(11) 
  pd(29,11) = pd(29,11) + rrt(185) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(185) * density(11) * 2.0d0
  pd(10,11) = pd(10,11) + rrt(186) * density(40) 
  pd(10,40) = pd(10,40) + rrt(186) * density(11) 
  pd(11,11) = pd(11,11) - rrt(186) * density(40) 
  pd(11,40) = pd(11,40) - rrt(186) * density(11) 
  pd(12,01) = pd(12,01) + rrt(187) * density(13) 
  pd(12,13) = pd(12,13) + rrt(187) * density(01) 
  pd(13,01) = pd(13,01) - rrt(187) * density(13) 
  pd(13,13) = pd(13,13) - rrt(187) * density(01) 
  pd(01,13) = pd(01,13) + rrt(188) * density(21) 
  pd(01,21) = pd(01,21) + rrt(188) * density(13) 
  pd(13,13) = pd(13,13) - rrt(188) * density(21) 
  pd(13,21) = pd(13,21) - rrt(188) * density(13) 
  pd(21,13) = pd(21,13) - rrt(188) * density(21) 
  pd(21,21) = pd(21,21) - rrt(188) * density(13) 
  pd(29,13) = pd(29,13) + rrt(188) * density(21) 
  pd(29,21) = pd(29,21) + rrt(188) * density(13) 
  pd(31,13) = pd(31,13) + rrt(188) * density(21) 
  pd(31,21) = pd(31,21) + rrt(188) * density(13) 
  pd(11,01) = pd(11,01) + rrt(189) * density(12) 
  pd(11,12) = pd(11,12) + rrt(189) * density(01) 
  pd(12,01) = pd(12,01) - rrt(189) * density(12) 
  pd(12,12) = pd(12,12) - rrt(189) * density(01) 
  pd(01,12) = pd(01,12) + rrt(190) * density(21) 
  pd(01,21) = pd(01,21) + rrt(190) * density(12) 
  pd(12,12) = pd(12,12) - rrt(190) * density(21) 
  pd(12,21) = pd(12,21) - rrt(190) * density(12) 
  pd(21,12) = pd(21,12) - rrt(190) * density(21) 
  pd(21,21) = pd(21,21) - rrt(190) * density(12) 
  pd(29,12) = pd(29,12) + rrt(190) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(190) * density(12) * 2.0d0
  pd(01,12) = pd(01,12) + rrt(191) * density(40) 
  pd(01,40) = pd(01,40) + rrt(191) * density(12) 
  pd(12,12) = pd(12,12) - rrt(191) * density(40) 
  pd(12,40) = pd(12,40) - rrt(191) * density(12) 
  pd(14,12) = pd(14,12) + rrt(191) * density(40) 
  pd(14,40) = pd(14,40) + rrt(191) * density(12) 
  pd(29,12) = pd(29,12) + rrt(191) * density(40) 
  pd(29,40) = pd(29,40) + rrt(191) * density(12) 
  pd(40,12) = pd(40,12) - rrt(191) * density(40) 
  pd(40,40) = pd(40,40) - rrt(191) * density(12) 
  pd(10,10) = pd(10,10) - rrt(192) * density(12) 
  pd(10,12) = pd(10,12) - rrt(192) * density(10) 
  pd(12,10) = pd(12,10) - rrt(192) * density(12) 
  pd(12,12) = pd(12,12) - rrt(192) * density(10) 
  pd(20,10) = pd(20,10) + rrt(192) * density(12) 
  pd(20,12) = pd(20,12) + rrt(192) * density(10) 
  pd(53,10) = pd(53,10) + rrt(192) * density(12) 
  pd(53,12) = pd(53,12) + rrt(192) * density(10) 
  pd(12,12) = pd(12,12) - rrt(193) * density(12) * 4.0d0
  pd(20,12) = pd(20,12) + rrt(193) * density(12) * 2.0d0
  pd(53,12) = pd(53,12) + rrt(193) * density(12) * 2.0d0
  pd(10,01) = pd(10,01) + rrt(194) * density(14)**2 
  pd(10,14) = pd(10,14) + rrt(194) * density(01) * density(14) * 2.0d0
  pd(14,01) = pd(14,01) - rrt(194) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(194) * density(01) * density(14) * 4.0d0
  pd(10,14) = pd(10,14) + rrt(195) * density(14) * density(21) * 2.0d0
  pd(10,21) = pd(10,21) + rrt(195) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(195) * density(14) * density(21) * 4.0d0
  pd(14,21) = pd(14,21) - rrt(195) * density(14)**2 * 2.0d0
  pd(10,14) = pd(10,14) + rrt(196) * density(14) * density(40) * 2.0d0
  pd(10,40) = pd(10,40) + rrt(196) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(196) * density(14) * density(40) * 4.0d0
  pd(14,40) = pd(14,40) - rrt(196) * density(14)**2 * 2.0d0
  pd(10,14) = pd(10,14) + rrt(197) * density(14)**2 * 3.0d0
  pd(14,14) = pd(14,14) - rrt(197) * density(14)**2 * 6.0d0
  pd(10,14) = pd(10,14) + rrt(198) * density(14) * density(29) * 2.0d0
  pd(10,29) = pd(10,29) + rrt(198) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(198) * density(14) * density(29) * 4.0d0
  pd(14,29) = pd(14,29) - rrt(198) * density(14)**2 * 2.0d0
  pd(11,01) = pd(11,01) + rrt(199) * density(14)**2 
  pd(11,14) = pd(11,14) + rrt(199) * density(01) * density(14) * 2.0d0
  pd(14,01) = pd(14,01) - rrt(199) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(199) * density(01) * density(14) * 4.0d0
  pd(11,14) = pd(11,14) + rrt(200) * density(14) * density(21) * 2.0d0
  pd(11,21) = pd(11,21) + rrt(200) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(200) * density(14) * density(21) * 4.0d0
  pd(14,21) = pd(14,21) - rrt(200) * density(14)**2 * 2.0d0
  pd(11,14) = pd(11,14) + rrt(201) * density(14) * density(40) * 2.0d0
  pd(11,40) = pd(11,40) + rrt(201) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(201) * density(14) * density(40) * 4.0d0
  pd(14,40) = pd(14,40) - rrt(201) * density(14)**2 * 2.0d0
  pd(11,14) = pd(11,14) + rrt(202) * density(14)**2 * 3.0d0
  pd(14,14) = pd(14,14) - rrt(202) * density(14)**2 * 6.0d0
  pd(11,14) = pd(11,14) + rrt(203) * density(14) * density(29) * 2.0d0
  pd(11,29) = pd(11,29) + rrt(203) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(203) * density(14) * density(29) * 4.0d0
  pd(14,29) = pd(14,29) - rrt(203) * density(14)**2 * 2.0d0
  pd(14,15) = pd(14,15) + rrt(204) * density(29) 
  pd(14,29) = pd(14,29) + rrt(204) * density(15) 
  pd(15,15) = pd(15,15) - rrt(204) * density(29) 
  pd(15,29) = pd(15,29) - rrt(204) * density(15) 
  pd(29,15) = pd(29,15) - rrt(204) * density(29) 
  pd(29,29) = pd(29,29) - rrt(204) * density(15) 
  pd(30,15) = pd(30,15) + rrt(204) * density(29) 
  pd(30,29) = pd(30,29) + rrt(204) * density(15) 
  pd(15,15) = pd(15,15) - rrt(205) * density(21) 
  pd(15,21) = pd(15,21) - rrt(205) * density(15) 
  pd(21,15) = pd(21,15) - rrt(205) * density(21) 
  pd(21,21) = pd(21,21) - rrt(205) * density(15) 
  pd(29,15) = pd(29,15) + rrt(205) * density(21) 
  pd(29,21) = pd(29,21) + rrt(205) * density(15) 
  pd(40,15) = pd(40,15) + rrt(205) * density(21) 
  pd(40,21) = pd(40,21) + rrt(205) * density(15) 
  pd(01,15) = pd(01,15) + rrt(206) * density(40) 
  pd(01,40) = pd(01,40) + rrt(206) * density(15) 
  pd(15,15) = pd(15,15) - rrt(206) * density(40) 
  pd(15,40) = pd(15,40) - rrt(206) * density(15) 
  pd(29,15) = pd(29,15) + rrt(206) * density(40) 
  pd(29,40) = pd(29,40) + rrt(206) * density(15) 
  pd(40,15) = pd(40,15) - rrt(206) * density(40) 
  pd(40,40) = pd(40,40) - rrt(206) * density(15) 
  pd(01,15) = pd(01,15) + rrt(207) * density(41) 
  pd(01,41) = pd(01,41) + rrt(207) * density(15) 
  pd(15,15) = pd(15,15) - rrt(207) * density(41) 
  pd(15,41) = pd(15,41) - rrt(207) * density(15) 
  pd(40,15) = pd(40,15) + rrt(207) * density(41) 
  pd(40,41) = pd(40,41) + rrt(207) * density(15) 
  pd(41,15) = pd(41,15) - rrt(207) * density(41) 
  pd(41,41) = pd(41,41) - rrt(207) * density(15) 
  pd(14,01) = pd(14,01) + rrt(208) * density(15) 
  pd(14,15) = pd(14,15) + rrt(208) * density(01) 
  pd(15,01) = pd(15,01) - rrt(208) * density(15) 
  pd(15,15) = pd(15,15) - rrt(208) * density(01) 
  pd(14,14) = pd(14,14) + rrt(209) * density(16) 
  pd(14,16) = pd(14,16) + rrt(209) * density(14) 
  pd(16,14) = pd(16,14) - rrt(209) * density(16) 
  pd(16,16) = pd(16,16) - rrt(209) * density(14) 
  pd(14,16) = pd(14,16) + rrt(210) * density(29) 
  pd(14,29) = pd(14,29) + rrt(210) * density(16) 
  pd(16,16) = pd(16,16) - rrt(210) * density(29) 
  pd(16,29) = pd(16,29) - rrt(210) * density(16) 
  pd(15,14) = pd(15,14) + rrt(211) * density(16) 
  pd(15,16) = pd(15,16) + rrt(211) * density(14) 
  pd(16,14) = pd(16,14) - rrt(211) * density(16) 
  pd(16,16) = pd(16,16) - rrt(211) * density(14) 
  pd(14,01) = pd(14,01) + rrt(212) * density(16) 
  pd(14,16) = pd(14,16) + rrt(212) * density(01) 
  pd(16,01) = pd(16,01) - rrt(212) * density(16) 
  pd(16,16) = pd(16,16) - rrt(212) * density(01) 
  pd(15,15) = pd(15,15) - rrt(213) * density(16) 
  pd(15,16) = pd(15,16) - rrt(213) * density(15) 
  pd(16,15) = pd(16,15) - rrt(213) * density(16) 
  pd(16,16) = pd(16,16) - rrt(213) * density(15) 
  pd(18,15) = pd(18,15) + rrt(213) * density(16) 
  pd(18,16) = pd(18,16) + rrt(213) * density(15) 
  pd(53,15) = pd(53,15) + rrt(213) * density(16) 
  pd(53,16) = pd(53,16) + rrt(213) * density(15) 
  pd(16,16) = pd(16,16) - rrt(214) * density(21) 
  pd(16,21) = pd(16,21) - rrt(214) * density(16) 
  pd(21,16) = pd(21,16) - rrt(214) * density(21) 
  pd(21,21) = pd(21,21) - rrt(214) * density(16) 
  pd(29,16) = pd(29,16) + rrt(214) * density(21) 
  pd(29,21) = pd(29,21) + rrt(214) * density(16) 
  pd(40,16) = pd(40,16) + rrt(214) * density(21) 
  pd(40,21) = pd(40,21) + rrt(214) * density(16) 
  pd(10,16) = pd(10,16) + rrt(215) * density(40) 
  pd(10,40) = pd(10,40) + rrt(215) * density(16) 
  pd(16,16) = pd(16,16) - rrt(215) * density(40) 
  pd(16,40) = pd(16,40) - rrt(215) * density(16) 
  pd(29,16) = pd(29,16) + rrt(215) * density(40) 
  pd(29,40) = pd(29,40) + rrt(215) * density(16) 
  pd(40,16) = pd(40,16) - rrt(215) * density(40) 
  pd(40,40) = pd(40,40) - rrt(215) * density(16) 
  pd(21,26) = pd(21,26) + rrt(216) * density(29) 
  pd(21,29) = pd(21,29) + rrt(216) * density(26) 
  pd(26,26) = pd(26,26) - rrt(216) * density(29) 
  pd(26,29) = pd(26,29) - rrt(216) * density(26) 
  pd(14,14) = pd(14,14) - rrt(217) * density(26) 
  pd(14,26) = pd(14,26) - rrt(217) * density(14) 
  pd(26,14) = pd(26,14) - rrt(217) * density(26) 
  pd(26,26) = pd(26,26) - rrt(217) * density(14) 
  pd(29,14) = pd(29,14) + rrt(217) * density(26) 
  pd(29,26) = pd(29,26) + rrt(217) * density(14) 
  pd(40,14) = pd(40,14) + rrt(217) * density(26) 
  pd(40,26) = pd(40,26) + rrt(217) * density(14) 
  pd(21,21) = pd(21,21) + rrt(218) * density(26) 
  pd(21,26) = pd(21,26) + rrt(218) * density(21) 
  pd(26,21) = pd(26,21) - rrt(218) * density(26) 
  pd(26,26) = pd(26,26) - rrt(218) * density(21) 
  pd(21,01) = pd(21,01) + rrt(219) * density(26) 
  pd(21,26) = pd(21,26) + rrt(219) * density(01) 
  pd(26,01) = pd(26,01) - rrt(219) * density(26) 
  pd(26,26) = pd(26,26) - rrt(219) * density(01) 
  pd(21,26) = pd(21,26) + rrt(220) * density(40) 
  pd(21,40) = pd(21,40) + rrt(220) * density(26) 
  pd(26,26) = pd(26,26) - rrt(220) * density(40) 
  pd(26,40) = pd(26,40) - rrt(220) * density(26) 
  pd(21,26) = pd(21,26) + rrt(221) * density(32) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(221) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(221) * density(32) 
  pd(26,32) = pd(26,32) - rrt(221) * density(26) 
  pd(30,26) = pd(30,26) + rrt(221) * density(32) 
  pd(30,32) = pd(30,32) + rrt(221) * density(26) 
  pd(32,26) = pd(32,26) - rrt(221) * density(32) 
  pd(32,32) = pd(32,32) - rrt(221) * density(26) 
  pd(21,26) = pd(21,26) + rrt(222) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(222) * density(26) * 4.0d0
  pd(27,26) = pd(27,26) + rrt(222) * density(26) * 2.0d0
  pd(21,29) = pd(21,29) + rrt(223) * density(32) 
  pd(21,32) = pd(21,32) + rrt(223) * density(29) 
  pd(26,29) = pd(26,29) + rrt(223) * density(32) 
  pd(26,32) = pd(26,32) + rrt(223) * density(29) 
  pd(29,29) = pd(29,29) - rrt(223) * density(32) 
  pd(29,32) = pd(29,32) - rrt(223) * density(29) 
  pd(32,29) = pd(32,29) - rrt(223) * density(32) 
  pd(32,32) = pd(32,32) - rrt(223) * density(29) 
  pd(26,27) = pd(26,27) + rrt(224) * density(29) 
  pd(26,29) = pd(26,29) + rrt(224) * density(27) 
  pd(27,27) = pd(27,27) - rrt(224) * density(29) 
  pd(27,29) = pd(27,29) - rrt(224) * density(27) 
  pd(21,27) = pd(21,27) + rrt(225) * density(29) 
  pd(21,29) = pd(21,29) + rrt(225) * density(27) 
  pd(27,27) = pd(27,27) - rrt(225) * density(29) 
  pd(27,29) = pd(27,29) - rrt(225) * density(27) 
  pd(29,27) = pd(29,27) - rrt(225) * density(29) 
  pd(29,29) = pd(29,29) - rrt(225) * density(27) 
  pd(30,27) = pd(30,27) + rrt(225) * density(29) 
  pd(30,29) = pd(30,29) + rrt(225) * density(27) 
  pd(26,21) = pd(26,21) + rrt(226) * density(27) 
  pd(26,27) = pd(26,27) + rrt(226) * density(21) 
  pd(27,21) = pd(27,21) - rrt(226) * density(27) 
  pd(27,27) = pd(27,27) - rrt(226) * density(21) 
  pd(26,01) = pd(26,01) + rrt(227) * density(27) 
  pd(26,27) = pd(26,27) + rrt(227) * density(01) 
  pd(27,01) = pd(27,01) - rrt(227) * density(27) 
  pd(27,27) = pd(27,27) - rrt(227) * density(01) 
  pd(26,27) = pd(26,27) + rrt(228) * density(40) 
  pd(26,40) = pd(26,40) + rrt(228) * density(27) 
  pd(27,27) = pd(27,27) - rrt(228) * density(40) 
  pd(27,40) = pd(27,40) - rrt(228) * density(27) 
  pd(21,27) = pd(21,27) + rrt(229) * density(32) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(229) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(229) * density(32) 
  pd(27,32) = pd(27,32) - rrt(229) * density(27) 
  pd(29,27) = pd(29,27) + rrt(229) * density(32) 
  pd(29,32) = pd(29,32) + rrt(229) * density(27) 
  pd(32,27) = pd(32,27) - rrt(229) * density(32) 
  pd(32,32) = pd(32,32) - rrt(229) * density(27) 
  pd(21,28) = pd(21,28) + rrt(230) * density(29) 
  pd(21,29) = pd(21,29) + rrt(230) * density(28) 
  pd(28,28) = pd(28,28) - rrt(230) * density(29) 
  pd(28,29) = pd(28,29) - rrt(230) * density(28) 
  pd(29,28) = pd(29,28) - rrt(230) * density(29) 
  pd(29,29) = pd(29,29) - rrt(230) * density(28) 
  pd(31,28) = pd(31,28) + rrt(230) * density(29) 
  pd(31,29) = pd(31,29) + rrt(230) * density(28) 
  pd(21,21) = pd(21,21) - rrt(231) * density(28) 
  pd(21,28) = pd(21,28) - rrt(231) * density(21) 
  pd(27,21) = pd(27,21) + rrt(231) * density(28) * 2.0d0
  pd(27,28) = pd(27,28) + rrt(231) * density(21) * 2.0d0
  pd(28,21) = pd(28,21) - rrt(231) * density(28) 
  pd(28,28) = pd(28,28) - rrt(231) * density(21) 
  pd(27,01) = pd(27,01) + rrt(232) * density(28) 
  pd(27,28) = pd(27,28) + rrt(232) * density(01) 
  pd(28,01) = pd(28,01) - rrt(232) * density(28) 
  pd(28,28) = pd(28,28) - rrt(232) * density(01) 
  pd(29,29) = pd(29,29) + rrt(233) * density(30) 
  pd(29,30) = pd(29,30) + rrt(233) * density(29) 
  pd(30,29) = pd(30,29) - rrt(233) * density(30) 
  pd(30,30) = pd(30,30) - rrt(233) * density(29) 
  pd(29,21) = pd(29,21) + rrt(234) * density(30) 
  pd(29,30) = pd(29,30) + rrt(234) * density(21) 
  pd(30,21) = pd(30,21) - rrt(234) * density(30) 
  pd(30,30) = pd(30,30) - rrt(234) * density(21) 
  pd(21,21) = pd(21,21) - rrt(235) * density(30) 
  pd(21,30) = pd(21,30) - rrt(235) * density(21) 
  pd(26,21) = pd(26,21) + rrt(235) * density(30) 
  pd(26,30) = pd(26,30) + rrt(235) * density(21) 
  pd(29,21) = pd(29,21) + rrt(235) * density(30) 
  pd(29,30) = pd(29,30) + rrt(235) * density(21) 
  pd(30,21) = pd(30,21) - rrt(235) * density(30) 
  pd(30,30) = pd(30,30) - rrt(235) * density(21) 
  pd(21,21) = pd(21,21) - rrt(236) * density(30) 
  pd(21,30) = pd(21,30) - rrt(236) * density(21) 
  pd(27,21) = pd(27,21) + rrt(236) * density(30) 
  pd(27,30) = pd(27,30) + rrt(236) * density(21) 
  pd(29,21) = pd(29,21) + rrt(236) * density(30) 
  pd(29,30) = pd(29,30) + rrt(236) * density(21) 
  pd(30,21) = pd(30,21) - rrt(236) * density(30) 
  pd(30,30) = pd(30,30) - rrt(236) * density(21) 
  pd(29,01) = pd(29,01) + rrt(237) * density(30) 
  pd(29,30) = pd(29,30) + rrt(237) * density(01) 
  pd(30,01) = pd(30,01) - rrt(237) * density(30) 
  pd(30,30) = pd(30,30) - rrt(237) * density(01) 
  pd(21,30) = pd(21,30) + rrt(238) * density(32) 
  pd(21,32) = pd(21,32) + rrt(238) * density(30) 
  pd(29,30) = pd(29,30) + rrt(238) * density(32) * 2.0d0
  pd(29,32) = pd(29,32) + rrt(238) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(238) * density(32) 
  pd(30,32) = pd(30,32) - rrt(238) * density(30) 
  pd(32,30) = pd(32,30) - rrt(238) * density(32) 
  pd(32,32) = pd(32,32) - rrt(238) * density(30) 
  pd(21,30) = pd(21,30) + rrt(239) * density(32) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(239) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(239) * density(32) 
  pd(30,32) = pd(30,32) - rrt(239) * density(30) 
  pd(32,30) = pd(32,30) - rrt(239) * density(32) 
  pd(32,32) = pd(32,32) - rrt(239) * density(30) 
  pd(14,30) = pd(14,30) + rrt(240) * density(40) 
  pd(14,40) = pd(14,40) + rrt(240) * density(30) 
  pd(21,30) = pd(21,30) + rrt(240) * density(40) 
  pd(21,40) = pd(21,40) + rrt(240) * density(30) 
  pd(30,30) = pd(30,30) - rrt(240) * density(40) 
  pd(30,40) = pd(30,40) - rrt(240) * density(30) 
  pd(40,30) = pd(40,30) - rrt(240) * density(40) 
  pd(40,40) = pd(40,40) - rrt(240) * density(30) 
  pd(30,30) = pd(30,30) - rrt(241) * density(41) 
  pd(30,41) = pd(30,41) - rrt(241) * density(30) 
  pd(40,30) = pd(40,30) + rrt(241) * density(41) * 2.0d0
  pd(40,41) = pd(40,41) + rrt(241) * density(30) * 2.0d0
  pd(41,30) = pd(41,30) - rrt(241) * density(41) 
  pd(41,41) = pd(41,41) - rrt(241) * density(30) 
  pd(01,30) = pd(01,30) + rrt(242) * density(41) 
  pd(01,41) = pd(01,41) + rrt(242) * density(30) 
  pd(21,30) = pd(21,30) + rrt(242) * density(41) 
  pd(21,41) = pd(21,41) + rrt(242) * density(30) 
  pd(30,30) = pd(30,30) - rrt(242) * density(41) 
  pd(30,41) = pd(30,41) - rrt(242) * density(30) 
  pd(41,30) = pd(41,30) - rrt(242) * density(41) 
  pd(41,41) = pd(41,41) - rrt(242) * density(30) 
  pd(30,29) = pd(30,29) + rrt(243) * density(31) 
  pd(30,31) = pd(30,31) + rrt(243) * density(29) 
  pd(31,29) = pd(31,29) - rrt(243) * density(31) 
  pd(31,31) = pd(31,31) - rrt(243) * density(29) 
  pd(29,14) = pd(29,14) + rrt(244) * density(31) 
  pd(29,31) = pd(29,31) + rrt(244) * density(14) 
  pd(31,14) = pd(31,14) - rrt(244) * density(31) 
  pd(31,31) = pd(31,31) - rrt(244) * density(14) 
  pd(30,21) = pd(30,21) + rrt(245) * density(31) 
  pd(30,31) = pd(30,31) + rrt(245) * density(21) 
  pd(31,21) = pd(31,21) - rrt(245) * density(31) 
  pd(31,31) = pd(31,31) - rrt(245) * density(21) 
  pd(21,21) = pd(21,21) - rrt(246) * density(31) 
  pd(21,31) = pd(21,31) - rrt(246) * density(21) 
  pd(29,21) = pd(29,21) + rrt(246) * density(31) * 3.0d0
  pd(29,31) = pd(29,31) + rrt(246) * density(21) * 3.0d0
  pd(31,21) = pd(31,21) - rrt(246) * density(31) 
  pd(31,31) = pd(31,31) - rrt(246) * density(21) 
  pd(29,01) = pd(29,01) + rrt(247) * density(31) 
  pd(29,31) = pd(29,31) + rrt(247) * density(01) 
  pd(31,01) = pd(31,01) - rrt(247) * density(31) 
  pd(31,31) = pd(31,31) - rrt(247) * density(01) 
  pd(26,26) = pd(26,26) - rrt(248) * density(31) 
  pd(26,31) = pd(26,31) - rrt(248) * density(26) 
  pd(28,26) = pd(28,26) + rrt(248) * density(31) 
  pd(28,31) = pd(28,31) + rrt(248) * density(26) 
  pd(29,26) = pd(29,26) + rrt(248) * density(31) 
  pd(29,31) = pd(29,31) + rrt(248) * density(26) 
  pd(31,26) = pd(31,26) - rrt(248) * density(31) 
  pd(31,31) = pd(31,31) - rrt(248) * density(26) 
  pd(26,26) = pd(26,26) - rrt(249) * density(31) 
  pd(26,31) = pd(26,31) - rrt(249) * density(26) 
  pd(27,26) = pd(27,26) + rrt(249) * density(31) 
  pd(27,31) = pd(27,31) + rrt(249) * density(26) 
  pd(30,26) = pd(30,26) + rrt(249) * density(31) 
  pd(30,31) = pd(30,31) + rrt(249) * density(26) 
  pd(31,26) = pd(31,26) - rrt(249) * density(31) 
  pd(31,31) = pd(31,31) - rrt(249) * density(26) 
  pd(26,26) = pd(26,26) - rrt(250) * density(31) 
  pd(26,31) = pd(26,31) - rrt(250) * density(26) 
  pd(29,26) = pd(29,26) + rrt(250) * density(31) * 3.0d0
  pd(29,31) = pd(29,31) + rrt(250) * density(26) * 3.0d0
  pd(31,26) = pd(31,26) - rrt(250) * density(31) 
  pd(31,31) = pd(31,31) - rrt(250) * density(26) 
  pd(29,31) = pd(29,31) + rrt(251) * density(40) 
  pd(29,40) = pd(29,40) + rrt(251) * density(31) 
  pd(31,31) = pd(31,31) - rrt(251) * density(40) 
  pd(31,40) = pd(31,40) - rrt(251) * density(31) 
  pd(30,31) = pd(30,31) + rrt(252) * density(40) 
  pd(30,40) = pd(30,40) + rrt(252) * density(31) 
  pd(31,31) = pd(31,31) - rrt(252) * density(40) 
  pd(31,40) = pd(31,40) - rrt(252) * density(31) 
  pd(21,31) = pd(21,31) + rrt(253) * density(32) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(253) * density(31) * 2.0d0
  pd(31,31) = pd(31,31) - rrt(253) * density(32) 
  pd(31,32) = pd(31,32) - rrt(253) * density(31) 
  pd(32,31) = pd(32,31) - rrt(253) * density(32) 
  pd(32,32) = pd(32,32) - rrt(253) * density(31) 
  pd(21,31) = pd(21,31) + rrt(254) * density(32) 
  pd(21,32) = pd(21,32) + rrt(254) * density(31) 
  pd(29,31) = pd(29,31) + rrt(254) * density(32) 
  pd(29,32) = pd(29,32) + rrt(254) * density(31) 
  pd(30,31) = pd(30,31) + rrt(254) * density(32) 
  pd(30,32) = pd(30,32) + rrt(254) * density(31) 
  pd(31,31) = pd(31,31) - rrt(254) * density(32) 
  pd(31,32) = pd(31,32) - rrt(254) * density(31) 
  pd(32,31) = pd(32,31) - rrt(254) * density(32) 
  pd(32,32) = pd(32,32) - rrt(254) * density(31) 
  pd(29,31) = pd(29,31) + rrt(255) * density(41) 
  pd(29,41) = pd(29,41) + rrt(255) * density(31) 
  pd(31,31) = pd(31,31) - rrt(255) * density(41) 
  pd(31,41) = pd(31,41) - rrt(255) * density(31) 
  pd(30,31) = pd(30,31) + rrt(256) * density(41) 
  pd(30,41) = pd(30,41) + rrt(256) * density(31) 
  pd(31,31) = pd(31,31) - rrt(256) * density(41) 
  pd(31,41) = pd(31,41) - rrt(256) * density(31) 
  pd(01,14) = pd(01,14) + rrt(257) * density(40) 
  pd(01,40) = pd(01,40) + rrt(257) * density(14) 
  pd(14,14) = pd(14,14) - rrt(257) * density(40) 
  pd(14,40) = pd(14,40) - rrt(257) * density(14) 
  pd(29,14) = pd(29,14) + rrt(257) * density(40) 
  pd(29,40) = pd(29,40) + rrt(257) * density(14) 
  pd(40,14) = pd(40,14) - rrt(257) * density(40) 
  pd(40,40) = pd(40,40) - rrt(257) * density(14) 
  pd(14,14) = pd(14,14) - rrt(258) * density(21) 
  pd(14,21) = pd(14,21) - rrt(258) * density(14) 
  pd(21,14) = pd(21,14) - rrt(258) * density(21) 
  pd(21,21) = pd(21,21) - rrt(258) * density(14) 
  pd(29,14) = pd(29,14) + rrt(258) * density(21) 
  pd(29,21) = pd(29,21) + rrt(258) * density(14) 
  pd(40,14) = pd(40,14) + rrt(258) * density(21) 
  pd(40,21) = pd(40,21) + rrt(258) * density(14) 
  pd(01,14) = pd(01,14) + rrt(259) * density(42) 
  pd(01,42) = pd(01,42) + rrt(259) * density(14) 
  pd(14,14) = pd(14,14) - rrt(259) * density(42) 
  pd(14,42) = pd(14,42) - rrt(259) * density(14) 
  pd(29,14) = pd(29,14) + rrt(259) * density(42) * 2.0d0
  pd(29,42) = pd(29,42) + rrt(259) * density(14) * 2.0d0
  pd(42,14) = pd(42,14) - rrt(259) * density(42) 
  pd(42,42) = pd(42,42) - rrt(259) * density(14) 
  pd(14,14) = pd(14,14) - rrt(260) * density(42) 
  pd(14,42) = pd(14,42) - rrt(260) * density(14) 
  pd(29,14) = pd(29,14) + rrt(260) * density(42) 
  pd(29,42) = pd(29,42) + rrt(260) * density(14) 
  pd(41,14) = pd(41,14) + rrt(260) * density(42) 
  pd(41,42) = pd(41,42) + rrt(260) * density(14) 
  pd(42,14) = pd(42,14) - rrt(260) * density(42) 
  pd(42,42) = pd(42,42) - rrt(260) * density(14) 
  pd(01,14) = pd(01,14) + rrt(261) * density(42) 
  pd(01,42) = pd(01,42) + rrt(261) * density(14) 
  pd(14,14) = pd(14,14) - rrt(261) * density(42) 
  pd(14,42) = pd(14,42) - rrt(261) * density(14) 
  pd(21,14) = pd(21,14) + rrt(261) * density(42) 
  pd(21,42) = pd(21,42) + rrt(261) * density(14) 
  pd(42,14) = pd(42,14) - rrt(261) * density(42) 
  pd(42,42) = pd(42,42) - rrt(261) * density(14) 
  pd(14,14) = pd(14,14) - rrt(262) * density(42) 
  pd(14,42) = pd(14,42) - rrt(262) * density(14) 
  pd(40,14) = pd(40,14) + rrt(262) * density(42) * 2.0d0
  pd(40,42) = pd(40,42) + rrt(262) * density(14) * 2.0d0
  pd(42,14) = pd(42,14) - rrt(262) * density(42) 
  pd(42,42) = pd(42,42) - rrt(262) * density(14) 
  pd(01,01) = pd(01,01) - rrt(263) * density(29) 
  pd(01,29) = pd(01,29) - rrt(263) * density(01) 
  pd(14,01) = pd(14,01) + rrt(263) * density(29) 
  pd(14,29) = pd(14,29) + rrt(263) * density(01) 
  pd(29,01) = pd(29,01) - rrt(263) * density(29) 
  pd(29,29) = pd(29,29) - rrt(263) * density(01) 
  pd(40,01) = pd(40,01) + rrt(263) * density(29) 
  pd(40,29) = pd(40,29) + rrt(263) * density(01) 
  pd(14,29) = pd(14,29) + rrt(264) * density(40) 
  pd(14,40) = pd(14,40) + rrt(264) * density(29) 
  pd(21,29) = pd(21,29) + rrt(264) * density(40) 
  pd(21,40) = pd(21,40) + rrt(264) * density(29) 
  pd(29,29) = pd(29,29) - rrt(264) * density(40) 
  pd(29,40) = pd(29,40) - rrt(264) * density(29) 
  pd(40,29) = pd(40,29) - rrt(264) * density(40) 
  pd(40,40) = pd(40,40) - rrt(264) * density(29) 
  pd(29,29) = pd(29,29) - rrt(265) * density(40) 
  pd(29,40) = pd(29,40) - rrt(265) * density(29) 
  pd(40,29) = pd(40,29) - rrt(265) * density(40) 
  pd(40,40) = pd(40,40) - rrt(265) * density(29) 
  pd(42,29) = pd(42,29) + rrt(265) * density(40) 
  pd(42,40) = pd(42,40) + rrt(265) * density(29) 
  pd(01,29) = pd(01,29) + rrt(266) * density(41) 
  pd(01,41) = pd(01,41) + rrt(266) * density(29) 
  pd(21,29) = pd(21,29) + rrt(266) * density(41) 
  pd(21,41) = pd(21,41) + rrt(266) * density(29) 
  pd(29,29) = pd(29,29) - rrt(266) * density(41) 
  pd(29,41) = pd(29,41) - rrt(266) * density(29) 
  pd(41,29) = pd(41,29) - rrt(266) * density(41) 
  pd(41,41) = pd(41,41) - rrt(266) * density(29) 
  pd(29,29) = pd(29,29) - rrt(267) * density(41) 
  pd(29,41) = pd(29,41) - rrt(267) * density(29) 
  pd(40,29) = pd(40,29) + rrt(267) * density(41) * 2.0d0
  pd(40,41) = pd(40,41) + rrt(267) * density(29) * 2.0d0
  pd(41,29) = pd(41,29) - rrt(267) * density(41) 
  pd(41,41) = pd(41,41) - rrt(267) * density(29) 
  pd(21,29) = pd(21,29) + rrt(268) * density(42) 
  pd(21,42) = pd(21,42) + rrt(268) * density(29) 
  pd(29,29) = pd(29,29) - rrt(268) * density(42) 
  pd(29,42) = pd(29,42) - rrt(268) * density(29) 
  pd(40,29) = pd(40,29) + rrt(268) * density(42) 
  pd(40,42) = pd(40,42) + rrt(268) * density(29) 
  pd(42,29) = pd(42,29) - rrt(268) * density(42) 
  pd(42,42) = pd(42,42) - rrt(268) * density(29) 
  pd(21,29) = pd(21,29) + rrt(269) * density(43) 
  pd(21,43) = pd(21,43) + rrt(269) * density(29) 
  pd(29,29) = pd(29,29) - rrt(269) * density(43) 
  pd(29,43) = pd(29,43) - rrt(269) * density(29) 
  pd(42,29) = pd(42,29) + rrt(269) * density(43) 
  pd(42,43) = pd(42,43) + rrt(269) * density(29) 
  pd(43,29) = pd(43,29) - rrt(269) * density(43) 
  pd(43,43) = pd(43,43) - rrt(269) * density(29) 
  pd(01,01) = pd(01,01) - rrt(270) * density(21) 
  pd(01,21) = pd(01,21) - rrt(270) * density(01) 
  pd(21,01) = pd(21,01) - rrt(270) * density(21) 
  pd(21,21) = pd(21,21) - rrt(270) * density(01) 
  pd(29,01) = pd(29,01) + rrt(270) * density(21) 
  pd(29,21) = pd(29,21) + rrt(270) * density(01) 
  pd(41,01) = pd(41,01) + rrt(270) * density(21) 
  pd(41,21) = pd(41,21) + rrt(270) * density(01) 
  pd(14,40) = pd(14,40) + rrt(271) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(271) * density(40) * 4.0d0
  pd(42,40) = pd(42,40) + rrt(271) * density(40) * 2.0d0
  pd(29,40) = pd(29,40) + rrt(272) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(272) * density(40) * 4.0d0
  pd(41,40) = pd(41,40) + rrt(272) * density(40) * 2.0d0
  pd(01,40) = pd(01,40) + rrt(273) * density(40) * 2.0d0
  pd(21,40) = pd(21,40) + rrt(273) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(273) * density(40) * 4.0d0
  pd(21,21) = pd(21,21) - rrt(274) * density(40) 
  pd(21,40) = pd(21,40) - rrt(274) * density(21) 
  pd(29,21) = pd(29,21) + rrt(274) * density(40) 
  pd(29,40) = pd(29,40) + rrt(274) * density(21) 
  pd(40,21) = pd(40,21) - rrt(274) * density(40) 
  pd(40,40) = pd(40,40) - rrt(274) * density(21) 
  pd(42,21) = pd(42,21) + rrt(274) * density(40) 
  pd(42,40) = pd(42,40) + rrt(274) * density(21) 
  pd(21,32) = pd(21,32) + rrt(275) * density(40) 
  pd(21,40) = pd(21,40) + rrt(275) * density(32) 
  pd(32,32) = pd(32,32) - rrt(275) * density(40) 
  pd(32,40) = pd(32,40) - rrt(275) * density(32) 
  pd(40,32) = pd(40,32) - rrt(275) * density(40) 
  pd(40,40) = pd(40,40) - rrt(275) * density(32) 
  pd(42,32) = pd(42,32) + rrt(275) * density(40) 
  pd(42,40) = pd(42,40) + rrt(275) * density(32) 
  pd(01,40) = pd(01,40) + rrt(276) * density(41) 
  pd(01,41) = pd(01,41) + rrt(276) * density(40) 
  pd(40,40) = pd(40,40) - rrt(276) * density(41) 
  pd(40,41) = pd(40,41) - rrt(276) * density(40) 
  pd(41,40) = pd(41,40) - rrt(276) * density(41) 
  pd(41,41) = pd(41,41) - rrt(276) * density(40) 
  pd(42,40) = pd(42,40) + rrt(276) * density(41) 
  pd(42,41) = pd(42,41) + rrt(276) * density(40) 
  pd(40,40) = pd(40,40) - rrt(277) * density(43) 
  pd(40,43) = pd(40,43) - rrt(277) * density(40) 
  pd(42,40) = pd(42,40) + rrt(277) * density(43) * 2.0d0
  pd(42,43) = pd(42,43) + rrt(277) * density(40) * 2.0d0
  pd(43,40) = pd(43,40) - rrt(277) * density(43) 
  pd(43,43) = pd(43,43) - rrt(277) * density(40) 
  pd(21,21) = pd(21,21) - rrt(278) * density(21) * 4.0d0
  pd(29,21) = pd(29,21) + rrt(278) * density(21) * 2.0d0
  pd(32,21) = pd(32,21) + rrt(278) * density(21) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(279) * density(42) 
  pd(21,42) = pd(21,42) - rrt(279) * density(21) 
  pd(32,21) = pd(32,21) + rrt(279) * density(42) 
  pd(32,42) = pd(32,42) + rrt(279) * density(21) 
  pd(40,21) = pd(40,21) + rrt(279) * density(42) 
  pd(40,42) = pd(40,42) + rrt(279) * density(21) 
  pd(42,21) = pd(42,21) - rrt(279) * density(42) 
  pd(42,42) = pd(42,42) - rrt(279) * density(21) 
  pd(21,42) = pd(21,42) + rrt(280) * density(42) * 2.0d0
  pd(40,42) = pd(40,42) + rrt(280) * density(42) * 4.0d0
  pd(42,42) = pd(42,42) - rrt(280) * density(42) * 4.0d0
  pd(40,42) = pd(40,42) + rrt(281) * density(42) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(281) * density(42) * 4.0d0
  pd(43,42) = pd(43,42) + rrt(281) * density(42) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(282) * density(42) 
  pd(21,42) = pd(21,42) + rrt(282) * density(32) 
  pd(32,32) = pd(32,32) - rrt(282) * density(42) 
  pd(32,42) = pd(32,42) - rrt(282) * density(32) 
  pd(42,32) = pd(42,32) - rrt(282) * density(42) 
  pd(42,42) = pd(42,42) - rrt(282) * density(32) 
  pd(43,32) = pd(43,32) + rrt(282) * density(42) 
  pd(43,42) = pd(43,42) + rrt(282) * density(32) 
  pd(21,42) = pd(21,42) + rrt(283) * density(43) 
  pd(21,43) = pd(21,43) + rrt(283) * density(42) 
  pd(40,42) = pd(40,42) + rrt(283) * density(43) 
  pd(40,43) = pd(40,43) + rrt(283) * density(42) 
  pd(43,42) = pd(43,42) - rrt(283) * density(43) 
  pd(43,43) = pd(43,43) - rrt(283) * density(42) 
  pd(21,21) = pd(21,21) - rrt(284) * density(43) 
  pd(21,43) = pd(21,43) - rrt(284) * density(21) 
  pd(32,21) = pd(32,21) + rrt(284) * density(43) 
  pd(32,43) = pd(32,43) + rrt(284) * density(21) 
  pd(42,21) = pd(42,21) + rrt(284) * density(43) 
  pd(42,43) = pd(42,43) + rrt(284) * density(21) 
  pd(43,21) = pd(43,21) - rrt(284) * density(43) 
  pd(43,43) = pd(43,43) - rrt(284) * density(21) 
  pd(21,43) = pd(21,43) + rrt(285) * density(43) * 2.0d0
  pd(42,43) = pd(42,43) + rrt(285) * density(43) * 4.0d0
  pd(43,43) = pd(43,43) - rrt(285) * density(43) * 4.0d0
  pd(14,14) = pd(14,14) - rrt(286) * density(14) * 4.0d0
  pd(18,14) = pd(18,14) + rrt(286) * density(14) * 2.0d0
  pd(53,14) = pd(53,14) + rrt(286) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(287) * density(29) 
  pd(14,29) = pd(14,29) - rrt(287) * density(14) 
  pd(29,14) = pd(29,14) - rrt(287) * density(29) 
  pd(29,29) = pd(29,29) - rrt(287) * density(14) 
  pd(45,14) = pd(45,14) + rrt(287) * density(29) 
  pd(45,29) = pd(45,29) + rrt(287) * density(14) 
  pd(53,14) = pd(53,14) + rrt(287) * density(29) 
  pd(53,29) = pd(53,29) + rrt(287) * density(14) 
  pd(01,01) = pd(01,01) - rrt(288) * density(01) * 2.0d0
  pd(14,01) = pd(14,01) + rrt(288) * density(01) * 4.0d0
  pd(01,01) = pd(01,01) - rrt(289) * density(21) 
  pd(01,21) = pd(01,21) - rrt(289) * density(01) 
  pd(14,01) = pd(14,01) + rrt(289) * density(21) * 2.0d0
  pd(14,21) = pd(14,21) + rrt(289) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(290) * density(40) 
  pd(01,40) = pd(01,40) - rrt(290) * density(01) 
  pd(14,01) = pd(14,01) + rrt(290) * density(40) * 2.0d0
  pd(14,40) = pd(14,40) + rrt(290) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(291) * density(29) 
  pd(01,29) = pd(01,29) - rrt(291) * density(01) 
  pd(14,01) = pd(14,01) + rrt(291) * density(29) * 2.0d0
  pd(14,29) = pd(14,29) + rrt(291) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(292) * density(14) 
  pd(01,14) = pd(01,14) - rrt(292) * density(01) 
  pd(14,01) = pd(14,01) + rrt(292) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) + rrt(292) * density(01) * 2.0d0
  pd(21,01) = pd(21,01) - rrt(293) * density(21) 
  pd(21,21) = pd(21,21) - rrt(293) * density(01) 
  pd(29,01) = pd(29,01) + rrt(293) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(293) * density(01) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(294) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(294) * density(21) * 4.0d0
  pd(21,21) = pd(21,21) - rrt(295) * density(29) 
  pd(21,29) = pd(21,29) - rrt(295) * density(21) 
  pd(29,21) = pd(29,21) + rrt(295) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) + rrt(295) * density(21) * 2.0d0
  pd(21,14) = pd(21,14) - rrt(296) * density(21) 
  pd(21,21) = pd(21,21) - rrt(296) * density(14) 
  pd(29,14) = pd(29,14) + rrt(296) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(296) * density(14) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(297) * density(40) 
  pd(21,40) = pd(21,40) - rrt(297) * density(21) 
  pd(29,21) = pd(29,21) + rrt(297) * density(40) * 2.0d0
  pd(29,40) = pd(29,40) + rrt(297) * density(21) * 2.0d0
  pd(14,01) = pd(14,01) + rrt(298) * density(40) 
  pd(14,40) = pd(14,40) + rrt(298) * density(01) 
  pd(29,01) = pd(29,01) + rrt(298) * density(40) 
  pd(29,40) = pd(29,40) + rrt(298) * density(01) 
  pd(40,01) = pd(40,01) - rrt(298) * density(40) 
  pd(40,40) = pd(40,40) - rrt(298) * density(01) 
  pd(14,21) = pd(14,21) + rrt(299) * density(40) 
  pd(14,40) = pd(14,40) + rrt(299) * density(21) 
  pd(29,21) = pd(29,21) + rrt(299) * density(40) 
  pd(29,40) = pd(29,40) + rrt(299) * density(21) 
  pd(40,21) = pd(40,21) - rrt(299) * density(40) 
  pd(40,40) = pd(40,40) - rrt(299) * density(21) 
  pd(14,29) = pd(14,29) + rrt(300) * density(40) 
  pd(14,40) = pd(14,40) + rrt(300) * density(29) 
  pd(29,29) = pd(29,29) + rrt(300) * density(40) 
  pd(29,40) = pd(29,40) + rrt(300) * density(29) 
  pd(40,29) = pd(40,29) - rrt(300) * density(40) 
  pd(40,40) = pd(40,40) - rrt(300) * density(29) 
  pd(14,14) = pd(14,14) + rrt(301) * density(40) 
  pd(14,40) = pd(14,40) + rrt(301) * density(14) 
  pd(29,14) = pd(29,14) + rrt(301) * density(40) 
  pd(29,40) = pd(29,40) + rrt(301) * density(14) 
  pd(40,14) = pd(40,14) - rrt(301) * density(40) 
  pd(40,40) = pd(40,40) - rrt(301) * density(14) 
  pd(14,40) = pd(14,40) + rrt(302) * density(40) * 2.0d0
  pd(29,40) = pd(29,40) + rrt(302) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(302) * density(40) * 2.0d0
  pd(21,01) = pd(21,01) + rrt(303) * density(32) 
  pd(21,32) = pd(21,32) + rrt(303) * density(01) 
  pd(29,01) = pd(29,01) + rrt(303) * density(32) 
  pd(29,32) = pd(29,32) + rrt(303) * density(01) 
  pd(32,01) = pd(32,01) - rrt(303) * density(32) 
  pd(32,32) = pd(32,32) - rrt(303) * density(01) 
  pd(21,21) = pd(21,21) + rrt(304) * density(32) 
  pd(21,32) = pd(21,32) + rrt(304) * density(21) 
  pd(29,21) = pd(29,21) + rrt(304) * density(32) 
  pd(29,32) = pd(29,32) + rrt(304) * density(21) 
  pd(32,21) = pd(32,21) - rrt(304) * density(32) 
  pd(32,32) = pd(32,32) - rrt(304) * density(21) 
  pd(21,14) = pd(21,14) + rrt(305) * density(32) 
  pd(21,32) = pd(21,32) + rrt(305) * density(14) 
  pd(29,14) = pd(29,14) + rrt(305) * density(32) 
  pd(29,32) = pd(29,32) + rrt(305) * density(14) 
  pd(32,14) = pd(32,14) - rrt(305) * density(32) 
  pd(32,32) = pd(32,32) - rrt(305) * density(14) 
  pd(21,29) = pd(21,29) + rrt(306) * density(32) 
  pd(21,32) = pd(21,32) + rrt(306) * density(29) 
  pd(29,29) = pd(29,29) + rrt(306) * density(32) 
  pd(29,32) = pd(29,32) + rrt(306) * density(29) 
  pd(32,29) = pd(32,29) - rrt(306) * density(32) 
  pd(32,32) = pd(32,32) - rrt(306) * density(29) 
  pd(01,01) = pd(01,01) + rrt(307) * density(41) 
  pd(01,41) = pd(01,41) + rrt(307) * density(01) 
  pd(29,01) = pd(29,01) + rrt(307) * density(41) 
  pd(29,41) = pd(29,41) + rrt(307) * density(01) 
  pd(41,01) = pd(41,01) - rrt(307) * density(41) 
  pd(41,41) = pd(41,41) - rrt(307) * density(01) 
  pd(01,21) = pd(01,21) + rrt(308) * density(41) 
  pd(01,41) = pd(01,41) + rrt(308) * density(21) 
  pd(29,21) = pd(29,21) + rrt(308) * density(41) 
  pd(29,41) = pd(29,41) + rrt(308) * density(21) 
  pd(41,21) = pd(41,21) - rrt(308) * density(41) 
  pd(41,41) = pd(41,41) - rrt(308) * density(21) 
  pd(01,40) = pd(01,40) + rrt(309) * density(41) 
  pd(01,41) = pd(01,41) + rrt(309) * density(40) 
  pd(29,40) = pd(29,40) + rrt(309) * density(41) 
  pd(29,41) = pd(29,41) + rrt(309) * density(40) 
  pd(41,40) = pd(41,40) - rrt(309) * density(41) 
  pd(41,41) = pd(41,41) - rrt(309) * density(40) 
  pd(01,41) = pd(01,41) + rrt(310) * density(41) * 2.0d0
  pd(29,41) = pd(29,41) + rrt(310) * density(41) * 2.0d0
  pd(41,41) = pd(41,41) - rrt(310) * density(41) * 2.0d0
  pd(29,01) = pd(29,01) + rrt(311) * density(42) 
  pd(29,42) = pd(29,42) + rrt(311) * density(01) 
  pd(40,01) = pd(40,01) + rrt(311) * density(42) 
  pd(40,42) = pd(40,42) + rrt(311) * density(01) 
  pd(42,01) = pd(42,01) - rrt(311) * density(42) 
  pd(42,42) = pd(42,42) - rrt(311) * density(01) 
  pd(29,21) = pd(29,21) + rrt(312) * density(42) 
  pd(29,42) = pd(29,42) + rrt(312) * density(21) 
  pd(40,21) = pd(40,21) + rrt(312) * density(42) 
  pd(40,42) = pd(40,42) + rrt(312) * density(21) 
  pd(42,21) = pd(42,21) - rrt(312) * density(42) 
  pd(42,42) = pd(42,42) - rrt(312) * density(21) 
  pd(29,40) = pd(29,40) + rrt(313) * density(42) 
  pd(29,42) = pd(29,42) + rrt(313) * density(40) 
  pd(40,40) = pd(40,40) + rrt(313) * density(42) 
  pd(40,42) = pd(40,42) + rrt(313) * density(40) 
  pd(42,40) = pd(42,40) - rrt(313) * density(42) 
  pd(42,42) = pd(42,42) - rrt(313) * density(40) 
  pd(29,42) = pd(29,42) + rrt(314) * density(42) * 2.0d0
  pd(40,42) = pd(40,42) + rrt(314) * density(42) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(314) * density(42) * 2.0d0
  pd(29,01) = pd(29,01) + rrt(315) * density(43) 
  pd(29,43) = pd(29,43) + rrt(315) * density(01) 
  pd(42,01) = pd(42,01) + rrt(315) * density(43) 
  pd(42,43) = pd(42,43) + rrt(315) * density(01) 
  pd(43,01) = pd(43,01) - rrt(315) * density(43) 
  pd(43,43) = pd(43,43) - rrt(315) * density(01) 
  pd(29,21) = pd(29,21) + rrt(316) * density(43) 
  pd(29,43) = pd(29,43) + rrt(316) * density(21) 
  pd(42,21) = pd(42,21) + rrt(316) * density(43) 
  pd(42,43) = pd(42,43) + rrt(316) * density(21) 
  pd(43,21) = pd(43,21) - rrt(316) * density(43) 
  pd(43,43) = pd(43,43) - rrt(316) * density(21) 
  pd(29,40) = pd(29,40) + rrt(317) * density(43) 
  pd(29,43) = pd(29,43) + rrt(317) * density(40) 
  pd(42,40) = pd(42,40) + rrt(317) * density(43) 
  pd(42,43) = pd(42,43) + rrt(317) * density(40) 
  pd(43,40) = pd(43,40) - rrt(317) * density(43) 
  pd(43,43) = pd(43,43) - rrt(317) * density(40) 
  pd(29,14) = pd(29,14) + rrt(318) * density(43) 
  pd(29,43) = pd(29,43) + rrt(318) * density(14) 
  pd(42,14) = pd(42,14) + rrt(318) * density(43) 
  pd(42,43) = pd(42,43) + rrt(318) * density(14) 
  pd(43,14) = pd(43,14) - rrt(318) * density(43) 
  pd(43,43) = pd(43,43) - rrt(318) * density(14) 
  pd(29,29) = pd(29,29) + rrt(319) * density(43) 
  pd(29,43) = pd(29,43) + rrt(319) * density(29) 
  pd(42,29) = pd(42,29) + rrt(319) * density(43) 
  pd(42,43) = pd(42,43) + rrt(319) * density(29) 
  pd(43,29) = pd(43,29) - rrt(319) * density(43) 
  pd(43,43) = pd(43,43) - rrt(319) * density(29) 
  pd(21,01) = pd(21,01) + rrt(320) * density(43) 
  pd(21,43) = pd(21,43) + rrt(320) * density(01) 
  pd(40,01) = pd(40,01) + rrt(320) * density(43) 
  pd(40,43) = pd(40,43) + rrt(320) * density(01) 
  pd(43,01) = pd(43,01) - rrt(320) * density(43) 
  pd(43,43) = pd(43,43) - rrt(320) * density(01) 
  pd(21,21) = pd(21,21) + rrt(321) * density(43) 
  pd(21,43) = pd(21,43) + rrt(321) * density(21) 
  pd(40,21) = pd(40,21) + rrt(321) * density(43) 
  pd(40,43) = pd(40,43) + rrt(321) * density(21) 
  pd(43,21) = pd(43,21) - rrt(321) * density(43) 
  pd(43,43) = pd(43,43) - rrt(321) * density(21) 
  pd(21,40) = pd(21,40) + rrt(322) * density(43) 
  pd(21,43) = pd(21,43) + rrt(322) * density(40) 
  pd(40,40) = pd(40,40) + rrt(322) * density(43) 
  pd(40,43) = pd(40,43) + rrt(322) * density(40) 
  pd(43,40) = pd(43,40) - rrt(322) * density(43) 
  pd(43,43) = pd(43,43) - rrt(322) * density(40) 
  pd(21,14) = pd(21,14) + rrt(323) * density(43) 
  pd(21,43) = pd(21,43) + rrt(323) * density(14) 
  pd(40,14) = pd(40,14) + rrt(323) * density(43) 
  pd(40,43) = pd(40,43) + rrt(323) * density(14) 
  pd(43,14) = pd(43,14) - rrt(323) * density(43) 
  pd(43,43) = pd(43,43) - rrt(323) * density(14) 
  pd(21,29) = pd(21,29) + rrt(324) * density(43) 
  pd(21,43) = pd(21,43) + rrt(324) * density(29) 
  pd(40,29) = pd(40,29) + rrt(324) * density(43) 
  pd(40,43) = pd(40,43) + rrt(324) * density(29) 
  pd(43,29) = pd(43,29) - rrt(324) * density(43) 
  pd(43,43) = pd(43,43) - rrt(324) * density(29) 
  pd(42,44) = pd(42,44) + rrt(325) 
  pd(43,44) = pd(43,44) + rrt(325) 
  pd(44,44) = pd(44,44) - rrt(325) 
  pd(01,01) = pd(01,01) + rrt(326) * density(14)**2 
  pd(01,14) = pd(01,14) + rrt(326) * density(01) * density(14) * 2.0d0
  pd(14,01) = pd(14,01) - rrt(326) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(326) * density(01) * density(14) * 4.0d0
  pd(01,14) = pd(01,14) + rrt(327) * density(14) * density(21) * 2.0d0
  pd(01,21) = pd(01,21) + rrt(327) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(327) * density(14) * density(21) * 4.0d0
  pd(14,21) = pd(14,21) - rrt(327) * density(14)**2 * 2.0d0
  pd(01,14) = pd(01,14) + rrt(328) * density(14) * density(40) * 2.0d0
  pd(01,40) = pd(01,40) + rrt(328) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(328) * density(14) * density(40) * 4.0d0
  pd(14,40) = pd(14,40) - rrt(328) * density(14)**2 * 2.0d0
  pd(01,14) = pd(01,14) + rrt(329) * density(14)**2 * 3.0d0
  pd(14,14) = pd(14,14) - rrt(329) * density(14)**2 * 6.0d0
  pd(01,14) = pd(01,14) + rrt(330) * density(14) * density(29) * 2.0d0
  pd(01,29) = pd(01,29) + rrt(330) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(330) * density(14) * density(29) * 4.0d0
  pd(14,29) = pd(14,29) - rrt(330) * density(14)**2 * 2.0d0
  pd(21,01) = pd(21,01) + rrt(331) * density(29)**2 
  pd(21,29) = pd(21,29) + rrt(331) * density(01) * density(29) * 2.0d0
  pd(29,01) = pd(29,01) - rrt(331) * density(29)**2 * 2.0d0
  pd(29,29) = pd(29,29) - rrt(331) * density(01) * density(29) * 4.0d0
  pd(21,21) = pd(21,21) + rrt(332) * density(29)**2 
  pd(21,29) = pd(21,29) + rrt(332) * density(21) * density(29) * 2.0d0
  pd(29,21) = pd(29,21) - rrt(332) * density(29)**2 * 2.0d0
  pd(29,29) = pd(29,29) - rrt(332) * density(21) * density(29) * 4.0d0
  pd(21,14) = pd(21,14) + rrt(333) * density(29)**2 
  pd(21,29) = pd(21,29) + rrt(333) * density(14) * density(29) * 2.0d0
  pd(29,14) = pd(29,14) - rrt(333) * density(29)**2 * 2.0d0
  pd(29,29) = pd(29,29) - rrt(333) * density(14) * density(29) * 4.0d0
  pd(21,29) = pd(21,29) + rrt(334) * density(29)**2 * 3.0d0
  pd(29,29) = pd(29,29) - rrt(334) * density(29)**2 * 6.0d0
  pd(21,29) = pd(21,29) + rrt(335) * density(29) * density(40) * 2.0d0
  pd(21,40) = pd(21,40) + rrt(335) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(335) * density(29) * density(40) * 4.0d0
  pd(29,40) = pd(29,40) - rrt(335) * density(29)**2 * 2.0d0
  pd(14,01) = pd(14,01) - rrt(336) * density(14) * density(29) 
  pd(14,14) = pd(14,14) - rrt(336) * density(01) * density(29) 
  pd(14,29) = pd(14,29) - rrt(336) * density(01) * density(14) 
  pd(29,01) = pd(29,01) - rrt(336) * density(14) * density(29) 
  pd(29,14) = pd(29,14) - rrt(336) * density(01) * density(29) 
  pd(29,29) = pd(29,29) - rrt(336) * density(01) * density(14) 
  pd(40,01) = pd(40,01) + rrt(336) * density(14) * density(29) 
  pd(40,14) = pd(40,14) + rrt(336) * density(01) * density(29) 
  pd(40,29) = pd(40,29) + rrt(336) * density(01) * density(14) 
  pd(14,14) = pd(14,14) - rrt(337) * density(21) * density(29) 
  pd(14,21) = pd(14,21) - rrt(337) * density(14) * density(29) 
  pd(14,29) = pd(14,29) - rrt(337) * density(14) * density(21) 
  pd(29,14) = pd(29,14) - rrt(337) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(337) * density(14) * density(29) 
  pd(29,29) = pd(29,29) - rrt(337) * density(14) * density(21) 
  pd(40,14) = pd(40,14) + rrt(337) * density(21) * density(29) 
  pd(40,21) = pd(40,21) + rrt(337) * density(14) * density(29) 
  pd(40,29) = pd(40,29) + rrt(337) * density(14) * density(21) 
  pd(14,14) = pd(14,14) - rrt(338) * density(14) * density(29) * 2.0d0
  pd(14,29) = pd(14,29) - rrt(338) * density(14)**2 
  pd(29,14) = pd(29,14) - rrt(338) * density(14) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(338) * density(14)**2 
  pd(40,14) = pd(40,14) + rrt(338) * density(14) * density(29) * 2.0d0
  pd(40,29) = pd(40,29) + rrt(338) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(339) * density(29)**2 
  pd(14,29) = pd(14,29) - rrt(339) * density(14) * density(29) * 2.0d0
  pd(29,14) = pd(29,14) - rrt(339) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(339) * density(14) * density(29) * 2.0d0
  pd(40,14) = pd(40,14) + rrt(339) * density(29)**2 
  pd(40,29) = pd(40,29) + rrt(339) * density(14) * density(29) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(340) * density(29) * density(40) 
  pd(14,29) = pd(14,29) - rrt(340) * density(14) * density(40) 
  pd(14,40) = pd(14,40) - rrt(340) * density(14) * density(29) 
  pd(29,14) = pd(29,14) - rrt(340) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(340) * density(14) * density(40) 
  pd(29,40) = pd(29,40) - rrt(340) * density(14) * density(29) 
  pd(40,14) = pd(40,14) + rrt(340) * density(29) * density(40) 
  pd(40,29) = pd(40,29) + rrt(340) * density(14) * density(40) 
  pd(40,40) = pd(40,40) + rrt(340) * density(14) * density(29) 
  pd(21,01) = pd(21,01) - rrt(341) * density(21) * density(29) 
  pd(21,21) = pd(21,21) - rrt(341) * density(01) * density(29) 
  pd(21,29) = pd(21,29) - rrt(341) * density(01) * density(21) 
  pd(29,01) = pd(29,01) - rrt(341) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(341) * density(01) * density(29) 
  pd(29,29) = pd(29,29) - rrt(341) * density(01) * density(21) 
  pd(32,01) = pd(32,01) + rrt(341) * density(21) * density(29) 
  pd(32,21) = pd(32,21) + rrt(341) * density(01) * density(29) 
  pd(32,29) = pd(32,29) + rrt(341) * density(01) * density(21) 
  pd(21,21) = pd(21,21) - rrt(342) * density(21) * density(29) * 2.0d0
  pd(21,29) = pd(21,29) - rrt(342) * density(21)**2 
  pd(29,21) = pd(29,21) - rrt(342) * density(21) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(342) * density(21)**2 
  pd(32,21) = pd(32,21) + rrt(342) * density(21) * density(29) * 2.0d0
  pd(32,29) = pd(32,29) + rrt(342) * density(21)**2 
  pd(21,21) = pd(21,21) - rrt(343) * density(29) * density(40) 
  pd(21,29) = pd(21,29) - rrt(343) * density(21) * density(40) 
  pd(21,40) = pd(21,40) - rrt(343) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(343) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(343) * density(21) * density(40) 
  pd(29,40) = pd(29,40) - rrt(343) * density(21) * density(29) 
  pd(32,21) = pd(32,21) + rrt(343) * density(29) * density(40) 
  pd(32,29) = pd(32,29) + rrt(343) * density(21) * density(40) 
  pd(32,40) = pd(32,40) + rrt(343) * density(21) * density(29) 
  pd(21,14) = pd(21,14) - rrt(344) * density(21) * density(29) 
  pd(21,21) = pd(21,21) - rrt(344) * density(14) * density(29) 
  pd(21,29) = pd(21,29) - rrt(344) * density(14) * density(21) 
  pd(29,14) = pd(29,14) - rrt(344) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(344) * density(14) * density(29) 
  pd(29,29) = pd(29,29) - rrt(344) * density(14) * density(21) 
  pd(32,14) = pd(32,14) + rrt(344) * density(21) * density(29) 
  pd(32,21) = pd(32,21) + rrt(344) * density(14) * density(29) 
  pd(32,29) = pd(32,29) + rrt(344) * density(14) * density(21) 
  pd(21,21) = pd(21,21) - rrt(345) * density(29)**2 
  pd(21,29) = pd(21,29) - rrt(345) * density(21) * density(29) * 2.0d0
  pd(29,21) = pd(29,21) - rrt(345) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(345) * density(21) * density(29) * 2.0d0
  pd(32,21) = pd(32,21) + rrt(345) * density(29)**2 
  pd(32,29) = pd(32,29) + rrt(345) * density(21) * density(29) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(346) * density(29) 
  pd(01,29) = pd(01,29) - rrt(346) * density(01) 
  pd(29,01) = pd(29,01) - rrt(346) * density(29) 
  pd(29,29) = pd(29,29) - rrt(346) * density(01) 
  pd(41,01) = pd(41,01) + rrt(346) * density(29) 
  pd(41,29) = pd(41,29) + rrt(346) * density(01) 
  pd(29,01) = pd(29,01) - rrt(347) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(347) * density(01) * density(40) 
  pd(29,40) = pd(29,40) - rrt(347) * density(01) * density(29) 
  pd(40,01) = pd(40,01) - rrt(347) * density(29) * density(40) 
  pd(40,29) = pd(40,29) - rrt(347) * density(01) * density(40) 
  pd(40,40) = pd(40,40) - rrt(347) * density(01) * density(29) 
  pd(42,01) = pd(42,01) + rrt(347) * density(29) * density(40) 
  pd(42,29) = pd(42,29) + rrt(347) * density(01) * density(40) 
  pd(42,40) = pd(42,40) + rrt(347) * density(01) * density(29) 
  pd(29,21) = pd(29,21) - rrt(348) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(348) * density(21) * density(40) 
  pd(29,40) = pd(29,40) - rrt(348) * density(21) * density(29) 
  pd(40,21) = pd(40,21) - rrt(348) * density(29) * density(40) 
  pd(40,29) = pd(40,29) - rrt(348) * density(21) * density(40) 
  pd(40,40) = pd(40,40) - rrt(348) * density(21) * density(29) 
  pd(42,21) = pd(42,21) + rrt(348) * density(29) * density(40) 
  pd(42,29) = pd(42,29) + rrt(348) * density(21) * density(40) 
  pd(42,40) = pd(42,40) + rrt(348) * density(21) * density(29) 
  pd(29,29) = pd(29,29) - rrt(349) * density(40)**2 
  pd(29,40) = pd(29,40) - rrt(349) * density(29) * density(40) * 2.0d0
  pd(40,29) = pd(40,29) - rrt(349) * density(40)**2 
  pd(40,40) = pd(40,40) - rrt(349) * density(29) * density(40) * 2.0d0
  pd(42,29) = pd(42,29) + rrt(349) * density(40)**2 
  pd(42,40) = pd(42,40) + rrt(349) * density(29) * density(40) * 2.0d0
  pd(29,01) = pd(29,01) - rrt(350) * density(29) * density(42) 
  pd(29,29) = pd(29,29) - rrt(350) * density(01) * density(42) 
  pd(29,42) = pd(29,42) - rrt(350) * density(01) * density(29) 
  pd(42,01) = pd(42,01) - rrt(350) * density(29) * density(42) 
  pd(42,29) = pd(42,29) - rrt(350) * density(01) * density(42) 
  pd(42,42) = pd(42,42) - rrt(350) * density(01) * density(29) 
  pd(43,01) = pd(43,01) + rrt(350) * density(29) * density(42) 
  pd(43,29) = pd(43,29) + rrt(350) * density(01) * density(42) 
  pd(43,42) = pd(43,42) + rrt(350) * density(01) * density(29) 
  pd(29,21) = pd(29,21) - rrt(351) * density(29) * density(42) 
  pd(29,29) = pd(29,29) - rrt(351) * density(21) * density(42) 
  pd(29,42) = pd(29,42) - rrt(351) * density(21) * density(29) 
  pd(42,21) = pd(42,21) - rrt(351) * density(29) * density(42) 
  pd(42,29) = pd(42,29) - rrt(351) * density(21) * density(42) 
  pd(42,42) = pd(42,42) - rrt(351) * density(21) * density(29) 
  pd(43,21) = pd(43,21) + rrt(351) * density(29) * density(42) 
  pd(43,29) = pd(43,29) + rrt(351) * density(21) * density(42) 
  pd(43,42) = pd(43,42) + rrt(351) * density(21) * density(29) 
  pd(29,14) = pd(29,14) - rrt(352) * density(29) * density(42) 
  pd(29,29) = pd(29,29) - rrt(352) * density(14) * density(42) 
  pd(29,42) = pd(29,42) - rrt(352) * density(14) * density(29) 
  pd(42,14) = pd(42,14) - rrt(352) * density(29) * density(42) 
  pd(42,29) = pd(42,29) - rrt(352) * density(14) * density(42) 
  pd(42,42) = pd(42,42) - rrt(352) * density(14) * density(29) 
  pd(43,14) = pd(43,14) + rrt(352) * density(29) * density(42) 
  pd(43,29) = pd(43,29) + rrt(352) * density(14) * density(42) 
  pd(43,42) = pd(43,42) + rrt(352) * density(14) * density(29) 
  pd(29,29) = pd(29,29) - rrt(353) * density(29) * density(42) * 2.0d0
  pd(29,42) = pd(29,42) - rrt(353) * density(29)**2 
  pd(42,29) = pd(42,29) - rrt(353) * density(29) * density(42) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(353) * density(29)**2 
  pd(43,29) = pd(43,29) + rrt(353) * density(29) * density(42) * 2.0d0
  pd(43,42) = pd(43,42) + rrt(353) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(354) * density(40) * density(42) 
  pd(29,40) = pd(29,40) - rrt(354) * density(29) * density(42) 
  pd(29,42) = pd(29,42) - rrt(354) * density(29) * density(40) 
  pd(42,29) = pd(42,29) - rrt(354) * density(40) * density(42) 
  pd(42,40) = pd(42,40) - rrt(354) * density(29) * density(42) 
  pd(42,42) = pd(42,42) - rrt(354) * density(29) * density(40) 
  pd(43,29) = pd(43,29) + rrt(354) * density(40) * density(42) 
  pd(43,40) = pd(43,40) + rrt(354) * density(29) * density(42) 
  pd(43,42) = pd(43,42) + rrt(354) * density(29) * density(40) 
  pd(42,42) = pd(42,42) - rrt(355) * density(43) 
  pd(42,43) = pd(42,43) - rrt(355) * density(42) 
  pd(43,42) = pd(43,42) - rrt(355) * density(43) 
  pd(43,43) = pd(43,43) - rrt(355) * density(42) 
  pd(44,42) = pd(44,42) + rrt(355) * density(43) 
  pd(44,43) = pd(44,43) + rrt(355) * density(42) 
  pd(14,17) = pd(14,17) + rrt(356) * density(29) 
  pd(14,29) = pd(14,29) + rrt(356) * density(17) 
  pd(17,17) = pd(17,17) - rrt(356) * density(29) 
  pd(17,29) = pd(17,29) - rrt(356) * density(17) 
  pd(29,17) = pd(29,17) - rrt(356) * density(29) 
  pd(29,29) = pd(29,29) - rrt(356) * density(17) 
  pd(33,17) = pd(33,17) + rrt(356) * density(29) 
  pd(33,29) = pd(33,29) + rrt(356) * density(17) 
  pd(14,17) = pd(14,17) + rrt(357) * density(21) 
  pd(14,21) = pd(14,21) + rrt(357) * density(17) 
  pd(17,17) = pd(17,17) - rrt(357) * density(21) 
  pd(17,21) = pd(17,21) - rrt(357) * density(17) 
  pd(21,17) = pd(21,17) - rrt(357) * density(21) 
  pd(21,21) = pd(21,21) - rrt(357) * density(17) 
  pd(34,17) = pd(34,17) + rrt(357) * density(21) 
  pd(34,21) = pd(34,21) + rrt(357) * density(17) 
  pd(17,17) = pd(17,17) - rrt(358) * density(21) 
  pd(17,21) = pd(17,21) - rrt(358) * density(17) 
  pd(21,17) = pd(21,17) - rrt(358) * density(21) 
  pd(21,21) = pd(21,21) - rrt(358) * density(17) 
  pd(29,17) = pd(29,17) + rrt(358) * density(21) 
  pd(29,21) = pd(29,21) + rrt(358) * density(17) 
  pd(45,17) = pd(45,17) + rrt(358) * density(21) 
  pd(45,21) = pd(45,21) + rrt(358) * density(17) 
  pd(17,17) = pd(17,17) - rrt(359) * density(21) 
  pd(17,21) = pd(17,21) - rrt(359) * density(17) 
  pd(21,17) = pd(21,17) - rrt(359) * density(21) 
  pd(21,21) = pd(21,21) - rrt(359) * density(17) 
  pd(33,17) = pd(33,17) + rrt(359) * density(21) 
  pd(33,21) = pd(33,21) + rrt(359) * density(17) 
  pd(40,17) = pd(40,17) + rrt(359) * density(21) 
  pd(40,21) = pd(40,21) + rrt(359) * density(17) 
  pd(17,17) = pd(17,17) - rrt(360) * density(32) 
  pd(17,32) = pd(17,32) - rrt(360) * density(17) 
  pd(21,17) = pd(21,17) + rrt(360) * density(32) 
  pd(21,32) = pd(21,32) + rrt(360) * density(17) 
  pd(32,17) = pd(32,17) - rrt(360) * density(32) 
  pd(32,32) = pd(32,32) - rrt(360) * density(17) 
  pd(45,17) = pd(45,17) + rrt(360) * density(32) 
  pd(45,32) = pd(45,32) + rrt(360) * density(17) 
  pd(14,17) = pd(14,17) + rrt(361) * density(40) 
  pd(14,40) = pd(14,40) + rrt(361) * density(17) 
  pd(17,17) = pd(17,17) - rrt(361) * density(40) 
  pd(17,40) = pd(17,40) - rrt(361) * density(17) 
  pd(40,17) = pd(40,17) - rrt(361) * density(40) 
  pd(40,40) = pd(40,40) - rrt(361) * density(17) 
  pd(45,17) = pd(45,17) + rrt(361) * density(40) 
  pd(45,40) = pd(45,40) + rrt(361) * density(17) 
  pd(17,17) = pd(17,17) - rrt(362) * density(40) 
  pd(17,40) = pd(17,40) - rrt(362) * density(17) 
  pd(18,17) = pd(18,17) + rrt(362) * density(40) 
  pd(18,40) = pd(18,40) + rrt(362) * density(17) 
  pd(29,17) = pd(29,17) + rrt(362) * density(40) 
  pd(29,40) = pd(29,40) + rrt(362) * density(17) 
  pd(40,17) = pd(40,17) - rrt(362) * density(40) 
  pd(40,40) = pd(40,40) - rrt(362) * density(17) 
  pd(01,17) = pd(01,17) + rrt(363) * density(40) 
  pd(01,40) = pd(01,40) + rrt(363) * density(17) 
  pd(17,17) = pd(17,17) - rrt(363) * density(40) 
  pd(17,40) = pd(17,40) - rrt(363) * density(17) 
  pd(33,17) = pd(33,17) + rrt(363) * density(40) 
  pd(33,40) = pd(33,40) + rrt(363) * density(17) 
  pd(40,17) = pd(40,17) - rrt(363) * density(40) 
  pd(40,40) = pd(40,40) - rrt(363) * density(17) 
  pd(01,17) = pd(01,17) + rrt(364) * density(41) 
  pd(01,41) = pd(01,41) + rrt(364) * density(17) 
  pd(17,17) = pd(17,17) - rrt(364) * density(41) 
  pd(17,41) = pd(17,41) - rrt(364) * density(17) 
  pd(41,17) = pd(41,17) - rrt(364) * density(41) 
  pd(41,41) = pd(41,41) - rrt(364) * density(17) 
  pd(45,17) = pd(45,17) + rrt(364) * density(41) 
  pd(45,41) = pd(45,41) + rrt(364) * density(17) 
  pd(01,01) = pd(01,01) - rrt(365) * density(33) 
  pd(01,33) = pd(01,33) - rrt(365) * density(01) 
  pd(14,01) = pd(14,01) + rrt(365) * density(33) 
  pd(14,33) = pd(14,33) + rrt(365) * density(01) 
  pd(33,01) = pd(33,01) - rrt(365) * density(33) 
  pd(33,33) = pd(33,33) - rrt(365) * density(01) 
  pd(45,01) = pd(45,01) + rrt(365) * density(33) 
  pd(45,33) = pd(45,33) + rrt(365) * density(01) 
  pd(21,21) = pd(21,21) - rrt(366) * density(33) 
  pd(21,33) = pd(21,33) - rrt(366) * density(21) 
  pd(29,21) = pd(29,21) + rrt(366) * density(33) 
  pd(29,33) = pd(29,33) + rrt(366) * density(21) 
  pd(33,21) = pd(33,21) - rrt(366) * density(33) 
  pd(33,33) = pd(33,33) - rrt(366) * density(21) 
  pd(34,21) = pd(34,21) + rrt(366) * density(33) 
  pd(34,33) = pd(34,33) + rrt(366) * density(21) 
  pd(21,32) = pd(21,32) + rrt(367) * density(33) 
  pd(21,33) = pd(21,33) + rrt(367) * density(32) 
  pd(32,32) = pd(32,32) - rrt(367) * density(33) 
  pd(32,33) = pd(32,33) - rrt(367) * density(32) 
  pd(33,32) = pd(33,32) - rrt(367) * density(33) 
  pd(33,33) = pd(33,33) - rrt(367) * density(32) 
  pd(34,32) = pd(34,32) + rrt(367) * density(33) 
  pd(34,33) = pd(34,33) + rrt(367) * density(32) 
  pd(29,33) = pd(29,33) + rrt(368) * density(40) 
  pd(29,40) = pd(29,40) + rrt(368) * density(33) 
  pd(33,33) = pd(33,33) - rrt(368) * density(40) 
  pd(33,40) = pd(33,40) - rrt(368) * density(33) 
  pd(40,33) = pd(40,33) - rrt(368) * density(40) 
  pd(40,40) = pd(40,40) - rrt(368) * density(33) 
  pd(45,33) = pd(45,33) + rrt(368) * density(40) 
  pd(45,40) = pd(45,40) + rrt(368) * density(33) 
  pd(14,33) = pd(14,33) + rrt(369) * density(40) 
  pd(14,40) = pd(14,40) + rrt(369) * density(33) 
  pd(33,33) = pd(33,33) - rrt(369) * density(40) 
  pd(33,40) = pd(33,40) - rrt(369) * density(33) 
  pd(34,33) = pd(34,33) + rrt(369) * density(40) 
  pd(34,40) = pd(34,40) + rrt(369) * density(33) 
  pd(40,33) = pd(40,33) - rrt(369) * density(40) 
  pd(40,40) = pd(40,40) - rrt(369) * density(33) 
  pd(15,15) = pd(15,15) - rrt(370) * density(33) 
  pd(15,33) = pd(15,33) - rrt(370) * density(15) 
  pd(17,15) = pd(17,15) + rrt(370) * density(33) 
  pd(17,33) = pd(17,33) + rrt(370) * density(15) 
  pd(29,15) = pd(29,15) + rrt(370) * density(33) 
  pd(29,33) = pd(29,33) + rrt(370) * density(15) 
  pd(33,15) = pd(33,15) - rrt(370) * density(33) 
  pd(33,33) = pd(33,33) - rrt(370) * density(15) 
  pd(33,33) = pd(33,33) - rrt(371) * density(41) 
  pd(33,41) = pd(33,41) - rrt(371) * density(33) 
  pd(40,33) = pd(40,33) + rrt(371) * density(41) 
  pd(40,41) = pd(40,41) + rrt(371) * density(33) 
  pd(41,33) = pd(41,33) - rrt(371) * density(41) 
  pd(41,41) = pd(41,41) - rrt(371) * density(33) 
  pd(45,33) = pd(45,33) + rrt(371) * density(41) 
  pd(45,41) = pd(45,41) + rrt(371) * density(33) 
  pd(29,33) = pd(29,33) + rrt(372) * density(41) 
  pd(29,41) = pd(29,41) + rrt(372) * density(33) 
  pd(33,33) = pd(33,33) - rrt(372) * density(41) 
  pd(33,41) = pd(33,41) - rrt(372) * density(33) 
  pd(41,33) = pd(41,33) - rrt(372) * density(41) 
  pd(41,41) = pd(41,41) - rrt(372) * density(33) 
  pd(46,33) = pd(46,33) + rrt(372) * density(41) 
  pd(46,41) = pd(46,41) + rrt(372) * density(33) 
  pd(01,33) = pd(01,33) + rrt(373) * density(41) 
  pd(01,41) = pd(01,41) + rrt(373) * density(33) 
  pd(33,33) = pd(33,33) - rrt(373) * density(41) 
  pd(33,41) = pd(33,41) - rrt(373) * density(33) 
  pd(34,33) = pd(34,33) + rrt(373) * density(41) 
  pd(34,41) = pd(34,41) + rrt(373) * density(33) 
  pd(41,33) = pd(41,33) - rrt(373) * density(41) 
  pd(41,41) = pd(41,41) - rrt(373) * density(33) 
  pd(29,33) = pd(29,33) + rrt(374) * density(42) 
  pd(29,42) = pd(29,42) + rrt(374) * density(33) 
  pd(33,33) = pd(33,33) - rrt(374) * density(42) 
  pd(33,42) = pd(33,42) - rrt(374) * density(33) 
  pd(42,33) = pd(42,33) - rrt(374) * density(42) 
  pd(42,42) = pd(42,42) - rrt(374) * density(33) 
  pd(47,33) = pd(47,33) + rrt(374) * density(42) 
  pd(47,42) = pd(47,42) + rrt(374) * density(33) 
  pd(01,18) = pd(01,18) + rrt(375) * density(21) 
  pd(01,21) = pd(01,21) + rrt(375) * density(18) 
  pd(18,18) = pd(18,18) - rrt(375) * density(21) 
  pd(18,21) = pd(18,21) - rrt(375) * density(18) 
  pd(21,18) = pd(21,18) - rrt(375) * density(21) 
  pd(21,21) = pd(21,21) - rrt(375) * density(18) 
  pd(34,18) = pd(34,18) + rrt(375) * density(21) 
  pd(34,21) = pd(34,21) + rrt(375) * density(18) 
  pd(14,18) = pd(14,18) + rrt(376) * density(29) 
  pd(14,29) = pd(14,29) + rrt(376) * density(18) 
  pd(18,18) = pd(18,18) - rrt(376) * density(29) 
  pd(18,29) = pd(18,29) - rrt(376) * density(18) 
  pd(29,18) = pd(29,18) - rrt(376) * density(29) 
  pd(29,29) = pd(29,29) - rrt(376) * density(18) 
  pd(45,18) = pd(45,18) + rrt(376) * density(29) 
  pd(45,29) = pd(45,29) + rrt(376) * density(18) 
  pd(01,18) = pd(01,18) + rrt(377) * density(32) 
  pd(01,32) = pd(01,32) + rrt(377) * density(18) 
  pd(18,18) = pd(18,18) - rrt(377) * density(32) 
  pd(18,32) = pd(18,32) - rrt(377) * density(18) 
  pd(29,18) = pd(29,18) + rrt(377) * density(32) 
  pd(29,32) = pd(29,32) + rrt(377) * density(18) 
  pd(32,18) = pd(32,18) - rrt(377) * density(32) 
  pd(32,32) = pd(32,32) - rrt(377) * density(18) 
  pd(34,18) = pd(34,18) + rrt(377) * density(32) 
  pd(34,32) = pd(34,32) + rrt(377) * density(18) 
  pd(01,14) = pd(01,14) + rrt(378) * density(18) 
  pd(01,18) = pd(01,18) + rrt(378) * density(14) 
  pd(14,14) = pd(14,14) - rrt(378) * density(18) 
  pd(14,18) = pd(14,18) - rrt(378) * density(14) 
  pd(17,14) = pd(17,14) + rrt(378) * density(18) 
  pd(17,18) = pd(17,18) + rrt(378) * density(14) 
  pd(18,14) = pd(18,14) - rrt(378) * density(18) 
  pd(18,18) = pd(18,18) - rrt(378) * density(14) 
  pd(01,18) = pd(01,18) + rrt(379) * density(40) 
  pd(01,40) = pd(01,40) + rrt(379) * density(18) 
  pd(18,18) = pd(18,18) - rrt(379) * density(40) 
  pd(18,40) = pd(18,40) - rrt(379) * density(18) 
  pd(40,18) = pd(40,18) - rrt(379) * density(40) 
  pd(40,40) = pd(40,40) - rrt(379) * density(18) 
  pd(45,18) = pd(45,18) + rrt(379) * density(40) 
  pd(45,40) = pd(45,40) + rrt(379) * density(18) 
  pd(01,18) = pd(01,18) + rrt(380) * density(41) 
  pd(01,41) = pd(01,41) + rrt(380) * density(18) 
  pd(18,18) = pd(18,18) - rrt(380) * density(41) 
  pd(18,41) = pd(18,41) - rrt(380) * density(18) 
  pd(41,18) = pd(41,18) - rrt(380) * density(41) 
  pd(41,41) = pd(41,41) - rrt(380) * density(18) 
  pd(46,18) = pd(46,18) + rrt(380) * density(41) 
  pd(46,41) = pd(46,41) + rrt(380) * density(18) 
  pd(01,18) = pd(01,18) + rrt(381) * density(41) 
  pd(01,41) = pd(01,41) + rrt(381) * density(18) 
  pd(14,18) = pd(14,18) + rrt(381) * density(41) 
  pd(14,41) = pd(14,41) + rrt(381) * density(18) 
  pd(18,18) = pd(18,18) - rrt(381) * density(41) 
  pd(18,41) = pd(18,41) - rrt(381) * density(18) 
  pd(41,18) = pd(41,18) - rrt(381) * density(41) 
  pd(41,41) = pd(41,41) - rrt(381) * density(18) 
  pd(45,18) = pd(45,18) + rrt(381) * density(41) 
  pd(45,41) = pd(45,41) + rrt(381) * density(18) 
  pd(01,01) = pd(01,01) - rrt(382) * density(34) 
  pd(01,34) = pd(01,34) - rrt(382) * density(01) 
  pd(34,01) = pd(34,01) - rrt(382) * density(34) 
  pd(34,34) = pd(34,34) - rrt(382) * density(01) 
  pd(40,01) = pd(40,01) + rrt(382) * density(34) 
  pd(40,34) = pd(40,34) + rrt(382) * density(01) 
  pd(45,01) = pd(45,01) + rrt(382) * density(34) 
  pd(45,34) = pd(45,34) + rrt(382) * density(01) 
  pd(14,14) = pd(14,14) - rrt(383) * density(34) 
  pd(14,34) = pd(14,34) - rrt(383) * density(14) 
  pd(29,14) = pd(29,14) + rrt(383) * density(34) 
  pd(29,34) = pd(29,34) + rrt(383) * density(14) 
  pd(34,14) = pd(34,14) - rrt(383) * density(34) 
  pd(34,34) = pd(34,34) - rrt(383) * density(14) 
  pd(45,14) = pd(45,14) + rrt(383) * density(34) 
  pd(45,34) = pd(45,34) + rrt(383) * density(14) 
  pd(21,34) = pd(21,34) + rrt(384) * density(40) 
  pd(21,40) = pd(21,40) + rrt(384) * density(34) 
  pd(34,34) = pd(34,34) - rrt(384) * density(40) 
  pd(34,40) = pd(34,40) - rrt(384) * density(34) 
  pd(40,34) = pd(40,34) - rrt(384) * density(40) 
  pd(40,40) = pd(40,40) - rrt(384) * density(34) 
  pd(45,34) = pd(45,34) + rrt(384) * density(40) 
  pd(45,40) = pd(45,40) + rrt(384) * density(34) 
  pd(32,34) = pd(32,34) + rrt(385) * density(42) 
  pd(32,42) = pd(32,42) + rrt(385) * density(34) 
  pd(34,34) = pd(34,34) - rrt(385) * density(42) 
  pd(34,42) = pd(34,42) - rrt(385) * density(34) 
  pd(42,34) = pd(42,34) - rrt(385) * density(42) 
  pd(42,42) = pd(42,42) - rrt(385) * density(34) 
  pd(45,34) = pd(45,34) + rrt(385) * density(42) 
  pd(45,42) = pd(45,42) + rrt(385) * density(34) 
  pd(21,34) = pd(21,34) + rrt(386) * density(42) 
  pd(21,42) = pd(21,42) + rrt(386) * density(34) 
  pd(34,34) = pd(34,34) - rrt(386) * density(42) 
  pd(34,42) = pd(34,42) - rrt(386) * density(34) 
  pd(42,34) = pd(42,34) - rrt(386) * density(42) 
  pd(42,42) = pd(42,42) - rrt(386) * density(34) 
  pd(47,34) = pd(47,34) + rrt(386) * density(42) 
  pd(47,42) = pd(47,42) + rrt(386) * density(34) 
  pd(01,19) = pd(01,19) + rrt(387) * density(21) 
  pd(01,21) = pd(01,21) + rrt(387) * density(19) 
  pd(14,19) = pd(14,19) + rrt(387) * density(21) 
  pd(14,21) = pd(14,21) + rrt(387) * density(19) 
  pd(19,19) = pd(19,19) - rrt(387) * density(21) 
  pd(19,21) = pd(19,21) - rrt(387) * density(19) 
  pd(21,19) = pd(21,19) - rrt(387) * density(21) 
  pd(21,21) = pd(21,21) - rrt(387) * density(19) 
  pd(34,19) = pd(34,19) + rrt(387) * density(21) 
  pd(34,21) = pd(34,21) + rrt(387) * density(19) 
  pd(01,19) = pd(01,19) + rrt(388) * density(21) 
  pd(01,21) = pd(01,21) + rrt(388) * density(19) 
  pd(19,19) = pd(19,19) - rrt(388) * density(21) 
  pd(19,21) = pd(19,21) - rrt(388) * density(19) 
  pd(21,19) = pd(21,19) - rrt(388) * density(21) 
  pd(21,21) = pd(21,21) - rrt(388) * density(19) 
  pd(47,19) = pd(47,19) + rrt(388) * density(21) 
  pd(47,21) = pd(47,21) + rrt(388) * density(19) 
  pd(01,14) = pd(01,14) + rrt(389) * density(19) 
  pd(01,19) = pd(01,19) + rrt(389) * density(14) 
  pd(14,14) = pd(14,14) - rrt(389) * density(19) 
  pd(14,19) = pd(14,19) - rrt(389) * density(14) 
  pd(18,14) = pd(18,14) + rrt(389) * density(19) 
  pd(18,19) = pd(18,19) + rrt(389) * density(14) 
  pd(19,14) = pd(19,14) - rrt(389) * density(19) 
  pd(19,19) = pd(19,19) - rrt(389) * density(14) 
  pd(01,19) = pd(01,19) + rrt(390) * density(40) 
  pd(01,40) = pd(01,40) + rrt(390) * density(19) 
  pd(14,19) = pd(14,19) + rrt(390) * density(40) 
  pd(14,40) = pd(14,40) + rrt(390) * density(19) 
  pd(19,19) = pd(19,19) - rrt(390) * density(40) 
  pd(19,40) = pd(19,40) - rrt(390) * density(19) 
  pd(40,19) = pd(40,19) - rrt(390) * density(40) 
  pd(40,40) = pd(40,40) - rrt(390) * density(19) 
  pd(45,19) = pd(45,19) + rrt(390) * density(40) 
  pd(45,40) = pd(45,40) + rrt(390) * density(19) 
  pd(01,19) = pd(01,19) + rrt(391) * density(40) 
  pd(01,40) = pd(01,40) + rrt(391) * density(19) 
  pd(19,19) = pd(19,19) - rrt(391) * density(40) 
  pd(19,40) = pd(19,40) - rrt(391) * density(19) 
  pd(40,19) = pd(40,19) - rrt(391) * density(40) 
  pd(40,40) = pd(40,40) - rrt(391) * density(19) 
  pd(46,19) = pd(46,19) + rrt(391) * density(40) 
  pd(46,40) = pd(46,40) + rrt(391) * density(19) 
  pd(40,40) = pd(40,40) - rrt(392) * density(47) 
  pd(40,47) = pd(40,47) - rrt(392) * density(40) 
  pd(42,40) = pd(42,40) + rrt(392) * density(47) 
  pd(42,47) = pd(42,47) + rrt(392) * density(40) 
  pd(45,40) = pd(45,40) + rrt(392) * density(47) 
  pd(45,47) = pd(45,47) + rrt(392) * density(40) 
  pd(47,40) = pd(47,40) - rrt(392) * density(47) 
  pd(47,47) = pd(47,47) - rrt(392) * density(40) 
  pd(40,40) = pd(40,40) - rrt(393) * density(46) 
  pd(40,46) = pd(40,46) - rrt(393) * density(40) 
  pd(41,40) = pd(41,40) + rrt(393) * density(46) 
  pd(41,46) = pd(41,46) + rrt(393) * density(40) 
  pd(45,40) = pd(45,40) + rrt(393) * density(46) 
  pd(45,46) = pd(45,46) + rrt(393) * density(40) 
  pd(46,40) = pd(46,40) - rrt(393) * density(46) 
  pd(46,46) = pd(46,46) - rrt(393) * density(40) 
  pd(01,01) = pd(01,01) + rrt(394) * density(20) 
  pd(01,20) = pd(01,20) + rrt(394) * density(01) 
  pd(18,01) = pd(18,01) + rrt(394) * density(20) 
  pd(18,20) = pd(18,20) + rrt(394) * density(01) 
  pd(20,01) = pd(20,01) - rrt(394) * density(20) 
  pd(20,20) = pd(20,20) - rrt(394) * density(01) 
  pd(01,20) = pd(01,20) + rrt(395) * density(21) * 2.0d0
  pd(01,21) = pd(01,21) + rrt(395) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(395) * density(21) 
  pd(20,21) = pd(20,21) - rrt(395) * density(20) 
  pd(21,20) = pd(21,20) - rrt(395) * density(21) 
  pd(21,21) = pd(21,21) - rrt(395) * density(20) 
  pd(34,20) = pd(34,20) + rrt(395) * density(21) 
  pd(34,21) = pd(34,21) + rrt(395) * density(20) 
  pd(01,20) = pd(01,20) + rrt(396) * density(29) * 2.0d0
  pd(01,29) = pd(01,29) + rrt(396) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(396) * density(29) 
  pd(20,29) = pd(20,29) - rrt(396) * density(20) 
  pd(29,20) = pd(29,20) - rrt(396) * density(29) 
  pd(29,29) = pd(29,29) - rrt(396) * density(20) 
  pd(33,20) = pd(33,20) + rrt(396) * density(29) 
  pd(33,29) = pd(33,29) + rrt(396) * density(20) 
  pd(01,14) = pd(01,14) + rrt(397) * density(20) * 2.0d0
  pd(01,20) = pd(01,20) + rrt(397) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(397) * density(20) 
  pd(14,20) = pd(14,20) - rrt(397) * density(14) 
  pd(17,14) = pd(17,14) + rrt(397) * density(20) 
  pd(17,20) = pd(17,20) + rrt(397) * density(14) 
  pd(20,14) = pd(20,14) - rrt(397) * density(20) 
  pd(20,20) = pd(20,20) - rrt(397) * density(14) 
  pd(01,20) = pd(01,20) + rrt(398) * density(40) * 2.0d0
  pd(01,40) = pd(01,40) + rrt(398) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(398) * density(40) 
  pd(20,40) = pd(20,40) - rrt(398) * density(20) 
  pd(40,20) = pd(40,20) - rrt(398) * density(40) 
  pd(40,40) = pd(40,40) - rrt(398) * density(20) 
  pd(45,20) = pd(45,20) + rrt(398) * density(40) 
  pd(45,40) = pd(45,40) + rrt(398) * density(20) 
  pd(01,01) = pd(01,01) - rrt(399) * density(35) 
  pd(01,35) = pd(01,35) - rrt(399) * density(01) 
  pd(21,01) = pd(21,01) + rrt(399) * density(35) 
  pd(21,35) = pd(21,35) + rrt(399) * density(01) 
  pd(35,01) = pd(35,01) - rrt(399) * density(35) 
  pd(35,35) = pd(35,35) - rrt(399) * density(01) 
  pd(52,01) = pd(52,01) + rrt(399) * density(35) 
  pd(52,35) = pd(52,35) + rrt(399) * density(01) 
  pd(21,21) = pd(21,21) + rrt(400) * density(35) 
  pd(21,35) = pd(21,35) + rrt(400) * density(21) 
  pd(34,21) = pd(34,21) + rrt(400) * density(35) 
  pd(34,35) = pd(34,35) + rrt(400) * density(21) 
  pd(35,21) = pd(35,21) - rrt(400) * density(35) 
  pd(35,35) = pd(35,35) - rrt(400) * density(21) 
  pd(21,26) = pd(21,26) + rrt(401) * density(35) * 2.0d0
  pd(21,35) = pd(21,35) + rrt(401) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(401) * density(35) 
  pd(26,35) = pd(26,35) - rrt(401) * density(26) 
  pd(34,26) = pd(34,26) + rrt(401) * density(35) 
  pd(34,35) = pd(34,35) + rrt(401) * density(26) 
  pd(35,26) = pd(35,26) - rrt(401) * density(35) 
  pd(35,35) = pd(35,35) - rrt(401) * density(26) 
  pd(21,27) = pd(21,27) + rrt(402) * density(35) * 2.0d0
  pd(21,35) = pd(21,35) + rrt(402) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(402) * density(35) 
  pd(27,35) = pd(27,35) - rrt(402) * density(27) 
  pd(34,27) = pd(34,27) + rrt(402) * density(35) 
  pd(34,35) = pd(34,35) + rrt(402) * density(27) 
  pd(35,27) = pd(35,27) - rrt(402) * density(35) 
  pd(35,35) = pd(35,35) - rrt(402) * density(27) 
  pd(29,29) = pd(29,29) - rrt(403) * density(35) 
  pd(29,35) = pd(29,35) - rrt(403) * density(29) 
  pd(32,29) = pd(32,29) + rrt(403) * density(35) 
  pd(32,35) = pd(32,35) + rrt(403) * density(29) 
  pd(34,29) = pd(34,29) + rrt(403) * density(35) 
  pd(34,35) = pd(34,35) + rrt(403) * density(29) 
  pd(35,29) = pd(35,29) - rrt(403) * density(35) 
  pd(35,35) = pd(35,35) - rrt(403) * density(29) 
  pd(21,35) = pd(21,35) + rrt(404) * density(40) * 2.0d0
  pd(21,40) = pd(21,40) + rrt(404) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(404) * density(40) 
  pd(35,40) = pd(35,40) - rrt(404) * density(35) 
  pd(40,35) = pd(40,35) - rrt(404) * density(40) 
  pd(40,40) = pd(40,40) - rrt(404) * density(35) 
  pd(45,35) = pd(45,35) + rrt(404) * density(40) 
  pd(45,40) = pd(45,40) + rrt(404) * density(35) 
  pd(01,01) = pd(01,01) + rrt(405) * density(52) 
  pd(01,52) = pd(01,52) + rrt(405) * density(01) 
  pd(34,01) = pd(34,01) + rrt(405) * density(52) 
  pd(34,52) = pd(34,52) + rrt(405) * density(01) 
  pd(52,01) = pd(52,01) - rrt(405) * density(52) 
  pd(52,52) = pd(52,52) - rrt(405) * density(01) 
  pd(01,21) = pd(01,21) + rrt(406) * density(52) 
  pd(01,52) = pd(01,52) + rrt(406) * density(21) 
  pd(21,21) = pd(21,21) - rrt(406) * density(52) 
  pd(21,52) = pd(21,52) - rrt(406) * density(21) 
  pd(35,21) = pd(35,21) + rrt(406) * density(52) 
  pd(35,52) = pd(35,52) + rrt(406) * density(21) 
  pd(52,21) = pd(52,21) - rrt(406) * density(52) 
  pd(52,52) = pd(52,52) - rrt(406) * density(21) 
  pd(01,01) = pd(01,01) - rrt(407) * density(01) * density(17) * 2.0d0
  pd(01,17) = pd(01,17) - rrt(407) * density(01)**2 
  pd(17,01) = pd(17,01) - rrt(407) * density(01) * density(17) * 2.0d0
  pd(17,17) = pd(17,17) - rrt(407) * density(01)**2 
  pd(19,01) = pd(19,01) + rrt(407) * density(01) * density(17) * 2.0d0
  pd(19,17) = pd(19,17) + rrt(407) * density(01)**2 
  pd(17,17) = pd(17,17) - rrt(408) * density(29) 
  pd(17,29) = pd(17,29) - rrt(408) * density(17) 
  pd(29,17) = pd(29,17) - rrt(408) * density(29) 
  pd(29,29) = pd(29,29) - rrt(408) * density(17) 
  pd(45,17) = pd(45,17) + rrt(408) * density(29) 
  pd(45,29) = pd(45,29) + rrt(408) * density(17) 
  pd(14,14) = pd(14,14) - rrt(409) * density(17) 
  pd(14,17) = pd(14,17) - rrt(409) * density(14) 
  pd(17,14) = pd(17,14) - rrt(409) * density(17) 
  pd(17,17) = pd(17,17) - rrt(409) * density(14) 
  pd(18,14) = pd(18,14) + rrt(409) * density(17) 
  pd(18,17) = pd(18,17) + rrt(409) * density(14) 
  pd(01,01) = pd(01,01) - rrt(410) * density(33) 
  pd(01,33) = pd(01,33) - rrt(410) * density(01) 
  pd(14,01) = pd(14,01) + rrt(410) * density(33) 
  pd(14,33) = pd(14,33) + rrt(410) * density(01) 
  pd(33,01) = pd(33,01) - rrt(410) * density(33) 
  pd(33,33) = pd(33,33) - rrt(410) * density(01) 
  pd(45,01) = pd(45,01) + rrt(410) * density(33) 
  pd(45,33) = pd(45,33) + rrt(410) * density(01) 
  pd(29,29) = pd(29,29) - rrt(411) * density(33) 
  pd(29,33) = pd(29,33) - rrt(411) * density(29) 
  pd(33,29) = pd(33,29) - rrt(411) * density(33) 
  pd(33,33) = pd(33,33) - rrt(411) * density(29) 
  pd(34,29) = pd(34,29) + rrt(411) * density(33) 
  pd(34,33) = pd(34,33) + rrt(411) * density(29) 
  pd(14,14) = pd(14,14) - rrt(412) * density(33) 
  pd(14,33) = pd(14,33) - rrt(412) * density(14) 
  pd(33,14) = pd(33,14) - rrt(412) * density(33) 
  pd(33,33) = pd(33,33) - rrt(412) * density(14) 
  pd(45,14) = pd(45,14) + rrt(412) * density(33) 
  pd(45,33) = pd(45,33) + rrt(412) * density(14) 
  pd(01,01) = pd(01,01) - rrt(413) * density(01) * density(18) * 2.0d0
  pd(01,18) = pd(01,18) - rrt(413) * density(01)**2 
  pd(18,01) = pd(18,01) - rrt(413) * density(01) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(413) * density(01)**2 
  pd(20,01) = pd(20,01) + rrt(413) * density(01) * density(18) * 2.0d0
  pd(20,18) = pd(20,18) + rrt(413) * density(01)**2 
  pd(14,01) = pd(14,01) - rrt(414) * density(14) * density(18) 
  pd(14,14) = pd(14,14) - rrt(414) * density(01) * density(18) 
  pd(14,18) = pd(14,18) - rrt(414) * density(01) * density(14) 
  pd(18,01) = pd(18,01) - rrt(414) * density(14) * density(18) 
  pd(18,14) = pd(18,14) - rrt(414) * density(01) * density(18) 
  pd(18,18) = pd(18,18) - rrt(414) * density(01) * density(14) 
  pd(19,01) = pd(19,01) + rrt(414) * density(14) * density(18) 
  pd(19,14) = pd(19,14) + rrt(414) * density(01) * density(18) 
  pd(19,18) = pd(19,18) + rrt(414) * density(01) * density(14) 
  pd(21,21) = pd(21,21) - rrt(415) * density(21) * density(34) * 2.0d0
  pd(21,34) = pd(21,34) - rrt(415) * density(21)**2 
  pd(34,21) = pd(34,21) - rrt(415) * density(21) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(415) * density(21)**2 
  pd(35,21) = pd(35,21) + rrt(415) * density(21) * density(34) * 2.0d0
  pd(35,34) = pd(35,34) + rrt(415) * density(21)**2 
  pd(01,01) = pd(01,01) - rrt(416) * density(01) * density(34) * 2.0d0
  pd(01,34) = pd(01,34) - rrt(416) * density(01)**2 
  pd(34,01) = pd(34,01) - rrt(416) * density(01) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(416) * density(01)**2 
  pd(52,01) = pd(52,01) + rrt(416) * density(01) * density(34) * 2.0d0
  pd(52,34) = pd(52,34) + rrt(416) * density(01)**2 
  pd(26,26) = pd(26,26) - rrt(417) * density(36) 
  pd(26,36) = pd(26,36) - rrt(417) * density(26) 
  pd(29,26) = pd(29,26) + rrt(417) * density(36) 
  pd(29,36) = pd(29,36) + rrt(417) * density(26) 
  pd(36,26) = pd(36,26) - rrt(417) * density(36) 
  pd(36,36) = pd(36,36) - rrt(417) * density(26) 
  pd(37,26) = pd(37,26) + rrt(417) * density(36) 
  pd(37,36) = pd(37,36) + rrt(417) * density(26) 
  pd(29,32) = pd(29,32) + rrt(418) * density(36) 
  pd(29,36) = pd(29,36) + rrt(418) * density(32) 
  pd(32,32) = pd(32,32) - rrt(418) * density(36) 
  pd(32,36) = pd(32,36) - rrt(418) * density(32) 
  pd(36,32) = pd(36,32) - rrt(418) * density(36) 
  pd(36,36) = pd(36,36) - rrt(418) * density(32) 
  pd(38,32) = pd(38,32) + rrt(418) * density(36) 
  pd(38,36) = pd(38,36) + rrt(418) * density(32) 
  pd(29,36) = pd(29,36) + rrt(419) * density(42) 
  pd(29,42) = pd(29,42) + rrt(419) * density(36) 
  pd(36,36) = pd(36,36) - rrt(419) * density(42) 
  pd(36,42) = pd(36,42) - rrt(419) * density(36) 
  pd(42,36) = pd(42,36) - rrt(419) * density(42) 
  pd(42,42) = pd(42,42) - rrt(419) * density(36) 
  pd(50,36) = pd(50,36) + rrt(419) * density(42) 
  pd(50,42) = pd(50,42) + rrt(419) * density(36) 
  pd(36,36) = pd(36,36) - rrt(420) * density(41) 
  pd(36,41) = pd(36,41) - rrt(420) * density(36) 
  pd(40,36) = pd(40,36) + rrt(420) * density(41) 
  pd(40,41) = pd(40,41) + rrt(420) * density(36) 
  pd(41,36) = pd(41,36) - rrt(420) * density(41) 
  pd(41,41) = pd(41,41) - rrt(420) * density(36) 
  pd(48,36) = pd(48,36) + rrt(420) * density(41) 
  pd(48,41) = pd(48,41) + rrt(420) * density(36) 
  pd(29,36) = pd(29,36) + rrt(421) * density(41) 
  pd(29,41) = pd(29,41) + rrt(421) * density(36) 
  pd(36,36) = pd(36,36) - rrt(421) * density(41) 
  pd(36,41) = pd(36,41) - rrt(421) * density(36) 
  pd(41,36) = pd(41,36) - rrt(421) * density(41) 
  pd(41,41) = pd(41,41) - rrt(421) * density(36) 
  pd(49,36) = pd(49,36) + rrt(421) * density(41) 
  pd(49,41) = pd(49,41) + rrt(421) * density(36) 
  pd(21,29) = pd(21,29) + rrt(422) * density(37) 
  pd(21,37) = pd(21,37) + rrt(422) * density(29) 
  pd(29,29) = pd(29,29) - rrt(422) * density(37) 
  pd(29,37) = pd(29,37) - rrt(422) * density(29) 
  pd(36,29) = pd(36,29) + rrt(422) * density(37) 
  pd(36,37) = pd(36,37) + rrt(422) * density(29) 
  pd(37,29) = pd(37,29) - rrt(422) * density(37) 
  pd(37,37) = pd(37,37) - rrt(422) * density(29) 
  pd(21,32) = pd(21,32) + rrt(423) * density(37) 
  pd(21,37) = pd(21,37) + rrt(423) * density(32) 
  pd(32,32) = pd(32,32) - rrt(423) * density(37) 
  pd(32,37) = pd(32,37) - rrt(423) * density(32) 
  pd(37,32) = pd(37,32) - rrt(423) * density(37) 
  pd(37,37) = pd(37,37) - rrt(423) * density(32) 
  pd(38,32) = pd(38,32) + rrt(423) * density(37) 
  pd(38,37) = pd(38,37) + rrt(423) * density(32) 
  pd(21,37) = pd(21,37) + rrt(424) * density(42) 
  pd(21,42) = pd(21,42) + rrt(424) * density(37) 
  pd(37,37) = pd(37,37) - rrt(424) * density(42) 
  pd(37,42) = pd(37,42) - rrt(424) * density(37) 
  pd(42,37) = pd(42,37) - rrt(424) * density(42) 
  pd(42,42) = pd(42,42) - rrt(424) * density(37) 
  pd(50,37) = pd(50,37) + rrt(424) * density(42) 
  pd(50,42) = pd(50,42) + rrt(424) * density(37) 
  pd(21,37) = pd(21,37) + rrt(425) * density(43) 
  pd(21,43) = pd(21,43) + rrt(425) * density(37) 
  pd(37,37) = pd(37,37) - rrt(425) * density(43) 
  pd(37,43) = pd(37,43) - rrt(425) * density(37) 
  pd(43,37) = pd(43,37) - rrt(425) * density(43) 
  pd(43,43) = pd(43,43) - rrt(425) * density(37) 
  pd(51,37) = pd(51,37) + rrt(425) * density(43) 
  pd(51,43) = pd(51,43) + rrt(425) * density(37) 
  pd(21,29) = pd(21,29) + rrt(426) * density(38) 
  pd(21,38) = pd(21,38) + rrt(426) * density(29) 
  pd(29,29) = pd(29,29) - rrt(426) * density(38) 
  pd(29,38) = pd(29,38) - rrt(426) * density(29) 
  pd(37,29) = pd(37,29) + rrt(426) * density(38) 
  pd(37,38) = pd(37,38) + rrt(426) * density(29) 
  pd(38,29) = pd(38,29) - rrt(426) * density(38) 
  pd(38,38) = pd(38,38) - rrt(426) * density(29) 
  pd(29,38) = pd(29,38) + rrt(427) * density(40) 
  pd(29,40) = pd(29,40) + rrt(427) * density(38) 
  pd(38,38) = pd(38,38) - rrt(427) * density(40) 
  pd(38,40) = pd(38,40) - rrt(427) * density(38) 
  pd(40,38) = pd(40,38) - rrt(427) * density(40) 
  pd(40,40) = pd(40,40) - rrt(427) * density(38) 
  pd(51,38) = pd(51,38) + rrt(427) * density(40) 
  pd(51,40) = pd(51,40) + rrt(427) * density(38) 
  pd(21,38) = pd(21,38) + rrt(428) * density(40) 
  pd(21,40) = pd(21,40) + rrt(428) * density(38) 
  pd(38,38) = pd(38,38) - rrt(428) * density(40) 
  pd(38,40) = pd(38,40) - rrt(428) * density(38) 
  pd(40,38) = pd(40,38) - rrt(428) * density(40) 
  pd(40,40) = pd(40,40) - rrt(428) * density(38) 
  pd(50,38) = pd(50,38) + rrt(428) * density(40) 
  pd(50,40) = pd(50,40) + rrt(428) * density(38) 
  pd(32,38) = pd(32,38) + rrt(429) * density(42) 
  pd(32,42) = pd(32,42) + rrt(429) * density(38) 
  pd(38,38) = pd(38,38) - rrt(429) * density(42) 
  pd(38,42) = pd(38,42) - rrt(429) * density(38) 
  pd(42,38) = pd(42,38) - rrt(429) * density(42) 
  pd(42,42) = pd(42,42) - rrt(429) * density(38) 
  pd(50,38) = pd(50,38) + rrt(429) * density(42) 
  pd(50,42) = pd(50,42) + rrt(429) * density(38) 
  pd(21,38) = pd(21,38) + rrt(430) * density(42) 
  pd(21,42) = pd(21,42) + rrt(430) * density(38) 
  pd(38,38) = pd(38,38) - rrt(430) * density(42) 
  pd(38,42) = pd(38,42) - rrt(430) * density(38) 
  pd(42,38) = pd(42,38) - rrt(430) * density(42) 
  pd(42,42) = pd(42,42) - rrt(430) * density(38) 
  pd(51,38) = pd(51,38) + rrt(430) * density(42) 
  pd(51,42) = pd(51,42) + rrt(430) * density(38) 
  pd(32,38) = pd(32,38) + rrt(431) * density(43) 
  pd(32,43) = pd(32,43) + rrt(431) * density(38) 
  pd(38,38) = pd(38,38) - rrt(431) * density(43) 
  pd(38,43) = pd(38,43) - rrt(431) * density(38) 
  pd(43,38) = pd(43,38) - rrt(431) * density(43) 
  pd(43,43) = pd(43,43) - rrt(431) * density(38) 
  pd(51,38) = pd(51,38) + rrt(431) * density(43) 
  pd(51,43) = pd(51,43) + rrt(431) * density(38) 
  pd(21,21) = pd(21,21) - rrt(432) * density(48) 
  pd(21,48) = pd(21,48) - rrt(432) * density(21) 
  pd(37,21) = pd(37,21) + rrt(432) * density(48) 
  pd(37,48) = pd(37,48) + rrt(432) * density(21) 
  pd(40,21) = pd(40,21) + rrt(432) * density(48) 
  pd(40,48) = pd(40,48) + rrt(432) * density(21) 
  pd(48,21) = pd(48,21) - rrt(432) * density(48) 
  pd(48,48) = pd(48,48) - rrt(432) * density(21) 
  pd(40,42) = pd(40,42) + rrt(433) * density(48) 
  pd(40,48) = pd(40,48) + rrt(433) * density(42) 
  pd(42,42) = pd(42,42) - rrt(433) * density(48) 
  pd(42,48) = pd(42,48) - rrt(433) * density(42) 
  pd(48,42) = pd(48,42) - rrt(433) * density(48) 
  pd(48,48) = pd(48,48) - rrt(433) * density(42) 
  pd(50,42) = pd(50,42) + rrt(433) * density(48) 
  pd(50,48) = pd(50,48) + rrt(433) * density(42) 
  pd(01,41) = pd(01,41) + rrt(434) * density(48) 
  pd(01,48) = pd(01,48) + rrt(434) * density(41) 
  pd(41,41) = pd(41,41) - rrt(434) * density(48) 
  pd(41,48) = pd(41,48) - rrt(434) * density(41) 
  pd(48,41) = pd(48,41) - rrt(434) * density(48) 
  pd(48,48) = pd(48,48) - rrt(434) * density(41) 
  pd(50,41) = pd(50,41) + rrt(434) * density(48) 
  pd(50,48) = pd(50,48) + rrt(434) * density(41) 
  pd(21,32) = pd(21,32) + rrt(435) * density(50) 
  pd(21,50) = pd(21,50) + rrt(435) * density(32) 
  pd(32,32) = pd(32,32) - rrt(435) * density(50) 
  pd(32,50) = pd(32,50) - rrt(435) * density(32) 
  pd(50,32) = pd(50,32) - rrt(435) * density(50) 
  pd(50,50) = pd(50,50) - rrt(435) * density(32) 
  pd(51,32) = pd(51,32) + rrt(435) * density(50) 
  pd(51,50) = pd(51,50) + rrt(435) * density(32) 
  pd(40,42) = pd(40,42) + rrt(436) * density(50) 
  pd(40,50) = pd(40,50) + rrt(436) * density(42) 
  pd(42,42) = pd(42,42) - rrt(436) * density(50) 
  pd(42,50) = pd(42,50) - rrt(436) * density(42) 
  pd(50,42) = pd(50,42) - rrt(436) * density(50) 
  pd(50,50) = pd(50,50) - rrt(436) * density(42) 
  pd(51,42) = pd(51,42) + rrt(436) * density(50) 
  pd(51,50) = pd(51,50) + rrt(436) * density(42) 
  pd(42,43) = pd(42,43) + rrt(437) * density(50) 
  pd(42,50) = pd(42,50) + rrt(437) * density(43) 
  pd(43,43) = pd(43,43) - rrt(437) * density(50) 
  pd(43,50) = pd(43,50) - rrt(437) * density(43) 
  pd(50,43) = pd(50,43) - rrt(437) * density(50) 
  pd(50,50) = pd(50,50) - rrt(437) * density(43) 
  pd(51,43) = pd(51,43) + rrt(437) * density(50) 
  pd(51,50) = pd(51,50) + rrt(437) * density(43) 
  pd(42,44) = pd(42,44) + rrt(438) * density(50) * 2.0d0
  pd(42,50) = pd(42,50) + rrt(438) * density(44) * 2.0d0
  pd(44,44) = pd(44,44) - rrt(438) * density(50) 
  pd(44,50) = pd(44,50) - rrt(438) * density(44) 
  pd(50,44) = pd(50,44) - rrt(438) * density(50) 
  pd(50,50) = pd(50,50) - rrt(438) * density(44) 
  pd(51,44) = pd(51,44) + rrt(438) * density(50) 
  pd(51,50) = pd(51,50) + rrt(438) * density(44) 
  pd(40,40) = pd(40,40) - rrt(439) * density(51) 
  pd(40,51) = pd(40,51) - rrt(439) * density(40) 
  pd(42,40) = pd(42,40) + rrt(439) * density(51) 
  pd(42,51) = pd(42,51) + rrt(439) * density(40) 
  pd(50,40) = pd(50,40) + rrt(439) * density(51) 
  pd(50,51) = pd(50,51) + rrt(439) * density(40) 
  pd(51,40) = pd(51,40) - rrt(439) * density(51) 
  pd(51,51) = pd(51,51) - rrt(439) * density(40) 
  pd(21,01) = pd(21,01) + rrt(440) * density(39) 
  pd(21,39) = pd(21,39) + rrt(440) * density(01) 
  pd(37,01) = pd(37,01) + rrt(440) * density(39) 
  pd(37,39) = pd(37,39) + rrt(440) * density(01) 
  pd(39,01) = pd(39,01) - rrt(440) * density(39) 
  pd(39,39) = pd(39,39) - rrt(440) * density(01) 
  pd(21,21) = pd(21,21) + rrt(441) * density(39) 
  pd(21,39) = pd(21,39) + rrt(441) * density(21) 
  pd(37,21) = pd(37,21) + rrt(441) * density(39) 
  pd(37,39) = pd(37,39) + rrt(441) * density(21) 
  pd(39,21) = pd(39,21) - rrt(441) * density(39) 
  pd(39,39) = pd(39,39) - rrt(441) * density(21) 
  pd(21,29) = pd(21,29) + rrt(442) * density(39) 
  pd(21,39) = pd(21,39) + rrt(442) * density(29) 
  pd(29,29) = pd(29,29) - rrt(442) * density(39) 
  pd(29,39) = pd(29,39) - rrt(442) * density(29) 
  pd(38,29) = pd(38,29) + rrt(442) * density(39) 
  pd(38,39) = pd(38,39) + rrt(442) * density(29) 
  pd(39,29) = pd(39,29) - rrt(442) * density(39) 
  pd(39,39) = pd(39,39) - rrt(442) * density(29) 
  pd(21,29) = pd(21,29) + rrt(443) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(443) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(443) * density(39) 
  pd(29,39) = pd(29,39) - rrt(443) * density(29) 
  pd(36,29) = pd(36,29) + rrt(443) * density(39) 
  pd(36,39) = pd(36,39) + rrt(443) * density(29) 
  pd(39,29) = pd(39,29) - rrt(443) * density(39) 
  pd(39,39) = pd(39,39) - rrt(443) * density(29) 
  pd(21,26) = pd(21,26) + rrt(444) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(444) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(444) * density(39) 
  pd(26,39) = pd(26,39) - rrt(444) * density(26) 
  pd(37,26) = pd(37,26) + rrt(444) * density(39) 
  pd(37,39) = pd(37,39) + rrt(444) * density(26) 
  pd(39,26) = pd(39,26) - rrt(444) * density(39) 
  pd(39,39) = pd(39,39) - rrt(444) * density(26) 
  pd(21,27) = pd(21,27) + rrt(445) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(445) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(445) * density(39) 
  pd(27,39) = pd(27,39) - rrt(445) * density(27) 
  pd(37,27) = pd(37,27) + rrt(445) * density(39) 
  pd(37,39) = pd(37,39) + rrt(445) * density(27) 
  pd(39,27) = pd(39,27) - rrt(445) * density(39) 
  pd(39,39) = pd(39,39) - rrt(445) * density(27) 
  pd(21,39) = pd(21,39) + rrt(446) * density(40) 
  pd(21,40) = pd(21,40) + rrt(446) * density(39) 
  pd(39,39) = pd(39,39) - rrt(446) * density(40) 
  pd(39,40) = pd(39,40) - rrt(446) * density(39) 
  pd(40,39) = pd(40,39) - rrt(446) * density(40) 
  pd(40,40) = pd(40,40) - rrt(446) * density(39) 
  pd(51,39) = pd(51,39) + rrt(446) * density(40) 
  pd(51,40) = pd(51,40) + rrt(446) * density(39) 
  pd(21,21) = pd(21,21) - rrt(447) * density(36) 
  pd(21,36) = pd(21,36) - rrt(447) * density(21) 
  pd(36,21) = pd(36,21) - rrt(447) * density(36) 
  pd(36,36) = pd(36,36) - rrt(447) * density(21) 
  pd(38,21) = pd(38,21) + rrt(447) * density(36) 
  pd(38,36) = pd(38,36) + rrt(447) * density(21) 
  pd(36,36) = pd(36,36) - rrt(448) * density(40) 
  pd(36,40) = pd(36,40) - rrt(448) * density(36) 
  pd(40,36) = pd(40,36) - rrt(448) * density(40) 
  pd(40,40) = pd(40,40) - rrt(448) * density(36) 
  pd(50,36) = pd(50,36) + rrt(448) * density(40) 
  pd(50,40) = pd(50,40) + rrt(448) * density(36) 
  pd(21,21) = pd(21,21) - rrt(449) * density(37) 
  pd(21,37) = pd(21,37) - rrt(449) * density(21) 
  pd(37,21) = pd(37,21) - rrt(449) * density(37) 
  pd(37,37) = pd(37,37) - rrt(449) * density(21) 
  pd(39,21) = pd(39,21) + rrt(449) * density(37) 
  pd(39,37) = pd(39,37) + rrt(449) * density(21) 
  pd(14,17) = pd(14,17) + rrt(450) * density(36) 
  pd(14,36) = pd(14,36) + rrt(450) * density(17) 
  pd(17,17) = pd(17,17) - rrt(450) * density(36) 
  pd(17,36) = pd(17,36) - rrt(450) * density(17) 
  pd(29,17) = pd(29,17) + rrt(450) * density(36) 
  pd(29,36) = pd(29,36) + rrt(450) * density(17) 
  pd(36,17) = pd(36,17) - rrt(450) * density(36) 
  pd(36,36) = pd(36,36) - rrt(450) * density(17) 
  pd(01,18) = pd(01,18) + rrt(451) * density(36) 
  pd(01,36) = pd(01,36) + rrt(451) * density(18) 
  pd(18,18) = pd(18,18) - rrt(451) * density(36) 
  pd(18,36) = pd(18,36) - rrt(451) * density(18) 
  pd(29,18) = pd(29,18) + rrt(451) * density(36) 
  pd(29,36) = pd(29,36) + rrt(451) * density(18) 
  pd(36,18) = pd(36,18) - rrt(451) * density(36) 
  pd(36,36) = pd(36,36) - rrt(451) * density(18) 
  pd(29,33) = pd(29,33) + rrt(452) * density(36) * 2.0d0
  pd(29,36) = pd(29,36) + rrt(452) * density(33) * 2.0d0
  pd(33,33) = pd(33,33) - rrt(452) * density(36) 
  pd(33,36) = pd(33,36) - rrt(452) * density(33) 
  pd(36,33) = pd(36,33) - rrt(452) * density(36) 
  pd(36,36) = pd(36,36) - rrt(452) * density(33) 
  pd(21,34) = pd(21,34) + rrt(453) * density(36) 
  pd(21,36) = pd(21,36) + rrt(453) * density(34) 
  pd(29,34) = pd(29,34) + rrt(453) * density(36) 
  pd(29,36) = pd(29,36) + rrt(453) * density(34) 
  pd(34,34) = pd(34,34) - rrt(453) * density(36) 
  pd(34,36) = pd(34,36) - rrt(453) * density(34) 
  pd(36,34) = pd(36,34) - rrt(453) * density(36) 
  pd(36,36) = pd(36,36) - rrt(453) * density(34) 
  pd(29,36) = pd(29,36) + rrt(454) * density(45) 
  pd(29,45) = pd(29,45) + rrt(454) * density(36) 
  pd(36,36) = pd(36,36) - rrt(454) * density(45) 
  pd(36,45) = pd(36,45) - rrt(454) * density(36) 
  pd(40,36) = pd(40,36) + rrt(454) * density(45) 
  pd(40,45) = pd(40,45) + rrt(454) * density(36) 
  pd(45,36) = pd(45,36) - rrt(454) * density(45) 
  pd(45,45) = pd(45,45) - rrt(454) * density(36) 
  pd(29,36) = pd(29,36) + rrt(455) * density(46) 
  pd(29,46) = pd(29,46) + rrt(455) * density(36) 
  pd(36,36) = pd(36,36) - rrt(455) * density(46) 
  pd(36,46) = pd(36,46) - rrt(455) * density(36) 
  pd(41,36) = pd(41,36) + rrt(455) * density(46) 
  pd(41,46) = pd(41,46) + rrt(455) * density(36) 
  pd(46,36) = pd(46,36) - rrt(455) * density(46) 
  pd(46,46) = pd(46,46) - rrt(455) * density(36) 
  pd(29,36) = pd(29,36) + rrt(456) * density(47) 
  pd(29,47) = pd(29,47) + rrt(456) * density(36) 
  pd(36,36) = pd(36,36) - rrt(456) * density(47) 
  pd(36,47) = pd(36,47) - rrt(456) * density(36) 
  pd(42,36) = pd(42,36) + rrt(456) * density(47) 
  pd(42,47) = pd(42,47) + rrt(456) * density(36) 
  pd(47,36) = pd(47,36) - rrt(456) * density(47) 
  pd(47,47) = pd(47,47) - rrt(456) * density(36) 
  pd(14,17) = pd(14,17) + rrt(457) * density(37) 
  pd(14,37) = pd(14,37) + rrt(457) * density(17) 
  pd(17,17) = pd(17,17) - rrt(457) * density(37) 
  pd(17,37) = pd(17,37) - rrt(457) * density(17) 
  pd(21,17) = pd(21,17) + rrt(457) * density(37) 
  pd(21,37) = pd(21,37) + rrt(457) * density(17) 
  pd(37,17) = pd(37,17) - rrt(457) * density(37) 
  pd(37,37) = pd(37,37) - rrt(457) * density(17) 
  pd(01,18) = pd(01,18) + rrt(458) * density(37) 
  pd(01,37) = pd(01,37) + rrt(458) * density(18) 
  pd(18,18) = pd(18,18) - rrt(458) * density(37) 
  pd(18,37) = pd(18,37) - rrt(458) * density(18) 
  pd(21,18) = pd(21,18) + rrt(458) * density(37) 
  pd(21,37) = pd(21,37) + rrt(458) * density(18) 
  pd(37,18) = pd(37,18) - rrt(458) * density(37) 
  pd(37,37) = pd(37,37) - rrt(458) * density(18) 
  pd(21,33) = pd(21,33) + rrt(459) * density(37) 
  pd(21,37) = pd(21,37) + rrt(459) * density(33) 
  pd(29,33) = pd(29,33) + rrt(459) * density(37) 
  pd(29,37) = pd(29,37) + rrt(459) * density(33) 
  pd(33,33) = pd(33,33) - rrt(459) * density(37) 
  pd(33,37) = pd(33,37) - rrt(459) * density(33) 
  pd(37,33) = pd(37,33) - rrt(459) * density(37) 
  pd(37,37) = pd(37,37) - rrt(459) * density(33) 
  pd(21,34) = pd(21,34) + rrt(460) * density(37) * 2.0d0
  pd(21,37) = pd(21,37) + rrt(460) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(460) * density(37) 
  pd(34,37) = pd(34,37) - rrt(460) * density(34) 
  pd(37,34) = pd(37,34) - rrt(460) * density(37) 
  pd(37,37) = pd(37,37) - rrt(460) * density(34) 
  pd(21,37) = pd(21,37) + rrt(461) * density(45) 
  pd(21,45) = pd(21,45) + rrt(461) * density(37) 
  pd(37,37) = pd(37,37) - rrt(461) * density(45) 
  pd(37,45) = pd(37,45) - rrt(461) * density(37) 
  pd(40,37) = pd(40,37) + rrt(461) * density(45) 
  pd(40,45) = pd(40,45) + rrt(461) * density(37) 
  pd(45,37) = pd(45,37) - rrt(461) * density(45) 
  pd(45,45) = pd(45,45) - rrt(461) * density(37) 
  pd(21,37) = pd(21,37) + rrt(462) * density(46) 
  pd(21,46) = pd(21,46) + rrt(462) * density(37) 
  pd(37,37) = pd(37,37) - rrt(462) * density(46) 
  pd(37,46) = pd(37,46) - rrt(462) * density(37) 
  pd(41,37) = pd(41,37) + rrt(462) * density(46) 
  pd(41,46) = pd(41,46) + rrt(462) * density(37) 
  pd(46,37) = pd(46,37) - rrt(462) * density(46) 
  pd(46,46) = pd(46,46) - rrt(462) * density(37) 
  pd(21,37) = pd(21,37) + rrt(463) * density(47) 
  pd(21,47) = pd(21,47) + rrt(463) * density(37) 
  pd(37,37) = pd(37,37) - rrt(463) * density(47) 
  pd(37,47) = pd(37,47) - rrt(463) * density(37) 
  pd(42,37) = pd(42,37) + rrt(463) * density(47) 
  pd(42,47) = pd(42,47) + rrt(463) * density(37) 
  pd(47,37) = pd(47,37) - rrt(463) * density(47) 
  pd(47,47) = pd(47,47) - rrt(463) * density(37) 
  pd(14,17) = pd(14,17) + rrt(464) * density(38) 
  pd(14,38) = pd(14,38) + rrt(464) * density(17) 
  pd(17,17) = pd(17,17) - rrt(464) * density(38) 
  pd(17,38) = pd(17,38) - rrt(464) * density(17) 
  pd(32,17) = pd(32,17) + rrt(464) * density(38) 
  pd(32,38) = pd(32,38) + rrt(464) * density(17) 
  pd(38,17) = pd(38,17) - rrt(464) * density(38) 
  pd(38,38) = pd(38,38) - rrt(464) * density(17) 
  pd(01,18) = pd(01,18) + rrt(465) * density(38) 
  pd(01,38) = pd(01,38) + rrt(465) * density(18) 
  pd(18,18) = pd(18,18) - rrt(465) * density(38) 
  pd(18,38) = pd(18,38) - rrt(465) * density(18) 
  pd(32,18) = pd(32,18) + rrt(465) * density(38) 
  pd(32,38) = pd(32,38) + rrt(465) * density(18) 
  pd(38,18) = pd(38,18) - rrt(465) * density(38) 
  pd(38,38) = pd(38,38) - rrt(465) * density(18) 
  pd(29,33) = pd(29,33) + rrt(466) * density(38) 
  pd(29,38) = pd(29,38) + rrt(466) * density(33) 
  pd(32,33) = pd(32,33) + rrt(466) * density(38) 
  pd(32,38) = pd(32,38) + rrt(466) * density(33) 
  pd(33,33) = pd(33,33) - rrt(466) * density(38) 
  pd(33,38) = pd(33,38) - rrt(466) * density(33) 
  pd(38,33) = pd(38,33) - rrt(466) * density(38) 
  pd(38,38) = pd(38,38) - rrt(466) * density(33) 
  pd(21,34) = pd(21,34) + rrt(467) * density(38) 
  pd(21,38) = pd(21,38) + rrt(467) * density(34) 
  pd(32,34) = pd(32,34) + rrt(467) * density(38) 
  pd(32,38) = pd(32,38) + rrt(467) * density(34) 
  pd(34,34) = pd(34,34) - rrt(467) * density(38) 
  pd(34,38) = pd(34,38) - rrt(467) * density(34) 
  pd(38,34) = pd(38,34) - rrt(467) * density(38) 
  pd(38,38) = pd(38,38) - rrt(467) * density(34) 
  pd(32,38) = pd(32,38) + rrt(468) * density(45) 
  pd(32,45) = pd(32,45) + rrt(468) * density(38) 
  pd(38,38) = pd(38,38) - rrt(468) * density(45) 
  pd(38,45) = pd(38,45) - rrt(468) * density(38) 
  pd(40,38) = pd(40,38) + rrt(468) * density(45) 
  pd(40,45) = pd(40,45) + rrt(468) * density(38) 
  pd(45,38) = pd(45,38) - rrt(468) * density(45) 
  pd(45,45) = pd(45,45) - rrt(468) * density(38) 
  pd(32,38) = pd(32,38) + rrt(469) * density(46) 
  pd(32,46) = pd(32,46) + rrt(469) * density(38) 
  pd(38,38) = pd(38,38) - rrt(469) * density(46) 
  pd(38,46) = pd(38,46) - rrt(469) * density(38) 
  pd(41,38) = pd(41,38) + rrt(469) * density(46) 
  pd(41,46) = pd(41,46) + rrt(469) * density(38) 
  pd(46,38) = pd(46,38) - rrt(469) * density(46) 
  pd(46,46) = pd(46,46) - rrt(469) * density(38) 
  pd(32,38) = pd(32,38) + rrt(470) * density(47) 
  pd(32,47) = pd(32,47) + rrt(470) * density(38) 
  pd(38,38) = pd(38,38) - rrt(470) * density(47) 
  pd(38,47) = pd(38,47) - rrt(470) * density(38) 
  pd(42,38) = pd(42,38) + rrt(470) * density(47) 
  pd(42,47) = pd(42,47) + rrt(470) * density(38) 
  pd(47,38) = pd(47,38) - rrt(470) * density(47) 
  pd(47,47) = pd(47,47) - rrt(470) * density(38) 
  pd(14,17) = pd(14,17) + rrt(471) * density(48) 
  pd(14,48) = pd(14,48) + rrt(471) * density(17) 
  pd(17,17) = pd(17,17) - rrt(471) * density(48) 
  pd(17,48) = pd(17,48) - rrt(471) * density(17) 
  pd(40,17) = pd(40,17) + rrt(471) * density(48) 
  pd(40,48) = pd(40,48) + rrt(471) * density(17) 
  pd(48,17) = pd(48,17) - rrt(471) * density(48) 
  pd(48,48) = pd(48,48) - rrt(471) * density(17) 
  pd(01,18) = pd(01,18) + rrt(472) * density(48) 
  pd(01,48) = pd(01,48) + rrt(472) * density(18) 
  pd(18,18) = pd(18,18) - rrt(472) * density(48) 
  pd(18,48) = pd(18,48) - rrt(472) * density(18) 
  pd(40,18) = pd(40,18) + rrt(472) * density(48) 
  pd(40,48) = pd(40,48) + rrt(472) * density(18) 
  pd(48,18) = pd(48,18) - rrt(472) * density(48) 
  pd(48,48) = pd(48,48) - rrt(472) * density(18) 
  pd(29,33) = pd(29,33) + rrt(473) * density(48) 
  pd(29,48) = pd(29,48) + rrt(473) * density(33) 
  pd(33,33) = pd(33,33) - rrt(473) * density(48) 
  pd(33,48) = pd(33,48) - rrt(473) * density(33) 
  pd(40,33) = pd(40,33) + rrt(473) * density(48) 
  pd(40,48) = pd(40,48) + rrt(473) * density(33) 
  pd(48,33) = pd(48,33) - rrt(473) * density(48) 
  pd(48,48) = pd(48,48) - rrt(473) * density(33) 
  pd(21,34) = pd(21,34) + rrt(474) * density(48) 
  pd(21,48) = pd(21,48) + rrt(474) * density(34) 
  pd(34,34) = pd(34,34) - rrt(474) * density(48) 
  pd(34,48) = pd(34,48) - rrt(474) * density(34) 
  pd(40,34) = pd(40,34) + rrt(474) * density(48) 
  pd(40,48) = pd(40,48) + rrt(474) * density(34) 
  pd(48,34) = pd(48,34) - rrt(474) * density(48) 
  pd(48,48) = pd(48,48) - rrt(474) * density(34) 
  pd(40,45) = pd(40,45) + rrt(475) * density(48) * 2.0d0
  pd(40,48) = pd(40,48) + rrt(475) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(475) * density(48) 
  pd(45,48) = pd(45,48) - rrt(475) * density(45) 
  pd(48,45) = pd(48,45) - rrt(475) * density(48) 
  pd(48,48) = pd(48,48) - rrt(475) * density(45) 
  pd(40,46) = pd(40,46) + rrt(476) * density(48) 
  pd(40,48) = pd(40,48) + rrt(476) * density(46) 
  pd(41,46) = pd(41,46) + rrt(476) * density(48) 
  pd(41,48) = pd(41,48) + rrt(476) * density(46) 
  pd(46,46) = pd(46,46) - rrt(476) * density(48) 
  pd(46,48) = pd(46,48) - rrt(476) * density(46) 
  pd(48,46) = pd(48,46) - rrt(476) * density(48) 
  pd(48,48) = pd(48,48) - rrt(476) * density(46) 
  pd(40,47) = pd(40,47) + rrt(477) * density(48) 
  pd(40,48) = pd(40,48) + rrt(477) * density(47) 
  pd(42,47) = pd(42,47) + rrt(477) * density(48) 
  pd(42,48) = pd(42,48) + rrt(477) * density(47) 
  pd(47,47) = pd(47,47) - rrt(477) * density(48) 
  pd(47,48) = pd(47,48) - rrt(477) * density(47) 
  pd(48,47) = pd(48,47) - rrt(477) * density(48) 
  pd(48,48) = pd(48,48) - rrt(477) * density(47) 
  pd(14,17) = pd(14,17) + rrt(478) * density(49) 
  pd(14,49) = pd(14,49) + rrt(478) * density(17) 
  pd(17,17) = pd(17,17) - rrt(478) * density(49) 
  pd(17,49) = pd(17,49) - rrt(478) * density(17) 
  pd(41,17) = pd(41,17) + rrt(478) * density(49) 
  pd(41,49) = pd(41,49) + rrt(478) * density(17) 
  pd(49,17) = pd(49,17) - rrt(478) * density(49) 
  pd(49,49) = pd(49,49) - rrt(478) * density(17) 
  pd(01,18) = pd(01,18) + rrt(479) * density(49) 
  pd(01,49) = pd(01,49) + rrt(479) * density(18) 
  pd(18,18) = pd(18,18) - rrt(479) * density(49) 
  pd(18,49) = pd(18,49) - rrt(479) * density(18) 
  pd(41,18) = pd(41,18) + rrt(479) * density(49) 
  pd(41,49) = pd(41,49) + rrt(479) * density(18) 
  pd(49,18) = pd(49,18) - rrt(479) * density(49) 
  pd(49,49) = pd(49,49) - rrt(479) * density(18) 
  pd(29,33) = pd(29,33) + rrt(480) * density(49) 
  pd(29,49) = pd(29,49) + rrt(480) * density(33) 
  pd(33,33) = pd(33,33) - rrt(480) * density(49) 
  pd(33,49) = pd(33,49) - rrt(480) * density(33) 
  pd(41,33) = pd(41,33) + rrt(480) * density(49) 
  pd(41,49) = pd(41,49) + rrt(480) * density(33) 
  pd(49,33) = pd(49,33) - rrt(480) * density(49) 
  pd(49,49) = pd(49,49) - rrt(480) * density(33) 
  pd(21,34) = pd(21,34) + rrt(481) * density(49) 
  pd(21,49) = pd(21,49) + rrt(481) * density(34) 
  pd(34,34) = pd(34,34) - rrt(481) * density(49) 
  pd(34,49) = pd(34,49) - rrt(481) * density(34) 
  pd(41,34) = pd(41,34) + rrt(481) * density(49) 
  pd(41,49) = pd(41,49) + rrt(481) * density(34) 
  pd(49,34) = pd(49,34) - rrt(481) * density(49) 
  pd(49,49) = pd(49,49) - rrt(481) * density(34) 
  pd(40,45) = pd(40,45) + rrt(482) * density(49) 
  pd(40,49) = pd(40,49) + rrt(482) * density(45) 
  pd(41,45) = pd(41,45) + rrt(482) * density(49) 
  pd(41,49) = pd(41,49) + rrt(482) * density(45) 
  pd(45,45) = pd(45,45) - rrt(482) * density(49) 
  pd(45,49) = pd(45,49) - rrt(482) * density(45) 
  pd(49,45) = pd(49,45) - rrt(482) * density(49) 
  pd(49,49) = pd(49,49) - rrt(482) * density(45) 
  pd(41,46) = pd(41,46) + rrt(483) * density(49) * 2.0d0
  pd(41,49) = pd(41,49) + rrt(483) * density(46) * 2.0d0
  pd(46,46) = pd(46,46) - rrt(483) * density(49) 
  pd(46,49) = pd(46,49) - rrt(483) * density(46) 
  pd(49,46) = pd(49,46) - rrt(483) * density(49) 
  pd(49,49) = pd(49,49) - rrt(483) * density(46) 
  pd(41,47) = pd(41,47) + rrt(484) * density(49) 
  pd(41,49) = pd(41,49) + rrt(484) * density(47) 
  pd(42,47) = pd(42,47) + rrt(484) * density(49) 
  pd(42,49) = pd(42,49) + rrt(484) * density(47) 
  pd(47,47) = pd(47,47) - rrt(484) * density(49) 
  pd(47,49) = pd(47,49) - rrt(484) * density(47) 
  pd(49,47) = pd(49,47) - rrt(484) * density(49) 
  pd(49,49) = pd(49,49) - rrt(484) * density(47) 
  pd(14,17) = pd(14,17) + rrt(485) * density(50) 
  pd(14,50) = pd(14,50) + rrt(485) * density(17) 
  pd(17,17) = pd(17,17) - rrt(485) * density(50) 
  pd(17,50) = pd(17,50) - rrt(485) * density(17) 
  pd(42,17) = pd(42,17) + rrt(485) * density(50) 
  pd(42,50) = pd(42,50) + rrt(485) * density(17) 
  pd(50,17) = pd(50,17) - rrt(485) * density(50) 
  pd(50,50) = pd(50,50) - rrt(485) * density(17) 
  pd(01,18) = pd(01,18) + rrt(486) * density(50) 
  pd(01,50) = pd(01,50) + rrt(486) * density(18) 
  pd(18,18) = pd(18,18) - rrt(486) * density(50) 
  pd(18,50) = pd(18,50) - rrt(486) * density(18) 
  pd(42,18) = pd(42,18) + rrt(486) * density(50) 
  pd(42,50) = pd(42,50) + rrt(486) * density(18) 
  pd(50,18) = pd(50,18) - rrt(486) * density(50) 
  pd(50,50) = pd(50,50) - rrt(486) * density(18) 
  pd(29,33) = pd(29,33) + rrt(487) * density(50) 
  pd(29,50) = pd(29,50) + rrt(487) * density(33) 
  pd(33,33) = pd(33,33) - rrt(487) * density(50) 
  pd(33,50) = pd(33,50) - rrt(487) * density(33) 
  pd(42,33) = pd(42,33) + rrt(487) * density(50) 
  pd(42,50) = pd(42,50) + rrt(487) * density(33) 
  pd(50,33) = pd(50,33) - rrt(487) * density(50) 
  pd(50,50) = pd(50,50) - rrt(487) * density(33) 
  pd(21,34) = pd(21,34) + rrt(488) * density(50) 
  pd(21,50) = pd(21,50) + rrt(488) * density(34) 
  pd(34,34) = pd(34,34) - rrt(488) * density(50) 
  pd(34,50) = pd(34,50) - rrt(488) * density(34) 
  pd(42,34) = pd(42,34) + rrt(488) * density(50) 
  pd(42,50) = pd(42,50) + rrt(488) * density(34) 
  pd(50,34) = pd(50,34) - rrt(488) * density(50) 
  pd(50,50) = pd(50,50) - rrt(488) * density(34) 
  pd(40,45) = pd(40,45) + rrt(489) * density(50) 
  pd(40,50) = pd(40,50) + rrt(489) * density(45) 
  pd(42,45) = pd(42,45) + rrt(489) * density(50) 
  pd(42,50) = pd(42,50) + rrt(489) * density(45) 
  pd(45,45) = pd(45,45) - rrt(489) * density(50) 
  pd(45,50) = pd(45,50) - rrt(489) * density(45) 
  pd(50,45) = pd(50,45) - rrt(489) * density(50) 
  pd(50,50) = pd(50,50) - rrt(489) * density(45) 
  pd(41,46) = pd(41,46) + rrt(490) * density(50) 
  pd(41,50) = pd(41,50) + rrt(490) * density(46) 
  pd(42,46) = pd(42,46) + rrt(490) * density(50) 
  pd(42,50) = pd(42,50) + rrt(490) * density(46) 
  pd(46,46) = pd(46,46) - rrt(490) * density(50) 
  pd(46,50) = pd(46,50) - rrt(490) * density(46) 
  pd(50,46) = pd(50,46) - rrt(490) * density(50) 
  pd(50,50) = pd(50,50) - rrt(490) * density(46) 
  pd(42,47) = pd(42,47) + rrt(491) * density(50) * 2.0d0
  pd(42,50) = pd(42,50) + rrt(491) * density(47) * 2.0d0
  pd(47,47) = pd(47,47) - rrt(491) * density(50) 
  pd(47,50) = pd(47,50) - rrt(491) * density(47) 
  pd(50,47) = pd(50,47) - rrt(491) * density(50) 
  pd(50,50) = pd(50,50) - rrt(491) * density(47) 
  pd(14,17) = pd(14,17) + rrt(492) * density(51) 
  pd(14,51) = pd(14,51) + rrt(492) * density(17) 
  pd(17,17) = pd(17,17) - rrt(492) * density(51) 
  pd(17,51) = pd(17,51) - rrt(492) * density(17) 
  pd(43,17) = pd(43,17) + rrt(492) * density(51) 
  pd(43,51) = pd(43,51) + rrt(492) * density(17) 
  pd(51,17) = pd(51,17) - rrt(492) * density(51) 
  pd(51,51) = pd(51,51) - rrt(492) * density(17) 
  pd(01,18) = pd(01,18) + rrt(493) * density(51) 
  pd(01,51) = pd(01,51) + rrt(493) * density(18) 
  pd(18,18) = pd(18,18) - rrt(493) * density(51) 
  pd(18,51) = pd(18,51) - rrt(493) * density(18) 
  pd(43,18) = pd(43,18) + rrt(493) * density(51) 
  pd(43,51) = pd(43,51) + rrt(493) * density(18) 
  pd(51,18) = pd(51,18) - rrt(493) * density(51) 
  pd(51,51) = pd(51,51) - rrt(493) * density(18) 
  pd(29,33) = pd(29,33) + rrt(494) * density(51) 
  pd(29,51) = pd(29,51) + rrt(494) * density(33) 
  pd(33,33) = pd(33,33) - rrt(494) * density(51) 
  pd(33,51) = pd(33,51) - rrt(494) * density(33) 
  pd(43,33) = pd(43,33) + rrt(494) * density(51) 
  pd(43,51) = pd(43,51) + rrt(494) * density(33) 
  pd(51,33) = pd(51,33) - rrt(494) * density(51) 
  pd(51,51) = pd(51,51) - rrt(494) * density(33) 
  pd(21,34) = pd(21,34) + rrt(495) * density(51) 
  pd(21,51) = pd(21,51) + rrt(495) * density(34) 
  pd(34,34) = pd(34,34) - rrt(495) * density(51) 
  pd(34,51) = pd(34,51) - rrt(495) * density(34) 
  pd(43,34) = pd(43,34) + rrt(495) * density(51) 
  pd(43,51) = pd(43,51) + rrt(495) * density(34) 
  pd(51,34) = pd(51,34) - rrt(495) * density(51) 
  pd(51,51) = pd(51,51) - rrt(495) * density(34) 
  pd(40,45) = pd(40,45) + rrt(496) * density(51) 
  pd(40,51) = pd(40,51) + rrt(496) * density(45) 
  pd(43,45) = pd(43,45) + rrt(496) * density(51) 
  pd(43,51) = pd(43,51) + rrt(496) * density(45) 
  pd(45,45) = pd(45,45) - rrt(496) * density(51) 
  pd(45,51) = pd(45,51) - rrt(496) * density(45) 
  pd(51,45) = pd(51,45) - rrt(496) * density(51) 
  pd(51,51) = pd(51,51) - rrt(496) * density(45) 
  pd(41,46) = pd(41,46) + rrt(497) * density(51) 
  pd(41,51) = pd(41,51) + rrt(497) * density(46) 
  pd(43,46) = pd(43,46) + rrt(497) * density(51) 
  pd(43,51) = pd(43,51) + rrt(497) * density(46) 
  pd(46,46) = pd(46,46) - rrt(497) * density(51) 
  pd(46,51) = pd(46,51) - rrt(497) * density(46) 
  pd(51,46) = pd(51,46) - rrt(497) * density(51) 
  pd(51,51) = pd(51,51) - rrt(497) * density(46) 
  pd(42,47) = pd(42,47) + rrt(498) * density(51) 
  pd(42,51) = pd(42,51) + rrt(498) * density(47) 
  pd(43,47) = pd(43,47) + rrt(498) * density(51) 
  pd(43,51) = pd(43,51) + rrt(498) * density(47) 
  pd(47,47) = pd(47,47) - rrt(498) * density(51) 
  pd(47,51) = pd(47,51) - rrt(498) * density(47) 
  pd(51,47) = pd(51,47) - rrt(498) * density(51) 
  pd(51,51) = pd(51,51) - rrt(498) * density(47) 
  pd(14,18) = pd(14,18) + rrt(499) * density(36) * 2.0d0
  pd(14,36) = pd(14,36) + rrt(499) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(499) * density(36) 
  pd(18,36) = pd(18,36) - rrt(499) * density(18) 
  pd(29,18) = pd(29,18) + rrt(499) * density(36) 
  pd(29,36) = pd(29,36) + rrt(499) * density(18) 
  pd(36,18) = pd(36,18) - rrt(499) * density(36) 
  pd(36,36) = pd(36,36) - rrt(499) * density(18) 
  pd(01,19) = pd(01,19) + rrt(500) * density(36) 
  pd(01,36) = pd(01,36) + rrt(500) * density(19) 
  pd(14,19) = pd(14,19) + rrt(500) * density(36) 
  pd(14,36) = pd(14,36) + rrt(500) * density(19) 
  pd(19,19) = pd(19,19) - rrt(500) * density(36) 
  pd(19,36) = pd(19,36) - rrt(500) * density(19) 
  pd(29,19) = pd(29,19) + rrt(500) * density(36) 
  pd(29,36) = pd(29,36) + rrt(500) * density(19) 
  pd(36,19) = pd(36,19) - rrt(500) * density(36) 
  pd(36,36) = pd(36,36) - rrt(500) * density(19) 
  pd(01,20) = pd(01,20) + rrt(501) * density(36) * 2.0d0
  pd(01,36) = pd(01,36) + rrt(501) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(501) * density(36) 
  pd(20,36) = pd(20,36) - rrt(501) * density(20) 
  pd(29,20) = pd(29,20) + rrt(501) * density(36) 
  pd(29,36) = pd(29,36) + rrt(501) * density(20) 
  pd(36,20) = pd(36,20) - rrt(501) * density(36) 
  pd(36,36) = pd(36,36) - rrt(501) * density(20) 
  pd(29,34) = pd(29,34) + rrt(502) * density(36) * 3.0d0
  pd(29,36) = pd(29,36) + rrt(502) * density(34) * 3.0d0
  pd(34,34) = pd(34,34) - rrt(502) * density(36) 
  pd(34,36) = pd(34,36) - rrt(502) * density(34) 
  pd(36,34) = pd(36,34) - rrt(502) * density(36) 
  pd(36,36) = pd(36,36) - rrt(502) * density(34) 
  pd(21,35) = pd(21,35) + rrt(503) * density(36) * 2.0d0
  pd(21,36) = pd(21,36) + rrt(503) * density(35) * 2.0d0
  pd(29,35) = pd(29,35) + rrt(503) * density(36) 
  pd(29,36) = pd(29,36) + rrt(503) * density(35) 
  pd(35,35) = pd(35,35) - rrt(503) * density(36) 
  pd(35,36) = pd(35,36) - rrt(503) * density(35) 
  pd(36,35) = pd(36,35) - rrt(503) * density(36) 
  pd(36,36) = pd(36,36) - rrt(503) * density(35) 
  pd(14,36) = pd(14,36) + rrt(504) * density(45) 
  pd(14,45) = pd(14,45) + rrt(504) * density(36) 
  pd(29,36) = pd(29,36) + rrt(504) * density(45) * 2.0d0
  pd(29,45) = pd(29,45) + rrt(504) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(504) * density(45) 
  pd(36,45) = pd(36,45) - rrt(504) * density(36) 
  pd(45,36) = pd(45,36) - rrt(504) * density(45) 
  pd(45,45) = pd(45,45) - rrt(504) * density(36) 
  pd(01,36) = pd(01,36) + rrt(505) * density(46) 
  pd(01,46) = pd(01,46) + rrt(505) * density(36) 
  pd(29,36) = pd(29,36) + rrt(505) * density(46) * 2.0d0
  pd(29,46) = pd(29,46) + rrt(505) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(505) * density(46) 
  pd(36,46) = pd(36,46) - rrt(505) * density(36) 
  pd(46,36) = pd(46,36) - rrt(505) * density(46) 
  pd(46,46) = pd(46,46) - rrt(505) * density(36) 
  pd(14,36) = pd(14,36) + rrt(506) * density(47) 
  pd(14,47) = pd(14,47) + rrt(506) * density(36) 
  pd(21,36) = pd(21,36) + rrt(506) * density(47) 
  pd(21,47) = pd(21,47) + rrt(506) * density(36) 
  pd(29,36) = pd(29,36) + rrt(506) * density(47) 
  pd(29,47) = pd(29,47) + rrt(506) * density(36) 
  pd(36,36) = pd(36,36) - rrt(506) * density(47) 
  pd(36,47) = pd(36,47) - rrt(506) * density(36) 
  pd(47,36) = pd(47,36) - rrt(506) * density(47) 
  pd(47,47) = pd(47,47) - rrt(506) * density(36) 
  pd(01,36) = pd(01,36) + rrt(507) * density(52) 
  pd(01,52) = pd(01,52) + rrt(507) * density(36) 
  pd(21,36) = pd(21,36) + rrt(507) * density(52) 
  pd(21,52) = pd(21,52) + rrt(507) * density(36) 
  pd(29,36) = pd(29,36) + rrt(507) * density(52) 
  pd(29,52) = pd(29,52) + rrt(507) * density(36) 
  pd(36,36) = pd(36,36) - rrt(507) * density(52) 
  pd(36,52) = pd(36,52) - rrt(507) * density(36) 
  pd(52,36) = pd(52,36) - rrt(507) * density(52) 
  pd(52,52) = pd(52,52) - rrt(507) * density(36) 
  pd(14,18) = pd(14,18) + rrt(508) * density(37) * 2.0d0
  pd(14,37) = pd(14,37) + rrt(508) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(508) * density(37) 
  pd(18,37) = pd(18,37) - rrt(508) * density(18) 
  pd(21,18) = pd(21,18) + rrt(508) * density(37) 
  pd(21,37) = pd(21,37) + rrt(508) * density(18) 
  pd(37,18) = pd(37,18) - rrt(508) * density(37) 
  pd(37,37) = pd(37,37) - rrt(508) * density(18) 
  pd(01,19) = pd(01,19) + rrt(509) * density(37) 
  pd(01,37) = pd(01,37) + rrt(509) * density(19) 
  pd(14,19) = pd(14,19) + rrt(509) * density(37) 
  pd(14,37) = pd(14,37) + rrt(509) * density(19) 
  pd(19,19) = pd(19,19) - rrt(509) * density(37) 
  pd(19,37) = pd(19,37) - rrt(509) * density(19) 
  pd(21,19) = pd(21,19) + rrt(509) * density(37) 
  pd(21,37) = pd(21,37) + rrt(509) * density(19) 
  pd(37,19) = pd(37,19) - rrt(509) * density(37) 
  pd(37,37) = pd(37,37) - rrt(509) * density(19) 
  pd(01,20) = pd(01,20) + rrt(510) * density(37) * 2.0d0
  pd(01,37) = pd(01,37) + rrt(510) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(510) * density(37) 
  pd(20,37) = pd(20,37) - rrt(510) * density(20) 
  pd(21,20) = pd(21,20) + rrt(510) * density(37) 
  pd(21,37) = pd(21,37) + rrt(510) * density(20) 
  pd(37,20) = pd(37,20) - rrt(510) * density(37) 
  pd(37,37) = pd(37,37) - rrt(510) * density(20) 
  pd(21,34) = pd(21,34) + rrt(511) * density(37) 
  pd(21,37) = pd(21,37) + rrt(511) * density(34) 
  pd(29,34) = pd(29,34) + rrt(511) * density(37) * 2.0d0
  pd(29,37) = pd(29,37) + rrt(511) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(511) * density(37) 
  pd(34,37) = pd(34,37) - rrt(511) * density(34) 
  pd(37,34) = pd(37,34) - rrt(511) * density(37) 
  pd(37,37) = pd(37,37) - rrt(511) * density(34) 
  pd(21,35) = pd(21,35) + rrt(512) * density(37) * 3.0d0
  pd(21,37) = pd(21,37) + rrt(512) * density(35) * 3.0d0
  pd(35,35) = pd(35,35) - rrt(512) * density(37) 
  pd(35,37) = pd(35,37) - rrt(512) * density(35) 
  pd(37,35) = pd(37,35) - rrt(512) * density(37) 
  pd(37,37) = pd(37,37) - rrt(512) * density(35) 
  pd(14,37) = pd(14,37) + rrt(513) * density(45) 
  pd(14,45) = pd(14,45) + rrt(513) * density(37) 
  pd(21,37) = pd(21,37) + rrt(513) * density(45) 
  pd(21,45) = pd(21,45) + rrt(513) * density(37) 
  pd(29,37) = pd(29,37) + rrt(513) * density(45) 
  pd(29,45) = pd(29,45) + rrt(513) * density(37) 
  pd(37,37) = pd(37,37) - rrt(513) * density(45) 
  pd(37,45) = pd(37,45) - rrt(513) * density(37) 
  pd(45,37) = pd(45,37) - rrt(513) * density(45) 
  pd(45,45) = pd(45,45) - rrt(513) * density(37) 
  pd(01,37) = pd(01,37) + rrt(514) * density(46) 
  pd(01,46) = pd(01,46) + rrt(514) * density(37) 
  pd(21,37) = pd(21,37) + rrt(514) * density(46) 
  pd(21,46) = pd(21,46) + rrt(514) * density(37) 
  pd(29,37) = pd(29,37) + rrt(514) * density(46) 
  pd(29,46) = pd(29,46) + rrt(514) * density(37) 
  pd(37,37) = pd(37,37) - rrt(514) * density(46) 
  pd(37,46) = pd(37,46) - rrt(514) * density(37) 
  pd(46,37) = pd(46,37) - rrt(514) * density(46) 
  pd(46,46) = pd(46,46) - rrt(514) * density(37) 
  pd(14,37) = pd(14,37) + rrt(515) * density(47) 
  pd(14,47) = pd(14,47) + rrt(515) * density(37) 
  pd(21,37) = pd(21,37) + rrt(515) * density(47) * 2.0d0
  pd(21,47) = pd(21,47) + rrt(515) * density(37) * 2.0d0
  pd(37,37) = pd(37,37) - rrt(515) * density(47) 
  pd(37,47) = pd(37,47) - rrt(515) * density(37) 
  pd(47,37) = pd(47,37) - rrt(515) * density(47) 
  pd(47,47) = pd(47,47) - rrt(515) * density(37) 
  pd(01,37) = pd(01,37) + rrt(516) * density(52) 
  pd(01,52) = pd(01,52) + rrt(516) * density(37) 
  pd(21,37) = pd(21,37) + rrt(516) * density(52) * 2.0d0
  pd(21,52) = pd(21,52) + rrt(516) * density(37) * 2.0d0
  pd(37,37) = pd(37,37) - rrt(516) * density(52) 
  pd(37,52) = pd(37,52) - rrt(516) * density(37) 
  pd(52,37) = pd(52,37) - rrt(516) * density(52) 
  pd(52,52) = pd(52,52) - rrt(516) * density(37) 
  pd(14,18) = pd(14,18) + rrt(517) * density(38) * 2.0d0
  pd(14,38) = pd(14,38) + rrt(517) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(517) * density(38) 
  pd(18,38) = pd(18,38) - rrt(517) * density(18) 
  pd(32,18) = pd(32,18) + rrt(517) * density(38) 
  pd(32,38) = pd(32,38) + rrt(517) * density(18) 
  pd(38,18) = pd(38,18) - rrt(517) * density(38) 
  pd(38,38) = pd(38,38) - rrt(517) * density(18) 
  pd(01,19) = pd(01,19) + rrt(518) * density(38) 
  pd(01,38) = pd(01,38) + rrt(518) * density(19) 
  pd(14,19) = pd(14,19) + rrt(518) * density(38) 
  pd(14,38) = pd(14,38) + rrt(518) * density(19) 
  pd(19,19) = pd(19,19) - rrt(518) * density(38) 
  pd(19,38) = pd(19,38) - rrt(518) * density(19) 
  pd(32,19) = pd(32,19) + rrt(518) * density(38) 
  pd(32,38) = pd(32,38) + rrt(518) * density(19) 
  pd(38,19) = pd(38,19) - rrt(518) * density(38) 
  pd(38,38) = pd(38,38) - rrt(518) * density(19) 
  pd(01,20) = pd(01,20) + rrt(519) * density(38) * 2.0d0
  pd(01,38) = pd(01,38) + rrt(519) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(519) * density(38) 
  pd(20,38) = pd(20,38) - rrt(519) * density(20) 
  pd(32,20) = pd(32,20) + rrt(519) * density(38) 
  pd(32,38) = pd(32,38) + rrt(519) * density(20) 
  pd(38,20) = pd(38,20) - rrt(519) * density(38) 
  pd(38,38) = pd(38,38) - rrt(519) * density(20) 
  pd(29,34) = pd(29,34) + rrt(520) * density(38) * 2.0d0
  pd(29,38) = pd(29,38) + rrt(520) * density(34) * 2.0d0
  pd(32,34) = pd(32,34) + rrt(520) * density(38) 
  pd(32,38) = pd(32,38) + rrt(520) * density(34) 
  pd(34,34) = pd(34,34) - rrt(520) * density(38) 
  pd(34,38) = pd(34,38) - rrt(520) * density(34) 
  pd(38,34) = pd(38,34) - rrt(520) * density(38) 
  pd(38,38) = pd(38,38) - rrt(520) * density(34) 
  pd(21,35) = pd(21,35) + rrt(521) * density(38) * 2.0d0
  pd(21,38) = pd(21,38) + rrt(521) * density(35) * 2.0d0
  pd(32,35) = pd(32,35) + rrt(521) * density(38) 
  pd(32,38) = pd(32,38) + rrt(521) * density(35) 
  pd(35,35) = pd(35,35) - rrt(521) * density(38) 
  pd(35,38) = pd(35,38) - rrt(521) * density(35) 
  pd(38,35) = pd(38,35) - rrt(521) * density(38) 
  pd(38,38) = pd(38,38) - rrt(521) * density(35) 
  pd(14,38) = pd(14,38) + rrt(522) * density(45) 
  pd(14,45) = pd(14,45) + rrt(522) * density(38) 
  pd(29,38) = pd(29,38) + rrt(522) * density(45) 
  pd(29,45) = pd(29,45) + rrt(522) * density(38) 
  pd(32,38) = pd(32,38) + rrt(522) * density(45) 
  pd(32,45) = pd(32,45) + rrt(522) * density(38) 
  pd(38,38) = pd(38,38) - rrt(522) * density(45) 
  pd(38,45) = pd(38,45) - rrt(522) * density(38) 
  pd(45,38) = pd(45,38) - rrt(522) * density(45) 
  pd(45,45) = pd(45,45) - rrt(522) * density(38) 
  pd(01,38) = pd(01,38) + rrt(523) * density(46) 
  pd(01,46) = pd(01,46) + rrt(523) * density(38) 
  pd(29,38) = pd(29,38) + rrt(523) * density(46) 
  pd(29,46) = pd(29,46) + rrt(523) * density(38) 
  pd(32,38) = pd(32,38) + rrt(523) * density(46) 
  pd(32,46) = pd(32,46) + rrt(523) * density(38) 
  pd(38,38) = pd(38,38) - rrt(523) * density(46) 
  pd(38,46) = pd(38,46) - rrt(523) * density(38) 
  pd(46,38) = pd(46,38) - rrt(523) * density(46) 
  pd(46,46) = pd(46,46) - rrt(523) * density(38) 
  pd(14,38) = pd(14,38) + rrt(524) * density(47) 
  pd(14,47) = pd(14,47) + rrt(524) * density(38) 
  pd(21,38) = pd(21,38) + rrt(524) * density(47) 
  pd(21,47) = pd(21,47) + rrt(524) * density(38) 
  pd(32,38) = pd(32,38) + rrt(524) * density(47) 
  pd(32,47) = pd(32,47) + rrt(524) * density(38) 
  pd(38,38) = pd(38,38) - rrt(524) * density(47) 
  pd(38,47) = pd(38,47) - rrt(524) * density(38) 
  pd(47,38) = pd(47,38) - rrt(524) * density(47) 
  pd(47,47) = pd(47,47) - rrt(524) * density(38) 
  pd(01,38) = pd(01,38) + rrt(525) * density(52) 
  pd(01,52) = pd(01,52) + rrt(525) * density(38) 
  pd(21,38) = pd(21,38) + rrt(525) * density(52) 
  pd(21,52) = pd(21,52) + rrt(525) * density(38) 
  pd(32,38) = pd(32,38) + rrt(525) * density(52) 
  pd(32,52) = pd(32,52) + rrt(525) * density(38) 
  pd(38,38) = pd(38,38) - rrt(525) * density(52) 
  pd(38,52) = pd(38,52) - rrt(525) * density(38) 
  pd(52,38) = pd(52,38) - rrt(525) * density(52) 
  pd(52,52) = pd(52,52) - rrt(525) * density(38) 
  pd(14,18) = pd(14,18) + rrt(526) * density(48) * 2.0d0
  pd(14,48) = pd(14,48) + rrt(526) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(526) * density(48) 
  pd(18,48) = pd(18,48) - rrt(526) * density(18) 
  pd(40,18) = pd(40,18) + rrt(526) * density(48) 
  pd(40,48) = pd(40,48) + rrt(526) * density(18) 
  pd(48,18) = pd(48,18) - rrt(526) * density(48) 
  pd(48,48) = pd(48,48) - rrt(526) * density(18) 
  pd(01,19) = pd(01,19) + rrt(527) * density(48) 
  pd(01,48) = pd(01,48) + rrt(527) * density(19) 
  pd(14,19) = pd(14,19) + rrt(527) * density(48) 
  pd(14,48) = pd(14,48) + rrt(527) * density(19) 
  pd(19,19) = pd(19,19) - rrt(527) * density(48) 
  pd(19,48) = pd(19,48) - rrt(527) * density(19) 
  pd(40,19) = pd(40,19) + rrt(527) * density(48) 
  pd(40,48) = pd(40,48) + rrt(527) * density(19) 
  pd(48,19) = pd(48,19) - rrt(527) * density(48) 
  pd(48,48) = pd(48,48) - rrt(527) * density(19) 
  pd(01,20) = pd(01,20) + rrt(528) * density(48) * 2.0d0
  pd(01,48) = pd(01,48) + rrt(528) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(528) * density(48) 
  pd(20,48) = pd(20,48) - rrt(528) * density(20) 
  pd(40,20) = pd(40,20) + rrt(528) * density(48) 
  pd(40,48) = pd(40,48) + rrt(528) * density(20) 
  pd(48,20) = pd(48,20) - rrt(528) * density(48) 
  pd(48,48) = pd(48,48) - rrt(528) * density(20) 
  pd(29,34) = pd(29,34) + rrt(529) * density(48) * 2.0d0
  pd(29,48) = pd(29,48) + rrt(529) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(529) * density(48) 
  pd(34,48) = pd(34,48) - rrt(529) * density(34) 
  pd(40,34) = pd(40,34) + rrt(529) * density(48) 
  pd(40,48) = pd(40,48) + rrt(529) * density(34) 
  pd(48,34) = pd(48,34) - rrt(529) * density(48) 
  pd(48,48) = pd(48,48) - rrt(529) * density(34) 
  pd(21,35) = pd(21,35) + rrt(530) * density(48) * 2.0d0
  pd(21,48) = pd(21,48) + rrt(530) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(530) * density(48) 
  pd(35,48) = pd(35,48) - rrt(530) * density(35) 
  pd(40,35) = pd(40,35) + rrt(530) * density(48) 
  pd(40,48) = pd(40,48) + rrt(530) * density(35) 
  pd(48,35) = pd(48,35) - rrt(530) * density(48) 
  pd(48,48) = pd(48,48) - rrt(530) * density(35) 
  pd(14,45) = pd(14,45) + rrt(531) * density(48) 
  pd(14,48) = pd(14,48) + rrt(531) * density(45) 
  pd(29,45) = pd(29,45) + rrt(531) * density(48) 
  pd(29,48) = pd(29,48) + rrt(531) * density(45) 
  pd(40,45) = pd(40,45) + rrt(531) * density(48) 
  pd(40,48) = pd(40,48) + rrt(531) * density(45) 
  pd(45,45) = pd(45,45) - rrt(531) * density(48) 
  pd(45,48) = pd(45,48) - rrt(531) * density(45) 
  pd(48,45) = pd(48,45) - rrt(531) * density(48) 
  pd(48,48) = pd(48,48) - rrt(531) * density(45) 
  pd(01,46) = pd(01,46) + rrt(532) * density(48) 
  pd(01,48) = pd(01,48) + rrt(532) * density(46) 
  pd(29,46) = pd(29,46) + rrt(532) * density(48) 
  pd(29,48) = pd(29,48) + rrt(532) * density(46) 
  pd(40,46) = pd(40,46) + rrt(532) * density(48) 
  pd(40,48) = pd(40,48) + rrt(532) * density(46) 
  pd(46,46) = pd(46,46) - rrt(532) * density(48) 
  pd(46,48) = pd(46,48) - rrt(532) * density(46) 
  pd(48,46) = pd(48,46) - rrt(532) * density(48) 
  pd(48,48) = pd(48,48) - rrt(532) * density(46) 
  pd(14,47) = pd(14,47) + rrt(533) * density(48) 
  pd(14,48) = pd(14,48) + rrt(533) * density(47) 
  pd(21,47) = pd(21,47) + rrt(533) * density(48) 
  pd(21,48) = pd(21,48) + rrt(533) * density(47) 
  pd(40,47) = pd(40,47) + rrt(533) * density(48) 
  pd(40,48) = pd(40,48) + rrt(533) * density(47) 
  pd(47,47) = pd(47,47) - rrt(533) * density(48) 
  pd(47,48) = pd(47,48) - rrt(533) * density(47) 
  pd(48,47) = pd(48,47) - rrt(533) * density(48) 
  pd(48,48) = pd(48,48) - rrt(533) * density(47) 
  pd(01,48) = pd(01,48) + rrt(534) * density(52) 
  pd(01,52) = pd(01,52) + rrt(534) * density(48) 
  pd(21,48) = pd(21,48) + rrt(534) * density(52) 
  pd(21,52) = pd(21,52) + rrt(534) * density(48) 
  pd(40,48) = pd(40,48) + rrt(534) * density(52) 
  pd(40,52) = pd(40,52) + rrt(534) * density(48) 
  pd(48,48) = pd(48,48) - rrt(534) * density(52) 
  pd(48,52) = pd(48,52) - rrt(534) * density(48) 
  pd(52,48) = pd(52,48) - rrt(534) * density(52) 
  pd(52,52) = pd(52,52) - rrt(534) * density(48) 
  pd(14,18) = pd(14,18) + rrt(535) * density(49) * 2.0d0
  pd(14,49) = pd(14,49) + rrt(535) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(535) * density(49) 
  pd(18,49) = pd(18,49) - rrt(535) * density(18) 
  pd(41,18) = pd(41,18) + rrt(535) * density(49) 
  pd(41,49) = pd(41,49) + rrt(535) * density(18) 
  pd(49,18) = pd(49,18) - rrt(535) * density(49) 
  pd(49,49) = pd(49,49) - rrt(535) * density(18) 
  pd(01,19) = pd(01,19) + rrt(536) * density(49) 
  pd(01,49) = pd(01,49) + rrt(536) * density(19) 
  pd(14,19) = pd(14,19) + rrt(536) * density(49) 
  pd(14,49) = pd(14,49) + rrt(536) * density(19) 
  pd(19,19) = pd(19,19) - rrt(536) * density(49) 
  pd(19,49) = pd(19,49) - rrt(536) * density(19) 
  pd(41,19) = pd(41,19) + rrt(536) * density(49) 
  pd(41,49) = pd(41,49) + rrt(536) * density(19) 
  pd(49,19) = pd(49,19) - rrt(536) * density(49) 
  pd(49,49) = pd(49,49) - rrt(536) * density(19) 
  pd(01,20) = pd(01,20) + rrt(537) * density(49) * 2.0d0
  pd(01,49) = pd(01,49) + rrt(537) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(537) * density(49) 
  pd(20,49) = pd(20,49) - rrt(537) * density(20) 
  pd(41,20) = pd(41,20) + rrt(537) * density(49) 
  pd(41,49) = pd(41,49) + rrt(537) * density(20) 
  pd(49,20) = pd(49,20) - rrt(537) * density(49) 
  pd(49,49) = pd(49,49) - rrt(537) * density(20) 
  pd(29,34) = pd(29,34) + rrt(538) * density(49) * 2.0d0
  pd(29,49) = pd(29,49) + rrt(538) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(538) * density(49) 
  pd(34,49) = pd(34,49) - rrt(538) * density(34) 
  pd(41,34) = pd(41,34) + rrt(538) * density(49) 
  pd(41,49) = pd(41,49) + rrt(538) * density(34) 
  pd(49,34) = pd(49,34) - rrt(538) * density(49) 
  pd(49,49) = pd(49,49) - rrt(538) * density(34) 
  pd(21,35) = pd(21,35) + rrt(539) * density(49) * 2.0d0
  pd(21,49) = pd(21,49) + rrt(539) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(539) * density(49) 
  pd(35,49) = pd(35,49) - rrt(539) * density(35) 
  pd(41,35) = pd(41,35) + rrt(539) * density(49) 
  pd(41,49) = pd(41,49) + rrt(539) * density(35) 
  pd(49,35) = pd(49,35) - rrt(539) * density(49) 
  pd(49,49) = pd(49,49) - rrt(539) * density(35) 
  pd(14,45) = pd(14,45) + rrt(540) * density(49) 
  pd(14,49) = pd(14,49) + rrt(540) * density(45) 
  pd(29,45) = pd(29,45) + rrt(540) * density(49) 
  pd(29,49) = pd(29,49) + rrt(540) * density(45) 
  pd(41,45) = pd(41,45) + rrt(540) * density(49) 
  pd(41,49) = pd(41,49) + rrt(540) * density(45) 
  pd(45,45) = pd(45,45) - rrt(540) * density(49) 
  pd(45,49) = pd(45,49) - rrt(540) * density(45) 
  pd(49,45) = pd(49,45) - rrt(540) * density(49) 
  pd(49,49) = pd(49,49) - rrt(540) * density(45) 
  pd(01,46) = pd(01,46) + rrt(541) * density(49) 
  pd(01,49) = pd(01,49) + rrt(541) * density(46) 
  pd(29,46) = pd(29,46) + rrt(541) * density(49) 
  pd(29,49) = pd(29,49) + rrt(541) * density(46) 
  pd(41,46) = pd(41,46) + rrt(541) * density(49) 
  pd(41,49) = pd(41,49) + rrt(541) * density(46) 
  pd(46,46) = pd(46,46) - rrt(541) * density(49) 
  pd(46,49) = pd(46,49) - rrt(541) * density(46) 
  pd(49,46) = pd(49,46) - rrt(541) * density(49) 
  pd(49,49) = pd(49,49) - rrt(541) * density(46) 
  pd(14,47) = pd(14,47) + rrt(542) * density(49) 
  pd(14,49) = pd(14,49) + rrt(542) * density(47) 
  pd(21,47) = pd(21,47) + rrt(542) * density(49) 
  pd(21,49) = pd(21,49) + rrt(542) * density(47) 
  pd(41,47) = pd(41,47) + rrt(542) * density(49) 
  pd(41,49) = pd(41,49) + rrt(542) * density(47) 
  pd(47,47) = pd(47,47) - rrt(542) * density(49) 
  pd(47,49) = pd(47,49) - rrt(542) * density(47) 
  pd(49,47) = pd(49,47) - rrt(542) * density(49) 
  pd(49,49) = pd(49,49) - rrt(542) * density(47) 
  pd(01,49) = pd(01,49) + rrt(543) * density(52) 
  pd(01,52) = pd(01,52) + rrt(543) * density(49) 
  pd(21,49) = pd(21,49) + rrt(543) * density(52) 
  pd(21,52) = pd(21,52) + rrt(543) * density(49) 
  pd(41,49) = pd(41,49) + rrt(543) * density(52) 
  pd(41,52) = pd(41,52) + rrt(543) * density(49) 
  pd(49,49) = pd(49,49) - rrt(543) * density(52) 
  pd(49,52) = pd(49,52) - rrt(543) * density(49) 
  pd(52,49) = pd(52,49) - rrt(543) * density(52) 
  pd(52,52) = pd(52,52) - rrt(543) * density(49) 
  pd(14,18) = pd(14,18) + rrt(544) * density(50) * 2.0d0
  pd(14,50) = pd(14,50) + rrt(544) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(544) * density(50) 
  pd(18,50) = pd(18,50) - rrt(544) * density(18) 
  pd(42,18) = pd(42,18) + rrt(544) * density(50) 
  pd(42,50) = pd(42,50) + rrt(544) * density(18) 
  pd(50,18) = pd(50,18) - rrt(544) * density(50) 
  pd(50,50) = pd(50,50) - rrt(544) * density(18) 
  pd(01,19) = pd(01,19) + rrt(545) * density(50) 
  pd(01,50) = pd(01,50) + rrt(545) * density(19) 
  pd(14,19) = pd(14,19) + rrt(545) * density(50) 
  pd(14,50) = pd(14,50) + rrt(545) * density(19) 
  pd(19,19) = pd(19,19) - rrt(545) * density(50) 
  pd(19,50) = pd(19,50) - rrt(545) * density(19) 
  pd(42,19) = pd(42,19) + rrt(545) * density(50) 
  pd(42,50) = pd(42,50) + rrt(545) * density(19) 
  pd(50,19) = pd(50,19) - rrt(545) * density(50) 
  pd(50,50) = pd(50,50) - rrt(545) * density(19) 
  pd(01,20) = pd(01,20) + rrt(546) * density(50) * 2.0d0
  pd(01,50) = pd(01,50) + rrt(546) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(546) * density(50) 
  pd(20,50) = pd(20,50) - rrt(546) * density(20) 
  pd(42,20) = pd(42,20) + rrt(546) * density(50) 
  pd(42,50) = pd(42,50) + rrt(546) * density(20) 
  pd(50,20) = pd(50,20) - rrt(546) * density(50) 
  pd(50,50) = pd(50,50) - rrt(546) * density(20) 
  pd(29,34) = pd(29,34) + rrt(547) * density(50) * 2.0d0
  pd(29,50) = pd(29,50) + rrt(547) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(547) * density(50) 
  pd(34,50) = pd(34,50) - rrt(547) * density(34) 
  pd(42,34) = pd(42,34) + rrt(547) * density(50) 
  pd(42,50) = pd(42,50) + rrt(547) * density(34) 
  pd(50,34) = pd(50,34) - rrt(547) * density(50) 
  pd(50,50) = pd(50,50) - rrt(547) * density(34) 
  pd(21,35) = pd(21,35) + rrt(548) * density(50) * 2.0d0
  pd(21,50) = pd(21,50) + rrt(548) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(548) * density(50) 
  pd(35,50) = pd(35,50) - rrt(548) * density(35) 
  pd(42,35) = pd(42,35) + rrt(548) * density(50) 
  pd(42,50) = pd(42,50) + rrt(548) * density(35) 
  pd(50,35) = pd(50,35) - rrt(548) * density(50) 
  pd(50,50) = pd(50,50) - rrt(548) * density(35) 
  pd(14,45) = pd(14,45) + rrt(549) * density(50) 
  pd(14,50) = pd(14,50) + rrt(549) * density(45) 
  pd(29,45) = pd(29,45) + rrt(549) * density(50) 
  pd(29,50) = pd(29,50) + rrt(549) * density(45) 
  pd(42,45) = pd(42,45) + rrt(549) * density(50) 
  pd(42,50) = pd(42,50) + rrt(549) * density(45) 
  pd(45,45) = pd(45,45) - rrt(549) * density(50) 
  pd(45,50) = pd(45,50) - rrt(549) * density(45) 
  pd(50,45) = pd(50,45) - rrt(549) * density(50) 
  pd(50,50) = pd(50,50) - rrt(549) * density(45) 
  pd(01,46) = pd(01,46) + rrt(550) * density(50) 
  pd(01,50) = pd(01,50) + rrt(550) * density(46) 
  pd(29,46) = pd(29,46) + rrt(550) * density(50) 
  pd(29,50) = pd(29,50) + rrt(550) * density(46) 
  pd(42,46) = pd(42,46) + rrt(550) * density(50) 
  pd(42,50) = pd(42,50) + rrt(550) * density(46) 
  pd(46,46) = pd(46,46) - rrt(550) * density(50) 
  pd(46,50) = pd(46,50) - rrt(550) * density(46) 
  pd(50,46) = pd(50,46) - rrt(550) * density(50) 
  pd(50,50) = pd(50,50) - rrt(550) * density(46) 
  pd(14,47) = pd(14,47) + rrt(551) * density(50) 
  pd(14,50) = pd(14,50) + rrt(551) * density(47) 
  pd(21,47) = pd(21,47) + rrt(551) * density(50) 
  pd(21,50) = pd(21,50) + rrt(551) * density(47) 
  pd(42,47) = pd(42,47) + rrt(551) * density(50) 
  pd(42,50) = pd(42,50) + rrt(551) * density(47) 
  pd(47,47) = pd(47,47) - rrt(551) * density(50) 
  pd(47,50) = pd(47,50) - rrt(551) * density(47) 
  pd(50,47) = pd(50,47) - rrt(551) * density(50) 
  pd(50,50) = pd(50,50) - rrt(551) * density(47) 
  pd(01,50) = pd(01,50) + rrt(552) * density(52) 
  pd(01,52) = pd(01,52) + rrt(552) * density(50) 
  pd(21,50) = pd(21,50) + rrt(552) * density(52) 
  pd(21,52) = pd(21,52) + rrt(552) * density(50) 
  pd(42,50) = pd(42,50) + rrt(552) * density(52) 
  pd(42,52) = pd(42,52) + rrt(552) * density(50) 
  pd(50,50) = pd(50,50) - rrt(552) * density(52) 
  pd(50,52) = pd(50,52) - rrt(552) * density(50) 
  pd(52,50) = pd(52,50) - rrt(552) * density(52) 
  pd(52,52) = pd(52,52) - rrt(552) * density(50) 
  pd(14,18) = pd(14,18) + rrt(553) * density(51) * 2.0d0
  pd(14,51) = pd(14,51) + rrt(553) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(553) * density(51) 
  pd(18,51) = pd(18,51) - rrt(553) * density(18) 
  pd(43,18) = pd(43,18) + rrt(553) * density(51) 
  pd(43,51) = pd(43,51) + rrt(553) * density(18) 
  pd(51,18) = pd(51,18) - rrt(553) * density(51) 
  pd(51,51) = pd(51,51) - rrt(553) * density(18) 
  pd(01,19) = pd(01,19) + rrt(554) * density(51) 
  pd(01,51) = pd(01,51) + rrt(554) * density(19) 
  pd(14,19) = pd(14,19) + rrt(554) * density(51) 
  pd(14,51) = pd(14,51) + rrt(554) * density(19) 
  pd(19,19) = pd(19,19) - rrt(554) * density(51) 
  pd(19,51) = pd(19,51) - rrt(554) * density(19) 
  pd(43,19) = pd(43,19) + rrt(554) * density(51) 
  pd(43,51) = pd(43,51) + rrt(554) * density(19) 
  pd(51,19) = pd(51,19) - rrt(554) * density(51) 
  pd(51,51) = pd(51,51) - rrt(554) * density(19) 
  pd(01,20) = pd(01,20) + rrt(555) * density(51) * 2.0d0
  pd(01,51) = pd(01,51) + rrt(555) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(555) * density(51) 
  pd(20,51) = pd(20,51) - rrt(555) * density(20) 
  pd(43,20) = pd(43,20) + rrt(555) * density(51) 
  pd(43,51) = pd(43,51) + rrt(555) * density(20) 
  pd(51,20) = pd(51,20) - rrt(555) * density(51) 
  pd(51,51) = pd(51,51) - rrt(555) * density(20) 
  pd(29,34) = pd(29,34) + rrt(556) * density(51) * 2.0d0
  pd(29,51) = pd(29,51) + rrt(556) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(556) * density(51) 
  pd(34,51) = pd(34,51) - rrt(556) * density(34) 
  pd(43,34) = pd(43,34) + rrt(556) * density(51) 
  pd(43,51) = pd(43,51) + rrt(556) * density(34) 
  pd(51,34) = pd(51,34) - rrt(556) * density(51) 
  pd(51,51) = pd(51,51) - rrt(556) * density(34) 
  pd(21,35) = pd(21,35) + rrt(557) * density(51) * 2.0d0
  pd(21,51) = pd(21,51) + rrt(557) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(557) * density(51) 
  pd(35,51) = pd(35,51) - rrt(557) * density(35) 
  pd(43,35) = pd(43,35) + rrt(557) * density(51) 
  pd(43,51) = pd(43,51) + rrt(557) * density(35) 
  pd(51,35) = pd(51,35) - rrt(557) * density(51) 
  pd(51,51) = pd(51,51) - rrt(557) * density(35) 
  pd(14,45) = pd(14,45) + rrt(558) * density(51) 
  pd(14,51) = pd(14,51) + rrt(558) * density(45) 
  pd(29,45) = pd(29,45) + rrt(558) * density(51) 
  pd(29,51) = pd(29,51) + rrt(558) * density(45) 
  pd(43,45) = pd(43,45) + rrt(558) * density(51) 
  pd(43,51) = pd(43,51) + rrt(558) * density(45) 
  pd(45,45) = pd(45,45) - rrt(558) * density(51) 
  pd(45,51) = pd(45,51) - rrt(558) * density(45) 
  pd(51,45) = pd(51,45) - rrt(558) * density(51) 
  pd(51,51) = pd(51,51) - rrt(558) * density(45) 
  pd(01,46) = pd(01,46) + rrt(559) * density(51) 
  pd(01,51) = pd(01,51) + rrt(559) * density(46) 
  pd(29,46) = pd(29,46) + rrt(559) * density(51) 
  pd(29,51) = pd(29,51) + rrt(559) * density(46) 
  pd(43,46) = pd(43,46) + rrt(559) * density(51) 
  pd(43,51) = pd(43,51) + rrt(559) * density(46) 
  pd(46,46) = pd(46,46) - rrt(559) * density(51) 
  pd(46,51) = pd(46,51) - rrt(559) * density(46) 
  pd(51,46) = pd(51,46) - rrt(559) * density(51) 
  pd(51,51) = pd(51,51) - rrt(559) * density(46) 
  pd(14,47) = pd(14,47) + rrt(560) * density(51) 
  pd(14,51) = pd(14,51) + rrt(560) * density(47) 
  pd(21,47) = pd(21,47) + rrt(560) * density(51) 
  pd(21,51) = pd(21,51) + rrt(560) * density(47) 
  pd(43,47) = pd(43,47) + rrt(560) * density(51) 
  pd(43,51) = pd(43,51) + rrt(560) * density(47) 
  pd(47,47) = pd(47,47) - rrt(560) * density(51) 
  pd(47,51) = pd(47,51) - rrt(560) * density(47) 
  pd(51,47) = pd(51,47) - rrt(560) * density(51) 
  pd(51,51) = pd(51,51) - rrt(560) * density(47) 
  pd(01,51) = pd(01,51) + rrt(561) * density(52) 
  pd(01,52) = pd(01,52) + rrt(561) * density(51) 
  pd(21,51) = pd(21,51) + rrt(561) * density(52) 
  pd(21,52) = pd(21,52) + rrt(561) * density(51) 
  pd(43,51) = pd(43,51) + rrt(561) * density(52) 
  pd(43,52) = pd(43,52) + rrt(561) * density(51) 
  pd(51,51) = pd(51,51) - rrt(561) * density(52) 
  pd(51,52) = pd(51,52) - rrt(561) * density(51) 
  pd(52,51) = pd(52,51) - rrt(561) * density(52) 
  pd(52,52) = pd(52,52) - rrt(561) * density(51) 
  pd(14,17) = pd(14,17) + rrt(562) * density(39) 
  pd(14,39) = pd(14,39) + rrt(562) * density(17) 
  pd(17,17) = pd(17,17) - rrt(562) * density(39) 
  pd(17,39) = pd(17,39) - rrt(562) * density(17) 
  pd(21,17) = pd(21,17) + rrt(562) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(562) * density(17) * 2.0d0
  pd(39,17) = pd(39,17) - rrt(562) * density(39) 
  pd(39,39) = pd(39,39) - rrt(562) * density(17) 
  pd(01,18) = pd(01,18) + rrt(563) * density(39) 
  pd(01,39) = pd(01,39) + rrt(563) * density(18) 
  pd(18,18) = pd(18,18) - rrt(563) * density(39) 
  pd(18,39) = pd(18,39) - rrt(563) * density(18) 
  pd(21,18) = pd(21,18) + rrt(563) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(563) * density(18) * 2.0d0
  pd(39,18) = pd(39,18) - rrt(563) * density(39) 
  pd(39,39) = pd(39,39) - rrt(563) * density(18) 
  pd(21,33) = pd(21,33) + rrt(564) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(564) * density(33) * 2.0d0
  pd(29,33) = pd(29,33) + rrt(564) * density(39) 
  pd(29,39) = pd(29,39) + rrt(564) * density(33) 
  pd(33,33) = pd(33,33) - rrt(564) * density(39) 
  pd(33,39) = pd(33,39) - rrt(564) * density(33) 
  pd(39,33) = pd(39,33) - rrt(564) * density(39) 
  pd(39,39) = pd(39,39) - rrt(564) * density(33) 
  pd(21,34) = pd(21,34) + rrt(565) * density(39) * 3.0d0
  pd(21,39) = pd(21,39) + rrt(565) * density(34) * 3.0d0
  pd(34,34) = pd(34,34) - rrt(565) * density(39) 
  pd(34,39) = pd(34,39) - rrt(565) * density(34) 
  pd(39,34) = pd(39,34) - rrt(565) * density(39) 
  pd(39,39) = pd(39,39) - rrt(565) * density(34) 
  pd(21,39) = pd(21,39) + rrt(566) * density(45) * 2.0d0
  pd(21,45) = pd(21,45) + rrt(566) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(566) * density(45) 
  pd(39,45) = pd(39,45) - rrt(566) * density(39) 
  pd(40,39) = pd(40,39) + rrt(566) * density(45) 
  pd(40,45) = pd(40,45) + rrt(566) * density(39) 
  pd(45,39) = pd(45,39) - rrt(566) * density(45) 
  pd(45,45) = pd(45,45) - rrt(566) * density(39) 
  pd(21,39) = pd(21,39) + rrt(567) * density(46) * 2.0d0
  pd(21,46) = pd(21,46) + rrt(567) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(567) * density(46) 
  pd(39,46) = pd(39,46) - rrt(567) * density(39) 
  pd(41,39) = pd(41,39) + rrt(567) * density(46) 
  pd(41,46) = pd(41,46) + rrt(567) * density(39) 
  pd(46,39) = pd(46,39) - rrt(567) * density(46) 
  pd(46,46) = pd(46,46) - rrt(567) * density(39) 
  pd(21,39) = pd(21,39) + rrt(568) * density(47) * 2.0d0
  pd(21,47) = pd(21,47) + rrt(568) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(568) * density(47) 
  pd(39,47) = pd(39,47) - rrt(568) * density(39) 
  pd(42,39) = pd(42,39) + rrt(568) * density(47) 
  pd(42,47) = pd(42,47) + rrt(568) * density(39) 
  pd(47,39) = pd(47,39) - rrt(568) * density(47) 
  pd(47,47) = pd(47,47) - rrt(568) * density(39) 
  pd(01,19) = pd(01,19) + rrt(569) * density(39) 
  pd(01,39) = pd(01,39) + rrt(569) * density(19) 
  pd(14,19) = pd(14,19) + rrt(569) * density(39) 
  pd(14,39) = pd(14,39) + rrt(569) * density(19) 
  pd(19,19) = pd(19,19) - rrt(569) * density(39) 
  pd(19,39) = pd(19,39) - rrt(569) * density(19) 
  pd(21,19) = pd(21,19) + rrt(569) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(569) * density(19) * 2.0d0
  pd(39,19) = pd(39,19) - rrt(569) * density(39) 
  pd(39,39) = pd(39,39) - rrt(569) * density(19) 
  pd(01,20) = pd(01,20) + rrt(570) * density(39) * 2.0d0
  pd(01,39) = pd(01,39) + rrt(570) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(570) * density(39) 
  pd(20,39) = pd(20,39) - rrt(570) * density(20) 
  pd(21,20) = pd(21,20) + rrt(570) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(570) * density(20) * 2.0d0
  pd(39,20) = pd(39,20) - rrt(570) * density(39) 
  pd(39,39) = pd(39,39) - rrt(570) * density(20) 
  pd(21,35) = pd(21,35) + rrt(571) * density(39) * 4.0d0
  pd(21,39) = pd(21,39) + rrt(571) * density(35) * 4.0d0
  pd(35,35) = pd(35,35) - rrt(571) * density(39) 
  pd(35,39) = pd(35,39) - rrt(571) * density(35) 
  pd(39,35) = pd(39,35) - rrt(571) * density(39) 
  pd(39,39) = pd(39,39) - rrt(571) * density(35) 
  pd(01,39) = pd(01,39) + rrt(572) * density(52) 
  pd(01,52) = pd(01,52) + rrt(572) * density(39) 
  pd(21,39) = pd(21,39) + rrt(572) * density(52) * 3.0d0
  pd(21,52) = pd(21,52) + rrt(572) * density(39) * 3.0d0
  pd(39,39) = pd(39,39) - rrt(572) * density(52) 
  pd(39,52) = pd(39,52) - rrt(572) * density(39) 
  pd(52,39) = pd(52,39) - rrt(572) * density(52) 
  pd(52,52) = pd(52,52) - rrt(572) * density(39) 
  pd(14,17) = pd(14,17) + rrt(573) * density(36) 
  pd(14,36) = pd(14,36) + rrt(573) * density(17) 
  pd(17,17) = pd(17,17) - rrt(573) * density(36) 
  pd(17,36) = pd(17,36) - rrt(573) * density(17) 
  pd(29,17) = pd(29,17) + rrt(573) * density(36) 
  pd(29,36) = pd(29,36) + rrt(573) * density(17) 
  pd(36,17) = pd(36,17) - rrt(573) * density(36) 
  pd(36,36) = pd(36,36) - rrt(573) * density(17) 
  pd(01,18) = pd(01,18) + rrt(574) * density(36) 
  pd(01,36) = pd(01,36) + rrt(574) * density(18) 
  pd(18,18) = pd(18,18) - rrt(574) * density(36) 
  pd(18,36) = pd(18,36) - rrt(574) * density(18) 
  pd(29,18) = pd(29,18) + rrt(574) * density(36) 
  pd(29,36) = pd(29,36) + rrt(574) * density(18) 
  pd(36,18) = pd(36,18) - rrt(574) * density(36) 
  pd(36,36) = pd(36,36) - rrt(574) * density(18) 
  pd(29,33) = pd(29,33) + rrt(575) * density(36) * 2.0d0
  pd(29,36) = pd(29,36) + rrt(575) * density(33) * 2.0d0
  pd(33,33) = pd(33,33) - rrt(575) * density(36) 
  pd(33,36) = pd(33,36) - rrt(575) * density(33) 
  pd(36,33) = pd(36,33) - rrt(575) * density(36) 
  pd(36,36) = pd(36,36) - rrt(575) * density(33) 
  pd(21,34) = pd(21,34) + rrt(576) * density(36) 
  pd(21,36) = pd(21,36) + rrt(576) * density(34) 
  pd(29,34) = pd(29,34) + rrt(576) * density(36) 
  pd(29,36) = pd(29,36) + rrt(576) * density(34) 
  pd(34,34) = pd(34,34) - rrt(576) * density(36) 
  pd(34,36) = pd(34,36) - rrt(576) * density(34) 
  pd(36,34) = pd(36,34) - rrt(576) * density(36) 
  pd(36,36) = pd(36,36) - rrt(576) * density(34) 
  pd(29,36) = pd(29,36) + rrt(577) * density(45) 
  pd(29,45) = pd(29,45) + rrt(577) * density(36) 
  pd(36,36) = pd(36,36) - rrt(577) * density(45) 
  pd(36,45) = pd(36,45) - rrt(577) * density(36) 
  pd(40,36) = pd(40,36) + rrt(577) * density(45) 
  pd(40,45) = pd(40,45) + rrt(577) * density(36) 
  pd(45,36) = pd(45,36) - rrt(577) * density(45) 
  pd(45,45) = pd(45,45) - rrt(577) * density(36) 
  pd(14,17) = pd(14,17) + rrt(578) * density(37) 
  pd(14,37) = pd(14,37) + rrt(578) * density(17) 
  pd(17,17) = pd(17,17) - rrt(578) * density(37) 
  pd(17,37) = pd(17,37) - rrt(578) * density(17) 
  pd(21,17) = pd(21,17) + rrt(578) * density(37) 
  pd(21,37) = pd(21,37) + rrt(578) * density(17) 
  pd(37,17) = pd(37,17) - rrt(578) * density(37) 
  pd(37,37) = pd(37,37) - rrt(578) * density(17) 
  pd(01,18) = pd(01,18) + rrt(579) * density(37) 
  pd(01,37) = pd(01,37) + rrt(579) * density(18) 
  pd(18,18) = pd(18,18) - rrt(579) * density(37) 
  pd(18,37) = pd(18,37) - rrt(579) * density(18) 
  pd(21,18) = pd(21,18) + rrt(579) * density(37) 
  pd(21,37) = pd(21,37) + rrt(579) * density(18) 
  pd(37,18) = pd(37,18) - rrt(579) * density(37) 
  pd(37,37) = pd(37,37) - rrt(579) * density(18) 
  pd(21,33) = pd(21,33) + rrt(580) * density(37) 
  pd(21,37) = pd(21,37) + rrt(580) * density(33) 
  pd(29,33) = pd(29,33) + rrt(580) * density(37) 
  pd(29,37) = pd(29,37) + rrt(580) * density(33) 
  pd(33,33) = pd(33,33) - rrt(580) * density(37) 
  pd(33,37) = pd(33,37) - rrt(580) * density(33) 
  pd(37,33) = pd(37,33) - rrt(580) * density(37) 
  pd(37,37) = pd(37,37) - rrt(580) * density(33) 
  pd(21,34) = pd(21,34) + rrt(581) * density(37) * 2.0d0
  pd(21,37) = pd(21,37) + rrt(581) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(581) * density(37) 
  pd(34,37) = pd(34,37) - rrt(581) * density(34) 
  pd(37,34) = pd(37,34) - rrt(581) * density(37) 
  pd(37,37) = pd(37,37) - rrt(581) * density(34) 
  pd(21,37) = pd(21,37) + rrt(582) * density(45) 
  pd(21,45) = pd(21,45) + rrt(582) * density(37) 
  pd(37,37) = pd(37,37) - rrt(582) * density(45) 
  pd(37,45) = pd(37,45) - rrt(582) * density(37) 
  pd(40,37) = pd(40,37) + rrt(582) * density(45) 
  pd(40,45) = pd(40,45) + rrt(582) * density(37) 
  pd(45,37) = pd(45,37) - rrt(582) * density(45) 
  pd(45,45) = pd(45,45) - rrt(582) * density(37) 
  pd(17,17) = pd(17,17) - rrt(583) * density(36) 
  pd(17,36) = pd(17,36) - rrt(583) * density(17) 
  pd(36,17) = pd(36,17) - rrt(583) * density(36) 
  pd(36,36) = pd(36,36) - rrt(583) * density(17) 
  pd(40,17) = pd(40,17) + rrt(583) * density(36) 
  pd(40,36) = pd(40,36) + rrt(583) * density(17) 
  pd(18,18) = pd(18,18) - rrt(584) * density(36) 
  pd(18,36) = pd(18,36) - rrt(584) * density(18) 
  pd(36,18) = pd(36,18) - rrt(584) * density(36) 
  pd(36,36) = pd(36,36) - rrt(584) * density(18) 
  pd(41,18) = pd(41,18) + rrt(584) * density(36) 
  pd(41,36) = pd(41,36) + rrt(584) * density(18) 
  pd(21,33) = pd(21,33) + rrt(585) * density(36) 
  pd(21,36) = pd(21,36) + rrt(585) * density(33) 
  pd(33,33) = pd(33,33) - rrt(585) * density(36) 
  pd(33,36) = pd(33,36) - rrt(585) * density(33) 
  pd(36,33) = pd(36,33) - rrt(585) * density(36) 
  pd(36,36) = pd(36,36) - rrt(585) * density(33) 
  pd(32,34) = pd(32,34) + rrt(586) * density(36) 
  pd(32,36) = pd(32,36) + rrt(586) * density(34) 
  pd(34,34) = pd(34,34) - rrt(586) * density(36) 
  pd(34,36) = pd(34,36) - rrt(586) * density(34) 
  pd(36,34) = pd(36,34) - rrt(586) * density(36) 
  pd(36,36) = pd(36,36) - rrt(586) * density(34) 
  pd(36,36) = pd(36,36) - rrt(587) * density(45) 
  pd(36,45) = pd(36,45) - rrt(587) * density(36) 
  pd(42,36) = pd(42,36) + rrt(587) * density(45) 
  pd(42,45) = pd(42,45) + rrt(587) * density(36) 
  pd(45,36) = pd(45,36) - rrt(587) * density(45) 
  pd(45,45) = pd(45,45) - rrt(587) * density(36) 
  pd(17,17) = pd(17,17) - rrt(588) * density(37) 
  pd(17,37) = pd(17,37) - rrt(588) * density(17) 
  pd(37,17) = pd(37,17) - rrt(588) * density(37) 
  pd(37,37) = pd(37,37) - rrt(588) * density(17) 
  pd(42,17) = pd(42,17) + rrt(588) * density(37) 
  pd(42,37) = pd(42,37) + rrt(588) * density(17) 
  pd(32,33) = pd(32,33) + rrt(589) * density(37) 
  pd(32,37) = pd(32,37) + rrt(589) * density(33) 
  pd(33,33) = pd(33,33) - rrt(589) * density(37) 
  pd(33,37) = pd(33,37) - rrt(589) * density(33) 
  pd(37,33) = pd(37,33) - rrt(589) * density(37) 
  pd(37,37) = pd(37,37) - rrt(589) * density(33) 
  pd(37,37) = pd(37,37) - rrt(590) * density(45) 
  pd(37,45) = pd(37,45) - rrt(590) * density(37) 
  pd(43,37) = pd(43,37) + rrt(590) * density(45) 
  pd(43,45) = pd(43,45) + rrt(590) * density(37) 
  pd(45,37) = pd(45,37) - rrt(590) * density(45) 
  pd(45,45) = pd(45,45) - rrt(590) * density(37) 
  pd(14,17) = pd(14,17) + rrt(591) * density(38) 
  pd(14,38) = pd(14,38) + rrt(591) * density(17) 
  pd(17,17) = pd(17,17) - rrt(591) * density(38) 
  pd(17,38) = pd(17,38) - rrt(591) * density(17) 
  pd(32,17) = pd(32,17) + rrt(591) * density(38) 
  pd(32,38) = pd(32,38) + rrt(591) * density(17) 
  pd(38,17) = pd(38,17) - rrt(591) * density(38) 
  pd(38,38) = pd(38,38) - rrt(591) * density(17) 
  pd(01,18) = pd(01,18) + rrt(592) * density(38) 
  pd(01,38) = pd(01,38) + rrt(592) * density(18) 
  pd(18,18) = pd(18,18) - rrt(592) * density(38) 
  pd(18,38) = pd(18,38) - rrt(592) * density(18) 
  pd(32,18) = pd(32,18) + rrt(592) * density(38) 
  pd(32,38) = pd(32,38) + rrt(592) * density(18) 
  pd(38,18) = pd(38,18) - rrt(592) * density(38) 
  pd(38,38) = pd(38,38) - rrt(592) * density(18) 
  pd(29,33) = pd(29,33) + rrt(593) * density(38) 
  pd(29,38) = pd(29,38) + rrt(593) * density(33) 
  pd(32,33) = pd(32,33) + rrt(593) * density(38) 
  pd(32,38) = pd(32,38) + rrt(593) * density(33) 
  pd(33,33) = pd(33,33) - rrt(593) * density(38) 
  pd(33,38) = pd(33,38) - rrt(593) * density(33) 
  pd(38,33) = pd(38,33) - rrt(593) * density(38) 
  pd(38,38) = pd(38,38) - rrt(593) * density(33) 
  pd(21,34) = pd(21,34) + rrt(594) * density(38) 
  pd(21,38) = pd(21,38) + rrt(594) * density(34) 
  pd(32,34) = pd(32,34) + rrt(594) * density(38) 
  pd(32,38) = pd(32,38) + rrt(594) * density(34) 
  pd(34,34) = pd(34,34) - rrt(594) * density(38) 
  pd(34,38) = pd(34,38) - rrt(594) * density(34) 
  pd(38,34) = pd(38,34) - rrt(594) * density(38) 
  pd(38,38) = pd(38,38) - rrt(594) * density(34) 
  pd(32,38) = pd(32,38) + rrt(595) * density(45) 
  pd(32,45) = pd(32,45) + rrt(595) * density(38) 
  pd(38,38) = pd(38,38) - rrt(595) * density(45) 
  pd(38,45) = pd(38,45) - rrt(595) * density(38) 
  pd(40,38) = pd(40,38) + rrt(595) * density(45) 
  pd(40,45) = pd(40,45) + rrt(595) * density(38) 
  pd(45,38) = pd(45,38) - rrt(595) * density(45) 
  pd(45,45) = pd(45,45) - rrt(595) * density(38) 
  pd(32,38) = pd(32,38) + rrt(596) * density(46) 
  pd(32,46) = pd(32,46) + rrt(596) * density(38) 
  pd(38,38) = pd(38,38) - rrt(596) * density(46) 
  pd(38,46) = pd(38,46) - rrt(596) * density(38) 
  pd(41,38) = pd(41,38) + rrt(596) * density(46) 
  pd(41,46) = pd(41,46) + rrt(596) * density(38) 
  pd(46,38) = pd(46,38) - rrt(596) * density(46) 
  pd(46,46) = pd(46,46) - rrt(596) * density(38) 
  pd(32,38) = pd(32,38) + rrt(597) * density(47) 
  pd(32,47) = pd(32,47) + rrt(597) * density(38) 
  pd(38,38) = pd(38,38) - rrt(597) * density(47) 
  pd(38,47) = pd(38,47) - rrt(597) * density(38) 
  pd(42,38) = pd(42,38) + rrt(597) * density(47) 
  pd(42,47) = pd(42,47) + rrt(597) * density(38) 
  pd(47,38) = pd(47,38) - rrt(597) * density(47) 
  pd(47,47) = pd(47,47) - rrt(597) * density(38) 
  pd(14,17) = pd(14,17) + rrt(598) * density(48) 
  pd(14,48) = pd(14,48) + rrt(598) * density(17) 
  pd(17,17) = pd(17,17) - rrt(598) * density(48) 
  pd(17,48) = pd(17,48) - rrt(598) * density(17) 
  pd(40,17) = pd(40,17) + rrt(598) * density(48) 
  pd(40,48) = pd(40,48) + rrt(598) * density(17) 
  pd(48,17) = pd(48,17) - rrt(598) * density(48) 
  pd(48,48) = pd(48,48) - rrt(598) * density(17) 
  pd(01,18) = pd(01,18) + rrt(599) * density(48) 
  pd(01,48) = pd(01,48) + rrt(599) * density(18) 
  pd(18,18) = pd(18,18) - rrt(599) * density(48) 
  pd(18,48) = pd(18,48) - rrt(599) * density(18) 
  pd(40,18) = pd(40,18) + rrt(599) * density(48) 
  pd(40,48) = pd(40,48) + rrt(599) * density(18) 
  pd(48,18) = pd(48,18) - rrt(599) * density(48) 
  pd(48,48) = pd(48,48) - rrt(599) * density(18) 
  pd(29,33) = pd(29,33) + rrt(600) * density(48) 
  pd(29,48) = pd(29,48) + rrt(600) * density(33) 
  pd(33,33) = pd(33,33) - rrt(600) * density(48) 
  pd(33,48) = pd(33,48) - rrt(600) * density(33) 
  pd(40,33) = pd(40,33) + rrt(600) * density(48) 
  pd(40,48) = pd(40,48) + rrt(600) * density(33) 
  pd(48,33) = pd(48,33) - rrt(600) * density(48) 
  pd(48,48) = pd(48,48) - rrt(600) * density(33) 
  pd(21,34) = pd(21,34) + rrt(601) * density(48) 
  pd(21,48) = pd(21,48) + rrt(601) * density(34) 
  pd(34,34) = pd(34,34) - rrt(601) * density(48) 
  pd(34,48) = pd(34,48) - rrt(601) * density(34) 
  pd(40,34) = pd(40,34) + rrt(601) * density(48) 
  pd(40,48) = pd(40,48) + rrt(601) * density(34) 
  pd(48,34) = pd(48,34) - rrt(601) * density(48) 
  pd(48,48) = pd(48,48) - rrt(601) * density(34) 
  pd(40,45) = pd(40,45) + rrt(602) * density(48) * 2.0d0
  pd(40,48) = pd(40,48) + rrt(602) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(602) * density(48) 
  pd(45,48) = pd(45,48) - rrt(602) * density(45) 
  pd(48,45) = pd(48,45) - rrt(602) * density(48) 
  pd(48,48) = pd(48,48) - rrt(602) * density(45) 
  pd(40,46) = pd(40,46) + rrt(603) * density(48) 
  pd(40,48) = pd(40,48) + rrt(603) * density(46) 
  pd(41,46) = pd(41,46) + rrt(603) * density(48) 
  pd(41,48) = pd(41,48) + rrt(603) * density(46) 
  pd(46,46) = pd(46,46) - rrt(603) * density(48) 
  pd(46,48) = pd(46,48) - rrt(603) * density(46) 
  pd(48,46) = pd(48,46) - rrt(603) * density(48) 
  pd(48,48) = pd(48,48) - rrt(603) * density(46) 
  pd(40,47) = pd(40,47) + rrt(604) * density(48) 
  pd(40,48) = pd(40,48) + rrt(604) * density(47) 
  pd(42,47) = pd(42,47) + rrt(604) * density(48) 
  pd(42,48) = pd(42,48) + rrt(604) * density(47) 
  pd(47,47) = pd(47,47) - rrt(604) * density(48) 
  pd(47,48) = pd(47,48) - rrt(604) * density(47) 
  pd(48,47) = pd(48,47) - rrt(604) * density(48) 
  pd(48,48) = pd(48,48) - rrt(604) * density(47) 
  pd(14,17) = pd(14,17) + rrt(605) * density(49) 
  pd(14,49) = pd(14,49) + rrt(605) * density(17) 
  pd(17,17) = pd(17,17) - rrt(605) * density(49) 
  pd(17,49) = pd(17,49) - rrt(605) * density(17) 
  pd(41,17) = pd(41,17) + rrt(605) * density(49) 
  pd(41,49) = pd(41,49) + rrt(605) * density(17) 
  pd(49,17) = pd(49,17) - rrt(605) * density(49) 
  pd(49,49) = pd(49,49) - rrt(605) * density(17) 
  pd(01,18) = pd(01,18) + rrt(606) * density(49) 
  pd(01,49) = pd(01,49) + rrt(606) * density(18) 
  pd(18,18) = pd(18,18) - rrt(606) * density(49) 
  pd(18,49) = pd(18,49) - rrt(606) * density(18) 
  pd(41,18) = pd(41,18) + rrt(606) * density(49) 
  pd(41,49) = pd(41,49) + rrt(606) * density(18) 
  pd(49,18) = pd(49,18) - rrt(606) * density(49) 
  pd(49,49) = pd(49,49) - rrt(606) * density(18) 
  pd(29,33) = pd(29,33) + rrt(607) * density(49) 
  pd(29,49) = pd(29,49) + rrt(607) * density(33) 
  pd(33,33) = pd(33,33) - rrt(607) * density(49) 
  pd(33,49) = pd(33,49) - rrt(607) * density(33) 
  pd(41,33) = pd(41,33) + rrt(607) * density(49) 
  pd(41,49) = pd(41,49) + rrt(607) * density(33) 
  pd(49,33) = pd(49,33) - rrt(607) * density(49) 
  pd(49,49) = pd(49,49) - rrt(607) * density(33) 
  pd(21,34) = pd(21,34) + rrt(608) * density(49) 
  pd(21,49) = pd(21,49) + rrt(608) * density(34) 
  pd(34,34) = pd(34,34) - rrt(608) * density(49) 
  pd(34,49) = pd(34,49) - rrt(608) * density(34) 
  pd(41,34) = pd(41,34) + rrt(608) * density(49) 
  pd(41,49) = pd(41,49) + rrt(608) * density(34) 
  pd(49,34) = pd(49,34) - rrt(608) * density(49) 
  pd(49,49) = pd(49,49) - rrt(608) * density(34) 
  pd(40,45) = pd(40,45) + rrt(609) * density(49) 
  pd(40,49) = pd(40,49) + rrt(609) * density(45) 
  pd(41,45) = pd(41,45) + rrt(609) * density(49) 
  pd(41,49) = pd(41,49) + rrt(609) * density(45) 
  pd(45,45) = pd(45,45) - rrt(609) * density(49) 
  pd(45,49) = pd(45,49) - rrt(609) * density(45) 
  pd(49,45) = pd(49,45) - rrt(609) * density(49) 
  pd(49,49) = pd(49,49) - rrt(609) * density(45) 
  pd(41,46) = pd(41,46) + rrt(610) * density(49) * 2.0d0
  pd(41,49) = pd(41,49) + rrt(610) * density(46) * 2.0d0
  pd(46,46) = pd(46,46) - rrt(610) * density(49) 
  pd(46,49) = pd(46,49) - rrt(610) * density(46) 
  pd(49,46) = pd(49,46) - rrt(610) * density(49) 
  pd(49,49) = pd(49,49) - rrt(610) * density(46) 
  pd(41,47) = pd(41,47) + rrt(611) * density(49) 
  pd(41,49) = pd(41,49) + rrt(611) * density(47) 
  pd(42,47) = pd(42,47) + rrt(611) * density(49) 
  pd(42,49) = pd(42,49) + rrt(611) * density(47) 
  pd(47,47) = pd(47,47) - rrt(611) * density(49) 
  pd(47,49) = pd(47,49) - rrt(611) * density(47) 
  pd(49,47) = pd(49,47) - rrt(611) * density(49) 
  pd(49,49) = pd(49,49) - rrt(611) * density(47) 
  pd(14,17) = pd(14,17) + rrt(612) * density(50) 
  pd(14,50) = pd(14,50) + rrt(612) * density(17) 
  pd(17,17) = pd(17,17) - rrt(612) * density(50) 
  pd(17,50) = pd(17,50) - rrt(612) * density(17) 
  pd(42,17) = pd(42,17) + rrt(612) * density(50) 
  pd(42,50) = pd(42,50) + rrt(612) * density(17) 
  pd(50,17) = pd(50,17) - rrt(612) * density(50) 
  pd(50,50) = pd(50,50) - rrt(612) * density(17) 
  pd(01,18) = pd(01,18) + rrt(613) * density(50) 
  pd(01,50) = pd(01,50) + rrt(613) * density(18) 
  pd(18,18) = pd(18,18) - rrt(613) * density(50) 
  pd(18,50) = pd(18,50) - rrt(613) * density(18) 
  pd(42,18) = pd(42,18) + rrt(613) * density(50) 
  pd(42,50) = pd(42,50) + rrt(613) * density(18) 
  pd(50,18) = pd(50,18) - rrt(613) * density(50) 
  pd(50,50) = pd(50,50) - rrt(613) * density(18) 
  pd(29,33) = pd(29,33) + rrt(614) * density(50) 
  pd(29,50) = pd(29,50) + rrt(614) * density(33) 
  pd(33,33) = pd(33,33) - rrt(614) * density(50) 
  pd(33,50) = pd(33,50) - rrt(614) * density(33) 
  pd(42,33) = pd(42,33) + rrt(614) * density(50) 
  pd(42,50) = pd(42,50) + rrt(614) * density(33) 
  pd(50,33) = pd(50,33) - rrt(614) * density(50) 
  pd(50,50) = pd(50,50) - rrt(614) * density(33) 
  pd(21,34) = pd(21,34) + rrt(615) * density(50) 
  pd(21,50) = pd(21,50) + rrt(615) * density(34) 
  pd(34,34) = pd(34,34) - rrt(615) * density(50) 
  pd(34,50) = pd(34,50) - rrt(615) * density(34) 
  pd(42,34) = pd(42,34) + rrt(615) * density(50) 
  pd(42,50) = pd(42,50) + rrt(615) * density(34) 
  pd(50,34) = pd(50,34) - rrt(615) * density(50) 
  pd(50,50) = pd(50,50) - rrt(615) * density(34) 
  pd(40,45) = pd(40,45) + rrt(616) * density(50) 
  pd(40,50) = pd(40,50) + rrt(616) * density(45) 
  pd(42,45) = pd(42,45) + rrt(616) * density(50) 
  pd(42,50) = pd(42,50) + rrt(616) * density(45) 
  pd(45,45) = pd(45,45) - rrt(616) * density(50) 
  pd(45,50) = pd(45,50) - rrt(616) * density(45) 
  pd(50,45) = pd(50,45) - rrt(616) * density(50) 
  pd(50,50) = pd(50,50) - rrt(616) * density(45) 
  pd(41,46) = pd(41,46) + rrt(617) * density(50) 
  pd(41,50) = pd(41,50) + rrt(617) * density(46) 
  pd(42,46) = pd(42,46) + rrt(617) * density(50) 
  pd(42,50) = pd(42,50) + rrt(617) * density(46) 
  pd(46,46) = pd(46,46) - rrt(617) * density(50) 
  pd(46,50) = pd(46,50) - rrt(617) * density(46) 
  pd(50,46) = pd(50,46) - rrt(617) * density(50) 
  pd(50,50) = pd(50,50) - rrt(617) * density(46) 
  pd(42,47) = pd(42,47) + rrt(618) * density(50) * 2.0d0
  pd(42,50) = pd(42,50) + rrt(618) * density(47) * 2.0d0
  pd(47,47) = pd(47,47) - rrt(618) * density(50) 
  pd(47,50) = pd(47,50) - rrt(618) * density(47) 
  pd(50,47) = pd(50,47) - rrt(618) * density(50) 
  pd(50,50) = pd(50,50) - rrt(618) * density(47) 
  pd(14,17) = pd(14,17) + rrt(619) * density(51) 
  pd(14,51) = pd(14,51) + rrt(619) * density(17) 
  pd(17,17) = pd(17,17) - rrt(619) * density(51) 
  pd(17,51) = pd(17,51) - rrt(619) * density(17) 
  pd(43,17) = pd(43,17) + rrt(619) * density(51) 
  pd(43,51) = pd(43,51) + rrt(619) * density(17) 
  pd(51,17) = pd(51,17) - rrt(619) * density(51) 
  pd(51,51) = pd(51,51) - rrt(619) * density(17) 
  pd(01,18) = pd(01,18) + rrt(620) * density(51) 
  pd(01,51) = pd(01,51) + rrt(620) * density(18) 
  pd(18,18) = pd(18,18) - rrt(620) * density(51) 
  pd(18,51) = pd(18,51) - rrt(620) * density(18) 
  pd(43,18) = pd(43,18) + rrt(620) * density(51) 
  pd(43,51) = pd(43,51) + rrt(620) * density(18) 
  pd(51,18) = pd(51,18) - rrt(620) * density(51) 
  pd(51,51) = pd(51,51) - rrt(620) * density(18) 
  pd(29,33) = pd(29,33) + rrt(621) * density(51) 
  pd(29,51) = pd(29,51) + rrt(621) * density(33) 
  pd(33,33) = pd(33,33) - rrt(621) * density(51) 
  pd(33,51) = pd(33,51) - rrt(621) * density(33) 
  pd(43,33) = pd(43,33) + rrt(621) * density(51) 
  pd(43,51) = pd(43,51) + rrt(621) * density(33) 
  pd(51,33) = pd(51,33) - rrt(621) * density(51) 
  pd(51,51) = pd(51,51) - rrt(621) * density(33) 
  pd(21,34) = pd(21,34) + rrt(622) * density(51) 
  pd(21,51) = pd(21,51) + rrt(622) * density(34) 
  pd(34,34) = pd(34,34) - rrt(622) * density(51) 
  pd(34,51) = pd(34,51) - rrt(622) * density(34) 
  pd(43,34) = pd(43,34) + rrt(622) * density(51) 
  pd(43,51) = pd(43,51) + rrt(622) * density(34) 
  pd(51,34) = pd(51,34) - rrt(622) * density(51) 
  pd(51,51) = pd(51,51) - rrt(622) * density(34) 
  pd(40,45) = pd(40,45) + rrt(623) * density(51) 
  pd(40,51) = pd(40,51) + rrt(623) * density(45) 
  pd(43,45) = pd(43,45) + rrt(623) * density(51) 
  pd(43,51) = pd(43,51) + rrt(623) * density(45) 
  pd(45,45) = pd(45,45) - rrt(623) * density(51) 
  pd(45,51) = pd(45,51) - rrt(623) * density(45) 
  pd(51,45) = pd(51,45) - rrt(623) * density(51) 
  pd(51,51) = pd(51,51) - rrt(623) * density(45) 
  pd(41,46) = pd(41,46) + rrt(624) * density(51) 
  pd(41,51) = pd(41,51) + rrt(624) * density(46) 
  pd(43,46) = pd(43,46) + rrt(624) * density(51) 
  pd(43,51) = pd(43,51) + rrt(624) * density(46) 
  pd(46,46) = pd(46,46) - rrt(624) * density(51) 
  pd(46,51) = pd(46,51) - rrt(624) * density(46) 
  pd(51,46) = pd(51,46) - rrt(624) * density(51) 
  pd(51,51) = pd(51,51) - rrt(624) * density(46) 
  pd(42,47) = pd(42,47) + rrt(625) * density(51) 
  pd(42,51) = pd(42,51) + rrt(625) * density(47) 
  pd(43,47) = pd(43,47) + rrt(625) * density(51) 
  pd(43,51) = pd(43,51) + rrt(625) * density(47) 
  pd(47,47) = pd(47,47) - rrt(625) * density(51) 
  pd(47,51) = pd(47,51) - rrt(625) * density(47) 
  pd(51,47) = pd(51,47) - rrt(625) * density(51) 
  pd(51,51) = pd(51,51) - rrt(625) * density(47) 
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(54,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(54,:) = pd(54,:) * ZDPlasKin_cfg(13)
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
  rrt(001) = KVT10_N2N2*1.0D0
  rrt(002) = KVT10_N2N2*2.0D0
  rrt(003) = KVT10_N2N2*3.0D0
  rrt(004) = KVT10_N2N2*4.0D0
  rrt(005) = KVT10_N2N2*5.0D0
  rrt(006) = KVT10_N2N2*6.0D0
  rrt(007) = KVT10_N2N2*7.0D0
  rrt(008) = KVT10_N2N2*8.0D0
  rrt(009) = KVT01_N2N2*1.0D0
  rrt(010) = KVT01_N2N2*2.0D0
  rrt(011) = KVT01_N2N2*3.0D0
  rrt(012) = KVT01_N2N2*4.0D0
  rrt(013) = KVT01_N2N2*5.0D0
  rrt(014) = KVT01_N2N2*6.0D0
  rrt(015) = KVT01_N2N2*7.0D0
  rrt(016) = KVT01_N2N2*8.0D0
  rrt(017) = KVT10_N2N*1.0D0
  rrt(018) = KVT10_N2N*2.0D0
  rrt(019) = KVT10_N2N*3.0D0
  rrt(020) = KVT10_N2N*4.0D0
  rrt(021) = KVT10_N2N*5.0D0
  rrt(022) = KVT10_N2N*6.0D0
  rrt(023) = KVT10_N2N*7.0D0
  rrt(024) = KVT10_N2N*8.0D0
  rrt(025) = KVT01_N2N*1.0D0
  rrt(026) = KVT01_N2N*2.0D0
  rrt(027) = KVT01_N2N*3.0D0
  rrt(028) = KVT01_N2N*4.0D0
  rrt(029) = KVT01_N2N*5.0D0
  rrt(030) = KVT01_N2N*6.0D0
  rrt(031) = KVT01_N2N*7.0D0
  rrt(032) = KVT01_N2N*8.0D0
  rrt(033) = KVT10_N2O*1.0D0
  rrt(034) = KVT10_N2O*2.0D0
  rrt(035) = KVT10_N2O*3.0D0
  rrt(036) = KVT10_N2O*4.0D0
  rrt(037) = KVT10_N2O*5.0D0
  rrt(038) = KVT10_N2O*6.0D0
  rrt(039) = KVT10_N2O*7.0D0
  rrt(040) = KVT10_N2O*8.0D0
  rrt(041) = KVT01_N2O*1.0D0
  rrt(042) = KVT01_N2O*2.0D0
  rrt(043) = KVT01_N2O*3.0D0
  rrt(044) = KVT01_N2O*4.0D0
  rrt(045) = KVT01_N2O*5.0D0
  rrt(046) = KVT01_N2O*6.0D0
  rrt(047) = KVT01_N2O*7.0D0
  rrt(048) = KVT01_N2O*8.0D0
  rrt(049) = KVT10_O2O2*1.0D0
  rrt(050) = KVT10_O2O2*2.0D0
  rrt(051) = KVT10_O2O2*3.0D0
  rrt(052) = KVT10_O2O2*4.0D0
  rrt(053) = KVT01_O2O2*1.0D0
  rrt(054) = KVT01_O2O2*2.0D0
  rrt(055) = KVT01_O2O2*3.0D0
  rrt(056) = KVT01_O2O2*4.0D0
  rrt(057) = KVT10_O2O*1.0D0
  rrt(058) = KVT10_O2O*2.0D0
  rrt(059) = KVT10_O2O*3.0D0
  rrt(060) = KVT10_O2O*4.0D0
  rrt(061) = KVT01_O2O*1.0D0
  rrt(062) = KVT01_O2O*2.0D0
  rrt(063) = KVT01_O2O*3.0D0
  rrt(064) = KVT01_O2O*4.0D0
  rrt(065) = bolsig_rates(bolsig_pointer(1))
  rrt(066) = bolsig_rates(bolsig_pointer(2))
  rrt(067) = bolsig_rates(bolsig_pointer(3))
  rrt(068) = bolsig_rates(bolsig_pointer(4))
  rrt(069) = bolsig_rates(bolsig_pointer(5))
  rrt(070) = bolsig_rates(bolsig_pointer(6))
  rrt(071) = bolsig_rates(bolsig_pointer(7))
  rrt(072) = bolsig_rates(bolsig_pointer(8))
  rrt(073) = bolsig_rates(bolsig_pointer(9))
  rrt(074) = bolsig_rates(bolsig_pointer(10))
  rrt(075) = bolsig_rates(bolsig_pointer(11))
  rrt(076) = bolsig_rates(bolsig_pointer(12))
  rrt(077) = bolsig_rates(bolsig_pointer(13))
  rrt(078) = bolsig_rates(bolsig_pointer(14))
  rrt(079) = bolsig_rates(bolsig_pointer(15))
  rrt(080) = bolsig_rates(bolsig_pointer(16))
  rrt(081) = bolsig_rates(bolsig_pointer(17))
  rrt(082) = bolsig_rates(bolsig_pointer(18))
  rrt(083) = bolsig_rates(bolsig_pointer(19))
  rrt(084) = bolsig_rates(bolsig_pointer(20))
  rrt(085) = bolsig_rates(bolsig_pointer(21))
  rrt(086) = bolsig_rates(bolsig_pointer(22))
  rrt(087) = bolsig_rates(bolsig_pointer(23))
  rrt(088) = bolsig_rates(bolsig_pointer(24))
  rrt(089) = bolsig_rates(bolsig_pointer(25))
  rrt(090) = bolsig_rates(bolsig_pointer(26))
  rrt(091) = bolsig_rates(bolsig_pointer(27))
  rrt(092) = 1.8D-7*(300.0D0/TE)**0.39*0.50D0
  rrt(093) = 1.8D-7*(300.0D0/TE)**0.39*0.45D0
  rrt(094) = 1.8D-7*(300.0D0/TE)**0.39*0.05D0
  rrt(095) = 2.7D-7*(300.0D0/TE)**0.7*0.55D0
  rrt(096) = 2.7D-7*(300.0D0/TE)**0.7*0.40D0
  rrt(097) = 2.7D-7*(300.0D0/TE)**0.7*0.05D0
  rrt(098) = 4.2D-7*(300.0D0/TE)**0.85*0.20D0
  rrt(099) = 4.2D-7*(300.0D0/TE)**0.85*0.80D0
  rrt(100) = 2.0D-7*(300.0D0/TE)**0.5
  rrt(101) = 2.3D-6*(300.0D0/TE)**0.53
  rrt(102) = rrt(100)
  rrt(103) = rrt(100)
  rrt(104) = 1.4D-6*(300.0D0/TE)**0.5
  rrt(105) = 1.3D-6*(300.0D0/TE)**0.5
  rrt(106) = 7.0D-20*(300.0D0/TE)**4.5
  rrt(107) = rrt(106)
  rrt(108) = 6.0D-27*(300.0D0/TE)**1.5*ANY_NEUTRAL
  rrt(109) = rrt(108)
  rrt(110) = bolsig_rates(bolsig_pointer(28))
  rrt(111) = bolsig_rates(bolsig_pointer(29))
  rrt(112) = bolsig_rates(bolsig_pointer(30))
  rrt(113) = bolsig_rates(bolsig_pointer(31))
  rrt(114) = bolsig_rates(bolsig_pointer(32))
  rrt(115) = bolsig_rates(bolsig_pointer(33))
  rrt(116) = 1.0D-11
  rrt(117) = 1.0D-31
  rrt(118) = 1.0D-31
  rrt(119) = 1.0D-31*ANY_NEUTRAL
  rrt(120) = 8.0D-31*ANY_NEUTRAL
  rrt(121) = 6.0D-33*ANY_NEUTRAL
  rrt(122) = 1.1D-31*(300.0D0/TE)**2*EXP(-70.0D0/TGAS)*EXP(1500.0D0*(TE-TGAS)/(TE*TGAS))
  rrt(123) = 1.4D-10
  rrt(124) = 2.6D-10
  rrt(125) = 2.6D-10
  rrt(126) = 5.0D-13
  rrt(127) = 5.0D-15
  rrt(128) = 3.0D-10
  rrt(129) = 6.9D-10
  rrt(130) = 2.2D-9
  rrt(131) = 1.9D-9
  rrt(132) = 3.0D-10
  rrt(133) = 1.5D-10
  rrt(134) = 5.0D-10
  rrt(135) = 2.7D-10*(TEFFN2/300.0D0)**0.5*EXP(-5590.0D0/TEFFN2)
  rrt(136) = 2.0D-10
  rrt(137) = 3.6D-10
  rrt(138) = 1.9D-12*(TEFFN2/300.0D0)**0.5*EXP(-4990.0D0/TEFFN2)
  rrt(139) = 2.1D-9
  rrt(140) = 2.5D-9
  rrt(141) = 3.0D-10
  rrt(142) = 5.0D-10
  rrt(143) = 5.0D-10
  rrt(144) = 5.0D-10
  rrt(145) = 5.0D-10
  rrt(146) = 5.0D-10
  rrt(147) = 1.5D-10
  rrt(148) = 1.5D-10
  rrt(149) = 1.5D-10
  rrt(150) = 1.5D-10
  rrt(151) = 2.1D-9
  rrt(152) = 2.1D-9
  rrt(153) = 2.1D-9
  rrt(154) = 2.1D-9
  rrt(155) = 2.1D-9
  rrt(156) = 2.5D-9
  rrt(157) = 2.5D-9
  rrt(158) = 2.5D-9
  rrt(159) = 2.5D-9
  rrt(160) = 2.5D-9
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
  rrt(176) = 2.0D-14*(TGAS/300.0D0)**0.55
  rrt(177) = 3.0D-16
  rrt(178) = 6.9D-11
  rrt(179) = 1.0D-11
  rrt(180) = 1.0D-12
  rrt(181) = 3.0D-10
  rrt(182) = 1.5D-10
  rrt(183) = 3.0D-11
  rrt(184) = 2.0D-12
  rrt(185) = 3.0D-10
  rrt(186) = 2.4D-10
  rrt(187) = 1.0D-11
  rrt(188) = 3.0D-10
  rrt(189) = 1.9D-13
  rrt(190) = 2.8D-11
  rrt(191) = 3.6D-10
  rrt(192) = 4.0D-12
  rrt(193) = 1.0D-11
  rrt(194) = 1.7D-33
  rrt(195) = 1.7D-33
  rrt(196) = 1.7D-33
  rrt(197) = 1.0D-32
  rrt(198) = 1.0D-32
  rrt(199) = 2.4D-33
  rrt(200) = 2.4D-33
  rrt(201) = 2.4D-33
  rrt(202) = 1.4D-32
  rrt(203) = 1.4D-32
  rrt(204) = 4.0D-13
  rrt(205) = 5.2D-12
  rrt(206) = 1.8D-10
  rrt(207) = 3.5D-12
  rrt(208) = 1.0D-13*EXP(-510.0D0/TGAS)
  rrt(209) = 1.8D-12
  rrt(210) = 1.0D-12
  rrt(211) = 6.0D-13
  rrt(212) = 6.0D-14
  rrt(213) = 1.0D-13
  rrt(214) = 2.6D-12
  rrt(215) = 3.0D-11
  rrt(216) = 7.0D-16
  rrt(217) = 2.0D-14*EXP(-600.0D0/TGAS)
  rrt(218) = 3.8D-18*EXP(-205.0D0/TGAS)
  rrt(219) = 3.0D-21
  rrt(220) = 2.5D-11
  rrt(221) = 5.2D-11*EXP(-2840.0D0/TGAS)
  rrt(222) = 7.0D-28*TGAS**3.8*EXP(700.0D0/TGAS)
  rrt(223) = 1.0D-11*EXP(-2300.0D0/TGAS)
  rrt(224) = 8.1D-14
  rrt(225) = 3.4D-11*(300.0D0/TGAS)**0.1*EXP(-4200.0D0/TGAS)
  rrt(226) = 4.3D-22*TGAS**2.4*EXP(-281.0D0/TGAS)
  rrt(227) = 1.7D-15*(TGAS/300.0D0)
  rrt(228) = 6.0D-14
  rrt(229) = 2.2D-11
  rrt(230) = 9.0D-12
  rrt(231) = 3.0D-13
  rrt(232) = 9.0D-15
  rrt(233) = 8.0D-12
  rrt(234) = 6.4D-12*EXP(67.0D0/TGAS)
  rrt(235) = 1.0D-12
  rrt(236) = 2.6D-11*EXP(67.0D0/TGAS)
  rrt(237) = 2.3D-11
  rrt(238) = 1.2D-10
  rrt(239) = 1.2D-10
  rrt(240) = 1.7D-10
  rrt(241) = 7.2D-11
  rrt(242) = 4.4D-11
  rrt(243) = 5.0D-11*EXP(-300.0D0/TGAS)
  rrt(244) = 1.0D-12
  rrt(245) = 1.3D-12*EXP(-850.0D0/TGAS)
  rrt(246) = 3.0D-12*EXP(-850.0D0/TGAS)
  rrt(247) = 1.0D-17
  rrt(248) = 1.1D-10
  rrt(249) = 2.9D-11
  rrt(250) = 3.2D-11
  rrt(251) = 2.9D-10
  rrt(252) = 5.1D-10
  rrt(253) = 2.9D-10
  rrt(254) = 2.9D-10
  rrt(255) = 6.3D-12
  rrt(256) = 3.1D-12
  rrt(257) = 1.8D-11*(TGAS/300.0)**0.5
  rrt(258) = 3.2D-12*(TGAS/300.0)*EXP(-3150.0D0/TGAS)
  rrt(259) = 9.1D-13
  rrt(260) = 3.0D-12
  rrt(261) = 7.0D-13
  rrt(262) = 2.3D-12
  rrt(263) = 3.0D-10*EXP(-38370.0D0/TGAS)
  rrt(264) = 7.5D-12*(TGAS/300.0)*EXP(-19500.0D0/TGAS)
  rrt(265) = 4.2D-18
  rrt(266) = 8.3D-12*EXP(-14000.0D0/TGAS)
  rrt(267) = 1.5D-10*EXP(-14090.0D0/TGAS)
  rrt(268) = 9.1D-12*(TGAS/300.0D0)**0.18
  rrt(269) = 1.0D-11
  rrt(270) = 2.5D-10*EXP(-50390.0D0/TGAS)
  rrt(271) = 3.3D-16*(300.0D0/TGAS)**0.5*EXP(-39200.0D0/TGAS)
  rrt(272) = 2.2D-12*EXP(-32100.0D0/TGAS)
  rrt(273) = 5.1D-13*EXP(-33660.0D0/TGAS)
  rrt(274) = 2.8D-12*EXP(-23400.0D0/TGAS)
  rrt(275) = 2.5D-13*EXP(-765.0D0/TGAS)
  rrt(276) = 4.6D-10*EXP(-25170.0D0/TGAS)
  rrt(277) = 1.7D-11
  rrt(278) = 2.0D-11*EXP(-49800.0D0/TGAS)
  rrt(279) = 2.8D-12*EXP(-25400.0D0/TGAS)
  rrt(280) = 3.3D-12*EXP(-13500.0D0/TGAS)
  rrt(281) = 4.5D-10*EXP(-18500.0D0/TGAS)
  rrt(282) = 1.2D-13*EXP(-2450.0D0/TGAS)
  rrt(283) = 2.3D-13*EXP(-1600.0D0/TGAS)
  rrt(284) = 1.5D-12*EXP(-15020.0D0/TGAS)
  rrt(285) = 4.3D-12*EXP(-3850.0D0/TGAS)
  rrt(286) = 2.7D-11*EXP(-6.74D4/TGAS)
  rrt(287) = 1.6D-12*(TGAS/300.0D0)**0.5*(0.19D0+8.6D0*TGAS)*EXP(-32000.0D0/TGAS)
  rrt(288) = 5.4D-8*(1.0D0-EXP(-3354.0D0/TGAS))*EXP(-113200.0D0/TGAS)*1.0D0
  rrt(289) = rrt(288)
  rrt(290) = rrt(288)
  rrt(291) = 5.4D-8*(1.0D0-EXP(-3354.0D0/TGAS))*EXP(-113200.0D0/TGAS)*6.6D0
  rrt(292) = rrt(291)
  rrt(293) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*1.0D0
  rrt(294) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*5.9D0
  rrt(295) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*21.D0
  rrt(296) = rrt(293)
  rrt(297) = rrt(293)
  rrt(298) = 8.7D-9*EXP(-75994.0D0/TGAS)*1.0D0
  rrt(299) = rrt(298)
  rrt(300) = 8.7D-9*EXP(-75994.0D0/TGAS)*20.D0
  rrt(301) = rrt(300)
  rrt(302) = rrt(300)
  rrt(303) = 6.6D-10*EXP(-11600.0D0/TGAS)*1.0D0
  rrt(304) = 6.6D-10*EXP(-11600.0D0/TGAS)*0.38D0
  rrt(305) = 6.6D-10*EXP(-11600.0D0/TGAS)*6.3D0*EXP(170.0D0/TGAS)
  rrt(306) = rrt(305)
  rrt(307) = 1.2D-8*(300.0D0/TGAS)*EXP(-29000.0D0/TGAS)*1.0D0
  rrt(308) = rrt(307)
  rrt(309) = 1.2D-8*(300.0D0/TGAS)*EXP(-29000.0D0/TGAS)*2.0D0
  rrt(310) = 1.2D-8*(300.0D0/TGAS)*EXP(-29000.0D0/TGAS)*4.0D0
  rrt(311) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*1.0D0
  rrt(312) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*0.78D0
  rrt(313) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*7.8D0
  rrt(314) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*5.9D0
  rrt(315) = 3.1D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*1.0D0
  rrt(316) = rrt(315)
  rrt(317) = rrt(315)
  rrt(318) = 3.1D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*10.D0
  rrt(319) = rrt(318)
  rrt(320) = 6.2D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*1.0D0
  rrt(321) = rrt(320)
  rrt(322) = rrt(320)
  rrt(323) = 6.2D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*12.D0
  rrt(324) = rrt(323)
  rrt(325) = 2.1D-11*(300.0D0/TGAS)**4.4*EXP(-11080.0D0/TGAS)*ANY_NEUTRAL
  rrt(326) = MAX(8.3D-34*EXP(500.0D0/TGAS),1.91D-33)
  rrt(327) = 1.8D-33*EXP(435.0D0/TGAS)*1.0D0
  rrt(328) = rrt(327)
  rrt(329) = 1.8D-33*EXP(435.0D0/TGAS)*3.0D0
  rrt(330) = rrt(329)
  rrt(331) = MAX(2.8D-34*EXP(720.0D0/TGAS),1.0D-33*(300.0D0/TGAS)**0.41)
  rrt(332) = 4.0D-33*(300.0D0/TGAS)**0.41*1.0D0
  rrt(333) = 4.0D-33*(300.0D0/TGAS)**0.41*0.8D0
  rrt(334) = 4.0D-33*(300.0D0/TGAS)**0.41*3.6D0
  rrt(335) = 4.0D-33*(300.0D0/TGAS)**0.41*0.17D0
  rrt(336) = 1.0D-32*(300.0D0/TGAS)**0.5
  rrt(337) = rrt(336)
  rrt(338) = 1.8D-31*(300.0D0/TGAS)
  rrt(339) = rrt(338)
  rrt(340) = rrt(338)
  rrt(341) = MAX(5.8D-34*(300.0D0/TGAS)**2.8,5.4D-34*(300.0D0/TGAS)**1.9)
  rrt(342) = 7.6D-34*(300.0D0/TGAS)**1.9
  rrt(343) = rrt(342)
  rrt(344) = MIN(3.9D-33*(300.0D0/TGAS)**1.9,1.1D-34*EXP(1060.0D0/TGAS))
  rrt(345) = rrt(344)
  rrt(346) = 3.9D-35*EXP(-10400.0D0/TGAS)*ANY_NEUTRAL
  rrt(347) = 1.2D-31*(300.0D0/TGAS)**1.8*1.0D0
  rrt(348) = 1.2D-31*(300.0D0/TGAS)**1.8*0.78D0
  rrt(349) = rrt(348)
  rrt(350) = 8.9D-32*(300.0D0/TGAS)**2*1.0D0
  rrt(351) = rrt(350)
  rrt(352) = 8.9D-32*(300.0D0/TGAS)**2*13.D0
  rrt(353) = rrt(352)
  rrt(354) = 8.9D-32*(300.0D0/TGAS)**2*2.4D0
  rrt(355) = 3.7D-30*(300.0D0/TGAS)**4.1*ANY_NEUTRAL
  rrt(356) = 1.0D-12
  rrt(357) = 2.8D-10
  rrt(358) = 2.5D-10
  rrt(359) = 2.8D-11
  rrt(360) = 5.0D-10
  rrt(361) = 8.0D-10
  rrt(362) = 3.0D-12
  rrt(363) = 1.0D-12
  rrt(364) = 5.5D-10
  rrt(365) = (1.5D0-2.0D-3*TEFFN+9.6D-7*TEFFN**2)*1.0D-12
  rrt(366) = 2.0D-11*(300.0D0/TEFFN)**0.5
  rrt(367) = 1.0D-10
  rrt(368) = 2.4D-11
  rrt(369) = 3.0D-12
  rrt(370) = 1.3D-10
  rrt(371) = 2.3D-10
  rrt(372) = 2.2D-10
  rrt(373) = 2.0D-11
  rrt(374) = 1.6D-9
  rrt(375) = 6.0D-11*(300.0D0/TEFFN2)**0.5
  rrt(376) = 1.3D-10*(300.0D0/TEFFN2)**0.5
  rrt(377) = 1.0D-10
  rrt(378) = 7.2D-13*(TEFFN2/300.0D0)
  rrt(379) = 3.3D-10
  rrt(380) = 5.0D-10
  rrt(381) = 4.0D-10
  rrt(382) = 1.0D-17
  rrt(383) = 1.2D-10
  rrt(384) = 6.3D-10
  rrt(385) = 1.0D-11
  rrt(386) = 6.6D-10
  rrt(387) = 2.3D-11
  rrt(388) = 4.4D-11
  rrt(389) = 6.6D-11
  rrt(390) = 7.0D-11
  rrt(391) = 7.0D-11
  rrt(392) = 2.9D-10
  rrt(393) = 2.9D-10
  rrt(394) = MIN(2.1D-16*EXP(TEFFN4/121.0D0),1.0D-10)
  rrt(395) = 2.5D-10
  rrt(396) = 2.5D-10
  rrt(397) = 1.0D-11
  rrt(398) = 4.0D-10
  rrt(399) = 4.6D-12*(TEFFN4/300.0D0)**2.5*EXP(-2650.0D0/TEFFN4)
  rrt(400) = 3.3D-6*(300.0D0/TEFFN4)**4*EXP(-5030.0D0/TEFFN4)
  rrt(401) = 1.0D-10
  rrt(402) = 1.0D-10
  rrt(403) = 3.0D-10
  rrt(404) = 1.0D-10
  rrt(405) = 1.1D-6*(300.0D0/TEFFN4)**5.3*EXP(-2360.0D0/TEFFN4)
  rrt(406) = 1.0D-9
  rrt(407) = 1.7D-29*(300.0D0/TEFFN)**2.1
  rrt(408) = 1.0D-29*ANY_NEUTRAL
  rrt(409) = rrt(408)
  rrt(410) = 6.0D-29*(300.0D0/TEFFN)**2*ANY_NEUTRAL
  rrt(411) = rrt(408)
  rrt(412) = rrt(408)
  rrt(413) = 5.2D-29*(300.0D0/TEFFN2)**2.2
  rrt(414) = 9.0D-30*EXP(400.0D0/TEFFN2)
  rrt(415) = 2.4D-30*(300.0D0/TEFFN2)**3.2
  rrt(416) = 9.0D-31*(300.0D0/TEFFN2)**2
  rrt(417) = 1.0D-10
  rrt(418) = 8.0D-10
  rrt(419) = 1.2D-9
  rrt(420) = 2.0D-10
  rrt(421) = 2.0D-12
  rrt(422) = 3.3D-10
  rrt(423) = 3.5D-10
  rrt(424) = 7.0D-10
  rrt(425) = 5.0D-10
  rrt(426) = 1.0D-11
  rrt(427) = 1.0D-11
  rrt(428) = 2.6D-12
  rrt(429) = 7.0D-11
  rrt(430) = 2.0D-11
  rrt(431) = 5.0D-10
  rrt(432) = 5.0D-10
  rrt(433) = 7.4D-10
  rrt(434) = 2.8D-14
  rrt(435) = 1.8D-11
  rrt(436) = 4.0D-12
  rrt(437) = 5.0D-10
  rrt(438) = 7.0D-10
  rrt(439) = 3.0D-15
  rrt(440) = 1.0D-10*EXP(-1044.0D0/TEFFN4)
  rrt(441) = rrt(440)
  rrt(442) = 4.0D-10
  rrt(443) = 3.0D-10
  rrt(444) = 1.0D-10
  rrt(445) = 1.0D-10
  rrt(446) = 2.5D-10
  rrt(447) = 1.1D-30*(300.0D0/TEFFN)*ANY_NEUTRAL
  rrt(448) = rrt(408)
  rrt(449) = 3.5D-31*(300.0D0/TEFFN2)*ANY_NEUTRAL
  rrt(450) = 2.0D-7*(300.0D0/TIONN)**0.5
  rrt(451) = rrt(450)
  rrt(452) = rrt(450)
  rrt(453) = rrt(450)
  rrt(454) = rrt(450)
  rrt(455) = rrt(450)
  rrt(456) = rrt(450)
  rrt(457) = rrt(450)
  rrt(458) = rrt(450)
  rrt(459) = rrt(450)
  rrt(460) = rrt(450)
  rrt(461) = rrt(450)
  rrt(462) = rrt(450)
  rrt(463) = rrt(450)
  rrt(464) = rrt(450)
  rrt(465) = rrt(450)
  rrt(466) = rrt(450)
  rrt(467) = rrt(450)
  rrt(468) = rrt(450)
  rrt(469) = rrt(450)
  rrt(470) = rrt(450)
  rrt(471) = rrt(450)
  rrt(472) = rrt(450)
  rrt(473) = rrt(450)
  rrt(474) = rrt(450)
  rrt(475) = rrt(450)
  rrt(476) = rrt(450)
  rrt(477) = rrt(450)
  rrt(478) = rrt(450)
  rrt(479) = rrt(450)
  rrt(480) = rrt(450)
  rrt(481) = rrt(450)
  rrt(482) = rrt(450)
  rrt(483) = rrt(450)
  rrt(484) = rrt(450)
  rrt(485) = rrt(450)
  rrt(486) = rrt(450)
  rrt(487) = rrt(450)
  rrt(488) = rrt(450)
  rrt(489) = rrt(450)
  rrt(490) = rrt(450)
  rrt(491) = rrt(450)
  rrt(492) = rrt(450)
  rrt(493) = rrt(450)
  rrt(494) = rrt(450)
  rrt(495) = rrt(450)
  rrt(496) = rrt(450)
  rrt(497) = rrt(450)
  rrt(498) = rrt(450)
  rrt(499) = 1.0D-7
  rrt(500) = 1.0D-7
  rrt(501) = 1.0D-7
  rrt(502) = 1.0D-7
  rrt(503) = 1.0D-7
  rrt(504) = 1.0D-7
  rrt(505) = 1.0D-7
  rrt(506) = 1.0D-7
  rrt(507) = 1.0D-7
  rrt(508) = 1.0D-7
  rrt(509) = 1.0D-7
  rrt(510) = 1.0D-7
  rrt(511) = 1.0D-7
  rrt(512) = 1.0D-7
  rrt(513) = 1.0D-7
  rrt(514) = 1.0D-7
  rrt(515) = 1.0D-7
  rrt(516) = 1.0D-7
  rrt(517) = 1.0D-7
  rrt(518) = 1.0D-7
  rrt(519) = 1.0D-7
  rrt(520) = 1.0D-7
  rrt(521) = 1.0D-7
  rrt(522) = 1.0D-7
  rrt(523) = 1.0D-7
  rrt(524) = 1.0D-7
  rrt(525) = 1.0D-7
  rrt(526) = 1.0D-7
  rrt(527) = 1.0D-7
  rrt(528) = 1.0D-7
  rrt(529) = 1.0D-7
  rrt(530) = 1.0D-7
  rrt(531) = 1.0D-7
  rrt(532) = 1.0D-7
  rrt(533) = 1.0D-7
  rrt(534) = 1.0D-7
  rrt(535) = 1.0D-7
  rrt(536) = 1.0D-7
  rrt(537) = 1.0D-7
  rrt(538) = 1.0D-7
  rrt(539) = 1.0D-7
  rrt(540) = 1.0D-7
  rrt(541) = 1.0D-7
  rrt(542) = 1.0D-7
  rrt(543) = 1.0D-7
  rrt(544) = 1.0D-7
  rrt(545) = 1.0D-7
  rrt(546) = 1.0D-7
  rrt(547) = 1.0D-7
  rrt(548) = 1.0D-7
  rrt(549) = 1.0D-7
  rrt(550) = 1.0D-7
  rrt(551) = 1.0D-7
  rrt(552) = 1.0D-7
  rrt(553) = 1.0D-7
  rrt(554) = 1.0D-7
  rrt(555) = 1.0D-7
  rrt(556) = 1.0D-7
  rrt(557) = 1.0D-7
  rrt(558) = 1.0D-7
  rrt(559) = 1.0D-7
  rrt(560) = 1.0D-7
  rrt(561) = 1.0D-7
  rrt(562) = 1.0D-7
  rrt(563) = 1.0D-7
  rrt(564) = 1.0D-7
  rrt(565) = 1.0D-7
  rrt(566) = 1.0D-7
  rrt(567) = 1.0D-7
  rrt(568) = 1.0D-7
  rrt(569) = 1.0D-7
  rrt(570) = 1.0D-7
  rrt(571) = 1.0D-7
  rrt(572) = 1.0D-7
  rrt(573) = 2.0D-25*(300.0D0/TIONN)**2.5*ANY_NEUTRAL
  rrt(574) = rrt(573)
  rrt(575) = rrt(573)
  rrt(576) = rrt(573)
  rrt(577) = rrt(573)
  rrt(578) = rrt(573)
  rrt(579) = rrt(573)
  rrt(580) = rrt(573)
  rrt(581) = rrt(573)
  rrt(582) = rrt(573)
  rrt(583) = rrt(573)
  rrt(584) = rrt(573)
  rrt(585) = rrt(573)
  rrt(586) = rrt(573)
  rrt(587) = rrt(573)
  rrt(588) = rrt(573)
  rrt(589) = rrt(573)
  rrt(590) = rrt(573)
  rrt(591) = 2.0D-25*(300.0D0/TIONN2)**2.5*ANY_NEUTRAL
  rrt(592) = rrt(591)
  rrt(593) = rrt(591)
  rrt(594) = rrt(591)
  rrt(595) = rrt(591)
  rrt(596) = rrt(591)
  rrt(597) = rrt(591)
  rrt(598) = rrt(591)
  rrt(599) = rrt(591)
  rrt(600) = rrt(591)
  rrt(601) = rrt(591)
  rrt(602) = rrt(591)
  rrt(603) = rrt(591)
  rrt(604) = rrt(591)
  rrt(605) = rrt(591)
  rrt(606) = rrt(591)
  rrt(607) = rrt(591)
  rrt(608) = rrt(591)
  rrt(609) = rrt(591)
  rrt(610) = rrt(591)
  rrt(611) = rrt(591)
  rrt(612) = rrt(591)
  rrt(613) = rrt(591)
  rrt(614) = rrt(591)
  rrt(615) = rrt(591)
  rrt(616) = rrt(591)
  rrt(617) = rrt(591)
  rrt(618) = rrt(591)
  rrt(619) = rrt(591)
  rrt(620) = rrt(591)
  rrt(621) = rrt(591)
  rrt(622) = rrt(591)
  rrt(623) = rrt(591)
  rrt(624) = rrt(591)
  rrt(625) = rrt(591)
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
