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
! Tue Dec 10 14:05:27 2019
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
  integer, parameter :: species_max = 53, species_electrons = 53, species_length = 9, reactions_max = 646, reactions_length = 44
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
  integer, parameter, private               :: bolsig_species_max = 21, bolsig_species_length = 6, bolsig_rates_max = 54, &
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
  /"bolsig:N2->N2(V1RES)                        ","bolsig:N2->N2(V1)                           ",&
   "bolsig:N2->N2(V2)                           ","bolsig:N2->N2(V3)                           ",&
   "bolsig:N2->N2(V4)                           ","bolsig:N2->N2(V5)                           ",&
   "bolsig:N2->N2(V6)                           ","bolsig:N2->N2(V7)                           ",&
   "bolsig:N2->N2(V8)                           ","bolsig:O2->O2(V1RES)                        ",&
   "bolsig:O2->O2(V1)                           ","bolsig:O2->O2(V2RES)                        ",&
   "bolsig:O2->O2(V2)                           ","bolsig:O2->O2(V3)                           ",&
   "bolsig:O2->O2(V4)                           ","bolsig:O2(V1)->O2                           ",&
   "bolsig:O2(V2)->O2                           ","bolsig:O2(V3)->O2                           ",&
   "bolsig:O2(V4)->O2                           ","N2(V1)+N2=>N2+N2                            ",&
   "N2(V2)+N2=>N2(V1)+N2                        ","N2(V3)+N2=>N2(V2)+N2                        ",&
   "N2(V4)+N2=>N2(V3)+N2                        ","N2(V5)+N2=>N2(V4)+N2                        ",&
   "N2(V6)+N2=>N2(V5)+N2                        ","N2(V7)+N2=>N2(V6)+N2                        ",&
   "N2(V8)+N2=>N2(V7)+N2                        ","N2+N2=>N2(V1)+N2                            ",&
   "N2(V1)+N2=>N2(V2)+N2                        ","N2(V2)+N2=>N2(V3)+N2                        ",&
   "N2(V3)+N2=>N2(V4)+N2                        ","N2(V4)+N2=>N2(V5)+N2                        ",&
   "N2(V5)+N2=>N2(V6)+N2                        ","N2(V6)+N2=>N2(V7)+N2                        ",&
   "N2(V7)+N2=>N2(V8)+N2                        ","N2(V1)+N=>N2+N                              "/
  data reaction_sign(37:72) &
  /"N2(V2)+N=>N2(V1)+N                          ","N2(V3)+N=>N2(V2)+N                          ",&
   "N2(V4)+N=>N2(V3)+N                          ","N2(V5)+N=>N2(V4)+N                          ",&
   "N2(V6)+N=>N2(V5)+N                          ","N2(V7)+N=>N2(V6)+N                          ",&
   "N2(V8)+N=>N2(V7)+N                          ","N2+N=>N2(V1)+N                              ",&
   "N2(V1)+N=>N2(V2)+N                          ","N2(V2)+N=>N2(V3)+N                          ",&
   "N2(V3)+N=>N2(V4)+N                          ","N2(V4)+N=>N2(V5)+N                          ",&
   "N2(V5)+N=>N2(V6)+N                          ","N2(V6)+N=>N2(V7)+N                          ",&
   "N2(V7)+N=>N2(V8)+N                          ","N2(V1)+O=>N2+O                              ",&
   "N2(V2)+O=>N2(V1)+O                          ","N2(V3)+O=>N2(V2)+O                          ",&
   "N2(V4)+O=>N2(V3)+O                          ","N2(V5)+O=>N2(V4)+O                          ",&
   "N2(V6)+O=>N2(V5)+O                          ","N2(V7)+O=>N2(V6)+O                          ",&
   "N2(V8)+O=>N2(V7)+O                          ","N2+O=>N2(V1)+O                              ",&
   "N2(V1)+O=>N2(V2)+O                          ","N2(V2)+O=>N2(V3)+O                          ",&
   "N2(V3)+O=>N2(V4)+O                          ","N2(V4)+O=>N2(V5)+O                          ",&
   "N2(V5)+O=>N2(V6)+O                          ","N2(V6)+O=>N2(V7)+O                          ",&
   "N2(V7)+O=>N2(V8)+O                          ","O2(V1)+O2=>O2+O2                            ",&
   "O2(V2)+O2=>O2(V1)+O2                        ","O2(V3)+O2=>O2(V2)+O2                        ",&
   "O2(V4)+O2=>O2(V3)+O2                        ","O2+O2=>O2(V1)+O2                            "/
  data reaction_sign(73:108) &
  /"O2(V1)+O2=>O2(V2)+O2                        ","O2(V2)+O2=>O2(V3)+O2                        ",&
   "O2(V3)+O2=>O2(V4)+O2                        ","O2(V1)+O=>O2+O                              ",&
   "O2(V2)+O=>O2(V1)+O                          ","O2(V3)+O=>O2(V2)+O                          ",&
   "O2(V4)+O=>O2(V3)+O                          ","O2+O=>O2(V1)+O                              ",&
   "O2(V1)+O=>O2(V2)+O                          ","O2(V2)+O=>O2(V3)+O                          ",&
   "O2(V3)+O=>O2(V4)+O                          ","bolsig:N2->N2(A3)                           ",&
   "bolsig:N2->N2(A3,V5-9)                      ","bolsig:N2->N2(A3,V10-)                      ",&
   "bolsig:N2->N2(B3)                           ","bolsig:N2->N2(W3)                           ",&
   "bolsig:N2->N2(B'3)                          ","bolsig:N2->N2(A'1)                          ",&
   "bolsig:N2->N2(A1)                           ","bolsig:N2->N2(W1)                           ",&
   "bolsig:N2->N2(C3)                           ","bolsig:N2->N2(E3)                           ",&
   "bolsig:N2->N2(A''1)                         ","bolsig:O2->O2(A1)                           ",&
   "bolsig:O2->O2(B1)                           ","bolsig:O2->O2(4.5EV)                        ",&
   "bolsig:O2->O2(6.0EV)                        ","bolsig:O2->O2(8.4EV)                        ",&
   "bolsig:O2->O2(9.97EV)                       ","bolsig:O->O(1D)                             ",&
   "bolsig:O->O(1S)                             ","bolsig:O2(A1)->O2                           ",&
   "bolsig:N->N^+                               ","bolsig:O->O^+                               ",&
   "bolsig:N2->N2^+                             ","bolsig:N2(A3)->N2^+                         "/
  data reaction_sign(109:144) &
  /"bolsig:O2->O2^+                             ","bolsig:O2(A1)->O2^+                         ",&
   "bolsig:NO->NO^+                             ","bolsig:N2O->N2O^+                           ",&
   "E+N2^+=>N+N                                 ","E+N2^+=>N+N(2D)                             ",&
   "E+N2^+=>N+N(2P)                             ","E+O2^+=>O+O                                 ",&
   "E+O2^+=>O+O(1D)                             ","E+O2^+=>O+O(1S)                             ",&
   "E+NO^+=>O+N                                 ","E+NO^+=>O+N(2D)                             ",&
   "E+N3^+=>N2+N                                ","E+N4^+=>N2+N2                               ",&
   "E+N2O^+=>N2+O                               ","E+NO2^+=>NO+O                               ",&
   "E+O4^+=>O2+O2                               ","E+O2^+N2=>O2+N2                             ",&
   "E+N^++E=>N+E                                ","E+O^++E=>O+E                                ",&
   "E+N^++ANY_NEUTRAL=>N+ANY_NEUTRAL            ","E+O^++ANY_NEUTRAL=>O+ANY_NEUTRAL            ",&
   "bolsig:O2->O^-+O                            ","bolsig:NO->O^-+N                            ",&
   "bolsig:O3->O^-+O2                           ","bolsig:O3->O2^-+O                           ",&
   "bolsig:N2O->NO^-+N                          ","bolsig:O2+O2->O2^-+O2                       ",&
   "E+NO2=>O^-+NO                               ","E+O+O2=>O^-+O2                              ",&
   "E+O+O2=>O2^-+O                              ","E+O3+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL          ",&
   "E+NO+ANY_NEUTRAL=>NO^-+ANY_NEUTRAL          ","E+N2O+ANY_NEUTRAL=>N2O^-+ANY_NEUTRAL        ",&
   "E+O2+N2=>O2^-+N2                            ","O^-+O=>O2+E                                 "/
  data reaction_sign(145:180) &
  /"O^-+N=>NO+E                                 ","O^-+NO=>NO2+E                               ",&
   "O^-+N2=>N2O+E                               ","O^-+O2=>O3+E                                ",&
   "O^-+O2(A1)=>O3+E                            ","O^-+O2(B1)=>O+O2+E                          ",&
   "O^-+N2(A3)=>O+N2+E                          ","O^-+N2(B3)=>O+N2+E                          ",&
   "O^-+O3=>O2+O2+E                             ","O2^-+O=>O3+E                                ",&
   "O2^-+N=>NO2+E                               ","O2^-+O2=>O2+O2+E                            ",&
   "O2^-+O2(A1)=>O2+O2+E                        ","O2^-+O2(B1)=>O2+O2+E                        ",&
   "O2^-+N2=>O2+N2+E                            ","O2^-+N2(A3)=>O2+N2+E                        ",&
   "O2^-+N2(B3)=>O2+N2+E                        ","O3^-+O=>O2+O2+E                             ",&
   "NO^-+N=>N2O+E                               ","O3^-+N=>NO+O2+E                             ",&
   "N2O^-+N=>NO+N2+E                            ","NO2^-+N=>NO+NO+E                            ",&
   "NO3^-+N=>NO+NO2+E                           ","NO^-+O=>NO2+E                               ",&
   "N2O^-+O=>NO+NO+E                            ","NO2^-+O=>NO+O2+E                            ",&
   "NO3^-+O=>NO+O3+E                            ","O3^-+N2(A3)=>O3+N2+E                        ",&
   "NO^-+N2(A3)=>NO+N2+E                        ","N2O^-+N2(A3)=>N2O+N2+E                      ",&
   "NO2^-+N2(A3)=>NO2+N2+E                      ","NO3^-+N2(A3)=>NO3+N2+E                      ",&
   "O3^-+N2(B3)=>O3+N2+E                        ","NO^-+N2(B3)=>NO+N2+E                        ",&
   "N2O^-+N2(B3)=>N2O+N2+E                      ","NO2^-+N2(B3)=>NO2+N2+E                      "/
  data reaction_sign(181:216) &
  /"NO3^-+N2(B3)=>NO3+N2+E                      ","N2(A3)=>N2                                  ",&
   "N2(B3)=>N2(A3)                              ","N2(A`1)=>N2                                 ",&
   "N2(C3)=>N2(B3)                              ","O2(A1)=>O2                                  ",&
   "O2(B1)=>O2(A1)                              ","O2(B1)=>O2                                  ",&
   "O2(4.5EV)=>O2                               ","N2(A3)+O=>NO+N(2D)                          ",&
   "N2(A3)+O=>N2+O(1S)                          ","N2(A3)+N=>N2+N                              ",&
   "N2(A3)+N=>N2+N(2P)                          ","N2(A3)+O2=>N2+O+O(1D)                       ",&
   "N2(A3)+O2=>N2+O2(A1)                        ","N2(A3)+O2=>N2+O2(B1)                        ",&
   "N2(A3)+O2=>N2O+O                            ","N2(A3)+N2=>N2+N2                            ",&
   "N2(A3)+NO=>N2+NO                            ","N2(A3)+N2O=>N2+N+NO                         ",&
   "N2(A3)+NO2=>N2+O+NO                         ","N2(A3)+N2(A3)=>N2+N2(B3)                    ",&
   "N2(A3)+N2(A3)=>N2+N2(C3)                    ","N2(B3)+N2=>N2(A3)+N2                        ",&
   "N2(B3)+N2=>N2+N2                            ","N2(B3)+O2=>N2+O+O                           ",&
   "N2(B3)+NO=>N2(A3)+NO                        ","N2(C3)+N2=>N2(A`1)+N2                       ",&
   "N2(C3)+O2=>N2+O+O(1S)                       ","N2(A`1)+N2=>N2(B3)+N2                       ",&
   "N2(A`1)+O2=>N2+O+O                          ","N2(A`1)+NO=>N2+N+O                          ",&
   "N2(A`1)+N2(A3)=>N4^++E                      ","N2(A`1)+N2(A`1)=>N4^++E                     ",&
   "N+N+N2=>N2(A3)+N2                           ","N+N+O2=>N2(A3)+O2                           "/
  data reaction_sign(217:252) &
  /"N+N+NO=>N2(A3)+NO                           ","N+N+N=>N2(A3)+N                             ",&
   "N+N+O=>N2(A3)+O                             ","N+N+N2=>N2(B3)+N2                           ",&
   "N+N+O2=>N2(B3)+O2                           ","N+N+NO=>N2(B3)+NO                           ",&
   "N+N+N=>N2(B3)+N                             ","N+N+O=>N2(B3)+O                             ",&
   "N(2D)+O=>N+O(1D)                            ","N(2D)+O2=>NO+O                              ",&
   "N(2D)+NO=>N2+O                              ","N(2D)+N2O=>NO+N2                            ",&
   "N(2D)+N2=>N+N2                              ","N(2P)+N=>N+N                                ",&
   "N(2P)+O=>N+O                                ","N(2P)+N=>N(2D)+N                            ",&
   "N(2P)+N2=>N+N2                              ","N(2P)+N(2D)=>N2^++E                         ",&
   "N(2P)+O2=>NO+O                              ","N(2P)+NO=>N2(A3)+O                          ",&
   "O2(A1)+O=>O2+O                              ","O2(A1)+N=>NO+O                              ",&
   "O2(A1)+O2=>O2+O2                            ","O2(A1)+N2=>O2+N2                            ",&
   "O2(A1)+NO=>O2+NO                            ","O2(A1)+O3=>O2+O2+O(1D)                      ",&
   "O2(A1)+O2(A1)=>O2+O2(B1)                    ","O+O3=>O2+O2(A1)                             ",&
   "O2(B1)+O=>O2(A1)+O                          ","O2(B1)+O=>O2+O(1D)                          ",&
   "O2(B1)+O2=>O2(A1)+O2                        ","O2(B1)+N2=>O2(A1)+N2                        ",&
   "O2(B1)+NO=>O2(A1)+NO                        ","O2(B1)+O3=>O2+O2+O                          ",&
   "O2(4.5EV)+O=>O2+O(1S)                       ","O2(4.5EV)+O2=>O2(B1)+O2(B1)                 "/
  data reaction_sign(253:288) &
  /"O2(4.5EV)+N2=>O2(B1)+N2                     ","O(1D)+O=>O+O                                ",&
   "O(1D)+O2=>O+O2                              ","O(1D)+O2=>O+O2(A1)                          ",&
   "O(1D)+O2=>O+O2(B1)                          ","O(1D)+N2=>O+N2                              ",&
   "O(1D)+O3=>O2+O+O                            ","O(1D)+O3=>O2+O2                             ",&
   "O(1D)+NO=>O2+N                              ","O(1D)+N2O=>NO+NO                            ",&
   "O(1D)+N2O=>O2+N2                            ","O(1S)+O=>O(1D)+O                            ",&
   "O(1S)+N=>O+N                                ","O(1S)+O2=>O(1D)+O2                          ",&
   "O(1S)+O2=>O+O+O                             ","O(1S)+N2=>O+N2                              ",&
   "O(1S)+O2(A1)=>O+O2(4.5EV)                   ","O(1S)+O2(A1)=>O(1D)+O2(B1)                  ",&
   "O(1S)+O2(A1)=>O+O+O                         ","O(1S)+NO=>O+NO                              ",&
   "O(1S)+NO=>O(1D)+NO                          ","O(1S)+O3=>O2+O2                             ",&
   "O(1S)+O3=>O2+O+O(1D)                        ","O(1S)+N2O=>O+N2O                            ",&
   "O(1S)+N2O=>O(1D)+N2O                        ","N+NO=>O+N2                                  ",&
   "N+O2=>O+NO                                  ","N+NO2=>O+O+N2                               ",&
   "N+NO2=>O+N2O                                ","N+NO2=>N2+O2                                ",&
   "N+NO2=>NO+NO                                ","O+N2=>N+NO                                  ",&
   "O+NO=>N+O2                                  ","O+NO=>NO2                                   ",&
   "O+N2O=>N2+O2                                ","O+N2O=>NO+NO                                "/
  data reaction_sign(289:324) &
  /"O+NO2=>NO+O2                                ","O+NO3=>O2+NO2                               ",&
   "N2+O2=>O+N2O                                ","NO+NO=>N+NO2                                ",&
   "NO+NO=>O+N2O                                ","NO+NO=>N2+O2                                ",&
   "NO+O2=>O+NO2                                ","NO+O3=>O2+NO2                               ",&
   "NO+N2O=>N2+NO2                              ","NO+NO3=>NO2+NO2                             ",&
   "O2+O2=>O+O3                                 ","O2+NO2=>NO+O3                               ",&
   "NO2+NO2=>NO+NO+O2                           ","NO2+NO2=>NO+NO3                             ",&
   "NO2+O3=>O2+NO3                              ","NO2+NO3=>NO+NO2+O2                          ",&
   "NO3+O2=>NO2+O3                              ","NO3+NO3=>O2+NO2+NO2                         ",&
   "N+N=>N2^++E                                 ","N+O=>NO^++E                                 ",&
   "N2+N2=>N+N+N2                               ","N2+O2=>N+N+O2                               ",&
   "N2+NO=>N+N+NO                               ","N2+O=>N+N+O                                 ",&
   "N2+N=>N+N+N                                 ","O2+N2=>O+O+N2                               ",&
   "O2+O2=>O+O+O2                               ","O2+O=>O+O+O                                 ",&
   "O2+N=>O+O+N                                 ","O2+NO=>O+O+NO                               ",&
   "NO+N2=>N+O+N2                               ","NO+O2=>N+O+O2                               ",&
   "NO+O=>N+O+O                                 ","NO+N=>N+O+N                                 ",&
   "NO+NO=>N+O+NO                               ","O3+N2=>O2+O+N2                              "/
  data reaction_sign(325:360) &
  /"O3+O2=>O2+O+O2                              ","O3+N=>O2+O+N                                ",&
   "O3+O=>O2+O+O                                ","N2O+N2=>N2+O+N2                             ",&
   "N2O+O2=>N2+O+O2                             ","N2O+NO=>N2+O+NO                             ",&
   "N2O+N2O=>N2+O+N2O                           ","NO2+N2=>NO+O+N2                             ",&
   "NO2+O2=>NO+O+O2                             ","NO2+NO=>NO+O+NO                             ",&
   "NO2+NO2=>NO+O+NO2                           ","NO3+N2=>NO2+O+N2                            ",&
   "NO3+O2=>NO2+O+O2                            ","NO3+NO=>NO2+O+NO                            ",&
   "NO3+N=>NO2+O+N                              ","NO3+O=>NO2+O+O                              ",&
   "NO3+N2=>NO+O2+N2                            ","NO3+O2=>NO+O2+O2                            ",&
   "NO3+NO=>NO+O2+NO                            ","NO3+N=>NO+O2+N                              ",&
   "NO3+O=>NO+O2+O                              ","N2O5+ANY_NEUTRAL=>NO2+NO3+ANY_NEUTRAL       ",&
   "N+N+N2=>N2+N2                               ","N+N+O2=>N2+O2                               ",&
   "N+N+NO=>N2+NO                               ","N+N+N=>N2+N                                 ",&
   "N+N+O=>N2+O                                 ","O+O+N2=>O2+N2                               ",&
   "O+O+O2=>O2+O2                               ","O+O+N=>O2+N                                 ",&
   "O+O+O=>O2+O                                 ","O+O+NO=>O2+NO                               ",&
   "N+O+N2=>NO+N2                               ","N+O+O2=>NO+O2                               ",&
   "N+O+N=>NO+N                                 ","N+O+O=>NO+O                                 "/
  data reaction_sign(361:396) &
  /"N+O+NO=>NO+NO                               ","O+O2+N2=>O3+N2                              ",&
   "O+O2+O2=>O3+O2                              ","O+O2+NO=>O3+NO                              ",&
   "O+O2+N=>O3+N                                ","O+O2+O=>O3+O                                ",&
   "O+N2+ANY_NEUTRAL=>N2O+ANY_NEUTRAL           ","O+NO+N2=>NO2+N2                             ",&
   "O+NO+O2=>NO2+O2                             ","O+NO+NO=>NO2+NO                             ",&
   "O+NO2+N2=>NO3+N2                            ","O+NO2+O2=>NO3+O2                            ",&
   "O+NO2+N=>NO3+N                              ","O+NO2+O=>NO3+O                              ",&
   "O+NO2+NO=>NO3+NO                            ","NO2+NO3+ANY_NEUTRAL=>N2O5+ANY_NEUTRAL       ",&
   "N^++O=>N+O^+                                ","N^++O2=>O2^++N                              ",&
   "N^++O2=>NO^++O                              ","N^++O2=>O^++NO                              ",&
   "N^++O3=>NO^++O2                             ","N^++NO=>NO^++N                              ",&
   "N^++NO=>N2^++O                              ","N^++NO=>O^++N2                              ",&
   "N^++N2O=>NO^++N2                            ","O^++N2=>NO^++N                              ",&
   "O^++O2=>O2^++O                              ","O^++O3=>O2^++O2                             ",&
   "O^++NO=>NO^++O                              ","O^++NO=>O2^++N                              ",&
   "O^++N(2D)=>N^++O                            ","O^++N2O=>NO^++NO                            ",&
   "O^++N2O=>N2O^++O                            ","O^++N2O=>O2^++N2                            ",&
   "O^++NO2=>NO2^++O                            ","N2^++O2=>O2^++N2                            "/
  data reaction_sign(397:432) &
  /"N2^++O=>NO^++N                              ","N2^++O3=>O2^++O+N2                          ",&
   "N2^++N=>N^++N2                              ","N2^++NO=>NO^++N2                            ",&
   "N2^++N2O=>N2O^++N2                          ","N2^++N2O=>NO^++N+N2                         ",&
   "O2^++N2=>NO^++NO                            ","O2^++N=>NO^++O                              ",&
   "O2^++NO=>NO^++O2                            ","O2^++NO2=>NO^++O3                           ",&
   "O2^++NO2=>NO2^++O2                          ","N3^++O2=>O2^++N+N2                          ",&
   "N3^++O2=>NO2^++N2                           ","N3^++N=>N2^++N2                             ",&
   "N3^++NO=>NO^++N+N2                          ","N3^++NO=>N2O^++N2                           ",&
   "NO2^++NO=>NO^++NO2                          ","N2O^++NO=>NO^++N2O                          ",&
   "N4^++N2=>N2^++N2+N2                         ","N4^++O2=>O2^++N2+N2                         ",&
   "N4^++O=>O^++N2+N2                           ","N4^++N=>N^++N2+N2                           ",&
   "N4^++NO=>NO^++N2+N2                         ","O4^++N2=>O2^+N2+O2                          ",&
   "O4^++O2=>O2^++O2+O2                         ","O4^++O2(A1)=>O2^++O2+O2                     ",&
   "O4^++O2(B1)=>O2^++O2+O2                     ","O4^++O=>O2^++O3                             ",&
   "O4^++NO=>NO^++O2+O2                         ","O2^+N2+N2=>O2^++N2+N2                       ",&
   "O2^+N2+O2=>O4^++N2                          ","N^++N2+N2=>N3^++N2                          ",&
   "N^++O+ANY_NEUTRAL=>NO^++ANY_NEUTRAL         ","N^++N+ANY_NEUTRAL=>N2^++ANY_NEUTRAL         ",&
   "O^++N2+ANY_NEUTRAL=>NO^++N+ANY_NEUTRAL      ","O^++O+ANY_NEUTRAL=>O2^++ANY_NEUTRAL         "/
  data reaction_sign(433:468) &
  /"O^++N+ANY_NEUTRAL=>NO^++ANY_NEUTRAL         ","N2^++N2+N2=>N4^++N2                         ",&
   "N2^++N+N2=>N3^++N2                          ","O2^++O2+O2=>O4^++O2                         ",&
   "O2^++N2+N2=>O2^+N2+N2                       ","O^-+O2(A1)=>O2^-+O                          ",&
   "O^-+O3=>O3^-+O                              ","O^-+NO2=>NO2^-+O                            ",&
   "O^-+N2O=>NO^-+NO                            ","O^-+N2O=>N2O^-+O                            ",&
   "O2^-+O=>O^-+O2                              ","O2^-+O3=>O3^-+O2                            ",&
   "O2^-+NO2=>NO2^-+O2                          ","O2^-+NO3=>NO3^-+O2                          ",&
   "O3^-+O=>O2^-+O2                             ","O3^-+NO=>NO3^-+O                            ",&
   "O3^-+NO=>NO2^-+O2                           ","O3^-+NO2=>NO2^-+O3                          ",&
   "O3^-+NO2=>NO3^-+O2                          ","O3^-+NO3=>NO3^-+O3                          ",&
   "NO^-+O2=>O2^-+NO                            ","NO^-+NO2=>NO2^-+NO                          ",&
   "NO^-+N2O=>NO2^-+N2                          ","NO2^-+O3=>NO3^-+O2                          ",&
   "NO2^-+NO2=>NO3^-+NO                         ","NO2^-+NO3=>NO3^-+NO2                        ",&
   "NO2^-+N2O5=>NO3^-+NO2+NO2                   ","NO3^-+NO=>NO2^-+NO2                         ",&
   "O4^-+N2=>O2^-+O2+N2                         ","O4^-+O2=>O2^-+O2+O2                         ",&
   "O4^-+O=>O3^-+O2                             ","O4^-+O=>O^-+O2+O2                           ",&
   "O4^-+O2(A1)=>O2^-+O2+O2                     ","O4^-+O2(B1)=>O2^-+O2+O2                     ",&
   "O4^-+NO=>NO3^-+O2                           ","O^-+O2+ANY_NEUTRAL=>O3^-+ANY_NEUTRAL        "/
  data reaction_sign(469:504) &
  /"O^-+NO+ANY_NEUTRAL=>NO2^-+ANY_NEUTRAL       ","O2^-+O2+ANY_NEUTRAL=>O4^-+ANY_NEUTRAL       ",&
   "O^-+N^+=>O+N                                ","O^-+N2^+=>O+N2                              ",&
   "O^-+O^+=>O+O                                ","O^-+O2^+=>O+O2                              ",&
   "O^-+NO^+=>O+NO                              ","O^-+N2O^+=>O+N2O                            ",&
   "O^-+NO2^+=>O+NO2                            ","O2^-+N^+=>O2+N                              ",&
   "O2^-+N2^+=>O2+N2                            ","O2^-+O^+=>O2+O                              ",&
   "O2^-+O2^+=>O2+O2                            ","O2^-+NO^+=>O2+NO                            ",&
   "O2^-+N2O^+=>O2+N2O                          ","O2^-+NO2^+=>O2+NO2                          ",&
   "O3^-+N^+=>O3+N                              ","O3^-+N2^+=>O3+N2                            ",&
   "O3^-+O^+=>O3+O                              ","O3^-+O2^+=>O3+O2                            ",&
   "O3^-+NO^+=>O3+NO                            ","O3^-+N2O^+=>O3+N2O                          ",&
   "O3^-+NO2^+=>O3+NO2                          ","NO^-+N^+=>NO+N                              ",&
   "NO^-+N2^+=>NO+N2                            ","NO^-+O^+=>NO+O                              ",&
   "NO^-+O2^+=>NO+O2                            ","NO^-+NO^+=>NO+NO                            ",&
   "NO^-+N2O^+=>NO+N2O                          ","NO^-+NO2^+=>NO+NO2                          ",&
   "N2O^-+N^+=>N2O+N                            ","N2O^-+N2^+=>N2O+N2                          ",&
   "N2O^-+O^+=>N2O+O                            ","N2O^-+O2^+=>N2O+O2                          ",&
   "N2O^-+NO^+=>N2O+NO                          ","N2O^-+N2O^+=>N2O+N2O                        "/
  data reaction_sign(505:540) &
  /"N2O^-+NO2^+=>N2O+NO2                        ","NO2^-+N^+=>NO2+N                            ",&
   "NO2^-+N2^+=>NO2+N2                          ","NO2^-+O^+=>NO2+O                            ",&
   "NO2^-+O2^+=>NO2+O2                          ","NO2^-+NO^+=>NO2+NO                          ",&
   "NO2^-+N2O^+=>NO2+N2O                        ","NO2^-+NO2^+=>NO2+NO2                        ",&
   "NO3^-+N^+=>NO3+N                            ","NO3^-+N2^+=>NO3+N2                          ",&
   "NO3^-+O^+=>NO3+O                            ","NO3^-+O2^+=>NO3+O2                          ",&
   "NO3^-+NO^+=>NO3+NO                          ","NO3^-+N2O^+=>NO3+N2O                        ",&
   "NO3^-+NO2^+=>NO3+NO2                        ","O^-+N2^+=>O+N+N                             ",&
   "O^-+N3^+=>O+N+N2                            ","O^-+N4^+=>O+N2+N2                           ",&
   "O^-+O2^+=>O+O+O                             ","O^-+O4^+=>O+O2+O2                           ",&
   "O^-+NO^+=>O+N+O                             ","O^-+N2O^+=>O+N2+O                           ",&
   "O^-+NO2^+=>O+N+O2                           ","O^-+O2^+N2=>O+O2+N2                         ",&
   "O2^-+N2^+=>O2+N+N                           ","O2^-+N3^+=>O2+N+N2                          ",&
   "O2^-+N4^+=>O2+N2+N2                         ","O2^-+O2^+=>O2+O+O                           ",&
   "O2^-+O4^+=>O2+O2+O2                         ","O2^-+NO^+=>O2+N+O                           ",&
   "O2^-+N2O^+=>O2+N2+O                         ","O2^-+NO2^+=>O2+N+O2                         ",&
   "O2^-+O2^+N2=>O2+O2+N2                       ","O3^-+N2^+=>O3+N+N                           ",&
   "O3^-+N3^+=>O3+N+N2                          ","O3^-+N4^+=>O3+N2+N2                         "/
  data reaction_sign(541:576) &
  /"O3^-+O2^+=>O3+O+O                           ","O3^-+O4^+=>O3+O2+O2                         ",&
   "O3^-+NO^+=>O3+N+O                           ","O3^-+N2O^+=>O3+N2+O                         ",&
   "O3^-+NO2^+=>O3+N+O2                         ","O3^-+O2^+N2=>O3+O2+N2                       ",&
   "NO^-+N2^+=>NO+N+N                           ","NO^-+N3^+=>NO+N+N2                          ",&
   "NO^-+N4^+=>NO+N2+N2                         ","NO^-+O2^+=>NO+O+O                           ",&
   "NO^-+O4^+=>NO+O2+O2                         ","NO^-+NO^+=>NO+N+O                           ",&
   "NO^-+N2O^+=>NO+N2+O                         ","NO^-+NO2^+=>NO+N+O2                         ",&
   "NO^-+O2^+N2=>NO+O2+N2                       ","N2O^-+N2^+=>N2O+N+N                         ",&
   "N2O^-+N3^+=>N2O+N+N2                        ","N2O^-+N4^+=>N2O+N2+N2                       ",&
   "N2O^-+O2^+=>N2O+O+O                         ","N2O^-+O4^+=>N2O+O2+O2                       ",&
   "N2O^-+NO^+=>N2O+N+O                         ","N2O^-+N2O^+=>N2O+N2+O                       ",&
   "N2O^-+NO2^+=>N2O+N+O2                       ","N2O^-+O2^+N2=>N2O+O2+N2                     ",&
   "NO2^-+N2^+=>NO2+N+N                         ","NO2^-+N3^+=>NO2+N+N2                        ",&
   "NO2^-+N4^+=>NO2+N2+N2                       ","NO2^-+O2^+=>NO2+O+O                         ",&
   "NO2^-+O4^+=>NO2+O2+O2                       ","NO2^-+NO^+=>NO2+N+O                         ",&
   "NO2^-+N2O^+=>NO2+N2+O                       ","NO2^-+NO2^+=>NO2+N+O2                       ",&
   "NO2^-+O2^+N2=>NO2+O2+N2                     ","NO3^-+N2^+=>NO3+N+N                         ",&
   "NO3^-+N3^+=>NO3+N+N2                        ","NO3^-+N4^+=>NO3+N2+N2                       "/
  data reaction_sign(577:612) &
  /"NO3^-+O2^+=>NO3+O+O                         ","NO3^-+O4^+=>NO3+O2+O2                       ",&
   "NO3^-+NO^+=>NO3+N+O                         ","NO3^-+N2O^+=>NO3+N2+O                       ",&
   "NO3^-+NO2^+=>NO3+N+O2                       ","NO3^-+O2^+N2=>NO3+O2+N2                     ",&
   "O4^-+N^+=>O2+O2+N                           ","O4^-+N2^+=>O2+O2+N2                         ",&
   "O4^-+O^+=>O2+O2+O                           ","O4^-+O2^+=>O2+O2+O2                         ",&
   "O4^-+NO^+=>O2+O2+NO                         ","O4^-+N2O^+=>O2+O2+N2O                       ",&
   "O4^-+NO2^+=>O2+O2+NO2                       ","O4^-+N3^+=>O2+O2+N2+N                       ",&
   "O4^-+N4^+=>O2+O2+N2+N2                      ","O4^-+O4^+=>O2+O2+O2+O2                      ",&
   "O4^-+O2^+N2=>O2+O2+O2+N2                    ","O^-+N^++ANY_NEUTRAL=>O+N+ANY_NEUTRAL        ",&
   "O^-+N2^++ANY_NEUTRAL=>O+N2+ANY_NEUTRAL      ","O^-+O^++ANY_NEUTRAL=>O+O+ANY_NEUTRAL        ",&
   "O^-+O2^++ANY_NEUTRAL=>O+O2+ANY_NEUTRAL      ","O^-+NO^++ANY_NEUTRAL=>O+NO+ANY_NEUTRAL      ",&
   "O2^-+N^++ANY_NEUTRAL=>O2+N+ANY_NEUTRAL      ","O2^-+N2^++ANY_NEUTRAL=>O2+N2+ANY_NEUTRAL    ",&
   "O2^-+O^++ANY_NEUTRAL=>O2+O+ANY_NEUTRAL      ","O2^-+O2^++ANY_NEUTRAL=>O2+O2+ANY_NEUTRAL    ",&
   "O2^-+NO^++ANY_NEUTRAL=>O2+NO+ANY_NEUTRAL    ","O^-+N^++ANY_NEUTRAL=>NO+ANY_NEUTRAL         ",&
   "O^-+N2^++ANY_NEUTRAL=>N2O+ANY_NEUTRAL       ","O^-+O^++ANY_NEUTRAL=>O2+ANY_NEUTRAL         ",&
   "O^-+O2^++ANY_NEUTRAL=>O3+ANY_NEUTRAL        ","O^-+NO^++ANY_NEUTRAL=>NO2+ANY_NEUTRAL       ",&
   "O2^-+N^++ANY_NEUTRAL=>NO2+ANY_NEUTRAL       ","O2^-+O^++ANY_NEUTRAL=>O3+ANY_NEUTRAL        ",&
   "O2^-+NO^++ANY_NEUTRAL=>NO3+ANY_NEUTRAL      ","O3^-+N^++ANY_NEUTRAL=>O3+N+ANY_NEUTRAL      "/
  data reaction_sign(613:646) &
  /"O3^-+N2^++ANY_NEUTRAL=>O3+N2+ANY_NEUTRAL    ","O3^-+O^++ANY_NEUTRAL=>O3+O+ANY_NEUTRAL      ",&
   "O3^-+O2^++ANY_NEUTRAL=>O3+O2+ANY_NEUTRAL    ","O3^-+NO^++ANY_NEUTRAL=>O3+NO+ANY_NEUTRAL    ",&
   "O3^-+N2O^++ANY_NEUTRAL=>O3+N2O+ANY_NEUTRAL  ","O3^-+NO2^++ANY_NEUTRAL=>O3+NO2+ANY_NEUTRAL  ",&
   "NO^-+N^++ANY_NEUTRAL=>NO+N+ANY_NEUTRAL      ","NO^-+N2^++ANY_NEUTRAL=>NO+N2+ANY_NEUTRAL    ",&
   "NO^-+O^++ANY_NEUTRAL=>NO+O+ANY_NEUTRAL      ","NO^-+O2^++ANY_NEUTRAL=>NO+O2+ANY_NEUTRAL    ",&
   "NO^-+NO^++ANY_NEUTRAL=>NO+NO+ANY_NEUTRAL    ","NO^-+N2O^++ANY_NEUTRAL=>NO+N2O+ANY_NEUTRAL  ",&
   "NO^-+NO2^++ANY_NEUTRAL=>NO+NO2+ANY_NEUTRAL  ","N2O^-+N^++ANY_NEUTRAL=>N2O+N+ANY_NEUTRAL    ",&
   "N2O^-+N2^++ANY_NEUTRAL=>N2O+N2+ANY_NEUTRAL  ","N2O^-+O^++ANY_NEUTRAL=>N2O+O+ANY_NEUTRAL    ",&
   "N2O^-+O2^++ANY_NEUTRAL=>N2O+O2+ANY_NEUTRAL  ","N2O^-+NO^++ANY_NEUTRAL=>N2O+NO+ANY_NEUTRAL  ",&
   "N2O^-+N2O^++ANY_NEUTRAL=>N2O+N2O+ANY_NEUTRAL","N2O^-+NO2^++ANY_NEUTRAL=>N2O+NO2+ANY_NEUTRAL",&
   "NO2^-+N^++ANY_NEUTRAL=>NO2+N+ANY_NEUTRAL    ","NO2^-+N2^++ANY_NEUTRAL=>NO2+N2+ANY_NEUTRAL  ",&
   "NO2^-+O^++ANY_NEUTRAL=>NO2+O+ANY_NEUTRAL    ","NO2^-+O2^++ANY_NEUTRAL=>NO2+O2+ANY_NEUTRAL  ",&
   "NO2^-+NO^++ANY_NEUTRAL=>NO2+NO+ANY_NEUTRAL  ","NO2^-+N2O^++ANY_NEUTRAL=>NO2+N2O+ANY_NEUTRAL",&
   "NO2^-+NO2^++ANY_NEUTRAL=>NO2+NO2+ANY_NEUTRAL","NO3^-+N^++ANY_NEUTRAL=>NO3+N+ANY_NEUTRAL    ",&
   "NO3^-+N2^++ANY_NEUTRAL=>NO3+N2+ANY_NEUTRAL  ","NO3^-+O^++ANY_NEUTRAL=>NO3+O+ANY_NEUTRAL    ",&
   "NO3^-+O2^++ANY_NEUTRAL=>NO3+O2+ANY_NEUTRAL  ","NO3^-+NO^++ANY_NEUTRAL=>NO3+NO+ANY_NEUTRAL  ",&
   "NO3^-+N2O^++ANY_NEUTRAL=>NO3+N2O+ANY_NEUTRAL","NO3^-+NO2^++ANY_NEUTRAL=>NO3+NO2+ANY_NEUTRAL"/
  data bolsig_species(1:bolsig_species_max) &
  /"N2    ","N2(V1)","N2(V2)","N2(V3)","N2(V4)","N2(V5)","N2(V6)","N2(V7)","N2(V8)","N2(A3)","O2    ","O2(V1)","O2(V2)","O2(V3)",&
   "O2(V4)","O2(A1)","N     ","O     ","NO    ","O3    ","N2O   "/
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
      write(ifile_unit,"(646(i3))",err=200) int(mrtm(i,:))
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
    write(ifile_unit,"(1x,A12,646(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
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
      write(ifile_unit,"(647(1pe13.4))") densav(0,2), rrt_loc(:)
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
  reac_source_local(21,010) = - reac_rate_local(010) 
  reac_source_local(22,010) = + reac_rate_local(010) 
  reac_source_local(21,011) = - reac_rate_local(011) 
  reac_source_local(22,011) = + reac_rate_local(011) 
  reac_source_local(21,012) = - reac_rate_local(012) 
  reac_source_local(23,012) = + reac_rate_local(012) 
  reac_source_local(21,013) = - reac_rate_local(013) 
  reac_source_local(23,013) = + reac_rate_local(013) 
  reac_source_local(21,014) = - reac_rate_local(014) 
  reac_source_local(24,014) = + reac_rate_local(014) 
  reac_source_local(21,015) = - reac_rate_local(015) 
  reac_source_local(25,015) = + reac_rate_local(015) 
  reac_source_local(21,016) = + reac_rate_local(016) 
  reac_source_local(22,016) = - reac_rate_local(016) 
  reac_source_local(21,017) = + reac_rate_local(017) 
  reac_source_local(23,017) = - reac_rate_local(017) 
  reac_source_local(21,018) = + reac_rate_local(018) 
  reac_source_local(24,018) = - reac_rate_local(018) 
  reac_source_local(21,019) = + reac_rate_local(019) 
  reac_source_local(25,019) = - reac_rate_local(019) 
  reac_source_local(01,020) = + reac_rate_local(020) 
  reac_source_local(02,020) = - reac_rate_local(020) 
  reac_source_local(02,021) = + reac_rate_local(021) 
  reac_source_local(03,021) = - reac_rate_local(021) 
  reac_source_local(03,022) = + reac_rate_local(022) 
  reac_source_local(04,022) = - reac_rate_local(022) 
  reac_source_local(04,023) = + reac_rate_local(023) 
  reac_source_local(05,023) = - reac_rate_local(023) 
  reac_source_local(05,024) = + reac_rate_local(024) 
  reac_source_local(06,024) = - reac_rate_local(024) 
  reac_source_local(06,025) = + reac_rate_local(025) 
  reac_source_local(07,025) = - reac_rate_local(025) 
  reac_source_local(07,026) = + reac_rate_local(026) 
  reac_source_local(08,026) = - reac_rate_local(026) 
  reac_source_local(08,027) = + reac_rate_local(027) 
  reac_source_local(09,027) = - reac_rate_local(027) 
  reac_source_local(01,028) = - reac_rate_local(028) 
  reac_source_local(02,028) = + reac_rate_local(028) 
  reac_source_local(02,029) = - reac_rate_local(029) 
  reac_source_local(03,029) = + reac_rate_local(029) 
  reac_source_local(03,030) = - reac_rate_local(030) 
  reac_source_local(04,030) = + reac_rate_local(030) 
  reac_source_local(04,031) = - reac_rate_local(031) 
  reac_source_local(05,031) = + reac_rate_local(031) 
  reac_source_local(05,032) = - reac_rate_local(032) 
  reac_source_local(06,032) = + reac_rate_local(032) 
  reac_source_local(06,033) = - reac_rate_local(033) 
  reac_source_local(07,033) = + reac_rate_local(033) 
  reac_source_local(07,034) = - reac_rate_local(034) 
  reac_source_local(08,034) = + reac_rate_local(034) 
  reac_source_local(08,035) = - reac_rate_local(035) 
  reac_source_local(09,035) = + reac_rate_local(035) 
  reac_source_local(01,036) = + reac_rate_local(036) 
  reac_source_local(02,036) = - reac_rate_local(036) 
  reac_source_local(02,037) = + reac_rate_local(037) 
  reac_source_local(03,037) = - reac_rate_local(037) 
  reac_source_local(03,038) = + reac_rate_local(038) 
  reac_source_local(04,038) = - reac_rate_local(038) 
  reac_source_local(04,039) = + reac_rate_local(039) 
  reac_source_local(05,039) = - reac_rate_local(039) 
  reac_source_local(05,040) = + reac_rate_local(040) 
  reac_source_local(06,040) = - reac_rate_local(040) 
  reac_source_local(06,041) = + reac_rate_local(041) 
  reac_source_local(07,041) = - reac_rate_local(041) 
  reac_source_local(07,042) = + reac_rate_local(042) 
  reac_source_local(08,042) = - reac_rate_local(042) 
  reac_source_local(08,043) = + reac_rate_local(043) 
  reac_source_local(09,043) = - reac_rate_local(043) 
  reac_source_local(01,044) = - reac_rate_local(044) 
  reac_source_local(02,044) = + reac_rate_local(044) 
  reac_source_local(02,045) = - reac_rate_local(045) 
  reac_source_local(03,045) = + reac_rate_local(045) 
  reac_source_local(03,046) = - reac_rate_local(046) 
  reac_source_local(04,046) = + reac_rate_local(046) 
  reac_source_local(04,047) = - reac_rate_local(047) 
  reac_source_local(05,047) = + reac_rate_local(047) 
  reac_source_local(05,048) = - reac_rate_local(048) 
  reac_source_local(06,048) = + reac_rate_local(048) 
  reac_source_local(06,049) = - reac_rate_local(049) 
  reac_source_local(07,049) = + reac_rate_local(049) 
  reac_source_local(07,050) = - reac_rate_local(050) 
  reac_source_local(08,050) = + reac_rate_local(050) 
  reac_source_local(08,051) = - reac_rate_local(051) 
  reac_source_local(09,051) = + reac_rate_local(051) 
  reac_source_local(01,052) = + reac_rate_local(052) 
  reac_source_local(02,052) = - reac_rate_local(052) 
  reac_source_local(02,053) = + reac_rate_local(053) 
  reac_source_local(03,053) = - reac_rate_local(053) 
  reac_source_local(03,054) = + reac_rate_local(054) 
  reac_source_local(04,054) = - reac_rate_local(054) 
  reac_source_local(04,055) = + reac_rate_local(055) 
  reac_source_local(05,055) = - reac_rate_local(055) 
  reac_source_local(05,056) = + reac_rate_local(056) 
  reac_source_local(06,056) = - reac_rate_local(056) 
  reac_source_local(06,057) = + reac_rate_local(057) 
  reac_source_local(07,057) = - reac_rate_local(057) 
  reac_source_local(07,058) = + reac_rate_local(058) 
  reac_source_local(08,058) = - reac_rate_local(058) 
  reac_source_local(08,059) = + reac_rate_local(059) 
  reac_source_local(09,059) = - reac_rate_local(059) 
  reac_source_local(01,060) = - reac_rate_local(060) 
  reac_source_local(02,060) = + reac_rate_local(060) 
  reac_source_local(02,061) = - reac_rate_local(061) 
  reac_source_local(03,061) = + reac_rate_local(061) 
  reac_source_local(03,062) = - reac_rate_local(062) 
  reac_source_local(04,062) = + reac_rate_local(062) 
  reac_source_local(04,063) = - reac_rate_local(063) 
  reac_source_local(05,063) = + reac_rate_local(063) 
  reac_source_local(05,064) = - reac_rate_local(064) 
  reac_source_local(06,064) = + reac_rate_local(064) 
  reac_source_local(06,065) = - reac_rate_local(065) 
  reac_source_local(07,065) = + reac_rate_local(065) 
  reac_source_local(07,066) = - reac_rate_local(066) 
  reac_source_local(08,066) = + reac_rate_local(066) 
  reac_source_local(08,067) = - reac_rate_local(067) 
  reac_source_local(09,067) = + reac_rate_local(067) 
  reac_source_local(21,068) = + reac_rate_local(068) 
  reac_source_local(22,068) = - reac_rate_local(068) 
  reac_source_local(22,069) = + reac_rate_local(069) 
  reac_source_local(23,069) = - reac_rate_local(069) 
  reac_source_local(23,070) = + reac_rate_local(070) 
  reac_source_local(24,070) = - reac_rate_local(070) 
  reac_source_local(24,071) = + reac_rate_local(071) 
  reac_source_local(25,071) = - reac_rate_local(071) 
  reac_source_local(21,072) = - reac_rate_local(072) 
  reac_source_local(22,072) = + reac_rate_local(072) 
  reac_source_local(22,073) = - reac_rate_local(073) 
  reac_source_local(23,073) = + reac_rate_local(073) 
  reac_source_local(23,074) = - reac_rate_local(074) 
  reac_source_local(24,074) = + reac_rate_local(074) 
  reac_source_local(24,075) = - reac_rate_local(075) 
  reac_source_local(25,075) = + reac_rate_local(075) 
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
  reac_source_local(01,084) = - reac_rate_local(084) 
  reac_source_local(10,084) = + reac_rate_local(084) 
  reac_source_local(01,085) = - reac_rate_local(085) 
  reac_source_local(10,085) = + reac_rate_local(085) 
  reac_source_local(01,086) = - reac_rate_local(086) 
  reac_source_local(10,086) = + reac_rate_local(086) 
  reac_source_local(01,087) = - reac_rate_local(087) 
  reac_source_local(11,087) = + reac_rate_local(087) 
  reac_source_local(01,088) = - reac_rate_local(088) 
  reac_source_local(11,088) = + reac_rate_local(088) 
  reac_source_local(01,089) = - reac_rate_local(089) 
  reac_source_local(11,089) = + reac_rate_local(089) 
  reac_source_local(01,090) = - reac_rate_local(090) 
  reac_source_local(12,090) = + reac_rate_local(090) 
  reac_source_local(01,091) = - reac_rate_local(091) 
  reac_source_local(12,091) = + reac_rate_local(091) 
  reac_source_local(01,092) = - reac_rate_local(092) 
  reac_source_local(12,092) = + reac_rate_local(092) 
  reac_source_local(01,093) = - reac_rate_local(093) 
  reac_source_local(13,093) = + reac_rate_local(093) 
  reac_source_local(01,094) = - reac_rate_local(094) 
  reac_source_local(13,094) = + reac_rate_local(094) 
  reac_source_local(01,095) = - reac_rate_local(095) 
  reac_source_local(13,095) = + reac_rate_local(095) 
  reac_source_local(21,096) = - reac_rate_local(096) 
  reac_source_local(26,096) = + reac_rate_local(096) 
  reac_source_local(21,097) = - reac_rate_local(097) 
  reac_source_local(27,097) = + reac_rate_local(097) 
  reac_source_local(21,098) = - reac_rate_local(098) 
  reac_source_local(28,098) = + reac_rate_local(098) 
  reac_source_local(21,099) = - reac_rate_local(099) 
  reac_source_local(29,099) = + reac_rate_local(099) * 2.d0
  reac_source_local(21,100) = - reac_rate_local(100) 
  reac_source_local(29,100) = + reac_rate_local(100) 
  reac_source_local(30,100) = + reac_rate_local(100) 
  reac_source_local(21,101) = - reac_rate_local(101) 
  reac_source_local(29,101) = + reac_rate_local(101) 
  reac_source_local(31,101) = + reac_rate_local(101) 
  reac_source_local(29,102) = - reac_rate_local(102) 
  reac_source_local(30,102) = + reac_rate_local(102) 
  reac_source_local(29,103) = - reac_rate_local(103) 
  reac_source_local(31,103) = + reac_rate_local(103) 
  reac_source_local(21,104) = + reac_rate_local(104) 
  reac_source_local(26,104) = - reac_rate_local(104) 
  reac_source_local(14,105) = - reac_rate_local(105) 
  reac_source_local(17,105) = + reac_rate_local(105) 
  reac_source_local(53,105) = + reac_rate_local(105) 
  reac_source_local(29,106) = - reac_rate_local(106) 
  reac_source_local(33,106) = + reac_rate_local(106) 
  reac_source_local(53,106) = + reac_rate_local(106) 
  reac_source_local(01,107) = - reac_rate_local(107) 
  reac_source_local(18,107) = + reac_rate_local(107) 
  reac_source_local(53,107) = + reac_rate_local(107) 
  reac_source_local(10,108) = - reac_rate_local(108) 
  reac_source_local(18,108) = + reac_rate_local(108) 
  reac_source_local(53,108) = + reac_rate_local(108) 
  reac_source_local(21,109) = - reac_rate_local(109) 
  reac_source_local(34,109) = + reac_rate_local(109) 
  reac_source_local(53,109) = + reac_rate_local(109) 
  reac_source_local(26,110) = - reac_rate_local(110) 
  reac_source_local(34,110) = + reac_rate_local(110) 
  reac_source_local(53,110) = + reac_rate_local(110) 
  reac_source_local(40,111) = - reac_rate_local(111) 
  reac_source_local(45,111) = + reac_rate_local(111) 
  reac_source_local(53,111) = + reac_rate_local(111) 
  reac_source_local(41,112) = - reac_rate_local(112) 
  reac_source_local(46,112) = + reac_rate_local(112) 
  reac_source_local(53,112) = + reac_rate_local(112) 
  reac_source_local(14,113) = + reac_rate_local(113) * 2.d0
  reac_source_local(18,113) = - reac_rate_local(113) 
  reac_source_local(53,113) = - reac_rate_local(113) 
  reac_source_local(14,114) = + reac_rate_local(114) 
  reac_source_local(15,114) = + reac_rate_local(114) 
  reac_source_local(18,114) = - reac_rate_local(114) 
  reac_source_local(53,114) = - reac_rate_local(114) 
  reac_source_local(14,115) = + reac_rate_local(115) 
  reac_source_local(16,115) = + reac_rate_local(115) 
  reac_source_local(18,115) = - reac_rate_local(115) 
  reac_source_local(53,115) = - reac_rate_local(115) 
  reac_source_local(29,116) = + reac_rate_local(116) * 2.d0
  reac_source_local(34,116) = - reac_rate_local(116) 
  reac_source_local(53,116) = - reac_rate_local(116) 
  reac_source_local(29,117) = + reac_rate_local(117) 
  reac_source_local(30,117) = + reac_rate_local(117) 
  reac_source_local(34,117) = - reac_rate_local(117) 
  reac_source_local(53,117) = - reac_rate_local(117) 
  reac_source_local(29,118) = + reac_rate_local(118) 
  reac_source_local(31,118) = + reac_rate_local(118) 
  reac_source_local(34,118) = - reac_rate_local(118) 
  reac_source_local(53,118) = - reac_rate_local(118) 
  reac_source_local(14,119) = + reac_rate_local(119) 
  reac_source_local(29,119) = + reac_rate_local(119) 
  reac_source_local(45,119) = - reac_rate_local(119) 
  reac_source_local(53,119) = - reac_rate_local(119) 
  reac_source_local(15,120) = + reac_rate_local(120) 
  reac_source_local(29,120) = + reac_rate_local(120) 
  reac_source_local(45,120) = - reac_rate_local(120) 
  reac_source_local(53,120) = - reac_rate_local(120) 
  reac_source_local(01,121) = + reac_rate_local(121) 
  reac_source_local(14,121) = + reac_rate_local(121) 
  reac_source_local(19,121) = - reac_rate_local(121) 
  reac_source_local(53,121) = - reac_rate_local(121) 
  reac_source_local(01,122) = + reac_rate_local(122) * 2.d0
  reac_source_local(20,122) = - reac_rate_local(122) 
  reac_source_local(53,122) = - reac_rate_local(122) 
  reac_source_local(01,123) = + reac_rate_local(123) 
  reac_source_local(29,123) = + reac_rate_local(123) 
  reac_source_local(46,123) = - reac_rate_local(123) 
  reac_source_local(53,123) = - reac_rate_local(123) 
  reac_source_local(29,124) = + reac_rate_local(124) 
  reac_source_local(40,124) = + reac_rate_local(124) 
  reac_source_local(47,124) = - reac_rate_local(124) 
  reac_source_local(53,124) = - reac_rate_local(124) 
  reac_source_local(21,125) = + reac_rate_local(125) * 2.d0
  reac_source_local(35,125) = - reac_rate_local(125) 
  reac_source_local(53,125) = - reac_rate_local(125) 
  reac_source_local(01,126) = + reac_rate_local(126) 
  reac_source_local(21,126) = + reac_rate_local(126) 
  reac_source_local(52,126) = - reac_rate_local(126) 
  reac_source_local(53,126) = - reac_rate_local(126) 
  reac_source_local(14,127) = + reac_rate_local(127) 
  reac_source_local(17,127) = - reac_rate_local(127) 
  reac_source_local(53,127) = - reac_rate_local(127) 
  reac_source_local(29,128) = + reac_rate_local(128) 
  reac_source_local(33,128) = - reac_rate_local(128) 
  reac_source_local(53,128) = - reac_rate_local(128) 
  reac_source_local(14,129) = + reac_rate_local(129) 
  reac_source_local(17,129) = - reac_rate_local(129) 
  reac_source_local(53,129) = - reac_rate_local(129) 
  reac_source_local(29,130) = + reac_rate_local(130) 
  reac_source_local(33,130) = - reac_rate_local(130) 
  reac_source_local(53,130) = - reac_rate_local(130) 
  reac_source_local(21,131) = - reac_rate_local(131) 
  reac_source_local(29,131) = + reac_rate_local(131) 
  reac_source_local(36,131) = + reac_rate_local(131) 
  reac_source_local(53,131) = - reac_rate_local(131) 
  reac_source_local(14,132) = + reac_rate_local(132) 
  reac_source_local(36,132) = + reac_rate_local(132) 
  reac_source_local(40,132) = - reac_rate_local(132) 
  reac_source_local(53,132) = - reac_rate_local(132) 
  reac_source_local(21,133) = + reac_rate_local(133) 
  reac_source_local(32,133) = - reac_rate_local(133) 
  reac_source_local(36,133) = + reac_rate_local(133) 
  reac_source_local(53,133) = - reac_rate_local(133) 
  reac_source_local(29,134) = + reac_rate_local(134) 
  reac_source_local(32,134) = - reac_rate_local(134) 
  reac_source_local(37,134) = + reac_rate_local(134) 
  reac_source_local(53,134) = - reac_rate_local(134) 
  reac_source_local(14,135) = + reac_rate_local(135) 
  reac_source_local(41,135) = - reac_rate_local(135) 
  reac_source_local(48,135) = + reac_rate_local(135) 
  reac_source_local(53,135) = - reac_rate_local(135) 
  reac_source_local(21,136) = - reac_rate_local(136) 
  reac_source_local(37,136) = + reac_rate_local(136) 
  reac_source_local(53,136) = - reac_rate_local(136) 
  reac_source_local(36,137) = + reac_rate_local(137) 
  reac_source_local(40,137) = + reac_rate_local(137) 
  reac_source_local(42,137) = - reac_rate_local(137) 
  reac_source_local(53,137) = - reac_rate_local(137) 
  reac_source_local(29,138) = - reac_rate_local(138) 
  reac_source_local(36,138) = + reac_rate_local(138) 
  reac_source_local(53,138) = - reac_rate_local(138) 
  reac_source_local(21,139) = - reac_rate_local(139) 
  reac_source_local(37,139) = + reac_rate_local(139) 
  reac_source_local(53,139) = - reac_rate_local(139) 
  reac_source_local(32,140) = - reac_rate_local(140) 
  reac_source_local(38,140) = + reac_rate_local(140) 
  reac_source_local(53,140) = - reac_rate_local(140) 
  reac_source_local(40,141) = - reac_rate_local(141) 
  reac_source_local(48,141) = + reac_rate_local(141) 
  reac_source_local(53,141) = - reac_rate_local(141) 
  reac_source_local(41,142) = - reac_rate_local(142) 
  reac_source_local(49,142) = + reac_rate_local(142) 
  reac_source_local(53,142) = - reac_rate_local(142) 
  reac_source_local(21,143) = - reac_rate_local(143) 
  reac_source_local(37,143) = + reac_rate_local(143) 
  reac_source_local(53,143) = - reac_rate_local(143) 
  reac_source_local(21,144) = + reac_rate_local(144) 
  reac_source_local(29,144) = - reac_rate_local(144) 
  reac_source_local(36,144) = - reac_rate_local(144) 
  reac_source_local(53,144) = + reac_rate_local(144) 
  reac_source_local(14,145) = - reac_rate_local(145) 
  reac_source_local(36,145) = - reac_rate_local(145) 
  reac_source_local(40,145) = + reac_rate_local(145) 
  reac_source_local(53,145) = + reac_rate_local(145) 
  reac_source_local(36,146) = - reac_rate_local(146) 
  reac_source_local(40,146) = - reac_rate_local(146) 
  reac_source_local(42,146) = + reac_rate_local(146) 
  reac_source_local(53,146) = + reac_rate_local(146) 
  reac_source_local(01,147) = - reac_rate_local(147) 
  reac_source_local(36,147) = - reac_rate_local(147) 
  reac_source_local(41,147) = + reac_rate_local(147) 
  reac_source_local(53,147) = + reac_rate_local(147) 
  reac_source_local(21,148) = - reac_rate_local(148) 
  reac_source_local(32,148) = + reac_rate_local(148) 
  reac_source_local(36,148) = - reac_rate_local(148) 
  reac_source_local(53,148) = + reac_rate_local(148) 
  reac_source_local(26,149) = - reac_rate_local(149) 
  reac_source_local(32,149) = + reac_rate_local(149) 
  reac_source_local(36,149) = - reac_rate_local(149) 
  reac_source_local(53,149) = + reac_rate_local(149) 
  reac_source_local(21,150) = + reac_rate_local(150) 
  reac_source_local(27,150) = - reac_rate_local(150) 
  reac_source_local(29,150) = + reac_rate_local(150) 
  reac_source_local(36,150) = - reac_rate_local(150) 
  reac_source_local(53,150) = + reac_rate_local(150) 
  reac_source_local(01,151) = + reac_rate_local(151) 
  reac_source_local(10,151) = - reac_rate_local(151) 
  reac_source_local(29,151) = + reac_rate_local(151) 
  reac_source_local(36,151) = - reac_rate_local(151) 
  reac_source_local(53,151) = + reac_rate_local(151) 
  reac_source_local(01,152) = + reac_rate_local(152) 
  reac_source_local(11,152) = - reac_rate_local(152) 
  reac_source_local(29,152) = + reac_rate_local(152) 
  reac_source_local(36,152) = - reac_rate_local(152) 
  reac_source_local(53,152) = + reac_rate_local(152) 
  reac_source_local(21,153) = + reac_rate_local(153) * 2.d0
  reac_source_local(32,153) = - reac_rate_local(153) 
  reac_source_local(36,153) = - reac_rate_local(153) 
  reac_source_local(53,153) = + reac_rate_local(153) 
  reac_source_local(29,154) = - reac_rate_local(154) 
  reac_source_local(32,154) = + reac_rate_local(154) 
  reac_source_local(37,154) = - reac_rate_local(154) 
  reac_source_local(53,154) = + reac_rate_local(154) 
  reac_source_local(14,155) = - reac_rate_local(155) 
  reac_source_local(37,155) = - reac_rate_local(155) 
  reac_source_local(42,155) = + reac_rate_local(155) 
  reac_source_local(53,155) = + reac_rate_local(155) 
  reac_source_local(21,156) = + reac_rate_local(156) 
  reac_source_local(37,156) = - reac_rate_local(156) 
  reac_source_local(53,156) = + reac_rate_local(156) 
  reac_source_local(21,157) = + reac_rate_local(157) * 2.d0
  reac_source_local(26,157) = - reac_rate_local(157) 
  reac_source_local(37,157) = - reac_rate_local(157) 
  reac_source_local(53,157) = + reac_rate_local(157) 
  reac_source_local(21,158) = + reac_rate_local(158) * 2.d0
  reac_source_local(27,158) = - reac_rate_local(158) 
  reac_source_local(37,158) = - reac_rate_local(158) 
  reac_source_local(53,158) = + reac_rate_local(158) 
  reac_source_local(21,159) = + reac_rate_local(159) 
  reac_source_local(37,159) = - reac_rate_local(159) 
  reac_source_local(53,159) = + reac_rate_local(159) 
  reac_source_local(01,160) = + reac_rate_local(160) 
  reac_source_local(10,160) = - reac_rate_local(160) 
  reac_source_local(21,160) = + reac_rate_local(160) 
  reac_source_local(37,160) = - reac_rate_local(160) 
  reac_source_local(53,160) = + reac_rate_local(160) 
  reac_source_local(01,161) = + reac_rate_local(161) 
  reac_source_local(11,161) = - reac_rate_local(161) 
  reac_source_local(21,161) = + reac_rate_local(161) 
  reac_source_local(37,161) = - reac_rate_local(161) 
  reac_source_local(53,161) = + reac_rate_local(161) 
  reac_source_local(21,162) = + reac_rate_local(162) * 2.d0
  reac_source_local(29,162) = - reac_rate_local(162) 
  reac_source_local(38,162) = - reac_rate_local(162) 
  reac_source_local(53,162) = + reac_rate_local(162) 
  reac_source_local(14,163) = - reac_rate_local(163) 
  reac_source_local(41,163) = + reac_rate_local(163) 
  reac_source_local(48,163) = - reac_rate_local(163) 
  reac_source_local(53,163) = + reac_rate_local(163) 
  reac_source_local(14,164) = - reac_rate_local(164) 
  reac_source_local(21,164) = + reac_rate_local(164) 
  reac_source_local(38,164) = - reac_rate_local(164) 
  reac_source_local(40,164) = + reac_rate_local(164) 
  reac_source_local(53,164) = + reac_rate_local(164) 
  reac_source_local(01,165) = + reac_rate_local(165) 
  reac_source_local(14,165) = - reac_rate_local(165) 
  reac_source_local(40,165) = + reac_rate_local(165) 
  reac_source_local(49,165) = - reac_rate_local(165) 
  reac_source_local(53,165) = + reac_rate_local(165) 
  reac_source_local(14,166) = - reac_rate_local(166) 
  reac_source_local(40,166) = + reac_rate_local(166) * 2.d0
  reac_source_local(50,166) = - reac_rate_local(166) 
  reac_source_local(53,166) = + reac_rate_local(166) 
  reac_source_local(14,167) = - reac_rate_local(167) 
  reac_source_local(40,167) = + reac_rate_local(167) 
  reac_source_local(42,167) = + reac_rate_local(167) 
  reac_source_local(51,167) = - reac_rate_local(167) 
  reac_source_local(53,167) = + reac_rate_local(167) 
  reac_source_local(29,168) = - reac_rate_local(168) 
  reac_source_local(42,168) = + reac_rate_local(168) 
  reac_source_local(48,168) = - reac_rate_local(168) 
  reac_source_local(53,168) = + reac_rate_local(168) 
  reac_source_local(29,169) = - reac_rate_local(169) 
  reac_source_local(40,169) = + reac_rate_local(169) * 2.d0
  reac_source_local(49,169) = - reac_rate_local(169) 
  reac_source_local(53,169) = + reac_rate_local(169) 
  reac_source_local(21,170) = + reac_rate_local(170) 
  reac_source_local(29,170) = - reac_rate_local(170) 
  reac_source_local(40,170) = + reac_rate_local(170) 
  reac_source_local(50,170) = - reac_rate_local(170) 
  reac_source_local(53,170) = + reac_rate_local(170) 
  reac_source_local(29,171) = - reac_rate_local(171) 
  reac_source_local(32,171) = + reac_rate_local(171) 
  reac_source_local(40,171) = + reac_rate_local(171) 
  reac_source_local(51,171) = - reac_rate_local(171) 
  reac_source_local(53,171) = + reac_rate_local(171) 
  reac_source_local(01,172) = + reac_rate_local(172) 
  reac_source_local(10,172) = - reac_rate_local(172) 
  reac_source_local(32,172) = + reac_rate_local(172) 
  reac_source_local(38,172) = - reac_rate_local(172) 
  reac_source_local(53,172) = + reac_rate_local(172) 
  reac_source_local(01,173) = + reac_rate_local(173) 
  reac_source_local(10,173) = - reac_rate_local(173) 
  reac_source_local(40,173) = + reac_rate_local(173) 
  reac_source_local(48,173) = - reac_rate_local(173) 
  reac_source_local(53,173) = + reac_rate_local(173) 
  reac_source_local(01,174) = + reac_rate_local(174) 
  reac_source_local(10,174) = - reac_rate_local(174) 
  reac_source_local(41,174) = + reac_rate_local(174) 
  reac_source_local(49,174) = - reac_rate_local(174) 
  reac_source_local(53,174) = + reac_rate_local(174) 
  reac_source_local(01,175) = + reac_rate_local(175) 
  reac_source_local(10,175) = - reac_rate_local(175) 
  reac_source_local(42,175) = + reac_rate_local(175) 
  reac_source_local(50,175) = - reac_rate_local(175) 
  reac_source_local(53,175) = + reac_rate_local(175) 
  reac_source_local(01,176) = + reac_rate_local(176) 
  reac_source_local(10,176) = - reac_rate_local(176) 
  reac_source_local(43,176) = + reac_rate_local(176) 
  reac_source_local(51,176) = - reac_rate_local(176) 
  reac_source_local(53,176) = + reac_rate_local(176) 
  reac_source_local(01,177) = + reac_rate_local(177) 
  reac_source_local(11,177) = - reac_rate_local(177) 
  reac_source_local(32,177) = + reac_rate_local(177) 
  reac_source_local(38,177) = - reac_rate_local(177) 
  reac_source_local(53,177) = + reac_rate_local(177) 
  reac_source_local(01,178) = + reac_rate_local(178) 
  reac_source_local(11,178) = - reac_rate_local(178) 
  reac_source_local(40,178) = + reac_rate_local(178) 
  reac_source_local(48,178) = - reac_rate_local(178) 
  reac_source_local(53,178) = + reac_rate_local(178) 
  reac_source_local(01,179) = + reac_rate_local(179) 
  reac_source_local(11,179) = - reac_rate_local(179) 
  reac_source_local(41,179) = + reac_rate_local(179) 
  reac_source_local(49,179) = - reac_rate_local(179) 
  reac_source_local(53,179) = + reac_rate_local(179) 
  reac_source_local(01,180) = + reac_rate_local(180) 
  reac_source_local(11,180) = - reac_rate_local(180) 
  reac_source_local(42,180) = + reac_rate_local(180) 
  reac_source_local(50,180) = - reac_rate_local(180) 
  reac_source_local(53,180) = + reac_rate_local(180) 
  reac_source_local(01,181) = + reac_rate_local(181) 
  reac_source_local(11,181) = - reac_rate_local(181) 
  reac_source_local(43,181) = + reac_rate_local(181) 
  reac_source_local(51,181) = - reac_rate_local(181) 
  reac_source_local(53,181) = + reac_rate_local(181) 
  reac_source_local(01,182) = + reac_rate_local(182) 
  reac_source_local(10,182) = - reac_rate_local(182) 
  reac_source_local(10,183) = + reac_rate_local(183) 
  reac_source_local(11,183) = - reac_rate_local(183) 
  reac_source_local(01,184) = + reac_rate_local(184) 
  reac_source_local(12,184) = - reac_rate_local(184) 
  reac_source_local(11,185) = + reac_rate_local(185) 
  reac_source_local(13,185) = - reac_rate_local(185) 
  reac_source_local(21,186) = + reac_rate_local(186) 
  reac_source_local(26,186) = - reac_rate_local(186) 
  reac_source_local(26,187) = + reac_rate_local(187) 
  reac_source_local(27,187) = - reac_rate_local(187) 
  reac_source_local(21,188) = + reac_rate_local(188) 
  reac_source_local(27,188) = - reac_rate_local(188) 
  reac_source_local(21,189) = + reac_rate_local(189) 
  reac_source_local(28,189) = - reac_rate_local(189) 
  reac_source_local(10,190) = - reac_rate_local(190) 
  reac_source_local(15,190) = + reac_rate_local(190) 
  reac_source_local(29,190) = - reac_rate_local(190) 
  reac_source_local(40,190) = + reac_rate_local(190) 
  reac_source_local(01,191) = + reac_rate_local(191) 
  reac_source_local(10,191) = - reac_rate_local(191) 
  reac_source_local(29,191) = - reac_rate_local(191) 
  reac_source_local(31,191) = + reac_rate_local(191) 
  reac_source_local(01,192) = + reac_rate_local(192) 
  reac_source_local(10,192) = - reac_rate_local(192) 
  reac_source_local(01,193) = + reac_rate_local(193) 
  reac_source_local(10,193) = - reac_rate_local(193) 
  reac_source_local(14,193) = - reac_rate_local(193) 
  reac_source_local(16,193) = + reac_rate_local(193) 
  reac_source_local(01,194) = + reac_rate_local(194) 
  reac_source_local(10,194) = - reac_rate_local(194) 
  reac_source_local(21,194) = - reac_rate_local(194) 
  reac_source_local(29,194) = + reac_rate_local(194) 
  reac_source_local(30,194) = + reac_rate_local(194) 
  reac_source_local(01,195) = + reac_rate_local(195) 
  reac_source_local(10,195) = - reac_rate_local(195) 
  reac_source_local(21,195) = - reac_rate_local(195) 
  reac_source_local(26,195) = + reac_rate_local(195) 
  reac_source_local(01,196) = + reac_rate_local(196) 
  reac_source_local(10,196) = - reac_rate_local(196) 
  reac_source_local(21,196) = - reac_rate_local(196) 
  reac_source_local(27,196) = + reac_rate_local(196) 
  reac_source_local(10,197) = - reac_rate_local(197) 
  reac_source_local(21,197) = - reac_rate_local(197) 
  reac_source_local(29,197) = + reac_rate_local(197) 
  reac_source_local(41,197) = + reac_rate_local(197) 
  reac_source_local(01,198) = + reac_rate_local(198) 
  reac_source_local(10,198) = - reac_rate_local(198) 
  reac_source_local(01,199) = + reac_rate_local(199) 
  reac_source_local(10,199) = - reac_rate_local(199) 
  reac_source_local(01,200) = + reac_rate_local(200) 
  reac_source_local(10,200) = - reac_rate_local(200) 
  reac_source_local(14,200) = + reac_rate_local(200) 
  reac_source_local(40,200) = + reac_rate_local(200) 
  reac_source_local(41,200) = - reac_rate_local(200) 
  reac_source_local(01,201) = + reac_rate_local(201) 
  reac_source_local(10,201) = - reac_rate_local(201) 
  reac_source_local(29,201) = + reac_rate_local(201) 
  reac_source_local(40,201) = + reac_rate_local(201) 
  reac_source_local(42,201) = - reac_rate_local(201) 
  reac_source_local(01,202) = + reac_rate_local(202) 
  reac_source_local(10,202) = - reac_rate_local(202) * 2.d0
  reac_source_local(11,202) = + reac_rate_local(202) 
  reac_source_local(01,203) = + reac_rate_local(203) 
  reac_source_local(10,203) = - reac_rate_local(203) * 2.d0
  reac_source_local(13,203) = + reac_rate_local(203) 
  reac_source_local(10,204) = + reac_rate_local(204) 
  reac_source_local(11,204) = - reac_rate_local(204) 
  reac_source_local(01,205) = + reac_rate_local(205) 
  reac_source_local(11,205) = - reac_rate_local(205) 
  reac_source_local(01,206) = + reac_rate_local(206) 
  reac_source_local(11,206) = - reac_rate_local(206) 
  reac_source_local(21,206) = - reac_rate_local(206) 
  reac_source_local(29,206) = + reac_rate_local(206) * 2.d0
  reac_source_local(10,207) = + reac_rate_local(207) 
  reac_source_local(11,207) = - reac_rate_local(207) 
  reac_source_local(12,208) = + reac_rate_local(208) 
  reac_source_local(13,208) = - reac_rate_local(208) 
  reac_source_local(01,209) = + reac_rate_local(209) 
  reac_source_local(13,209) = - reac_rate_local(209) 
  reac_source_local(21,209) = - reac_rate_local(209) 
  reac_source_local(29,209) = + reac_rate_local(209) 
  reac_source_local(31,209) = + reac_rate_local(209) 
  reac_source_local(11,210) = + reac_rate_local(210) 
  reac_source_local(12,210) = - reac_rate_local(210) 
  reac_source_local(01,211) = + reac_rate_local(211) 
  reac_source_local(12,211) = - reac_rate_local(211) 
  reac_source_local(21,211) = - reac_rate_local(211) 
  reac_source_local(29,211) = + reac_rate_local(211) * 2.d0
  reac_source_local(01,212) = + reac_rate_local(212) 
  reac_source_local(12,212) = - reac_rate_local(212) 
  reac_source_local(14,212) = + reac_rate_local(212) 
  reac_source_local(29,212) = + reac_rate_local(212) 
  reac_source_local(40,212) = - reac_rate_local(212) 
  reac_source_local(10,213) = - reac_rate_local(213) 
  reac_source_local(12,213) = - reac_rate_local(213) 
  reac_source_local(20,213) = + reac_rate_local(213) 
  reac_source_local(53,213) = + reac_rate_local(213) 
  reac_source_local(12,214) = - reac_rate_local(214) * 2.d0
  reac_source_local(20,214) = + reac_rate_local(214) 
  reac_source_local(53,214) = + reac_rate_local(214) 
  reac_source_local(10,215) = + reac_rate_local(215) 
  reac_source_local(14,215) = - reac_rate_local(215) * 2.d0
  reac_source_local(10,216) = + reac_rate_local(216) 
  reac_source_local(14,216) = - reac_rate_local(216) * 2.d0
  reac_source_local(10,217) = + reac_rate_local(217) 
  reac_source_local(14,217) = - reac_rate_local(217) * 2.d0
  reac_source_local(10,218) = + reac_rate_local(218) 
  reac_source_local(14,218) = - reac_rate_local(218) * 2.d0
  reac_source_local(10,219) = + reac_rate_local(219) 
  reac_source_local(14,219) = - reac_rate_local(219) * 2.d0
  reac_source_local(11,220) = + reac_rate_local(220) 
  reac_source_local(14,220) = - reac_rate_local(220) * 2.d0
  reac_source_local(11,221) = + reac_rate_local(221) 
  reac_source_local(14,221) = - reac_rate_local(221) * 2.d0
  reac_source_local(11,222) = + reac_rate_local(222) 
  reac_source_local(14,222) = - reac_rate_local(222) * 2.d0
  reac_source_local(11,223) = + reac_rate_local(223) 
  reac_source_local(14,223) = - reac_rate_local(223) * 2.d0
  reac_source_local(11,224) = + reac_rate_local(224) 
  reac_source_local(14,224) = - reac_rate_local(224) * 2.d0
  reac_source_local(14,225) = + reac_rate_local(225) 
  reac_source_local(15,225) = - reac_rate_local(225) 
  reac_source_local(29,225) = - reac_rate_local(225) 
  reac_source_local(30,225) = + reac_rate_local(225) 
  reac_source_local(15,226) = - reac_rate_local(226) 
  reac_source_local(21,226) = - reac_rate_local(226) 
  reac_source_local(29,226) = + reac_rate_local(226) 
  reac_source_local(40,226) = + reac_rate_local(226) 
  reac_source_local(01,227) = + reac_rate_local(227) 
  reac_source_local(15,227) = - reac_rate_local(227) 
  reac_source_local(29,227) = + reac_rate_local(227) 
  reac_source_local(40,227) = - reac_rate_local(227) 
  reac_source_local(01,228) = + reac_rate_local(228) 
  reac_source_local(15,228) = - reac_rate_local(228) 
  reac_source_local(40,228) = + reac_rate_local(228) 
  reac_source_local(41,228) = - reac_rate_local(228) 
  reac_source_local(14,229) = + reac_rate_local(229) 
  reac_source_local(15,229) = - reac_rate_local(229) 
  reac_source_local(14,230) = + reac_rate_local(230) 
  reac_source_local(16,230) = - reac_rate_local(230) 
  reac_source_local(14,231) = + reac_rate_local(231) 
  reac_source_local(16,231) = - reac_rate_local(231) 
  reac_source_local(15,232) = + reac_rate_local(232) 
  reac_source_local(16,232) = - reac_rate_local(232) 
  reac_source_local(14,233) = + reac_rate_local(233) 
  reac_source_local(16,233) = - reac_rate_local(233) 
  reac_source_local(15,234) = - reac_rate_local(234) 
  reac_source_local(16,234) = - reac_rate_local(234) 
  reac_source_local(18,234) = + reac_rate_local(234) 
  reac_source_local(53,234) = + reac_rate_local(234) 
  reac_source_local(16,235) = - reac_rate_local(235) 
  reac_source_local(21,235) = - reac_rate_local(235) 
  reac_source_local(29,235) = + reac_rate_local(235) 
  reac_source_local(40,235) = + reac_rate_local(235) 
  reac_source_local(10,236) = + reac_rate_local(236) 
  reac_source_local(16,236) = - reac_rate_local(236) 
  reac_source_local(29,236) = + reac_rate_local(236) 
  reac_source_local(40,236) = - reac_rate_local(236) 
  reac_source_local(21,237) = + reac_rate_local(237) 
  reac_source_local(26,237) = - reac_rate_local(237) 
  reac_source_local(14,238) = - reac_rate_local(238) 
  reac_source_local(26,238) = - reac_rate_local(238) 
  reac_source_local(29,238) = + reac_rate_local(238) 
  reac_source_local(40,238) = + reac_rate_local(238) 
  reac_source_local(21,239) = + reac_rate_local(239) 
  reac_source_local(26,239) = - reac_rate_local(239) 
  reac_source_local(21,240) = + reac_rate_local(240) 
  reac_source_local(26,240) = - reac_rate_local(240) 
  reac_source_local(21,241) = + reac_rate_local(241) 
  reac_source_local(26,241) = - reac_rate_local(241) 
  reac_source_local(21,242) = + reac_rate_local(242) * 2.d0
  reac_source_local(26,242) = - reac_rate_local(242) 
  reac_source_local(30,242) = + reac_rate_local(242) 
  reac_source_local(32,242) = - reac_rate_local(242) 
  reac_source_local(21,243) = + reac_rate_local(243) 
  reac_source_local(26,243) = - reac_rate_local(243) * 2.d0
  reac_source_local(27,243) = + reac_rate_local(243) 
  reac_source_local(21,244) = + reac_rate_local(244) 
  reac_source_local(26,244) = + reac_rate_local(244) 
  reac_source_local(29,244) = - reac_rate_local(244) 
  reac_source_local(32,244) = - reac_rate_local(244) 
  reac_source_local(26,245) = + reac_rate_local(245) 
  reac_source_local(27,245) = - reac_rate_local(245) 
  reac_source_local(21,246) = + reac_rate_local(246) 
  reac_source_local(27,246) = - reac_rate_local(246) 
  reac_source_local(29,246) = - reac_rate_local(246) 
  reac_source_local(30,246) = + reac_rate_local(246) 
  reac_source_local(26,247) = + reac_rate_local(247) 
  reac_source_local(27,247) = - reac_rate_local(247) 
  reac_source_local(26,248) = + reac_rate_local(248) 
  reac_source_local(27,248) = - reac_rate_local(248) 
  reac_source_local(26,249) = + reac_rate_local(249) 
  reac_source_local(27,249) = - reac_rate_local(249) 
  reac_source_local(21,250) = + reac_rate_local(250) * 2.d0
  reac_source_local(27,250) = - reac_rate_local(250) 
  reac_source_local(29,250) = + reac_rate_local(250) 
  reac_source_local(32,250) = - reac_rate_local(250) 
  reac_source_local(21,251) = + reac_rate_local(251) 
  reac_source_local(28,251) = - reac_rate_local(251) 
  reac_source_local(29,251) = - reac_rate_local(251) 
  reac_source_local(31,251) = + reac_rate_local(251) 
  reac_source_local(21,252) = - reac_rate_local(252) 
  reac_source_local(27,252) = + reac_rate_local(252) * 2.d0
  reac_source_local(28,252) = - reac_rate_local(252) 
  reac_source_local(27,253) = + reac_rate_local(253) 
  reac_source_local(28,253) = - reac_rate_local(253) 
  reac_source_local(29,254) = + reac_rate_local(254) 
  reac_source_local(30,254) = - reac_rate_local(254) 
  reac_source_local(29,255) = + reac_rate_local(255) 
  reac_source_local(30,255) = - reac_rate_local(255) 
  reac_source_local(21,256) = - reac_rate_local(256) 
  reac_source_local(26,256) = + reac_rate_local(256) 
  reac_source_local(29,256) = + reac_rate_local(256) 
  reac_source_local(30,256) = - reac_rate_local(256) 
  reac_source_local(21,257) = - reac_rate_local(257) 
  reac_source_local(27,257) = + reac_rate_local(257) 
  reac_source_local(29,257) = + reac_rate_local(257) 
  reac_source_local(30,257) = - reac_rate_local(257) 
  reac_source_local(29,258) = + reac_rate_local(258) 
  reac_source_local(30,258) = - reac_rate_local(258) 
  reac_source_local(21,259) = + reac_rate_local(259) 
  reac_source_local(29,259) = + reac_rate_local(259) * 2.d0
  reac_source_local(30,259) = - reac_rate_local(259) 
  reac_source_local(32,259) = - reac_rate_local(259) 
  reac_source_local(21,260) = + reac_rate_local(260) * 2.d0
  reac_source_local(30,260) = - reac_rate_local(260) 
  reac_source_local(32,260) = - reac_rate_local(260) 
  reac_source_local(14,261) = + reac_rate_local(261) 
  reac_source_local(21,261) = + reac_rate_local(261) 
  reac_source_local(30,261) = - reac_rate_local(261) 
  reac_source_local(40,261) = - reac_rate_local(261) 
  reac_source_local(30,262) = - reac_rate_local(262) 
  reac_source_local(40,262) = + reac_rate_local(262) * 2.d0
  reac_source_local(41,262) = - reac_rate_local(262) 
  reac_source_local(01,263) = + reac_rate_local(263) 
  reac_source_local(21,263) = + reac_rate_local(263) 
  reac_source_local(30,263) = - reac_rate_local(263) 
  reac_source_local(41,263) = - reac_rate_local(263) 
  reac_source_local(30,264) = + reac_rate_local(264) 
  reac_source_local(31,264) = - reac_rate_local(264) 
  reac_source_local(29,265) = + reac_rate_local(265) 
  reac_source_local(31,265) = - reac_rate_local(265) 
  reac_source_local(30,266) = + reac_rate_local(266) 
  reac_source_local(31,266) = - reac_rate_local(266) 
  reac_source_local(21,267) = - reac_rate_local(267) 
  reac_source_local(29,267) = + reac_rate_local(267) * 3.d0
  reac_source_local(31,267) = - reac_rate_local(267) 
  reac_source_local(29,268) = + reac_rate_local(268) 
  reac_source_local(31,268) = - reac_rate_local(268) 
  reac_source_local(26,269) = - reac_rate_local(269) 
  reac_source_local(28,269) = + reac_rate_local(269) 
  reac_source_local(29,269) = + reac_rate_local(269) 
  reac_source_local(31,269) = - reac_rate_local(269) 
  reac_source_local(26,270) = - reac_rate_local(270) 
  reac_source_local(27,270) = + reac_rate_local(270) 
  reac_source_local(30,270) = + reac_rate_local(270) 
  reac_source_local(31,270) = - reac_rate_local(270) 
  reac_source_local(26,271) = - reac_rate_local(271) 
  reac_source_local(29,271) = + reac_rate_local(271) * 3.d0
  reac_source_local(31,271) = - reac_rate_local(271) 
  reac_source_local(29,272) = + reac_rate_local(272) 
  reac_source_local(31,272) = - reac_rate_local(272) 
  reac_source_local(30,273) = + reac_rate_local(273) 
  reac_source_local(31,273) = - reac_rate_local(273) 
  reac_source_local(21,274) = + reac_rate_local(274) * 2.d0
  reac_source_local(31,274) = - reac_rate_local(274) 
  reac_source_local(32,274) = - reac_rate_local(274) 
  reac_source_local(21,275) = + reac_rate_local(275) 
  reac_source_local(29,275) = + reac_rate_local(275) 
  reac_source_local(30,275) = + reac_rate_local(275) 
  reac_source_local(31,275) = - reac_rate_local(275) 
  reac_source_local(32,275) = - reac_rate_local(275) 
  reac_source_local(29,276) = + reac_rate_local(276) 
  reac_source_local(31,276) = - reac_rate_local(276) 
  reac_source_local(30,277) = + reac_rate_local(277) 
  reac_source_local(31,277) = - reac_rate_local(277) 
  reac_source_local(01,278) = + reac_rate_local(278) 
  reac_source_local(14,278) = - reac_rate_local(278) 
  reac_source_local(29,278) = + reac_rate_local(278) 
  reac_source_local(40,278) = - reac_rate_local(278) 
  reac_source_local(14,279) = - reac_rate_local(279) 
  reac_source_local(21,279) = - reac_rate_local(279) 
  reac_source_local(29,279) = + reac_rate_local(279) 
  reac_source_local(40,279) = + reac_rate_local(279) 
  reac_source_local(01,280) = + reac_rate_local(280) 
  reac_source_local(14,280) = - reac_rate_local(280) 
  reac_source_local(29,280) = + reac_rate_local(280) * 2.d0
  reac_source_local(42,280) = - reac_rate_local(280) 
  reac_source_local(14,281) = - reac_rate_local(281) 
  reac_source_local(29,281) = + reac_rate_local(281) 
  reac_source_local(41,281) = + reac_rate_local(281) 
  reac_source_local(42,281) = - reac_rate_local(281) 
  reac_source_local(01,282) = + reac_rate_local(282) 
  reac_source_local(14,282) = - reac_rate_local(282) 
  reac_source_local(21,282) = + reac_rate_local(282) 
  reac_source_local(42,282) = - reac_rate_local(282) 
  reac_source_local(14,283) = - reac_rate_local(283) 
  reac_source_local(40,283) = + reac_rate_local(283) * 2.d0
  reac_source_local(42,283) = - reac_rate_local(283) 
  reac_source_local(01,284) = - reac_rate_local(284) 
  reac_source_local(14,284) = + reac_rate_local(284) 
  reac_source_local(29,284) = - reac_rate_local(284) 
  reac_source_local(40,284) = + reac_rate_local(284) 
  reac_source_local(14,285) = + reac_rate_local(285) 
  reac_source_local(21,285) = + reac_rate_local(285) 
  reac_source_local(29,285) = - reac_rate_local(285) 
  reac_source_local(40,285) = - reac_rate_local(285) 
  reac_source_local(29,286) = - reac_rate_local(286) 
  reac_source_local(40,286) = - reac_rate_local(286) 
  reac_source_local(42,286) = + reac_rate_local(286) 
  reac_source_local(01,287) = + reac_rate_local(287) 
  reac_source_local(21,287) = + reac_rate_local(287) 
  reac_source_local(29,287) = - reac_rate_local(287) 
  reac_source_local(41,287) = - reac_rate_local(287) 
  reac_source_local(29,288) = - reac_rate_local(288) 
  reac_source_local(40,288) = + reac_rate_local(288) * 2.d0
  reac_source_local(41,288) = - reac_rate_local(288) 
  reac_source_local(21,289) = + reac_rate_local(289) 
  reac_source_local(29,289) = - reac_rate_local(289) 
  reac_source_local(40,289) = + reac_rate_local(289) 
  reac_source_local(42,289) = - reac_rate_local(289) 
  reac_source_local(21,290) = + reac_rate_local(290) 
  reac_source_local(29,290) = - reac_rate_local(290) 
  reac_source_local(42,290) = + reac_rate_local(290) 
  reac_source_local(43,290) = - reac_rate_local(290) 
  reac_source_local(01,291) = - reac_rate_local(291) 
  reac_source_local(21,291) = - reac_rate_local(291) 
  reac_source_local(29,291) = + reac_rate_local(291) 
  reac_source_local(41,291) = + reac_rate_local(291) 
  reac_source_local(14,292) = + reac_rate_local(292) 
  reac_source_local(40,292) = - reac_rate_local(292) * 2.d0
  reac_source_local(42,292) = + reac_rate_local(292) 
  reac_source_local(29,293) = + reac_rate_local(293) 
  reac_source_local(40,293) = - reac_rate_local(293) * 2.d0
  reac_source_local(41,293) = + reac_rate_local(293) 
  reac_source_local(01,294) = + reac_rate_local(294) 
  reac_source_local(21,294) = + reac_rate_local(294) 
  reac_source_local(40,294) = - reac_rate_local(294) * 2.d0
  reac_source_local(21,295) = - reac_rate_local(295) 
  reac_source_local(29,295) = + reac_rate_local(295) 
  reac_source_local(40,295) = - reac_rate_local(295) 
  reac_source_local(42,295) = + reac_rate_local(295) 
  reac_source_local(21,296) = + reac_rate_local(296) 
  reac_source_local(32,296) = - reac_rate_local(296) 
  reac_source_local(40,296) = - reac_rate_local(296) 
  reac_source_local(42,296) = + reac_rate_local(296) 
  reac_source_local(01,297) = + reac_rate_local(297) 
  reac_source_local(40,297) = - reac_rate_local(297) 
  reac_source_local(41,297) = - reac_rate_local(297) 
  reac_source_local(42,297) = + reac_rate_local(297) 
  reac_source_local(40,298) = - reac_rate_local(298) 
  reac_source_local(42,298) = + reac_rate_local(298) * 2.d0
  reac_source_local(43,298) = - reac_rate_local(298) 
  reac_source_local(21,299) = - reac_rate_local(299) * 2.d0
  reac_source_local(29,299) = + reac_rate_local(299) 
  reac_source_local(32,299) = + reac_rate_local(299) 
  reac_source_local(21,300) = - reac_rate_local(300) 
  reac_source_local(32,300) = + reac_rate_local(300) 
  reac_source_local(40,300) = + reac_rate_local(300) 
  reac_source_local(42,300) = - reac_rate_local(300) 
  reac_source_local(21,301) = + reac_rate_local(301) 
  reac_source_local(40,301) = + reac_rate_local(301) * 2.d0
  reac_source_local(42,301) = - reac_rate_local(301) * 2.d0
  reac_source_local(40,302) = + reac_rate_local(302) 
  reac_source_local(42,302) = - reac_rate_local(302) * 2.d0
  reac_source_local(43,302) = + reac_rate_local(302) 
  reac_source_local(21,303) = + reac_rate_local(303) 
  reac_source_local(32,303) = - reac_rate_local(303) 
  reac_source_local(42,303) = - reac_rate_local(303) 
  reac_source_local(43,303) = + reac_rate_local(303) 
  reac_source_local(21,304) = + reac_rate_local(304) 
  reac_source_local(40,304) = + reac_rate_local(304) 
  reac_source_local(43,304) = - reac_rate_local(304) 
  reac_source_local(21,305) = - reac_rate_local(305) 
  reac_source_local(32,305) = + reac_rate_local(305) 
  reac_source_local(42,305) = + reac_rate_local(305) 
  reac_source_local(43,305) = - reac_rate_local(305) 
  reac_source_local(21,306) = + reac_rate_local(306) 
  reac_source_local(42,306) = + reac_rate_local(306) * 2.d0
  reac_source_local(43,306) = - reac_rate_local(306) * 2.d0
  reac_source_local(14,307) = - reac_rate_local(307) * 2.d0
  reac_source_local(18,307) = + reac_rate_local(307) 
  reac_source_local(53,307) = + reac_rate_local(307) 
  reac_source_local(14,308) = - reac_rate_local(308) 
  reac_source_local(29,308) = - reac_rate_local(308) 
  reac_source_local(45,308) = + reac_rate_local(308) 
  reac_source_local(53,308) = + reac_rate_local(308) 
  reac_source_local(01,309) = - reac_rate_local(309) 
  reac_source_local(14,309) = + reac_rate_local(309) * 2.d0
  reac_source_local(01,310) = - reac_rate_local(310) 
  reac_source_local(14,310) = + reac_rate_local(310) * 2.d0
  reac_source_local(01,311) = - reac_rate_local(311) 
  reac_source_local(14,311) = + reac_rate_local(311) * 2.d0
  reac_source_local(01,312) = - reac_rate_local(312) 
  reac_source_local(14,312) = + reac_rate_local(312) * 2.d0
  reac_source_local(01,313) = - reac_rate_local(313) 
  reac_source_local(14,313) = + reac_rate_local(313) * 2.d0
  reac_source_local(21,314) = - reac_rate_local(314) 
  reac_source_local(29,314) = + reac_rate_local(314) * 2.d0
  reac_source_local(21,315) = - reac_rate_local(315) 
  reac_source_local(29,315) = + reac_rate_local(315) * 2.d0
  reac_source_local(21,316) = - reac_rate_local(316) 
  reac_source_local(29,316) = + reac_rate_local(316) * 2.d0
  reac_source_local(21,317) = - reac_rate_local(317) 
  reac_source_local(29,317) = + reac_rate_local(317) * 2.d0
  reac_source_local(21,318) = - reac_rate_local(318) 
  reac_source_local(29,318) = + reac_rate_local(318) * 2.d0
  reac_source_local(14,319) = + reac_rate_local(319) 
  reac_source_local(29,319) = + reac_rate_local(319) 
  reac_source_local(40,319) = - reac_rate_local(319) 
  reac_source_local(14,320) = + reac_rate_local(320) 
  reac_source_local(29,320) = + reac_rate_local(320) 
  reac_source_local(40,320) = - reac_rate_local(320) 
  reac_source_local(14,321) = + reac_rate_local(321) 
  reac_source_local(29,321) = + reac_rate_local(321) 
  reac_source_local(40,321) = - reac_rate_local(321) 
  reac_source_local(14,322) = + reac_rate_local(322) 
  reac_source_local(29,322) = + reac_rate_local(322) 
  reac_source_local(40,322) = - reac_rate_local(322) 
  reac_source_local(14,323) = + reac_rate_local(323) 
  reac_source_local(29,323) = + reac_rate_local(323) 
  reac_source_local(40,323) = - reac_rate_local(323) 
  reac_source_local(21,324) = + reac_rate_local(324) 
  reac_source_local(29,324) = + reac_rate_local(324) 
  reac_source_local(32,324) = - reac_rate_local(324) 
  reac_source_local(21,325) = + reac_rate_local(325) 
  reac_source_local(29,325) = + reac_rate_local(325) 
  reac_source_local(32,325) = - reac_rate_local(325) 
  reac_source_local(21,326) = + reac_rate_local(326) 
  reac_source_local(29,326) = + reac_rate_local(326) 
  reac_source_local(32,326) = - reac_rate_local(326) 
  reac_source_local(21,327) = + reac_rate_local(327) 
  reac_source_local(29,327) = + reac_rate_local(327) 
  reac_source_local(32,327) = - reac_rate_local(327) 
  reac_source_local(01,328) = + reac_rate_local(328) 
  reac_source_local(29,328) = + reac_rate_local(328) 
  reac_source_local(41,328) = - reac_rate_local(328) 
  reac_source_local(01,329) = + reac_rate_local(329) 
  reac_source_local(29,329) = + reac_rate_local(329) 
  reac_source_local(41,329) = - reac_rate_local(329) 
  reac_source_local(01,330) = + reac_rate_local(330) 
  reac_source_local(29,330) = + reac_rate_local(330) 
  reac_source_local(41,330) = - reac_rate_local(330) 
  reac_source_local(01,331) = + reac_rate_local(331) 
  reac_source_local(29,331) = + reac_rate_local(331) 
  reac_source_local(41,331) = - reac_rate_local(331) 
  reac_source_local(29,332) = + reac_rate_local(332) 
  reac_source_local(40,332) = + reac_rate_local(332) 
  reac_source_local(42,332) = - reac_rate_local(332) 
  reac_source_local(29,333) = + reac_rate_local(333) 
  reac_source_local(40,333) = + reac_rate_local(333) 
  reac_source_local(42,333) = - reac_rate_local(333) 
  reac_source_local(29,334) = + reac_rate_local(334) 
  reac_source_local(40,334) = + reac_rate_local(334) 
  reac_source_local(42,334) = - reac_rate_local(334) 
  reac_source_local(29,335) = + reac_rate_local(335) 
  reac_source_local(40,335) = + reac_rate_local(335) 
  reac_source_local(42,335) = - reac_rate_local(335) 
  reac_source_local(29,336) = + reac_rate_local(336) 
  reac_source_local(42,336) = + reac_rate_local(336) 
  reac_source_local(43,336) = - reac_rate_local(336) 
  reac_source_local(29,337) = + reac_rate_local(337) 
  reac_source_local(42,337) = + reac_rate_local(337) 
  reac_source_local(43,337) = - reac_rate_local(337) 
  reac_source_local(29,338) = + reac_rate_local(338) 
  reac_source_local(42,338) = + reac_rate_local(338) 
  reac_source_local(43,338) = - reac_rate_local(338) 
  reac_source_local(29,339) = + reac_rate_local(339) 
  reac_source_local(42,339) = + reac_rate_local(339) 
  reac_source_local(43,339) = - reac_rate_local(339) 
  reac_source_local(29,340) = + reac_rate_local(340) 
  reac_source_local(42,340) = + reac_rate_local(340) 
  reac_source_local(43,340) = - reac_rate_local(340) 
  reac_source_local(21,341) = + reac_rate_local(341) 
  reac_source_local(40,341) = + reac_rate_local(341) 
  reac_source_local(43,341) = - reac_rate_local(341) 
  reac_source_local(21,342) = + reac_rate_local(342) 
  reac_source_local(40,342) = + reac_rate_local(342) 
  reac_source_local(43,342) = - reac_rate_local(342) 
  reac_source_local(21,343) = + reac_rate_local(343) 
  reac_source_local(40,343) = + reac_rate_local(343) 
  reac_source_local(43,343) = - reac_rate_local(343) 
  reac_source_local(21,344) = + reac_rate_local(344) 
  reac_source_local(40,344) = + reac_rate_local(344) 
  reac_source_local(43,344) = - reac_rate_local(344) 
  reac_source_local(21,345) = + reac_rate_local(345) 
  reac_source_local(40,345) = + reac_rate_local(345) 
  reac_source_local(43,345) = - reac_rate_local(345) 
  reac_source_local(42,346) = + reac_rate_local(346) 
  reac_source_local(43,346) = + reac_rate_local(346) 
  reac_source_local(44,346) = - reac_rate_local(346) 
  reac_source_local(01,347) = + reac_rate_local(347) 
  reac_source_local(14,347) = - reac_rate_local(347) * 2.d0
  reac_source_local(01,348) = + reac_rate_local(348) 
  reac_source_local(14,348) = - reac_rate_local(348) * 2.d0
  reac_source_local(01,349) = + reac_rate_local(349) 
  reac_source_local(14,349) = - reac_rate_local(349) * 2.d0
  reac_source_local(01,350) = + reac_rate_local(350) 
  reac_source_local(14,350) = - reac_rate_local(350) * 2.d0
  reac_source_local(01,351) = + reac_rate_local(351) 
  reac_source_local(14,351) = - reac_rate_local(351) * 2.d0
  reac_source_local(21,352) = + reac_rate_local(352) 
  reac_source_local(29,352) = - reac_rate_local(352) * 2.d0
  reac_source_local(21,353) = + reac_rate_local(353) 
  reac_source_local(29,353) = - reac_rate_local(353) * 2.d0
  reac_source_local(21,354) = + reac_rate_local(354) 
  reac_source_local(29,354) = - reac_rate_local(354) * 2.d0
  reac_source_local(21,355) = + reac_rate_local(355) 
  reac_source_local(29,355) = - reac_rate_local(355) * 2.d0
  reac_source_local(21,356) = + reac_rate_local(356) 
  reac_source_local(29,356) = - reac_rate_local(356) * 2.d0
  reac_source_local(14,357) = - reac_rate_local(357) 
  reac_source_local(29,357) = - reac_rate_local(357) 
  reac_source_local(40,357) = + reac_rate_local(357) 
  reac_source_local(14,358) = - reac_rate_local(358) 
  reac_source_local(29,358) = - reac_rate_local(358) 
  reac_source_local(40,358) = + reac_rate_local(358) 
  reac_source_local(14,359) = - reac_rate_local(359) 
  reac_source_local(29,359) = - reac_rate_local(359) 
  reac_source_local(40,359) = + reac_rate_local(359) 
  reac_source_local(14,360) = - reac_rate_local(360) 
  reac_source_local(29,360) = - reac_rate_local(360) 
  reac_source_local(40,360) = + reac_rate_local(360) 
  reac_source_local(14,361) = - reac_rate_local(361) 
  reac_source_local(29,361) = - reac_rate_local(361) 
  reac_source_local(40,361) = + reac_rate_local(361) 
  reac_source_local(21,362) = - reac_rate_local(362) 
  reac_source_local(29,362) = - reac_rate_local(362) 
  reac_source_local(32,362) = + reac_rate_local(362) 
  reac_source_local(21,363) = - reac_rate_local(363) 
  reac_source_local(29,363) = - reac_rate_local(363) 
  reac_source_local(32,363) = + reac_rate_local(363) 
  reac_source_local(21,364) = - reac_rate_local(364) 
  reac_source_local(29,364) = - reac_rate_local(364) 
  reac_source_local(32,364) = + reac_rate_local(364) 
  reac_source_local(21,365) = - reac_rate_local(365) 
  reac_source_local(29,365) = - reac_rate_local(365) 
  reac_source_local(32,365) = + reac_rate_local(365) 
  reac_source_local(21,366) = - reac_rate_local(366) 
  reac_source_local(29,366) = - reac_rate_local(366) 
  reac_source_local(32,366) = + reac_rate_local(366) 
  reac_source_local(01,367) = - reac_rate_local(367) 
  reac_source_local(29,367) = - reac_rate_local(367) 
  reac_source_local(41,367) = + reac_rate_local(367) 
  reac_source_local(29,368) = - reac_rate_local(368) 
  reac_source_local(40,368) = - reac_rate_local(368) 
  reac_source_local(42,368) = + reac_rate_local(368) 
  reac_source_local(29,369) = - reac_rate_local(369) 
  reac_source_local(40,369) = - reac_rate_local(369) 
  reac_source_local(42,369) = + reac_rate_local(369) 
  reac_source_local(29,370) = - reac_rate_local(370) 
  reac_source_local(40,370) = - reac_rate_local(370) 
  reac_source_local(42,370) = + reac_rate_local(370) 
  reac_source_local(29,371) = - reac_rate_local(371) 
  reac_source_local(42,371) = - reac_rate_local(371) 
  reac_source_local(43,371) = + reac_rate_local(371) 
  reac_source_local(29,372) = - reac_rate_local(372) 
  reac_source_local(42,372) = - reac_rate_local(372) 
  reac_source_local(43,372) = + reac_rate_local(372) 
  reac_source_local(29,373) = - reac_rate_local(373) 
  reac_source_local(42,373) = - reac_rate_local(373) 
  reac_source_local(43,373) = + reac_rate_local(373) 
  reac_source_local(29,374) = - reac_rate_local(374) 
  reac_source_local(42,374) = - reac_rate_local(374) 
  reac_source_local(43,374) = + reac_rate_local(374) 
  reac_source_local(29,375) = - reac_rate_local(375) 
  reac_source_local(42,375) = - reac_rate_local(375) 
  reac_source_local(43,375) = + reac_rate_local(375) 
  reac_source_local(42,376) = - reac_rate_local(376) 
  reac_source_local(43,376) = - reac_rate_local(376) 
  reac_source_local(44,376) = + reac_rate_local(376) 
  reac_source_local(14,377) = + reac_rate_local(377) 
  reac_source_local(17,377) = - reac_rate_local(377) 
  reac_source_local(29,377) = - reac_rate_local(377) 
  reac_source_local(33,377) = + reac_rate_local(377) 
  reac_source_local(14,378) = + reac_rate_local(378) 
  reac_source_local(17,378) = - reac_rate_local(378) 
  reac_source_local(21,378) = - reac_rate_local(378) 
  reac_source_local(34,378) = + reac_rate_local(378) 
  reac_source_local(17,379) = - reac_rate_local(379) 
  reac_source_local(21,379) = - reac_rate_local(379) 
  reac_source_local(29,379) = + reac_rate_local(379) 
  reac_source_local(45,379) = + reac_rate_local(379) 
  reac_source_local(17,380) = - reac_rate_local(380) 
  reac_source_local(21,380) = - reac_rate_local(380) 
  reac_source_local(33,380) = + reac_rate_local(380) 
  reac_source_local(40,380) = + reac_rate_local(380) 
  reac_source_local(17,381) = - reac_rate_local(381) 
  reac_source_local(21,381) = + reac_rate_local(381) 
  reac_source_local(32,381) = - reac_rate_local(381) 
  reac_source_local(45,381) = + reac_rate_local(381) 
  reac_source_local(14,382) = + reac_rate_local(382) 
  reac_source_local(17,382) = - reac_rate_local(382) 
  reac_source_local(40,382) = - reac_rate_local(382) 
  reac_source_local(45,382) = + reac_rate_local(382) 
  reac_source_local(17,383) = - reac_rate_local(383) 
  reac_source_local(18,383) = + reac_rate_local(383) 
  reac_source_local(29,383) = + reac_rate_local(383) 
  reac_source_local(40,383) = - reac_rate_local(383) 
  reac_source_local(01,384) = + reac_rate_local(384) 
  reac_source_local(17,384) = - reac_rate_local(384) 
  reac_source_local(33,384) = + reac_rate_local(384) 
  reac_source_local(40,384) = - reac_rate_local(384) 
  reac_source_local(01,385) = + reac_rate_local(385) 
  reac_source_local(17,385) = - reac_rate_local(385) 
  reac_source_local(41,385) = - reac_rate_local(385) 
  reac_source_local(45,385) = + reac_rate_local(385) 
  reac_source_local(01,386) = - reac_rate_local(386) 
  reac_source_local(14,386) = + reac_rate_local(386) 
  reac_source_local(33,386) = - reac_rate_local(386) 
  reac_source_local(45,386) = + reac_rate_local(386) 
  reac_source_local(21,387) = - reac_rate_local(387) 
  reac_source_local(29,387) = + reac_rate_local(387) 
  reac_source_local(33,387) = - reac_rate_local(387) 
  reac_source_local(34,387) = + reac_rate_local(387) 
  reac_source_local(21,388) = + reac_rate_local(388) 
  reac_source_local(32,388) = - reac_rate_local(388) 
  reac_source_local(33,388) = - reac_rate_local(388) 
  reac_source_local(34,388) = + reac_rate_local(388) 
  reac_source_local(29,389) = + reac_rate_local(389) 
  reac_source_local(33,389) = - reac_rate_local(389) 
  reac_source_local(40,389) = - reac_rate_local(389) 
  reac_source_local(45,389) = + reac_rate_local(389) 
  reac_source_local(14,390) = + reac_rate_local(390) 
  reac_source_local(33,390) = - reac_rate_local(390) 
  reac_source_local(34,390) = + reac_rate_local(390) 
  reac_source_local(40,390) = - reac_rate_local(390) 
  reac_source_local(15,391) = - reac_rate_local(391) 
  reac_source_local(17,391) = + reac_rate_local(391) 
  reac_source_local(29,391) = + reac_rate_local(391) 
  reac_source_local(33,391) = - reac_rate_local(391) 
  reac_source_local(33,392) = - reac_rate_local(392) 
  reac_source_local(40,392) = + reac_rate_local(392) 
  reac_source_local(41,392) = - reac_rate_local(392) 
  reac_source_local(45,392) = + reac_rate_local(392) 
  reac_source_local(29,393) = + reac_rate_local(393) 
  reac_source_local(33,393) = - reac_rate_local(393) 
  reac_source_local(41,393) = - reac_rate_local(393) 
  reac_source_local(46,393) = + reac_rate_local(393) 
  reac_source_local(01,394) = + reac_rate_local(394) 
  reac_source_local(33,394) = - reac_rate_local(394) 
  reac_source_local(34,394) = + reac_rate_local(394) 
  reac_source_local(41,394) = - reac_rate_local(394) 
  reac_source_local(29,395) = + reac_rate_local(395) 
  reac_source_local(33,395) = - reac_rate_local(395) 
  reac_source_local(42,395) = - reac_rate_local(395) 
  reac_source_local(47,395) = + reac_rate_local(395) 
  reac_source_local(01,396) = + reac_rate_local(396) 
  reac_source_local(18,396) = - reac_rate_local(396) 
  reac_source_local(21,396) = - reac_rate_local(396) 
  reac_source_local(34,396) = + reac_rate_local(396) 
  reac_source_local(14,397) = + reac_rate_local(397) 
  reac_source_local(18,397) = - reac_rate_local(397) 
  reac_source_local(29,397) = - reac_rate_local(397) 
  reac_source_local(45,397) = + reac_rate_local(397) 
  reac_source_local(01,398) = + reac_rate_local(398) 
  reac_source_local(18,398) = - reac_rate_local(398) 
  reac_source_local(29,398) = + reac_rate_local(398) 
  reac_source_local(32,398) = - reac_rate_local(398) 
  reac_source_local(34,398) = + reac_rate_local(398) 
  reac_source_local(01,399) = + reac_rate_local(399) 
  reac_source_local(14,399) = - reac_rate_local(399) 
  reac_source_local(17,399) = + reac_rate_local(399) 
  reac_source_local(18,399) = - reac_rate_local(399) 
  reac_source_local(01,400) = + reac_rate_local(400) 
  reac_source_local(18,400) = - reac_rate_local(400) 
  reac_source_local(40,400) = - reac_rate_local(400) 
  reac_source_local(45,400) = + reac_rate_local(400) 
  reac_source_local(01,401) = + reac_rate_local(401) 
  reac_source_local(18,401) = - reac_rate_local(401) 
  reac_source_local(41,401) = - reac_rate_local(401) 
  reac_source_local(46,401) = + reac_rate_local(401) 
  reac_source_local(01,402) = + reac_rate_local(402) 
  reac_source_local(14,402) = + reac_rate_local(402) 
  reac_source_local(18,402) = - reac_rate_local(402) 
  reac_source_local(41,402) = - reac_rate_local(402) 
  reac_source_local(45,402) = + reac_rate_local(402) 
  reac_source_local(01,403) = - reac_rate_local(403) 
  reac_source_local(34,403) = - reac_rate_local(403) 
  reac_source_local(40,403) = + reac_rate_local(403) 
  reac_source_local(45,403) = + reac_rate_local(403) 
  reac_source_local(14,404) = - reac_rate_local(404) 
  reac_source_local(29,404) = + reac_rate_local(404) 
  reac_source_local(34,404) = - reac_rate_local(404) 
  reac_source_local(45,404) = + reac_rate_local(404) 
  reac_source_local(21,405) = + reac_rate_local(405) 
  reac_source_local(34,405) = - reac_rate_local(405) 
  reac_source_local(40,405) = - reac_rate_local(405) 
  reac_source_local(45,405) = + reac_rate_local(405) 
  reac_source_local(32,406) = + reac_rate_local(406) 
  reac_source_local(34,406) = - reac_rate_local(406) 
  reac_source_local(42,406) = - reac_rate_local(406) 
  reac_source_local(45,406) = + reac_rate_local(406) 
  reac_source_local(21,407) = + reac_rate_local(407) 
  reac_source_local(34,407) = - reac_rate_local(407) 
  reac_source_local(42,407) = - reac_rate_local(407) 
  reac_source_local(47,407) = + reac_rate_local(407) 
  reac_source_local(01,408) = + reac_rate_local(408) 
  reac_source_local(14,408) = + reac_rate_local(408) 
  reac_source_local(19,408) = - reac_rate_local(408) 
  reac_source_local(21,408) = - reac_rate_local(408) 
  reac_source_local(34,408) = + reac_rate_local(408) 
  reac_source_local(01,409) = + reac_rate_local(409) 
  reac_source_local(19,409) = - reac_rate_local(409) 
  reac_source_local(21,409) = - reac_rate_local(409) 
  reac_source_local(47,409) = + reac_rate_local(409) 
  reac_source_local(01,410) = + reac_rate_local(410) 
  reac_source_local(14,410) = - reac_rate_local(410) 
  reac_source_local(18,410) = + reac_rate_local(410) 
  reac_source_local(19,410) = - reac_rate_local(410) 
  reac_source_local(01,411) = + reac_rate_local(411) 
  reac_source_local(14,411) = + reac_rate_local(411) 
  reac_source_local(19,411) = - reac_rate_local(411) 
  reac_source_local(40,411) = - reac_rate_local(411) 
  reac_source_local(45,411) = + reac_rate_local(411) 
  reac_source_local(01,412) = + reac_rate_local(412) 
  reac_source_local(19,412) = - reac_rate_local(412) 
  reac_source_local(40,412) = - reac_rate_local(412) 
  reac_source_local(46,412) = + reac_rate_local(412) 
  reac_source_local(40,413) = - reac_rate_local(413) 
  reac_source_local(42,413) = + reac_rate_local(413) 
  reac_source_local(45,413) = + reac_rate_local(413) 
  reac_source_local(47,413) = - reac_rate_local(413) 
  reac_source_local(40,414) = - reac_rate_local(414) 
  reac_source_local(41,414) = + reac_rate_local(414) 
  reac_source_local(45,414) = + reac_rate_local(414) 
  reac_source_local(46,414) = - reac_rate_local(414) 
  reac_source_local(01,415) = + reac_rate_local(415) 
  reac_source_local(18,415) = + reac_rate_local(415) 
  reac_source_local(20,415) = - reac_rate_local(415) 
  reac_source_local(01,416) = + reac_rate_local(416) * 2.d0
  reac_source_local(20,416) = - reac_rate_local(416) 
  reac_source_local(21,416) = - reac_rate_local(416) 
  reac_source_local(34,416) = + reac_rate_local(416) 
  reac_source_local(01,417) = + reac_rate_local(417) * 2.d0
  reac_source_local(20,417) = - reac_rate_local(417) 
  reac_source_local(29,417) = - reac_rate_local(417) 
  reac_source_local(33,417) = + reac_rate_local(417) 
  reac_source_local(01,418) = + reac_rate_local(418) * 2.d0
  reac_source_local(14,418) = - reac_rate_local(418) 
  reac_source_local(17,418) = + reac_rate_local(418) 
  reac_source_local(20,418) = - reac_rate_local(418) 
  reac_source_local(01,419) = + reac_rate_local(419) * 2.d0
  reac_source_local(20,419) = - reac_rate_local(419) 
  reac_source_local(40,419) = - reac_rate_local(419) 
  reac_source_local(45,419) = + reac_rate_local(419) 
  reac_source_local(01,420) = - reac_rate_local(420) 
  reac_source_local(21,420) = + reac_rate_local(420) 
  reac_source_local(35,420) = - reac_rate_local(420) 
  reac_source_local(52,420) = + reac_rate_local(420) 
  reac_source_local(21,421) = + reac_rate_local(421) 
  reac_source_local(34,421) = + reac_rate_local(421) 
  reac_source_local(35,421) = - reac_rate_local(421) 
  reac_source_local(21,422) = + reac_rate_local(422) * 2.d0
  reac_source_local(26,422) = - reac_rate_local(422) 
  reac_source_local(34,422) = + reac_rate_local(422) 
  reac_source_local(35,422) = - reac_rate_local(422) 
  reac_source_local(21,423) = + reac_rate_local(423) * 2.d0
  reac_source_local(27,423) = - reac_rate_local(423) 
  reac_source_local(34,423) = + reac_rate_local(423) 
  reac_source_local(35,423) = - reac_rate_local(423) 
  reac_source_local(29,424) = - reac_rate_local(424) 
  reac_source_local(32,424) = + reac_rate_local(424) 
  reac_source_local(34,424) = + reac_rate_local(424) 
  reac_source_local(35,424) = - reac_rate_local(424) 
  reac_source_local(21,425) = + reac_rate_local(425) * 2.d0
  reac_source_local(35,425) = - reac_rate_local(425) 
  reac_source_local(40,425) = - reac_rate_local(425) 
  reac_source_local(45,425) = + reac_rate_local(425) 
  reac_source_local(01,426) = + reac_rate_local(426) 
  reac_source_local(34,426) = + reac_rate_local(426) 
  reac_source_local(52,426) = - reac_rate_local(426) 
  reac_source_local(01,427) = + reac_rate_local(427) 
  reac_source_local(21,427) = - reac_rate_local(427) 
  reac_source_local(35,427) = + reac_rate_local(427) 
  reac_source_local(52,427) = - reac_rate_local(427) 
  reac_source_local(01,428) = - reac_rate_local(428) 
  reac_source_local(17,428) = - reac_rate_local(428) 
  reac_source_local(19,428) = + reac_rate_local(428) 
  reac_source_local(17,429) = - reac_rate_local(429) 
  reac_source_local(29,429) = - reac_rate_local(429) 
  reac_source_local(45,429) = + reac_rate_local(429) 
  reac_source_local(14,430) = - reac_rate_local(430) 
  reac_source_local(17,430) = - reac_rate_local(430) 
  reac_source_local(18,430) = + reac_rate_local(430) 
  reac_source_local(01,431) = - reac_rate_local(431) 
  reac_source_local(14,431) = + reac_rate_local(431) 
  reac_source_local(33,431) = - reac_rate_local(431) 
  reac_source_local(45,431) = + reac_rate_local(431) 
  reac_source_local(29,432) = - reac_rate_local(432) 
  reac_source_local(33,432) = - reac_rate_local(432) 
  reac_source_local(34,432) = + reac_rate_local(432) 
  reac_source_local(14,433) = - reac_rate_local(433) 
  reac_source_local(33,433) = - reac_rate_local(433) 
  reac_source_local(45,433) = + reac_rate_local(433) 
  reac_source_local(01,434) = - reac_rate_local(434) 
  reac_source_local(18,434) = - reac_rate_local(434) 
  reac_source_local(20,434) = + reac_rate_local(434) 
  reac_source_local(14,435) = - reac_rate_local(435) 
  reac_source_local(18,435) = - reac_rate_local(435) 
  reac_source_local(19,435) = + reac_rate_local(435) 
  reac_source_local(21,436) = - reac_rate_local(436) 
  reac_source_local(34,436) = - reac_rate_local(436) 
  reac_source_local(35,436) = + reac_rate_local(436) 
  reac_source_local(01,437) = - reac_rate_local(437) 
  reac_source_local(34,437) = - reac_rate_local(437) 
  reac_source_local(52,437) = + reac_rate_local(437) 
  reac_source_local(26,438) = - reac_rate_local(438) 
  reac_source_local(29,438) = + reac_rate_local(438) 
  reac_source_local(36,438) = - reac_rate_local(438) 
  reac_source_local(37,438) = + reac_rate_local(438) 
  reac_source_local(29,439) = + reac_rate_local(439) 
  reac_source_local(32,439) = - reac_rate_local(439) 
  reac_source_local(36,439) = - reac_rate_local(439) 
  reac_source_local(38,439) = + reac_rate_local(439) 
  reac_source_local(29,440) = + reac_rate_local(440) 
  reac_source_local(36,440) = - reac_rate_local(440) 
  reac_source_local(42,440) = - reac_rate_local(440) 
  reac_source_local(50,440) = + reac_rate_local(440) 
  reac_source_local(36,441) = - reac_rate_local(441) 
  reac_source_local(40,441) = + reac_rate_local(441) 
  reac_source_local(41,441) = - reac_rate_local(441) 
  reac_source_local(48,441) = + reac_rate_local(441) 
  reac_source_local(29,442) = + reac_rate_local(442) 
  reac_source_local(36,442) = - reac_rate_local(442) 
  reac_source_local(41,442) = - reac_rate_local(442) 
  reac_source_local(49,442) = + reac_rate_local(442) 
  reac_source_local(21,443) = + reac_rate_local(443) 
  reac_source_local(29,443) = - reac_rate_local(443) 
  reac_source_local(36,443) = + reac_rate_local(443) 
  reac_source_local(37,443) = - reac_rate_local(443) 
  reac_source_local(21,444) = + reac_rate_local(444) 
  reac_source_local(32,444) = - reac_rate_local(444) 
  reac_source_local(37,444) = - reac_rate_local(444) 
  reac_source_local(38,444) = + reac_rate_local(444) 
  reac_source_local(21,445) = + reac_rate_local(445) 
  reac_source_local(37,445) = - reac_rate_local(445) 
  reac_source_local(42,445) = - reac_rate_local(445) 
  reac_source_local(50,445) = + reac_rate_local(445) 
  reac_source_local(21,446) = + reac_rate_local(446) 
  reac_source_local(37,446) = - reac_rate_local(446) 
  reac_source_local(43,446) = - reac_rate_local(446) 
  reac_source_local(51,446) = + reac_rate_local(446) 
  reac_source_local(21,447) = + reac_rate_local(447) 
  reac_source_local(29,447) = - reac_rate_local(447) 
  reac_source_local(37,447) = + reac_rate_local(447) 
  reac_source_local(38,447) = - reac_rate_local(447) 
  reac_source_local(29,448) = + reac_rate_local(448) 
  reac_source_local(38,448) = - reac_rate_local(448) 
  reac_source_local(40,448) = - reac_rate_local(448) 
  reac_source_local(51,448) = + reac_rate_local(448) 
  reac_source_local(21,449) = + reac_rate_local(449) 
  reac_source_local(38,449) = - reac_rate_local(449) 
  reac_source_local(40,449) = - reac_rate_local(449) 
  reac_source_local(50,449) = + reac_rate_local(449) 
  reac_source_local(32,450) = + reac_rate_local(450) 
  reac_source_local(38,450) = - reac_rate_local(450) 
  reac_source_local(42,450) = - reac_rate_local(450) 
  reac_source_local(50,450) = + reac_rate_local(450) 
  reac_source_local(21,451) = + reac_rate_local(451) 
  reac_source_local(38,451) = - reac_rate_local(451) 
  reac_source_local(42,451) = - reac_rate_local(451) 
  reac_source_local(51,451) = + reac_rate_local(451) 
  reac_source_local(32,452) = + reac_rate_local(452) 
  reac_source_local(38,452) = - reac_rate_local(452) 
  reac_source_local(43,452) = - reac_rate_local(452) 
  reac_source_local(51,452) = + reac_rate_local(452) 
  reac_source_local(21,453) = - reac_rate_local(453) 
  reac_source_local(37,453) = + reac_rate_local(453) 
  reac_source_local(40,453) = + reac_rate_local(453) 
  reac_source_local(48,453) = - reac_rate_local(453) 
  reac_source_local(40,454) = + reac_rate_local(454) 
  reac_source_local(42,454) = - reac_rate_local(454) 
  reac_source_local(48,454) = - reac_rate_local(454) 
  reac_source_local(50,454) = + reac_rate_local(454) 
  reac_source_local(01,455) = + reac_rate_local(455) 
  reac_source_local(41,455) = - reac_rate_local(455) 
  reac_source_local(48,455) = - reac_rate_local(455) 
  reac_source_local(50,455) = + reac_rate_local(455) 
  reac_source_local(21,456) = + reac_rate_local(456) 
  reac_source_local(32,456) = - reac_rate_local(456) 
  reac_source_local(50,456) = - reac_rate_local(456) 
  reac_source_local(51,456) = + reac_rate_local(456) 
  reac_source_local(40,457) = + reac_rate_local(457) 
  reac_source_local(42,457) = - reac_rate_local(457) 
  reac_source_local(50,457) = - reac_rate_local(457) 
  reac_source_local(51,457) = + reac_rate_local(457) 
  reac_source_local(42,458) = + reac_rate_local(458) 
  reac_source_local(43,458) = - reac_rate_local(458) 
  reac_source_local(50,458) = - reac_rate_local(458) 
  reac_source_local(51,458) = + reac_rate_local(458) 
  reac_source_local(42,459) = + reac_rate_local(459) * 2.d0
  reac_source_local(44,459) = - reac_rate_local(459) 
  reac_source_local(50,459) = - reac_rate_local(459) 
  reac_source_local(51,459) = + reac_rate_local(459) 
  reac_source_local(40,460) = - reac_rate_local(460) 
  reac_source_local(42,460) = + reac_rate_local(460) 
  reac_source_local(50,460) = + reac_rate_local(460) 
  reac_source_local(51,460) = - reac_rate_local(460) 
  reac_source_local(21,461) = + reac_rate_local(461) 
  reac_source_local(37,461) = + reac_rate_local(461) 
  reac_source_local(39,461) = - reac_rate_local(461) 
  reac_source_local(21,462) = + reac_rate_local(462) 
  reac_source_local(37,462) = + reac_rate_local(462) 
  reac_source_local(39,462) = - reac_rate_local(462) 
  reac_source_local(21,463) = + reac_rate_local(463) 
  reac_source_local(29,463) = - reac_rate_local(463) 
  reac_source_local(38,463) = + reac_rate_local(463) 
  reac_source_local(39,463) = - reac_rate_local(463) 
  reac_source_local(21,464) = + reac_rate_local(464) * 2.d0
  reac_source_local(29,464) = - reac_rate_local(464) 
  reac_source_local(36,464) = + reac_rate_local(464) 
  reac_source_local(39,464) = - reac_rate_local(464) 
  reac_source_local(21,465) = + reac_rate_local(465) * 2.d0
  reac_source_local(26,465) = - reac_rate_local(465) 
  reac_source_local(37,465) = + reac_rate_local(465) 
  reac_source_local(39,465) = - reac_rate_local(465) 
  reac_source_local(21,466) = + reac_rate_local(466) * 2.d0
  reac_source_local(27,466) = - reac_rate_local(466) 
  reac_source_local(37,466) = + reac_rate_local(466) 
  reac_source_local(39,466) = - reac_rate_local(466) 
  reac_source_local(21,467) = + reac_rate_local(467) 
  reac_source_local(39,467) = - reac_rate_local(467) 
  reac_source_local(40,467) = - reac_rate_local(467) 
  reac_source_local(51,467) = + reac_rate_local(467) 
  reac_source_local(21,468) = - reac_rate_local(468) 
  reac_source_local(36,468) = - reac_rate_local(468) 
  reac_source_local(38,468) = + reac_rate_local(468) 
  reac_source_local(36,469) = - reac_rate_local(469) 
  reac_source_local(40,469) = - reac_rate_local(469) 
  reac_source_local(50,469) = + reac_rate_local(469) 
  reac_source_local(21,470) = - reac_rate_local(470) 
  reac_source_local(37,470) = - reac_rate_local(470) 
  reac_source_local(39,470) = + reac_rate_local(470) 
  reac_source_local(14,471) = + reac_rate_local(471) 
  reac_source_local(17,471) = - reac_rate_local(471) 
  reac_source_local(29,471) = + reac_rate_local(471) 
  reac_source_local(36,471) = - reac_rate_local(471) 
  reac_source_local(01,472) = + reac_rate_local(472) 
  reac_source_local(18,472) = - reac_rate_local(472) 
  reac_source_local(29,472) = + reac_rate_local(472) 
  reac_source_local(36,472) = - reac_rate_local(472) 
  reac_source_local(29,473) = + reac_rate_local(473) * 2.d0
  reac_source_local(33,473) = - reac_rate_local(473) 
  reac_source_local(36,473) = - reac_rate_local(473) 
  reac_source_local(21,474) = + reac_rate_local(474) 
  reac_source_local(29,474) = + reac_rate_local(474) 
  reac_source_local(34,474) = - reac_rate_local(474) 
  reac_source_local(36,474) = - reac_rate_local(474) 
  reac_source_local(29,475) = + reac_rate_local(475) 
  reac_source_local(36,475) = - reac_rate_local(475) 
  reac_source_local(40,475) = + reac_rate_local(475) 
  reac_source_local(45,475) = - reac_rate_local(475) 
  reac_source_local(29,476) = + reac_rate_local(476) 
  reac_source_local(36,476) = - reac_rate_local(476) 
  reac_source_local(41,476) = + reac_rate_local(476) 
  reac_source_local(46,476) = - reac_rate_local(476) 
  reac_source_local(29,477) = + reac_rate_local(477) 
  reac_source_local(36,477) = - reac_rate_local(477) 
  reac_source_local(42,477) = + reac_rate_local(477) 
  reac_source_local(47,477) = - reac_rate_local(477) 
  reac_source_local(14,478) = + reac_rate_local(478) 
  reac_source_local(17,478) = - reac_rate_local(478) 
  reac_source_local(21,478) = + reac_rate_local(478) 
  reac_source_local(37,478) = - reac_rate_local(478) 
  reac_source_local(01,479) = + reac_rate_local(479) 
  reac_source_local(18,479) = - reac_rate_local(479) 
  reac_source_local(21,479) = + reac_rate_local(479) 
  reac_source_local(37,479) = - reac_rate_local(479) 
  reac_source_local(21,480) = + reac_rate_local(480) 
  reac_source_local(29,480) = + reac_rate_local(480) 
  reac_source_local(33,480) = - reac_rate_local(480) 
  reac_source_local(37,480) = - reac_rate_local(480) 
  reac_source_local(21,481) = + reac_rate_local(481) * 2.d0
  reac_source_local(34,481) = - reac_rate_local(481) 
  reac_source_local(37,481) = - reac_rate_local(481) 
  reac_source_local(21,482) = + reac_rate_local(482) 
  reac_source_local(37,482) = - reac_rate_local(482) 
  reac_source_local(40,482) = + reac_rate_local(482) 
  reac_source_local(45,482) = - reac_rate_local(482) 
  reac_source_local(21,483) = + reac_rate_local(483) 
  reac_source_local(37,483) = - reac_rate_local(483) 
  reac_source_local(41,483) = + reac_rate_local(483) 
  reac_source_local(46,483) = - reac_rate_local(483) 
  reac_source_local(21,484) = + reac_rate_local(484) 
  reac_source_local(37,484) = - reac_rate_local(484) 
  reac_source_local(42,484) = + reac_rate_local(484) 
  reac_source_local(47,484) = - reac_rate_local(484) 
  reac_source_local(14,485) = + reac_rate_local(485) 
  reac_source_local(17,485) = - reac_rate_local(485) 
  reac_source_local(32,485) = + reac_rate_local(485) 
  reac_source_local(38,485) = - reac_rate_local(485) 
  reac_source_local(01,486) = + reac_rate_local(486) 
  reac_source_local(18,486) = - reac_rate_local(486) 
  reac_source_local(32,486) = + reac_rate_local(486) 
  reac_source_local(38,486) = - reac_rate_local(486) 
  reac_source_local(29,487) = + reac_rate_local(487) 
  reac_source_local(32,487) = + reac_rate_local(487) 
  reac_source_local(33,487) = - reac_rate_local(487) 
  reac_source_local(38,487) = - reac_rate_local(487) 
  reac_source_local(21,488) = + reac_rate_local(488) 
  reac_source_local(32,488) = + reac_rate_local(488) 
  reac_source_local(34,488) = - reac_rate_local(488) 
  reac_source_local(38,488) = - reac_rate_local(488) 
  reac_source_local(32,489) = + reac_rate_local(489) 
  reac_source_local(38,489) = - reac_rate_local(489) 
  reac_source_local(40,489) = + reac_rate_local(489) 
  reac_source_local(45,489) = - reac_rate_local(489) 
  reac_source_local(32,490) = + reac_rate_local(490) 
  reac_source_local(38,490) = - reac_rate_local(490) 
  reac_source_local(41,490) = + reac_rate_local(490) 
  reac_source_local(46,490) = - reac_rate_local(490) 
  reac_source_local(32,491) = + reac_rate_local(491) 
  reac_source_local(38,491) = - reac_rate_local(491) 
  reac_source_local(42,491) = + reac_rate_local(491) 
  reac_source_local(47,491) = - reac_rate_local(491) 
  reac_source_local(14,492) = + reac_rate_local(492) 
  reac_source_local(17,492) = - reac_rate_local(492) 
  reac_source_local(40,492) = + reac_rate_local(492) 
  reac_source_local(48,492) = - reac_rate_local(492) 
  reac_source_local(01,493) = + reac_rate_local(493) 
  reac_source_local(18,493) = - reac_rate_local(493) 
  reac_source_local(40,493) = + reac_rate_local(493) 
  reac_source_local(48,493) = - reac_rate_local(493) 
  reac_source_local(29,494) = + reac_rate_local(494) 
  reac_source_local(33,494) = - reac_rate_local(494) 
  reac_source_local(40,494) = + reac_rate_local(494) 
  reac_source_local(48,494) = - reac_rate_local(494) 
  reac_source_local(21,495) = + reac_rate_local(495) 
  reac_source_local(34,495) = - reac_rate_local(495) 
  reac_source_local(40,495) = + reac_rate_local(495) 
  reac_source_local(48,495) = - reac_rate_local(495) 
  reac_source_local(40,496) = + reac_rate_local(496) * 2.d0
  reac_source_local(45,496) = - reac_rate_local(496) 
  reac_source_local(48,496) = - reac_rate_local(496) 
  reac_source_local(40,497) = + reac_rate_local(497) 
  reac_source_local(41,497) = + reac_rate_local(497) 
  reac_source_local(46,497) = - reac_rate_local(497) 
  reac_source_local(48,497) = - reac_rate_local(497) 
  reac_source_local(40,498) = + reac_rate_local(498) 
  reac_source_local(42,498) = + reac_rate_local(498) 
  reac_source_local(47,498) = - reac_rate_local(498) 
  reac_source_local(48,498) = - reac_rate_local(498) 
  reac_source_local(14,499) = + reac_rate_local(499) 
  reac_source_local(17,499) = - reac_rate_local(499) 
  reac_source_local(41,499) = + reac_rate_local(499) 
  reac_source_local(49,499) = - reac_rate_local(499) 
  reac_source_local(01,500) = + reac_rate_local(500) 
  reac_source_local(18,500) = - reac_rate_local(500) 
  reac_source_local(41,500) = + reac_rate_local(500) 
  reac_source_local(49,500) = - reac_rate_local(500) 
  reac_source_local(29,501) = + reac_rate_local(501) 
  reac_source_local(33,501) = - reac_rate_local(501) 
  reac_source_local(41,501) = + reac_rate_local(501) 
  reac_source_local(49,501) = - reac_rate_local(501) 
  reac_source_local(21,502) = + reac_rate_local(502) 
  reac_source_local(34,502) = - reac_rate_local(502) 
  reac_source_local(41,502) = + reac_rate_local(502) 
  reac_source_local(49,502) = - reac_rate_local(502) 
  reac_source_local(40,503) = + reac_rate_local(503) 
  reac_source_local(41,503) = + reac_rate_local(503) 
  reac_source_local(45,503) = - reac_rate_local(503) 
  reac_source_local(49,503) = - reac_rate_local(503) 
  reac_source_local(41,504) = + reac_rate_local(504) * 2.d0
  reac_source_local(46,504) = - reac_rate_local(504) 
  reac_source_local(49,504) = - reac_rate_local(504) 
  reac_source_local(41,505) = + reac_rate_local(505) 
  reac_source_local(42,505) = + reac_rate_local(505) 
  reac_source_local(47,505) = - reac_rate_local(505) 
  reac_source_local(49,505) = - reac_rate_local(505) 
  reac_source_local(14,506) = + reac_rate_local(506) 
  reac_source_local(17,506) = - reac_rate_local(506) 
  reac_source_local(42,506) = + reac_rate_local(506) 
  reac_source_local(50,506) = - reac_rate_local(506) 
  reac_source_local(01,507) = + reac_rate_local(507) 
  reac_source_local(18,507) = - reac_rate_local(507) 
  reac_source_local(42,507) = + reac_rate_local(507) 
  reac_source_local(50,507) = - reac_rate_local(507) 
  reac_source_local(29,508) = + reac_rate_local(508) 
  reac_source_local(33,508) = - reac_rate_local(508) 
  reac_source_local(42,508) = + reac_rate_local(508) 
  reac_source_local(50,508) = - reac_rate_local(508) 
  reac_source_local(21,509) = + reac_rate_local(509) 
  reac_source_local(34,509) = - reac_rate_local(509) 
  reac_source_local(42,509) = + reac_rate_local(509) 
  reac_source_local(50,509) = - reac_rate_local(509) 
  reac_source_local(40,510) = + reac_rate_local(510) 
  reac_source_local(42,510) = + reac_rate_local(510) 
  reac_source_local(45,510) = - reac_rate_local(510) 
  reac_source_local(50,510) = - reac_rate_local(510) 
  reac_source_local(41,511) = + reac_rate_local(511) 
  reac_source_local(42,511) = + reac_rate_local(511) 
  reac_source_local(46,511) = - reac_rate_local(511) 
  reac_source_local(50,511) = - reac_rate_local(511) 
  reac_source_local(42,512) = + reac_rate_local(512) * 2.d0
  reac_source_local(47,512) = - reac_rate_local(512) 
  reac_source_local(50,512) = - reac_rate_local(512) 
  reac_source_local(14,513) = + reac_rate_local(513) 
  reac_source_local(17,513) = - reac_rate_local(513) 
  reac_source_local(43,513) = + reac_rate_local(513) 
  reac_source_local(51,513) = - reac_rate_local(513) 
  reac_source_local(01,514) = + reac_rate_local(514) 
  reac_source_local(18,514) = - reac_rate_local(514) 
  reac_source_local(43,514) = + reac_rate_local(514) 
  reac_source_local(51,514) = - reac_rate_local(514) 
  reac_source_local(29,515) = + reac_rate_local(515) 
  reac_source_local(33,515) = - reac_rate_local(515) 
  reac_source_local(43,515) = + reac_rate_local(515) 
  reac_source_local(51,515) = - reac_rate_local(515) 
  reac_source_local(21,516) = + reac_rate_local(516) 
  reac_source_local(34,516) = - reac_rate_local(516) 
  reac_source_local(43,516) = + reac_rate_local(516) 
  reac_source_local(51,516) = - reac_rate_local(516) 
  reac_source_local(40,517) = + reac_rate_local(517) 
  reac_source_local(43,517) = + reac_rate_local(517) 
  reac_source_local(45,517) = - reac_rate_local(517) 
  reac_source_local(51,517) = - reac_rate_local(517) 
  reac_source_local(41,518) = + reac_rate_local(518) 
  reac_source_local(43,518) = + reac_rate_local(518) 
  reac_source_local(46,518) = - reac_rate_local(518) 
  reac_source_local(51,518) = - reac_rate_local(518) 
  reac_source_local(42,519) = + reac_rate_local(519) 
  reac_source_local(43,519) = + reac_rate_local(519) 
  reac_source_local(47,519) = - reac_rate_local(519) 
  reac_source_local(51,519) = - reac_rate_local(519) 
  reac_source_local(14,520) = + reac_rate_local(520) * 2.d0
  reac_source_local(18,520) = - reac_rate_local(520) 
  reac_source_local(29,520) = + reac_rate_local(520) 
  reac_source_local(36,520) = - reac_rate_local(520) 
  reac_source_local(01,521) = + reac_rate_local(521) 
  reac_source_local(14,521) = + reac_rate_local(521) 
  reac_source_local(19,521) = - reac_rate_local(521) 
  reac_source_local(29,521) = + reac_rate_local(521) 
  reac_source_local(36,521) = - reac_rate_local(521) 
  reac_source_local(01,522) = + reac_rate_local(522) * 2.d0
  reac_source_local(20,522) = - reac_rate_local(522) 
  reac_source_local(29,522) = + reac_rate_local(522) 
  reac_source_local(36,522) = - reac_rate_local(522) 
  reac_source_local(29,523) = + reac_rate_local(523) * 3.d0
  reac_source_local(34,523) = - reac_rate_local(523) 
  reac_source_local(36,523) = - reac_rate_local(523) 
  reac_source_local(21,524) = + reac_rate_local(524) * 2.d0
  reac_source_local(29,524) = + reac_rate_local(524) 
  reac_source_local(35,524) = - reac_rate_local(524) 
  reac_source_local(36,524) = - reac_rate_local(524) 
  reac_source_local(14,525) = + reac_rate_local(525) 
  reac_source_local(29,525) = + reac_rate_local(525) * 2.d0
  reac_source_local(36,525) = - reac_rate_local(525) 
  reac_source_local(45,525) = - reac_rate_local(525) 
  reac_source_local(01,526) = + reac_rate_local(526) 
  reac_source_local(29,526) = + reac_rate_local(526) * 2.d0
  reac_source_local(36,526) = - reac_rate_local(526) 
  reac_source_local(46,526) = - reac_rate_local(526) 
  reac_source_local(14,527) = + reac_rate_local(527) 
  reac_source_local(21,527) = + reac_rate_local(527) 
  reac_source_local(29,527) = + reac_rate_local(527) 
  reac_source_local(36,527) = - reac_rate_local(527) 
  reac_source_local(47,527) = - reac_rate_local(527) 
  reac_source_local(01,528) = + reac_rate_local(528) 
  reac_source_local(21,528) = + reac_rate_local(528) 
  reac_source_local(29,528) = + reac_rate_local(528) 
  reac_source_local(36,528) = - reac_rate_local(528) 
  reac_source_local(52,528) = - reac_rate_local(528) 
  reac_source_local(14,529) = + reac_rate_local(529) * 2.d0
  reac_source_local(18,529) = - reac_rate_local(529) 
  reac_source_local(21,529) = + reac_rate_local(529) 
  reac_source_local(37,529) = - reac_rate_local(529) 
  reac_source_local(01,530) = + reac_rate_local(530) 
  reac_source_local(14,530) = + reac_rate_local(530) 
  reac_source_local(19,530) = - reac_rate_local(530) 
  reac_source_local(21,530) = + reac_rate_local(530) 
  reac_source_local(37,530) = - reac_rate_local(530) 
  reac_source_local(01,531) = + reac_rate_local(531) * 2.d0
  reac_source_local(20,531) = - reac_rate_local(531) 
  reac_source_local(21,531) = + reac_rate_local(531) 
  reac_source_local(37,531) = - reac_rate_local(531) 
  reac_source_local(21,532) = + reac_rate_local(532) 
  reac_source_local(29,532) = + reac_rate_local(532) * 2.d0
  reac_source_local(34,532) = - reac_rate_local(532) 
  reac_source_local(37,532) = - reac_rate_local(532) 
  reac_source_local(21,533) = + reac_rate_local(533) * 3.d0
  reac_source_local(35,533) = - reac_rate_local(533) 
  reac_source_local(37,533) = - reac_rate_local(533) 
  reac_source_local(14,534) = + reac_rate_local(534) 
  reac_source_local(21,534) = + reac_rate_local(534) 
  reac_source_local(29,534) = + reac_rate_local(534) 
  reac_source_local(37,534) = - reac_rate_local(534) 
  reac_source_local(45,534) = - reac_rate_local(534) 
  reac_source_local(01,535) = + reac_rate_local(535) 
  reac_source_local(21,535) = + reac_rate_local(535) 
  reac_source_local(29,535) = + reac_rate_local(535) 
  reac_source_local(37,535) = - reac_rate_local(535) 
  reac_source_local(46,535) = - reac_rate_local(535) 
  reac_source_local(14,536) = + reac_rate_local(536) 
  reac_source_local(21,536) = + reac_rate_local(536) * 2.d0
  reac_source_local(37,536) = - reac_rate_local(536) 
  reac_source_local(47,536) = - reac_rate_local(536) 
  reac_source_local(01,537) = + reac_rate_local(537) 
  reac_source_local(21,537) = + reac_rate_local(537) * 2.d0
  reac_source_local(37,537) = - reac_rate_local(537) 
  reac_source_local(52,537) = - reac_rate_local(537) 
  reac_source_local(14,538) = + reac_rate_local(538) * 2.d0
  reac_source_local(18,538) = - reac_rate_local(538) 
  reac_source_local(32,538) = + reac_rate_local(538) 
  reac_source_local(38,538) = - reac_rate_local(538) 
  reac_source_local(01,539) = + reac_rate_local(539) 
  reac_source_local(14,539) = + reac_rate_local(539) 
  reac_source_local(19,539) = - reac_rate_local(539) 
  reac_source_local(32,539) = + reac_rate_local(539) 
  reac_source_local(38,539) = - reac_rate_local(539) 
  reac_source_local(01,540) = + reac_rate_local(540) * 2.d0
  reac_source_local(20,540) = - reac_rate_local(540) 
  reac_source_local(32,540) = + reac_rate_local(540) 
  reac_source_local(38,540) = - reac_rate_local(540) 
  reac_source_local(29,541) = + reac_rate_local(541) * 2.d0
  reac_source_local(32,541) = + reac_rate_local(541) 
  reac_source_local(34,541) = - reac_rate_local(541) 
  reac_source_local(38,541) = - reac_rate_local(541) 
  reac_source_local(21,542) = + reac_rate_local(542) * 2.d0
  reac_source_local(32,542) = + reac_rate_local(542) 
  reac_source_local(35,542) = - reac_rate_local(542) 
  reac_source_local(38,542) = - reac_rate_local(542) 
  reac_source_local(14,543) = + reac_rate_local(543) 
  reac_source_local(29,543) = + reac_rate_local(543) 
  reac_source_local(32,543) = + reac_rate_local(543) 
  reac_source_local(38,543) = - reac_rate_local(543) 
  reac_source_local(45,543) = - reac_rate_local(543) 
  reac_source_local(01,544) = + reac_rate_local(544) 
  reac_source_local(29,544) = + reac_rate_local(544) 
  reac_source_local(32,544) = + reac_rate_local(544) 
  reac_source_local(38,544) = - reac_rate_local(544) 
  reac_source_local(46,544) = - reac_rate_local(544) 
  reac_source_local(14,545) = + reac_rate_local(545) 
  reac_source_local(21,545) = + reac_rate_local(545) 
  reac_source_local(32,545) = + reac_rate_local(545) 
  reac_source_local(38,545) = - reac_rate_local(545) 
  reac_source_local(47,545) = - reac_rate_local(545) 
  reac_source_local(01,546) = + reac_rate_local(546) 
  reac_source_local(21,546) = + reac_rate_local(546) 
  reac_source_local(32,546) = + reac_rate_local(546) 
  reac_source_local(38,546) = - reac_rate_local(546) 
  reac_source_local(52,546) = - reac_rate_local(546) 
  reac_source_local(14,547) = + reac_rate_local(547) * 2.d0
  reac_source_local(18,547) = - reac_rate_local(547) 
  reac_source_local(40,547) = + reac_rate_local(547) 
  reac_source_local(48,547) = - reac_rate_local(547) 
  reac_source_local(01,548) = + reac_rate_local(548) 
  reac_source_local(14,548) = + reac_rate_local(548) 
  reac_source_local(19,548) = - reac_rate_local(548) 
  reac_source_local(40,548) = + reac_rate_local(548) 
  reac_source_local(48,548) = - reac_rate_local(548) 
  reac_source_local(01,549) = + reac_rate_local(549) * 2.d0
  reac_source_local(20,549) = - reac_rate_local(549) 
  reac_source_local(40,549) = + reac_rate_local(549) 
  reac_source_local(48,549) = - reac_rate_local(549) 
  reac_source_local(29,550) = + reac_rate_local(550) * 2.d0
  reac_source_local(34,550) = - reac_rate_local(550) 
  reac_source_local(40,550) = + reac_rate_local(550) 
  reac_source_local(48,550) = - reac_rate_local(550) 
  reac_source_local(21,551) = + reac_rate_local(551) * 2.d0
  reac_source_local(35,551) = - reac_rate_local(551) 
  reac_source_local(40,551) = + reac_rate_local(551) 
  reac_source_local(48,551) = - reac_rate_local(551) 
  reac_source_local(14,552) = + reac_rate_local(552) 
  reac_source_local(29,552) = + reac_rate_local(552) 
  reac_source_local(40,552) = + reac_rate_local(552) 
  reac_source_local(45,552) = - reac_rate_local(552) 
  reac_source_local(48,552) = - reac_rate_local(552) 
  reac_source_local(01,553) = + reac_rate_local(553) 
  reac_source_local(29,553) = + reac_rate_local(553) 
  reac_source_local(40,553) = + reac_rate_local(553) 
  reac_source_local(46,553) = - reac_rate_local(553) 
  reac_source_local(48,553) = - reac_rate_local(553) 
  reac_source_local(14,554) = + reac_rate_local(554) 
  reac_source_local(21,554) = + reac_rate_local(554) 
  reac_source_local(40,554) = + reac_rate_local(554) 
  reac_source_local(47,554) = - reac_rate_local(554) 
  reac_source_local(48,554) = - reac_rate_local(554) 
  reac_source_local(01,555) = + reac_rate_local(555) 
  reac_source_local(21,555) = + reac_rate_local(555) 
  reac_source_local(40,555) = + reac_rate_local(555) 
  reac_source_local(48,555) = - reac_rate_local(555) 
  reac_source_local(52,555) = - reac_rate_local(555) 
  reac_source_local(14,556) = + reac_rate_local(556) * 2.d0
  reac_source_local(18,556) = - reac_rate_local(556) 
  reac_source_local(41,556) = + reac_rate_local(556) 
  reac_source_local(49,556) = - reac_rate_local(556) 
  reac_source_local(01,557) = + reac_rate_local(557) 
  reac_source_local(14,557) = + reac_rate_local(557) 
  reac_source_local(19,557) = - reac_rate_local(557) 
  reac_source_local(41,557) = + reac_rate_local(557) 
  reac_source_local(49,557) = - reac_rate_local(557) 
  reac_source_local(01,558) = + reac_rate_local(558) * 2.d0
  reac_source_local(20,558) = - reac_rate_local(558) 
  reac_source_local(41,558) = + reac_rate_local(558) 
  reac_source_local(49,558) = - reac_rate_local(558) 
  reac_source_local(29,559) = + reac_rate_local(559) * 2.d0
  reac_source_local(34,559) = - reac_rate_local(559) 
  reac_source_local(41,559) = + reac_rate_local(559) 
  reac_source_local(49,559) = - reac_rate_local(559) 
  reac_source_local(21,560) = + reac_rate_local(560) * 2.d0
  reac_source_local(35,560) = - reac_rate_local(560) 
  reac_source_local(41,560) = + reac_rate_local(560) 
  reac_source_local(49,560) = - reac_rate_local(560) 
  reac_source_local(14,561) = + reac_rate_local(561) 
  reac_source_local(29,561) = + reac_rate_local(561) 
  reac_source_local(41,561) = + reac_rate_local(561) 
  reac_source_local(45,561) = - reac_rate_local(561) 
  reac_source_local(49,561) = - reac_rate_local(561) 
  reac_source_local(01,562) = + reac_rate_local(562) 
  reac_source_local(29,562) = + reac_rate_local(562) 
  reac_source_local(41,562) = + reac_rate_local(562) 
  reac_source_local(46,562) = - reac_rate_local(562) 
  reac_source_local(49,562) = - reac_rate_local(562) 
  reac_source_local(14,563) = + reac_rate_local(563) 
  reac_source_local(21,563) = + reac_rate_local(563) 
  reac_source_local(41,563) = + reac_rate_local(563) 
  reac_source_local(47,563) = - reac_rate_local(563) 
  reac_source_local(49,563) = - reac_rate_local(563) 
  reac_source_local(01,564) = + reac_rate_local(564) 
  reac_source_local(21,564) = + reac_rate_local(564) 
  reac_source_local(41,564) = + reac_rate_local(564) 
  reac_source_local(49,564) = - reac_rate_local(564) 
  reac_source_local(52,564) = - reac_rate_local(564) 
  reac_source_local(14,565) = + reac_rate_local(565) * 2.d0
  reac_source_local(18,565) = - reac_rate_local(565) 
  reac_source_local(42,565) = + reac_rate_local(565) 
  reac_source_local(50,565) = - reac_rate_local(565) 
  reac_source_local(01,566) = + reac_rate_local(566) 
  reac_source_local(14,566) = + reac_rate_local(566) 
  reac_source_local(19,566) = - reac_rate_local(566) 
  reac_source_local(42,566) = + reac_rate_local(566) 
  reac_source_local(50,566) = - reac_rate_local(566) 
  reac_source_local(01,567) = + reac_rate_local(567) * 2.d0
  reac_source_local(20,567) = - reac_rate_local(567) 
  reac_source_local(42,567) = + reac_rate_local(567) 
  reac_source_local(50,567) = - reac_rate_local(567) 
  reac_source_local(29,568) = + reac_rate_local(568) * 2.d0
  reac_source_local(34,568) = - reac_rate_local(568) 
  reac_source_local(42,568) = + reac_rate_local(568) 
  reac_source_local(50,568) = - reac_rate_local(568) 
  reac_source_local(21,569) = + reac_rate_local(569) * 2.d0
  reac_source_local(35,569) = - reac_rate_local(569) 
  reac_source_local(42,569) = + reac_rate_local(569) 
  reac_source_local(50,569) = - reac_rate_local(569) 
  reac_source_local(14,570) = + reac_rate_local(570) 
  reac_source_local(29,570) = + reac_rate_local(570) 
  reac_source_local(42,570) = + reac_rate_local(570) 
  reac_source_local(45,570) = - reac_rate_local(570) 
  reac_source_local(50,570) = - reac_rate_local(570) 
  reac_source_local(01,571) = + reac_rate_local(571) 
  reac_source_local(29,571) = + reac_rate_local(571) 
  reac_source_local(42,571) = + reac_rate_local(571) 
  reac_source_local(46,571) = - reac_rate_local(571) 
  reac_source_local(50,571) = - reac_rate_local(571) 
  reac_source_local(14,572) = + reac_rate_local(572) 
  reac_source_local(21,572) = + reac_rate_local(572) 
  reac_source_local(42,572) = + reac_rate_local(572) 
  reac_source_local(47,572) = - reac_rate_local(572) 
  reac_source_local(50,572) = - reac_rate_local(572) 
  reac_source_local(01,573) = + reac_rate_local(573) 
  reac_source_local(21,573) = + reac_rate_local(573) 
  reac_source_local(42,573) = + reac_rate_local(573) 
  reac_source_local(50,573) = - reac_rate_local(573) 
  reac_source_local(52,573) = - reac_rate_local(573) 
  reac_source_local(14,574) = + reac_rate_local(574) * 2.d0
  reac_source_local(18,574) = - reac_rate_local(574) 
  reac_source_local(43,574) = + reac_rate_local(574) 
  reac_source_local(51,574) = - reac_rate_local(574) 
  reac_source_local(01,575) = + reac_rate_local(575) 
  reac_source_local(14,575) = + reac_rate_local(575) 
  reac_source_local(19,575) = - reac_rate_local(575) 
  reac_source_local(43,575) = + reac_rate_local(575) 
  reac_source_local(51,575) = - reac_rate_local(575) 
  reac_source_local(01,576) = + reac_rate_local(576) * 2.d0
  reac_source_local(20,576) = - reac_rate_local(576) 
  reac_source_local(43,576) = + reac_rate_local(576) 
  reac_source_local(51,576) = - reac_rate_local(576) 
  reac_source_local(29,577) = + reac_rate_local(577) * 2.d0
  reac_source_local(34,577) = - reac_rate_local(577) 
  reac_source_local(43,577) = + reac_rate_local(577) 
  reac_source_local(51,577) = - reac_rate_local(577) 
  reac_source_local(21,578) = + reac_rate_local(578) * 2.d0
  reac_source_local(35,578) = - reac_rate_local(578) 
  reac_source_local(43,578) = + reac_rate_local(578) 
  reac_source_local(51,578) = - reac_rate_local(578) 
  reac_source_local(14,579) = + reac_rate_local(579) 
  reac_source_local(29,579) = + reac_rate_local(579) 
  reac_source_local(43,579) = + reac_rate_local(579) 
  reac_source_local(45,579) = - reac_rate_local(579) 
  reac_source_local(51,579) = - reac_rate_local(579) 
  reac_source_local(01,580) = + reac_rate_local(580) 
  reac_source_local(29,580) = + reac_rate_local(580) 
  reac_source_local(43,580) = + reac_rate_local(580) 
  reac_source_local(46,580) = - reac_rate_local(580) 
  reac_source_local(51,580) = - reac_rate_local(580) 
  reac_source_local(14,581) = + reac_rate_local(581) 
  reac_source_local(21,581) = + reac_rate_local(581) 
  reac_source_local(43,581) = + reac_rate_local(581) 
  reac_source_local(47,581) = - reac_rate_local(581) 
  reac_source_local(51,581) = - reac_rate_local(581) 
  reac_source_local(01,582) = + reac_rate_local(582) 
  reac_source_local(21,582) = + reac_rate_local(582) 
  reac_source_local(43,582) = + reac_rate_local(582) 
  reac_source_local(51,582) = - reac_rate_local(582) 
  reac_source_local(52,582) = - reac_rate_local(582) 
  reac_source_local(14,583) = + reac_rate_local(583) 
  reac_source_local(17,583) = - reac_rate_local(583) 
  reac_source_local(21,583) = + reac_rate_local(583) * 2.d0
  reac_source_local(39,583) = - reac_rate_local(583) 
  reac_source_local(01,584) = + reac_rate_local(584) 
  reac_source_local(18,584) = - reac_rate_local(584) 
  reac_source_local(21,584) = + reac_rate_local(584) * 2.d0
  reac_source_local(39,584) = - reac_rate_local(584) 
  reac_source_local(21,585) = + reac_rate_local(585) * 2.d0
  reac_source_local(29,585) = + reac_rate_local(585) 
  reac_source_local(33,585) = - reac_rate_local(585) 
  reac_source_local(39,585) = - reac_rate_local(585) 
  reac_source_local(21,586) = + reac_rate_local(586) * 3.d0
  reac_source_local(34,586) = - reac_rate_local(586) 
  reac_source_local(39,586) = - reac_rate_local(586) 
  reac_source_local(21,587) = + reac_rate_local(587) * 2.d0
  reac_source_local(39,587) = - reac_rate_local(587) 
  reac_source_local(40,587) = + reac_rate_local(587) 
  reac_source_local(45,587) = - reac_rate_local(587) 
  reac_source_local(21,588) = + reac_rate_local(588) * 2.d0
  reac_source_local(39,588) = - reac_rate_local(588) 
  reac_source_local(41,588) = + reac_rate_local(588) 
  reac_source_local(46,588) = - reac_rate_local(588) 
  reac_source_local(21,589) = + reac_rate_local(589) * 2.d0
  reac_source_local(39,589) = - reac_rate_local(589) 
  reac_source_local(42,589) = + reac_rate_local(589) 
  reac_source_local(47,589) = - reac_rate_local(589) 
  reac_source_local(01,590) = + reac_rate_local(590) 
  reac_source_local(14,590) = + reac_rate_local(590) 
  reac_source_local(19,590) = - reac_rate_local(590) 
  reac_source_local(21,590) = + reac_rate_local(590) * 2.d0
  reac_source_local(39,590) = - reac_rate_local(590) 
  reac_source_local(01,591) = + reac_rate_local(591) * 2.d0
  reac_source_local(20,591) = - reac_rate_local(591) 
  reac_source_local(21,591) = + reac_rate_local(591) * 2.d0
  reac_source_local(39,591) = - reac_rate_local(591) 
  reac_source_local(21,592) = + reac_rate_local(592) * 4.d0
  reac_source_local(35,592) = - reac_rate_local(592) 
  reac_source_local(39,592) = - reac_rate_local(592) 
  reac_source_local(01,593) = + reac_rate_local(593) 
  reac_source_local(21,593) = + reac_rate_local(593) * 3.d0
  reac_source_local(39,593) = - reac_rate_local(593) 
  reac_source_local(52,593) = - reac_rate_local(593) 
  reac_source_local(14,594) = + reac_rate_local(594) 
  reac_source_local(17,594) = - reac_rate_local(594) 
  reac_source_local(29,594) = + reac_rate_local(594) 
  reac_source_local(36,594) = - reac_rate_local(594) 
  reac_source_local(01,595) = + reac_rate_local(595) 
  reac_source_local(18,595) = - reac_rate_local(595) 
  reac_source_local(29,595) = + reac_rate_local(595) 
  reac_source_local(36,595) = - reac_rate_local(595) 
  reac_source_local(29,596) = + reac_rate_local(596) * 2.d0
  reac_source_local(33,596) = - reac_rate_local(596) 
  reac_source_local(36,596) = - reac_rate_local(596) 
  reac_source_local(21,597) = + reac_rate_local(597) 
  reac_source_local(29,597) = + reac_rate_local(597) 
  reac_source_local(34,597) = - reac_rate_local(597) 
  reac_source_local(36,597) = - reac_rate_local(597) 
  reac_source_local(29,598) = + reac_rate_local(598) 
  reac_source_local(36,598) = - reac_rate_local(598) 
  reac_source_local(40,598) = + reac_rate_local(598) 
  reac_source_local(45,598) = - reac_rate_local(598) 
  reac_source_local(14,599) = + reac_rate_local(599) 
  reac_source_local(17,599) = - reac_rate_local(599) 
  reac_source_local(21,599) = + reac_rate_local(599) 
  reac_source_local(37,599) = - reac_rate_local(599) 
  reac_source_local(01,600) = + reac_rate_local(600) 
  reac_source_local(18,600) = - reac_rate_local(600) 
  reac_source_local(21,600) = + reac_rate_local(600) 
  reac_source_local(37,600) = - reac_rate_local(600) 
  reac_source_local(21,601) = + reac_rate_local(601) 
  reac_source_local(29,601) = + reac_rate_local(601) 
  reac_source_local(33,601) = - reac_rate_local(601) 
  reac_source_local(37,601) = - reac_rate_local(601) 
  reac_source_local(21,602) = + reac_rate_local(602) * 2.d0
  reac_source_local(34,602) = - reac_rate_local(602) 
  reac_source_local(37,602) = - reac_rate_local(602) 
  reac_source_local(21,603) = + reac_rate_local(603) 
  reac_source_local(37,603) = - reac_rate_local(603) 
  reac_source_local(40,603) = + reac_rate_local(603) 
  reac_source_local(45,603) = - reac_rate_local(603) 
  reac_source_local(17,604) = - reac_rate_local(604) 
  reac_source_local(36,604) = - reac_rate_local(604) 
  reac_source_local(40,604) = + reac_rate_local(604) 
  reac_source_local(18,605) = - reac_rate_local(605) 
  reac_source_local(36,605) = - reac_rate_local(605) 
  reac_source_local(41,605) = + reac_rate_local(605) 
  reac_source_local(21,606) = + reac_rate_local(606) 
  reac_source_local(33,606) = - reac_rate_local(606) 
  reac_source_local(36,606) = - reac_rate_local(606) 
  reac_source_local(32,607) = + reac_rate_local(607) 
  reac_source_local(34,607) = - reac_rate_local(607) 
  reac_source_local(36,607) = - reac_rate_local(607) 
  reac_source_local(36,608) = - reac_rate_local(608) 
  reac_source_local(42,608) = + reac_rate_local(608) 
  reac_source_local(45,608) = - reac_rate_local(608) 
  reac_source_local(17,609) = - reac_rate_local(609) 
  reac_source_local(37,609) = - reac_rate_local(609) 
  reac_source_local(42,609) = + reac_rate_local(609) 
  reac_source_local(32,610) = + reac_rate_local(610) 
  reac_source_local(33,610) = - reac_rate_local(610) 
  reac_source_local(37,610) = - reac_rate_local(610) 
  reac_source_local(37,611) = - reac_rate_local(611) 
  reac_source_local(43,611) = + reac_rate_local(611) 
  reac_source_local(45,611) = - reac_rate_local(611) 
  reac_source_local(14,612) = + reac_rate_local(612) 
  reac_source_local(17,612) = - reac_rate_local(612) 
  reac_source_local(32,612) = + reac_rate_local(612) 
  reac_source_local(38,612) = - reac_rate_local(612) 
  reac_source_local(01,613) = + reac_rate_local(613) 
  reac_source_local(18,613) = - reac_rate_local(613) 
  reac_source_local(32,613) = + reac_rate_local(613) 
  reac_source_local(38,613) = - reac_rate_local(613) 
  reac_source_local(29,614) = + reac_rate_local(614) 
  reac_source_local(32,614) = + reac_rate_local(614) 
  reac_source_local(33,614) = - reac_rate_local(614) 
  reac_source_local(38,614) = - reac_rate_local(614) 
  reac_source_local(21,615) = + reac_rate_local(615) 
  reac_source_local(32,615) = + reac_rate_local(615) 
  reac_source_local(34,615) = - reac_rate_local(615) 
  reac_source_local(38,615) = - reac_rate_local(615) 
  reac_source_local(32,616) = + reac_rate_local(616) 
  reac_source_local(38,616) = - reac_rate_local(616) 
  reac_source_local(40,616) = + reac_rate_local(616) 
  reac_source_local(45,616) = - reac_rate_local(616) 
  reac_source_local(32,617) = + reac_rate_local(617) 
  reac_source_local(38,617) = - reac_rate_local(617) 
  reac_source_local(41,617) = + reac_rate_local(617) 
  reac_source_local(46,617) = - reac_rate_local(617) 
  reac_source_local(32,618) = + reac_rate_local(618) 
  reac_source_local(38,618) = - reac_rate_local(618) 
  reac_source_local(42,618) = + reac_rate_local(618) 
  reac_source_local(47,618) = - reac_rate_local(618) 
  reac_source_local(14,619) = + reac_rate_local(619) 
  reac_source_local(17,619) = - reac_rate_local(619) 
  reac_source_local(40,619) = + reac_rate_local(619) 
  reac_source_local(48,619) = - reac_rate_local(619) 
  reac_source_local(01,620) = + reac_rate_local(620) 
  reac_source_local(18,620) = - reac_rate_local(620) 
  reac_source_local(40,620) = + reac_rate_local(620) 
  reac_source_local(48,620) = - reac_rate_local(620) 
  reac_source_local(29,621) = + reac_rate_local(621) 
  reac_source_local(33,621) = - reac_rate_local(621) 
  reac_source_local(40,621) = + reac_rate_local(621) 
  reac_source_local(48,621) = - reac_rate_local(621) 
  reac_source_local(21,622) = + reac_rate_local(622) 
  reac_source_local(34,622) = - reac_rate_local(622) 
  reac_source_local(40,622) = + reac_rate_local(622) 
  reac_source_local(48,622) = - reac_rate_local(622) 
  reac_source_local(40,623) = + reac_rate_local(623) * 2.d0
  reac_source_local(45,623) = - reac_rate_local(623) 
  reac_source_local(48,623) = - reac_rate_local(623) 
  reac_source_local(40,624) = + reac_rate_local(624) 
  reac_source_local(41,624) = + reac_rate_local(624) 
  reac_source_local(46,624) = - reac_rate_local(624) 
  reac_source_local(48,624) = - reac_rate_local(624) 
  reac_source_local(40,625) = + reac_rate_local(625) 
  reac_source_local(42,625) = + reac_rate_local(625) 
  reac_source_local(47,625) = - reac_rate_local(625) 
  reac_source_local(48,625) = - reac_rate_local(625) 
  reac_source_local(14,626) = + reac_rate_local(626) 
  reac_source_local(17,626) = - reac_rate_local(626) 
  reac_source_local(41,626) = + reac_rate_local(626) 
  reac_source_local(49,626) = - reac_rate_local(626) 
  reac_source_local(01,627) = + reac_rate_local(627) 
  reac_source_local(18,627) = - reac_rate_local(627) 
  reac_source_local(41,627) = + reac_rate_local(627) 
  reac_source_local(49,627) = - reac_rate_local(627) 
  reac_source_local(29,628) = + reac_rate_local(628) 
  reac_source_local(33,628) = - reac_rate_local(628) 
  reac_source_local(41,628) = + reac_rate_local(628) 
  reac_source_local(49,628) = - reac_rate_local(628) 
  reac_source_local(21,629) = + reac_rate_local(629) 
  reac_source_local(34,629) = - reac_rate_local(629) 
  reac_source_local(41,629) = + reac_rate_local(629) 
  reac_source_local(49,629) = - reac_rate_local(629) 
  reac_source_local(40,630) = + reac_rate_local(630) 
  reac_source_local(41,630) = + reac_rate_local(630) 
  reac_source_local(45,630) = - reac_rate_local(630) 
  reac_source_local(49,630) = - reac_rate_local(630) 
  reac_source_local(41,631) = + reac_rate_local(631) * 2.d0
  reac_source_local(46,631) = - reac_rate_local(631) 
  reac_source_local(49,631) = - reac_rate_local(631) 
  reac_source_local(41,632) = + reac_rate_local(632) 
  reac_source_local(42,632) = + reac_rate_local(632) 
  reac_source_local(47,632) = - reac_rate_local(632) 
  reac_source_local(49,632) = - reac_rate_local(632) 
  reac_source_local(14,633) = + reac_rate_local(633) 
  reac_source_local(17,633) = - reac_rate_local(633) 
  reac_source_local(42,633) = + reac_rate_local(633) 
  reac_source_local(50,633) = - reac_rate_local(633) 
  reac_source_local(01,634) = + reac_rate_local(634) 
  reac_source_local(18,634) = - reac_rate_local(634) 
  reac_source_local(42,634) = + reac_rate_local(634) 
  reac_source_local(50,634) = - reac_rate_local(634) 
  reac_source_local(29,635) = + reac_rate_local(635) 
  reac_source_local(33,635) = - reac_rate_local(635) 
  reac_source_local(42,635) = + reac_rate_local(635) 
  reac_source_local(50,635) = - reac_rate_local(635) 
  reac_source_local(21,636) = + reac_rate_local(636) 
  reac_source_local(34,636) = - reac_rate_local(636) 
  reac_source_local(42,636) = + reac_rate_local(636) 
  reac_source_local(50,636) = - reac_rate_local(636) 
  reac_source_local(40,637) = + reac_rate_local(637) 
  reac_source_local(42,637) = + reac_rate_local(637) 
  reac_source_local(45,637) = - reac_rate_local(637) 
  reac_source_local(50,637) = - reac_rate_local(637) 
  reac_source_local(41,638) = + reac_rate_local(638) 
  reac_source_local(42,638) = + reac_rate_local(638) 
  reac_source_local(46,638) = - reac_rate_local(638) 
  reac_source_local(50,638) = - reac_rate_local(638) 
  reac_source_local(42,639) = + reac_rate_local(639) * 2.d0
  reac_source_local(47,639) = - reac_rate_local(639) 
  reac_source_local(50,639) = - reac_rate_local(639) 
  reac_source_local(14,640) = + reac_rate_local(640) 
  reac_source_local(17,640) = - reac_rate_local(640) 
  reac_source_local(43,640) = + reac_rate_local(640) 
  reac_source_local(51,640) = - reac_rate_local(640) 
  reac_source_local(01,641) = + reac_rate_local(641) 
  reac_source_local(18,641) = - reac_rate_local(641) 
  reac_source_local(43,641) = + reac_rate_local(641) 
  reac_source_local(51,641) = - reac_rate_local(641) 
  reac_source_local(29,642) = + reac_rate_local(642) 
  reac_source_local(33,642) = - reac_rate_local(642) 
  reac_source_local(43,642) = + reac_rate_local(642) 
  reac_source_local(51,642) = - reac_rate_local(642) 
  reac_source_local(21,643) = + reac_rate_local(643) 
  reac_source_local(34,643) = - reac_rate_local(643) 
  reac_source_local(43,643) = + reac_rate_local(643) 
  reac_source_local(51,643) = - reac_rate_local(643) 
  reac_source_local(40,644) = + reac_rate_local(644) 
  reac_source_local(43,644) = + reac_rate_local(644) 
  reac_source_local(45,644) = - reac_rate_local(644) 
  reac_source_local(51,644) = - reac_rate_local(644) 
  reac_source_local(41,645) = + reac_rate_local(645) 
  reac_source_local(43,645) = + reac_rate_local(645) 
  reac_source_local(46,645) = - reac_rate_local(645) 
  reac_source_local(51,645) = - reac_rate_local(645) 
  reac_source_local(42,646) = + reac_rate_local(646) 
  reac_source_local(43,646) = + reac_rate_local(646) 
  reac_source_local(47,646) = - reac_rate_local(646) 
  reac_source_local(51,646) = - reac_rate_local(646) 
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
  rrt(001) = rrt(001) * density(01) * density(53) 
  rrt(002) = rrt(002) * density(01) * density(53) 
  rrt(003) = rrt(003) * density(01) * density(53) 
  rrt(004) = rrt(004) * density(01) * density(53) 
  rrt(005) = rrt(005) * density(01) * density(53) 
  rrt(006) = rrt(006) * density(01) * density(53) 
  rrt(007) = rrt(007) * density(01) * density(53) 
  rrt(008) = rrt(008) * density(01) * density(53) 
  rrt(009) = rrt(009) * density(01) * density(53) 
  rrt(010) = rrt(010) * density(21) * density(53) 
  rrt(011) = rrt(011) * density(21) * density(53) 
  rrt(012) = rrt(012) * density(21) * density(53) 
  rrt(013) = rrt(013) * density(21) * density(53) 
  rrt(014) = rrt(014) * density(21) * density(53) 
  rrt(015) = rrt(015) * density(21) * density(53) 
  rrt(016) = rrt(016) * density(22) * density(53) 
  rrt(017) = rrt(017) * density(23) * density(53) 
  rrt(018) = rrt(018) * density(24) * density(53) 
  rrt(019) = rrt(019) * density(25) * density(53) 
  rrt(020) = rrt(020) * density(01) * density(02) 
  rrt(021) = rrt(021) * density(01) * density(03) 
  rrt(022) = rrt(022) * density(01) * density(04) 
  rrt(023) = rrt(023) * density(01) * density(05) 
  rrt(024) = rrt(024) * density(01) * density(06) 
  rrt(025) = rrt(025) * density(01) * density(07) 
  rrt(026) = rrt(026) * density(01) * density(08) 
  rrt(027) = rrt(027) * density(01) * density(09) 
  rrt(028) = rrt(028) * density(01)**2 
  rrt(029) = rrt(029) * density(01) * density(02) 
  rrt(030) = rrt(030) * density(01) * density(03) 
  rrt(031) = rrt(031) * density(01) * density(04) 
  rrt(032) = rrt(032) * density(01) * density(05) 
  rrt(033) = rrt(033) * density(01) * density(06) 
  rrt(034) = rrt(034) * density(01) * density(07) 
  rrt(035) = rrt(035) * density(01) * density(08) 
  rrt(036) = rrt(036) * density(02) * density(14) 
  rrt(037) = rrt(037) * density(03) * density(14) 
  rrt(038) = rrt(038) * density(04) * density(14) 
  rrt(039) = rrt(039) * density(05) * density(14) 
  rrt(040) = rrt(040) * density(06) * density(14) 
  rrt(041) = rrt(041) * density(07) * density(14) 
  rrt(042) = rrt(042) * density(08) * density(14) 
  rrt(043) = rrt(043) * density(09) * density(14) 
  rrt(044) = rrt(044) * density(01) * density(14) 
  rrt(045) = rrt(045) * density(02) * density(14) 
  rrt(046) = rrt(046) * density(03) * density(14) 
  rrt(047) = rrt(047) * density(04) * density(14) 
  rrt(048) = rrt(048) * density(05) * density(14) 
  rrt(049) = rrt(049) * density(06) * density(14) 
  rrt(050) = rrt(050) * density(07) * density(14) 
  rrt(051) = rrt(051) * density(08) * density(14) 
  rrt(052) = rrt(052) * density(02) * density(29) 
  rrt(053) = rrt(053) * density(03) * density(29) 
  rrt(054) = rrt(054) * density(04) * density(29) 
  rrt(055) = rrt(055) * density(05) * density(29) 
  rrt(056) = rrt(056) * density(06) * density(29) 
  rrt(057) = rrt(057) * density(07) * density(29) 
  rrt(058) = rrt(058) * density(08) * density(29) 
  rrt(059) = rrt(059) * density(09) * density(29) 
  rrt(060) = rrt(060) * density(01) * density(29) 
  rrt(061) = rrt(061) * density(02) * density(29) 
  rrt(062) = rrt(062) * density(03) * density(29) 
  rrt(063) = rrt(063) * density(04) * density(29) 
  rrt(064) = rrt(064) * density(05) * density(29) 
  rrt(065) = rrt(065) * density(06) * density(29) 
  rrt(066) = rrt(066) * density(07) * density(29) 
  rrt(067) = rrt(067) * density(08) * density(29) 
  rrt(068) = rrt(068) * density(21) * density(22) 
  rrt(069) = rrt(069) * density(21) * density(23) 
  rrt(070) = rrt(070) * density(21) * density(24) 
  rrt(071) = rrt(071) * density(21) * density(25) 
  rrt(072) = rrt(072) * density(21)**2 
  rrt(073) = rrt(073) * density(21) * density(22) 
  rrt(074) = rrt(074) * density(21) * density(23) 
  rrt(075) = rrt(075) * density(21) * density(24) 
  rrt(076) = rrt(076) * density(22) * density(29) 
  rrt(077) = rrt(077) * density(23) * density(29) 
  rrt(078) = rrt(078) * density(24) * density(29) 
  rrt(079) = rrt(079) * density(25) * density(29) 
  rrt(080) = rrt(080) * density(21) * density(29) 
  rrt(081) = rrt(081) * density(22) * density(29) 
  rrt(082) = rrt(082) * density(23) * density(29) 
  rrt(083) = rrt(083) * density(24) * density(29) 
  rrt(084) = rrt(084) * density(01) * density(53) 
  rrt(085) = rrt(085) * density(01) * density(53) 
  rrt(086) = rrt(086) * density(01) * density(53) 
  rrt(087) = rrt(087) * density(01) * density(53) 
  rrt(088) = rrt(088) * density(01) * density(53) 
  rrt(089) = rrt(089) * density(01) * density(53) 
  rrt(090) = rrt(090) * density(01) * density(53) 
  rrt(091) = rrt(091) * density(01) * density(53) 
  rrt(092) = rrt(092) * density(01) * density(53) 
  rrt(093) = rrt(093) * density(01) * density(53) 
  rrt(094) = rrt(094) * density(01) * density(53) 
  rrt(095) = rrt(095) * density(01) * density(53) 
  rrt(096) = rrt(096) * density(21) * density(53) 
  rrt(097) = rrt(097) * density(21) * density(53) 
  rrt(098) = rrt(098) * density(21) * density(53) 
  rrt(099) = rrt(099) * density(21) * density(53) 
  rrt(100) = rrt(100) * density(21) * density(53) 
  rrt(101) = rrt(101) * density(21) * density(53) 
  rrt(102) = rrt(102) * density(29) * density(53) 
  rrt(103) = rrt(103) * density(29) * density(53) 
  rrt(104) = rrt(104) * density(26) * density(53) 
  rrt(105) = rrt(105) * density(14) * density(53) 
  rrt(106) = rrt(106) * density(29) * density(53) 
  rrt(107) = rrt(107) * density(01) * density(53) 
  rrt(108) = rrt(108) * density(10) * density(53) 
  rrt(109) = rrt(109) * density(21) * density(53) 
  rrt(110) = rrt(110) * density(26) * density(53) 
  rrt(111) = rrt(111) * density(40) * density(53) 
  rrt(112) = rrt(112) * density(41) * density(53) 
  rrt(113) = rrt(113) * density(18) * density(53) 
  rrt(114) = rrt(114) * density(18) * density(53) 
  rrt(115) = rrt(115) * density(18) * density(53) 
  rrt(116) = rrt(116) * density(34) * density(53) 
  rrt(117) = rrt(117) * density(34) * density(53) 
  rrt(118) = rrt(118) * density(34) * density(53) 
  rrt(119) = rrt(119) * density(45) * density(53) 
  rrt(120) = rrt(120) * density(45) * density(53) 
  rrt(121) = rrt(121) * density(19) * density(53) 
  rrt(122) = rrt(122) * density(20) * density(53) 
  rrt(123) = rrt(123) * density(46) * density(53) 
  rrt(124) = rrt(124) * density(47) * density(53) 
  rrt(125) = rrt(125) * density(35) * density(53) 
  rrt(126) = rrt(126) * density(52) * density(53) 
  rrt(127) = rrt(127) * density(17) * density(53)**2 
  rrt(128) = rrt(128) * density(33) * density(53)**2 
  rrt(129) = rrt(129) * density(17) * density(53) 
  rrt(130) = rrt(130) * density(33) * density(53) 
  rrt(131) = rrt(131) * density(21) * density(53) 
  rrt(132) = rrt(132) * density(40) * density(53) 
  rrt(133) = rrt(133) * density(32) * density(53) 
  rrt(134) = rrt(134) * density(32) * density(53) 
  rrt(135) = rrt(135) * density(41) * density(53) 
  rrt(136) = rrt(136) * density(21)**2 * density(53) 
  rrt(137) = rrt(137) * density(42) * density(53) 
  rrt(138) = rrt(138) * density(21) * density(29) * density(53) 
  rrt(139) = rrt(139) * density(21) * density(29) * density(53) 
  rrt(140) = rrt(140) * density(32) * density(53) 
  rrt(141) = rrt(141) * density(40) * density(53) 
  rrt(142) = rrt(142) * density(41) * density(53) 
  rrt(143) = rrt(143) * density(01) * density(21) * density(53) 
  rrt(144) = rrt(144) * density(29) * density(36) 
  rrt(145) = rrt(145) * density(14) * density(36) 
  rrt(146) = rrt(146) * density(36) * density(40) 
  rrt(147) = rrt(147) * density(01) * density(36) 
  rrt(148) = rrt(148) * density(21) * density(36) 
  rrt(149) = rrt(149) * density(26) * density(36) 
  rrt(150) = rrt(150) * density(27) * density(36) 
  rrt(151) = rrt(151) * density(10) * density(36) 
  rrt(152) = rrt(152) * density(11) * density(36) 
  rrt(153) = rrt(153) * density(32) * density(36) 
  rrt(154) = rrt(154) * density(29) * density(37) 
  rrt(155) = rrt(155) * density(14) * density(37) 
  rrt(156) = rrt(156) * density(21) * density(37) 
  rrt(157) = rrt(157) * density(26) * density(37) 
  rrt(158) = rrt(158) * density(27) * density(37) 
  rrt(159) = rrt(159) * density(01) * density(37) 
  rrt(160) = rrt(160) * density(10) * density(37) 
  rrt(161) = rrt(161) * density(11) * density(37) 
  rrt(162) = rrt(162) * density(29) * density(38) 
  rrt(163) = rrt(163) * density(14) * density(48) 
  rrt(164) = rrt(164) * density(14) * density(38) 
  rrt(165) = rrt(165) * density(14) * density(49) 
  rrt(166) = rrt(166) * density(14) * density(50) 
  rrt(167) = rrt(167) * density(14) * density(51) 
  rrt(168) = rrt(168) * density(29) * density(48) 
  rrt(169) = rrt(169) * density(29) * density(49) 
  rrt(170) = rrt(170) * density(29) * density(50) 
  rrt(171) = rrt(171) * density(29) * density(51) 
  rrt(172) = rrt(172) * density(10) * density(38) 
  rrt(173) = rrt(173) * density(10) * density(48) 
  rrt(174) = rrt(174) * density(10) * density(49) 
  rrt(175) = rrt(175) * density(10) * density(50) 
  rrt(176) = rrt(176) * density(10) * density(51) 
  rrt(177) = rrt(177) * density(11) * density(38) 
  rrt(178) = rrt(178) * density(11) * density(48) 
  rrt(179) = rrt(179) * density(11) * density(49) 
  rrt(180) = rrt(180) * density(11) * density(50) 
  rrt(181) = rrt(181) * density(11) * density(51) 
  rrt(182) = rrt(182) * density(10) 
  rrt(183) = rrt(183) * density(11) 
  rrt(184) = rrt(184) * density(12) 
  rrt(185) = rrt(185) * density(13) 
  rrt(186) = rrt(186) * density(26) 
  rrt(187) = rrt(187) * density(27) 
  rrt(188) = rrt(188) * density(27) 
  rrt(189) = rrt(189) * density(28) 
  rrt(190) = rrt(190) * density(10) * density(29) 
  rrt(191) = rrt(191) * density(10) * density(29) 
  rrt(192) = rrt(192) * density(10) * density(14) 
  rrt(193) = rrt(193) * density(10) * density(14) 
  rrt(194) = rrt(194) * density(10) * density(21) 
  rrt(195) = rrt(195) * density(10) * density(21) 
  rrt(196) = rrt(196) * density(10) * density(21) 
  rrt(197) = rrt(197) * density(10) * density(21) 
  rrt(198) = rrt(198) * density(01) * density(10) 
  rrt(199) = rrt(199) * density(10) * density(40) 
  rrt(200) = rrt(200) * density(10) * density(41) 
  rrt(201) = rrt(201) * density(10) * density(42) 
  rrt(202) = rrt(202) * density(10)**2 
  rrt(203) = rrt(203) * density(10)**2 
  rrt(204) = rrt(204) * density(01) * density(11) 
  rrt(205) = rrt(205) * density(01) * density(11) 
  rrt(206) = rrt(206) * density(11) * density(21) 
  rrt(207) = rrt(207) * density(11) * density(40) 
  rrt(208) = rrt(208) * density(01) * density(13) 
  rrt(209) = rrt(209) * density(13) * density(21) 
  rrt(210) = rrt(210) * density(01) * density(12) 
  rrt(211) = rrt(211) * density(12) * density(21) 
  rrt(212) = rrt(212) * density(12) * density(40) 
  rrt(213) = rrt(213) * density(10) * density(12) 
  rrt(214) = rrt(214) * density(12)**2 
  rrt(215) = rrt(215) * density(01) * density(14)**2 
  rrt(216) = rrt(216) * density(14)**2 * density(21) 
  rrt(217) = rrt(217) * density(14)**2 * density(40) 
  rrt(218) = rrt(218) * density(14)**3 
  rrt(219) = rrt(219) * density(14)**2 * density(29) 
  rrt(220) = rrt(220) * density(01) * density(14)**2 
  rrt(221) = rrt(221) * density(14)**2 * density(21) 
  rrt(222) = rrt(222) * density(14)**2 * density(40) 
  rrt(223) = rrt(223) * density(14)**3 
  rrt(224) = rrt(224) * density(14)**2 * density(29) 
  rrt(225) = rrt(225) * density(15) * density(29) 
  rrt(226) = rrt(226) * density(15) * density(21) 
  rrt(227) = rrt(227) * density(15) * density(40) 
  rrt(228) = rrt(228) * density(15) * density(41) 
  rrt(229) = rrt(229) * density(01) * density(15) 
  rrt(230) = rrt(230) * density(14) * density(16) 
  rrt(231) = rrt(231) * density(16) * density(29) 
  rrt(232) = rrt(232) * density(14) * density(16) 
  rrt(233) = rrt(233) * density(01) * density(16) 
  rrt(234) = rrt(234) * density(15) * density(16) 
  rrt(235) = rrt(235) * density(16) * density(21) 
  rrt(236) = rrt(236) * density(16) * density(40) 
  rrt(237) = rrt(237) * density(26) * density(29) 
  rrt(238) = rrt(238) * density(14) * density(26) 
  rrt(239) = rrt(239) * density(21) * density(26) 
  rrt(240) = rrt(240) * density(01) * density(26) 
  rrt(241) = rrt(241) * density(26) * density(40) 
  rrt(242) = rrt(242) * density(26) * density(32) 
  rrt(243) = rrt(243) * density(26)**2 
  rrt(244) = rrt(244) * density(29) * density(32) 
  rrt(245) = rrt(245) * density(27) * density(29) 
  rrt(246) = rrt(246) * density(27) * density(29) 
  rrt(247) = rrt(247) * density(21) * density(27) 
  rrt(248) = rrt(248) * density(01) * density(27) 
  rrt(249) = rrt(249) * density(27) * density(40) 
  rrt(250) = rrt(250) * density(27) * density(32) 
  rrt(251) = rrt(251) * density(28) * density(29) 
  rrt(252) = rrt(252) * density(21) * density(28) 
  rrt(253) = rrt(253) * density(01) * density(28) 
  rrt(254) = rrt(254) * density(29) * density(30) 
  rrt(255) = rrt(255) * density(21) * density(30) 
  rrt(256) = rrt(256) * density(21) * density(30) 
  rrt(257) = rrt(257) * density(21) * density(30) 
  rrt(258) = rrt(258) * density(01) * density(30) 
  rrt(259) = rrt(259) * density(30) * density(32) 
  rrt(260) = rrt(260) * density(30) * density(32) 
  rrt(261) = rrt(261) * density(30) * density(40) 
  rrt(262) = rrt(262) * density(30) * density(41) 
  rrt(263) = rrt(263) * density(30) * density(41) 
  rrt(264) = rrt(264) * density(29) * density(31) 
  rrt(265) = rrt(265) * density(14) * density(31) 
  rrt(266) = rrt(266) * density(21) * density(31) 
  rrt(267) = rrt(267) * density(21) * density(31) 
  rrt(268) = rrt(268) * density(01) * density(31) 
  rrt(269) = rrt(269) * density(26) * density(31) 
  rrt(270) = rrt(270) * density(26) * density(31) 
  rrt(271) = rrt(271) * density(26) * density(31) 
  rrt(272) = rrt(272) * density(31) * density(40) 
  rrt(273) = rrt(273) * density(31) * density(40) 
  rrt(274) = rrt(274) * density(31) * density(32) 
  rrt(275) = rrt(275) * density(31) * density(32) 
  rrt(276) = rrt(276) * density(31) * density(41) 
  rrt(277) = rrt(277) * density(31) * density(41) 
  rrt(278) = rrt(278) * density(14) * density(40) 
  rrt(279) = rrt(279) * density(14) * density(21) 
  rrt(280) = rrt(280) * density(14) * density(42) 
  rrt(281) = rrt(281) * density(14) * density(42) 
  rrt(282) = rrt(282) * density(14) * density(42) 
  rrt(283) = rrt(283) * density(14) * density(42) 
  rrt(284) = rrt(284) * density(01) * density(29) 
  rrt(285) = rrt(285) * density(29) * density(40) 
  rrt(286) = rrt(286) * density(29) * density(40) 
  rrt(287) = rrt(287) * density(29) * density(41) 
  rrt(288) = rrt(288) * density(29) * density(41) 
  rrt(289) = rrt(289) * density(29) * density(42) 
  rrt(290) = rrt(290) * density(29) * density(43) 
  rrt(291) = rrt(291) * density(01) * density(21) 
  rrt(292) = rrt(292) * density(40)**2 
  rrt(293) = rrt(293) * density(40)**2 
  rrt(294) = rrt(294) * density(40)**2 
  rrt(295) = rrt(295) * density(21) * density(40) 
  rrt(296) = rrt(296) * density(32) * density(40) 
  rrt(297) = rrt(297) * density(40) * density(41) 
  rrt(298) = rrt(298) * density(40) * density(43) 
  rrt(299) = rrt(299) * density(21)**2 
  rrt(300) = rrt(300) * density(21) * density(42) 
  rrt(301) = rrt(301) * density(42)**2 
  rrt(302) = rrt(302) * density(42)**2 
  rrt(303) = rrt(303) * density(32) * density(42) 
  rrt(304) = rrt(304) * density(42) * density(43) 
  rrt(305) = rrt(305) * density(21) * density(43) 
  rrt(306) = rrt(306) * density(43)**2 
  rrt(307) = rrt(307) * density(14)**2 
  rrt(308) = rrt(308) * density(14) * density(29) 
  rrt(309) = rrt(309) * density(01)**2 
  rrt(310) = rrt(310) * density(01) * density(21) 
  rrt(311) = rrt(311) * density(01) * density(40) 
  rrt(312) = rrt(312) * density(01) * density(29) 
  rrt(313) = rrt(313) * density(01) * density(14) 
  rrt(314) = rrt(314) * density(01) * density(21) 
  rrt(315) = rrt(315) * density(21)**2 
  rrt(316) = rrt(316) * density(21) * density(29) 
  rrt(317) = rrt(317) * density(14) * density(21) 
  rrt(318) = rrt(318) * density(21) * density(40) 
  rrt(319) = rrt(319) * density(01) * density(40) 
  rrt(320) = rrt(320) * density(21) * density(40) 
  rrt(321) = rrt(321) * density(29) * density(40) 
  rrt(322) = rrt(322) * density(14) * density(40) 
  rrt(323) = rrt(323) * density(40)**2 
  rrt(324) = rrt(324) * density(01) * density(32) 
  rrt(325) = rrt(325) * density(21) * density(32) 
  rrt(326) = rrt(326) * density(14) * density(32) 
  rrt(327) = rrt(327) * density(29) * density(32) 
  rrt(328) = rrt(328) * density(01) * density(41) 
  rrt(329) = rrt(329) * density(21) * density(41) 
  rrt(330) = rrt(330) * density(40) * density(41) 
  rrt(331) = rrt(331) * density(41)**2 
  rrt(332) = rrt(332) * density(01) * density(42) 
  rrt(333) = rrt(333) * density(21) * density(42) 
  rrt(334) = rrt(334) * density(40) * density(42) 
  rrt(335) = rrt(335) * density(42)**2 
  rrt(336) = rrt(336) * density(01) * density(43) 
  rrt(337) = rrt(337) * density(21) * density(43) 
  rrt(338) = rrt(338) * density(40) * density(43) 
  rrt(339) = rrt(339) * density(14) * density(43) 
  rrt(340) = rrt(340) * density(29) * density(43) 
  rrt(341) = rrt(341) * density(01) * density(43) 
  rrt(342) = rrt(342) * density(21) * density(43) 
  rrt(343) = rrt(343) * density(40) * density(43) 
  rrt(344) = rrt(344) * density(14) * density(43) 
  rrt(345) = rrt(345) * density(29) * density(43) 
  rrt(346) = rrt(346) * density(44) 
  rrt(347) = rrt(347) * density(01) * density(14)**2 
  rrt(348) = rrt(348) * density(14)**2 * density(21) 
  rrt(349) = rrt(349) * density(14)**2 * density(40) 
  rrt(350) = rrt(350) * density(14)**3 
  rrt(351) = rrt(351) * density(14)**2 * density(29) 
  rrt(352) = rrt(352) * density(01) * density(29)**2 
  rrt(353) = rrt(353) * density(21) * density(29)**2 
  rrt(354) = rrt(354) * density(14) * density(29)**2 
  rrt(355) = rrt(355) * density(29)**3 
  rrt(356) = rrt(356) * density(29)**2 * density(40) 
  rrt(357) = rrt(357) * density(01) * density(14) * density(29) 
  rrt(358) = rrt(358) * density(14) * density(21) * density(29) 
  rrt(359) = rrt(359) * density(14)**2 * density(29) 
  rrt(360) = rrt(360) * density(14) * density(29)**2 
  rrt(361) = rrt(361) * density(14) * density(29) * density(40) 
  rrt(362) = rrt(362) * density(01) * density(21) * density(29) 
  rrt(363) = rrt(363) * density(21)**2 * density(29) 
  rrt(364) = rrt(364) * density(21) * density(29) * density(40) 
  rrt(365) = rrt(365) * density(14) * density(21) * density(29) 
  rrt(366) = rrt(366) * density(21) * density(29)**2 
  rrt(367) = rrt(367) * density(01) * density(29) 
  rrt(368) = rrt(368) * density(01) * density(29) * density(40) 
  rrt(369) = rrt(369) * density(21) * density(29) * density(40) 
  rrt(370) = rrt(370) * density(29) * density(40)**2 
  rrt(371) = rrt(371) * density(01) * density(29) * density(42) 
  rrt(372) = rrt(372) * density(21) * density(29) * density(42) 
  rrt(373) = rrt(373) * density(14) * density(29) * density(42) 
  rrt(374) = rrt(374) * density(29)**2 * density(42) 
  rrt(375) = rrt(375) * density(29) * density(40) * density(42) 
  rrt(376) = rrt(376) * density(42) * density(43) 
  rrt(377) = rrt(377) * density(17) * density(29) 
  rrt(378) = rrt(378) * density(17) * density(21) 
  rrt(379) = rrt(379) * density(17) * density(21) 
  rrt(380) = rrt(380) * density(17) * density(21) 
  rrt(381) = rrt(381) * density(17) * density(32) 
  rrt(382) = rrt(382) * density(17) * density(40) 
  rrt(383) = rrt(383) * density(17) * density(40) 
  rrt(384) = rrt(384) * density(17) * density(40) 
  rrt(385) = rrt(385) * density(17) * density(41) 
  rrt(386) = rrt(386) * density(01) * density(33) 
  rrt(387) = rrt(387) * density(21) * density(33) 
  rrt(388) = rrt(388) * density(32) * density(33) 
  rrt(389) = rrt(389) * density(33) * density(40) 
  rrt(390) = rrt(390) * density(33) * density(40) 
  rrt(391) = rrt(391) * density(15) * density(33) 
  rrt(392) = rrt(392) * density(33) * density(41) 
  rrt(393) = rrt(393) * density(33) * density(41) 
  rrt(394) = rrt(394) * density(33) * density(41) 
  rrt(395) = rrt(395) * density(33) * density(42) 
  rrt(396) = rrt(396) * density(18) * density(21) 
  rrt(397) = rrt(397) * density(18) * density(29) 
  rrt(398) = rrt(398) * density(18) * density(32) 
  rrt(399) = rrt(399) * density(14) * density(18) 
  rrt(400) = rrt(400) * density(18) * density(40) 
  rrt(401) = rrt(401) * density(18) * density(41) 
  rrt(402) = rrt(402) * density(18) * density(41) 
  rrt(403) = rrt(403) * density(01) * density(34) 
  rrt(404) = rrt(404) * density(14) * density(34) 
  rrt(405) = rrt(405) * density(34) * density(40) 
  rrt(406) = rrt(406) * density(34) * density(42) 
  rrt(407) = rrt(407) * density(34) * density(42) 
  rrt(408) = rrt(408) * density(19) * density(21) 
  rrt(409) = rrt(409) * density(19) * density(21) 
  rrt(410) = rrt(410) * density(14) * density(19) 
  rrt(411) = rrt(411) * density(19) * density(40) 
  rrt(412) = rrt(412) * density(19) * density(40) 
  rrt(413) = rrt(413) * density(40) * density(47) 
  rrt(414) = rrt(414) * density(40) * density(46) 
  rrt(415) = rrt(415) * density(01) * density(20) 
  rrt(416) = rrt(416) * density(20) * density(21) 
  rrt(417) = rrt(417) * density(20) * density(29) 
  rrt(418) = rrt(418) * density(14) * density(20) 
  rrt(419) = rrt(419) * density(20) * density(40) 
  rrt(420) = rrt(420) * density(01) * density(35) 
  rrt(421) = rrt(421) * density(21) * density(35) 
  rrt(422) = rrt(422) * density(26) * density(35) 
  rrt(423) = rrt(423) * density(27) * density(35) 
  rrt(424) = rrt(424) * density(29) * density(35) 
  rrt(425) = rrt(425) * density(35) * density(40) 
  rrt(426) = rrt(426) * density(01) * density(52) 
  rrt(427) = rrt(427) * density(21) * density(52) 
  rrt(428) = rrt(428) * density(01)**2 * density(17) 
  rrt(429) = rrt(429) * density(17) * density(29) 
  rrt(430) = rrt(430) * density(14) * density(17) 
  rrt(431) = rrt(431) * density(01) * density(33) 
  rrt(432) = rrt(432) * density(29) * density(33) 
  rrt(433) = rrt(433) * density(14) * density(33) 
  rrt(434) = rrt(434) * density(01)**2 * density(18) 
  rrt(435) = rrt(435) * density(01) * density(14) * density(18) 
  rrt(436) = rrt(436) * density(21)**2 * density(34) 
  rrt(437) = rrt(437) * density(01)**2 * density(34) 
  rrt(438) = rrt(438) * density(26) * density(36) 
  rrt(439) = rrt(439) * density(32) * density(36) 
  rrt(440) = rrt(440) * density(36) * density(42) 
  rrt(441) = rrt(441) * density(36) * density(41) 
  rrt(442) = rrt(442) * density(36) * density(41) 
  rrt(443) = rrt(443) * density(29) * density(37) 
  rrt(444) = rrt(444) * density(32) * density(37) 
  rrt(445) = rrt(445) * density(37) * density(42) 
  rrt(446) = rrt(446) * density(37) * density(43) 
  rrt(447) = rrt(447) * density(29) * density(38) 
  rrt(448) = rrt(448) * density(38) * density(40) 
  rrt(449) = rrt(449) * density(38) * density(40) 
  rrt(450) = rrt(450) * density(38) * density(42) 
  rrt(451) = rrt(451) * density(38) * density(42) 
  rrt(452) = rrt(452) * density(38) * density(43) 
  rrt(453) = rrt(453) * density(21) * density(48) 
  rrt(454) = rrt(454) * density(42) * density(48) 
  rrt(455) = rrt(455) * density(41) * density(48) 
  rrt(456) = rrt(456) * density(32) * density(50) 
  rrt(457) = rrt(457) * density(42) * density(50) 
  rrt(458) = rrt(458) * density(43) * density(50) 
  rrt(459) = rrt(459) * density(44) * density(50) 
  rrt(460) = rrt(460) * density(40) * density(51) 
  rrt(461) = rrt(461) * density(01) * density(39) 
  rrt(462) = rrt(462) * density(21) * density(39) 
  rrt(463) = rrt(463) * density(29) * density(39) 
  rrt(464) = rrt(464) * density(29) * density(39) 
  rrt(465) = rrt(465) * density(26) * density(39) 
  rrt(466) = rrt(466) * density(27) * density(39) 
  rrt(467) = rrt(467) * density(39) * density(40) 
  rrt(468) = rrt(468) * density(21) * density(36) 
  rrt(469) = rrt(469) * density(36) * density(40) 
  rrt(470) = rrt(470) * density(21) * density(37) 
  rrt(471) = rrt(471) * density(17) * density(36) 
  rrt(472) = rrt(472) * density(18) * density(36) 
  rrt(473) = rrt(473) * density(33) * density(36) 
  rrt(474) = rrt(474) * density(34) * density(36) 
  rrt(475) = rrt(475) * density(36) * density(45) 
  rrt(476) = rrt(476) * density(36) * density(46) 
  rrt(477) = rrt(477) * density(36) * density(47) 
  rrt(478) = rrt(478) * density(17) * density(37) 
  rrt(479) = rrt(479) * density(18) * density(37) 
  rrt(480) = rrt(480) * density(33) * density(37) 
  rrt(481) = rrt(481) * density(34) * density(37) 
  rrt(482) = rrt(482) * density(37) * density(45) 
  rrt(483) = rrt(483) * density(37) * density(46) 
  rrt(484) = rrt(484) * density(37) * density(47) 
  rrt(485) = rrt(485) * density(17) * density(38) 
  rrt(486) = rrt(486) * density(18) * density(38) 
  rrt(487) = rrt(487) * density(33) * density(38) 
  rrt(488) = rrt(488) * density(34) * density(38) 
  rrt(489) = rrt(489) * density(38) * density(45) 
  rrt(490) = rrt(490) * density(38) * density(46) 
  rrt(491) = rrt(491) * density(38) * density(47) 
  rrt(492) = rrt(492) * density(17) * density(48) 
  rrt(493) = rrt(493) * density(18) * density(48) 
  rrt(494) = rrt(494) * density(33) * density(48) 
  rrt(495) = rrt(495) * density(34) * density(48) 
  rrt(496) = rrt(496) * density(45) * density(48) 
  rrt(497) = rrt(497) * density(46) * density(48) 
  rrt(498) = rrt(498) * density(47) * density(48) 
  rrt(499) = rrt(499) * density(17) * density(49) 
  rrt(500) = rrt(500) * density(18) * density(49) 
  rrt(501) = rrt(501) * density(33) * density(49) 
  rrt(502) = rrt(502) * density(34) * density(49) 
  rrt(503) = rrt(503) * density(45) * density(49) 
  rrt(504) = rrt(504) * density(46) * density(49) 
  rrt(505) = rrt(505) * density(47) * density(49) 
  rrt(506) = rrt(506) * density(17) * density(50) 
  rrt(507) = rrt(507) * density(18) * density(50) 
  rrt(508) = rrt(508) * density(33) * density(50) 
  rrt(509) = rrt(509) * density(34) * density(50) 
  rrt(510) = rrt(510) * density(45) * density(50) 
  rrt(511) = rrt(511) * density(46) * density(50) 
  rrt(512) = rrt(512) * density(47) * density(50) 
  rrt(513) = rrt(513) * density(17) * density(51) 
  rrt(514) = rrt(514) * density(18) * density(51) 
  rrt(515) = rrt(515) * density(33) * density(51) 
  rrt(516) = rrt(516) * density(34) * density(51) 
  rrt(517) = rrt(517) * density(45) * density(51) 
  rrt(518) = rrt(518) * density(46) * density(51) 
  rrt(519) = rrt(519) * density(47) * density(51) 
  rrt(520) = rrt(520) * density(18) * density(36) 
  rrt(521) = rrt(521) * density(19) * density(36) 
  rrt(522) = rrt(522) * density(20) * density(36) 
  rrt(523) = rrt(523) * density(34) * density(36) 
  rrt(524) = rrt(524) * density(35) * density(36) 
  rrt(525) = rrt(525) * density(36) * density(45) 
  rrt(526) = rrt(526) * density(36) * density(46) 
  rrt(527) = rrt(527) * density(36) * density(47) 
  rrt(528) = rrt(528) * density(36) * density(52) 
  rrt(529) = rrt(529) * density(18) * density(37) 
  rrt(530) = rrt(530) * density(19) * density(37) 
  rrt(531) = rrt(531) * density(20) * density(37) 
  rrt(532) = rrt(532) * density(34) * density(37) 
  rrt(533) = rrt(533) * density(35) * density(37) 
  rrt(534) = rrt(534) * density(37) * density(45) 
  rrt(535) = rrt(535) * density(37) * density(46) 
  rrt(536) = rrt(536) * density(37) * density(47) 
  rrt(537) = rrt(537) * density(37) * density(52) 
  rrt(538) = rrt(538) * density(18) * density(38) 
  rrt(539) = rrt(539) * density(19) * density(38) 
  rrt(540) = rrt(540) * density(20) * density(38) 
  rrt(541) = rrt(541) * density(34) * density(38) 
  rrt(542) = rrt(542) * density(35) * density(38) 
  rrt(543) = rrt(543) * density(38) * density(45) 
  rrt(544) = rrt(544) * density(38) * density(46) 
  rrt(545) = rrt(545) * density(38) * density(47) 
  rrt(546) = rrt(546) * density(38) * density(52) 
  rrt(547) = rrt(547) * density(18) * density(48) 
  rrt(548) = rrt(548) * density(19) * density(48) 
  rrt(549) = rrt(549) * density(20) * density(48) 
  rrt(550) = rrt(550) * density(34) * density(48) 
  rrt(551) = rrt(551) * density(35) * density(48) 
  rrt(552) = rrt(552) * density(45) * density(48) 
  rrt(553) = rrt(553) * density(46) * density(48) 
  rrt(554) = rrt(554) * density(47) * density(48) 
  rrt(555) = rrt(555) * density(48) * density(52) 
  rrt(556) = rrt(556) * density(18) * density(49) 
  rrt(557) = rrt(557) * density(19) * density(49) 
  rrt(558) = rrt(558) * density(20) * density(49) 
  rrt(559) = rrt(559) * density(34) * density(49) 
  rrt(560) = rrt(560) * density(35) * density(49) 
  rrt(561) = rrt(561) * density(45) * density(49) 
  rrt(562) = rrt(562) * density(46) * density(49) 
  rrt(563) = rrt(563) * density(47) * density(49) 
  rrt(564) = rrt(564) * density(49) * density(52) 
  rrt(565) = rrt(565) * density(18) * density(50) 
  rrt(566) = rrt(566) * density(19) * density(50) 
  rrt(567) = rrt(567) * density(20) * density(50) 
  rrt(568) = rrt(568) * density(34) * density(50) 
  rrt(569) = rrt(569) * density(35) * density(50) 
  rrt(570) = rrt(570) * density(45) * density(50) 
  rrt(571) = rrt(571) * density(46) * density(50) 
  rrt(572) = rrt(572) * density(47) * density(50) 
  rrt(573) = rrt(573) * density(50) * density(52) 
  rrt(574) = rrt(574) * density(18) * density(51) 
  rrt(575) = rrt(575) * density(19) * density(51) 
  rrt(576) = rrt(576) * density(20) * density(51) 
  rrt(577) = rrt(577) * density(34) * density(51) 
  rrt(578) = rrt(578) * density(35) * density(51) 
  rrt(579) = rrt(579) * density(45) * density(51) 
  rrt(580) = rrt(580) * density(46) * density(51) 
  rrt(581) = rrt(581) * density(47) * density(51) 
  rrt(582) = rrt(582) * density(51) * density(52) 
  rrt(583) = rrt(583) * density(17) * density(39) 
  rrt(584) = rrt(584) * density(18) * density(39) 
  rrt(585) = rrt(585) * density(33) * density(39) 
  rrt(586) = rrt(586) * density(34) * density(39) 
  rrt(587) = rrt(587) * density(39) * density(45) 
  rrt(588) = rrt(588) * density(39) * density(46) 
  rrt(589) = rrt(589) * density(39) * density(47) 
  rrt(590) = rrt(590) * density(19) * density(39) 
  rrt(591) = rrt(591) * density(20) * density(39) 
  rrt(592) = rrt(592) * density(35) * density(39) 
  rrt(593) = rrt(593) * density(39) * density(52) 
  rrt(594) = rrt(594) * density(17) * density(36) 
  rrt(595) = rrt(595) * density(18) * density(36) 
  rrt(596) = rrt(596) * density(33) * density(36) 
  rrt(597) = rrt(597) * density(34) * density(36) 
  rrt(598) = rrt(598) * density(36) * density(45) 
  rrt(599) = rrt(599) * density(17) * density(37) 
  rrt(600) = rrt(600) * density(18) * density(37) 
  rrt(601) = rrt(601) * density(33) * density(37) 
  rrt(602) = rrt(602) * density(34) * density(37) 
  rrt(603) = rrt(603) * density(37) * density(45) 
  rrt(604) = rrt(604) * density(17) * density(36) 
  rrt(605) = rrt(605) * density(18) * density(36) 
  rrt(606) = rrt(606) * density(33) * density(36) 
  rrt(607) = rrt(607) * density(34) * density(36) 
  rrt(608) = rrt(608) * density(36) * density(45) 
  rrt(609) = rrt(609) * density(17) * density(37) 
  rrt(610) = rrt(610) * density(33) * density(37) 
  rrt(611) = rrt(611) * density(37) * density(45) 
  rrt(612) = rrt(612) * density(17) * density(38) 
  rrt(613) = rrt(613) * density(18) * density(38) 
  rrt(614) = rrt(614) * density(33) * density(38) 
  rrt(615) = rrt(615) * density(34) * density(38) 
  rrt(616) = rrt(616) * density(38) * density(45) 
  rrt(617) = rrt(617) * density(38) * density(46) 
  rrt(618) = rrt(618) * density(38) * density(47) 
  rrt(619) = rrt(619) * density(17) * density(48) 
  rrt(620) = rrt(620) * density(18) * density(48) 
  rrt(621) = rrt(621) * density(33) * density(48) 
  rrt(622) = rrt(622) * density(34) * density(48) 
  rrt(623) = rrt(623) * density(45) * density(48) 
  rrt(624) = rrt(624) * density(46) * density(48) 
  rrt(625) = rrt(625) * density(47) * density(48) 
  rrt(626) = rrt(626) * density(17) * density(49) 
  rrt(627) = rrt(627) * density(18) * density(49) 
  rrt(628) = rrt(628) * density(33) * density(49) 
  rrt(629) = rrt(629) * density(34) * density(49) 
  rrt(630) = rrt(630) * density(45) * density(49) 
  rrt(631) = rrt(631) * density(46) * density(49) 
  rrt(632) = rrt(632) * density(47) * density(49) 
  rrt(633) = rrt(633) * density(17) * density(50) 
  rrt(634) = rrt(634) * density(18) * density(50) 
  rrt(635) = rrt(635) * density(33) * density(50) 
  rrt(636) = rrt(636) * density(34) * density(50) 
  rrt(637) = rrt(637) * density(45) * density(50) 
  rrt(638) = rrt(638) * density(46) * density(50) 
  rrt(639) = rrt(639) * density(47) * density(50) 
  rrt(640) = rrt(640) * density(17) * density(51) 
  rrt(641) = rrt(641) * density(18) * density(51) 
  rrt(642) = rrt(642) * density(33) * density(51) 
  rrt(643) = rrt(643) * density(34) * density(51) 
  rrt(644) = rrt(644) * density(45) * density(51) 
  rrt(645) = rrt(645) * density(46) * density(51) 
  rrt(646) = rrt(646) * density(47) * density(51) 
  ydot(01) = -rrt(001)-rrt(002)-rrt(003)-rrt(004)-rrt(005)-rrt(006)-rrt(007)-rrt(008)-rrt(009)+rrt(020)-rrt(028)+rrt(036)-rrt(044)&
             +rrt(052)-rrt(060)-rrt(084)-rrt(085)-rrt(086)-rrt(087)-rrt(088)-rrt(089)-rrt(090)-rrt(091)-rrt(092)-rrt(093)-rrt(094)&
             -rrt(095)-rrt(107)+rrt(121)+  2.d0 * rrt(122)+rrt(123)+rrt(126)-rrt(147)+rrt(151)+rrt(152)+rrt(160)+rrt(161)+rrt(165)&
             +rrt(172)+rrt(173)+rrt(174)+rrt(175)+rrt(176)+rrt(177)+rrt(178)+rrt(179)+rrt(180)+rrt(181)+rrt(182)+rrt(184)+rrt(191)&
             +rrt(192)+rrt(193)+rrt(194)+rrt(195)+rrt(196)+rrt(198)+rrt(199)+rrt(200)+rrt(201)+rrt(202)+rrt(203)+rrt(205)+rrt(206)&
             +rrt(209)+rrt(211)+rrt(212)+rrt(227)+rrt(228)+rrt(263)+rrt(278)+rrt(280)+rrt(282)-rrt(284)+rrt(287)-rrt(291)+rrt(294)&
             +rrt(297)-rrt(309)-rrt(310)-rrt(311)-rrt(312)-rrt(313)+rrt(328)+rrt(329)+rrt(330)+rrt(331)+rrt(347)+rrt(348)+rrt(349)&
             +rrt(350)+rrt(351)-rrt(367)+rrt(384)+rrt(385)-rrt(386)+rrt(394)+rrt(396)+rrt(398)+rrt(399)+rrt(400)+rrt(401)+rrt(402)&
             -rrt(403)+rrt(408)+rrt(409)+rrt(410)+rrt(411)+rrt(412)+rrt(415)+  2.d0 * rrt(416)+  2.d0 * rrt(417)+  2.d0 * rrt(418)&
             +  2.d0 * rrt(419)-rrt(420)+rrt(426)+rrt(427)-rrt(428)-rrt(431)-rrt(434)-rrt(437)+rrt(455)+rrt(472)+rrt(479)+rrt(486)&
             +rrt(493)+rrt(500)+rrt(507)+rrt(514)+rrt(521)+  2.d0 * rrt(522)+rrt(526)+rrt(528)+rrt(530)+  2.d0 * rrt(531)+rrt(535)&
             +rrt(537)+rrt(539)+  2.d0 * rrt(540)+rrt(544)+rrt(546)+rrt(548)+  2.d0 * rrt(549)+rrt(553)+rrt(555)+rrt(557)&
             +  2.d0 * rrt(558)+rrt(562)+rrt(564)+rrt(566)+  2.d0 * rrt(567)+rrt(571)+rrt(573)+rrt(575)+  2.d0 * rrt(576)+rrt(580)&
             +rrt(582)+rrt(584)+rrt(590)+  2.d0 * rrt(591)+rrt(593)+rrt(595)+rrt(600)+rrt(613)+rrt(620)+rrt(627)+rrt(634)+rrt(641) 
  ydot(02) = +rrt(001)+rrt(002)-rrt(020)+rrt(021)+rrt(028)-rrt(029)-rrt(036)+rrt(037)+rrt(044)-rrt(045)-rrt(052)+rrt(053)+rrt(060)&
             -rrt(061) 
  ydot(03) = +rrt(003)-rrt(021)+rrt(022)+rrt(029)-rrt(030)-rrt(037)+rrt(038)+rrt(045)-rrt(046)-rrt(053)+rrt(054)+rrt(061)-rrt(062) 
  ydot(04) = +rrt(004)-rrt(022)+rrt(023)+rrt(030)-rrt(031)-rrt(038)+rrt(039)+rrt(046)-rrt(047)-rrt(054)+rrt(055)+rrt(062)-rrt(063) 
  ydot(05) = +rrt(005)-rrt(023)+rrt(024)+rrt(031)-rrt(032)-rrt(039)+rrt(040)+rrt(047)-rrt(048)-rrt(055)+rrt(056)+rrt(063)-rrt(064) 
  ydot(06) = +rrt(006)-rrt(024)+rrt(025)+rrt(032)-rrt(033)-rrt(040)+rrt(041)+rrt(048)-rrt(049)-rrt(056)+rrt(057)+rrt(064)-rrt(065) 
  ydot(07) = +rrt(007)-rrt(025)+rrt(026)+rrt(033)-rrt(034)-rrt(041)+rrt(042)+rrt(049)-rrt(050)-rrt(057)+rrt(058)+rrt(065)-rrt(066) 
  ydot(08) = +rrt(008)-rrt(026)+rrt(027)+rrt(034)-rrt(035)-rrt(042)+rrt(043)+rrt(050)-rrt(051)-rrt(058)+rrt(059)+rrt(066)-rrt(067) 
  ydot(09) = +rrt(009)-rrt(027)+rrt(035)-rrt(043)+rrt(051)-rrt(059)+rrt(067) 
  ydot(10) = +rrt(084)+rrt(085)+rrt(086)-rrt(108)-rrt(151)-rrt(160)-rrt(172)-rrt(173)-rrt(174)-rrt(175)-rrt(176)-rrt(182)+rrt(183)&
             -rrt(190)-rrt(191)-rrt(192)-rrt(193)-rrt(194)-rrt(195)-rrt(196)-rrt(197)-rrt(198)-rrt(199)-rrt(200)-rrt(201)&
             -  2.d0 * rrt(202)-  2.d0 * rrt(203)+rrt(204)+rrt(207)-rrt(213)+rrt(215)+rrt(216)+rrt(217)+rrt(218)+rrt(219)+rrt(236) 
  ydot(11) = +rrt(087)+rrt(088)+rrt(089)-rrt(152)-rrt(161)-rrt(177)-rrt(178)-rrt(179)-rrt(180)-rrt(181)-rrt(183)+rrt(185)+rrt(202)&
             -rrt(204)-rrt(205)-rrt(206)-rrt(207)+rrt(210)+rrt(220)+rrt(221)+rrt(222)+rrt(223)+rrt(224) 
  ydot(12) = +rrt(090)+rrt(091)+rrt(092)-rrt(184)+rrt(208)-rrt(210)-rrt(211)-rrt(212)-rrt(213)-  2.d0 * rrt(214) 
  ydot(13) = +rrt(093)+rrt(094)+rrt(095)-rrt(185)+rrt(203)-rrt(208)-rrt(209) 
  ydot(14) = -rrt(105)+  2.d0 * rrt(113)+rrt(114)+rrt(115)+rrt(119)+rrt(121)+rrt(127)+rrt(129)+rrt(132)+rrt(135)-rrt(145)-rrt(155)&
             -rrt(163)-rrt(164)-rrt(165)-rrt(166)-rrt(167)-rrt(193)+rrt(200)+rrt(212)-  2.d0 * rrt(215)-  2.d0 * rrt(216)&
             -  2.d0 * rrt(217)-  2.d0 * rrt(218)-  2.d0 * rrt(219)-  2.d0 * rrt(220)-  2.d0 * rrt(221)-  2.d0 * rrt(222)&
             -  2.d0 * rrt(223)-  2.d0 * rrt(224)+rrt(225)+rrt(229)+rrt(230)+rrt(231)+rrt(233)-rrt(238)+rrt(261)-rrt(278)-rrt(279)&
             -rrt(280)-rrt(281)-rrt(282)-rrt(283)+rrt(284)+rrt(285)+rrt(292)-  2.d0 * rrt(307)-rrt(308)+  2.d0 * rrt(309)&
             +  2.d0 * rrt(310)+  2.d0 * rrt(311)+  2.d0 * rrt(312)+  2.d0 * rrt(313)+rrt(319)+rrt(320)+rrt(321)+rrt(322)+rrt(323)&
             -  2.d0 * rrt(347)-  2.d0 * rrt(348)-  2.d0 * rrt(349)-  2.d0 * rrt(350)-  2.d0 * rrt(351)-rrt(357)-rrt(358)-rrt(359)&
             -rrt(360)-rrt(361)+rrt(377)+rrt(378)+rrt(382)+rrt(386)+rrt(390)+rrt(397)-rrt(399)+rrt(402)-rrt(404)+rrt(408)-rrt(410)&
             +rrt(411)-rrt(418)-rrt(430)+rrt(431)-rrt(433)-rrt(435)+rrt(471)+rrt(478)+rrt(485)+rrt(492)+rrt(499)+rrt(506)+rrt(513)&
             +  2.d0 * rrt(520)+rrt(521)+rrt(525)+rrt(527)+  2.d0 * rrt(529)+rrt(530)+rrt(534)+rrt(536)+  2.d0 * rrt(538)+rrt(539)&
             +rrt(543)+rrt(545)+  2.d0 * rrt(547)+rrt(548)+rrt(552)+rrt(554)+  2.d0 * rrt(556)+rrt(557)+rrt(561)+rrt(563)&
             +  2.d0 * rrt(565)+rrt(566)+rrt(570)+rrt(572)+  2.d0 * rrt(574)+rrt(575)+rrt(579)+rrt(581)+rrt(583)+rrt(590)+rrt(594)&
             +rrt(599)+rrt(612)+rrt(619)+rrt(626)+rrt(633)+rrt(640) 
  ydot(15) = +rrt(114)+rrt(120)+rrt(190)-rrt(225)-rrt(226)-rrt(227)-rrt(228)-rrt(229)+rrt(232)-rrt(234)-rrt(391) 
  ydot(16) = +rrt(115)+rrt(193)-rrt(230)-rrt(231)-rrt(232)-rrt(233)-rrt(234)-rrt(235)-rrt(236) 
  ydot(17) = +rrt(105)-rrt(127)-rrt(129)-rrt(377)-rrt(378)-rrt(379)-rrt(380)-rrt(381)-rrt(382)-rrt(383)-rrt(384)-rrt(385)+rrt(391)&
             +rrt(399)+rrt(418)-rrt(428)-rrt(429)-rrt(430)-rrt(471)-rrt(478)-rrt(485)-rrt(492)-rrt(499)-rrt(506)-rrt(513)-rrt(583)&
             -rrt(594)-rrt(599)-rrt(604)-rrt(609)-rrt(612)-rrt(619)-rrt(626)-rrt(633)-rrt(640) 
  ydot(18) = +rrt(107)+rrt(108)-rrt(113)-rrt(114)-rrt(115)+rrt(234)+rrt(307)+rrt(383)-rrt(396)-rrt(397)-rrt(398)-rrt(399)-rrt(400)&
             -rrt(401)-rrt(402)+rrt(410)+rrt(415)+rrt(430)-rrt(434)-rrt(435)-rrt(472)-rrt(479)-rrt(486)-rrt(493)-rrt(500)-rrt(507)&
             -rrt(514)-rrt(520)-rrt(529)-rrt(538)-rrt(547)-rrt(556)-rrt(565)-rrt(574)-rrt(584)-rrt(595)-rrt(600)-rrt(605)-rrt(613)&
             -rrt(620)-rrt(627)-rrt(634)-rrt(641) 
  ydot(19) = -rrt(121)-rrt(408)-rrt(409)-rrt(410)-rrt(411)-rrt(412)+rrt(428)+rrt(435)-rrt(521)-rrt(530)-rrt(539)-rrt(548)-rrt(557)&
             -rrt(566)-rrt(575)-rrt(590) 
  ydot(20) = -rrt(122)+rrt(213)+rrt(214)-rrt(415)-rrt(416)-rrt(417)-rrt(418)-rrt(419)+rrt(434)-rrt(522)-rrt(531)-rrt(540)-rrt(549)&
             -rrt(558)-rrt(567)-rrt(576)-rrt(591) 
  ydot(21) = -rrt(010)-rrt(011)-rrt(012)-rrt(013)-rrt(014)-rrt(015)+rrt(016)+rrt(017)+rrt(018)+rrt(019)+rrt(068)-rrt(072)+rrt(076)&
             -rrt(080)-rrt(096)-rrt(097)-rrt(098)-rrt(099)-rrt(100)-rrt(101)+rrt(104)-rrt(109)+  2.d0 * rrt(125)+rrt(126)-rrt(131)&
             +rrt(133)-rrt(136)-rrt(139)-rrt(143)+rrt(144)-rrt(148)+rrt(150)+  2.d0 * rrt(153)+rrt(156)+  2.d0 * rrt(157)&
             +  2.d0 * rrt(158)+rrt(159)+rrt(160)+rrt(161)+  2.d0 * rrt(162)+rrt(164)+rrt(170)+rrt(186)+rrt(188)+rrt(189)-rrt(194)&
             -rrt(195)-rrt(196)-rrt(197)-rrt(206)-rrt(209)-rrt(211)-rrt(226)-rrt(235)+rrt(237)+rrt(239)+rrt(240)+rrt(241)&
             +  2.d0 * rrt(242)+rrt(243)+rrt(244)+rrt(246)+  2.d0 * rrt(250)+rrt(251)-rrt(252)-rrt(256)-rrt(257)+rrt(259)&
             +  2.d0 * rrt(260)+rrt(261)+rrt(263)-rrt(267)+  2.d0 * rrt(274)+rrt(275)-rrt(279)+rrt(282)+rrt(285)+rrt(287)+rrt(289)&
             +rrt(290)-rrt(291)+rrt(294)-rrt(295)+rrt(296)-  2.d0 * rrt(299)-rrt(300)+rrt(301)+rrt(303)+rrt(304)-rrt(305)+rrt(306)&
             -rrt(314)-rrt(315)-rrt(316)-rrt(317)-rrt(318)+rrt(324)+rrt(325)+rrt(326)+rrt(327)+rrt(341)+rrt(342)+rrt(343)+rrt(344)&
             +rrt(345)+rrt(352)+rrt(353)+rrt(354)+rrt(355)+rrt(356)-rrt(362)-rrt(363)-rrt(364)-rrt(365)-rrt(366)-rrt(378)-rrt(379)&
             -rrt(380)+rrt(381)-rrt(387)+rrt(388)-rrt(396)+rrt(405)+rrt(407)-rrt(408)-rrt(409)-rrt(416)+rrt(420)+rrt(421)&
             +  2.d0 * rrt(422)+  2.d0 * rrt(423)+  2.d0 * rrt(425)-rrt(427)-rrt(436)+rrt(443)+rrt(444)+rrt(445)+rrt(446)+rrt(447)&
             +rrt(449)+rrt(451)-rrt(453)+rrt(456)+rrt(461)+rrt(462)+rrt(463)+  2.d0 * rrt(464)+  2.d0 * rrt(465)+  2.d0 * rrt(466)&
             +rrt(467)-rrt(468)-rrt(470)+rrt(474)+rrt(478)+rrt(479)+rrt(480)+  2.d0 * rrt(481)+rrt(482)+rrt(483)+rrt(484)+rrt(488)&
             +rrt(495)+rrt(502)+rrt(509)+rrt(516)+  2.d0 * rrt(524)+rrt(527)+rrt(528)+rrt(529)+rrt(530)+rrt(531)+rrt(532)&
             +  3.d0 * rrt(533)+rrt(534)+rrt(535)+  2.d0 * rrt(536)+  2.d0 * rrt(537)+  2.d0 * rrt(542)+rrt(545)+rrt(546)&
             +  2.d0 * rrt(551)+rrt(554)+rrt(555)+  2.d0 * rrt(560)+rrt(563)+rrt(564)+  2.d0 * rrt(569)+rrt(572)+rrt(573)&
             +  2.d0 * rrt(578)+rrt(581)+rrt(582)+  2.d0 * rrt(583)+  2.d0 * rrt(584)+  2.d0 * rrt(585)+  3.d0 * rrt(586)&
             +  2.d0 * rrt(587)+  2.d0 * rrt(588)+  2.d0 * rrt(589)+  2.d0 * rrt(590)+  2.d0 * rrt(591)+  4.d0 * rrt(592)
  ydot(21) = ydot(21) &
             +  3.d0 * rrt(593)+rrt(597)+rrt(599)+rrt(600)+rrt(601)+  2.d0 * rrt(602)+rrt(603)+rrt(606)+rrt(615)+rrt(622)+rrt(629)&
             +rrt(636)+rrt(643) 
  ydot(22) = +rrt(010)+rrt(011)-rrt(016)-rrt(068)+rrt(069)+rrt(072)-rrt(073)-rrt(076)+rrt(077)+rrt(080)-rrt(081) 
  ydot(23) = +rrt(012)+rrt(013)-rrt(017)-rrt(069)+rrt(070)+rrt(073)-rrt(074)-rrt(077)+rrt(078)+rrt(081)-rrt(082) 
  ydot(24) = +rrt(014)-rrt(018)-rrt(070)+rrt(071)+rrt(074)-rrt(075)-rrt(078)+rrt(079)+rrt(082)-rrt(083) 
  ydot(25) = +rrt(015)-rrt(019)-rrt(071)+rrt(075)-rrt(079)+rrt(083) 
  ydot(26) = +rrt(096)-rrt(104)-rrt(110)-rrt(149)-rrt(157)-rrt(186)+rrt(187)+rrt(195)-rrt(237)-rrt(238)-rrt(239)-rrt(240)-rrt(241)&
             -rrt(242)-  2.d0 * rrt(243)+rrt(244)+rrt(245)+rrt(247)+rrt(248)+rrt(249)+rrt(256)-rrt(269)-rrt(270)-rrt(271)-rrt(422)&
             -rrt(438)-rrt(465) 
  ydot(27) = +rrt(097)-rrt(150)-rrt(158)-rrt(187)-rrt(188)+rrt(196)+rrt(243)-rrt(245)-rrt(246)-rrt(247)-rrt(248)-rrt(249)-rrt(250)&
             +  2.d0 * rrt(252)+rrt(253)+rrt(257)+rrt(270)-rrt(423)-rrt(466) 
  ydot(28) = +rrt(098)-rrt(189)-rrt(251)-rrt(252)-rrt(253)+rrt(269) 
  ydot(29) = +  2.d0 * rrt(099)+rrt(100)+rrt(101)-rrt(102)-rrt(103)-rrt(106)+  2.d0 * rrt(116)+rrt(117)+rrt(118)+rrt(119)+rrt(120)&
             +rrt(123)+rrt(124)+rrt(128)+rrt(130)+rrt(131)+rrt(134)-rrt(138)-rrt(144)+rrt(150)+rrt(151)+rrt(152)-rrt(154)-rrt(162)&
             -rrt(168)-rrt(169)-rrt(170)-rrt(171)-rrt(190)-rrt(191)+rrt(194)+rrt(197)+rrt(201)+  2.d0 * rrt(206)+rrt(209)&
             +  2.d0 * rrt(211)+rrt(212)-rrt(225)+rrt(226)+rrt(227)+rrt(235)+rrt(236)+rrt(238)-rrt(244)-rrt(246)+rrt(250)-rrt(251)&
             +rrt(254)+rrt(255)+rrt(256)+rrt(257)+rrt(258)+  2.d0 * rrt(259)+rrt(265)+  3.d0 * rrt(267)+rrt(268)+rrt(269)&
             +  3.d0 * rrt(271)+rrt(272)+rrt(275)+rrt(276)+rrt(278)+rrt(279)+  2.d0 * rrt(280)+rrt(281)-rrt(284)-rrt(285)-rrt(286)&
             -rrt(287)-rrt(288)-rrt(289)-rrt(290)+rrt(291)+rrt(293)+rrt(295)+rrt(299)-rrt(308)+  2.d0 * rrt(314)+  2.d0 * rrt(315)&
             +  2.d0 * rrt(316)+  2.d0 * rrt(317)+  2.d0 * rrt(318)+rrt(319)+rrt(320)+rrt(321)+rrt(322)+rrt(323)+rrt(324)+rrt(325)&
             +rrt(326)+rrt(327)+rrt(328)+rrt(329)+rrt(330)+rrt(331)+rrt(332)+rrt(333)+rrt(334)+rrt(335)+rrt(336)+rrt(337)+rrt(338)&
             +rrt(339)+rrt(340)-  2.d0 * rrt(352)-  2.d0 * rrt(353)-  2.d0 * rrt(354)-  2.d0 * rrt(355)-  2.d0 * rrt(356)-rrt(357)&
             -rrt(358)-rrt(359)-rrt(360)-rrt(361)-rrt(362)-rrt(363)-rrt(364)-rrt(365)-rrt(366)-rrt(367)-rrt(368)-rrt(369)-rrt(370)&
             -rrt(371)-rrt(372)-rrt(373)-rrt(374)-rrt(375)-rrt(377)+rrt(379)+rrt(383)+rrt(387)+rrt(389)+rrt(391)+rrt(393)+rrt(395)&
             -rrt(397)+rrt(398)+rrt(404)-rrt(417)-rrt(424)-rrt(429)-rrt(432)+rrt(438)+rrt(439)+rrt(440)+rrt(442)-rrt(443)-rrt(447)&
             +rrt(448)-rrt(463)-rrt(464)+rrt(471)+rrt(472)+  2.d0 * rrt(473)+rrt(474)+rrt(475)+rrt(476)+rrt(477)+rrt(480)+rrt(487)&
             +rrt(494)+rrt(501)+rrt(508)+rrt(515)+rrt(520)+rrt(521)+rrt(522)+  3.d0 * rrt(523)+rrt(524)+  2.d0 * rrt(525)&
             +  2.d0 * rrt(526)+rrt(527)+rrt(528)+  2.d0 * rrt(532)+rrt(534)+rrt(535)+  2.d0 * rrt(541)+rrt(543)+rrt(544)&
             +  2.d0 * rrt(550)+rrt(552)+rrt(553)+  2.d0 * rrt(559)+rrt(561)+rrt(562)+  2.d0 * rrt(568)+rrt(570)+rrt(571)&
             +  2.d0 * rrt(577)+rrt(579)+rrt(580)+rrt(585)+rrt(594)+rrt(595)+  2.d0 * rrt(596)+rrt(597)+rrt(598)+rrt(601)+rrt(614)&
             +rrt(621)+rrt(628)+rrt(635)+rrt(642) 
  ydot(30) = +rrt(100)+rrt(102)+rrt(117)+rrt(194)+rrt(225)+rrt(242)+rrt(246)-rrt(254)-rrt(255)-rrt(256)-rrt(257)-rrt(258)-rrt(259)&
             -rrt(260)-rrt(261)-rrt(262)-rrt(263)+rrt(264)+rrt(266)+rrt(270)+rrt(273)+rrt(275)+rrt(277) 
  ydot(31) = +rrt(101)+rrt(103)+rrt(118)+rrt(191)+rrt(209)+rrt(251)-rrt(264)-rrt(265)-rrt(266)-rrt(267)-rrt(268)-rrt(269)-rrt(270)&
             -rrt(271)-rrt(272)-rrt(273)-rrt(274)-rrt(275)-rrt(276)-rrt(277) 
  ydot(32) = -rrt(133)-rrt(134)-rrt(140)+rrt(148)+rrt(149)-rrt(153)+rrt(154)+rrt(171)+rrt(172)+rrt(177)-rrt(242)-rrt(244)-rrt(250)&
             -rrt(259)-rrt(260)-rrt(274)-rrt(275)-rrt(296)+rrt(299)+rrt(300)-rrt(303)+rrt(305)-rrt(324)-rrt(325)-rrt(326)-rrt(327)&
             +rrt(362)+rrt(363)+rrt(364)+rrt(365)+rrt(366)-rrt(381)-rrt(388)-rrt(398)+rrt(406)+rrt(424)-rrt(439)-rrt(444)+rrt(450)&
             +rrt(452)-rrt(456)+rrt(485)+rrt(486)+rrt(487)+rrt(488)+rrt(489)+rrt(490)+rrt(491)+rrt(538)+rrt(539)+rrt(540)+rrt(541)&
             +rrt(542)+rrt(543)+rrt(544)+rrt(545)+rrt(546)+rrt(607)+rrt(610)+rrt(612)+rrt(613)+rrt(614)+rrt(615)+rrt(616)+rrt(617)&
             +rrt(618) 
  ydot(33) = +rrt(106)-rrt(128)-rrt(130)+rrt(377)+rrt(380)+rrt(384)-rrt(386)-rrt(387)-rrt(388)-rrt(389)-rrt(390)-rrt(391)-rrt(392)&
             -rrt(393)-rrt(394)-rrt(395)+rrt(417)-rrt(431)-rrt(432)-rrt(433)-rrt(473)-rrt(480)-rrt(487)-rrt(494)-rrt(501)-rrt(508)&
             -rrt(515)-rrt(585)-rrt(596)-rrt(601)-rrt(606)-rrt(610)-rrt(614)-rrt(621)-rrt(628)-rrt(635)-rrt(642) 
  ydot(34) = +rrt(109)+rrt(110)-rrt(116)-rrt(117)-rrt(118)+rrt(378)+rrt(387)+rrt(388)+rrt(390)+rrt(394)+rrt(396)+rrt(398)-rrt(403)&
             -rrt(404)-rrt(405)-rrt(406)-rrt(407)+rrt(408)+rrt(416)+rrt(421)+rrt(422)+rrt(423)+rrt(424)+rrt(426)+rrt(432)-rrt(436)&
             -rrt(437)-rrt(474)-rrt(481)-rrt(488)-rrt(495)-rrt(502)-rrt(509)-rrt(516)-rrt(523)-rrt(532)-rrt(541)-rrt(550)-rrt(559)&
             -rrt(568)-rrt(577)-rrt(586)-rrt(597)-rrt(602)-rrt(607)-rrt(615)-rrt(622)-rrt(629)-rrt(636)-rrt(643) 
  ydot(35) = -rrt(125)-rrt(420)-rrt(421)-rrt(422)-rrt(423)-rrt(424)-rrt(425)+rrt(427)+rrt(436)-rrt(524)-rrt(533)-rrt(542)-rrt(551)&
             -rrt(560)-rrt(569)-rrt(578)-rrt(592) 
  ydot(36) = +rrt(131)+rrt(132)+rrt(133)+rrt(137)+rrt(138)-rrt(144)-rrt(145)-rrt(146)-rrt(147)-rrt(148)-rrt(149)-rrt(150)-rrt(151)&
             -rrt(152)-rrt(153)-rrt(438)-rrt(439)-rrt(440)-rrt(441)-rrt(442)+rrt(443)+rrt(464)-rrt(468)-rrt(469)-rrt(471)-rrt(472)&
             -rrt(473)-rrt(474)-rrt(475)-rrt(476)-rrt(477)-rrt(520)-rrt(521)-rrt(522)-rrt(523)-rrt(524)-rrt(525)-rrt(526)-rrt(527)&
             -rrt(528)-rrt(594)-rrt(595)-rrt(596)-rrt(597)-rrt(598)-rrt(604)-rrt(605)-rrt(606)-rrt(607)-rrt(608) 
  ydot(37) = +rrt(134)+rrt(136)+rrt(139)+rrt(143)-rrt(154)-rrt(155)-rrt(156)-rrt(157)-rrt(158)-rrt(159)-rrt(160)-rrt(161)+rrt(438)&
             -rrt(443)-rrt(444)-rrt(445)-rrt(446)+rrt(447)+rrt(453)+rrt(461)+rrt(462)+rrt(465)+rrt(466)-rrt(470)-rrt(478)-rrt(479)&
             -rrt(480)-rrt(481)-rrt(482)-rrt(483)-rrt(484)-rrt(529)-rrt(530)-rrt(531)-rrt(532)-rrt(533)-rrt(534)-rrt(535)-rrt(536)&
             -rrt(537)-rrt(599)-rrt(600)-rrt(601)-rrt(602)-rrt(603)-rrt(609)-rrt(610)-rrt(611) 
  ydot(38) = +rrt(140)-rrt(162)-rrt(164)-rrt(172)-rrt(177)+rrt(439)+rrt(444)-rrt(447)-rrt(448)-rrt(449)-rrt(450)-rrt(451)-rrt(452)&
             +rrt(463)+rrt(468)-rrt(485)-rrt(486)-rrt(487)-rrt(488)-rrt(489)-rrt(490)-rrt(491)-rrt(538)-rrt(539)-rrt(540)-rrt(541)&
             -rrt(542)-rrt(543)-rrt(544)-rrt(545)-rrt(546)-rrt(612)-rrt(613)-rrt(614)-rrt(615)-rrt(616)-rrt(617)-rrt(618) 
  ydot(39) = -rrt(461)-rrt(462)-rrt(463)-rrt(464)-rrt(465)-rrt(466)-rrt(467)+rrt(470)-rrt(583)-rrt(584)-rrt(585)-rrt(586)-rrt(587)&
             -rrt(588)-rrt(589)-rrt(590)-rrt(591)-rrt(592)-rrt(593) 
  ydot(40) = -rrt(111)+rrt(124)-rrt(132)+rrt(137)-rrt(141)+rrt(145)-rrt(146)+rrt(164)+rrt(165)+  2.d0 * rrt(166)+rrt(167)&
             +  2.d0 * rrt(169)+rrt(170)+rrt(171)+rrt(173)+rrt(178)+rrt(190)+rrt(200)+rrt(201)-rrt(212)+rrt(226)-rrt(227)+rrt(228)&
             +rrt(235)-rrt(236)+rrt(238)-rrt(261)+  2.d0 * rrt(262)-rrt(278)+rrt(279)+  2.d0 * rrt(283)+rrt(284)-rrt(285)-rrt(286)&
             +  2.d0 * rrt(288)+rrt(289)-  2.d0 * rrt(292)-  2.d0 * rrt(293)-  2.d0 * rrt(294)-rrt(295)-rrt(296)-rrt(297)-rrt(298)&
             +rrt(300)+  2.d0 * rrt(301)+rrt(302)+rrt(304)-rrt(319)-rrt(320)-rrt(321)-rrt(322)-rrt(323)+rrt(332)+rrt(333)+rrt(334)&
             +rrt(335)+rrt(341)+rrt(342)+rrt(343)+rrt(344)+rrt(345)+rrt(357)+rrt(358)+rrt(359)+rrt(360)+rrt(361)-rrt(368)-rrt(369)&
             -rrt(370)+rrt(380)-rrt(382)-rrt(383)-rrt(384)-rrt(389)-rrt(390)+rrt(392)-rrt(400)+rrt(403)-rrt(405)-rrt(411)-rrt(412)&
             -rrt(413)-rrt(414)-rrt(419)-rrt(425)+rrt(441)-rrt(448)-rrt(449)+rrt(453)+rrt(454)+rrt(457)-rrt(460)-rrt(467)-rrt(469)&
             +rrt(475)+rrt(482)+rrt(489)+rrt(492)+rrt(493)+rrt(494)+rrt(495)+  2.d0 * rrt(496)+rrt(497)+rrt(498)+rrt(503)+rrt(510)&
             +rrt(517)+rrt(547)+rrt(548)+rrt(549)+rrt(550)+rrt(551)+rrt(552)+rrt(553)+rrt(554)+rrt(555)+rrt(587)+rrt(598)+rrt(603)&
             +rrt(604)+rrt(616)+rrt(619)+rrt(620)+rrt(621)+rrt(622)+  2.d0 * rrt(623)+rrt(624)+rrt(625)+rrt(630)+rrt(637)+rrt(644) 
  ydot(41) = -rrt(112)-rrt(135)-rrt(142)+rrt(147)+rrt(163)+rrt(174)+rrt(179)+rrt(197)-rrt(200)-rrt(228)-rrt(262)-rrt(263)+rrt(281)&
             -rrt(287)-rrt(288)+rrt(291)+rrt(293)-rrt(297)-rrt(328)-rrt(329)-rrt(330)-rrt(331)+rrt(367)-rrt(385)-rrt(392)-rrt(393)&
             -rrt(394)-rrt(401)-rrt(402)+rrt(414)-rrt(441)-rrt(442)-rrt(455)+rrt(476)+rrt(483)+rrt(490)+rrt(497)+rrt(499)+rrt(500)&
             +rrt(501)+rrt(502)+rrt(503)+  2.d0 * rrt(504)+rrt(505)+rrt(511)+rrt(518)+rrt(556)+rrt(557)+rrt(558)+rrt(559)+rrt(560)&
             +rrt(561)+rrt(562)+rrt(563)+rrt(564)+rrt(588)+rrt(605)+rrt(617)+rrt(624)+rrt(626)+rrt(627)+rrt(628)+rrt(629)+rrt(630)&
             +  2.d0 * rrt(631)+rrt(632)+rrt(638)+rrt(645) 
  ydot(42) = -rrt(137)+rrt(146)+rrt(155)+rrt(167)+rrt(168)+rrt(175)+rrt(180)-rrt(201)-rrt(280)-rrt(281)-rrt(282)-rrt(283)+rrt(286)&
             -rrt(289)+rrt(290)+rrt(292)+rrt(295)+rrt(296)+rrt(297)+  2.d0 * rrt(298)-rrt(300)-  2.d0 * rrt(301)-  2.d0 * rrt(302)&
             -rrt(303)+rrt(305)+  2.d0 * rrt(306)-rrt(332)-rrt(333)-rrt(334)-rrt(335)+rrt(336)+rrt(337)+rrt(338)+rrt(339)+rrt(340)&
             +rrt(346)+rrt(368)+rrt(369)+rrt(370)-rrt(371)-rrt(372)-rrt(373)-rrt(374)-rrt(375)-rrt(376)-rrt(395)-rrt(406)-rrt(407)&
             +rrt(413)-rrt(440)-rrt(445)-rrt(450)-rrt(451)-rrt(454)-rrt(457)+rrt(458)+  2.d0 * rrt(459)+rrt(460)+rrt(477)+rrt(484)&
             +rrt(491)+rrt(498)+rrt(505)+rrt(506)+rrt(507)+rrt(508)+rrt(509)+rrt(510)+rrt(511)+  2.d0 * rrt(512)+rrt(519)+rrt(565)&
             +rrt(566)+rrt(567)+rrt(568)+rrt(569)+rrt(570)+rrt(571)+rrt(572)+rrt(573)+rrt(589)+rrt(608)+rrt(609)+rrt(618)+rrt(625)&
             +rrt(632)+rrt(633)+rrt(634)+rrt(635)+rrt(636)+rrt(637)+rrt(638)+  2.d0 * rrt(639)+rrt(646) 
  ydot(43) = +rrt(176)+rrt(181)-rrt(290)-rrt(298)+rrt(302)+rrt(303)-rrt(304)-rrt(305)-  2.d0 * rrt(306)-rrt(336)-rrt(337)-rrt(338)&
             -rrt(339)-rrt(340)-rrt(341)-rrt(342)-rrt(343)-rrt(344)-rrt(345)+rrt(346)+rrt(371)+rrt(372)+rrt(373)+rrt(374)+rrt(375)&
             -rrt(376)-rrt(446)-rrt(452)-rrt(458)+rrt(513)+rrt(514)+rrt(515)+rrt(516)+rrt(517)+rrt(518)+rrt(519)+rrt(574)+rrt(575)&
             +rrt(576)+rrt(577)+rrt(578)+rrt(579)+rrt(580)+rrt(581)+rrt(582)+rrt(611)+rrt(640)+rrt(641)+rrt(642)+rrt(643)+rrt(644)&
             +rrt(645)+rrt(646) 
  ydot(44) = -rrt(346)+rrt(376)-rrt(459) 
  ydot(45) = +rrt(111)-rrt(119)-rrt(120)+rrt(308)+rrt(379)+rrt(381)+rrt(382)+rrt(385)+rrt(386)+rrt(389)+rrt(392)+rrt(397)+rrt(400)&
             +rrt(402)+rrt(403)+rrt(404)+rrt(405)+rrt(406)+rrt(411)+rrt(413)+rrt(414)+rrt(419)+rrt(425)+rrt(429)+rrt(431)+rrt(433)&
             -rrt(475)-rrt(482)-rrt(489)-rrt(496)-rrt(503)-rrt(510)-rrt(517)-rrt(525)-rrt(534)-rrt(543)-rrt(552)-rrt(561)-rrt(570)&
             -rrt(579)-rrt(587)-rrt(598)-rrt(603)-rrt(608)-rrt(611)-rrt(616)-rrt(623)-rrt(630)-rrt(637)-rrt(644) 
  ydot(46) = +rrt(112)-rrt(123)+rrt(393)+rrt(401)+rrt(412)-rrt(414)-rrt(476)-rrt(483)-rrt(490)-rrt(497)-rrt(504)-rrt(511)-rrt(518)&
             -rrt(526)-rrt(535)-rrt(544)-rrt(553)-rrt(562)-rrt(571)-rrt(580)-rrt(588)-rrt(617)-rrt(624)-rrt(631)-rrt(638)-rrt(645) 
  ydot(47) = -rrt(124)+rrt(395)+rrt(407)+rrt(409)-rrt(413)-rrt(477)-rrt(484)-rrt(491)-rrt(498)-rrt(505)-rrt(512)-rrt(519)-rrt(527)&
             -rrt(536)-rrt(545)-rrt(554)-rrt(563)-rrt(572)-rrt(581)-rrt(589)-rrt(618)-rrt(625)-rrt(632)-rrt(639)-rrt(646) 
  ydot(48) = +rrt(135)+rrt(141)-rrt(163)-rrt(168)-rrt(173)-rrt(178)+rrt(441)-rrt(453)-rrt(454)-rrt(455)-rrt(492)-rrt(493)-rrt(494)&
             -rrt(495)-rrt(496)-rrt(497)-rrt(498)-rrt(547)-rrt(548)-rrt(549)-rrt(550)-rrt(551)-rrt(552)-rrt(553)-rrt(554)-rrt(555)&
             -rrt(619)-rrt(620)-rrt(621)-rrt(622)-rrt(623)-rrt(624)-rrt(625) 
  ydot(49) = +rrt(142)-rrt(165)-rrt(169)-rrt(174)-rrt(179)+rrt(442)-rrt(499)-rrt(500)-rrt(501)-rrt(502)-rrt(503)-rrt(504)-rrt(505)&
             -rrt(556)-rrt(557)-rrt(558)-rrt(559)-rrt(560)-rrt(561)-rrt(562)-rrt(563)-rrt(564)-rrt(626)-rrt(627)-rrt(628)-rrt(629)&
             -rrt(630)-rrt(631)-rrt(632) 
  ydot(50) = -rrt(166)-rrt(170)-rrt(175)-rrt(180)+rrt(440)+rrt(445)+rrt(449)+rrt(450)+rrt(454)+rrt(455)-rrt(456)-rrt(457)-rrt(458)&
             -rrt(459)+rrt(460)+rrt(469)-rrt(506)-rrt(507)-rrt(508)-rrt(509)-rrt(510)-rrt(511)-rrt(512)-rrt(565)-rrt(566)-rrt(567)&
             -rrt(568)-rrt(569)-rrt(570)-rrt(571)-rrt(572)-rrt(573)-rrt(633)-rrt(634)-rrt(635)-rrt(636)-rrt(637)-rrt(638)-rrt(639) 
  ydot(51) = -rrt(167)-rrt(171)-rrt(176)-rrt(181)+rrt(446)+rrt(448)+rrt(451)+rrt(452)+rrt(456)+rrt(457)+rrt(458)+rrt(459)-rrt(460)&
             +rrt(467)-rrt(513)-rrt(514)-rrt(515)-rrt(516)-rrt(517)-rrt(518)-rrt(519)-rrt(574)-rrt(575)-rrt(576)-rrt(577)-rrt(578)&
             -rrt(579)-rrt(580)-rrt(581)-rrt(582)-rrt(640)-rrt(641)-rrt(642)-rrt(643)-rrt(644)-rrt(645)-rrt(646) 
  ydot(52) = -rrt(126)+rrt(420)-rrt(426)-rrt(427)+rrt(437)-rrt(528)-rrt(537)-rrt(546)-rrt(555)-rrt(564)-rrt(573)-rrt(582)-rrt(593) 
  ydot(53) = +rrt(105)+rrt(106)+rrt(107)+rrt(108)+rrt(109)+rrt(110)+rrt(111)+rrt(112)-rrt(113)-rrt(114)-rrt(115)-rrt(116)-rrt(117)&
             -rrt(118)-rrt(119)-rrt(120)-rrt(121)-rrt(122)-rrt(123)-rrt(124)-rrt(125)-rrt(126)-rrt(127)-rrt(128)-rrt(129)-rrt(130)&
             -rrt(131)-rrt(132)-rrt(133)-rrt(134)-rrt(135)-rrt(136)-rrt(137)-rrt(138)-rrt(139)-rrt(140)-rrt(141)-rrt(142)-rrt(143)&
             +rrt(144)+rrt(145)+rrt(146)+rrt(147)+rrt(148)+rrt(149)+rrt(150)+rrt(151)+rrt(152)+rrt(153)+rrt(154)+rrt(155)+rrt(156)&
             +rrt(157)+rrt(158)+rrt(159)+rrt(160)+rrt(161)+rrt(162)+rrt(163)+rrt(164)+rrt(165)+rrt(166)+rrt(167)+rrt(168)+rrt(169)&
             +rrt(170)+rrt(171)+rrt(172)+rrt(173)+rrt(174)+rrt(175)+rrt(176)+rrt(177)+rrt(178)+rrt(179)+rrt(180)+rrt(181)+rrt(213)&
             +rrt(214)+rrt(234)+rrt(307)+rrt(308) 
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
  pd(01,01) = pd(01,01) - rrt(001) * density(53) 
  pd(01,53) = pd(01,53) - rrt(001) * density(01) 
  pd(02,01) = pd(02,01) + rrt(001) * density(53) 
  pd(02,53) = pd(02,53) + rrt(001) * density(01) 
  pd(01,01) = pd(01,01) - rrt(002) * density(53) 
  pd(01,53) = pd(01,53) - rrt(002) * density(01) 
  pd(02,01) = pd(02,01) + rrt(002) * density(53) 
  pd(02,53) = pd(02,53) + rrt(002) * density(01) 
  pd(01,01) = pd(01,01) - rrt(003) * density(53) 
  pd(01,53) = pd(01,53) - rrt(003) * density(01) 
  pd(03,01) = pd(03,01) + rrt(003) * density(53) 
  pd(03,53) = pd(03,53) + rrt(003) * density(01) 
  pd(01,01) = pd(01,01) - rrt(004) * density(53) 
  pd(01,53) = pd(01,53) - rrt(004) * density(01) 
  pd(04,01) = pd(04,01) + rrt(004) * density(53) 
  pd(04,53) = pd(04,53) + rrt(004) * density(01) 
  pd(01,01) = pd(01,01) - rrt(005) * density(53) 
  pd(01,53) = pd(01,53) - rrt(005) * density(01) 
  pd(05,01) = pd(05,01) + rrt(005) * density(53) 
  pd(05,53) = pd(05,53) + rrt(005) * density(01) 
  pd(01,01) = pd(01,01) - rrt(006) * density(53) 
  pd(01,53) = pd(01,53) - rrt(006) * density(01) 
  pd(06,01) = pd(06,01) + rrt(006) * density(53) 
  pd(06,53) = pd(06,53) + rrt(006) * density(01) 
  pd(01,01) = pd(01,01) - rrt(007) * density(53) 
  pd(01,53) = pd(01,53) - rrt(007) * density(01) 
  pd(07,01) = pd(07,01) + rrt(007) * density(53) 
  pd(07,53) = pd(07,53) + rrt(007) * density(01) 
  pd(01,01) = pd(01,01) - rrt(008) * density(53) 
  pd(01,53) = pd(01,53) - rrt(008) * density(01) 
  pd(08,01) = pd(08,01) + rrt(008) * density(53) 
  pd(08,53) = pd(08,53) + rrt(008) * density(01) 
  pd(01,01) = pd(01,01) - rrt(009) * density(53) 
  pd(01,53) = pd(01,53) - rrt(009) * density(01) 
  pd(09,01) = pd(09,01) + rrt(009) * density(53) 
  pd(09,53) = pd(09,53) + rrt(009) * density(01) 
  pd(21,21) = pd(21,21) - rrt(010) * density(53) 
  pd(21,53) = pd(21,53) - rrt(010) * density(21) 
  pd(22,21) = pd(22,21) + rrt(010) * density(53) 
  pd(22,53) = pd(22,53) + rrt(010) * density(21) 
  pd(21,21) = pd(21,21) - rrt(011) * density(53) 
  pd(21,53) = pd(21,53) - rrt(011) * density(21) 
  pd(22,21) = pd(22,21) + rrt(011) * density(53) 
  pd(22,53) = pd(22,53) + rrt(011) * density(21) 
  pd(21,21) = pd(21,21) - rrt(012) * density(53) 
  pd(21,53) = pd(21,53) - rrt(012) * density(21) 
  pd(23,21) = pd(23,21) + rrt(012) * density(53) 
  pd(23,53) = pd(23,53) + rrt(012) * density(21) 
  pd(21,21) = pd(21,21) - rrt(013) * density(53) 
  pd(21,53) = pd(21,53) - rrt(013) * density(21) 
  pd(23,21) = pd(23,21) + rrt(013) * density(53) 
  pd(23,53) = pd(23,53) + rrt(013) * density(21) 
  pd(21,21) = pd(21,21) - rrt(014) * density(53) 
  pd(21,53) = pd(21,53) - rrt(014) * density(21) 
  pd(24,21) = pd(24,21) + rrt(014) * density(53) 
  pd(24,53) = pd(24,53) + rrt(014) * density(21) 
  pd(21,21) = pd(21,21) - rrt(015) * density(53) 
  pd(21,53) = pd(21,53) - rrt(015) * density(21) 
  pd(25,21) = pd(25,21) + rrt(015) * density(53) 
  pd(25,53) = pd(25,53) + rrt(015) * density(21) 
  pd(21,22) = pd(21,22) + rrt(016) * density(53) 
  pd(21,53) = pd(21,53) + rrt(016) * density(22) 
  pd(22,22) = pd(22,22) - rrt(016) * density(53) 
  pd(22,53) = pd(22,53) - rrt(016) * density(22) 
  pd(21,23) = pd(21,23) + rrt(017) * density(53) 
  pd(21,53) = pd(21,53) + rrt(017) * density(23) 
  pd(23,23) = pd(23,23) - rrt(017) * density(53) 
  pd(23,53) = pd(23,53) - rrt(017) * density(23) 
  pd(21,24) = pd(21,24) + rrt(018) * density(53) 
  pd(21,53) = pd(21,53) + rrt(018) * density(24) 
  pd(24,24) = pd(24,24) - rrt(018) * density(53) 
  pd(24,53) = pd(24,53) - rrt(018) * density(24) 
  pd(21,25) = pd(21,25) + rrt(019) * density(53) 
  pd(21,53) = pd(21,53) + rrt(019) * density(25) 
  pd(25,25) = pd(25,25) - rrt(019) * density(53) 
  pd(25,53) = pd(25,53) - rrt(019) * density(25) 
  pd(01,01) = pd(01,01) + rrt(020) * density(02) 
  pd(01,02) = pd(01,02) + rrt(020) * density(01) 
  pd(02,01) = pd(02,01) - rrt(020) * density(02) 
  pd(02,02) = pd(02,02) - rrt(020) * density(01) 
  pd(02,01) = pd(02,01) + rrt(021) * density(03) 
  pd(02,03) = pd(02,03) + rrt(021) * density(01) 
  pd(03,01) = pd(03,01) - rrt(021) * density(03) 
  pd(03,03) = pd(03,03) - rrt(021) * density(01) 
  pd(03,01) = pd(03,01) + rrt(022) * density(04) 
  pd(03,04) = pd(03,04) + rrt(022) * density(01) 
  pd(04,01) = pd(04,01) - rrt(022) * density(04) 
  pd(04,04) = pd(04,04) - rrt(022) * density(01) 
  pd(04,01) = pd(04,01) + rrt(023) * density(05) 
  pd(04,05) = pd(04,05) + rrt(023) * density(01) 
  pd(05,01) = pd(05,01) - rrt(023) * density(05) 
  pd(05,05) = pd(05,05) - rrt(023) * density(01) 
  pd(05,01) = pd(05,01) + rrt(024) * density(06) 
  pd(05,06) = pd(05,06) + rrt(024) * density(01) 
  pd(06,01) = pd(06,01) - rrt(024) * density(06) 
  pd(06,06) = pd(06,06) - rrt(024) * density(01) 
  pd(06,01) = pd(06,01) + rrt(025) * density(07) 
  pd(06,07) = pd(06,07) + rrt(025) * density(01) 
  pd(07,01) = pd(07,01) - rrt(025) * density(07) 
  pd(07,07) = pd(07,07) - rrt(025) * density(01) 
  pd(07,01) = pd(07,01) + rrt(026) * density(08) 
  pd(07,08) = pd(07,08) + rrt(026) * density(01) 
  pd(08,01) = pd(08,01) - rrt(026) * density(08) 
  pd(08,08) = pd(08,08) - rrt(026) * density(01) 
  pd(08,01) = pd(08,01) + rrt(027) * density(09) 
  pd(08,09) = pd(08,09) + rrt(027) * density(01) 
  pd(09,01) = pd(09,01) - rrt(027) * density(09) 
  pd(09,09) = pd(09,09) - rrt(027) * density(01) 
  pd(01,01) = pd(01,01) - rrt(028) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) + rrt(028) * density(01) * 2.0d0
  pd(02,01) = pd(02,01) - rrt(029) * density(02) 
  pd(02,02) = pd(02,02) - rrt(029) * density(01) 
  pd(03,01) = pd(03,01) + rrt(029) * density(02) 
  pd(03,02) = pd(03,02) + rrt(029) * density(01) 
  pd(03,01) = pd(03,01) - rrt(030) * density(03) 
  pd(03,03) = pd(03,03) - rrt(030) * density(01) 
  pd(04,01) = pd(04,01) + rrt(030) * density(03) 
  pd(04,03) = pd(04,03) + rrt(030) * density(01) 
  pd(04,01) = pd(04,01) - rrt(031) * density(04) 
  pd(04,04) = pd(04,04) - rrt(031) * density(01) 
  pd(05,01) = pd(05,01) + rrt(031) * density(04) 
  pd(05,04) = pd(05,04) + rrt(031) * density(01) 
  pd(05,01) = pd(05,01) - rrt(032) * density(05) 
  pd(05,05) = pd(05,05) - rrt(032) * density(01) 
  pd(06,01) = pd(06,01) + rrt(032) * density(05) 
  pd(06,05) = pd(06,05) + rrt(032) * density(01) 
  pd(06,01) = pd(06,01) - rrt(033) * density(06) 
  pd(06,06) = pd(06,06) - rrt(033) * density(01) 
  pd(07,01) = pd(07,01) + rrt(033) * density(06) 
  pd(07,06) = pd(07,06) + rrt(033) * density(01) 
  pd(07,01) = pd(07,01) - rrt(034) * density(07) 
  pd(07,07) = pd(07,07) - rrt(034) * density(01) 
  pd(08,01) = pd(08,01) + rrt(034) * density(07) 
  pd(08,07) = pd(08,07) + rrt(034) * density(01) 
  pd(08,01) = pd(08,01) - rrt(035) * density(08) 
  pd(08,08) = pd(08,08) - rrt(035) * density(01) 
  pd(09,01) = pd(09,01) + rrt(035) * density(08) 
  pd(09,08) = pd(09,08) + rrt(035) * density(01) 
  pd(01,02) = pd(01,02) + rrt(036) * density(14) 
  pd(01,14) = pd(01,14) + rrt(036) * density(02) 
  pd(02,02) = pd(02,02) - rrt(036) * density(14) 
  pd(02,14) = pd(02,14) - rrt(036) * density(02) 
  pd(02,03) = pd(02,03) + rrt(037) * density(14) 
  pd(02,14) = pd(02,14) + rrt(037) * density(03) 
  pd(03,03) = pd(03,03) - rrt(037) * density(14) 
  pd(03,14) = pd(03,14) - rrt(037) * density(03) 
  pd(03,04) = pd(03,04) + rrt(038) * density(14) 
  pd(03,14) = pd(03,14) + rrt(038) * density(04) 
  pd(04,04) = pd(04,04) - rrt(038) * density(14) 
  pd(04,14) = pd(04,14) - rrt(038) * density(04) 
  pd(04,05) = pd(04,05) + rrt(039) * density(14) 
  pd(04,14) = pd(04,14) + rrt(039) * density(05) 
  pd(05,05) = pd(05,05) - rrt(039) * density(14) 
  pd(05,14) = pd(05,14) - rrt(039) * density(05) 
  pd(05,06) = pd(05,06) + rrt(040) * density(14) 
  pd(05,14) = pd(05,14) + rrt(040) * density(06) 
  pd(06,06) = pd(06,06) - rrt(040) * density(14) 
  pd(06,14) = pd(06,14) - rrt(040) * density(06) 
  pd(06,07) = pd(06,07) + rrt(041) * density(14) 
  pd(06,14) = pd(06,14) + rrt(041) * density(07) 
  pd(07,07) = pd(07,07) - rrt(041) * density(14) 
  pd(07,14) = pd(07,14) - rrt(041) * density(07) 
  pd(07,08) = pd(07,08) + rrt(042) * density(14) 
  pd(07,14) = pd(07,14) + rrt(042) * density(08) 
  pd(08,08) = pd(08,08) - rrt(042) * density(14) 
  pd(08,14) = pd(08,14) - rrt(042) * density(08) 
  pd(08,09) = pd(08,09) + rrt(043) * density(14) 
  pd(08,14) = pd(08,14) + rrt(043) * density(09) 
  pd(09,09) = pd(09,09) - rrt(043) * density(14) 
  pd(09,14) = pd(09,14) - rrt(043) * density(09) 
  pd(01,01) = pd(01,01) - rrt(044) * density(14) 
  pd(01,14) = pd(01,14) - rrt(044) * density(01) 
  pd(02,01) = pd(02,01) + rrt(044) * density(14) 
  pd(02,14) = pd(02,14) + rrt(044) * density(01) 
  pd(02,02) = pd(02,02) - rrt(045) * density(14) 
  pd(02,14) = pd(02,14) - rrt(045) * density(02) 
  pd(03,02) = pd(03,02) + rrt(045) * density(14) 
  pd(03,14) = pd(03,14) + rrt(045) * density(02) 
  pd(03,03) = pd(03,03) - rrt(046) * density(14) 
  pd(03,14) = pd(03,14) - rrt(046) * density(03) 
  pd(04,03) = pd(04,03) + rrt(046) * density(14) 
  pd(04,14) = pd(04,14) + rrt(046) * density(03) 
  pd(04,04) = pd(04,04) - rrt(047) * density(14) 
  pd(04,14) = pd(04,14) - rrt(047) * density(04) 
  pd(05,04) = pd(05,04) + rrt(047) * density(14) 
  pd(05,14) = pd(05,14) + rrt(047) * density(04) 
  pd(05,05) = pd(05,05) - rrt(048) * density(14) 
  pd(05,14) = pd(05,14) - rrt(048) * density(05) 
  pd(06,05) = pd(06,05) + rrt(048) * density(14) 
  pd(06,14) = pd(06,14) + rrt(048) * density(05) 
  pd(06,06) = pd(06,06) - rrt(049) * density(14) 
  pd(06,14) = pd(06,14) - rrt(049) * density(06) 
  pd(07,06) = pd(07,06) + rrt(049) * density(14) 
  pd(07,14) = pd(07,14) + rrt(049) * density(06) 
  pd(07,07) = pd(07,07) - rrt(050) * density(14) 
  pd(07,14) = pd(07,14) - rrt(050) * density(07) 
  pd(08,07) = pd(08,07) + rrt(050) * density(14) 
  pd(08,14) = pd(08,14) + rrt(050) * density(07) 
  pd(08,08) = pd(08,08) - rrt(051) * density(14) 
  pd(08,14) = pd(08,14) - rrt(051) * density(08) 
  pd(09,08) = pd(09,08) + rrt(051) * density(14) 
  pd(09,14) = pd(09,14) + rrt(051) * density(08) 
  pd(01,02) = pd(01,02) + rrt(052) * density(29) 
  pd(01,29) = pd(01,29) + rrt(052) * density(02) 
  pd(02,02) = pd(02,02) - rrt(052) * density(29) 
  pd(02,29) = pd(02,29) - rrt(052) * density(02) 
  pd(02,03) = pd(02,03) + rrt(053) * density(29) 
  pd(02,29) = pd(02,29) + rrt(053) * density(03) 
  pd(03,03) = pd(03,03) - rrt(053) * density(29) 
  pd(03,29) = pd(03,29) - rrt(053) * density(03) 
  pd(03,04) = pd(03,04) + rrt(054) * density(29) 
  pd(03,29) = pd(03,29) + rrt(054) * density(04) 
  pd(04,04) = pd(04,04) - rrt(054) * density(29) 
  pd(04,29) = pd(04,29) - rrt(054) * density(04) 
  pd(04,05) = pd(04,05) + rrt(055) * density(29) 
  pd(04,29) = pd(04,29) + rrt(055) * density(05) 
  pd(05,05) = pd(05,05) - rrt(055) * density(29) 
  pd(05,29) = pd(05,29) - rrt(055) * density(05) 
  pd(05,06) = pd(05,06) + rrt(056) * density(29) 
  pd(05,29) = pd(05,29) + rrt(056) * density(06) 
  pd(06,06) = pd(06,06) - rrt(056) * density(29) 
  pd(06,29) = pd(06,29) - rrt(056) * density(06) 
  pd(06,07) = pd(06,07) + rrt(057) * density(29) 
  pd(06,29) = pd(06,29) + rrt(057) * density(07) 
  pd(07,07) = pd(07,07) - rrt(057) * density(29) 
  pd(07,29) = pd(07,29) - rrt(057) * density(07) 
  pd(07,08) = pd(07,08) + rrt(058) * density(29) 
  pd(07,29) = pd(07,29) + rrt(058) * density(08) 
  pd(08,08) = pd(08,08) - rrt(058) * density(29) 
  pd(08,29) = pd(08,29) - rrt(058) * density(08) 
  pd(08,09) = pd(08,09) + rrt(059) * density(29) 
  pd(08,29) = pd(08,29) + rrt(059) * density(09) 
  pd(09,09) = pd(09,09) - rrt(059) * density(29) 
  pd(09,29) = pd(09,29) - rrt(059) * density(09) 
  pd(01,01) = pd(01,01) - rrt(060) * density(29) 
  pd(01,29) = pd(01,29) - rrt(060) * density(01) 
  pd(02,01) = pd(02,01) + rrt(060) * density(29) 
  pd(02,29) = pd(02,29) + rrt(060) * density(01) 
  pd(02,02) = pd(02,02) - rrt(061) * density(29) 
  pd(02,29) = pd(02,29) - rrt(061) * density(02) 
  pd(03,02) = pd(03,02) + rrt(061) * density(29) 
  pd(03,29) = pd(03,29) + rrt(061) * density(02) 
  pd(03,03) = pd(03,03) - rrt(062) * density(29) 
  pd(03,29) = pd(03,29) - rrt(062) * density(03) 
  pd(04,03) = pd(04,03) + rrt(062) * density(29) 
  pd(04,29) = pd(04,29) + rrt(062) * density(03) 
  pd(04,04) = pd(04,04) - rrt(063) * density(29) 
  pd(04,29) = pd(04,29) - rrt(063) * density(04) 
  pd(05,04) = pd(05,04) + rrt(063) * density(29) 
  pd(05,29) = pd(05,29) + rrt(063) * density(04) 
  pd(05,05) = pd(05,05) - rrt(064) * density(29) 
  pd(05,29) = pd(05,29) - rrt(064) * density(05) 
  pd(06,05) = pd(06,05) + rrt(064) * density(29) 
  pd(06,29) = pd(06,29) + rrt(064) * density(05) 
  pd(06,06) = pd(06,06) - rrt(065) * density(29) 
  pd(06,29) = pd(06,29) - rrt(065) * density(06) 
  pd(07,06) = pd(07,06) + rrt(065) * density(29) 
  pd(07,29) = pd(07,29) + rrt(065) * density(06) 
  pd(07,07) = pd(07,07) - rrt(066) * density(29) 
  pd(07,29) = pd(07,29) - rrt(066) * density(07) 
  pd(08,07) = pd(08,07) + rrt(066) * density(29) 
  pd(08,29) = pd(08,29) + rrt(066) * density(07) 
  pd(08,08) = pd(08,08) - rrt(067) * density(29) 
  pd(08,29) = pd(08,29) - rrt(067) * density(08) 
  pd(09,08) = pd(09,08) + rrt(067) * density(29) 
  pd(09,29) = pd(09,29) + rrt(067) * density(08) 
  pd(21,21) = pd(21,21) + rrt(068) * density(22) 
  pd(21,22) = pd(21,22) + rrt(068) * density(21) 
  pd(22,21) = pd(22,21) - rrt(068) * density(22) 
  pd(22,22) = pd(22,22) - rrt(068) * density(21) 
  pd(22,21) = pd(22,21) + rrt(069) * density(23) 
  pd(22,23) = pd(22,23) + rrt(069) * density(21) 
  pd(23,21) = pd(23,21) - rrt(069) * density(23) 
  pd(23,23) = pd(23,23) - rrt(069) * density(21) 
  pd(23,21) = pd(23,21) + rrt(070) * density(24) 
  pd(23,24) = pd(23,24) + rrt(070) * density(21) 
  pd(24,21) = pd(24,21) - rrt(070) * density(24) 
  pd(24,24) = pd(24,24) - rrt(070) * density(21) 
  pd(24,21) = pd(24,21) + rrt(071) * density(25) 
  pd(24,25) = pd(24,25) + rrt(071) * density(21) 
  pd(25,21) = pd(25,21) - rrt(071) * density(25) 
  pd(25,25) = pd(25,25) - rrt(071) * density(21) 
  pd(21,21) = pd(21,21) - rrt(072) * density(21) * 2.0d0
  pd(22,21) = pd(22,21) + rrt(072) * density(21) * 2.0d0
  pd(22,21) = pd(22,21) - rrt(073) * density(22) 
  pd(22,22) = pd(22,22) - rrt(073) * density(21) 
  pd(23,21) = pd(23,21) + rrt(073) * density(22) 
  pd(23,22) = pd(23,22) + rrt(073) * density(21) 
  pd(23,21) = pd(23,21) - rrt(074) * density(23) 
  pd(23,23) = pd(23,23) - rrt(074) * density(21) 
  pd(24,21) = pd(24,21) + rrt(074) * density(23) 
  pd(24,23) = pd(24,23) + rrt(074) * density(21) 
  pd(24,21) = pd(24,21) - rrt(075) * density(24) 
  pd(24,24) = pd(24,24) - rrt(075) * density(21) 
  pd(25,21) = pd(25,21) + rrt(075) * density(24) 
  pd(25,24) = pd(25,24) + rrt(075) * density(21) 
  pd(21,22) = pd(21,22) + rrt(076) * density(29) 
  pd(21,29) = pd(21,29) + rrt(076) * density(22) 
  pd(22,22) = pd(22,22) - rrt(076) * density(29) 
  pd(22,29) = pd(22,29) - rrt(076) * density(22) 
  pd(22,23) = pd(22,23) + rrt(077) * density(29) 
  pd(22,29) = pd(22,29) + rrt(077) * density(23) 
  pd(23,23) = pd(23,23) - rrt(077) * density(29) 
  pd(23,29) = pd(23,29) - rrt(077) * density(23) 
  pd(23,24) = pd(23,24) + rrt(078) * density(29) 
  pd(23,29) = pd(23,29) + rrt(078) * density(24) 
  pd(24,24) = pd(24,24) - rrt(078) * density(29) 
  pd(24,29) = pd(24,29) - rrt(078) * density(24) 
  pd(24,25) = pd(24,25) + rrt(079) * density(29) 
  pd(24,29) = pd(24,29) + rrt(079) * density(25) 
  pd(25,25) = pd(25,25) - rrt(079) * density(29) 
  pd(25,29) = pd(25,29) - rrt(079) * density(25) 
  pd(21,21) = pd(21,21) - rrt(080) * density(29) 
  pd(21,29) = pd(21,29) - rrt(080) * density(21) 
  pd(22,21) = pd(22,21) + rrt(080) * density(29) 
  pd(22,29) = pd(22,29) + rrt(080) * density(21) 
  pd(22,22) = pd(22,22) - rrt(081) * density(29) 
  pd(22,29) = pd(22,29) - rrt(081) * density(22) 
  pd(23,22) = pd(23,22) + rrt(081) * density(29) 
  pd(23,29) = pd(23,29) + rrt(081) * density(22) 
  pd(23,23) = pd(23,23) - rrt(082) * density(29) 
  pd(23,29) = pd(23,29) - rrt(082) * density(23) 
  pd(24,23) = pd(24,23) + rrt(082) * density(29) 
  pd(24,29) = pd(24,29) + rrt(082) * density(23) 
  pd(24,24) = pd(24,24) - rrt(083) * density(29) 
  pd(24,29) = pd(24,29) - rrt(083) * density(24) 
  pd(25,24) = pd(25,24) + rrt(083) * density(29) 
  pd(25,29) = pd(25,29) + rrt(083) * density(24) 
  pd(01,01) = pd(01,01) - rrt(084) * density(53) 
  pd(01,53) = pd(01,53) - rrt(084) * density(01) 
  pd(10,01) = pd(10,01) + rrt(084) * density(53) 
  pd(10,53) = pd(10,53) + rrt(084) * density(01) 
  pd(01,01) = pd(01,01) - rrt(085) * density(53) 
  pd(01,53) = pd(01,53) - rrt(085) * density(01) 
  pd(10,01) = pd(10,01) + rrt(085) * density(53) 
  pd(10,53) = pd(10,53) + rrt(085) * density(01) 
  pd(01,01) = pd(01,01) - rrt(086) * density(53) 
  pd(01,53) = pd(01,53) - rrt(086) * density(01) 
  pd(10,01) = pd(10,01) + rrt(086) * density(53) 
  pd(10,53) = pd(10,53) + rrt(086) * density(01) 
  pd(01,01) = pd(01,01) - rrt(087) * density(53) 
  pd(01,53) = pd(01,53) - rrt(087) * density(01) 
  pd(11,01) = pd(11,01) + rrt(087) * density(53) 
  pd(11,53) = pd(11,53) + rrt(087) * density(01) 
  pd(01,01) = pd(01,01) - rrt(088) * density(53) 
  pd(01,53) = pd(01,53) - rrt(088) * density(01) 
  pd(11,01) = pd(11,01) + rrt(088) * density(53) 
  pd(11,53) = pd(11,53) + rrt(088) * density(01) 
  pd(01,01) = pd(01,01) - rrt(089) * density(53) 
  pd(01,53) = pd(01,53) - rrt(089) * density(01) 
  pd(11,01) = pd(11,01) + rrt(089) * density(53) 
  pd(11,53) = pd(11,53) + rrt(089) * density(01) 
  pd(01,01) = pd(01,01) - rrt(090) * density(53) 
  pd(01,53) = pd(01,53) - rrt(090) * density(01) 
  pd(12,01) = pd(12,01) + rrt(090) * density(53) 
  pd(12,53) = pd(12,53) + rrt(090) * density(01) 
  pd(01,01) = pd(01,01) - rrt(091) * density(53) 
  pd(01,53) = pd(01,53) - rrt(091) * density(01) 
  pd(12,01) = pd(12,01) + rrt(091) * density(53) 
  pd(12,53) = pd(12,53) + rrt(091) * density(01) 
  pd(01,01) = pd(01,01) - rrt(092) * density(53) 
  pd(01,53) = pd(01,53) - rrt(092) * density(01) 
  pd(12,01) = pd(12,01) + rrt(092) * density(53) 
  pd(12,53) = pd(12,53) + rrt(092) * density(01) 
  pd(01,01) = pd(01,01) - rrt(093) * density(53) 
  pd(01,53) = pd(01,53) - rrt(093) * density(01) 
  pd(13,01) = pd(13,01) + rrt(093) * density(53) 
  pd(13,53) = pd(13,53) + rrt(093) * density(01) 
  pd(01,01) = pd(01,01) - rrt(094) * density(53) 
  pd(01,53) = pd(01,53) - rrt(094) * density(01) 
  pd(13,01) = pd(13,01) + rrt(094) * density(53) 
  pd(13,53) = pd(13,53) + rrt(094) * density(01) 
  pd(01,01) = pd(01,01) - rrt(095) * density(53) 
  pd(01,53) = pd(01,53) - rrt(095) * density(01) 
  pd(13,01) = pd(13,01) + rrt(095) * density(53) 
  pd(13,53) = pd(13,53) + rrt(095) * density(01) 
  pd(21,21) = pd(21,21) - rrt(096) * density(53) 
  pd(21,53) = pd(21,53) - rrt(096) * density(21) 
  pd(26,21) = pd(26,21) + rrt(096) * density(53) 
  pd(26,53) = pd(26,53) + rrt(096) * density(21) 
  pd(21,21) = pd(21,21) - rrt(097) * density(53) 
  pd(21,53) = pd(21,53) - rrt(097) * density(21) 
  pd(27,21) = pd(27,21) + rrt(097) * density(53) 
  pd(27,53) = pd(27,53) + rrt(097) * density(21) 
  pd(21,21) = pd(21,21) - rrt(098) * density(53) 
  pd(21,53) = pd(21,53) - rrt(098) * density(21) 
  pd(28,21) = pd(28,21) + rrt(098) * density(53) 
  pd(28,53) = pd(28,53) + rrt(098) * density(21) 
  pd(21,21) = pd(21,21) - rrt(099) * density(53) 
  pd(21,53) = pd(21,53) - rrt(099) * density(21) 
  pd(29,21) = pd(29,21) + rrt(099) * density(53) * 2.0d0
  pd(29,53) = pd(29,53) + rrt(099) * density(21) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(100) * density(53) 
  pd(21,53) = pd(21,53) - rrt(100) * density(21) 
  pd(29,21) = pd(29,21) + rrt(100) * density(53) 
  pd(29,53) = pd(29,53) + rrt(100) * density(21) 
  pd(30,21) = pd(30,21) + rrt(100) * density(53) 
  pd(30,53) = pd(30,53) + rrt(100) * density(21) 
  pd(21,21) = pd(21,21) - rrt(101) * density(53) 
  pd(21,53) = pd(21,53) - rrt(101) * density(21) 
  pd(29,21) = pd(29,21) + rrt(101) * density(53) 
  pd(29,53) = pd(29,53) + rrt(101) * density(21) 
  pd(31,21) = pd(31,21) + rrt(101) * density(53) 
  pd(31,53) = pd(31,53) + rrt(101) * density(21) 
  pd(29,29) = pd(29,29) - rrt(102) * density(53) 
  pd(29,53) = pd(29,53) - rrt(102) * density(29) 
  pd(30,29) = pd(30,29) + rrt(102) * density(53) 
  pd(30,53) = pd(30,53) + rrt(102) * density(29) 
  pd(29,29) = pd(29,29) - rrt(103) * density(53) 
  pd(29,53) = pd(29,53) - rrt(103) * density(29) 
  pd(31,29) = pd(31,29) + rrt(103) * density(53) 
  pd(31,53) = pd(31,53) + rrt(103) * density(29) 
  pd(21,26) = pd(21,26) + rrt(104) * density(53) 
  pd(21,53) = pd(21,53) + rrt(104) * density(26) 
  pd(26,26) = pd(26,26) - rrt(104) * density(53) 
  pd(26,53) = pd(26,53) - rrt(104) * density(26) 
  pd(14,14) = pd(14,14) - rrt(105) * density(53) 
  pd(14,53) = pd(14,53) - rrt(105) * density(14) 
  pd(17,14) = pd(17,14) + rrt(105) * density(53) 
  pd(17,53) = pd(17,53) + rrt(105) * density(14) 
  pd(53,14) = pd(53,14) + rrt(105) * density(53) 
  pd(53,53) = pd(53,53) + rrt(105) * density(14) 
  pd(29,29) = pd(29,29) - rrt(106) * density(53) 
  pd(29,53) = pd(29,53) - rrt(106) * density(29) 
  pd(33,29) = pd(33,29) + rrt(106) * density(53) 
  pd(33,53) = pd(33,53) + rrt(106) * density(29) 
  pd(53,29) = pd(53,29) + rrt(106) * density(53) 
  pd(53,53) = pd(53,53) + rrt(106) * density(29) 
  pd(01,01) = pd(01,01) - rrt(107) * density(53) 
  pd(01,53) = pd(01,53) - rrt(107) * density(01) 
  pd(18,01) = pd(18,01) + rrt(107) * density(53) 
  pd(18,53) = pd(18,53) + rrt(107) * density(01) 
  pd(53,01) = pd(53,01) + rrt(107) * density(53) 
  pd(53,53) = pd(53,53) + rrt(107) * density(01) 
  pd(10,10) = pd(10,10) - rrt(108) * density(53) 
  pd(10,53) = pd(10,53) - rrt(108) * density(10) 
  pd(18,10) = pd(18,10) + rrt(108) * density(53) 
  pd(18,53) = pd(18,53) + rrt(108) * density(10) 
  pd(53,10) = pd(53,10) + rrt(108) * density(53) 
  pd(53,53) = pd(53,53) + rrt(108) * density(10) 
  pd(21,21) = pd(21,21) - rrt(109) * density(53) 
  pd(21,53) = pd(21,53) - rrt(109) * density(21) 
  pd(34,21) = pd(34,21) + rrt(109) * density(53) 
  pd(34,53) = pd(34,53) + rrt(109) * density(21) 
  pd(53,21) = pd(53,21) + rrt(109) * density(53) 
  pd(53,53) = pd(53,53) + rrt(109) * density(21) 
  pd(26,26) = pd(26,26) - rrt(110) * density(53) 
  pd(26,53) = pd(26,53) - rrt(110) * density(26) 
  pd(34,26) = pd(34,26) + rrt(110) * density(53) 
  pd(34,53) = pd(34,53) + rrt(110) * density(26) 
  pd(53,26) = pd(53,26) + rrt(110) * density(53) 
  pd(53,53) = pd(53,53) + rrt(110) * density(26) 
  pd(40,40) = pd(40,40) - rrt(111) * density(53) 
  pd(40,53) = pd(40,53) - rrt(111) * density(40) 
  pd(45,40) = pd(45,40) + rrt(111) * density(53) 
  pd(45,53) = pd(45,53) + rrt(111) * density(40) 
  pd(53,40) = pd(53,40) + rrt(111) * density(53) 
  pd(53,53) = pd(53,53) + rrt(111) * density(40) 
  pd(41,41) = pd(41,41) - rrt(112) * density(53) 
  pd(41,53) = pd(41,53) - rrt(112) * density(41) 
  pd(46,41) = pd(46,41) + rrt(112) * density(53) 
  pd(46,53) = pd(46,53) + rrt(112) * density(41) 
  pd(53,41) = pd(53,41) + rrt(112) * density(53) 
  pd(53,53) = pd(53,53) + rrt(112) * density(41) 
  pd(14,18) = pd(14,18) + rrt(113) * density(53) * 2.0d0
  pd(14,53) = pd(14,53) + rrt(113) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(113) * density(53) 
  pd(18,53) = pd(18,53) - rrt(113) * density(18) 
  pd(53,18) = pd(53,18) - rrt(113) * density(53) 
  pd(53,53) = pd(53,53) - rrt(113) * density(18) 
  pd(14,18) = pd(14,18) + rrt(114) * density(53) 
  pd(14,53) = pd(14,53) + rrt(114) * density(18) 
  pd(15,18) = pd(15,18) + rrt(114) * density(53) 
  pd(15,53) = pd(15,53) + rrt(114) * density(18) 
  pd(18,18) = pd(18,18) - rrt(114) * density(53) 
  pd(18,53) = pd(18,53) - rrt(114) * density(18) 
  pd(53,18) = pd(53,18) - rrt(114) * density(53) 
  pd(53,53) = pd(53,53) - rrt(114) * density(18) 
  pd(14,18) = pd(14,18) + rrt(115) * density(53) 
  pd(14,53) = pd(14,53) + rrt(115) * density(18) 
  pd(16,18) = pd(16,18) + rrt(115) * density(53) 
  pd(16,53) = pd(16,53) + rrt(115) * density(18) 
  pd(18,18) = pd(18,18) - rrt(115) * density(53) 
  pd(18,53) = pd(18,53) - rrt(115) * density(18) 
  pd(53,18) = pd(53,18) - rrt(115) * density(53) 
  pd(53,53) = pd(53,53) - rrt(115) * density(18) 
  pd(29,34) = pd(29,34) + rrt(116) * density(53) * 2.0d0
  pd(29,53) = pd(29,53) + rrt(116) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(116) * density(53) 
  pd(34,53) = pd(34,53) - rrt(116) * density(34) 
  pd(53,34) = pd(53,34) - rrt(116) * density(53) 
  pd(53,53) = pd(53,53) - rrt(116) * density(34) 
  pd(29,34) = pd(29,34) + rrt(117) * density(53) 
  pd(29,53) = pd(29,53) + rrt(117) * density(34) 
  pd(30,34) = pd(30,34) + rrt(117) * density(53) 
  pd(30,53) = pd(30,53) + rrt(117) * density(34) 
  pd(34,34) = pd(34,34) - rrt(117) * density(53) 
  pd(34,53) = pd(34,53) - rrt(117) * density(34) 
  pd(53,34) = pd(53,34) - rrt(117) * density(53) 
  pd(53,53) = pd(53,53) - rrt(117) * density(34) 
  pd(29,34) = pd(29,34) + rrt(118) * density(53) 
  pd(29,53) = pd(29,53) + rrt(118) * density(34) 
  pd(31,34) = pd(31,34) + rrt(118) * density(53) 
  pd(31,53) = pd(31,53) + rrt(118) * density(34) 
  pd(34,34) = pd(34,34) - rrt(118) * density(53) 
  pd(34,53) = pd(34,53) - rrt(118) * density(34) 
  pd(53,34) = pd(53,34) - rrt(118) * density(53) 
  pd(53,53) = pd(53,53) - rrt(118) * density(34) 
  pd(14,45) = pd(14,45) + rrt(119) * density(53) 
  pd(14,53) = pd(14,53) + rrt(119) * density(45) 
  pd(29,45) = pd(29,45) + rrt(119) * density(53) 
  pd(29,53) = pd(29,53) + rrt(119) * density(45) 
  pd(45,45) = pd(45,45) - rrt(119) * density(53) 
  pd(45,53) = pd(45,53) - rrt(119) * density(45) 
  pd(53,45) = pd(53,45) - rrt(119) * density(53) 
  pd(53,53) = pd(53,53) - rrt(119) * density(45) 
  pd(15,45) = pd(15,45) + rrt(120) * density(53) 
  pd(15,53) = pd(15,53) + rrt(120) * density(45) 
  pd(29,45) = pd(29,45) + rrt(120) * density(53) 
  pd(29,53) = pd(29,53) + rrt(120) * density(45) 
  pd(45,45) = pd(45,45) - rrt(120) * density(53) 
  pd(45,53) = pd(45,53) - rrt(120) * density(45) 
  pd(53,45) = pd(53,45) - rrt(120) * density(53) 
  pd(53,53) = pd(53,53) - rrt(120) * density(45) 
  pd(01,19) = pd(01,19) + rrt(121) * density(53) 
  pd(01,53) = pd(01,53) + rrt(121) * density(19) 
  pd(14,19) = pd(14,19) + rrt(121) * density(53) 
  pd(14,53) = pd(14,53) + rrt(121) * density(19) 
  pd(19,19) = pd(19,19) - rrt(121) * density(53) 
  pd(19,53) = pd(19,53) - rrt(121) * density(19) 
  pd(53,19) = pd(53,19) - rrt(121) * density(53) 
  pd(53,53) = pd(53,53) - rrt(121) * density(19) 
  pd(01,20) = pd(01,20) + rrt(122) * density(53) * 2.0d0
  pd(01,53) = pd(01,53) + rrt(122) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(122) * density(53) 
  pd(20,53) = pd(20,53) - rrt(122) * density(20) 
  pd(53,20) = pd(53,20) - rrt(122) * density(53) 
  pd(53,53) = pd(53,53) - rrt(122) * density(20) 
  pd(01,46) = pd(01,46) + rrt(123) * density(53) 
  pd(01,53) = pd(01,53) + rrt(123) * density(46) 
  pd(29,46) = pd(29,46) + rrt(123) * density(53) 
  pd(29,53) = pd(29,53) + rrt(123) * density(46) 
  pd(46,46) = pd(46,46) - rrt(123) * density(53) 
  pd(46,53) = pd(46,53) - rrt(123) * density(46) 
  pd(53,46) = pd(53,46) - rrt(123) * density(53) 
  pd(53,53) = pd(53,53) - rrt(123) * density(46) 
  pd(29,47) = pd(29,47) + rrt(124) * density(53) 
  pd(29,53) = pd(29,53) + rrt(124) * density(47) 
  pd(40,47) = pd(40,47) + rrt(124) * density(53) 
  pd(40,53) = pd(40,53) + rrt(124) * density(47) 
  pd(47,47) = pd(47,47) - rrt(124) * density(53) 
  pd(47,53) = pd(47,53) - rrt(124) * density(47) 
  pd(53,47) = pd(53,47) - rrt(124) * density(53) 
  pd(53,53) = pd(53,53) - rrt(124) * density(47) 
  pd(21,35) = pd(21,35) + rrt(125) * density(53) * 2.0d0
  pd(21,53) = pd(21,53) + rrt(125) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(125) * density(53) 
  pd(35,53) = pd(35,53) - rrt(125) * density(35) 
  pd(53,35) = pd(53,35) - rrt(125) * density(53) 
  pd(53,53) = pd(53,53) - rrt(125) * density(35) 
  pd(01,52) = pd(01,52) + rrt(126) * density(53) 
  pd(01,53) = pd(01,53) + rrt(126) * density(52) 
  pd(21,52) = pd(21,52) + rrt(126) * density(53) 
  pd(21,53) = pd(21,53) + rrt(126) * density(52) 
  pd(52,52) = pd(52,52) - rrt(126) * density(53) 
  pd(52,53) = pd(52,53) - rrt(126) * density(52) 
  pd(53,52) = pd(53,52) - rrt(126) * density(53) 
  pd(53,53) = pd(53,53) - rrt(126) * density(52) 
  pd(14,17) = pd(14,17) + rrt(127) * density(53)**2 
  pd(14,53) = pd(14,53) + rrt(127) * density(17) * density(53) * 2.0d0
  pd(17,17) = pd(17,17) - rrt(127) * density(53)**2 
  pd(17,53) = pd(17,53) - rrt(127) * density(17) * density(53) * 2.0d0
  pd(53,17) = pd(53,17) - rrt(127) * density(53)**2 
  pd(53,53) = pd(53,53) - rrt(127) * density(17) * density(53) * 2.0d0
  pd(29,33) = pd(29,33) + rrt(128) * density(53)**2 
  pd(29,53) = pd(29,53) + rrt(128) * density(33) * density(53) * 2.0d0
  pd(33,33) = pd(33,33) - rrt(128) * density(53)**2 
  pd(33,53) = pd(33,53) - rrt(128) * density(33) * density(53) * 2.0d0
  pd(53,33) = pd(53,33) - rrt(128) * density(53)**2 
  pd(53,53) = pd(53,53) - rrt(128) * density(33) * density(53) * 2.0d0
  pd(14,17) = pd(14,17) + rrt(129) * density(53) 
  pd(14,53) = pd(14,53) + rrt(129) * density(17) 
  pd(17,17) = pd(17,17) - rrt(129) * density(53) 
  pd(17,53) = pd(17,53) - rrt(129) * density(17) 
  pd(53,17) = pd(53,17) - rrt(129) * density(53) 
  pd(53,53) = pd(53,53) - rrt(129) * density(17) 
  pd(29,33) = pd(29,33) + rrt(130) * density(53) 
  pd(29,53) = pd(29,53) + rrt(130) * density(33) 
  pd(33,33) = pd(33,33) - rrt(130) * density(53) 
  pd(33,53) = pd(33,53) - rrt(130) * density(33) 
  pd(53,33) = pd(53,33) - rrt(130) * density(53) 
  pd(53,53) = pd(53,53) - rrt(130) * density(33) 
  pd(21,21) = pd(21,21) - rrt(131) * density(53) 
  pd(21,53) = pd(21,53) - rrt(131) * density(21) 
  pd(29,21) = pd(29,21) + rrt(131) * density(53) 
  pd(29,53) = pd(29,53) + rrt(131) * density(21) 
  pd(36,21) = pd(36,21) + rrt(131) * density(53) 
  pd(36,53) = pd(36,53) + rrt(131) * density(21) 
  pd(53,21) = pd(53,21) - rrt(131) * density(53) 
  pd(53,53) = pd(53,53) - rrt(131) * density(21) 
  pd(14,40) = pd(14,40) + rrt(132) * density(53) 
  pd(14,53) = pd(14,53) + rrt(132) * density(40) 
  pd(36,40) = pd(36,40) + rrt(132) * density(53) 
  pd(36,53) = pd(36,53) + rrt(132) * density(40) 
  pd(40,40) = pd(40,40) - rrt(132) * density(53) 
  pd(40,53) = pd(40,53) - rrt(132) * density(40) 
  pd(53,40) = pd(53,40) - rrt(132) * density(53) 
  pd(53,53) = pd(53,53) - rrt(132) * density(40) 
  pd(21,32) = pd(21,32) + rrt(133) * density(53) 
  pd(21,53) = pd(21,53) + rrt(133) * density(32) 
  pd(32,32) = pd(32,32) - rrt(133) * density(53) 
  pd(32,53) = pd(32,53) - rrt(133) * density(32) 
  pd(36,32) = pd(36,32) + rrt(133) * density(53) 
  pd(36,53) = pd(36,53) + rrt(133) * density(32) 
  pd(53,32) = pd(53,32) - rrt(133) * density(53) 
  pd(53,53) = pd(53,53) - rrt(133) * density(32) 
  pd(29,32) = pd(29,32) + rrt(134) * density(53) 
  pd(29,53) = pd(29,53) + rrt(134) * density(32) 
  pd(32,32) = pd(32,32) - rrt(134) * density(53) 
  pd(32,53) = pd(32,53) - rrt(134) * density(32) 
  pd(37,32) = pd(37,32) + rrt(134) * density(53) 
  pd(37,53) = pd(37,53) + rrt(134) * density(32) 
  pd(53,32) = pd(53,32) - rrt(134) * density(53) 
  pd(53,53) = pd(53,53) - rrt(134) * density(32) 
  pd(14,41) = pd(14,41) + rrt(135) * density(53) 
  pd(14,53) = pd(14,53) + rrt(135) * density(41) 
  pd(41,41) = pd(41,41) - rrt(135) * density(53) 
  pd(41,53) = pd(41,53) - rrt(135) * density(41) 
  pd(48,41) = pd(48,41) + rrt(135) * density(53) 
  pd(48,53) = pd(48,53) + rrt(135) * density(41) 
  pd(53,41) = pd(53,41) - rrt(135) * density(53) 
  pd(53,53) = pd(53,53) - rrt(135) * density(41) 
  pd(21,21) = pd(21,21) - rrt(136) * density(21) * density(53) * 2.0d0
  pd(21,53) = pd(21,53) - rrt(136) * density(21)**2 
  pd(37,21) = pd(37,21) + rrt(136) * density(21) * density(53) * 2.0d0
  pd(37,53) = pd(37,53) + rrt(136) * density(21)**2 
  pd(53,21) = pd(53,21) - rrt(136) * density(21) * density(53) * 2.0d0
  pd(53,53) = pd(53,53) - rrt(136) * density(21)**2 
  pd(36,42) = pd(36,42) + rrt(137) * density(53) 
  pd(36,53) = pd(36,53) + rrt(137) * density(42) 
  pd(40,42) = pd(40,42) + rrt(137) * density(53) 
  pd(40,53) = pd(40,53) + rrt(137) * density(42) 
  pd(42,42) = pd(42,42) - rrt(137) * density(53) 
  pd(42,53) = pd(42,53) - rrt(137) * density(42) 
  pd(53,42) = pd(53,42) - rrt(137) * density(53) 
  pd(53,53) = pd(53,53) - rrt(137) * density(42) 
  pd(29,21) = pd(29,21) - rrt(138) * density(29) * density(53) 
  pd(29,29) = pd(29,29) - rrt(138) * density(21) * density(53) 
  pd(29,53) = pd(29,53) - rrt(138) * density(21) * density(29) 
  pd(36,21) = pd(36,21) + rrt(138) * density(29) * density(53) 
  pd(36,29) = pd(36,29) + rrt(138) * density(21) * density(53) 
  pd(36,53) = pd(36,53) + rrt(138) * density(21) * density(29) 
  pd(53,21) = pd(53,21) - rrt(138) * density(29) * density(53) 
  pd(53,29) = pd(53,29) - rrt(138) * density(21) * density(53) 
  pd(53,53) = pd(53,53) - rrt(138) * density(21) * density(29) 
  pd(21,21) = pd(21,21) - rrt(139) * density(29) * density(53) 
  pd(21,29) = pd(21,29) - rrt(139) * density(21) * density(53) 
  pd(21,53) = pd(21,53) - rrt(139) * density(21) * density(29) 
  pd(37,21) = pd(37,21) + rrt(139) * density(29) * density(53) 
  pd(37,29) = pd(37,29) + rrt(139) * density(21) * density(53) 
  pd(37,53) = pd(37,53) + rrt(139) * density(21) * density(29) 
  pd(53,21) = pd(53,21) - rrt(139) * density(29) * density(53) 
  pd(53,29) = pd(53,29) - rrt(139) * density(21) * density(53) 
  pd(53,53) = pd(53,53) - rrt(139) * density(21) * density(29) 
  pd(32,32) = pd(32,32) - rrt(140) * density(53) 
  pd(32,53) = pd(32,53) - rrt(140) * density(32) 
  pd(38,32) = pd(38,32) + rrt(140) * density(53) 
  pd(38,53) = pd(38,53) + rrt(140) * density(32) 
  pd(53,32) = pd(53,32) - rrt(140) * density(53) 
  pd(53,53) = pd(53,53) - rrt(140) * density(32) 
  pd(40,40) = pd(40,40) - rrt(141) * density(53) 
  pd(40,53) = pd(40,53) - rrt(141) * density(40) 
  pd(48,40) = pd(48,40) + rrt(141) * density(53) 
  pd(48,53) = pd(48,53) + rrt(141) * density(40) 
  pd(53,40) = pd(53,40) - rrt(141) * density(53) 
  pd(53,53) = pd(53,53) - rrt(141) * density(40) 
  pd(41,41) = pd(41,41) - rrt(142) * density(53) 
  pd(41,53) = pd(41,53) - rrt(142) * density(41) 
  pd(49,41) = pd(49,41) + rrt(142) * density(53) 
  pd(49,53) = pd(49,53) + rrt(142) * density(41) 
  pd(53,41) = pd(53,41) - rrt(142) * density(53) 
  pd(53,53) = pd(53,53) - rrt(142) * density(41) 
  pd(21,01) = pd(21,01) - rrt(143) * density(21) * density(53) 
  pd(21,21) = pd(21,21) - rrt(143) * density(01) * density(53) 
  pd(21,53) = pd(21,53) - rrt(143) * density(01) * density(21) 
  pd(37,01) = pd(37,01) + rrt(143) * density(21) * density(53) 
  pd(37,21) = pd(37,21) + rrt(143) * density(01) * density(53) 
  pd(37,53) = pd(37,53) + rrt(143) * density(01) * density(21) 
  pd(53,01) = pd(53,01) - rrt(143) * density(21) * density(53) 
  pd(53,21) = pd(53,21) - rrt(143) * density(01) * density(53) 
  pd(53,53) = pd(53,53) - rrt(143) * density(01) * density(21) 
  pd(21,29) = pd(21,29) + rrt(144) * density(36) 
  pd(21,36) = pd(21,36) + rrt(144) * density(29) 
  pd(29,29) = pd(29,29) - rrt(144) * density(36) 
  pd(29,36) = pd(29,36) - rrt(144) * density(29) 
  pd(36,29) = pd(36,29) - rrt(144) * density(36) 
  pd(36,36) = pd(36,36) - rrt(144) * density(29) 
  pd(53,29) = pd(53,29) + rrt(144) * density(36) 
  pd(53,36) = pd(53,36) + rrt(144) * density(29) 
  pd(14,14) = pd(14,14) - rrt(145) * density(36) 
  pd(14,36) = pd(14,36) - rrt(145) * density(14) 
  pd(36,14) = pd(36,14) - rrt(145) * density(36) 
  pd(36,36) = pd(36,36) - rrt(145) * density(14) 
  pd(40,14) = pd(40,14) + rrt(145) * density(36) 
  pd(40,36) = pd(40,36) + rrt(145) * density(14) 
  pd(53,14) = pd(53,14) + rrt(145) * density(36) 
  pd(53,36) = pd(53,36) + rrt(145) * density(14) 
  pd(36,36) = pd(36,36) - rrt(146) * density(40) 
  pd(36,40) = pd(36,40) - rrt(146) * density(36) 
  pd(40,36) = pd(40,36) - rrt(146) * density(40) 
  pd(40,40) = pd(40,40) - rrt(146) * density(36) 
  pd(42,36) = pd(42,36) + rrt(146) * density(40) 
  pd(42,40) = pd(42,40) + rrt(146) * density(36) 
  pd(53,36) = pd(53,36) + rrt(146) * density(40) 
  pd(53,40) = pd(53,40) + rrt(146) * density(36) 
  pd(01,01) = pd(01,01) - rrt(147) * density(36) 
  pd(01,36) = pd(01,36) - rrt(147) * density(01) 
  pd(36,01) = pd(36,01) - rrt(147) * density(36) 
  pd(36,36) = pd(36,36) - rrt(147) * density(01) 
  pd(41,01) = pd(41,01) + rrt(147) * density(36) 
  pd(41,36) = pd(41,36) + rrt(147) * density(01) 
  pd(53,01) = pd(53,01) + rrt(147) * density(36) 
  pd(53,36) = pd(53,36) + rrt(147) * density(01) 
  pd(21,21) = pd(21,21) - rrt(148) * density(36) 
  pd(21,36) = pd(21,36) - rrt(148) * density(21) 
  pd(32,21) = pd(32,21) + rrt(148) * density(36) 
  pd(32,36) = pd(32,36) + rrt(148) * density(21) 
  pd(36,21) = pd(36,21) - rrt(148) * density(36) 
  pd(36,36) = pd(36,36) - rrt(148) * density(21) 
  pd(53,21) = pd(53,21) + rrt(148) * density(36) 
  pd(53,36) = pd(53,36) + rrt(148) * density(21) 
  pd(26,26) = pd(26,26) - rrt(149) * density(36) 
  pd(26,36) = pd(26,36) - rrt(149) * density(26) 
  pd(32,26) = pd(32,26) + rrt(149) * density(36) 
  pd(32,36) = pd(32,36) + rrt(149) * density(26) 
  pd(36,26) = pd(36,26) - rrt(149) * density(36) 
  pd(36,36) = pd(36,36) - rrt(149) * density(26) 
  pd(53,26) = pd(53,26) + rrt(149) * density(36) 
  pd(53,36) = pd(53,36) + rrt(149) * density(26) 
  pd(21,27) = pd(21,27) + rrt(150) * density(36) 
  pd(21,36) = pd(21,36) + rrt(150) * density(27) 
  pd(27,27) = pd(27,27) - rrt(150) * density(36) 
  pd(27,36) = pd(27,36) - rrt(150) * density(27) 
  pd(29,27) = pd(29,27) + rrt(150) * density(36) 
  pd(29,36) = pd(29,36) + rrt(150) * density(27) 
  pd(36,27) = pd(36,27) - rrt(150) * density(36) 
  pd(36,36) = pd(36,36) - rrt(150) * density(27) 
  pd(53,27) = pd(53,27) + rrt(150) * density(36) 
  pd(53,36) = pd(53,36) + rrt(150) * density(27) 
  pd(01,10) = pd(01,10) + rrt(151) * density(36) 
  pd(01,36) = pd(01,36) + rrt(151) * density(10) 
  pd(10,10) = pd(10,10) - rrt(151) * density(36) 
  pd(10,36) = pd(10,36) - rrt(151) * density(10) 
  pd(29,10) = pd(29,10) + rrt(151) * density(36) 
  pd(29,36) = pd(29,36) + rrt(151) * density(10) 
  pd(36,10) = pd(36,10) - rrt(151) * density(36) 
  pd(36,36) = pd(36,36) - rrt(151) * density(10) 
  pd(53,10) = pd(53,10) + rrt(151) * density(36) 
  pd(53,36) = pd(53,36) + rrt(151) * density(10) 
  pd(01,11) = pd(01,11) + rrt(152) * density(36) 
  pd(01,36) = pd(01,36) + rrt(152) * density(11) 
  pd(11,11) = pd(11,11) - rrt(152) * density(36) 
  pd(11,36) = pd(11,36) - rrt(152) * density(11) 
  pd(29,11) = pd(29,11) + rrt(152) * density(36) 
  pd(29,36) = pd(29,36) + rrt(152) * density(11) 
  pd(36,11) = pd(36,11) - rrt(152) * density(36) 
  pd(36,36) = pd(36,36) - rrt(152) * density(11) 
  pd(53,11) = pd(53,11) + rrt(152) * density(36) 
  pd(53,36) = pd(53,36) + rrt(152) * density(11) 
  pd(21,32) = pd(21,32) + rrt(153) * density(36) * 2.0d0
  pd(21,36) = pd(21,36) + rrt(153) * density(32) * 2.0d0
  pd(32,32) = pd(32,32) - rrt(153) * density(36) 
  pd(32,36) = pd(32,36) - rrt(153) * density(32) 
  pd(36,32) = pd(36,32) - rrt(153) * density(36) 
  pd(36,36) = pd(36,36) - rrt(153) * density(32) 
  pd(53,32) = pd(53,32) + rrt(153) * density(36) 
  pd(53,36) = pd(53,36) + rrt(153) * density(32) 
  pd(29,29) = pd(29,29) - rrt(154) * density(37) 
  pd(29,37) = pd(29,37) - rrt(154) * density(29) 
  pd(32,29) = pd(32,29) + rrt(154) * density(37) 
  pd(32,37) = pd(32,37) + rrt(154) * density(29) 
  pd(37,29) = pd(37,29) - rrt(154) * density(37) 
  pd(37,37) = pd(37,37) - rrt(154) * density(29) 
  pd(53,29) = pd(53,29) + rrt(154) * density(37) 
  pd(53,37) = pd(53,37) + rrt(154) * density(29) 
  pd(14,14) = pd(14,14) - rrt(155) * density(37) 
  pd(14,37) = pd(14,37) - rrt(155) * density(14) 
  pd(37,14) = pd(37,14) - rrt(155) * density(37) 
  pd(37,37) = pd(37,37) - rrt(155) * density(14) 
  pd(42,14) = pd(42,14) + rrt(155) * density(37) 
  pd(42,37) = pd(42,37) + rrt(155) * density(14) 
  pd(53,14) = pd(53,14) + rrt(155) * density(37) 
  pd(53,37) = pd(53,37) + rrt(155) * density(14) 
  pd(21,21) = pd(21,21) + rrt(156) * density(37) 
  pd(21,37) = pd(21,37) + rrt(156) * density(21) 
  pd(37,21) = pd(37,21) - rrt(156) * density(37) 
  pd(37,37) = pd(37,37) - rrt(156) * density(21) 
  pd(53,21) = pd(53,21) + rrt(156) * density(37) 
  pd(53,37) = pd(53,37) + rrt(156) * density(21) 
  pd(21,26) = pd(21,26) + rrt(157) * density(37) * 2.0d0
  pd(21,37) = pd(21,37) + rrt(157) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(157) * density(37) 
  pd(26,37) = pd(26,37) - rrt(157) * density(26) 
  pd(37,26) = pd(37,26) - rrt(157) * density(37) 
  pd(37,37) = pd(37,37) - rrt(157) * density(26) 
  pd(53,26) = pd(53,26) + rrt(157) * density(37) 
  pd(53,37) = pd(53,37) + rrt(157) * density(26) 
  pd(21,27) = pd(21,27) + rrt(158) * density(37) * 2.0d0
  pd(21,37) = pd(21,37) + rrt(158) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(158) * density(37) 
  pd(27,37) = pd(27,37) - rrt(158) * density(27) 
  pd(37,27) = pd(37,27) - rrt(158) * density(37) 
  pd(37,37) = pd(37,37) - rrt(158) * density(27) 
  pd(53,27) = pd(53,27) + rrt(158) * density(37) 
  pd(53,37) = pd(53,37) + rrt(158) * density(27) 
  pd(21,01) = pd(21,01) + rrt(159) * density(37) 
  pd(21,37) = pd(21,37) + rrt(159) * density(01) 
  pd(37,01) = pd(37,01) - rrt(159) * density(37) 
  pd(37,37) = pd(37,37) - rrt(159) * density(01) 
  pd(53,01) = pd(53,01) + rrt(159) * density(37) 
  pd(53,37) = pd(53,37) + rrt(159) * density(01) 
  pd(01,10) = pd(01,10) + rrt(160) * density(37) 
  pd(01,37) = pd(01,37) + rrt(160) * density(10) 
  pd(10,10) = pd(10,10) - rrt(160) * density(37) 
  pd(10,37) = pd(10,37) - rrt(160) * density(10) 
  pd(21,10) = pd(21,10) + rrt(160) * density(37) 
  pd(21,37) = pd(21,37) + rrt(160) * density(10) 
  pd(37,10) = pd(37,10) - rrt(160) * density(37) 
  pd(37,37) = pd(37,37) - rrt(160) * density(10) 
  pd(53,10) = pd(53,10) + rrt(160) * density(37) 
  pd(53,37) = pd(53,37) + rrt(160) * density(10) 
  pd(01,11) = pd(01,11) + rrt(161) * density(37) 
  pd(01,37) = pd(01,37) + rrt(161) * density(11) 
  pd(11,11) = pd(11,11) - rrt(161) * density(37) 
  pd(11,37) = pd(11,37) - rrt(161) * density(11) 
  pd(21,11) = pd(21,11) + rrt(161) * density(37) 
  pd(21,37) = pd(21,37) + rrt(161) * density(11) 
  pd(37,11) = pd(37,11) - rrt(161) * density(37) 
  pd(37,37) = pd(37,37) - rrt(161) * density(11) 
  pd(53,11) = pd(53,11) + rrt(161) * density(37) 
  pd(53,37) = pd(53,37) + rrt(161) * density(11) 
  pd(21,29) = pd(21,29) + rrt(162) * density(38) * 2.0d0
  pd(21,38) = pd(21,38) + rrt(162) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(162) * density(38) 
  pd(29,38) = pd(29,38) - rrt(162) * density(29) 
  pd(38,29) = pd(38,29) - rrt(162) * density(38) 
  pd(38,38) = pd(38,38) - rrt(162) * density(29) 
  pd(53,29) = pd(53,29) + rrt(162) * density(38) 
  pd(53,38) = pd(53,38) + rrt(162) * density(29) 
  pd(14,14) = pd(14,14) - rrt(163) * density(48) 
  pd(14,48) = pd(14,48) - rrt(163) * density(14) 
  pd(41,14) = pd(41,14) + rrt(163) * density(48) 
  pd(41,48) = pd(41,48) + rrt(163) * density(14) 
  pd(48,14) = pd(48,14) - rrt(163) * density(48) 
  pd(48,48) = pd(48,48) - rrt(163) * density(14) 
  pd(53,14) = pd(53,14) + rrt(163) * density(48) 
  pd(53,48) = pd(53,48) + rrt(163) * density(14) 
  pd(14,14) = pd(14,14) - rrt(164) * density(38) 
  pd(14,38) = pd(14,38) - rrt(164) * density(14) 
  pd(21,14) = pd(21,14) + rrt(164) * density(38) 
  pd(21,38) = pd(21,38) + rrt(164) * density(14) 
  pd(38,14) = pd(38,14) - rrt(164) * density(38) 
  pd(38,38) = pd(38,38) - rrt(164) * density(14) 
  pd(40,14) = pd(40,14) + rrt(164) * density(38) 
  pd(40,38) = pd(40,38) + rrt(164) * density(14) 
  pd(53,14) = pd(53,14) + rrt(164) * density(38) 
  pd(53,38) = pd(53,38) + rrt(164) * density(14) 
  pd(01,14) = pd(01,14) + rrt(165) * density(49) 
  pd(01,49) = pd(01,49) + rrt(165) * density(14) 
  pd(14,14) = pd(14,14) - rrt(165) * density(49) 
  pd(14,49) = pd(14,49) - rrt(165) * density(14) 
  pd(40,14) = pd(40,14) + rrt(165) * density(49) 
  pd(40,49) = pd(40,49) + rrt(165) * density(14) 
  pd(49,14) = pd(49,14) - rrt(165) * density(49) 
  pd(49,49) = pd(49,49) - rrt(165) * density(14) 
  pd(53,14) = pd(53,14) + rrt(165) * density(49) 
  pd(53,49) = pd(53,49) + rrt(165) * density(14) 
  pd(14,14) = pd(14,14) - rrt(166) * density(50) 
  pd(14,50) = pd(14,50) - rrt(166) * density(14) 
  pd(40,14) = pd(40,14) + rrt(166) * density(50) * 2.0d0
  pd(40,50) = pd(40,50) + rrt(166) * density(14) * 2.0d0
  pd(50,14) = pd(50,14) - rrt(166) * density(50) 
  pd(50,50) = pd(50,50) - rrt(166) * density(14) 
  pd(53,14) = pd(53,14) + rrt(166) * density(50) 
  pd(53,50) = pd(53,50) + rrt(166) * density(14) 
  pd(14,14) = pd(14,14) - rrt(167) * density(51) 
  pd(14,51) = pd(14,51) - rrt(167) * density(14) 
  pd(40,14) = pd(40,14) + rrt(167) * density(51) 
  pd(40,51) = pd(40,51) + rrt(167) * density(14) 
  pd(42,14) = pd(42,14) + rrt(167) * density(51) 
  pd(42,51) = pd(42,51) + rrt(167) * density(14) 
  pd(51,14) = pd(51,14) - rrt(167) * density(51) 
  pd(51,51) = pd(51,51) - rrt(167) * density(14) 
  pd(53,14) = pd(53,14) + rrt(167) * density(51) 
  pd(53,51) = pd(53,51) + rrt(167) * density(14) 
  pd(29,29) = pd(29,29) - rrt(168) * density(48) 
  pd(29,48) = pd(29,48) - rrt(168) * density(29) 
  pd(42,29) = pd(42,29) + rrt(168) * density(48) 
  pd(42,48) = pd(42,48) + rrt(168) * density(29) 
  pd(48,29) = pd(48,29) - rrt(168) * density(48) 
  pd(48,48) = pd(48,48) - rrt(168) * density(29) 
  pd(53,29) = pd(53,29) + rrt(168) * density(48) 
  pd(53,48) = pd(53,48) + rrt(168) * density(29) 
  pd(29,29) = pd(29,29) - rrt(169) * density(49) 
  pd(29,49) = pd(29,49) - rrt(169) * density(29) 
  pd(40,29) = pd(40,29) + rrt(169) * density(49) * 2.0d0
  pd(40,49) = pd(40,49) + rrt(169) * density(29) * 2.0d0
  pd(49,29) = pd(49,29) - rrt(169) * density(49) 
  pd(49,49) = pd(49,49) - rrt(169) * density(29) 
  pd(53,29) = pd(53,29) + rrt(169) * density(49) 
  pd(53,49) = pd(53,49) + rrt(169) * density(29) 
  pd(21,29) = pd(21,29) + rrt(170) * density(50) 
  pd(21,50) = pd(21,50) + rrt(170) * density(29) 
  pd(29,29) = pd(29,29) - rrt(170) * density(50) 
  pd(29,50) = pd(29,50) - rrt(170) * density(29) 
  pd(40,29) = pd(40,29) + rrt(170) * density(50) 
  pd(40,50) = pd(40,50) + rrt(170) * density(29) 
  pd(50,29) = pd(50,29) - rrt(170) * density(50) 
  pd(50,50) = pd(50,50) - rrt(170) * density(29) 
  pd(53,29) = pd(53,29) + rrt(170) * density(50) 
  pd(53,50) = pd(53,50) + rrt(170) * density(29) 
  pd(29,29) = pd(29,29) - rrt(171) * density(51) 
  pd(29,51) = pd(29,51) - rrt(171) * density(29) 
  pd(32,29) = pd(32,29) + rrt(171) * density(51) 
  pd(32,51) = pd(32,51) + rrt(171) * density(29) 
  pd(40,29) = pd(40,29) + rrt(171) * density(51) 
  pd(40,51) = pd(40,51) + rrt(171) * density(29) 
  pd(51,29) = pd(51,29) - rrt(171) * density(51) 
  pd(51,51) = pd(51,51) - rrt(171) * density(29) 
  pd(53,29) = pd(53,29) + rrt(171) * density(51) 
  pd(53,51) = pd(53,51) + rrt(171) * density(29) 
  pd(01,10) = pd(01,10) + rrt(172) * density(38) 
  pd(01,38) = pd(01,38) + rrt(172) * density(10) 
  pd(10,10) = pd(10,10) - rrt(172) * density(38) 
  pd(10,38) = pd(10,38) - rrt(172) * density(10) 
  pd(32,10) = pd(32,10) + rrt(172) * density(38) 
  pd(32,38) = pd(32,38) + rrt(172) * density(10) 
  pd(38,10) = pd(38,10) - rrt(172) * density(38) 
  pd(38,38) = pd(38,38) - rrt(172) * density(10) 
  pd(53,10) = pd(53,10) + rrt(172) * density(38) 
  pd(53,38) = pd(53,38) + rrt(172) * density(10) 
  pd(01,10) = pd(01,10) + rrt(173) * density(48) 
  pd(01,48) = pd(01,48) + rrt(173) * density(10) 
  pd(10,10) = pd(10,10) - rrt(173) * density(48) 
  pd(10,48) = pd(10,48) - rrt(173) * density(10) 
  pd(40,10) = pd(40,10) + rrt(173) * density(48) 
  pd(40,48) = pd(40,48) + rrt(173) * density(10) 
  pd(48,10) = pd(48,10) - rrt(173) * density(48) 
  pd(48,48) = pd(48,48) - rrt(173) * density(10) 
  pd(53,10) = pd(53,10) + rrt(173) * density(48) 
  pd(53,48) = pd(53,48) + rrt(173) * density(10) 
  pd(01,10) = pd(01,10) + rrt(174) * density(49) 
  pd(01,49) = pd(01,49) + rrt(174) * density(10) 
  pd(10,10) = pd(10,10) - rrt(174) * density(49) 
  pd(10,49) = pd(10,49) - rrt(174) * density(10) 
  pd(41,10) = pd(41,10) + rrt(174) * density(49) 
  pd(41,49) = pd(41,49) + rrt(174) * density(10) 
  pd(49,10) = pd(49,10) - rrt(174) * density(49) 
  pd(49,49) = pd(49,49) - rrt(174) * density(10) 
  pd(53,10) = pd(53,10) + rrt(174) * density(49) 
  pd(53,49) = pd(53,49) + rrt(174) * density(10) 
  pd(01,10) = pd(01,10) + rrt(175) * density(50) 
  pd(01,50) = pd(01,50) + rrt(175) * density(10) 
  pd(10,10) = pd(10,10) - rrt(175) * density(50) 
  pd(10,50) = pd(10,50) - rrt(175) * density(10) 
  pd(42,10) = pd(42,10) + rrt(175) * density(50) 
  pd(42,50) = pd(42,50) + rrt(175) * density(10) 
  pd(50,10) = pd(50,10) - rrt(175) * density(50) 
  pd(50,50) = pd(50,50) - rrt(175) * density(10) 
  pd(53,10) = pd(53,10) + rrt(175) * density(50) 
  pd(53,50) = pd(53,50) + rrt(175) * density(10) 
  pd(01,10) = pd(01,10) + rrt(176) * density(51) 
  pd(01,51) = pd(01,51) + rrt(176) * density(10) 
  pd(10,10) = pd(10,10) - rrt(176) * density(51) 
  pd(10,51) = pd(10,51) - rrt(176) * density(10) 
  pd(43,10) = pd(43,10) + rrt(176) * density(51) 
  pd(43,51) = pd(43,51) + rrt(176) * density(10) 
  pd(51,10) = pd(51,10) - rrt(176) * density(51) 
  pd(51,51) = pd(51,51) - rrt(176) * density(10) 
  pd(53,10) = pd(53,10) + rrt(176) * density(51) 
  pd(53,51) = pd(53,51) + rrt(176) * density(10) 
  pd(01,11) = pd(01,11) + rrt(177) * density(38) 
  pd(01,38) = pd(01,38) + rrt(177) * density(11) 
  pd(11,11) = pd(11,11) - rrt(177) * density(38) 
  pd(11,38) = pd(11,38) - rrt(177) * density(11) 
  pd(32,11) = pd(32,11) + rrt(177) * density(38) 
  pd(32,38) = pd(32,38) + rrt(177) * density(11) 
  pd(38,11) = pd(38,11) - rrt(177) * density(38) 
  pd(38,38) = pd(38,38) - rrt(177) * density(11) 
  pd(53,11) = pd(53,11) + rrt(177) * density(38) 
  pd(53,38) = pd(53,38) + rrt(177) * density(11) 
  pd(01,11) = pd(01,11) + rrt(178) * density(48) 
  pd(01,48) = pd(01,48) + rrt(178) * density(11) 
  pd(11,11) = pd(11,11) - rrt(178) * density(48) 
  pd(11,48) = pd(11,48) - rrt(178) * density(11) 
  pd(40,11) = pd(40,11) + rrt(178) * density(48) 
  pd(40,48) = pd(40,48) + rrt(178) * density(11) 
  pd(48,11) = pd(48,11) - rrt(178) * density(48) 
  pd(48,48) = pd(48,48) - rrt(178) * density(11) 
  pd(53,11) = pd(53,11) + rrt(178) * density(48) 
  pd(53,48) = pd(53,48) + rrt(178) * density(11) 
  pd(01,11) = pd(01,11) + rrt(179) * density(49) 
  pd(01,49) = pd(01,49) + rrt(179) * density(11) 
  pd(11,11) = pd(11,11) - rrt(179) * density(49) 
  pd(11,49) = pd(11,49) - rrt(179) * density(11) 
  pd(41,11) = pd(41,11) + rrt(179) * density(49) 
  pd(41,49) = pd(41,49) + rrt(179) * density(11) 
  pd(49,11) = pd(49,11) - rrt(179) * density(49) 
  pd(49,49) = pd(49,49) - rrt(179) * density(11) 
  pd(53,11) = pd(53,11) + rrt(179) * density(49) 
  pd(53,49) = pd(53,49) + rrt(179) * density(11) 
  pd(01,11) = pd(01,11) + rrt(180) * density(50) 
  pd(01,50) = pd(01,50) + rrt(180) * density(11) 
  pd(11,11) = pd(11,11) - rrt(180) * density(50) 
  pd(11,50) = pd(11,50) - rrt(180) * density(11) 
  pd(42,11) = pd(42,11) + rrt(180) * density(50) 
  pd(42,50) = pd(42,50) + rrt(180) * density(11) 
  pd(50,11) = pd(50,11) - rrt(180) * density(50) 
  pd(50,50) = pd(50,50) - rrt(180) * density(11) 
  pd(53,11) = pd(53,11) + rrt(180) * density(50) 
  pd(53,50) = pd(53,50) + rrt(180) * density(11) 
  pd(01,11) = pd(01,11) + rrt(181) * density(51) 
  pd(01,51) = pd(01,51) + rrt(181) * density(11) 
  pd(11,11) = pd(11,11) - rrt(181) * density(51) 
  pd(11,51) = pd(11,51) - rrt(181) * density(11) 
  pd(43,11) = pd(43,11) + rrt(181) * density(51) 
  pd(43,51) = pd(43,51) + rrt(181) * density(11) 
  pd(51,11) = pd(51,11) - rrt(181) * density(51) 
  pd(51,51) = pd(51,51) - rrt(181) * density(11) 
  pd(53,11) = pd(53,11) + rrt(181) * density(51) 
  pd(53,51) = pd(53,51) + rrt(181) * density(11) 
  pd(01,10) = pd(01,10) + rrt(182) 
  pd(10,10) = pd(10,10) - rrt(182) 
  pd(10,11) = pd(10,11) + rrt(183) 
  pd(11,11) = pd(11,11) - rrt(183) 
  pd(01,12) = pd(01,12) + rrt(184) 
  pd(12,12) = pd(12,12) - rrt(184) 
  pd(11,13) = pd(11,13) + rrt(185) 
  pd(13,13) = pd(13,13) - rrt(185) 
  pd(21,26) = pd(21,26) + rrt(186) 
  pd(26,26) = pd(26,26) - rrt(186) 
  pd(26,27) = pd(26,27) + rrt(187) 
  pd(27,27) = pd(27,27) - rrt(187) 
  pd(21,27) = pd(21,27) + rrt(188) 
  pd(27,27) = pd(27,27) - rrt(188) 
  pd(21,28) = pd(21,28) + rrt(189) 
  pd(28,28) = pd(28,28) - rrt(189) 
  pd(10,10) = pd(10,10) - rrt(190) * density(29) 
  pd(10,29) = pd(10,29) - rrt(190) * density(10) 
  pd(15,10) = pd(15,10) + rrt(190) * density(29) 
  pd(15,29) = pd(15,29) + rrt(190) * density(10) 
  pd(29,10) = pd(29,10) - rrt(190) * density(29) 
  pd(29,29) = pd(29,29) - rrt(190) * density(10) 
  pd(40,10) = pd(40,10) + rrt(190) * density(29) 
  pd(40,29) = pd(40,29) + rrt(190) * density(10) 
  pd(01,10) = pd(01,10) + rrt(191) * density(29) 
  pd(01,29) = pd(01,29) + rrt(191) * density(10) 
  pd(10,10) = pd(10,10) - rrt(191) * density(29) 
  pd(10,29) = pd(10,29) - rrt(191) * density(10) 
  pd(29,10) = pd(29,10) - rrt(191) * density(29) 
  pd(29,29) = pd(29,29) - rrt(191) * density(10) 
  pd(31,10) = pd(31,10) + rrt(191) * density(29) 
  pd(31,29) = pd(31,29) + rrt(191) * density(10) 
  pd(01,10) = pd(01,10) + rrt(192) * density(14) 
  pd(01,14) = pd(01,14) + rrt(192) * density(10) 
  pd(10,10) = pd(10,10) - rrt(192) * density(14) 
  pd(10,14) = pd(10,14) - rrt(192) * density(10) 
  pd(01,10) = pd(01,10) + rrt(193) * density(14) 
  pd(01,14) = pd(01,14) + rrt(193) * density(10) 
  pd(10,10) = pd(10,10) - rrt(193) * density(14) 
  pd(10,14) = pd(10,14) - rrt(193) * density(10) 
  pd(14,10) = pd(14,10) - rrt(193) * density(14) 
  pd(14,14) = pd(14,14) - rrt(193) * density(10) 
  pd(16,10) = pd(16,10) + rrt(193) * density(14) 
  pd(16,14) = pd(16,14) + rrt(193) * density(10) 
  pd(01,10) = pd(01,10) + rrt(194) * density(21) 
  pd(01,21) = pd(01,21) + rrt(194) * density(10) 
  pd(10,10) = pd(10,10) - rrt(194) * density(21) 
  pd(10,21) = pd(10,21) - rrt(194) * density(10) 
  pd(21,10) = pd(21,10) - rrt(194) * density(21) 
  pd(21,21) = pd(21,21) - rrt(194) * density(10) 
  pd(29,10) = pd(29,10) + rrt(194) * density(21) 
  pd(29,21) = pd(29,21) + rrt(194) * density(10) 
  pd(30,10) = pd(30,10) + rrt(194) * density(21) 
  pd(30,21) = pd(30,21) + rrt(194) * density(10) 
  pd(01,10) = pd(01,10) + rrt(195) * density(21) 
  pd(01,21) = pd(01,21) + rrt(195) * density(10) 
  pd(10,10) = pd(10,10) - rrt(195) * density(21) 
  pd(10,21) = pd(10,21) - rrt(195) * density(10) 
  pd(21,10) = pd(21,10) - rrt(195) * density(21) 
  pd(21,21) = pd(21,21) - rrt(195) * density(10) 
  pd(26,10) = pd(26,10) + rrt(195) * density(21) 
  pd(26,21) = pd(26,21) + rrt(195) * density(10) 
  pd(01,10) = pd(01,10) + rrt(196) * density(21) 
  pd(01,21) = pd(01,21) + rrt(196) * density(10) 
  pd(10,10) = pd(10,10) - rrt(196) * density(21) 
  pd(10,21) = pd(10,21) - rrt(196) * density(10) 
  pd(21,10) = pd(21,10) - rrt(196) * density(21) 
  pd(21,21) = pd(21,21) - rrt(196) * density(10) 
  pd(27,10) = pd(27,10) + rrt(196) * density(21) 
  pd(27,21) = pd(27,21) + rrt(196) * density(10) 
  pd(10,10) = pd(10,10) - rrt(197) * density(21) 
  pd(10,21) = pd(10,21) - rrt(197) * density(10) 
  pd(21,10) = pd(21,10) - rrt(197) * density(21) 
  pd(21,21) = pd(21,21) - rrt(197) * density(10) 
  pd(29,10) = pd(29,10) + rrt(197) * density(21) 
  pd(29,21) = pd(29,21) + rrt(197) * density(10) 
  pd(41,10) = pd(41,10) + rrt(197) * density(21) 
  pd(41,21) = pd(41,21) + rrt(197) * density(10) 
  pd(01,01) = pd(01,01) + rrt(198) * density(10) 
  pd(01,10) = pd(01,10) + rrt(198) * density(01) 
  pd(10,01) = pd(10,01) - rrt(198) * density(10) 
  pd(10,10) = pd(10,10) - rrt(198) * density(01) 
  pd(01,10) = pd(01,10) + rrt(199) * density(40) 
  pd(01,40) = pd(01,40) + rrt(199) * density(10) 
  pd(10,10) = pd(10,10) - rrt(199) * density(40) 
  pd(10,40) = pd(10,40) - rrt(199) * density(10) 
  pd(01,10) = pd(01,10) + rrt(200) * density(41) 
  pd(01,41) = pd(01,41) + rrt(200) * density(10) 
  pd(10,10) = pd(10,10) - rrt(200) * density(41) 
  pd(10,41) = pd(10,41) - rrt(200) * density(10) 
  pd(14,10) = pd(14,10) + rrt(200) * density(41) 
  pd(14,41) = pd(14,41) + rrt(200) * density(10) 
  pd(40,10) = pd(40,10) + rrt(200) * density(41) 
  pd(40,41) = pd(40,41) + rrt(200) * density(10) 
  pd(41,10) = pd(41,10) - rrt(200) * density(41) 
  pd(41,41) = pd(41,41) - rrt(200) * density(10) 
  pd(01,10) = pd(01,10) + rrt(201) * density(42) 
  pd(01,42) = pd(01,42) + rrt(201) * density(10) 
  pd(10,10) = pd(10,10) - rrt(201) * density(42) 
  pd(10,42) = pd(10,42) - rrt(201) * density(10) 
  pd(29,10) = pd(29,10) + rrt(201) * density(42) 
  pd(29,42) = pd(29,42) + rrt(201) * density(10) 
  pd(40,10) = pd(40,10) + rrt(201) * density(42) 
  pd(40,42) = pd(40,42) + rrt(201) * density(10) 
  pd(42,10) = pd(42,10) - rrt(201) * density(42) 
  pd(42,42) = pd(42,42) - rrt(201) * density(10) 
  pd(01,10) = pd(01,10) + rrt(202) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(202) * density(10) * 4.0d0
  pd(11,10) = pd(11,10) + rrt(202) * density(10) * 2.0d0
  pd(01,10) = pd(01,10) + rrt(203) * density(10) * 2.0d0
  pd(10,10) = pd(10,10) - rrt(203) * density(10) * 4.0d0
  pd(13,10) = pd(13,10) + rrt(203) * density(10) * 2.0d0
  pd(10,01) = pd(10,01) + rrt(204) * density(11) 
  pd(10,11) = pd(10,11) + rrt(204) * density(01) 
  pd(11,01) = pd(11,01) - rrt(204) * density(11) 
  pd(11,11) = pd(11,11) - rrt(204) * density(01) 
  pd(01,01) = pd(01,01) + rrt(205) * density(11) 
  pd(01,11) = pd(01,11) + rrt(205) * density(01) 
  pd(11,01) = pd(11,01) - rrt(205) * density(11) 
  pd(11,11) = pd(11,11) - rrt(205) * density(01) 
  pd(01,11) = pd(01,11) + rrt(206) * density(21) 
  pd(01,21) = pd(01,21) + rrt(206) * density(11) 
  pd(11,11) = pd(11,11) - rrt(206) * density(21) 
  pd(11,21) = pd(11,21) - rrt(206) * density(11) 
  pd(21,11) = pd(21,11) - rrt(206) * density(21) 
  pd(21,21) = pd(21,21) - rrt(206) * density(11) 
  pd(29,11) = pd(29,11) + rrt(206) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(206) * density(11) * 2.0d0
  pd(10,11) = pd(10,11) + rrt(207) * density(40) 
  pd(10,40) = pd(10,40) + rrt(207) * density(11) 
  pd(11,11) = pd(11,11) - rrt(207) * density(40) 
  pd(11,40) = pd(11,40) - rrt(207) * density(11) 
  pd(12,01) = pd(12,01) + rrt(208) * density(13) 
  pd(12,13) = pd(12,13) + rrt(208) * density(01) 
  pd(13,01) = pd(13,01) - rrt(208) * density(13) 
  pd(13,13) = pd(13,13) - rrt(208) * density(01) 
  pd(01,13) = pd(01,13) + rrt(209) * density(21) 
  pd(01,21) = pd(01,21) + rrt(209) * density(13) 
  pd(13,13) = pd(13,13) - rrt(209) * density(21) 
  pd(13,21) = pd(13,21) - rrt(209) * density(13) 
  pd(21,13) = pd(21,13) - rrt(209) * density(21) 
  pd(21,21) = pd(21,21) - rrt(209) * density(13) 
  pd(29,13) = pd(29,13) + rrt(209) * density(21) 
  pd(29,21) = pd(29,21) + rrt(209) * density(13) 
  pd(31,13) = pd(31,13) + rrt(209) * density(21) 
  pd(31,21) = pd(31,21) + rrt(209) * density(13) 
  pd(11,01) = pd(11,01) + rrt(210) * density(12) 
  pd(11,12) = pd(11,12) + rrt(210) * density(01) 
  pd(12,01) = pd(12,01) - rrt(210) * density(12) 
  pd(12,12) = pd(12,12) - rrt(210) * density(01) 
  pd(01,12) = pd(01,12) + rrt(211) * density(21) 
  pd(01,21) = pd(01,21) + rrt(211) * density(12) 
  pd(12,12) = pd(12,12) - rrt(211) * density(21) 
  pd(12,21) = pd(12,21) - rrt(211) * density(12) 
  pd(21,12) = pd(21,12) - rrt(211) * density(21) 
  pd(21,21) = pd(21,21) - rrt(211) * density(12) 
  pd(29,12) = pd(29,12) + rrt(211) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(211) * density(12) * 2.0d0
  pd(01,12) = pd(01,12) + rrt(212) * density(40) 
  pd(01,40) = pd(01,40) + rrt(212) * density(12) 
  pd(12,12) = pd(12,12) - rrt(212) * density(40) 
  pd(12,40) = pd(12,40) - rrt(212) * density(12) 
  pd(14,12) = pd(14,12) + rrt(212) * density(40) 
  pd(14,40) = pd(14,40) + rrt(212) * density(12) 
  pd(29,12) = pd(29,12) + rrt(212) * density(40) 
  pd(29,40) = pd(29,40) + rrt(212) * density(12) 
  pd(40,12) = pd(40,12) - rrt(212) * density(40) 
  pd(40,40) = pd(40,40) - rrt(212) * density(12) 
  pd(10,10) = pd(10,10) - rrt(213) * density(12) 
  pd(10,12) = pd(10,12) - rrt(213) * density(10) 
  pd(12,10) = pd(12,10) - rrt(213) * density(12) 
  pd(12,12) = pd(12,12) - rrt(213) * density(10) 
  pd(20,10) = pd(20,10) + rrt(213) * density(12) 
  pd(20,12) = pd(20,12) + rrt(213) * density(10) 
  pd(53,10) = pd(53,10) + rrt(213) * density(12) 
  pd(53,12) = pd(53,12) + rrt(213) * density(10) 
  pd(12,12) = pd(12,12) - rrt(214) * density(12) * 4.0d0
  pd(20,12) = pd(20,12) + rrt(214) * density(12) * 2.0d0
  pd(53,12) = pd(53,12) + rrt(214) * density(12) * 2.0d0
  pd(10,01) = pd(10,01) + rrt(215) * density(14)**2 
  pd(10,14) = pd(10,14) + rrt(215) * density(01) * density(14) * 2.0d0
  pd(14,01) = pd(14,01) - rrt(215) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(215) * density(01) * density(14) * 4.0d0
  pd(10,14) = pd(10,14) + rrt(216) * density(14) * density(21) * 2.0d0
  pd(10,21) = pd(10,21) + rrt(216) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(216) * density(14) * density(21) * 4.0d0
  pd(14,21) = pd(14,21) - rrt(216) * density(14)**2 * 2.0d0
  pd(10,14) = pd(10,14) + rrt(217) * density(14) * density(40) * 2.0d0
  pd(10,40) = pd(10,40) + rrt(217) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(217) * density(14) * density(40) * 4.0d0
  pd(14,40) = pd(14,40) - rrt(217) * density(14)**2 * 2.0d0
  pd(10,14) = pd(10,14) + rrt(218) * density(14)**2 * 3.0d0
  pd(14,14) = pd(14,14) - rrt(218) * density(14)**2 * 6.0d0
  pd(10,14) = pd(10,14) + rrt(219) * density(14) * density(29) * 2.0d0
  pd(10,29) = pd(10,29) + rrt(219) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(219) * density(14) * density(29) * 4.0d0
  pd(14,29) = pd(14,29) - rrt(219) * density(14)**2 * 2.0d0
  pd(11,01) = pd(11,01) + rrt(220) * density(14)**2 
  pd(11,14) = pd(11,14) + rrt(220) * density(01) * density(14) * 2.0d0
  pd(14,01) = pd(14,01) - rrt(220) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(220) * density(01) * density(14) * 4.0d0
  pd(11,14) = pd(11,14) + rrt(221) * density(14) * density(21) * 2.0d0
  pd(11,21) = pd(11,21) + rrt(221) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(221) * density(14) * density(21) * 4.0d0
  pd(14,21) = pd(14,21) - rrt(221) * density(14)**2 * 2.0d0
  pd(11,14) = pd(11,14) + rrt(222) * density(14) * density(40) * 2.0d0
  pd(11,40) = pd(11,40) + rrt(222) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(222) * density(14) * density(40) * 4.0d0
  pd(14,40) = pd(14,40) - rrt(222) * density(14)**2 * 2.0d0
  pd(11,14) = pd(11,14) + rrt(223) * density(14)**2 * 3.0d0
  pd(14,14) = pd(14,14) - rrt(223) * density(14)**2 * 6.0d0
  pd(11,14) = pd(11,14) + rrt(224) * density(14) * density(29) * 2.0d0
  pd(11,29) = pd(11,29) + rrt(224) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(224) * density(14) * density(29) * 4.0d0
  pd(14,29) = pd(14,29) - rrt(224) * density(14)**2 * 2.0d0
  pd(14,15) = pd(14,15) + rrt(225) * density(29) 
  pd(14,29) = pd(14,29) + rrt(225) * density(15) 
  pd(15,15) = pd(15,15) - rrt(225) * density(29) 
  pd(15,29) = pd(15,29) - rrt(225) * density(15) 
  pd(29,15) = pd(29,15) - rrt(225) * density(29) 
  pd(29,29) = pd(29,29) - rrt(225) * density(15) 
  pd(30,15) = pd(30,15) + rrt(225) * density(29) 
  pd(30,29) = pd(30,29) + rrt(225) * density(15) 
  pd(15,15) = pd(15,15) - rrt(226) * density(21) 
  pd(15,21) = pd(15,21) - rrt(226) * density(15) 
  pd(21,15) = pd(21,15) - rrt(226) * density(21) 
  pd(21,21) = pd(21,21) - rrt(226) * density(15) 
  pd(29,15) = pd(29,15) + rrt(226) * density(21) 
  pd(29,21) = pd(29,21) + rrt(226) * density(15) 
  pd(40,15) = pd(40,15) + rrt(226) * density(21) 
  pd(40,21) = pd(40,21) + rrt(226) * density(15) 
  pd(01,15) = pd(01,15) + rrt(227) * density(40) 
  pd(01,40) = pd(01,40) + rrt(227) * density(15) 
  pd(15,15) = pd(15,15) - rrt(227) * density(40) 
  pd(15,40) = pd(15,40) - rrt(227) * density(15) 
  pd(29,15) = pd(29,15) + rrt(227) * density(40) 
  pd(29,40) = pd(29,40) + rrt(227) * density(15) 
  pd(40,15) = pd(40,15) - rrt(227) * density(40) 
  pd(40,40) = pd(40,40) - rrt(227) * density(15) 
  pd(01,15) = pd(01,15) + rrt(228) * density(41) 
  pd(01,41) = pd(01,41) + rrt(228) * density(15) 
  pd(15,15) = pd(15,15) - rrt(228) * density(41) 
  pd(15,41) = pd(15,41) - rrt(228) * density(15) 
  pd(40,15) = pd(40,15) + rrt(228) * density(41) 
  pd(40,41) = pd(40,41) + rrt(228) * density(15) 
  pd(41,15) = pd(41,15) - rrt(228) * density(41) 
  pd(41,41) = pd(41,41) - rrt(228) * density(15) 
  pd(14,01) = pd(14,01) + rrt(229) * density(15) 
  pd(14,15) = pd(14,15) + rrt(229) * density(01) 
  pd(15,01) = pd(15,01) - rrt(229) * density(15) 
  pd(15,15) = pd(15,15) - rrt(229) * density(01) 
  pd(14,14) = pd(14,14) + rrt(230) * density(16) 
  pd(14,16) = pd(14,16) + rrt(230) * density(14) 
  pd(16,14) = pd(16,14) - rrt(230) * density(16) 
  pd(16,16) = pd(16,16) - rrt(230) * density(14) 
  pd(14,16) = pd(14,16) + rrt(231) * density(29) 
  pd(14,29) = pd(14,29) + rrt(231) * density(16) 
  pd(16,16) = pd(16,16) - rrt(231) * density(29) 
  pd(16,29) = pd(16,29) - rrt(231) * density(16) 
  pd(15,14) = pd(15,14) + rrt(232) * density(16) 
  pd(15,16) = pd(15,16) + rrt(232) * density(14) 
  pd(16,14) = pd(16,14) - rrt(232) * density(16) 
  pd(16,16) = pd(16,16) - rrt(232) * density(14) 
  pd(14,01) = pd(14,01) + rrt(233) * density(16) 
  pd(14,16) = pd(14,16) + rrt(233) * density(01) 
  pd(16,01) = pd(16,01) - rrt(233) * density(16) 
  pd(16,16) = pd(16,16) - rrt(233) * density(01) 
  pd(15,15) = pd(15,15) - rrt(234) * density(16) 
  pd(15,16) = pd(15,16) - rrt(234) * density(15) 
  pd(16,15) = pd(16,15) - rrt(234) * density(16) 
  pd(16,16) = pd(16,16) - rrt(234) * density(15) 
  pd(18,15) = pd(18,15) + rrt(234) * density(16) 
  pd(18,16) = pd(18,16) + rrt(234) * density(15) 
  pd(53,15) = pd(53,15) + rrt(234) * density(16) 
  pd(53,16) = pd(53,16) + rrt(234) * density(15) 
  pd(16,16) = pd(16,16) - rrt(235) * density(21) 
  pd(16,21) = pd(16,21) - rrt(235) * density(16) 
  pd(21,16) = pd(21,16) - rrt(235) * density(21) 
  pd(21,21) = pd(21,21) - rrt(235) * density(16) 
  pd(29,16) = pd(29,16) + rrt(235) * density(21) 
  pd(29,21) = pd(29,21) + rrt(235) * density(16) 
  pd(40,16) = pd(40,16) + rrt(235) * density(21) 
  pd(40,21) = pd(40,21) + rrt(235) * density(16) 
  pd(10,16) = pd(10,16) + rrt(236) * density(40) 
  pd(10,40) = pd(10,40) + rrt(236) * density(16) 
  pd(16,16) = pd(16,16) - rrt(236) * density(40) 
  pd(16,40) = pd(16,40) - rrt(236) * density(16) 
  pd(29,16) = pd(29,16) + rrt(236) * density(40) 
  pd(29,40) = pd(29,40) + rrt(236) * density(16) 
  pd(40,16) = pd(40,16) - rrt(236) * density(40) 
  pd(40,40) = pd(40,40) - rrt(236) * density(16) 
  pd(21,26) = pd(21,26) + rrt(237) * density(29) 
  pd(21,29) = pd(21,29) + rrt(237) * density(26) 
  pd(26,26) = pd(26,26) - rrt(237) * density(29) 
  pd(26,29) = pd(26,29) - rrt(237) * density(26) 
  pd(14,14) = pd(14,14) - rrt(238) * density(26) 
  pd(14,26) = pd(14,26) - rrt(238) * density(14) 
  pd(26,14) = pd(26,14) - rrt(238) * density(26) 
  pd(26,26) = pd(26,26) - rrt(238) * density(14) 
  pd(29,14) = pd(29,14) + rrt(238) * density(26) 
  pd(29,26) = pd(29,26) + rrt(238) * density(14) 
  pd(40,14) = pd(40,14) + rrt(238) * density(26) 
  pd(40,26) = pd(40,26) + rrt(238) * density(14) 
  pd(21,21) = pd(21,21) + rrt(239) * density(26) 
  pd(21,26) = pd(21,26) + rrt(239) * density(21) 
  pd(26,21) = pd(26,21) - rrt(239) * density(26) 
  pd(26,26) = pd(26,26) - rrt(239) * density(21) 
  pd(21,01) = pd(21,01) + rrt(240) * density(26) 
  pd(21,26) = pd(21,26) + rrt(240) * density(01) 
  pd(26,01) = pd(26,01) - rrt(240) * density(26) 
  pd(26,26) = pd(26,26) - rrt(240) * density(01) 
  pd(21,26) = pd(21,26) + rrt(241) * density(40) 
  pd(21,40) = pd(21,40) + rrt(241) * density(26) 
  pd(26,26) = pd(26,26) - rrt(241) * density(40) 
  pd(26,40) = pd(26,40) - rrt(241) * density(26) 
  pd(21,26) = pd(21,26) + rrt(242) * density(32) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(242) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(242) * density(32) 
  pd(26,32) = pd(26,32) - rrt(242) * density(26) 
  pd(30,26) = pd(30,26) + rrt(242) * density(32) 
  pd(30,32) = pd(30,32) + rrt(242) * density(26) 
  pd(32,26) = pd(32,26) - rrt(242) * density(32) 
  pd(32,32) = pd(32,32) - rrt(242) * density(26) 
  pd(21,26) = pd(21,26) + rrt(243) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(243) * density(26) * 4.0d0
  pd(27,26) = pd(27,26) + rrt(243) * density(26) * 2.0d0
  pd(21,29) = pd(21,29) + rrt(244) * density(32) 
  pd(21,32) = pd(21,32) + rrt(244) * density(29) 
  pd(26,29) = pd(26,29) + rrt(244) * density(32) 
  pd(26,32) = pd(26,32) + rrt(244) * density(29) 
  pd(29,29) = pd(29,29) - rrt(244) * density(32) 
  pd(29,32) = pd(29,32) - rrt(244) * density(29) 
  pd(32,29) = pd(32,29) - rrt(244) * density(32) 
  pd(32,32) = pd(32,32) - rrt(244) * density(29) 
  pd(26,27) = pd(26,27) + rrt(245) * density(29) 
  pd(26,29) = pd(26,29) + rrt(245) * density(27) 
  pd(27,27) = pd(27,27) - rrt(245) * density(29) 
  pd(27,29) = pd(27,29) - rrt(245) * density(27) 
  pd(21,27) = pd(21,27) + rrt(246) * density(29) 
  pd(21,29) = pd(21,29) + rrt(246) * density(27) 
  pd(27,27) = pd(27,27) - rrt(246) * density(29) 
  pd(27,29) = pd(27,29) - rrt(246) * density(27) 
  pd(29,27) = pd(29,27) - rrt(246) * density(29) 
  pd(29,29) = pd(29,29) - rrt(246) * density(27) 
  pd(30,27) = pd(30,27) + rrt(246) * density(29) 
  pd(30,29) = pd(30,29) + rrt(246) * density(27) 
  pd(26,21) = pd(26,21) + rrt(247) * density(27) 
  pd(26,27) = pd(26,27) + rrt(247) * density(21) 
  pd(27,21) = pd(27,21) - rrt(247) * density(27) 
  pd(27,27) = pd(27,27) - rrt(247) * density(21) 
  pd(26,01) = pd(26,01) + rrt(248) * density(27) 
  pd(26,27) = pd(26,27) + rrt(248) * density(01) 
  pd(27,01) = pd(27,01) - rrt(248) * density(27) 
  pd(27,27) = pd(27,27) - rrt(248) * density(01) 
  pd(26,27) = pd(26,27) + rrt(249) * density(40) 
  pd(26,40) = pd(26,40) + rrt(249) * density(27) 
  pd(27,27) = pd(27,27) - rrt(249) * density(40) 
  pd(27,40) = pd(27,40) - rrt(249) * density(27) 
  pd(21,27) = pd(21,27) + rrt(250) * density(32) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(250) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(250) * density(32) 
  pd(27,32) = pd(27,32) - rrt(250) * density(27) 
  pd(29,27) = pd(29,27) + rrt(250) * density(32) 
  pd(29,32) = pd(29,32) + rrt(250) * density(27) 
  pd(32,27) = pd(32,27) - rrt(250) * density(32) 
  pd(32,32) = pd(32,32) - rrt(250) * density(27) 
  pd(21,28) = pd(21,28) + rrt(251) * density(29) 
  pd(21,29) = pd(21,29) + rrt(251) * density(28) 
  pd(28,28) = pd(28,28) - rrt(251) * density(29) 
  pd(28,29) = pd(28,29) - rrt(251) * density(28) 
  pd(29,28) = pd(29,28) - rrt(251) * density(29) 
  pd(29,29) = pd(29,29) - rrt(251) * density(28) 
  pd(31,28) = pd(31,28) + rrt(251) * density(29) 
  pd(31,29) = pd(31,29) + rrt(251) * density(28) 
  pd(21,21) = pd(21,21) - rrt(252) * density(28) 
  pd(21,28) = pd(21,28) - rrt(252) * density(21) 
  pd(27,21) = pd(27,21) + rrt(252) * density(28) * 2.0d0
  pd(27,28) = pd(27,28) + rrt(252) * density(21) * 2.0d0
  pd(28,21) = pd(28,21) - rrt(252) * density(28) 
  pd(28,28) = pd(28,28) - rrt(252) * density(21) 
  pd(27,01) = pd(27,01) + rrt(253) * density(28) 
  pd(27,28) = pd(27,28) + rrt(253) * density(01) 
  pd(28,01) = pd(28,01) - rrt(253) * density(28) 
  pd(28,28) = pd(28,28) - rrt(253) * density(01) 
  pd(29,29) = pd(29,29) + rrt(254) * density(30) 
  pd(29,30) = pd(29,30) + rrt(254) * density(29) 
  pd(30,29) = pd(30,29) - rrt(254) * density(30) 
  pd(30,30) = pd(30,30) - rrt(254) * density(29) 
  pd(29,21) = pd(29,21) + rrt(255) * density(30) 
  pd(29,30) = pd(29,30) + rrt(255) * density(21) 
  pd(30,21) = pd(30,21) - rrt(255) * density(30) 
  pd(30,30) = pd(30,30) - rrt(255) * density(21) 
  pd(21,21) = pd(21,21) - rrt(256) * density(30) 
  pd(21,30) = pd(21,30) - rrt(256) * density(21) 
  pd(26,21) = pd(26,21) + rrt(256) * density(30) 
  pd(26,30) = pd(26,30) + rrt(256) * density(21) 
  pd(29,21) = pd(29,21) + rrt(256) * density(30) 
  pd(29,30) = pd(29,30) + rrt(256) * density(21) 
  pd(30,21) = pd(30,21) - rrt(256) * density(30) 
  pd(30,30) = pd(30,30) - rrt(256) * density(21) 
  pd(21,21) = pd(21,21) - rrt(257) * density(30) 
  pd(21,30) = pd(21,30) - rrt(257) * density(21) 
  pd(27,21) = pd(27,21) + rrt(257) * density(30) 
  pd(27,30) = pd(27,30) + rrt(257) * density(21) 
  pd(29,21) = pd(29,21) + rrt(257) * density(30) 
  pd(29,30) = pd(29,30) + rrt(257) * density(21) 
  pd(30,21) = pd(30,21) - rrt(257) * density(30) 
  pd(30,30) = pd(30,30) - rrt(257) * density(21) 
  pd(29,01) = pd(29,01) + rrt(258) * density(30) 
  pd(29,30) = pd(29,30) + rrt(258) * density(01) 
  pd(30,01) = pd(30,01) - rrt(258) * density(30) 
  pd(30,30) = pd(30,30) - rrt(258) * density(01) 
  pd(21,30) = pd(21,30) + rrt(259) * density(32) 
  pd(21,32) = pd(21,32) + rrt(259) * density(30) 
  pd(29,30) = pd(29,30) + rrt(259) * density(32) * 2.0d0
  pd(29,32) = pd(29,32) + rrt(259) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(259) * density(32) 
  pd(30,32) = pd(30,32) - rrt(259) * density(30) 
  pd(32,30) = pd(32,30) - rrt(259) * density(32) 
  pd(32,32) = pd(32,32) - rrt(259) * density(30) 
  pd(21,30) = pd(21,30) + rrt(260) * density(32) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(260) * density(30) * 2.0d0
  pd(30,30) = pd(30,30) - rrt(260) * density(32) 
  pd(30,32) = pd(30,32) - rrt(260) * density(30) 
  pd(32,30) = pd(32,30) - rrt(260) * density(32) 
  pd(32,32) = pd(32,32) - rrt(260) * density(30) 
  pd(14,30) = pd(14,30) + rrt(261) * density(40) 
  pd(14,40) = pd(14,40) + rrt(261) * density(30) 
  pd(21,30) = pd(21,30) + rrt(261) * density(40) 
  pd(21,40) = pd(21,40) + rrt(261) * density(30) 
  pd(30,30) = pd(30,30) - rrt(261) * density(40) 
  pd(30,40) = pd(30,40) - rrt(261) * density(30) 
  pd(40,30) = pd(40,30) - rrt(261) * density(40) 
  pd(40,40) = pd(40,40) - rrt(261) * density(30) 
  pd(30,30) = pd(30,30) - rrt(262) * density(41) 
  pd(30,41) = pd(30,41) - rrt(262) * density(30) 
  pd(40,30) = pd(40,30) + rrt(262) * density(41) * 2.0d0
  pd(40,41) = pd(40,41) + rrt(262) * density(30) * 2.0d0
  pd(41,30) = pd(41,30) - rrt(262) * density(41) 
  pd(41,41) = pd(41,41) - rrt(262) * density(30) 
  pd(01,30) = pd(01,30) + rrt(263) * density(41) 
  pd(01,41) = pd(01,41) + rrt(263) * density(30) 
  pd(21,30) = pd(21,30) + rrt(263) * density(41) 
  pd(21,41) = pd(21,41) + rrt(263) * density(30) 
  pd(30,30) = pd(30,30) - rrt(263) * density(41) 
  pd(30,41) = pd(30,41) - rrt(263) * density(30) 
  pd(41,30) = pd(41,30) - rrt(263) * density(41) 
  pd(41,41) = pd(41,41) - rrt(263) * density(30) 
  pd(30,29) = pd(30,29) + rrt(264) * density(31) 
  pd(30,31) = pd(30,31) + rrt(264) * density(29) 
  pd(31,29) = pd(31,29) - rrt(264) * density(31) 
  pd(31,31) = pd(31,31) - rrt(264) * density(29) 
  pd(29,14) = pd(29,14) + rrt(265) * density(31) 
  pd(29,31) = pd(29,31) + rrt(265) * density(14) 
  pd(31,14) = pd(31,14) - rrt(265) * density(31) 
  pd(31,31) = pd(31,31) - rrt(265) * density(14) 
  pd(30,21) = pd(30,21) + rrt(266) * density(31) 
  pd(30,31) = pd(30,31) + rrt(266) * density(21) 
  pd(31,21) = pd(31,21) - rrt(266) * density(31) 
  pd(31,31) = pd(31,31) - rrt(266) * density(21) 
  pd(21,21) = pd(21,21) - rrt(267) * density(31) 
  pd(21,31) = pd(21,31) - rrt(267) * density(21) 
  pd(29,21) = pd(29,21) + rrt(267) * density(31) * 3.0d0
  pd(29,31) = pd(29,31) + rrt(267) * density(21) * 3.0d0
  pd(31,21) = pd(31,21) - rrt(267) * density(31) 
  pd(31,31) = pd(31,31) - rrt(267) * density(21) 
  pd(29,01) = pd(29,01) + rrt(268) * density(31) 
  pd(29,31) = pd(29,31) + rrt(268) * density(01) 
  pd(31,01) = pd(31,01) - rrt(268) * density(31) 
  pd(31,31) = pd(31,31) - rrt(268) * density(01) 
  pd(26,26) = pd(26,26) - rrt(269) * density(31) 
  pd(26,31) = pd(26,31) - rrt(269) * density(26) 
  pd(28,26) = pd(28,26) + rrt(269) * density(31) 
  pd(28,31) = pd(28,31) + rrt(269) * density(26) 
  pd(29,26) = pd(29,26) + rrt(269) * density(31) 
  pd(29,31) = pd(29,31) + rrt(269) * density(26) 
  pd(31,26) = pd(31,26) - rrt(269) * density(31) 
  pd(31,31) = pd(31,31) - rrt(269) * density(26) 
  pd(26,26) = pd(26,26) - rrt(270) * density(31) 
  pd(26,31) = pd(26,31) - rrt(270) * density(26) 
  pd(27,26) = pd(27,26) + rrt(270) * density(31) 
  pd(27,31) = pd(27,31) + rrt(270) * density(26) 
  pd(30,26) = pd(30,26) + rrt(270) * density(31) 
  pd(30,31) = pd(30,31) + rrt(270) * density(26) 
  pd(31,26) = pd(31,26) - rrt(270) * density(31) 
  pd(31,31) = pd(31,31) - rrt(270) * density(26) 
  pd(26,26) = pd(26,26) - rrt(271) * density(31) 
  pd(26,31) = pd(26,31) - rrt(271) * density(26) 
  pd(29,26) = pd(29,26) + rrt(271) * density(31) * 3.0d0
  pd(29,31) = pd(29,31) + rrt(271) * density(26) * 3.0d0
  pd(31,26) = pd(31,26) - rrt(271) * density(31) 
  pd(31,31) = pd(31,31) - rrt(271) * density(26) 
  pd(29,31) = pd(29,31) + rrt(272) * density(40) 
  pd(29,40) = pd(29,40) + rrt(272) * density(31) 
  pd(31,31) = pd(31,31) - rrt(272) * density(40) 
  pd(31,40) = pd(31,40) - rrt(272) * density(31) 
  pd(30,31) = pd(30,31) + rrt(273) * density(40) 
  pd(30,40) = pd(30,40) + rrt(273) * density(31) 
  pd(31,31) = pd(31,31) - rrt(273) * density(40) 
  pd(31,40) = pd(31,40) - rrt(273) * density(31) 
  pd(21,31) = pd(21,31) + rrt(274) * density(32) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(274) * density(31) * 2.0d0
  pd(31,31) = pd(31,31) - rrt(274) * density(32) 
  pd(31,32) = pd(31,32) - rrt(274) * density(31) 
  pd(32,31) = pd(32,31) - rrt(274) * density(32) 
  pd(32,32) = pd(32,32) - rrt(274) * density(31) 
  pd(21,31) = pd(21,31) + rrt(275) * density(32) 
  pd(21,32) = pd(21,32) + rrt(275) * density(31) 
  pd(29,31) = pd(29,31) + rrt(275) * density(32) 
  pd(29,32) = pd(29,32) + rrt(275) * density(31) 
  pd(30,31) = pd(30,31) + rrt(275) * density(32) 
  pd(30,32) = pd(30,32) + rrt(275) * density(31) 
  pd(31,31) = pd(31,31) - rrt(275) * density(32) 
  pd(31,32) = pd(31,32) - rrt(275) * density(31) 
  pd(32,31) = pd(32,31) - rrt(275) * density(32) 
  pd(32,32) = pd(32,32) - rrt(275) * density(31) 
  pd(29,31) = pd(29,31) + rrt(276) * density(41) 
  pd(29,41) = pd(29,41) + rrt(276) * density(31) 
  pd(31,31) = pd(31,31) - rrt(276) * density(41) 
  pd(31,41) = pd(31,41) - rrt(276) * density(31) 
  pd(30,31) = pd(30,31) + rrt(277) * density(41) 
  pd(30,41) = pd(30,41) + rrt(277) * density(31) 
  pd(31,31) = pd(31,31) - rrt(277) * density(41) 
  pd(31,41) = pd(31,41) - rrt(277) * density(31) 
  pd(01,14) = pd(01,14) + rrt(278) * density(40) 
  pd(01,40) = pd(01,40) + rrt(278) * density(14) 
  pd(14,14) = pd(14,14) - rrt(278) * density(40) 
  pd(14,40) = pd(14,40) - rrt(278) * density(14) 
  pd(29,14) = pd(29,14) + rrt(278) * density(40) 
  pd(29,40) = pd(29,40) + rrt(278) * density(14) 
  pd(40,14) = pd(40,14) - rrt(278) * density(40) 
  pd(40,40) = pd(40,40) - rrt(278) * density(14) 
  pd(14,14) = pd(14,14) - rrt(279) * density(21) 
  pd(14,21) = pd(14,21) - rrt(279) * density(14) 
  pd(21,14) = pd(21,14) - rrt(279) * density(21) 
  pd(21,21) = pd(21,21) - rrt(279) * density(14) 
  pd(29,14) = pd(29,14) + rrt(279) * density(21) 
  pd(29,21) = pd(29,21) + rrt(279) * density(14) 
  pd(40,14) = pd(40,14) + rrt(279) * density(21) 
  pd(40,21) = pd(40,21) + rrt(279) * density(14) 
  pd(01,14) = pd(01,14) + rrt(280) * density(42) 
  pd(01,42) = pd(01,42) + rrt(280) * density(14) 
  pd(14,14) = pd(14,14) - rrt(280) * density(42) 
  pd(14,42) = pd(14,42) - rrt(280) * density(14) 
  pd(29,14) = pd(29,14) + rrt(280) * density(42) * 2.0d0
  pd(29,42) = pd(29,42) + rrt(280) * density(14) * 2.0d0
  pd(42,14) = pd(42,14) - rrt(280) * density(42) 
  pd(42,42) = pd(42,42) - rrt(280) * density(14) 
  pd(14,14) = pd(14,14) - rrt(281) * density(42) 
  pd(14,42) = pd(14,42) - rrt(281) * density(14) 
  pd(29,14) = pd(29,14) + rrt(281) * density(42) 
  pd(29,42) = pd(29,42) + rrt(281) * density(14) 
  pd(41,14) = pd(41,14) + rrt(281) * density(42) 
  pd(41,42) = pd(41,42) + rrt(281) * density(14) 
  pd(42,14) = pd(42,14) - rrt(281) * density(42) 
  pd(42,42) = pd(42,42) - rrt(281) * density(14) 
  pd(01,14) = pd(01,14) + rrt(282) * density(42) 
  pd(01,42) = pd(01,42) + rrt(282) * density(14) 
  pd(14,14) = pd(14,14) - rrt(282) * density(42) 
  pd(14,42) = pd(14,42) - rrt(282) * density(14) 
  pd(21,14) = pd(21,14) + rrt(282) * density(42) 
  pd(21,42) = pd(21,42) + rrt(282) * density(14) 
  pd(42,14) = pd(42,14) - rrt(282) * density(42) 
  pd(42,42) = pd(42,42) - rrt(282) * density(14) 
  pd(14,14) = pd(14,14) - rrt(283) * density(42) 
  pd(14,42) = pd(14,42) - rrt(283) * density(14) 
  pd(40,14) = pd(40,14) + rrt(283) * density(42) * 2.0d0
  pd(40,42) = pd(40,42) + rrt(283) * density(14) * 2.0d0
  pd(42,14) = pd(42,14) - rrt(283) * density(42) 
  pd(42,42) = pd(42,42) - rrt(283) * density(14) 
  pd(01,01) = pd(01,01) - rrt(284) * density(29) 
  pd(01,29) = pd(01,29) - rrt(284) * density(01) 
  pd(14,01) = pd(14,01) + rrt(284) * density(29) 
  pd(14,29) = pd(14,29) + rrt(284) * density(01) 
  pd(29,01) = pd(29,01) - rrt(284) * density(29) 
  pd(29,29) = pd(29,29) - rrt(284) * density(01) 
  pd(40,01) = pd(40,01) + rrt(284) * density(29) 
  pd(40,29) = pd(40,29) + rrt(284) * density(01) 
  pd(14,29) = pd(14,29) + rrt(285) * density(40) 
  pd(14,40) = pd(14,40) + rrt(285) * density(29) 
  pd(21,29) = pd(21,29) + rrt(285) * density(40) 
  pd(21,40) = pd(21,40) + rrt(285) * density(29) 
  pd(29,29) = pd(29,29) - rrt(285) * density(40) 
  pd(29,40) = pd(29,40) - rrt(285) * density(29) 
  pd(40,29) = pd(40,29) - rrt(285) * density(40) 
  pd(40,40) = pd(40,40) - rrt(285) * density(29) 
  pd(29,29) = pd(29,29) - rrt(286) * density(40) 
  pd(29,40) = pd(29,40) - rrt(286) * density(29) 
  pd(40,29) = pd(40,29) - rrt(286) * density(40) 
  pd(40,40) = pd(40,40) - rrt(286) * density(29) 
  pd(42,29) = pd(42,29) + rrt(286) * density(40) 
  pd(42,40) = pd(42,40) + rrt(286) * density(29) 
  pd(01,29) = pd(01,29) + rrt(287) * density(41) 
  pd(01,41) = pd(01,41) + rrt(287) * density(29) 
  pd(21,29) = pd(21,29) + rrt(287) * density(41) 
  pd(21,41) = pd(21,41) + rrt(287) * density(29) 
  pd(29,29) = pd(29,29) - rrt(287) * density(41) 
  pd(29,41) = pd(29,41) - rrt(287) * density(29) 
  pd(41,29) = pd(41,29) - rrt(287) * density(41) 
  pd(41,41) = pd(41,41) - rrt(287) * density(29) 
  pd(29,29) = pd(29,29) - rrt(288) * density(41) 
  pd(29,41) = pd(29,41) - rrt(288) * density(29) 
  pd(40,29) = pd(40,29) + rrt(288) * density(41) * 2.0d0
  pd(40,41) = pd(40,41) + rrt(288) * density(29) * 2.0d0
  pd(41,29) = pd(41,29) - rrt(288) * density(41) 
  pd(41,41) = pd(41,41) - rrt(288) * density(29) 
  pd(21,29) = pd(21,29) + rrt(289) * density(42) 
  pd(21,42) = pd(21,42) + rrt(289) * density(29) 
  pd(29,29) = pd(29,29) - rrt(289) * density(42) 
  pd(29,42) = pd(29,42) - rrt(289) * density(29) 
  pd(40,29) = pd(40,29) + rrt(289) * density(42) 
  pd(40,42) = pd(40,42) + rrt(289) * density(29) 
  pd(42,29) = pd(42,29) - rrt(289) * density(42) 
  pd(42,42) = pd(42,42) - rrt(289) * density(29) 
  pd(21,29) = pd(21,29) + rrt(290) * density(43) 
  pd(21,43) = pd(21,43) + rrt(290) * density(29) 
  pd(29,29) = pd(29,29) - rrt(290) * density(43) 
  pd(29,43) = pd(29,43) - rrt(290) * density(29) 
  pd(42,29) = pd(42,29) + rrt(290) * density(43) 
  pd(42,43) = pd(42,43) + rrt(290) * density(29) 
  pd(43,29) = pd(43,29) - rrt(290) * density(43) 
  pd(43,43) = pd(43,43) - rrt(290) * density(29) 
  pd(01,01) = pd(01,01) - rrt(291) * density(21) 
  pd(01,21) = pd(01,21) - rrt(291) * density(01) 
  pd(21,01) = pd(21,01) - rrt(291) * density(21) 
  pd(21,21) = pd(21,21) - rrt(291) * density(01) 
  pd(29,01) = pd(29,01) + rrt(291) * density(21) 
  pd(29,21) = pd(29,21) + rrt(291) * density(01) 
  pd(41,01) = pd(41,01) + rrt(291) * density(21) 
  pd(41,21) = pd(41,21) + rrt(291) * density(01) 
  pd(14,40) = pd(14,40) + rrt(292) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(292) * density(40) * 4.0d0
  pd(42,40) = pd(42,40) + rrt(292) * density(40) * 2.0d0
  pd(29,40) = pd(29,40) + rrt(293) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(293) * density(40) * 4.0d0
  pd(41,40) = pd(41,40) + rrt(293) * density(40) * 2.0d0
  pd(01,40) = pd(01,40) + rrt(294) * density(40) * 2.0d0
  pd(21,40) = pd(21,40) + rrt(294) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(294) * density(40) * 4.0d0
  pd(21,21) = pd(21,21) - rrt(295) * density(40) 
  pd(21,40) = pd(21,40) - rrt(295) * density(21) 
  pd(29,21) = pd(29,21) + rrt(295) * density(40) 
  pd(29,40) = pd(29,40) + rrt(295) * density(21) 
  pd(40,21) = pd(40,21) - rrt(295) * density(40) 
  pd(40,40) = pd(40,40) - rrt(295) * density(21) 
  pd(42,21) = pd(42,21) + rrt(295) * density(40) 
  pd(42,40) = pd(42,40) + rrt(295) * density(21) 
  pd(21,32) = pd(21,32) + rrt(296) * density(40) 
  pd(21,40) = pd(21,40) + rrt(296) * density(32) 
  pd(32,32) = pd(32,32) - rrt(296) * density(40) 
  pd(32,40) = pd(32,40) - rrt(296) * density(32) 
  pd(40,32) = pd(40,32) - rrt(296) * density(40) 
  pd(40,40) = pd(40,40) - rrt(296) * density(32) 
  pd(42,32) = pd(42,32) + rrt(296) * density(40) 
  pd(42,40) = pd(42,40) + rrt(296) * density(32) 
  pd(01,40) = pd(01,40) + rrt(297) * density(41) 
  pd(01,41) = pd(01,41) + rrt(297) * density(40) 
  pd(40,40) = pd(40,40) - rrt(297) * density(41) 
  pd(40,41) = pd(40,41) - rrt(297) * density(40) 
  pd(41,40) = pd(41,40) - rrt(297) * density(41) 
  pd(41,41) = pd(41,41) - rrt(297) * density(40) 
  pd(42,40) = pd(42,40) + rrt(297) * density(41) 
  pd(42,41) = pd(42,41) + rrt(297) * density(40) 
  pd(40,40) = pd(40,40) - rrt(298) * density(43) 
  pd(40,43) = pd(40,43) - rrt(298) * density(40) 
  pd(42,40) = pd(42,40) + rrt(298) * density(43) * 2.0d0
  pd(42,43) = pd(42,43) + rrt(298) * density(40) * 2.0d0
  pd(43,40) = pd(43,40) - rrt(298) * density(43) 
  pd(43,43) = pd(43,43) - rrt(298) * density(40) 
  pd(21,21) = pd(21,21) - rrt(299) * density(21) * 4.0d0
  pd(29,21) = pd(29,21) + rrt(299) * density(21) * 2.0d0
  pd(32,21) = pd(32,21) + rrt(299) * density(21) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(300) * density(42) 
  pd(21,42) = pd(21,42) - rrt(300) * density(21) 
  pd(32,21) = pd(32,21) + rrt(300) * density(42) 
  pd(32,42) = pd(32,42) + rrt(300) * density(21) 
  pd(40,21) = pd(40,21) + rrt(300) * density(42) 
  pd(40,42) = pd(40,42) + rrt(300) * density(21) 
  pd(42,21) = pd(42,21) - rrt(300) * density(42) 
  pd(42,42) = pd(42,42) - rrt(300) * density(21) 
  pd(21,42) = pd(21,42) + rrt(301) * density(42) * 2.0d0
  pd(40,42) = pd(40,42) + rrt(301) * density(42) * 4.0d0
  pd(42,42) = pd(42,42) - rrt(301) * density(42) * 4.0d0
  pd(40,42) = pd(40,42) + rrt(302) * density(42) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(302) * density(42) * 4.0d0
  pd(43,42) = pd(43,42) + rrt(302) * density(42) * 2.0d0
  pd(21,32) = pd(21,32) + rrt(303) * density(42) 
  pd(21,42) = pd(21,42) + rrt(303) * density(32) 
  pd(32,32) = pd(32,32) - rrt(303) * density(42) 
  pd(32,42) = pd(32,42) - rrt(303) * density(32) 
  pd(42,32) = pd(42,32) - rrt(303) * density(42) 
  pd(42,42) = pd(42,42) - rrt(303) * density(32) 
  pd(43,32) = pd(43,32) + rrt(303) * density(42) 
  pd(43,42) = pd(43,42) + rrt(303) * density(32) 
  pd(21,42) = pd(21,42) + rrt(304) * density(43) 
  pd(21,43) = pd(21,43) + rrt(304) * density(42) 
  pd(40,42) = pd(40,42) + rrt(304) * density(43) 
  pd(40,43) = pd(40,43) + rrt(304) * density(42) 
  pd(43,42) = pd(43,42) - rrt(304) * density(43) 
  pd(43,43) = pd(43,43) - rrt(304) * density(42) 
  pd(21,21) = pd(21,21) - rrt(305) * density(43) 
  pd(21,43) = pd(21,43) - rrt(305) * density(21) 
  pd(32,21) = pd(32,21) + rrt(305) * density(43) 
  pd(32,43) = pd(32,43) + rrt(305) * density(21) 
  pd(42,21) = pd(42,21) + rrt(305) * density(43) 
  pd(42,43) = pd(42,43) + rrt(305) * density(21) 
  pd(43,21) = pd(43,21) - rrt(305) * density(43) 
  pd(43,43) = pd(43,43) - rrt(305) * density(21) 
  pd(21,43) = pd(21,43) + rrt(306) * density(43) * 2.0d0
  pd(42,43) = pd(42,43) + rrt(306) * density(43) * 4.0d0
  pd(43,43) = pd(43,43) - rrt(306) * density(43) * 4.0d0
  pd(14,14) = pd(14,14) - rrt(307) * density(14) * 4.0d0
  pd(18,14) = pd(18,14) + rrt(307) * density(14) * 2.0d0
  pd(53,14) = pd(53,14) + rrt(307) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(308) * density(29) 
  pd(14,29) = pd(14,29) - rrt(308) * density(14) 
  pd(29,14) = pd(29,14) - rrt(308) * density(29) 
  pd(29,29) = pd(29,29) - rrt(308) * density(14) 
  pd(45,14) = pd(45,14) + rrt(308) * density(29) 
  pd(45,29) = pd(45,29) + rrt(308) * density(14) 
  pd(53,14) = pd(53,14) + rrt(308) * density(29) 
  pd(53,29) = pd(53,29) + rrt(308) * density(14) 
  pd(01,01) = pd(01,01) - rrt(309) * density(01) * 2.0d0
  pd(14,01) = pd(14,01) + rrt(309) * density(01) * 4.0d0
  pd(01,01) = pd(01,01) - rrt(310) * density(21) 
  pd(01,21) = pd(01,21) - rrt(310) * density(01) 
  pd(14,01) = pd(14,01) + rrt(310) * density(21) * 2.0d0
  pd(14,21) = pd(14,21) + rrt(310) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(311) * density(40) 
  pd(01,40) = pd(01,40) - rrt(311) * density(01) 
  pd(14,01) = pd(14,01) + rrt(311) * density(40) * 2.0d0
  pd(14,40) = pd(14,40) + rrt(311) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(312) * density(29) 
  pd(01,29) = pd(01,29) - rrt(312) * density(01) 
  pd(14,01) = pd(14,01) + rrt(312) * density(29) * 2.0d0
  pd(14,29) = pd(14,29) + rrt(312) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(313) * density(14) 
  pd(01,14) = pd(01,14) - rrt(313) * density(01) 
  pd(14,01) = pd(14,01) + rrt(313) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) + rrt(313) * density(01) * 2.0d0
  pd(21,01) = pd(21,01) - rrt(314) * density(21) 
  pd(21,21) = pd(21,21) - rrt(314) * density(01) 
  pd(29,01) = pd(29,01) + rrt(314) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(314) * density(01) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(315) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(315) * density(21) * 4.0d0
  pd(21,21) = pd(21,21) - rrt(316) * density(29) 
  pd(21,29) = pd(21,29) - rrt(316) * density(21) 
  pd(29,21) = pd(29,21) + rrt(316) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) + rrt(316) * density(21) * 2.0d0
  pd(21,14) = pd(21,14) - rrt(317) * density(21) 
  pd(21,21) = pd(21,21) - rrt(317) * density(14) 
  pd(29,14) = pd(29,14) + rrt(317) * density(21) * 2.0d0
  pd(29,21) = pd(29,21) + rrt(317) * density(14) * 2.0d0
  pd(21,21) = pd(21,21) - rrt(318) * density(40) 
  pd(21,40) = pd(21,40) - rrt(318) * density(21) 
  pd(29,21) = pd(29,21) + rrt(318) * density(40) * 2.0d0
  pd(29,40) = pd(29,40) + rrt(318) * density(21) * 2.0d0
  pd(14,01) = pd(14,01) + rrt(319) * density(40) 
  pd(14,40) = pd(14,40) + rrt(319) * density(01) 
  pd(29,01) = pd(29,01) + rrt(319) * density(40) 
  pd(29,40) = pd(29,40) + rrt(319) * density(01) 
  pd(40,01) = pd(40,01) - rrt(319) * density(40) 
  pd(40,40) = pd(40,40) - rrt(319) * density(01) 
  pd(14,21) = pd(14,21) + rrt(320) * density(40) 
  pd(14,40) = pd(14,40) + rrt(320) * density(21) 
  pd(29,21) = pd(29,21) + rrt(320) * density(40) 
  pd(29,40) = pd(29,40) + rrt(320) * density(21) 
  pd(40,21) = pd(40,21) - rrt(320) * density(40) 
  pd(40,40) = pd(40,40) - rrt(320) * density(21) 
  pd(14,29) = pd(14,29) + rrt(321) * density(40) 
  pd(14,40) = pd(14,40) + rrt(321) * density(29) 
  pd(29,29) = pd(29,29) + rrt(321) * density(40) 
  pd(29,40) = pd(29,40) + rrt(321) * density(29) 
  pd(40,29) = pd(40,29) - rrt(321) * density(40) 
  pd(40,40) = pd(40,40) - rrt(321) * density(29) 
  pd(14,14) = pd(14,14) + rrt(322) * density(40) 
  pd(14,40) = pd(14,40) + rrt(322) * density(14) 
  pd(29,14) = pd(29,14) + rrt(322) * density(40) 
  pd(29,40) = pd(29,40) + rrt(322) * density(14) 
  pd(40,14) = pd(40,14) - rrt(322) * density(40) 
  pd(40,40) = pd(40,40) - rrt(322) * density(14) 
  pd(14,40) = pd(14,40) + rrt(323) * density(40) * 2.0d0
  pd(29,40) = pd(29,40) + rrt(323) * density(40) * 2.0d0
  pd(40,40) = pd(40,40) - rrt(323) * density(40) * 2.0d0
  pd(21,01) = pd(21,01) + rrt(324) * density(32) 
  pd(21,32) = pd(21,32) + rrt(324) * density(01) 
  pd(29,01) = pd(29,01) + rrt(324) * density(32) 
  pd(29,32) = pd(29,32) + rrt(324) * density(01) 
  pd(32,01) = pd(32,01) - rrt(324) * density(32) 
  pd(32,32) = pd(32,32) - rrt(324) * density(01) 
  pd(21,21) = pd(21,21) + rrt(325) * density(32) 
  pd(21,32) = pd(21,32) + rrt(325) * density(21) 
  pd(29,21) = pd(29,21) + rrt(325) * density(32) 
  pd(29,32) = pd(29,32) + rrt(325) * density(21) 
  pd(32,21) = pd(32,21) - rrt(325) * density(32) 
  pd(32,32) = pd(32,32) - rrt(325) * density(21) 
  pd(21,14) = pd(21,14) + rrt(326) * density(32) 
  pd(21,32) = pd(21,32) + rrt(326) * density(14) 
  pd(29,14) = pd(29,14) + rrt(326) * density(32) 
  pd(29,32) = pd(29,32) + rrt(326) * density(14) 
  pd(32,14) = pd(32,14) - rrt(326) * density(32) 
  pd(32,32) = pd(32,32) - rrt(326) * density(14) 
  pd(21,29) = pd(21,29) + rrt(327) * density(32) 
  pd(21,32) = pd(21,32) + rrt(327) * density(29) 
  pd(29,29) = pd(29,29) + rrt(327) * density(32) 
  pd(29,32) = pd(29,32) + rrt(327) * density(29) 
  pd(32,29) = pd(32,29) - rrt(327) * density(32) 
  pd(32,32) = pd(32,32) - rrt(327) * density(29) 
  pd(01,01) = pd(01,01) + rrt(328) * density(41) 
  pd(01,41) = pd(01,41) + rrt(328) * density(01) 
  pd(29,01) = pd(29,01) + rrt(328) * density(41) 
  pd(29,41) = pd(29,41) + rrt(328) * density(01) 
  pd(41,01) = pd(41,01) - rrt(328) * density(41) 
  pd(41,41) = pd(41,41) - rrt(328) * density(01) 
  pd(01,21) = pd(01,21) + rrt(329) * density(41) 
  pd(01,41) = pd(01,41) + rrt(329) * density(21) 
  pd(29,21) = pd(29,21) + rrt(329) * density(41) 
  pd(29,41) = pd(29,41) + rrt(329) * density(21) 
  pd(41,21) = pd(41,21) - rrt(329) * density(41) 
  pd(41,41) = pd(41,41) - rrt(329) * density(21) 
  pd(01,40) = pd(01,40) + rrt(330) * density(41) 
  pd(01,41) = pd(01,41) + rrt(330) * density(40) 
  pd(29,40) = pd(29,40) + rrt(330) * density(41) 
  pd(29,41) = pd(29,41) + rrt(330) * density(40) 
  pd(41,40) = pd(41,40) - rrt(330) * density(41) 
  pd(41,41) = pd(41,41) - rrt(330) * density(40) 
  pd(01,41) = pd(01,41) + rrt(331) * density(41) * 2.0d0
  pd(29,41) = pd(29,41) + rrt(331) * density(41) * 2.0d0
  pd(41,41) = pd(41,41) - rrt(331) * density(41) * 2.0d0
  pd(29,01) = pd(29,01) + rrt(332) * density(42) 
  pd(29,42) = pd(29,42) + rrt(332) * density(01) 
  pd(40,01) = pd(40,01) + rrt(332) * density(42) 
  pd(40,42) = pd(40,42) + rrt(332) * density(01) 
  pd(42,01) = pd(42,01) - rrt(332) * density(42) 
  pd(42,42) = pd(42,42) - rrt(332) * density(01) 
  pd(29,21) = pd(29,21) + rrt(333) * density(42) 
  pd(29,42) = pd(29,42) + rrt(333) * density(21) 
  pd(40,21) = pd(40,21) + rrt(333) * density(42) 
  pd(40,42) = pd(40,42) + rrt(333) * density(21) 
  pd(42,21) = pd(42,21) - rrt(333) * density(42) 
  pd(42,42) = pd(42,42) - rrt(333) * density(21) 
  pd(29,40) = pd(29,40) + rrt(334) * density(42) 
  pd(29,42) = pd(29,42) + rrt(334) * density(40) 
  pd(40,40) = pd(40,40) + rrt(334) * density(42) 
  pd(40,42) = pd(40,42) + rrt(334) * density(40) 
  pd(42,40) = pd(42,40) - rrt(334) * density(42) 
  pd(42,42) = pd(42,42) - rrt(334) * density(40) 
  pd(29,42) = pd(29,42) + rrt(335) * density(42) * 2.0d0
  pd(40,42) = pd(40,42) + rrt(335) * density(42) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(335) * density(42) * 2.0d0
  pd(29,01) = pd(29,01) + rrt(336) * density(43) 
  pd(29,43) = pd(29,43) + rrt(336) * density(01) 
  pd(42,01) = pd(42,01) + rrt(336) * density(43) 
  pd(42,43) = pd(42,43) + rrt(336) * density(01) 
  pd(43,01) = pd(43,01) - rrt(336) * density(43) 
  pd(43,43) = pd(43,43) - rrt(336) * density(01) 
  pd(29,21) = pd(29,21) + rrt(337) * density(43) 
  pd(29,43) = pd(29,43) + rrt(337) * density(21) 
  pd(42,21) = pd(42,21) + rrt(337) * density(43) 
  pd(42,43) = pd(42,43) + rrt(337) * density(21) 
  pd(43,21) = pd(43,21) - rrt(337) * density(43) 
  pd(43,43) = pd(43,43) - rrt(337) * density(21) 
  pd(29,40) = pd(29,40) + rrt(338) * density(43) 
  pd(29,43) = pd(29,43) + rrt(338) * density(40) 
  pd(42,40) = pd(42,40) + rrt(338) * density(43) 
  pd(42,43) = pd(42,43) + rrt(338) * density(40) 
  pd(43,40) = pd(43,40) - rrt(338) * density(43) 
  pd(43,43) = pd(43,43) - rrt(338) * density(40) 
  pd(29,14) = pd(29,14) + rrt(339) * density(43) 
  pd(29,43) = pd(29,43) + rrt(339) * density(14) 
  pd(42,14) = pd(42,14) + rrt(339) * density(43) 
  pd(42,43) = pd(42,43) + rrt(339) * density(14) 
  pd(43,14) = pd(43,14) - rrt(339) * density(43) 
  pd(43,43) = pd(43,43) - rrt(339) * density(14) 
  pd(29,29) = pd(29,29) + rrt(340) * density(43) 
  pd(29,43) = pd(29,43) + rrt(340) * density(29) 
  pd(42,29) = pd(42,29) + rrt(340) * density(43) 
  pd(42,43) = pd(42,43) + rrt(340) * density(29) 
  pd(43,29) = pd(43,29) - rrt(340) * density(43) 
  pd(43,43) = pd(43,43) - rrt(340) * density(29) 
  pd(21,01) = pd(21,01) + rrt(341) * density(43) 
  pd(21,43) = pd(21,43) + rrt(341) * density(01) 
  pd(40,01) = pd(40,01) + rrt(341) * density(43) 
  pd(40,43) = pd(40,43) + rrt(341) * density(01) 
  pd(43,01) = pd(43,01) - rrt(341) * density(43) 
  pd(43,43) = pd(43,43) - rrt(341) * density(01) 
  pd(21,21) = pd(21,21) + rrt(342) * density(43) 
  pd(21,43) = pd(21,43) + rrt(342) * density(21) 
  pd(40,21) = pd(40,21) + rrt(342) * density(43) 
  pd(40,43) = pd(40,43) + rrt(342) * density(21) 
  pd(43,21) = pd(43,21) - rrt(342) * density(43) 
  pd(43,43) = pd(43,43) - rrt(342) * density(21) 
  pd(21,40) = pd(21,40) + rrt(343) * density(43) 
  pd(21,43) = pd(21,43) + rrt(343) * density(40) 
  pd(40,40) = pd(40,40) + rrt(343) * density(43) 
  pd(40,43) = pd(40,43) + rrt(343) * density(40) 
  pd(43,40) = pd(43,40) - rrt(343) * density(43) 
  pd(43,43) = pd(43,43) - rrt(343) * density(40) 
  pd(21,14) = pd(21,14) + rrt(344) * density(43) 
  pd(21,43) = pd(21,43) + rrt(344) * density(14) 
  pd(40,14) = pd(40,14) + rrt(344) * density(43) 
  pd(40,43) = pd(40,43) + rrt(344) * density(14) 
  pd(43,14) = pd(43,14) - rrt(344) * density(43) 
  pd(43,43) = pd(43,43) - rrt(344) * density(14) 
  pd(21,29) = pd(21,29) + rrt(345) * density(43) 
  pd(21,43) = pd(21,43) + rrt(345) * density(29) 
  pd(40,29) = pd(40,29) + rrt(345) * density(43) 
  pd(40,43) = pd(40,43) + rrt(345) * density(29) 
  pd(43,29) = pd(43,29) - rrt(345) * density(43) 
  pd(43,43) = pd(43,43) - rrt(345) * density(29) 
  pd(42,44) = pd(42,44) + rrt(346) 
  pd(43,44) = pd(43,44) + rrt(346) 
  pd(44,44) = pd(44,44) - rrt(346) 
  pd(01,01) = pd(01,01) + rrt(347) * density(14)**2 
  pd(01,14) = pd(01,14) + rrt(347) * density(01) * density(14) * 2.0d0
  pd(14,01) = pd(14,01) - rrt(347) * density(14)**2 * 2.0d0
  pd(14,14) = pd(14,14) - rrt(347) * density(01) * density(14) * 4.0d0
  pd(01,14) = pd(01,14) + rrt(348) * density(14) * density(21) * 2.0d0
  pd(01,21) = pd(01,21) + rrt(348) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(348) * density(14) * density(21) * 4.0d0
  pd(14,21) = pd(14,21) - rrt(348) * density(14)**2 * 2.0d0
  pd(01,14) = pd(01,14) + rrt(349) * density(14) * density(40) * 2.0d0
  pd(01,40) = pd(01,40) + rrt(349) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(349) * density(14) * density(40) * 4.0d0
  pd(14,40) = pd(14,40) - rrt(349) * density(14)**2 * 2.0d0
  pd(01,14) = pd(01,14) + rrt(350) * density(14)**2 * 3.0d0
  pd(14,14) = pd(14,14) - rrt(350) * density(14)**2 * 6.0d0
  pd(01,14) = pd(01,14) + rrt(351) * density(14) * density(29) * 2.0d0
  pd(01,29) = pd(01,29) + rrt(351) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(351) * density(14) * density(29) * 4.0d0
  pd(14,29) = pd(14,29) - rrt(351) * density(14)**2 * 2.0d0
  pd(21,01) = pd(21,01) + rrt(352) * density(29)**2 
  pd(21,29) = pd(21,29) + rrt(352) * density(01) * density(29) * 2.0d0
  pd(29,01) = pd(29,01) - rrt(352) * density(29)**2 * 2.0d0
  pd(29,29) = pd(29,29) - rrt(352) * density(01) * density(29) * 4.0d0
  pd(21,21) = pd(21,21) + rrt(353) * density(29)**2 
  pd(21,29) = pd(21,29) + rrt(353) * density(21) * density(29) * 2.0d0
  pd(29,21) = pd(29,21) - rrt(353) * density(29)**2 * 2.0d0
  pd(29,29) = pd(29,29) - rrt(353) * density(21) * density(29) * 4.0d0
  pd(21,14) = pd(21,14) + rrt(354) * density(29)**2 
  pd(21,29) = pd(21,29) + rrt(354) * density(14) * density(29) * 2.0d0
  pd(29,14) = pd(29,14) - rrt(354) * density(29)**2 * 2.0d0
  pd(29,29) = pd(29,29) - rrt(354) * density(14) * density(29) * 4.0d0
  pd(21,29) = pd(21,29) + rrt(355) * density(29)**2 * 3.0d0
  pd(29,29) = pd(29,29) - rrt(355) * density(29)**2 * 6.0d0
  pd(21,29) = pd(21,29) + rrt(356) * density(29) * density(40) * 2.0d0
  pd(21,40) = pd(21,40) + rrt(356) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(356) * density(29) * density(40) * 4.0d0
  pd(29,40) = pd(29,40) - rrt(356) * density(29)**2 * 2.0d0
  pd(14,01) = pd(14,01) - rrt(357) * density(14) * density(29) 
  pd(14,14) = pd(14,14) - rrt(357) * density(01) * density(29) 
  pd(14,29) = pd(14,29) - rrt(357) * density(01) * density(14) 
  pd(29,01) = pd(29,01) - rrt(357) * density(14) * density(29) 
  pd(29,14) = pd(29,14) - rrt(357) * density(01) * density(29) 
  pd(29,29) = pd(29,29) - rrt(357) * density(01) * density(14) 
  pd(40,01) = pd(40,01) + rrt(357) * density(14) * density(29) 
  pd(40,14) = pd(40,14) + rrt(357) * density(01) * density(29) 
  pd(40,29) = pd(40,29) + rrt(357) * density(01) * density(14) 
  pd(14,14) = pd(14,14) - rrt(358) * density(21) * density(29) 
  pd(14,21) = pd(14,21) - rrt(358) * density(14) * density(29) 
  pd(14,29) = pd(14,29) - rrt(358) * density(14) * density(21) 
  pd(29,14) = pd(29,14) - rrt(358) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(358) * density(14) * density(29) 
  pd(29,29) = pd(29,29) - rrt(358) * density(14) * density(21) 
  pd(40,14) = pd(40,14) + rrt(358) * density(21) * density(29) 
  pd(40,21) = pd(40,21) + rrt(358) * density(14) * density(29) 
  pd(40,29) = pd(40,29) + rrt(358) * density(14) * density(21) 
  pd(14,14) = pd(14,14) - rrt(359) * density(14) * density(29) * 2.0d0
  pd(14,29) = pd(14,29) - rrt(359) * density(14)**2 
  pd(29,14) = pd(29,14) - rrt(359) * density(14) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(359) * density(14)**2 
  pd(40,14) = pd(40,14) + rrt(359) * density(14) * density(29) * 2.0d0
  pd(40,29) = pd(40,29) + rrt(359) * density(14)**2 
  pd(14,14) = pd(14,14) - rrt(360) * density(29)**2 
  pd(14,29) = pd(14,29) - rrt(360) * density(14) * density(29) * 2.0d0
  pd(29,14) = pd(29,14) - rrt(360) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(360) * density(14) * density(29) * 2.0d0
  pd(40,14) = pd(40,14) + rrt(360) * density(29)**2 
  pd(40,29) = pd(40,29) + rrt(360) * density(14) * density(29) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(361) * density(29) * density(40) 
  pd(14,29) = pd(14,29) - rrt(361) * density(14) * density(40) 
  pd(14,40) = pd(14,40) - rrt(361) * density(14) * density(29) 
  pd(29,14) = pd(29,14) - rrt(361) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(361) * density(14) * density(40) 
  pd(29,40) = pd(29,40) - rrt(361) * density(14) * density(29) 
  pd(40,14) = pd(40,14) + rrt(361) * density(29) * density(40) 
  pd(40,29) = pd(40,29) + rrt(361) * density(14) * density(40) 
  pd(40,40) = pd(40,40) + rrt(361) * density(14) * density(29) 
  pd(21,01) = pd(21,01) - rrt(362) * density(21) * density(29) 
  pd(21,21) = pd(21,21) - rrt(362) * density(01) * density(29) 
  pd(21,29) = pd(21,29) - rrt(362) * density(01) * density(21) 
  pd(29,01) = pd(29,01) - rrt(362) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(362) * density(01) * density(29) 
  pd(29,29) = pd(29,29) - rrt(362) * density(01) * density(21) 
  pd(32,01) = pd(32,01) + rrt(362) * density(21) * density(29) 
  pd(32,21) = pd(32,21) + rrt(362) * density(01) * density(29) 
  pd(32,29) = pd(32,29) + rrt(362) * density(01) * density(21) 
  pd(21,21) = pd(21,21) - rrt(363) * density(21) * density(29) * 2.0d0
  pd(21,29) = pd(21,29) - rrt(363) * density(21)**2 
  pd(29,21) = pd(29,21) - rrt(363) * density(21) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(363) * density(21)**2 
  pd(32,21) = pd(32,21) + rrt(363) * density(21) * density(29) * 2.0d0
  pd(32,29) = pd(32,29) + rrt(363) * density(21)**2 
  pd(21,21) = pd(21,21) - rrt(364) * density(29) * density(40) 
  pd(21,29) = pd(21,29) - rrt(364) * density(21) * density(40) 
  pd(21,40) = pd(21,40) - rrt(364) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(364) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(364) * density(21) * density(40) 
  pd(29,40) = pd(29,40) - rrt(364) * density(21) * density(29) 
  pd(32,21) = pd(32,21) + rrt(364) * density(29) * density(40) 
  pd(32,29) = pd(32,29) + rrt(364) * density(21) * density(40) 
  pd(32,40) = pd(32,40) + rrt(364) * density(21) * density(29) 
  pd(21,14) = pd(21,14) - rrt(365) * density(21) * density(29) 
  pd(21,21) = pd(21,21) - rrt(365) * density(14) * density(29) 
  pd(21,29) = pd(21,29) - rrt(365) * density(14) * density(21) 
  pd(29,14) = pd(29,14) - rrt(365) * density(21) * density(29) 
  pd(29,21) = pd(29,21) - rrt(365) * density(14) * density(29) 
  pd(29,29) = pd(29,29) - rrt(365) * density(14) * density(21) 
  pd(32,14) = pd(32,14) + rrt(365) * density(21) * density(29) 
  pd(32,21) = pd(32,21) + rrt(365) * density(14) * density(29) 
  pd(32,29) = pd(32,29) + rrt(365) * density(14) * density(21) 
  pd(21,21) = pd(21,21) - rrt(366) * density(29)**2 
  pd(21,29) = pd(21,29) - rrt(366) * density(21) * density(29) * 2.0d0
  pd(29,21) = pd(29,21) - rrt(366) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(366) * density(21) * density(29) * 2.0d0
  pd(32,21) = pd(32,21) + rrt(366) * density(29)**2 
  pd(32,29) = pd(32,29) + rrt(366) * density(21) * density(29) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(367) * density(29) 
  pd(01,29) = pd(01,29) - rrt(367) * density(01) 
  pd(29,01) = pd(29,01) - rrt(367) * density(29) 
  pd(29,29) = pd(29,29) - rrt(367) * density(01) 
  pd(41,01) = pd(41,01) + rrt(367) * density(29) 
  pd(41,29) = pd(41,29) + rrt(367) * density(01) 
  pd(29,01) = pd(29,01) - rrt(368) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(368) * density(01) * density(40) 
  pd(29,40) = pd(29,40) - rrt(368) * density(01) * density(29) 
  pd(40,01) = pd(40,01) - rrt(368) * density(29) * density(40) 
  pd(40,29) = pd(40,29) - rrt(368) * density(01) * density(40) 
  pd(40,40) = pd(40,40) - rrt(368) * density(01) * density(29) 
  pd(42,01) = pd(42,01) + rrt(368) * density(29) * density(40) 
  pd(42,29) = pd(42,29) + rrt(368) * density(01) * density(40) 
  pd(42,40) = pd(42,40) + rrt(368) * density(01) * density(29) 
  pd(29,21) = pd(29,21) - rrt(369) * density(29) * density(40) 
  pd(29,29) = pd(29,29) - rrt(369) * density(21) * density(40) 
  pd(29,40) = pd(29,40) - rrt(369) * density(21) * density(29) 
  pd(40,21) = pd(40,21) - rrt(369) * density(29) * density(40) 
  pd(40,29) = pd(40,29) - rrt(369) * density(21) * density(40) 
  pd(40,40) = pd(40,40) - rrt(369) * density(21) * density(29) 
  pd(42,21) = pd(42,21) + rrt(369) * density(29) * density(40) 
  pd(42,29) = pd(42,29) + rrt(369) * density(21) * density(40) 
  pd(42,40) = pd(42,40) + rrt(369) * density(21) * density(29) 
  pd(29,29) = pd(29,29) - rrt(370) * density(40)**2 
  pd(29,40) = pd(29,40) - rrt(370) * density(29) * density(40) * 2.0d0
  pd(40,29) = pd(40,29) - rrt(370) * density(40)**2 
  pd(40,40) = pd(40,40) - rrt(370) * density(29) * density(40) * 2.0d0
  pd(42,29) = pd(42,29) + rrt(370) * density(40)**2 
  pd(42,40) = pd(42,40) + rrt(370) * density(29) * density(40) * 2.0d0
  pd(29,01) = pd(29,01) - rrt(371) * density(29) * density(42) 
  pd(29,29) = pd(29,29) - rrt(371) * density(01) * density(42) 
  pd(29,42) = pd(29,42) - rrt(371) * density(01) * density(29) 
  pd(42,01) = pd(42,01) - rrt(371) * density(29) * density(42) 
  pd(42,29) = pd(42,29) - rrt(371) * density(01) * density(42) 
  pd(42,42) = pd(42,42) - rrt(371) * density(01) * density(29) 
  pd(43,01) = pd(43,01) + rrt(371) * density(29) * density(42) 
  pd(43,29) = pd(43,29) + rrt(371) * density(01) * density(42) 
  pd(43,42) = pd(43,42) + rrt(371) * density(01) * density(29) 
  pd(29,21) = pd(29,21) - rrt(372) * density(29) * density(42) 
  pd(29,29) = pd(29,29) - rrt(372) * density(21) * density(42) 
  pd(29,42) = pd(29,42) - rrt(372) * density(21) * density(29) 
  pd(42,21) = pd(42,21) - rrt(372) * density(29) * density(42) 
  pd(42,29) = pd(42,29) - rrt(372) * density(21) * density(42) 
  pd(42,42) = pd(42,42) - rrt(372) * density(21) * density(29) 
  pd(43,21) = pd(43,21) + rrt(372) * density(29) * density(42) 
  pd(43,29) = pd(43,29) + rrt(372) * density(21) * density(42) 
  pd(43,42) = pd(43,42) + rrt(372) * density(21) * density(29) 
  pd(29,14) = pd(29,14) - rrt(373) * density(29) * density(42) 
  pd(29,29) = pd(29,29) - rrt(373) * density(14) * density(42) 
  pd(29,42) = pd(29,42) - rrt(373) * density(14) * density(29) 
  pd(42,14) = pd(42,14) - rrt(373) * density(29) * density(42) 
  pd(42,29) = pd(42,29) - rrt(373) * density(14) * density(42) 
  pd(42,42) = pd(42,42) - rrt(373) * density(14) * density(29) 
  pd(43,14) = pd(43,14) + rrt(373) * density(29) * density(42) 
  pd(43,29) = pd(43,29) + rrt(373) * density(14) * density(42) 
  pd(43,42) = pd(43,42) + rrt(373) * density(14) * density(29) 
  pd(29,29) = pd(29,29) - rrt(374) * density(29) * density(42) * 2.0d0
  pd(29,42) = pd(29,42) - rrt(374) * density(29)**2 
  pd(42,29) = pd(42,29) - rrt(374) * density(29) * density(42) * 2.0d0
  pd(42,42) = pd(42,42) - rrt(374) * density(29)**2 
  pd(43,29) = pd(43,29) + rrt(374) * density(29) * density(42) * 2.0d0
  pd(43,42) = pd(43,42) + rrt(374) * density(29)**2 
  pd(29,29) = pd(29,29) - rrt(375) * density(40) * density(42) 
  pd(29,40) = pd(29,40) - rrt(375) * density(29) * density(42) 
  pd(29,42) = pd(29,42) - rrt(375) * density(29) * density(40) 
  pd(42,29) = pd(42,29) - rrt(375) * density(40) * density(42) 
  pd(42,40) = pd(42,40) - rrt(375) * density(29) * density(42) 
  pd(42,42) = pd(42,42) - rrt(375) * density(29) * density(40) 
  pd(43,29) = pd(43,29) + rrt(375) * density(40) * density(42) 
  pd(43,40) = pd(43,40) + rrt(375) * density(29) * density(42) 
  pd(43,42) = pd(43,42) + rrt(375) * density(29) * density(40) 
  pd(42,42) = pd(42,42) - rrt(376) * density(43) 
  pd(42,43) = pd(42,43) - rrt(376) * density(42) 
  pd(43,42) = pd(43,42) - rrt(376) * density(43) 
  pd(43,43) = pd(43,43) - rrt(376) * density(42) 
  pd(44,42) = pd(44,42) + rrt(376) * density(43) 
  pd(44,43) = pd(44,43) + rrt(376) * density(42) 
  pd(14,17) = pd(14,17) + rrt(377) * density(29) 
  pd(14,29) = pd(14,29) + rrt(377) * density(17) 
  pd(17,17) = pd(17,17) - rrt(377) * density(29) 
  pd(17,29) = pd(17,29) - rrt(377) * density(17) 
  pd(29,17) = pd(29,17) - rrt(377) * density(29) 
  pd(29,29) = pd(29,29) - rrt(377) * density(17) 
  pd(33,17) = pd(33,17) + rrt(377) * density(29) 
  pd(33,29) = pd(33,29) + rrt(377) * density(17) 
  pd(14,17) = pd(14,17) + rrt(378) * density(21) 
  pd(14,21) = pd(14,21) + rrt(378) * density(17) 
  pd(17,17) = pd(17,17) - rrt(378) * density(21) 
  pd(17,21) = pd(17,21) - rrt(378) * density(17) 
  pd(21,17) = pd(21,17) - rrt(378) * density(21) 
  pd(21,21) = pd(21,21) - rrt(378) * density(17) 
  pd(34,17) = pd(34,17) + rrt(378) * density(21) 
  pd(34,21) = pd(34,21) + rrt(378) * density(17) 
  pd(17,17) = pd(17,17) - rrt(379) * density(21) 
  pd(17,21) = pd(17,21) - rrt(379) * density(17) 
  pd(21,17) = pd(21,17) - rrt(379) * density(21) 
  pd(21,21) = pd(21,21) - rrt(379) * density(17) 
  pd(29,17) = pd(29,17) + rrt(379) * density(21) 
  pd(29,21) = pd(29,21) + rrt(379) * density(17) 
  pd(45,17) = pd(45,17) + rrt(379) * density(21) 
  pd(45,21) = pd(45,21) + rrt(379) * density(17) 
  pd(17,17) = pd(17,17) - rrt(380) * density(21) 
  pd(17,21) = pd(17,21) - rrt(380) * density(17) 
  pd(21,17) = pd(21,17) - rrt(380) * density(21) 
  pd(21,21) = pd(21,21) - rrt(380) * density(17) 
  pd(33,17) = pd(33,17) + rrt(380) * density(21) 
  pd(33,21) = pd(33,21) + rrt(380) * density(17) 
  pd(40,17) = pd(40,17) + rrt(380) * density(21) 
  pd(40,21) = pd(40,21) + rrt(380) * density(17) 
  pd(17,17) = pd(17,17) - rrt(381) * density(32) 
  pd(17,32) = pd(17,32) - rrt(381) * density(17) 
  pd(21,17) = pd(21,17) + rrt(381) * density(32) 
  pd(21,32) = pd(21,32) + rrt(381) * density(17) 
  pd(32,17) = pd(32,17) - rrt(381) * density(32) 
  pd(32,32) = pd(32,32) - rrt(381) * density(17) 
  pd(45,17) = pd(45,17) + rrt(381) * density(32) 
  pd(45,32) = pd(45,32) + rrt(381) * density(17) 
  pd(14,17) = pd(14,17) + rrt(382) * density(40) 
  pd(14,40) = pd(14,40) + rrt(382) * density(17) 
  pd(17,17) = pd(17,17) - rrt(382) * density(40) 
  pd(17,40) = pd(17,40) - rrt(382) * density(17) 
  pd(40,17) = pd(40,17) - rrt(382) * density(40) 
  pd(40,40) = pd(40,40) - rrt(382) * density(17) 
  pd(45,17) = pd(45,17) + rrt(382) * density(40) 
  pd(45,40) = pd(45,40) + rrt(382) * density(17) 
  pd(17,17) = pd(17,17) - rrt(383) * density(40) 
  pd(17,40) = pd(17,40) - rrt(383) * density(17) 
  pd(18,17) = pd(18,17) + rrt(383) * density(40) 
  pd(18,40) = pd(18,40) + rrt(383) * density(17) 
  pd(29,17) = pd(29,17) + rrt(383) * density(40) 
  pd(29,40) = pd(29,40) + rrt(383) * density(17) 
  pd(40,17) = pd(40,17) - rrt(383) * density(40) 
  pd(40,40) = pd(40,40) - rrt(383) * density(17) 
  pd(01,17) = pd(01,17) + rrt(384) * density(40) 
  pd(01,40) = pd(01,40) + rrt(384) * density(17) 
  pd(17,17) = pd(17,17) - rrt(384) * density(40) 
  pd(17,40) = pd(17,40) - rrt(384) * density(17) 
  pd(33,17) = pd(33,17) + rrt(384) * density(40) 
  pd(33,40) = pd(33,40) + rrt(384) * density(17) 
  pd(40,17) = pd(40,17) - rrt(384) * density(40) 
  pd(40,40) = pd(40,40) - rrt(384) * density(17) 
  pd(01,17) = pd(01,17) + rrt(385) * density(41) 
  pd(01,41) = pd(01,41) + rrt(385) * density(17) 
  pd(17,17) = pd(17,17) - rrt(385) * density(41) 
  pd(17,41) = pd(17,41) - rrt(385) * density(17) 
  pd(41,17) = pd(41,17) - rrt(385) * density(41) 
  pd(41,41) = pd(41,41) - rrt(385) * density(17) 
  pd(45,17) = pd(45,17) + rrt(385) * density(41) 
  pd(45,41) = pd(45,41) + rrt(385) * density(17) 
  pd(01,01) = pd(01,01) - rrt(386) * density(33) 
  pd(01,33) = pd(01,33) - rrt(386) * density(01) 
  pd(14,01) = pd(14,01) + rrt(386) * density(33) 
  pd(14,33) = pd(14,33) + rrt(386) * density(01) 
  pd(33,01) = pd(33,01) - rrt(386) * density(33) 
  pd(33,33) = pd(33,33) - rrt(386) * density(01) 
  pd(45,01) = pd(45,01) + rrt(386) * density(33) 
  pd(45,33) = pd(45,33) + rrt(386) * density(01) 
  pd(21,21) = pd(21,21) - rrt(387) * density(33) 
  pd(21,33) = pd(21,33) - rrt(387) * density(21) 
  pd(29,21) = pd(29,21) + rrt(387) * density(33) 
  pd(29,33) = pd(29,33) + rrt(387) * density(21) 
  pd(33,21) = pd(33,21) - rrt(387) * density(33) 
  pd(33,33) = pd(33,33) - rrt(387) * density(21) 
  pd(34,21) = pd(34,21) + rrt(387) * density(33) 
  pd(34,33) = pd(34,33) + rrt(387) * density(21) 
  pd(21,32) = pd(21,32) + rrt(388) * density(33) 
  pd(21,33) = pd(21,33) + rrt(388) * density(32) 
  pd(32,32) = pd(32,32) - rrt(388) * density(33) 
  pd(32,33) = pd(32,33) - rrt(388) * density(32) 
  pd(33,32) = pd(33,32) - rrt(388) * density(33) 
  pd(33,33) = pd(33,33) - rrt(388) * density(32) 
  pd(34,32) = pd(34,32) + rrt(388) * density(33) 
  pd(34,33) = pd(34,33) + rrt(388) * density(32) 
  pd(29,33) = pd(29,33) + rrt(389) * density(40) 
  pd(29,40) = pd(29,40) + rrt(389) * density(33) 
  pd(33,33) = pd(33,33) - rrt(389) * density(40) 
  pd(33,40) = pd(33,40) - rrt(389) * density(33) 
  pd(40,33) = pd(40,33) - rrt(389) * density(40) 
  pd(40,40) = pd(40,40) - rrt(389) * density(33) 
  pd(45,33) = pd(45,33) + rrt(389) * density(40) 
  pd(45,40) = pd(45,40) + rrt(389) * density(33) 
  pd(14,33) = pd(14,33) + rrt(390) * density(40) 
  pd(14,40) = pd(14,40) + rrt(390) * density(33) 
  pd(33,33) = pd(33,33) - rrt(390) * density(40) 
  pd(33,40) = pd(33,40) - rrt(390) * density(33) 
  pd(34,33) = pd(34,33) + rrt(390) * density(40) 
  pd(34,40) = pd(34,40) + rrt(390) * density(33) 
  pd(40,33) = pd(40,33) - rrt(390) * density(40) 
  pd(40,40) = pd(40,40) - rrt(390) * density(33) 
  pd(15,15) = pd(15,15) - rrt(391) * density(33) 
  pd(15,33) = pd(15,33) - rrt(391) * density(15) 
  pd(17,15) = pd(17,15) + rrt(391) * density(33) 
  pd(17,33) = pd(17,33) + rrt(391) * density(15) 
  pd(29,15) = pd(29,15) + rrt(391) * density(33) 
  pd(29,33) = pd(29,33) + rrt(391) * density(15) 
  pd(33,15) = pd(33,15) - rrt(391) * density(33) 
  pd(33,33) = pd(33,33) - rrt(391) * density(15) 
  pd(33,33) = pd(33,33) - rrt(392) * density(41) 
  pd(33,41) = pd(33,41) - rrt(392) * density(33) 
  pd(40,33) = pd(40,33) + rrt(392) * density(41) 
  pd(40,41) = pd(40,41) + rrt(392) * density(33) 
  pd(41,33) = pd(41,33) - rrt(392) * density(41) 
  pd(41,41) = pd(41,41) - rrt(392) * density(33) 
  pd(45,33) = pd(45,33) + rrt(392) * density(41) 
  pd(45,41) = pd(45,41) + rrt(392) * density(33) 
  pd(29,33) = pd(29,33) + rrt(393) * density(41) 
  pd(29,41) = pd(29,41) + rrt(393) * density(33) 
  pd(33,33) = pd(33,33) - rrt(393) * density(41) 
  pd(33,41) = pd(33,41) - rrt(393) * density(33) 
  pd(41,33) = pd(41,33) - rrt(393) * density(41) 
  pd(41,41) = pd(41,41) - rrt(393) * density(33) 
  pd(46,33) = pd(46,33) + rrt(393) * density(41) 
  pd(46,41) = pd(46,41) + rrt(393) * density(33) 
  pd(01,33) = pd(01,33) + rrt(394) * density(41) 
  pd(01,41) = pd(01,41) + rrt(394) * density(33) 
  pd(33,33) = pd(33,33) - rrt(394) * density(41) 
  pd(33,41) = pd(33,41) - rrt(394) * density(33) 
  pd(34,33) = pd(34,33) + rrt(394) * density(41) 
  pd(34,41) = pd(34,41) + rrt(394) * density(33) 
  pd(41,33) = pd(41,33) - rrt(394) * density(41) 
  pd(41,41) = pd(41,41) - rrt(394) * density(33) 
  pd(29,33) = pd(29,33) + rrt(395) * density(42) 
  pd(29,42) = pd(29,42) + rrt(395) * density(33) 
  pd(33,33) = pd(33,33) - rrt(395) * density(42) 
  pd(33,42) = pd(33,42) - rrt(395) * density(33) 
  pd(42,33) = pd(42,33) - rrt(395) * density(42) 
  pd(42,42) = pd(42,42) - rrt(395) * density(33) 
  pd(47,33) = pd(47,33) + rrt(395) * density(42) 
  pd(47,42) = pd(47,42) + rrt(395) * density(33) 
  pd(01,18) = pd(01,18) + rrt(396) * density(21) 
  pd(01,21) = pd(01,21) + rrt(396) * density(18) 
  pd(18,18) = pd(18,18) - rrt(396) * density(21) 
  pd(18,21) = pd(18,21) - rrt(396) * density(18) 
  pd(21,18) = pd(21,18) - rrt(396) * density(21) 
  pd(21,21) = pd(21,21) - rrt(396) * density(18) 
  pd(34,18) = pd(34,18) + rrt(396) * density(21) 
  pd(34,21) = pd(34,21) + rrt(396) * density(18) 
  pd(14,18) = pd(14,18) + rrt(397) * density(29) 
  pd(14,29) = pd(14,29) + rrt(397) * density(18) 
  pd(18,18) = pd(18,18) - rrt(397) * density(29) 
  pd(18,29) = pd(18,29) - rrt(397) * density(18) 
  pd(29,18) = pd(29,18) - rrt(397) * density(29) 
  pd(29,29) = pd(29,29) - rrt(397) * density(18) 
  pd(45,18) = pd(45,18) + rrt(397) * density(29) 
  pd(45,29) = pd(45,29) + rrt(397) * density(18) 
  pd(01,18) = pd(01,18) + rrt(398) * density(32) 
  pd(01,32) = pd(01,32) + rrt(398) * density(18) 
  pd(18,18) = pd(18,18) - rrt(398) * density(32) 
  pd(18,32) = pd(18,32) - rrt(398) * density(18) 
  pd(29,18) = pd(29,18) + rrt(398) * density(32) 
  pd(29,32) = pd(29,32) + rrt(398) * density(18) 
  pd(32,18) = pd(32,18) - rrt(398) * density(32) 
  pd(32,32) = pd(32,32) - rrt(398) * density(18) 
  pd(34,18) = pd(34,18) + rrt(398) * density(32) 
  pd(34,32) = pd(34,32) + rrt(398) * density(18) 
  pd(01,14) = pd(01,14) + rrt(399) * density(18) 
  pd(01,18) = pd(01,18) + rrt(399) * density(14) 
  pd(14,14) = pd(14,14) - rrt(399) * density(18) 
  pd(14,18) = pd(14,18) - rrt(399) * density(14) 
  pd(17,14) = pd(17,14) + rrt(399) * density(18) 
  pd(17,18) = pd(17,18) + rrt(399) * density(14) 
  pd(18,14) = pd(18,14) - rrt(399) * density(18) 
  pd(18,18) = pd(18,18) - rrt(399) * density(14) 
  pd(01,18) = pd(01,18) + rrt(400) * density(40) 
  pd(01,40) = pd(01,40) + rrt(400) * density(18) 
  pd(18,18) = pd(18,18) - rrt(400) * density(40) 
  pd(18,40) = pd(18,40) - rrt(400) * density(18) 
  pd(40,18) = pd(40,18) - rrt(400) * density(40) 
  pd(40,40) = pd(40,40) - rrt(400) * density(18) 
  pd(45,18) = pd(45,18) + rrt(400) * density(40) 
  pd(45,40) = pd(45,40) + rrt(400) * density(18) 
  pd(01,18) = pd(01,18) + rrt(401) * density(41) 
  pd(01,41) = pd(01,41) + rrt(401) * density(18) 
  pd(18,18) = pd(18,18) - rrt(401) * density(41) 
  pd(18,41) = pd(18,41) - rrt(401) * density(18) 
  pd(41,18) = pd(41,18) - rrt(401) * density(41) 
  pd(41,41) = pd(41,41) - rrt(401) * density(18) 
  pd(46,18) = pd(46,18) + rrt(401) * density(41) 
  pd(46,41) = pd(46,41) + rrt(401) * density(18) 
  pd(01,18) = pd(01,18) + rrt(402) * density(41) 
  pd(01,41) = pd(01,41) + rrt(402) * density(18) 
  pd(14,18) = pd(14,18) + rrt(402) * density(41) 
  pd(14,41) = pd(14,41) + rrt(402) * density(18) 
  pd(18,18) = pd(18,18) - rrt(402) * density(41) 
  pd(18,41) = pd(18,41) - rrt(402) * density(18) 
  pd(41,18) = pd(41,18) - rrt(402) * density(41) 
  pd(41,41) = pd(41,41) - rrt(402) * density(18) 
  pd(45,18) = pd(45,18) + rrt(402) * density(41) 
  pd(45,41) = pd(45,41) + rrt(402) * density(18) 
  pd(01,01) = pd(01,01) - rrt(403) * density(34) 
  pd(01,34) = pd(01,34) - rrt(403) * density(01) 
  pd(34,01) = pd(34,01) - rrt(403) * density(34) 
  pd(34,34) = pd(34,34) - rrt(403) * density(01) 
  pd(40,01) = pd(40,01) + rrt(403) * density(34) 
  pd(40,34) = pd(40,34) + rrt(403) * density(01) 
  pd(45,01) = pd(45,01) + rrt(403) * density(34) 
  pd(45,34) = pd(45,34) + rrt(403) * density(01) 
  pd(14,14) = pd(14,14) - rrt(404) * density(34) 
  pd(14,34) = pd(14,34) - rrt(404) * density(14) 
  pd(29,14) = pd(29,14) + rrt(404) * density(34) 
  pd(29,34) = pd(29,34) + rrt(404) * density(14) 
  pd(34,14) = pd(34,14) - rrt(404) * density(34) 
  pd(34,34) = pd(34,34) - rrt(404) * density(14) 
  pd(45,14) = pd(45,14) + rrt(404) * density(34) 
  pd(45,34) = pd(45,34) + rrt(404) * density(14) 
  pd(21,34) = pd(21,34) + rrt(405) * density(40) 
  pd(21,40) = pd(21,40) + rrt(405) * density(34) 
  pd(34,34) = pd(34,34) - rrt(405) * density(40) 
  pd(34,40) = pd(34,40) - rrt(405) * density(34) 
  pd(40,34) = pd(40,34) - rrt(405) * density(40) 
  pd(40,40) = pd(40,40) - rrt(405) * density(34) 
  pd(45,34) = pd(45,34) + rrt(405) * density(40) 
  pd(45,40) = pd(45,40) + rrt(405) * density(34) 
  pd(32,34) = pd(32,34) + rrt(406) * density(42) 
  pd(32,42) = pd(32,42) + rrt(406) * density(34) 
  pd(34,34) = pd(34,34) - rrt(406) * density(42) 
  pd(34,42) = pd(34,42) - rrt(406) * density(34) 
  pd(42,34) = pd(42,34) - rrt(406) * density(42) 
  pd(42,42) = pd(42,42) - rrt(406) * density(34) 
  pd(45,34) = pd(45,34) + rrt(406) * density(42) 
  pd(45,42) = pd(45,42) + rrt(406) * density(34) 
  pd(21,34) = pd(21,34) + rrt(407) * density(42) 
  pd(21,42) = pd(21,42) + rrt(407) * density(34) 
  pd(34,34) = pd(34,34) - rrt(407) * density(42) 
  pd(34,42) = pd(34,42) - rrt(407) * density(34) 
  pd(42,34) = pd(42,34) - rrt(407) * density(42) 
  pd(42,42) = pd(42,42) - rrt(407) * density(34) 
  pd(47,34) = pd(47,34) + rrt(407) * density(42) 
  pd(47,42) = pd(47,42) + rrt(407) * density(34) 
  pd(01,19) = pd(01,19) + rrt(408) * density(21) 
  pd(01,21) = pd(01,21) + rrt(408) * density(19) 
  pd(14,19) = pd(14,19) + rrt(408) * density(21) 
  pd(14,21) = pd(14,21) + rrt(408) * density(19) 
  pd(19,19) = pd(19,19) - rrt(408) * density(21) 
  pd(19,21) = pd(19,21) - rrt(408) * density(19) 
  pd(21,19) = pd(21,19) - rrt(408) * density(21) 
  pd(21,21) = pd(21,21) - rrt(408) * density(19) 
  pd(34,19) = pd(34,19) + rrt(408) * density(21) 
  pd(34,21) = pd(34,21) + rrt(408) * density(19) 
  pd(01,19) = pd(01,19) + rrt(409) * density(21) 
  pd(01,21) = pd(01,21) + rrt(409) * density(19) 
  pd(19,19) = pd(19,19) - rrt(409) * density(21) 
  pd(19,21) = pd(19,21) - rrt(409) * density(19) 
  pd(21,19) = pd(21,19) - rrt(409) * density(21) 
  pd(21,21) = pd(21,21) - rrt(409) * density(19) 
  pd(47,19) = pd(47,19) + rrt(409) * density(21) 
  pd(47,21) = pd(47,21) + rrt(409) * density(19) 
  pd(01,14) = pd(01,14) + rrt(410) * density(19) 
  pd(01,19) = pd(01,19) + rrt(410) * density(14) 
  pd(14,14) = pd(14,14) - rrt(410) * density(19) 
  pd(14,19) = pd(14,19) - rrt(410) * density(14) 
  pd(18,14) = pd(18,14) + rrt(410) * density(19) 
  pd(18,19) = pd(18,19) + rrt(410) * density(14) 
  pd(19,14) = pd(19,14) - rrt(410) * density(19) 
  pd(19,19) = pd(19,19) - rrt(410) * density(14) 
  pd(01,19) = pd(01,19) + rrt(411) * density(40) 
  pd(01,40) = pd(01,40) + rrt(411) * density(19) 
  pd(14,19) = pd(14,19) + rrt(411) * density(40) 
  pd(14,40) = pd(14,40) + rrt(411) * density(19) 
  pd(19,19) = pd(19,19) - rrt(411) * density(40) 
  pd(19,40) = pd(19,40) - rrt(411) * density(19) 
  pd(40,19) = pd(40,19) - rrt(411) * density(40) 
  pd(40,40) = pd(40,40) - rrt(411) * density(19) 
  pd(45,19) = pd(45,19) + rrt(411) * density(40) 
  pd(45,40) = pd(45,40) + rrt(411) * density(19) 
  pd(01,19) = pd(01,19) + rrt(412) * density(40) 
  pd(01,40) = pd(01,40) + rrt(412) * density(19) 
  pd(19,19) = pd(19,19) - rrt(412) * density(40) 
  pd(19,40) = pd(19,40) - rrt(412) * density(19) 
  pd(40,19) = pd(40,19) - rrt(412) * density(40) 
  pd(40,40) = pd(40,40) - rrt(412) * density(19) 
  pd(46,19) = pd(46,19) + rrt(412) * density(40) 
  pd(46,40) = pd(46,40) + rrt(412) * density(19) 
  pd(40,40) = pd(40,40) - rrt(413) * density(47) 
  pd(40,47) = pd(40,47) - rrt(413) * density(40) 
  pd(42,40) = pd(42,40) + rrt(413) * density(47) 
  pd(42,47) = pd(42,47) + rrt(413) * density(40) 
  pd(45,40) = pd(45,40) + rrt(413) * density(47) 
  pd(45,47) = pd(45,47) + rrt(413) * density(40) 
  pd(47,40) = pd(47,40) - rrt(413) * density(47) 
  pd(47,47) = pd(47,47) - rrt(413) * density(40) 
  pd(40,40) = pd(40,40) - rrt(414) * density(46) 
  pd(40,46) = pd(40,46) - rrt(414) * density(40) 
  pd(41,40) = pd(41,40) + rrt(414) * density(46) 
  pd(41,46) = pd(41,46) + rrt(414) * density(40) 
  pd(45,40) = pd(45,40) + rrt(414) * density(46) 
  pd(45,46) = pd(45,46) + rrt(414) * density(40) 
  pd(46,40) = pd(46,40) - rrt(414) * density(46) 
  pd(46,46) = pd(46,46) - rrt(414) * density(40) 
  pd(01,01) = pd(01,01) + rrt(415) * density(20) 
  pd(01,20) = pd(01,20) + rrt(415) * density(01) 
  pd(18,01) = pd(18,01) + rrt(415) * density(20) 
  pd(18,20) = pd(18,20) + rrt(415) * density(01) 
  pd(20,01) = pd(20,01) - rrt(415) * density(20) 
  pd(20,20) = pd(20,20) - rrt(415) * density(01) 
  pd(01,20) = pd(01,20) + rrt(416) * density(21) * 2.0d0
  pd(01,21) = pd(01,21) + rrt(416) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(416) * density(21) 
  pd(20,21) = pd(20,21) - rrt(416) * density(20) 
  pd(21,20) = pd(21,20) - rrt(416) * density(21) 
  pd(21,21) = pd(21,21) - rrt(416) * density(20) 
  pd(34,20) = pd(34,20) + rrt(416) * density(21) 
  pd(34,21) = pd(34,21) + rrt(416) * density(20) 
  pd(01,20) = pd(01,20) + rrt(417) * density(29) * 2.0d0
  pd(01,29) = pd(01,29) + rrt(417) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(417) * density(29) 
  pd(20,29) = pd(20,29) - rrt(417) * density(20) 
  pd(29,20) = pd(29,20) - rrt(417) * density(29) 
  pd(29,29) = pd(29,29) - rrt(417) * density(20) 
  pd(33,20) = pd(33,20) + rrt(417) * density(29) 
  pd(33,29) = pd(33,29) + rrt(417) * density(20) 
  pd(01,14) = pd(01,14) + rrt(418) * density(20) * 2.0d0
  pd(01,20) = pd(01,20) + rrt(418) * density(14) * 2.0d0
  pd(14,14) = pd(14,14) - rrt(418) * density(20) 
  pd(14,20) = pd(14,20) - rrt(418) * density(14) 
  pd(17,14) = pd(17,14) + rrt(418) * density(20) 
  pd(17,20) = pd(17,20) + rrt(418) * density(14) 
  pd(20,14) = pd(20,14) - rrt(418) * density(20) 
  pd(20,20) = pd(20,20) - rrt(418) * density(14) 
  pd(01,20) = pd(01,20) + rrt(419) * density(40) * 2.0d0
  pd(01,40) = pd(01,40) + rrt(419) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(419) * density(40) 
  pd(20,40) = pd(20,40) - rrt(419) * density(20) 
  pd(40,20) = pd(40,20) - rrt(419) * density(40) 
  pd(40,40) = pd(40,40) - rrt(419) * density(20) 
  pd(45,20) = pd(45,20) + rrt(419) * density(40) 
  pd(45,40) = pd(45,40) + rrt(419) * density(20) 
  pd(01,01) = pd(01,01) - rrt(420) * density(35) 
  pd(01,35) = pd(01,35) - rrt(420) * density(01) 
  pd(21,01) = pd(21,01) + rrt(420) * density(35) 
  pd(21,35) = pd(21,35) + rrt(420) * density(01) 
  pd(35,01) = pd(35,01) - rrt(420) * density(35) 
  pd(35,35) = pd(35,35) - rrt(420) * density(01) 
  pd(52,01) = pd(52,01) + rrt(420) * density(35) 
  pd(52,35) = pd(52,35) + rrt(420) * density(01) 
  pd(21,21) = pd(21,21) + rrt(421) * density(35) 
  pd(21,35) = pd(21,35) + rrt(421) * density(21) 
  pd(34,21) = pd(34,21) + rrt(421) * density(35) 
  pd(34,35) = pd(34,35) + rrt(421) * density(21) 
  pd(35,21) = pd(35,21) - rrt(421) * density(35) 
  pd(35,35) = pd(35,35) - rrt(421) * density(21) 
  pd(21,26) = pd(21,26) + rrt(422) * density(35) * 2.0d0
  pd(21,35) = pd(21,35) + rrt(422) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(422) * density(35) 
  pd(26,35) = pd(26,35) - rrt(422) * density(26) 
  pd(34,26) = pd(34,26) + rrt(422) * density(35) 
  pd(34,35) = pd(34,35) + rrt(422) * density(26) 
  pd(35,26) = pd(35,26) - rrt(422) * density(35) 
  pd(35,35) = pd(35,35) - rrt(422) * density(26) 
  pd(21,27) = pd(21,27) + rrt(423) * density(35) * 2.0d0
  pd(21,35) = pd(21,35) + rrt(423) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(423) * density(35) 
  pd(27,35) = pd(27,35) - rrt(423) * density(27) 
  pd(34,27) = pd(34,27) + rrt(423) * density(35) 
  pd(34,35) = pd(34,35) + rrt(423) * density(27) 
  pd(35,27) = pd(35,27) - rrt(423) * density(35) 
  pd(35,35) = pd(35,35) - rrt(423) * density(27) 
  pd(29,29) = pd(29,29) - rrt(424) * density(35) 
  pd(29,35) = pd(29,35) - rrt(424) * density(29) 
  pd(32,29) = pd(32,29) + rrt(424) * density(35) 
  pd(32,35) = pd(32,35) + rrt(424) * density(29) 
  pd(34,29) = pd(34,29) + rrt(424) * density(35) 
  pd(34,35) = pd(34,35) + rrt(424) * density(29) 
  pd(35,29) = pd(35,29) - rrt(424) * density(35) 
  pd(35,35) = pd(35,35) - rrt(424) * density(29) 
  pd(21,35) = pd(21,35) + rrt(425) * density(40) * 2.0d0
  pd(21,40) = pd(21,40) + rrt(425) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(425) * density(40) 
  pd(35,40) = pd(35,40) - rrt(425) * density(35) 
  pd(40,35) = pd(40,35) - rrt(425) * density(40) 
  pd(40,40) = pd(40,40) - rrt(425) * density(35) 
  pd(45,35) = pd(45,35) + rrt(425) * density(40) 
  pd(45,40) = pd(45,40) + rrt(425) * density(35) 
  pd(01,01) = pd(01,01) + rrt(426) * density(52) 
  pd(01,52) = pd(01,52) + rrt(426) * density(01) 
  pd(34,01) = pd(34,01) + rrt(426) * density(52) 
  pd(34,52) = pd(34,52) + rrt(426) * density(01) 
  pd(52,01) = pd(52,01) - rrt(426) * density(52) 
  pd(52,52) = pd(52,52) - rrt(426) * density(01) 
  pd(01,21) = pd(01,21) + rrt(427) * density(52) 
  pd(01,52) = pd(01,52) + rrt(427) * density(21) 
  pd(21,21) = pd(21,21) - rrt(427) * density(52) 
  pd(21,52) = pd(21,52) - rrt(427) * density(21) 
  pd(35,21) = pd(35,21) + rrt(427) * density(52) 
  pd(35,52) = pd(35,52) + rrt(427) * density(21) 
  pd(52,21) = pd(52,21) - rrt(427) * density(52) 
  pd(52,52) = pd(52,52) - rrt(427) * density(21) 
  pd(01,01) = pd(01,01) - rrt(428) * density(01) * density(17) * 2.0d0
  pd(01,17) = pd(01,17) - rrt(428) * density(01)**2 
  pd(17,01) = pd(17,01) - rrt(428) * density(01) * density(17) * 2.0d0
  pd(17,17) = pd(17,17) - rrt(428) * density(01)**2 
  pd(19,01) = pd(19,01) + rrt(428) * density(01) * density(17) * 2.0d0
  pd(19,17) = pd(19,17) + rrt(428) * density(01)**2 
  pd(17,17) = pd(17,17) - rrt(429) * density(29) 
  pd(17,29) = pd(17,29) - rrt(429) * density(17) 
  pd(29,17) = pd(29,17) - rrt(429) * density(29) 
  pd(29,29) = pd(29,29) - rrt(429) * density(17) 
  pd(45,17) = pd(45,17) + rrt(429) * density(29) 
  pd(45,29) = pd(45,29) + rrt(429) * density(17) 
  pd(14,14) = pd(14,14) - rrt(430) * density(17) 
  pd(14,17) = pd(14,17) - rrt(430) * density(14) 
  pd(17,14) = pd(17,14) - rrt(430) * density(17) 
  pd(17,17) = pd(17,17) - rrt(430) * density(14) 
  pd(18,14) = pd(18,14) + rrt(430) * density(17) 
  pd(18,17) = pd(18,17) + rrt(430) * density(14) 
  pd(01,01) = pd(01,01) - rrt(431) * density(33) 
  pd(01,33) = pd(01,33) - rrt(431) * density(01) 
  pd(14,01) = pd(14,01) + rrt(431) * density(33) 
  pd(14,33) = pd(14,33) + rrt(431) * density(01) 
  pd(33,01) = pd(33,01) - rrt(431) * density(33) 
  pd(33,33) = pd(33,33) - rrt(431) * density(01) 
  pd(45,01) = pd(45,01) + rrt(431) * density(33) 
  pd(45,33) = pd(45,33) + rrt(431) * density(01) 
  pd(29,29) = pd(29,29) - rrt(432) * density(33) 
  pd(29,33) = pd(29,33) - rrt(432) * density(29) 
  pd(33,29) = pd(33,29) - rrt(432) * density(33) 
  pd(33,33) = pd(33,33) - rrt(432) * density(29) 
  pd(34,29) = pd(34,29) + rrt(432) * density(33) 
  pd(34,33) = pd(34,33) + rrt(432) * density(29) 
  pd(14,14) = pd(14,14) - rrt(433) * density(33) 
  pd(14,33) = pd(14,33) - rrt(433) * density(14) 
  pd(33,14) = pd(33,14) - rrt(433) * density(33) 
  pd(33,33) = pd(33,33) - rrt(433) * density(14) 
  pd(45,14) = pd(45,14) + rrt(433) * density(33) 
  pd(45,33) = pd(45,33) + rrt(433) * density(14) 
  pd(01,01) = pd(01,01) - rrt(434) * density(01) * density(18) * 2.0d0
  pd(01,18) = pd(01,18) - rrt(434) * density(01)**2 
  pd(18,01) = pd(18,01) - rrt(434) * density(01) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(434) * density(01)**2 
  pd(20,01) = pd(20,01) + rrt(434) * density(01) * density(18) * 2.0d0
  pd(20,18) = pd(20,18) + rrt(434) * density(01)**2 
  pd(14,01) = pd(14,01) - rrt(435) * density(14) * density(18) 
  pd(14,14) = pd(14,14) - rrt(435) * density(01) * density(18) 
  pd(14,18) = pd(14,18) - rrt(435) * density(01) * density(14) 
  pd(18,01) = pd(18,01) - rrt(435) * density(14) * density(18) 
  pd(18,14) = pd(18,14) - rrt(435) * density(01) * density(18) 
  pd(18,18) = pd(18,18) - rrt(435) * density(01) * density(14) 
  pd(19,01) = pd(19,01) + rrt(435) * density(14) * density(18) 
  pd(19,14) = pd(19,14) + rrt(435) * density(01) * density(18) 
  pd(19,18) = pd(19,18) + rrt(435) * density(01) * density(14) 
  pd(21,21) = pd(21,21) - rrt(436) * density(21) * density(34) * 2.0d0
  pd(21,34) = pd(21,34) - rrt(436) * density(21)**2 
  pd(34,21) = pd(34,21) - rrt(436) * density(21) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(436) * density(21)**2 
  pd(35,21) = pd(35,21) + rrt(436) * density(21) * density(34) * 2.0d0
  pd(35,34) = pd(35,34) + rrt(436) * density(21)**2 
  pd(01,01) = pd(01,01) - rrt(437) * density(01) * density(34) * 2.0d0
  pd(01,34) = pd(01,34) - rrt(437) * density(01)**2 
  pd(34,01) = pd(34,01) - rrt(437) * density(01) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(437) * density(01)**2 
  pd(52,01) = pd(52,01) + rrt(437) * density(01) * density(34) * 2.0d0
  pd(52,34) = pd(52,34) + rrt(437) * density(01)**2 
  pd(26,26) = pd(26,26) - rrt(438) * density(36) 
  pd(26,36) = pd(26,36) - rrt(438) * density(26) 
  pd(29,26) = pd(29,26) + rrt(438) * density(36) 
  pd(29,36) = pd(29,36) + rrt(438) * density(26) 
  pd(36,26) = pd(36,26) - rrt(438) * density(36) 
  pd(36,36) = pd(36,36) - rrt(438) * density(26) 
  pd(37,26) = pd(37,26) + rrt(438) * density(36) 
  pd(37,36) = pd(37,36) + rrt(438) * density(26) 
  pd(29,32) = pd(29,32) + rrt(439) * density(36) 
  pd(29,36) = pd(29,36) + rrt(439) * density(32) 
  pd(32,32) = pd(32,32) - rrt(439) * density(36) 
  pd(32,36) = pd(32,36) - rrt(439) * density(32) 
  pd(36,32) = pd(36,32) - rrt(439) * density(36) 
  pd(36,36) = pd(36,36) - rrt(439) * density(32) 
  pd(38,32) = pd(38,32) + rrt(439) * density(36) 
  pd(38,36) = pd(38,36) + rrt(439) * density(32) 
  pd(29,36) = pd(29,36) + rrt(440) * density(42) 
  pd(29,42) = pd(29,42) + rrt(440) * density(36) 
  pd(36,36) = pd(36,36) - rrt(440) * density(42) 
  pd(36,42) = pd(36,42) - rrt(440) * density(36) 
  pd(42,36) = pd(42,36) - rrt(440) * density(42) 
  pd(42,42) = pd(42,42) - rrt(440) * density(36) 
  pd(50,36) = pd(50,36) + rrt(440) * density(42) 
  pd(50,42) = pd(50,42) + rrt(440) * density(36) 
  pd(36,36) = pd(36,36) - rrt(441) * density(41) 
  pd(36,41) = pd(36,41) - rrt(441) * density(36) 
  pd(40,36) = pd(40,36) + rrt(441) * density(41) 
  pd(40,41) = pd(40,41) + rrt(441) * density(36) 
  pd(41,36) = pd(41,36) - rrt(441) * density(41) 
  pd(41,41) = pd(41,41) - rrt(441) * density(36) 
  pd(48,36) = pd(48,36) + rrt(441) * density(41) 
  pd(48,41) = pd(48,41) + rrt(441) * density(36) 
  pd(29,36) = pd(29,36) + rrt(442) * density(41) 
  pd(29,41) = pd(29,41) + rrt(442) * density(36) 
  pd(36,36) = pd(36,36) - rrt(442) * density(41) 
  pd(36,41) = pd(36,41) - rrt(442) * density(36) 
  pd(41,36) = pd(41,36) - rrt(442) * density(41) 
  pd(41,41) = pd(41,41) - rrt(442) * density(36) 
  pd(49,36) = pd(49,36) + rrt(442) * density(41) 
  pd(49,41) = pd(49,41) + rrt(442) * density(36) 
  pd(21,29) = pd(21,29) + rrt(443) * density(37) 
  pd(21,37) = pd(21,37) + rrt(443) * density(29) 
  pd(29,29) = pd(29,29) - rrt(443) * density(37) 
  pd(29,37) = pd(29,37) - rrt(443) * density(29) 
  pd(36,29) = pd(36,29) + rrt(443) * density(37) 
  pd(36,37) = pd(36,37) + rrt(443) * density(29) 
  pd(37,29) = pd(37,29) - rrt(443) * density(37) 
  pd(37,37) = pd(37,37) - rrt(443) * density(29) 
  pd(21,32) = pd(21,32) + rrt(444) * density(37) 
  pd(21,37) = pd(21,37) + rrt(444) * density(32) 
  pd(32,32) = pd(32,32) - rrt(444) * density(37) 
  pd(32,37) = pd(32,37) - rrt(444) * density(32) 
  pd(37,32) = pd(37,32) - rrt(444) * density(37) 
  pd(37,37) = pd(37,37) - rrt(444) * density(32) 
  pd(38,32) = pd(38,32) + rrt(444) * density(37) 
  pd(38,37) = pd(38,37) + rrt(444) * density(32) 
  pd(21,37) = pd(21,37) + rrt(445) * density(42) 
  pd(21,42) = pd(21,42) + rrt(445) * density(37) 
  pd(37,37) = pd(37,37) - rrt(445) * density(42) 
  pd(37,42) = pd(37,42) - rrt(445) * density(37) 
  pd(42,37) = pd(42,37) - rrt(445) * density(42) 
  pd(42,42) = pd(42,42) - rrt(445) * density(37) 
  pd(50,37) = pd(50,37) + rrt(445) * density(42) 
  pd(50,42) = pd(50,42) + rrt(445) * density(37) 
  pd(21,37) = pd(21,37) + rrt(446) * density(43) 
  pd(21,43) = pd(21,43) + rrt(446) * density(37) 
  pd(37,37) = pd(37,37) - rrt(446) * density(43) 
  pd(37,43) = pd(37,43) - rrt(446) * density(37) 
  pd(43,37) = pd(43,37) - rrt(446) * density(43) 
  pd(43,43) = pd(43,43) - rrt(446) * density(37) 
  pd(51,37) = pd(51,37) + rrt(446) * density(43) 
  pd(51,43) = pd(51,43) + rrt(446) * density(37) 
  pd(21,29) = pd(21,29) + rrt(447) * density(38) 
  pd(21,38) = pd(21,38) + rrt(447) * density(29) 
  pd(29,29) = pd(29,29) - rrt(447) * density(38) 
  pd(29,38) = pd(29,38) - rrt(447) * density(29) 
  pd(37,29) = pd(37,29) + rrt(447) * density(38) 
  pd(37,38) = pd(37,38) + rrt(447) * density(29) 
  pd(38,29) = pd(38,29) - rrt(447) * density(38) 
  pd(38,38) = pd(38,38) - rrt(447) * density(29) 
  pd(29,38) = pd(29,38) + rrt(448) * density(40) 
  pd(29,40) = pd(29,40) + rrt(448) * density(38) 
  pd(38,38) = pd(38,38) - rrt(448) * density(40) 
  pd(38,40) = pd(38,40) - rrt(448) * density(38) 
  pd(40,38) = pd(40,38) - rrt(448) * density(40) 
  pd(40,40) = pd(40,40) - rrt(448) * density(38) 
  pd(51,38) = pd(51,38) + rrt(448) * density(40) 
  pd(51,40) = pd(51,40) + rrt(448) * density(38) 
  pd(21,38) = pd(21,38) + rrt(449) * density(40) 
  pd(21,40) = pd(21,40) + rrt(449) * density(38) 
  pd(38,38) = pd(38,38) - rrt(449) * density(40) 
  pd(38,40) = pd(38,40) - rrt(449) * density(38) 
  pd(40,38) = pd(40,38) - rrt(449) * density(40) 
  pd(40,40) = pd(40,40) - rrt(449) * density(38) 
  pd(50,38) = pd(50,38) + rrt(449) * density(40) 
  pd(50,40) = pd(50,40) + rrt(449) * density(38) 
  pd(32,38) = pd(32,38) + rrt(450) * density(42) 
  pd(32,42) = pd(32,42) + rrt(450) * density(38) 
  pd(38,38) = pd(38,38) - rrt(450) * density(42) 
  pd(38,42) = pd(38,42) - rrt(450) * density(38) 
  pd(42,38) = pd(42,38) - rrt(450) * density(42) 
  pd(42,42) = pd(42,42) - rrt(450) * density(38) 
  pd(50,38) = pd(50,38) + rrt(450) * density(42) 
  pd(50,42) = pd(50,42) + rrt(450) * density(38) 
  pd(21,38) = pd(21,38) + rrt(451) * density(42) 
  pd(21,42) = pd(21,42) + rrt(451) * density(38) 
  pd(38,38) = pd(38,38) - rrt(451) * density(42) 
  pd(38,42) = pd(38,42) - rrt(451) * density(38) 
  pd(42,38) = pd(42,38) - rrt(451) * density(42) 
  pd(42,42) = pd(42,42) - rrt(451) * density(38) 
  pd(51,38) = pd(51,38) + rrt(451) * density(42) 
  pd(51,42) = pd(51,42) + rrt(451) * density(38) 
  pd(32,38) = pd(32,38) + rrt(452) * density(43) 
  pd(32,43) = pd(32,43) + rrt(452) * density(38) 
  pd(38,38) = pd(38,38) - rrt(452) * density(43) 
  pd(38,43) = pd(38,43) - rrt(452) * density(38) 
  pd(43,38) = pd(43,38) - rrt(452) * density(43) 
  pd(43,43) = pd(43,43) - rrt(452) * density(38) 
  pd(51,38) = pd(51,38) + rrt(452) * density(43) 
  pd(51,43) = pd(51,43) + rrt(452) * density(38) 
  pd(21,21) = pd(21,21) - rrt(453) * density(48) 
  pd(21,48) = pd(21,48) - rrt(453) * density(21) 
  pd(37,21) = pd(37,21) + rrt(453) * density(48) 
  pd(37,48) = pd(37,48) + rrt(453) * density(21) 
  pd(40,21) = pd(40,21) + rrt(453) * density(48) 
  pd(40,48) = pd(40,48) + rrt(453) * density(21) 
  pd(48,21) = pd(48,21) - rrt(453) * density(48) 
  pd(48,48) = pd(48,48) - rrt(453) * density(21) 
  pd(40,42) = pd(40,42) + rrt(454) * density(48) 
  pd(40,48) = pd(40,48) + rrt(454) * density(42) 
  pd(42,42) = pd(42,42) - rrt(454) * density(48) 
  pd(42,48) = pd(42,48) - rrt(454) * density(42) 
  pd(48,42) = pd(48,42) - rrt(454) * density(48) 
  pd(48,48) = pd(48,48) - rrt(454) * density(42) 
  pd(50,42) = pd(50,42) + rrt(454) * density(48) 
  pd(50,48) = pd(50,48) + rrt(454) * density(42) 
  pd(01,41) = pd(01,41) + rrt(455) * density(48) 
  pd(01,48) = pd(01,48) + rrt(455) * density(41) 
  pd(41,41) = pd(41,41) - rrt(455) * density(48) 
  pd(41,48) = pd(41,48) - rrt(455) * density(41) 
  pd(48,41) = pd(48,41) - rrt(455) * density(48) 
  pd(48,48) = pd(48,48) - rrt(455) * density(41) 
  pd(50,41) = pd(50,41) + rrt(455) * density(48) 
  pd(50,48) = pd(50,48) + rrt(455) * density(41) 
  pd(21,32) = pd(21,32) + rrt(456) * density(50) 
  pd(21,50) = pd(21,50) + rrt(456) * density(32) 
  pd(32,32) = pd(32,32) - rrt(456) * density(50) 
  pd(32,50) = pd(32,50) - rrt(456) * density(32) 
  pd(50,32) = pd(50,32) - rrt(456) * density(50) 
  pd(50,50) = pd(50,50) - rrt(456) * density(32) 
  pd(51,32) = pd(51,32) + rrt(456) * density(50) 
  pd(51,50) = pd(51,50) + rrt(456) * density(32) 
  pd(40,42) = pd(40,42) + rrt(457) * density(50) 
  pd(40,50) = pd(40,50) + rrt(457) * density(42) 
  pd(42,42) = pd(42,42) - rrt(457) * density(50) 
  pd(42,50) = pd(42,50) - rrt(457) * density(42) 
  pd(50,42) = pd(50,42) - rrt(457) * density(50) 
  pd(50,50) = pd(50,50) - rrt(457) * density(42) 
  pd(51,42) = pd(51,42) + rrt(457) * density(50) 
  pd(51,50) = pd(51,50) + rrt(457) * density(42) 
  pd(42,43) = pd(42,43) + rrt(458) * density(50) 
  pd(42,50) = pd(42,50) + rrt(458) * density(43) 
  pd(43,43) = pd(43,43) - rrt(458) * density(50) 
  pd(43,50) = pd(43,50) - rrt(458) * density(43) 
  pd(50,43) = pd(50,43) - rrt(458) * density(50) 
  pd(50,50) = pd(50,50) - rrt(458) * density(43) 
  pd(51,43) = pd(51,43) + rrt(458) * density(50) 
  pd(51,50) = pd(51,50) + rrt(458) * density(43) 
  pd(42,44) = pd(42,44) + rrt(459) * density(50) * 2.0d0
  pd(42,50) = pd(42,50) + rrt(459) * density(44) * 2.0d0
  pd(44,44) = pd(44,44) - rrt(459) * density(50) 
  pd(44,50) = pd(44,50) - rrt(459) * density(44) 
  pd(50,44) = pd(50,44) - rrt(459) * density(50) 
  pd(50,50) = pd(50,50) - rrt(459) * density(44) 
  pd(51,44) = pd(51,44) + rrt(459) * density(50) 
  pd(51,50) = pd(51,50) + rrt(459) * density(44) 
  pd(40,40) = pd(40,40) - rrt(460) * density(51) 
  pd(40,51) = pd(40,51) - rrt(460) * density(40) 
  pd(42,40) = pd(42,40) + rrt(460) * density(51) 
  pd(42,51) = pd(42,51) + rrt(460) * density(40) 
  pd(50,40) = pd(50,40) + rrt(460) * density(51) 
  pd(50,51) = pd(50,51) + rrt(460) * density(40) 
  pd(51,40) = pd(51,40) - rrt(460) * density(51) 
  pd(51,51) = pd(51,51) - rrt(460) * density(40) 
  pd(21,01) = pd(21,01) + rrt(461) * density(39) 
  pd(21,39) = pd(21,39) + rrt(461) * density(01) 
  pd(37,01) = pd(37,01) + rrt(461) * density(39) 
  pd(37,39) = pd(37,39) + rrt(461) * density(01) 
  pd(39,01) = pd(39,01) - rrt(461) * density(39) 
  pd(39,39) = pd(39,39) - rrt(461) * density(01) 
  pd(21,21) = pd(21,21) + rrt(462) * density(39) 
  pd(21,39) = pd(21,39) + rrt(462) * density(21) 
  pd(37,21) = pd(37,21) + rrt(462) * density(39) 
  pd(37,39) = pd(37,39) + rrt(462) * density(21) 
  pd(39,21) = pd(39,21) - rrt(462) * density(39) 
  pd(39,39) = pd(39,39) - rrt(462) * density(21) 
  pd(21,29) = pd(21,29) + rrt(463) * density(39) 
  pd(21,39) = pd(21,39) + rrt(463) * density(29) 
  pd(29,29) = pd(29,29) - rrt(463) * density(39) 
  pd(29,39) = pd(29,39) - rrt(463) * density(29) 
  pd(38,29) = pd(38,29) + rrt(463) * density(39) 
  pd(38,39) = pd(38,39) + rrt(463) * density(29) 
  pd(39,29) = pd(39,29) - rrt(463) * density(39) 
  pd(39,39) = pd(39,39) - rrt(463) * density(29) 
  pd(21,29) = pd(21,29) + rrt(464) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(464) * density(29) * 2.0d0
  pd(29,29) = pd(29,29) - rrt(464) * density(39) 
  pd(29,39) = pd(29,39) - rrt(464) * density(29) 
  pd(36,29) = pd(36,29) + rrt(464) * density(39) 
  pd(36,39) = pd(36,39) + rrt(464) * density(29) 
  pd(39,29) = pd(39,29) - rrt(464) * density(39) 
  pd(39,39) = pd(39,39) - rrt(464) * density(29) 
  pd(21,26) = pd(21,26) + rrt(465) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(465) * density(26) * 2.0d0
  pd(26,26) = pd(26,26) - rrt(465) * density(39) 
  pd(26,39) = pd(26,39) - rrt(465) * density(26) 
  pd(37,26) = pd(37,26) + rrt(465) * density(39) 
  pd(37,39) = pd(37,39) + rrt(465) * density(26) 
  pd(39,26) = pd(39,26) - rrt(465) * density(39) 
  pd(39,39) = pd(39,39) - rrt(465) * density(26) 
  pd(21,27) = pd(21,27) + rrt(466) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(466) * density(27) * 2.0d0
  pd(27,27) = pd(27,27) - rrt(466) * density(39) 
  pd(27,39) = pd(27,39) - rrt(466) * density(27) 
  pd(37,27) = pd(37,27) + rrt(466) * density(39) 
  pd(37,39) = pd(37,39) + rrt(466) * density(27) 
  pd(39,27) = pd(39,27) - rrt(466) * density(39) 
  pd(39,39) = pd(39,39) - rrt(466) * density(27) 
  pd(21,39) = pd(21,39) + rrt(467) * density(40) 
  pd(21,40) = pd(21,40) + rrt(467) * density(39) 
  pd(39,39) = pd(39,39) - rrt(467) * density(40) 
  pd(39,40) = pd(39,40) - rrt(467) * density(39) 
  pd(40,39) = pd(40,39) - rrt(467) * density(40) 
  pd(40,40) = pd(40,40) - rrt(467) * density(39) 
  pd(51,39) = pd(51,39) + rrt(467) * density(40) 
  pd(51,40) = pd(51,40) + rrt(467) * density(39) 
  pd(21,21) = pd(21,21) - rrt(468) * density(36) 
  pd(21,36) = pd(21,36) - rrt(468) * density(21) 
  pd(36,21) = pd(36,21) - rrt(468) * density(36) 
  pd(36,36) = pd(36,36) - rrt(468) * density(21) 
  pd(38,21) = pd(38,21) + rrt(468) * density(36) 
  pd(38,36) = pd(38,36) + rrt(468) * density(21) 
  pd(36,36) = pd(36,36) - rrt(469) * density(40) 
  pd(36,40) = pd(36,40) - rrt(469) * density(36) 
  pd(40,36) = pd(40,36) - rrt(469) * density(40) 
  pd(40,40) = pd(40,40) - rrt(469) * density(36) 
  pd(50,36) = pd(50,36) + rrt(469) * density(40) 
  pd(50,40) = pd(50,40) + rrt(469) * density(36) 
  pd(21,21) = pd(21,21) - rrt(470) * density(37) 
  pd(21,37) = pd(21,37) - rrt(470) * density(21) 
  pd(37,21) = pd(37,21) - rrt(470) * density(37) 
  pd(37,37) = pd(37,37) - rrt(470) * density(21) 
  pd(39,21) = pd(39,21) + rrt(470) * density(37) 
  pd(39,37) = pd(39,37) + rrt(470) * density(21) 
  pd(14,17) = pd(14,17) + rrt(471) * density(36) 
  pd(14,36) = pd(14,36) + rrt(471) * density(17) 
  pd(17,17) = pd(17,17) - rrt(471) * density(36) 
  pd(17,36) = pd(17,36) - rrt(471) * density(17) 
  pd(29,17) = pd(29,17) + rrt(471) * density(36) 
  pd(29,36) = pd(29,36) + rrt(471) * density(17) 
  pd(36,17) = pd(36,17) - rrt(471) * density(36) 
  pd(36,36) = pd(36,36) - rrt(471) * density(17) 
  pd(01,18) = pd(01,18) + rrt(472) * density(36) 
  pd(01,36) = pd(01,36) + rrt(472) * density(18) 
  pd(18,18) = pd(18,18) - rrt(472) * density(36) 
  pd(18,36) = pd(18,36) - rrt(472) * density(18) 
  pd(29,18) = pd(29,18) + rrt(472) * density(36) 
  pd(29,36) = pd(29,36) + rrt(472) * density(18) 
  pd(36,18) = pd(36,18) - rrt(472) * density(36) 
  pd(36,36) = pd(36,36) - rrt(472) * density(18) 
  pd(29,33) = pd(29,33) + rrt(473) * density(36) * 2.0d0
  pd(29,36) = pd(29,36) + rrt(473) * density(33) * 2.0d0
  pd(33,33) = pd(33,33) - rrt(473) * density(36) 
  pd(33,36) = pd(33,36) - rrt(473) * density(33) 
  pd(36,33) = pd(36,33) - rrt(473) * density(36) 
  pd(36,36) = pd(36,36) - rrt(473) * density(33) 
  pd(21,34) = pd(21,34) + rrt(474) * density(36) 
  pd(21,36) = pd(21,36) + rrt(474) * density(34) 
  pd(29,34) = pd(29,34) + rrt(474) * density(36) 
  pd(29,36) = pd(29,36) + rrt(474) * density(34) 
  pd(34,34) = pd(34,34) - rrt(474) * density(36) 
  pd(34,36) = pd(34,36) - rrt(474) * density(34) 
  pd(36,34) = pd(36,34) - rrt(474) * density(36) 
  pd(36,36) = pd(36,36) - rrt(474) * density(34) 
  pd(29,36) = pd(29,36) + rrt(475) * density(45) 
  pd(29,45) = pd(29,45) + rrt(475) * density(36) 
  pd(36,36) = pd(36,36) - rrt(475) * density(45) 
  pd(36,45) = pd(36,45) - rrt(475) * density(36) 
  pd(40,36) = pd(40,36) + rrt(475) * density(45) 
  pd(40,45) = pd(40,45) + rrt(475) * density(36) 
  pd(45,36) = pd(45,36) - rrt(475) * density(45) 
  pd(45,45) = pd(45,45) - rrt(475) * density(36) 
  pd(29,36) = pd(29,36) + rrt(476) * density(46) 
  pd(29,46) = pd(29,46) + rrt(476) * density(36) 
  pd(36,36) = pd(36,36) - rrt(476) * density(46) 
  pd(36,46) = pd(36,46) - rrt(476) * density(36) 
  pd(41,36) = pd(41,36) + rrt(476) * density(46) 
  pd(41,46) = pd(41,46) + rrt(476) * density(36) 
  pd(46,36) = pd(46,36) - rrt(476) * density(46) 
  pd(46,46) = pd(46,46) - rrt(476) * density(36) 
  pd(29,36) = pd(29,36) + rrt(477) * density(47) 
  pd(29,47) = pd(29,47) + rrt(477) * density(36) 
  pd(36,36) = pd(36,36) - rrt(477) * density(47) 
  pd(36,47) = pd(36,47) - rrt(477) * density(36) 
  pd(42,36) = pd(42,36) + rrt(477) * density(47) 
  pd(42,47) = pd(42,47) + rrt(477) * density(36) 
  pd(47,36) = pd(47,36) - rrt(477) * density(47) 
  pd(47,47) = pd(47,47) - rrt(477) * density(36) 
  pd(14,17) = pd(14,17) + rrt(478) * density(37) 
  pd(14,37) = pd(14,37) + rrt(478) * density(17) 
  pd(17,17) = pd(17,17) - rrt(478) * density(37) 
  pd(17,37) = pd(17,37) - rrt(478) * density(17) 
  pd(21,17) = pd(21,17) + rrt(478) * density(37) 
  pd(21,37) = pd(21,37) + rrt(478) * density(17) 
  pd(37,17) = pd(37,17) - rrt(478) * density(37) 
  pd(37,37) = pd(37,37) - rrt(478) * density(17) 
  pd(01,18) = pd(01,18) + rrt(479) * density(37) 
  pd(01,37) = pd(01,37) + rrt(479) * density(18) 
  pd(18,18) = pd(18,18) - rrt(479) * density(37) 
  pd(18,37) = pd(18,37) - rrt(479) * density(18) 
  pd(21,18) = pd(21,18) + rrt(479) * density(37) 
  pd(21,37) = pd(21,37) + rrt(479) * density(18) 
  pd(37,18) = pd(37,18) - rrt(479) * density(37) 
  pd(37,37) = pd(37,37) - rrt(479) * density(18) 
  pd(21,33) = pd(21,33) + rrt(480) * density(37) 
  pd(21,37) = pd(21,37) + rrt(480) * density(33) 
  pd(29,33) = pd(29,33) + rrt(480) * density(37) 
  pd(29,37) = pd(29,37) + rrt(480) * density(33) 
  pd(33,33) = pd(33,33) - rrt(480) * density(37) 
  pd(33,37) = pd(33,37) - rrt(480) * density(33) 
  pd(37,33) = pd(37,33) - rrt(480) * density(37) 
  pd(37,37) = pd(37,37) - rrt(480) * density(33) 
  pd(21,34) = pd(21,34) + rrt(481) * density(37) * 2.0d0
  pd(21,37) = pd(21,37) + rrt(481) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(481) * density(37) 
  pd(34,37) = pd(34,37) - rrt(481) * density(34) 
  pd(37,34) = pd(37,34) - rrt(481) * density(37) 
  pd(37,37) = pd(37,37) - rrt(481) * density(34) 
  pd(21,37) = pd(21,37) + rrt(482) * density(45) 
  pd(21,45) = pd(21,45) + rrt(482) * density(37) 
  pd(37,37) = pd(37,37) - rrt(482) * density(45) 
  pd(37,45) = pd(37,45) - rrt(482) * density(37) 
  pd(40,37) = pd(40,37) + rrt(482) * density(45) 
  pd(40,45) = pd(40,45) + rrt(482) * density(37) 
  pd(45,37) = pd(45,37) - rrt(482) * density(45) 
  pd(45,45) = pd(45,45) - rrt(482) * density(37) 
  pd(21,37) = pd(21,37) + rrt(483) * density(46) 
  pd(21,46) = pd(21,46) + rrt(483) * density(37) 
  pd(37,37) = pd(37,37) - rrt(483) * density(46) 
  pd(37,46) = pd(37,46) - rrt(483) * density(37) 
  pd(41,37) = pd(41,37) + rrt(483) * density(46) 
  pd(41,46) = pd(41,46) + rrt(483) * density(37) 
  pd(46,37) = pd(46,37) - rrt(483) * density(46) 
  pd(46,46) = pd(46,46) - rrt(483) * density(37) 
  pd(21,37) = pd(21,37) + rrt(484) * density(47) 
  pd(21,47) = pd(21,47) + rrt(484) * density(37) 
  pd(37,37) = pd(37,37) - rrt(484) * density(47) 
  pd(37,47) = pd(37,47) - rrt(484) * density(37) 
  pd(42,37) = pd(42,37) + rrt(484) * density(47) 
  pd(42,47) = pd(42,47) + rrt(484) * density(37) 
  pd(47,37) = pd(47,37) - rrt(484) * density(47) 
  pd(47,47) = pd(47,47) - rrt(484) * density(37) 
  pd(14,17) = pd(14,17) + rrt(485) * density(38) 
  pd(14,38) = pd(14,38) + rrt(485) * density(17) 
  pd(17,17) = pd(17,17) - rrt(485) * density(38) 
  pd(17,38) = pd(17,38) - rrt(485) * density(17) 
  pd(32,17) = pd(32,17) + rrt(485) * density(38) 
  pd(32,38) = pd(32,38) + rrt(485) * density(17) 
  pd(38,17) = pd(38,17) - rrt(485) * density(38) 
  pd(38,38) = pd(38,38) - rrt(485) * density(17) 
  pd(01,18) = pd(01,18) + rrt(486) * density(38) 
  pd(01,38) = pd(01,38) + rrt(486) * density(18) 
  pd(18,18) = pd(18,18) - rrt(486) * density(38) 
  pd(18,38) = pd(18,38) - rrt(486) * density(18) 
  pd(32,18) = pd(32,18) + rrt(486) * density(38) 
  pd(32,38) = pd(32,38) + rrt(486) * density(18) 
  pd(38,18) = pd(38,18) - rrt(486) * density(38) 
  pd(38,38) = pd(38,38) - rrt(486) * density(18) 
  pd(29,33) = pd(29,33) + rrt(487) * density(38) 
  pd(29,38) = pd(29,38) + rrt(487) * density(33) 
  pd(32,33) = pd(32,33) + rrt(487) * density(38) 
  pd(32,38) = pd(32,38) + rrt(487) * density(33) 
  pd(33,33) = pd(33,33) - rrt(487) * density(38) 
  pd(33,38) = pd(33,38) - rrt(487) * density(33) 
  pd(38,33) = pd(38,33) - rrt(487) * density(38) 
  pd(38,38) = pd(38,38) - rrt(487) * density(33) 
  pd(21,34) = pd(21,34) + rrt(488) * density(38) 
  pd(21,38) = pd(21,38) + rrt(488) * density(34) 
  pd(32,34) = pd(32,34) + rrt(488) * density(38) 
  pd(32,38) = pd(32,38) + rrt(488) * density(34) 
  pd(34,34) = pd(34,34) - rrt(488) * density(38) 
  pd(34,38) = pd(34,38) - rrt(488) * density(34) 
  pd(38,34) = pd(38,34) - rrt(488) * density(38) 
  pd(38,38) = pd(38,38) - rrt(488) * density(34) 
  pd(32,38) = pd(32,38) + rrt(489) * density(45) 
  pd(32,45) = pd(32,45) + rrt(489) * density(38) 
  pd(38,38) = pd(38,38) - rrt(489) * density(45) 
  pd(38,45) = pd(38,45) - rrt(489) * density(38) 
  pd(40,38) = pd(40,38) + rrt(489) * density(45) 
  pd(40,45) = pd(40,45) + rrt(489) * density(38) 
  pd(45,38) = pd(45,38) - rrt(489) * density(45) 
  pd(45,45) = pd(45,45) - rrt(489) * density(38) 
  pd(32,38) = pd(32,38) + rrt(490) * density(46) 
  pd(32,46) = pd(32,46) + rrt(490) * density(38) 
  pd(38,38) = pd(38,38) - rrt(490) * density(46) 
  pd(38,46) = pd(38,46) - rrt(490) * density(38) 
  pd(41,38) = pd(41,38) + rrt(490) * density(46) 
  pd(41,46) = pd(41,46) + rrt(490) * density(38) 
  pd(46,38) = pd(46,38) - rrt(490) * density(46) 
  pd(46,46) = pd(46,46) - rrt(490) * density(38) 
  pd(32,38) = pd(32,38) + rrt(491) * density(47) 
  pd(32,47) = pd(32,47) + rrt(491) * density(38) 
  pd(38,38) = pd(38,38) - rrt(491) * density(47) 
  pd(38,47) = pd(38,47) - rrt(491) * density(38) 
  pd(42,38) = pd(42,38) + rrt(491) * density(47) 
  pd(42,47) = pd(42,47) + rrt(491) * density(38) 
  pd(47,38) = pd(47,38) - rrt(491) * density(47) 
  pd(47,47) = pd(47,47) - rrt(491) * density(38) 
  pd(14,17) = pd(14,17) + rrt(492) * density(48) 
  pd(14,48) = pd(14,48) + rrt(492) * density(17) 
  pd(17,17) = pd(17,17) - rrt(492) * density(48) 
  pd(17,48) = pd(17,48) - rrt(492) * density(17) 
  pd(40,17) = pd(40,17) + rrt(492) * density(48) 
  pd(40,48) = pd(40,48) + rrt(492) * density(17) 
  pd(48,17) = pd(48,17) - rrt(492) * density(48) 
  pd(48,48) = pd(48,48) - rrt(492) * density(17) 
  pd(01,18) = pd(01,18) + rrt(493) * density(48) 
  pd(01,48) = pd(01,48) + rrt(493) * density(18) 
  pd(18,18) = pd(18,18) - rrt(493) * density(48) 
  pd(18,48) = pd(18,48) - rrt(493) * density(18) 
  pd(40,18) = pd(40,18) + rrt(493) * density(48) 
  pd(40,48) = pd(40,48) + rrt(493) * density(18) 
  pd(48,18) = pd(48,18) - rrt(493) * density(48) 
  pd(48,48) = pd(48,48) - rrt(493) * density(18) 
  pd(29,33) = pd(29,33) + rrt(494) * density(48) 
  pd(29,48) = pd(29,48) + rrt(494) * density(33) 
  pd(33,33) = pd(33,33) - rrt(494) * density(48) 
  pd(33,48) = pd(33,48) - rrt(494) * density(33) 
  pd(40,33) = pd(40,33) + rrt(494) * density(48) 
  pd(40,48) = pd(40,48) + rrt(494) * density(33) 
  pd(48,33) = pd(48,33) - rrt(494) * density(48) 
  pd(48,48) = pd(48,48) - rrt(494) * density(33) 
  pd(21,34) = pd(21,34) + rrt(495) * density(48) 
  pd(21,48) = pd(21,48) + rrt(495) * density(34) 
  pd(34,34) = pd(34,34) - rrt(495) * density(48) 
  pd(34,48) = pd(34,48) - rrt(495) * density(34) 
  pd(40,34) = pd(40,34) + rrt(495) * density(48) 
  pd(40,48) = pd(40,48) + rrt(495) * density(34) 
  pd(48,34) = pd(48,34) - rrt(495) * density(48) 
  pd(48,48) = pd(48,48) - rrt(495) * density(34) 
  pd(40,45) = pd(40,45) + rrt(496) * density(48) * 2.0d0
  pd(40,48) = pd(40,48) + rrt(496) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(496) * density(48) 
  pd(45,48) = pd(45,48) - rrt(496) * density(45) 
  pd(48,45) = pd(48,45) - rrt(496) * density(48) 
  pd(48,48) = pd(48,48) - rrt(496) * density(45) 
  pd(40,46) = pd(40,46) + rrt(497) * density(48) 
  pd(40,48) = pd(40,48) + rrt(497) * density(46) 
  pd(41,46) = pd(41,46) + rrt(497) * density(48) 
  pd(41,48) = pd(41,48) + rrt(497) * density(46) 
  pd(46,46) = pd(46,46) - rrt(497) * density(48) 
  pd(46,48) = pd(46,48) - rrt(497) * density(46) 
  pd(48,46) = pd(48,46) - rrt(497) * density(48) 
  pd(48,48) = pd(48,48) - rrt(497) * density(46) 
  pd(40,47) = pd(40,47) + rrt(498) * density(48) 
  pd(40,48) = pd(40,48) + rrt(498) * density(47) 
  pd(42,47) = pd(42,47) + rrt(498) * density(48) 
  pd(42,48) = pd(42,48) + rrt(498) * density(47) 
  pd(47,47) = pd(47,47) - rrt(498) * density(48) 
  pd(47,48) = pd(47,48) - rrt(498) * density(47) 
  pd(48,47) = pd(48,47) - rrt(498) * density(48) 
  pd(48,48) = pd(48,48) - rrt(498) * density(47) 
  pd(14,17) = pd(14,17) + rrt(499) * density(49) 
  pd(14,49) = pd(14,49) + rrt(499) * density(17) 
  pd(17,17) = pd(17,17) - rrt(499) * density(49) 
  pd(17,49) = pd(17,49) - rrt(499) * density(17) 
  pd(41,17) = pd(41,17) + rrt(499) * density(49) 
  pd(41,49) = pd(41,49) + rrt(499) * density(17) 
  pd(49,17) = pd(49,17) - rrt(499) * density(49) 
  pd(49,49) = pd(49,49) - rrt(499) * density(17) 
  pd(01,18) = pd(01,18) + rrt(500) * density(49) 
  pd(01,49) = pd(01,49) + rrt(500) * density(18) 
  pd(18,18) = pd(18,18) - rrt(500) * density(49) 
  pd(18,49) = pd(18,49) - rrt(500) * density(18) 
  pd(41,18) = pd(41,18) + rrt(500) * density(49) 
  pd(41,49) = pd(41,49) + rrt(500) * density(18) 
  pd(49,18) = pd(49,18) - rrt(500) * density(49) 
  pd(49,49) = pd(49,49) - rrt(500) * density(18) 
  pd(29,33) = pd(29,33) + rrt(501) * density(49) 
  pd(29,49) = pd(29,49) + rrt(501) * density(33) 
  pd(33,33) = pd(33,33) - rrt(501) * density(49) 
  pd(33,49) = pd(33,49) - rrt(501) * density(33) 
  pd(41,33) = pd(41,33) + rrt(501) * density(49) 
  pd(41,49) = pd(41,49) + rrt(501) * density(33) 
  pd(49,33) = pd(49,33) - rrt(501) * density(49) 
  pd(49,49) = pd(49,49) - rrt(501) * density(33) 
  pd(21,34) = pd(21,34) + rrt(502) * density(49) 
  pd(21,49) = pd(21,49) + rrt(502) * density(34) 
  pd(34,34) = pd(34,34) - rrt(502) * density(49) 
  pd(34,49) = pd(34,49) - rrt(502) * density(34) 
  pd(41,34) = pd(41,34) + rrt(502) * density(49) 
  pd(41,49) = pd(41,49) + rrt(502) * density(34) 
  pd(49,34) = pd(49,34) - rrt(502) * density(49) 
  pd(49,49) = pd(49,49) - rrt(502) * density(34) 
  pd(40,45) = pd(40,45) + rrt(503) * density(49) 
  pd(40,49) = pd(40,49) + rrt(503) * density(45) 
  pd(41,45) = pd(41,45) + rrt(503) * density(49) 
  pd(41,49) = pd(41,49) + rrt(503) * density(45) 
  pd(45,45) = pd(45,45) - rrt(503) * density(49) 
  pd(45,49) = pd(45,49) - rrt(503) * density(45) 
  pd(49,45) = pd(49,45) - rrt(503) * density(49) 
  pd(49,49) = pd(49,49) - rrt(503) * density(45) 
  pd(41,46) = pd(41,46) + rrt(504) * density(49) * 2.0d0
  pd(41,49) = pd(41,49) + rrt(504) * density(46) * 2.0d0
  pd(46,46) = pd(46,46) - rrt(504) * density(49) 
  pd(46,49) = pd(46,49) - rrt(504) * density(46) 
  pd(49,46) = pd(49,46) - rrt(504) * density(49) 
  pd(49,49) = pd(49,49) - rrt(504) * density(46) 
  pd(41,47) = pd(41,47) + rrt(505) * density(49) 
  pd(41,49) = pd(41,49) + rrt(505) * density(47) 
  pd(42,47) = pd(42,47) + rrt(505) * density(49) 
  pd(42,49) = pd(42,49) + rrt(505) * density(47) 
  pd(47,47) = pd(47,47) - rrt(505) * density(49) 
  pd(47,49) = pd(47,49) - rrt(505) * density(47) 
  pd(49,47) = pd(49,47) - rrt(505) * density(49) 
  pd(49,49) = pd(49,49) - rrt(505) * density(47) 
  pd(14,17) = pd(14,17) + rrt(506) * density(50) 
  pd(14,50) = pd(14,50) + rrt(506) * density(17) 
  pd(17,17) = pd(17,17) - rrt(506) * density(50) 
  pd(17,50) = pd(17,50) - rrt(506) * density(17) 
  pd(42,17) = pd(42,17) + rrt(506) * density(50) 
  pd(42,50) = pd(42,50) + rrt(506) * density(17) 
  pd(50,17) = pd(50,17) - rrt(506) * density(50) 
  pd(50,50) = pd(50,50) - rrt(506) * density(17) 
  pd(01,18) = pd(01,18) + rrt(507) * density(50) 
  pd(01,50) = pd(01,50) + rrt(507) * density(18) 
  pd(18,18) = pd(18,18) - rrt(507) * density(50) 
  pd(18,50) = pd(18,50) - rrt(507) * density(18) 
  pd(42,18) = pd(42,18) + rrt(507) * density(50) 
  pd(42,50) = pd(42,50) + rrt(507) * density(18) 
  pd(50,18) = pd(50,18) - rrt(507) * density(50) 
  pd(50,50) = pd(50,50) - rrt(507) * density(18) 
  pd(29,33) = pd(29,33) + rrt(508) * density(50) 
  pd(29,50) = pd(29,50) + rrt(508) * density(33) 
  pd(33,33) = pd(33,33) - rrt(508) * density(50) 
  pd(33,50) = pd(33,50) - rrt(508) * density(33) 
  pd(42,33) = pd(42,33) + rrt(508) * density(50) 
  pd(42,50) = pd(42,50) + rrt(508) * density(33) 
  pd(50,33) = pd(50,33) - rrt(508) * density(50) 
  pd(50,50) = pd(50,50) - rrt(508) * density(33) 
  pd(21,34) = pd(21,34) + rrt(509) * density(50) 
  pd(21,50) = pd(21,50) + rrt(509) * density(34) 
  pd(34,34) = pd(34,34) - rrt(509) * density(50) 
  pd(34,50) = pd(34,50) - rrt(509) * density(34) 
  pd(42,34) = pd(42,34) + rrt(509) * density(50) 
  pd(42,50) = pd(42,50) + rrt(509) * density(34) 
  pd(50,34) = pd(50,34) - rrt(509) * density(50) 
  pd(50,50) = pd(50,50) - rrt(509) * density(34) 
  pd(40,45) = pd(40,45) + rrt(510) * density(50) 
  pd(40,50) = pd(40,50) + rrt(510) * density(45) 
  pd(42,45) = pd(42,45) + rrt(510) * density(50) 
  pd(42,50) = pd(42,50) + rrt(510) * density(45) 
  pd(45,45) = pd(45,45) - rrt(510) * density(50) 
  pd(45,50) = pd(45,50) - rrt(510) * density(45) 
  pd(50,45) = pd(50,45) - rrt(510) * density(50) 
  pd(50,50) = pd(50,50) - rrt(510) * density(45) 
  pd(41,46) = pd(41,46) + rrt(511) * density(50) 
  pd(41,50) = pd(41,50) + rrt(511) * density(46) 
  pd(42,46) = pd(42,46) + rrt(511) * density(50) 
  pd(42,50) = pd(42,50) + rrt(511) * density(46) 
  pd(46,46) = pd(46,46) - rrt(511) * density(50) 
  pd(46,50) = pd(46,50) - rrt(511) * density(46) 
  pd(50,46) = pd(50,46) - rrt(511) * density(50) 
  pd(50,50) = pd(50,50) - rrt(511) * density(46) 
  pd(42,47) = pd(42,47) + rrt(512) * density(50) * 2.0d0
  pd(42,50) = pd(42,50) + rrt(512) * density(47) * 2.0d0
  pd(47,47) = pd(47,47) - rrt(512) * density(50) 
  pd(47,50) = pd(47,50) - rrt(512) * density(47) 
  pd(50,47) = pd(50,47) - rrt(512) * density(50) 
  pd(50,50) = pd(50,50) - rrt(512) * density(47) 
  pd(14,17) = pd(14,17) + rrt(513) * density(51) 
  pd(14,51) = pd(14,51) + rrt(513) * density(17) 
  pd(17,17) = pd(17,17) - rrt(513) * density(51) 
  pd(17,51) = pd(17,51) - rrt(513) * density(17) 
  pd(43,17) = pd(43,17) + rrt(513) * density(51) 
  pd(43,51) = pd(43,51) + rrt(513) * density(17) 
  pd(51,17) = pd(51,17) - rrt(513) * density(51) 
  pd(51,51) = pd(51,51) - rrt(513) * density(17) 
  pd(01,18) = pd(01,18) + rrt(514) * density(51) 
  pd(01,51) = pd(01,51) + rrt(514) * density(18) 
  pd(18,18) = pd(18,18) - rrt(514) * density(51) 
  pd(18,51) = pd(18,51) - rrt(514) * density(18) 
  pd(43,18) = pd(43,18) + rrt(514) * density(51) 
  pd(43,51) = pd(43,51) + rrt(514) * density(18) 
  pd(51,18) = pd(51,18) - rrt(514) * density(51) 
  pd(51,51) = pd(51,51) - rrt(514) * density(18) 
  pd(29,33) = pd(29,33) + rrt(515) * density(51) 
  pd(29,51) = pd(29,51) + rrt(515) * density(33) 
  pd(33,33) = pd(33,33) - rrt(515) * density(51) 
  pd(33,51) = pd(33,51) - rrt(515) * density(33) 
  pd(43,33) = pd(43,33) + rrt(515) * density(51) 
  pd(43,51) = pd(43,51) + rrt(515) * density(33) 
  pd(51,33) = pd(51,33) - rrt(515) * density(51) 
  pd(51,51) = pd(51,51) - rrt(515) * density(33) 
  pd(21,34) = pd(21,34) + rrt(516) * density(51) 
  pd(21,51) = pd(21,51) + rrt(516) * density(34) 
  pd(34,34) = pd(34,34) - rrt(516) * density(51) 
  pd(34,51) = pd(34,51) - rrt(516) * density(34) 
  pd(43,34) = pd(43,34) + rrt(516) * density(51) 
  pd(43,51) = pd(43,51) + rrt(516) * density(34) 
  pd(51,34) = pd(51,34) - rrt(516) * density(51) 
  pd(51,51) = pd(51,51) - rrt(516) * density(34) 
  pd(40,45) = pd(40,45) + rrt(517) * density(51) 
  pd(40,51) = pd(40,51) + rrt(517) * density(45) 
  pd(43,45) = pd(43,45) + rrt(517) * density(51) 
  pd(43,51) = pd(43,51) + rrt(517) * density(45) 
  pd(45,45) = pd(45,45) - rrt(517) * density(51) 
  pd(45,51) = pd(45,51) - rrt(517) * density(45) 
  pd(51,45) = pd(51,45) - rrt(517) * density(51) 
  pd(51,51) = pd(51,51) - rrt(517) * density(45) 
  pd(41,46) = pd(41,46) + rrt(518) * density(51) 
  pd(41,51) = pd(41,51) + rrt(518) * density(46) 
  pd(43,46) = pd(43,46) + rrt(518) * density(51) 
  pd(43,51) = pd(43,51) + rrt(518) * density(46) 
  pd(46,46) = pd(46,46) - rrt(518) * density(51) 
  pd(46,51) = pd(46,51) - rrt(518) * density(46) 
  pd(51,46) = pd(51,46) - rrt(518) * density(51) 
  pd(51,51) = pd(51,51) - rrt(518) * density(46) 
  pd(42,47) = pd(42,47) + rrt(519) * density(51) 
  pd(42,51) = pd(42,51) + rrt(519) * density(47) 
  pd(43,47) = pd(43,47) + rrt(519) * density(51) 
  pd(43,51) = pd(43,51) + rrt(519) * density(47) 
  pd(47,47) = pd(47,47) - rrt(519) * density(51) 
  pd(47,51) = pd(47,51) - rrt(519) * density(47) 
  pd(51,47) = pd(51,47) - rrt(519) * density(51) 
  pd(51,51) = pd(51,51) - rrt(519) * density(47) 
  pd(14,18) = pd(14,18) + rrt(520) * density(36) * 2.0d0
  pd(14,36) = pd(14,36) + rrt(520) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(520) * density(36) 
  pd(18,36) = pd(18,36) - rrt(520) * density(18) 
  pd(29,18) = pd(29,18) + rrt(520) * density(36) 
  pd(29,36) = pd(29,36) + rrt(520) * density(18) 
  pd(36,18) = pd(36,18) - rrt(520) * density(36) 
  pd(36,36) = pd(36,36) - rrt(520) * density(18) 
  pd(01,19) = pd(01,19) + rrt(521) * density(36) 
  pd(01,36) = pd(01,36) + rrt(521) * density(19) 
  pd(14,19) = pd(14,19) + rrt(521) * density(36) 
  pd(14,36) = pd(14,36) + rrt(521) * density(19) 
  pd(19,19) = pd(19,19) - rrt(521) * density(36) 
  pd(19,36) = pd(19,36) - rrt(521) * density(19) 
  pd(29,19) = pd(29,19) + rrt(521) * density(36) 
  pd(29,36) = pd(29,36) + rrt(521) * density(19) 
  pd(36,19) = pd(36,19) - rrt(521) * density(36) 
  pd(36,36) = pd(36,36) - rrt(521) * density(19) 
  pd(01,20) = pd(01,20) + rrt(522) * density(36) * 2.0d0
  pd(01,36) = pd(01,36) + rrt(522) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(522) * density(36) 
  pd(20,36) = pd(20,36) - rrt(522) * density(20) 
  pd(29,20) = pd(29,20) + rrt(522) * density(36) 
  pd(29,36) = pd(29,36) + rrt(522) * density(20) 
  pd(36,20) = pd(36,20) - rrt(522) * density(36) 
  pd(36,36) = pd(36,36) - rrt(522) * density(20) 
  pd(29,34) = pd(29,34) + rrt(523) * density(36) * 3.0d0
  pd(29,36) = pd(29,36) + rrt(523) * density(34) * 3.0d0
  pd(34,34) = pd(34,34) - rrt(523) * density(36) 
  pd(34,36) = pd(34,36) - rrt(523) * density(34) 
  pd(36,34) = pd(36,34) - rrt(523) * density(36) 
  pd(36,36) = pd(36,36) - rrt(523) * density(34) 
  pd(21,35) = pd(21,35) + rrt(524) * density(36) * 2.0d0
  pd(21,36) = pd(21,36) + rrt(524) * density(35) * 2.0d0
  pd(29,35) = pd(29,35) + rrt(524) * density(36) 
  pd(29,36) = pd(29,36) + rrt(524) * density(35) 
  pd(35,35) = pd(35,35) - rrt(524) * density(36) 
  pd(35,36) = pd(35,36) - rrt(524) * density(35) 
  pd(36,35) = pd(36,35) - rrt(524) * density(36) 
  pd(36,36) = pd(36,36) - rrt(524) * density(35) 
  pd(14,36) = pd(14,36) + rrt(525) * density(45) 
  pd(14,45) = pd(14,45) + rrt(525) * density(36) 
  pd(29,36) = pd(29,36) + rrt(525) * density(45) * 2.0d0
  pd(29,45) = pd(29,45) + rrt(525) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(525) * density(45) 
  pd(36,45) = pd(36,45) - rrt(525) * density(36) 
  pd(45,36) = pd(45,36) - rrt(525) * density(45) 
  pd(45,45) = pd(45,45) - rrt(525) * density(36) 
  pd(01,36) = pd(01,36) + rrt(526) * density(46) 
  pd(01,46) = pd(01,46) + rrt(526) * density(36) 
  pd(29,36) = pd(29,36) + rrt(526) * density(46) * 2.0d0
  pd(29,46) = pd(29,46) + rrt(526) * density(36) * 2.0d0
  pd(36,36) = pd(36,36) - rrt(526) * density(46) 
  pd(36,46) = pd(36,46) - rrt(526) * density(36) 
  pd(46,36) = pd(46,36) - rrt(526) * density(46) 
  pd(46,46) = pd(46,46) - rrt(526) * density(36) 
  pd(14,36) = pd(14,36) + rrt(527) * density(47) 
  pd(14,47) = pd(14,47) + rrt(527) * density(36) 
  pd(21,36) = pd(21,36) + rrt(527) * density(47) 
  pd(21,47) = pd(21,47) + rrt(527) * density(36) 
  pd(29,36) = pd(29,36) + rrt(527) * density(47) 
  pd(29,47) = pd(29,47) + rrt(527) * density(36) 
  pd(36,36) = pd(36,36) - rrt(527) * density(47) 
  pd(36,47) = pd(36,47) - rrt(527) * density(36) 
  pd(47,36) = pd(47,36) - rrt(527) * density(47) 
  pd(47,47) = pd(47,47) - rrt(527) * density(36) 
  pd(01,36) = pd(01,36) + rrt(528) * density(52) 
  pd(01,52) = pd(01,52) + rrt(528) * density(36) 
  pd(21,36) = pd(21,36) + rrt(528) * density(52) 
  pd(21,52) = pd(21,52) + rrt(528) * density(36) 
  pd(29,36) = pd(29,36) + rrt(528) * density(52) 
  pd(29,52) = pd(29,52) + rrt(528) * density(36) 
  pd(36,36) = pd(36,36) - rrt(528) * density(52) 
  pd(36,52) = pd(36,52) - rrt(528) * density(36) 
  pd(52,36) = pd(52,36) - rrt(528) * density(52) 
  pd(52,52) = pd(52,52) - rrt(528) * density(36) 
  pd(14,18) = pd(14,18) + rrt(529) * density(37) * 2.0d0
  pd(14,37) = pd(14,37) + rrt(529) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(529) * density(37) 
  pd(18,37) = pd(18,37) - rrt(529) * density(18) 
  pd(21,18) = pd(21,18) + rrt(529) * density(37) 
  pd(21,37) = pd(21,37) + rrt(529) * density(18) 
  pd(37,18) = pd(37,18) - rrt(529) * density(37) 
  pd(37,37) = pd(37,37) - rrt(529) * density(18) 
  pd(01,19) = pd(01,19) + rrt(530) * density(37) 
  pd(01,37) = pd(01,37) + rrt(530) * density(19) 
  pd(14,19) = pd(14,19) + rrt(530) * density(37) 
  pd(14,37) = pd(14,37) + rrt(530) * density(19) 
  pd(19,19) = pd(19,19) - rrt(530) * density(37) 
  pd(19,37) = pd(19,37) - rrt(530) * density(19) 
  pd(21,19) = pd(21,19) + rrt(530) * density(37) 
  pd(21,37) = pd(21,37) + rrt(530) * density(19) 
  pd(37,19) = pd(37,19) - rrt(530) * density(37) 
  pd(37,37) = pd(37,37) - rrt(530) * density(19) 
  pd(01,20) = pd(01,20) + rrt(531) * density(37) * 2.0d0
  pd(01,37) = pd(01,37) + rrt(531) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(531) * density(37) 
  pd(20,37) = pd(20,37) - rrt(531) * density(20) 
  pd(21,20) = pd(21,20) + rrt(531) * density(37) 
  pd(21,37) = pd(21,37) + rrt(531) * density(20) 
  pd(37,20) = pd(37,20) - rrt(531) * density(37) 
  pd(37,37) = pd(37,37) - rrt(531) * density(20) 
  pd(21,34) = pd(21,34) + rrt(532) * density(37) 
  pd(21,37) = pd(21,37) + rrt(532) * density(34) 
  pd(29,34) = pd(29,34) + rrt(532) * density(37) * 2.0d0
  pd(29,37) = pd(29,37) + rrt(532) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(532) * density(37) 
  pd(34,37) = pd(34,37) - rrt(532) * density(34) 
  pd(37,34) = pd(37,34) - rrt(532) * density(37) 
  pd(37,37) = pd(37,37) - rrt(532) * density(34) 
  pd(21,35) = pd(21,35) + rrt(533) * density(37) * 3.0d0
  pd(21,37) = pd(21,37) + rrt(533) * density(35) * 3.0d0
  pd(35,35) = pd(35,35) - rrt(533) * density(37) 
  pd(35,37) = pd(35,37) - rrt(533) * density(35) 
  pd(37,35) = pd(37,35) - rrt(533) * density(37) 
  pd(37,37) = pd(37,37) - rrt(533) * density(35) 
  pd(14,37) = pd(14,37) + rrt(534) * density(45) 
  pd(14,45) = pd(14,45) + rrt(534) * density(37) 
  pd(21,37) = pd(21,37) + rrt(534) * density(45) 
  pd(21,45) = pd(21,45) + rrt(534) * density(37) 
  pd(29,37) = pd(29,37) + rrt(534) * density(45) 
  pd(29,45) = pd(29,45) + rrt(534) * density(37) 
  pd(37,37) = pd(37,37) - rrt(534) * density(45) 
  pd(37,45) = pd(37,45) - rrt(534) * density(37) 
  pd(45,37) = pd(45,37) - rrt(534) * density(45) 
  pd(45,45) = pd(45,45) - rrt(534) * density(37) 
  pd(01,37) = pd(01,37) + rrt(535) * density(46) 
  pd(01,46) = pd(01,46) + rrt(535) * density(37) 
  pd(21,37) = pd(21,37) + rrt(535) * density(46) 
  pd(21,46) = pd(21,46) + rrt(535) * density(37) 
  pd(29,37) = pd(29,37) + rrt(535) * density(46) 
  pd(29,46) = pd(29,46) + rrt(535) * density(37) 
  pd(37,37) = pd(37,37) - rrt(535) * density(46) 
  pd(37,46) = pd(37,46) - rrt(535) * density(37) 
  pd(46,37) = pd(46,37) - rrt(535) * density(46) 
  pd(46,46) = pd(46,46) - rrt(535) * density(37) 
  pd(14,37) = pd(14,37) + rrt(536) * density(47) 
  pd(14,47) = pd(14,47) + rrt(536) * density(37) 
  pd(21,37) = pd(21,37) + rrt(536) * density(47) * 2.0d0
  pd(21,47) = pd(21,47) + rrt(536) * density(37) * 2.0d0
  pd(37,37) = pd(37,37) - rrt(536) * density(47) 
  pd(37,47) = pd(37,47) - rrt(536) * density(37) 
  pd(47,37) = pd(47,37) - rrt(536) * density(47) 
  pd(47,47) = pd(47,47) - rrt(536) * density(37) 
  pd(01,37) = pd(01,37) + rrt(537) * density(52) 
  pd(01,52) = pd(01,52) + rrt(537) * density(37) 
  pd(21,37) = pd(21,37) + rrt(537) * density(52) * 2.0d0
  pd(21,52) = pd(21,52) + rrt(537) * density(37) * 2.0d0
  pd(37,37) = pd(37,37) - rrt(537) * density(52) 
  pd(37,52) = pd(37,52) - rrt(537) * density(37) 
  pd(52,37) = pd(52,37) - rrt(537) * density(52) 
  pd(52,52) = pd(52,52) - rrt(537) * density(37) 
  pd(14,18) = pd(14,18) + rrt(538) * density(38) * 2.0d0
  pd(14,38) = pd(14,38) + rrt(538) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(538) * density(38) 
  pd(18,38) = pd(18,38) - rrt(538) * density(18) 
  pd(32,18) = pd(32,18) + rrt(538) * density(38) 
  pd(32,38) = pd(32,38) + rrt(538) * density(18) 
  pd(38,18) = pd(38,18) - rrt(538) * density(38) 
  pd(38,38) = pd(38,38) - rrt(538) * density(18) 
  pd(01,19) = pd(01,19) + rrt(539) * density(38) 
  pd(01,38) = pd(01,38) + rrt(539) * density(19) 
  pd(14,19) = pd(14,19) + rrt(539) * density(38) 
  pd(14,38) = pd(14,38) + rrt(539) * density(19) 
  pd(19,19) = pd(19,19) - rrt(539) * density(38) 
  pd(19,38) = pd(19,38) - rrt(539) * density(19) 
  pd(32,19) = pd(32,19) + rrt(539) * density(38) 
  pd(32,38) = pd(32,38) + rrt(539) * density(19) 
  pd(38,19) = pd(38,19) - rrt(539) * density(38) 
  pd(38,38) = pd(38,38) - rrt(539) * density(19) 
  pd(01,20) = pd(01,20) + rrt(540) * density(38) * 2.0d0
  pd(01,38) = pd(01,38) + rrt(540) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(540) * density(38) 
  pd(20,38) = pd(20,38) - rrt(540) * density(20) 
  pd(32,20) = pd(32,20) + rrt(540) * density(38) 
  pd(32,38) = pd(32,38) + rrt(540) * density(20) 
  pd(38,20) = pd(38,20) - rrt(540) * density(38) 
  pd(38,38) = pd(38,38) - rrt(540) * density(20) 
  pd(29,34) = pd(29,34) + rrt(541) * density(38) * 2.0d0
  pd(29,38) = pd(29,38) + rrt(541) * density(34) * 2.0d0
  pd(32,34) = pd(32,34) + rrt(541) * density(38) 
  pd(32,38) = pd(32,38) + rrt(541) * density(34) 
  pd(34,34) = pd(34,34) - rrt(541) * density(38) 
  pd(34,38) = pd(34,38) - rrt(541) * density(34) 
  pd(38,34) = pd(38,34) - rrt(541) * density(38) 
  pd(38,38) = pd(38,38) - rrt(541) * density(34) 
  pd(21,35) = pd(21,35) + rrt(542) * density(38) * 2.0d0
  pd(21,38) = pd(21,38) + rrt(542) * density(35) * 2.0d0
  pd(32,35) = pd(32,35) + rrt(542) * density(38) 
  pd(32,38) = pd(32,38) + rrt(542) * density(35) 
  pd(35,35) = pd(35,35) - rrt(542) * density(38) 
  pd(35,38) = pd(35,38) - rrt(542) * density(35) 
  pd(38,35) = pd(38,35) - rrt(542) * density(38) 
  pd(38,38) = pd(38,38) - rrt(542) * density(35) 
  pd(14,38) = pd(14,38) + rrt(543) * density(45) 
  pd(14,45) = pd(14,45) + rrt(543) * density(38) 
  pd(29,38) = pd(29,38) + rrt(543) * density(45) 
  pd(29,45) = pd(29,45) + rrt(543) * density(38) 
  pd(32,38) = pd(32,38) + rrt(543) * density(45) 
  pd(32,45) = pd(32,45) + rrt(543) * density(38) 
  pd(38,38) = pd(38,38) - rrt(543) * density(45) 
  pd(38,45) = pd(38,45) - rrt(543) * density(38) 
  pd(45,38) = pd(45,38) - rrt(543) * density(45) 
  pd(45,45) = pd(45,45) - rrt(543) * density(38) 
  pd(01,38) = pd(01,38) + rrt(544) * density(46) 
  pd(01,46) = pd(01,46) + rrt(544) * density(38) 
  pd(29,38) = pd(29,38) + rrt(544) * density(46) 
  pd(29,46) = pd(29,46) + rrt(544) * density(38) 
  pd(32,38) = pd(32,38) + rrt(544) * density(46) 
  pd(32,46) = pd(32,46) + rrt(544) * density(38) 
  pd(38,38) = pd(38,38) - rrt(544) * density(46) 
  pd(38,46) = pd(38,46) - rrt(544) * density(38) 
  pd(46,38) = pd(46,38) - rrt(544) * density(46) 
  pd(46,46) = pd(46,46) - rrt(544) * density(38) 
  pd(14,38) = pd(14,38) + rrt(545) * density(47) 
  pd(14,47) = pd(14,47) + rrt(545) * density(38) 
  pd(21,38) = pd(21,38) + rrt(545) * density(47) 
  pd(21,47) = pd(21,47) + rrt(545) * density(38) 
  pd(32,38) = pd(32,38) + rrt(545) * density(47) 
  pd(32,47) = pd(32,47) + rrt(545) * density(38) 
  pd(38,38) = pd(38,38) - rrt(545) * density(47) 
  pd(38,47) = pd(38,47) - rrt(545) * density(38) 
  pd(47,38) = pd(47,38) - rrt(545) * density(47) 
  pd(47,47) = pd(47,47) - rrt(545) * density(38) 
  pd(01,38) = pd(01,38) + rrt(546) * density(52) 
  pd(01,52) = pd(01,52) + rrt(546) * density(38) 
  pd(21,38) = pd(21,38) + rrt(546) * density(52) 
  pd(21,52) = pd(21,52) + rrt(546) * density(38) 
  pd(32,38) = pd(32,38) + rrt(546) * density(52) 
  pd(32,52) = pd(32,52) + rrt(546) * density(38) 
  pd(38,38) = pd(38,38) - rrt(546) * density(52) 
  pd(38,52) = pd(38,52) - rrt(546) * density(38) 
  pd(52,38) = pd(52,38) - rrt(546) * density(52) 
  pd(52,52) = pd(52,52) - rrt(546) * density(38) 
  pd(14,18) = pd(14,18) + rrt(547) * density(48) * 2.0d0
  pd(14,48) = pd(14,48) + rrt(547) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(547) * density(48) 
  pd(18,48) = pd(18,48) - rrt(547) * density(18) 
  pd(40,18) = pd(40,18) + rrt(547) * density(48) 
  pd(40,48) = pd(40,48) + rrt(547) * density(18) 
  pd(48,18) = pd(48,18) - rrt(547) * density(48) 
  pd(48,48) = pd(48,48) - rrt(547) * density(18) 
  pd(01,19) = pd(01,19) + rrt(548) * density(48) 
  pd(01,48) = pd(01,48) + rrt(548) * density(19) 
  pd(14,19) = pd(14,19) + rrt(548) * density(48) 
  pd(14,48) = pd(14,48) + rrt(548) * density(19) 
  pd(19,19) = pd(19,19) - rrt(548) * density(48) 
  pd(19,48) = pd(19,48) - rrt(548) * density(19) 
  pd(40,19) = pd(40,19) + rrt(548) * density(48) 
  pd(40,48) = pd(40,48) + rrt(548) * density(19) 
  pd(48,19) = pd(48,19) - rrt(548) * density(48) 
  pd(48,48) = pd(48,48) - rrt(548) * density(19) 
  pd(01,20) = pd(01,20) + rrt(549) * density(48) * 2.0d0
  pd(01,48) = pd(01,48) + rrt(549) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(549) * density(48) 
  pd(20,48) = pd(20,48) - rrt(549) * density(20) 
  pd(40,20) = pd(40,20) + rrt(549) * density(48) 
  pd(40,48) = pd(40,48) + rrt(549) * density(20) 
  pd(48,20) = pd(48,20) - rrt(549) * density(48) 
  pd(48,48) = pd(48,48) - rrt(549) * density(20) 
  pd(29,34) = pd(29,34) + rrt(550) * density(48) * 2.0d0
  pd(29,48) = pd(29,48) + rrt(550) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(550) * density(48) 
  pd(34,48) = pd(34,48) - rrt(550) * density(34) 
  pd(40,34) = pd(40,34) + rrt(550) * density(48) 
  pd(40,48) = pd(40,48) + rrt(550) * density(34) 
  pd(48,34) = pd(48,34) - rrt(550) * density(48) 
  pd(48,48) = pd(48,48) - rrt(550) * density(34) 
  pd(21,35) = pd(21,35) + rrt(551) * density(48) * 2.0d0
  pd(21,48) = pd(21,48) + rrt(551) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(551) * density(48) 
  pd(35,48) = pd(35,48) - rrt(551) * density(35) 
  pd(40,35) = pd(40,35) + rrt(551) * density(48) 
  pd(40,48) = pd(40,48) + rrt(551) * density(35) 
  pd(48,35) = pd(48,35) - rrt(551) * density(48) 
  pd(48,48) = pd(48,48) - rrt(551) * density(35) 
  pd(14,45) = pd(14,45) + rrt(552) * density(48) 
  pd(14,48) = pd(14,48) + rrt(552) * density(45) 
  pd(29,45) = pd(29,45) + rrt(552) * density(48) 
  pd(29,48) = pd(29,48) + rrt(552) * density(45) 
  pd(40,45) = pd(40,45) + rrt(552) * density(48) 
  pd(40,48) = pd(40,48) + rrt(552) * density(45) 
  pd(45,45) = pd(45,45) - rrt(552) * density(48) 
  pd(45,48) = pd(45,48) - rrt(552) * density(45) 
  pd(48,45) = pd(48,45) - rrt(552) * density(48) 
  pd(48,48) = pd(48,48) - rrt(552) * density(45) 
  pd(01,46) = pd(01,46) + rrt(553) * density(48) 
  pd(01,48) = pd(01,48) + rrt(553) * density(46) 
  pd(29,46) = pd(29,46) + rrt(553) * density(48) 
  pd(29,48) = pd(29,48) + rrt(553) * density(46) 
  pd(40,46) = pd(40,46) + rrt(553) * density(48) 
  pd(40,48) = pd(40,48) + rrt(553) * density(46) 
  pd(46,46) = pd(46,46) - rrt(553) * density(48) 
  pd(46,48) = pd(46,48) - rrt(553) * density(46) 
  pd(48,46) = pd(48,46) - rrt(553) * density(48) 
  pd(48,48) = pd(48,48) - rrt(553) * density(46) 
  pd(14,47) = pd(14,47) + rrt(554) * density(48) 
  pd(14,48) = pd(14,48) + rrt(554) * density(47) 
  pd(21,47) = pd(21,47) + rrt(554) * density(48) 
  pd(21,48) = pd(21,48) + rrt(554) * density(47) 
  pd(40,47) = pd(40,47) + rrt(554) * density(48) 
  pd(40,48) = pd(40,48) + rrt(554) * density(47) 
  pd(47,47) = pd(47,47) - rrt(554) * density(48) 
  pd(47,48) = pd(47,48) - rrt(554) * density(47) 
  pd(48,47) = pd(48,47) - rrt(554) * density(48) 
  pd(48,48) = pd(48,48) - rrt(554) * density(47) 
  pd(01,48) = pd(01,48) + rrt(555) * density(52) 
  pd(01,52) = pd(01,52) + rrt(555) * density(48) 
  pd(21,48) = pd(21,48) + rrt(555) * density(52) 
  pd(21,52) = pd(21,52) + rrt(555) * density(48) 
  pd(40,48) = pd(40,48) + rrt(555) * density(52) 
  pd(40,52) = pd(40,52) + rrt(555) * density(48) 
  pd(48,48) = pd(48,48) - rrt(555) * density(52) 
  pd(48,52) = pd(48,52) - rrt(555) * density(48) 
  pd(52,48) = pd(52,48) - rrt(555) * density(52) 
  pd(52,52) = pd(52,52) - rrt(555) * density(48) 
  pd(14,18) = pd(14,18) + rrt(556) * density(49) * 2.0d0
  pd(14,49) = pd(14,49) + rrt(556) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(556) * density(49) 
  pd(18,49) = pd(18,49) - rrt(556) * density(18) 
  pd(41,18) = pd(41,18) + rrt(556) * density(49) 
  pd(41,49) = pd(41,49) + rrt(556) * density(18) 
  pd(49,18) = pd(49,18) - rrt(556) * density(49) 
  pd(49,49) = pd(49,49) - rrt(556) * density(18) 
  pd(01,19) = pd(01,19) + rrt(557) * density(49) 
  pd(01,49) = pd(01,49) + rrt(557) * density(19) 
  pd(14,19) = pd(14,19) + rrt(557) * density(49) 
  pd(14,49) = pd(14,49) + rrt(557) * density(19) 
  pd(19,19) = pd(19,19) - rrt(557) * density(49) 
  pd(19,49) = pd(19,49) - rrt(557) * density(19) 
  pd(41,19) = pd(41,19) + rrt(557) * density(49) 
  pd(41,49) = pd(41,49) + rrt(557) * density(19) 
  pd(49,19) = pd(49,19) - rrt(557) * density(49) 
  pd(49,49) = pd(49,49) - rrt(557) * density(19) 
  pd(01,20) = pd(01,20) + rrt(558) * density(49) * 2.0d0
  pd(01,49) = pd(01,49) + rrt(558) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(558) * density(49) 
  pd(20,49) = pd(20,49) - rrt(558) * density(20) 
  pd(41,20) = pd(41,20) + rrt(558) * density(49) 
  pd(41,49) = pd(41,49) + rrt(558) * density(20) 
  pd(49,20) = pd(49,20) - rrt(558) * density(49) 
  pd(49,49) = pd(49,49) - rrt(558) * density(20) 
  pd(29,34) = pd(29,34) + rrt(559) * density(49) * 2.0d0
  pd(29,49) = pd(29,49) + rrt(559) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(559) * density(49) 
  pd(34,49) = pd(34,49) - rrt(559) * density(34) 
  pd(41,34) = pd(41,34) + rrt(559) * density(49) 
  pd(41,49) = pd(41,49) + rrt(559) * density(34) 
  pd(49,34) = pd(49,34) - rrt(559) * density(49) 
  pd(49,49) = pd(49,49) - rrt(559) * density(34) 
  pd(21,35) = pd(21,35) + rrt(560) * density(49) * 2.0d0
  pd(21,49) = pd(21,49) + rrt(560) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(560) * density(49) 
  pd(35,49) = pd(35,49) - rrt(560) * density(35) 
  pd(41,35) = pd(41,35) + rrt(560) * density(49) 
  pd(41,49) = pd(41,49) + rrt(560) * density(35) 
  pd(49,35) = pd(49,35) - rrt(560) * density(49) 
  pd(49,49) = pd(49,49) - rrt(560) * density(35) 
  pd(14,45) = pd(14,45) + rrt(561) * density(49) 
  pd(14,49) = pd(14,49) + rrt(561) * density(45) 
  pd(29,45) = pd(29,45) + rrt(561) * density(49) 
  pd(29,49) = pd(29,49) + rrt(561) * density(45) 
  pd(41,45) = pd(41,45) + rrt(561) * density(49) 
  pd(41,49) = pd(41,49) + rrt(561) * density(45) 
  pd(45,45) = pd(45,45) - rrt(561) * density(49) 
  pd(45,49) = pd(45,49) - rrt(561) * density(45) 
  pd(49,45) = pd(49,45) - rrt(561) * density(49) 
  pd(49,49) = pd(49,49) - rrt(561) * density(45) 
  pd(01,46) = pd(01,46) + rrt(562) * density(49) 
  pd(01,49) = pd(01,49) + rrt(562) * density(46) 
  pd(29,46) = pd(29,46) + rrt(562) * density(49) 
  pd(29,49) = pd(29,49) + rrt(562) * density(46) 
  pd(41,46) = pd(41,46) + rrt(562) * density(49) 
  pd(41,49) = pd(41,49) + rrt(562) * density(46) 
  pd(46,46) = pd(46,46) - rrt(562) * density(49) 
  pd(46,49) = pd(46,49) - rrt(562) * density(46) 
  pd(49,46) = pd(49,46) - rrt(562) * density(49) 
  pd(49,49) = pd(49,49) - rrt(562) * density(46) 
  pd(14,47) = pd(14,47) + rrt(563) * density(49) 
  pd(14,49) = pd(14,49) + rrt(563) * density(47) 
  pd(21,47) = pd(21,47) + rrt(563) * density(49) 
  pd(21,49) = pd(21,49) + rrt(563) * density(47) 
  pd(41,47) = pd(41,47) + rrt(563) * density(49) 
  pd(41,49) = pd(41,49) + rrt(563) * density(47) 
  pd(47,47) = pd(47,47) - rrt(563) * density(49) 
  pd(47,49) = pd(47,49) - rrt(563) * density(47) 
  pd(49,47) = pd(49,47) - rrt(563) * density(49) 
  pd(49,49) = pd(49,49) - rrt(563) * density(47) 
  pd(01,49) = pd(01,49) + rrt(564) * density(52) 
  pd(01,52) = pd(01,52) + rrt(564) * density(49) 
  pd(21,49) = pd(21,49) + rrt(564) * density(52) 
  pd(21,52) = pd(21,52) + rrt(564) * density(49) 
  pd(41,49) = pd(41,49) + rrt(564) * density(52) 
  pd(41,52) = pd(41,52) + rrt(564) * density(49) 
  pd(49,49) = pd(49,49) - rrt(564) * density(52) 
  pd(49,52) = pd(49,52) - rrt(564) * density(49) 
  pd(52,49) = pd(52,49) - rrt(564) * density(52) 
  pd(52,52) = pd(52,52) - rrt(564) * density(49) 
  pd(14,18) = pd(14,18) + rrt(565) * density(50) * 2.0d0
  pd(14,50) = pd(14,50) + rrt(565) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(565) * density(50) 
  pd(18,50) = pd(18,50) - rrt(565) * density(18) 
  pd(42,18) = pd(42,18) + rrt(565) * density(50) 
  pd(42,50) = pd(42,50) + rrt(565) * density(18) 
  pd(50,18) = pd(50,18) - rrt(565) * density(50) 
  pd(50,50) = pd(50,50) - rrt(565) * density(18) 
  pd(01,19) = pd(01,19) + rrt(566) * density(50) 
  pd(01,50) = pd(01,50) + rrt(566) * density(19) 
  pd(14,19) = pd(14,19) + rrt(566) * density(50) 
  pd(14,50) = pd(14,50) + rrt(566) * density(19) 
  pd(19,19) = pd(19,19) - rrt(566) * density(50) 
  pd(19,50) = pd(19,50) - rrt(566) * density(19) 
  pd(42,19) = pd(42,19) + rrt(566) * density(50) 
  pd(42,50) = pd(42,50) + rrt(566) * density(19) 
  pd(50,19) = pd(50,19) - rrt(566) * density(50) 
  pd(50,50) = pd(50,50) - rrt(566) * density(19) 
  pd(01,20) = pd(01,20) + rrt(567) * density(50) * 2.0d0
  pd(01,50) = pd(01,50) + rrt(567) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(567) * density(50) 
  pd(20,50) = pd(20,50) - rrt(567) * density(20) 
  pd(42,20) = pd(42,20) + rrt(567) * density(50) 
  pd(42,50) = pd(42,50) + rrt(567) * density(20) 
  pd(50,20) = pd(50,20) - rrt(567) * density(50) 
  pd(50,50) = pd(50,50) - rrt(567) * density(20) 
  pd(29,34) = pd(29,34) + rrt(568) * density(50) * 2.0d0
  pd(29,50) = pd(29,50) + rrt(568) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(568) * density(50) 
  pd(34,50) = pd(34,50) - rrt(568) * density(34) 
  pd(42,34) = pd(42,34) + rrt(568) * density(50) 
  pd(42,50) = pd(42,50) + rrt(568) * density(34) 
  pd(50,34) = pd(50,34) - rrt(568) * density(50) 
  pd(50,50) = pd(50,50) - rrt(568) * density(34) 
  pd(21,35) = pd(21,35) + rrt(569) * density(50) * 2.0d0
  pd(21,50) = pd(21,50) + rrt(569) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(569) * density(50) 
  pd(35,50) = pd(35,50) - rrt(569) * density(35) 
  pd(42,35) = pd(42,35) + rrt(569) * density(50) 
  pd(42,50) = pd(42,50) + rrt(569) * density(35) 
  pd(50,35) = pd(50,35) - rrt(569) * density(50) 
  pd(50,50) = pd(50,50) - rrt(569) * density(35) 
  pd(14,45) = pd(14,45) + rrt(570) * density(50) 
  pd(14,50) = pd(14,50) + rrt(570) * density(45) 
  pd(29,45) = pd(29,45) + rrt(570) * density(50) 
  pd(29,50) = pd(29,50) + rrt(570) * density(45) 
  pd(42,45) = pd(42,45) + rrt(570) * density(50) 
  pd(42,50) = pd(42,50) + rrt(570) * density(45) 
  pd(45,45) = pd(45,45) - rrt(570) * density(50) 
  pd(45,50) = pd(45,50) - rrt(570) * density(45) 
  pd(50,45) = pd(50,45) - rrt(570) * density(50) 
  pd(50,50) = pd(50,50) - rrt(570) * density(45) 
  pd(01,46) = pd(01,46) + rrt(571) * density(50) 
  pd(01,50) = pd(01,50) + rrt(571) * density(46) 
  pd(29,46) = pd(29,46) + rrt(571) * density(50) 
  pd(29,50) = pd(29,50) + rrt(571) * density(46) 
  pd(42,46) = pd(42,46) + rrt(571) * density(50) 
  pd(42,50) = pd(42,50) + rrt(571) * density(46) 
  pd(46,46) = pd(46,46) - rrt(571) * density(50) 
  pd(46,50) = pd(46,50) - rrt(571) * density(46) 
  pd(50,46) = pd(50,46) - rrt(571) * density(50) 
  pd(50,50) = pd(50,50) - rrt(571) * density(46) 
  pd(14,47) = pd(14,47) + rrt(572) * density(50) 
  pd(14,50) = pd(14,50) + rrt(572) * density(47) 
  pd(21,47) = pd(21,47) + rrt(572) * density(50) 
  pd(21,50) = pd(21,50) + rrt(572) * density(47) 
  pd(42,47) = pd(42,47) + rrt(572) * density(50) 
  pd(42,50) = pd(42,50) + rrt(572) * density(47) 
  pd(47,47) = pd(47,47) - rrt(572) * density(50) 
  pd(47,50) = pd(47,50) - rrt(572) * density(47) 
  pd(50,47) = pd(50,47) - rrt(572) * density(50) 
  pd(50,50) = pd(50,50) - rrt(572) * density(47) 
  pd(01,50) = pd(01,50) + rrt(573) * density(52) 
  pd(01,52) = pd(01,52) + rrt(573) * density(50) 
  pd(21,50) = pd(21,50) + rrt(573) * density(52) 
  pd(21,52) = pd(21,52) + rrt(573) * density(50) 
  pd(42,50) = pd(42,50) + rrt(573) * density(52) 
  pd(42,52) = pd(42,52) + rrt(573) * density(50) 
  pd(50,50) = pd(50,50) - rrt(573) * density(52) 
  pd(50,52) = pd(50,52) - rrt(573) * density(50) 
  pd(52,50) = pd(52,50) - rrt(573) * density(52) 
  pd(52,52) = pd(52,52) - rrt(573) * density(50) 
  pd(14,18) = pd(14,18) + rrt(574) * density(51) * 2.0d0
  pd(14,51) = pd(14,51) + rrt(574) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(574) * density(51) 
  pd(18,51) = pd(18,51) - rrt(574) * density(18) 
  pd(43,18) = pd(43,18) + rrt(574) * density(51) 
  pd(43,51) = pd(43,51) + rrt(574) * density(18) 
  pd(51,18) = pd(51,18) - rrt(574) * density(51) 
  pd(51,51) = pd(51,51) - rrt(574) * density(18) 
  pd(01,19) = pd(01,19) + rrt(575) * density(51) 
  pd(01,51) = pd(01,51) + rrt(575) * density(19) 
  pd(14,19) = pd(14,19) + rrt(575) * density(51) 
  pd(14,51) = pd(14,51) + rrt(575) * density(19) 
  pd(19,19) = pd(19,19) - rrt(575) * density(51) 
  pd(19,51) = pd(19,51) - rrt(575) * density(19) 
  pd(43,19) = pd(43,19) + rrt(575) * density(51) 
  pd(43,51) = pd(43,51) + rrt(575) * density(19) 
  pd(51,19) = pd(51,19) - rrt(575) * density(51) 
  pd(51,51) = pd(51,51) - rrt(575) * density(19) 
  pd(01,20) = pd(01,20) + rrt(576) * density(51) * 2.0d0
  pd(01,51) = pd(01,51) + rrt(576) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(576) * density(51) 
  pd(20,51) = pd(20,51) - rrt(576) * density(20) 
  pd(43,20) = pd(43,20) + rrt(576) * density(51) 
  pd(43,51) = pd(43,51) + rrt(576) * density(20) 
  pd(51,20) = pd(51,20) - rrt(576) * density(51) 
  pd(51,51) = pd(51,51) - rrt(576) * density(20) 
  pd(29,34) = pd(29,34) + rrt(577) * density(51) * 2.0d0
  pd(29,51) = pd(29,51) + rrt(577) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(577) * density(51) 
  pd(34,51) = pd(34,51) - rrt(577) * density(34) 
  pd(43,34) = pd(43,34) + rrt(577) * density(51) 
  pd(43,51) = pd(43,51) + rrt(577) * density(34) 
  pd(51,34) = pd(51,34) - rrt(577) * density(51) 
  pd(51,51) = pd(51,51) - rrt(577) * density(34) 
  pd(21,35) = pd(21,35) + rrt(578) * density(51) * 2.0d0
  pd(21,51) = pd(21,51) + rrt(578) * density(35) * 2.0d0
  pd(35,35) = pd(35,35) - rrt(578) * density(51) 
  pd(35,51) = pd(35,51) - rrt(578) * density(35) 
  pd(43,35) = pd(43,35) + rrt(578) * density(51) 
  pd(43,51) = pd(43,51) + rrt(578) * density(35) 
  pd(51,35) = pd(51,35) - rrt(578) * density(51) 
  pd(51,51) = pd(51,51) - rrt(578) * density(35) 
  pd(14,45) = pd(14,45) + rrt(579) * density(51) 
  pd(14,51) = pd(14,51) + rrt(579) * density(45) 
  pd(29,45) = pd(29,45) + rrt(579) * density(51) 
  pd(29,51) = pd(29,51) + rrt(579) * density(45) 
  pd(43,45) = pd(43,45) + rrt(579) * density(51) 
  pd(43,51) = pd(43,51) + rrt(579) * density(45) 
  pd(45,45) = pd(45,45) - rrt(579) * density(51) 
  pd(45,51) = pd(45,51) - rrt(579) * density(45) 
  pd(51,45) = pd(51,45) - rrt(579) * density(51) 
  pd(51,51) = pd(51,51) - rrt(579) * density(45) 
  pd(01,46) = pd(01,46) + rrt(580) * density(51) 
  pd(01,51) = pd(01,51) + rrt(580) * density(46) 
  pd(29,46) = pd(29,46) + rrt(580) * density(51) 
  pd(29,51) = pd(29,51) + rrt(580) * density(46) 
  pd(43,46) = pd(43,46) + rrt(580) * density(51) 
  pd(43,51) = pd(43,51) + rrt(580) * density(46) 
  pd(46,46) = pd(46,46) - rrt(580) * density(51) 
  pd(46,51) = pd(46,51) - rrt(580) * density(46) 
  pd(51,46) = pd(51,46) - rrt(580) * density(51) 
  pd(51,51) = pd(51,51) - rrt(580) * density(46) 
  pd(14,47) = pd(14,47) + rrt(581) * density(51) 
  pd(14,51) = pd(14,51) + rrt(581) * density(47) 
  pd(21,47) = pd(21,47) + rrt(581) * density(51) 
  pd(21,51) = pd(21,51) + rrt(581) * density(47) 
  pd(43,47) = pd(43,47) + rrt(581) * density(51) 
  pd(43,51) = pd(43,51) + rrt(581) * density(47) 
  pd(47,47) = pd(47,47) - rrt(581) * density(51) 
  pd(47,51) = pd(47,51) - rrt(581) * density(47) 
  pd(51,47) = pd(51,47) - rrt(581) * density(51) 
  pd(51,51) = pd(51,51) - rrt(581) * density(47) 
  pd(01,51) = pd(01,51) + rrt(582) * density(52) 
  pd(01,52) = pd(01,52) + rrt(582) * density(51) 
  pd(21,51) = pd(21,51) + rrt(582) * density(52) 
  pd(21,52) = pd(21,52) + rrt(582) * density(51) 
  pd(43,51) = pd(43,51) + rrt(582) * density(52) 
  pd(43,52) = pd(43,52) + rrt(582) * density(51) 
  pd(51,51) = pd(51,51) - rrt(582) * density(52) 
  pd(51,52) = pd(51,52) - rrt(582) * density(51) 
  pd(52,51) = pd(52,51) - rrt(582) * density(52) 
  pd(52,52) = pd(52,52) - rrt(582) * density(51) 
  pd(14,17) = pd(14,17) + rrt(583) * density(39) 
  pd(14,39) = pd(14,39) + rrt(583) * density(17) 
  pd(17,17) = pd(17,17) - rrt(583) * density(39) 
  pd(17,39) = pd(17,39) - rrt(583) * density(17) 
  pd(21,17) = pd(21,17) + rrt(583) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(583) * density(17) * 2.0d0
  pd(39,17) = pd(39,17) - rrt(583) * density(39) 
  pd(39,39) = pd(39,39) - rrt(583) * density(17) 
  pd(01,18) = pd(01,18) + rrt(584) * density(39) 
  pd(01,39) = pd(01,39) + rrt(584) * density(18) 
  pd(18,18) = pd(18,18) - rrt(584) * density(39) 
  pd(18,39) = pd(18,39) - rrt(584) * density(18) 
  pd(21,18) = pd(21,18) + rrt(584) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(584) * density(18) * 2.0d0
  pd(39,18) = pd(39,18) - rrt(584) * density(39) 
  pd(39,39) = pd(39,39) - rrt(584) * density(18) 
  pd(21,33) = pd(21,33) + rrt(585) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(585) * density(33) * 2.0d0
  pd(29,33) = pd(29,33) + rrt(585) * density(39) 
  pd(29,39) = pd(29,39) + rrt(585) * density(33) 
  pd(33,33) = pd(33,33) - rrt(585) * density(39) 
  pd(33,39) = pd(33,39) - rrt(585) * density(33) 
  pd(39,33) = pd(39,33) - rrt(585) * density(39) 
  pd(39,39) = pd(39,39) - rrt(585) * density(33) 
  pd(21,34) = pd(21,34) + rrt(586) * density(39) * 3.0d0
  pd(21,39) = pd(21,39) + rrt(586) * density(34) * 3.0d0
  pd(34,34) = pd(34,34) - rrt(586) * density(39) 
  pd(34,39) = pd(34,39) - rrt(586) * density(34) 
  pd(39,34) = pd(39,34) - rrt(586) * density(39) 
  pd(39,39) = pd(39,39) - rrt(586) * density(34) 
  pd(21,39) = pd(21,39) + rrt(587) * density(45) * 2.0d0
  pd(21,45) = pd(21,45) + rrt(587) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(587) * density(45) 
  pd(39,45) = pd(39,45) - rrt(587) * density(39) 
  pd(40,39) = pd(40,39) + rrt(587) * density(45) 
  pd(40,45) = pd(40,45) + rrt(587) * density(39) 
  pd(45,39) = pd(45,39) - rrt(587) * density(45) 
  pd(45,45) = pd(45,45) - rrt(587) * density(39) 
  pd(21,39) = pd(21,39) + rrt(588) * density(46) * 2.0d0
  pd(21,46) = pd(21,46) + rrt(588) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(588) * density(46) 
  pd(39,46) = pd(39,46) - rrt(588) * density(39) 
  pd(41,39) = pd(41,39) + rrt(588) * density(46) 
  pd(41,46) = pd(41,46) + rrt(588) * density(39) 
  pd(46,39) = pd(46,39) - rrt(588) * density(46) 
  pd(46,46) = pd(46,46) - rrt(588) * density(39) 
  pd(21,39) = pd(21,39) + rrt(589) * density(47) * 2.0d0
  pd(21,47) = pd(21,47) + rrt(589) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(589) * density(47) 
  pd(39,47) = pd(39,47) - rrt(589) * density(39) 
  pd(42,39) = pd(42,39) + rrt(589) * density(47) 
  pd(42,47) = pd(42,47) + rrt(589) * density(39) 
  pd(47,39) = pd(47,39) - rrt(589) * density(47) 
  pd(47,47) = pd(47,47) - rrt(589) * density(39) 
  pd(01,19) = pd(01,19) + rrt(590) * density(39) 
  pd(01,39) = pd(01,39) + rrt(590) * density(19) 
  pd(14,19) = pd(14,19) + rrt(590) * density(39) 
  pd(14,39) = pd(14,39) + rrt(590) * density(19) 
  pd(19,19) = pd(19,19) - rrt(590) * density(39) 
  pd(19,39) = pd(19,39) - rrt(590) * density(19) 
  pd(21,19) = pd(21,19) + rrt(590) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(590) * density(19) * 2.0d0
  pd(39,19) = pd(39,19) - rrt(590) * density(39) 
  pd(39,39) = pd(39,39) - rrt(590) * density(19) 
  pd(01,20) = pd(01,20) + rrt(591) * density(39) * 2.0d0
  pd(01,39) = pd(01,39) + rrt(591) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(591) * density(39) 
  pd(20,39) = pd(20,39) - rrt(591) * density(20) 
  pd(21,20) = pd(21,20) + rrt(591) * density(39) * 2.0d0
  pd(21,39) = pd(21,39) + rrt(591) * density(20) * 2.0d0
  pd(39,20) = pd(39,20) - rrt(591) * density(39) 
  pd(39,39) = pd(39,39) - rrt(591) * density(20) 
  pd(21,35) = pd(21,35) + rrt(592) * density(39) * 4.0d0
  pd(21,39) = pd(21,39) + rrt(592) * density(35) * 4.0d0
  pd(35,35) = pd(35,35) - rrt(592) * density(39) 
  pd(35,39) = pd(35,39) - rrt(592) * density(35) 
  pd(39,35) = pd(39,35) - rrt(592) * density(39) 
  pd(39,39) = pd(39,39) - rrt(592) * density(35) 
  pd(01,39) = pd(01,39) + rrt(593) * density(52) 
  pd(01,52) = pd(01,52) + rrt(593) * density(39) 
  pd(21,39) = pd(21,39) + rrt(593) * density(52) * 3.0d0
  pd(21,52) = pd(21,52) + rrt(593) * density(39) * 3.0d0
  pd(39,39) = pd(39,39) - rrt(593) * density(52) 
  pd(39,52) = pd(39,52) - rrt(593) * density(39) 
  pd(52,39) = pd(52,39) - rrt(593) * density(52) 
  pd(52,52) = pd(52,52) - rrt(593) * density(39) 
  pd(14,17) = pd(14,17) + rrt(594) * density(36) 
  pd(14,36) = pd(14,36) + rrt(594) * density(17) 
  pd(17,17) = pd(17,17) - rrt(594) * density(36) 
  pd(17,36) = pd(17,36) - rrt(594) * density(17) 
  pd(29,17) = pd(29,17) + rrt(594) * density(36) 
  pd(29,36) = pd(29,36) + rrt(594) * density(17) 
  pd(36,17) = pd(36,17) - rrt(594) * density(36) 
  pd(36,36) = pd(36,36) - rrt(594) * density(17) 
  pd(01,18) = pd(01,18) + rrt(595) * density(36) 
  pd(01,36) = pd(01,36) + rrt(595) * density(18) 
  pd(18,18) = pd(18,18) - rrt(595) * density(36) 
  pd(18,36) = pd(18,36) - rrt(595) * density(18) 
  pd(29,18) = pd(29,18) + rrt(595) * density(36) 
  pd(29,36) = pd(29,36) + rrt(595) * density(18) 
  pd(36,18) = pd(36,18) - rrt(595) * density(36) 
  pd(36,36) = pd(36,36) - rrt(595) * density(18) 
  pd(29,33) = pd(29,33) + rrt(596) * density(36) * 2.0d0
  pd(29,36) = pd(29,36) + rrt(596) * density(33) * 2.0d0
  pd(33,33) = pd(33,33) - rrt(596) * density(36) 
  pd(33,36) = pd(33,36) - rrt(596) * density(33) 
  pd(36,33) = pd(36,33) - rrt(596) * density(36) 
  pd(36,36) = pd(36,36) - rrt(596) * density(33) 
  pd(21,34) = pd(21,34) + rrt(597) * density(36) 
  pd(21,36) = pd(21,36) + rrt(597) * density(34) 
  pd(29,34) = pd(29,34) + rrt(597) * density(36) 
  pd(29,36) = pd(29,36) + rrt(597) * density(34) 
  pd(34,34) = pd(34,34) - rrt(597) * density(36) 
  pd(34,36) = pd(34,36) - rrt(597) * density(34) 
  pd(36,34) = pd(36,34) - rrt(597) * density(36) 
  pd(36,36) = pd(36,36) - rrt(597) * density(34) 
  pd(29,36) = pd(29,36) + rrt(598) * density(45) 
  pd(29,45) = pd(29,45) + rrt(598) * density(36) 
  pd(36,36) = pd(36,36) - rrt(598) * density(45) 
  pd(36,45) = pd(36,45) - rrt(598) * density(36) 
  pd(40,36) = pd(40,36) + rrt(598) * density(45) 
  pd(40,45) = pd(40,45) + rrt(598) * density(36) 
  pd(45,36) = pd(45,36) - rrt(598) * density(45) 
  pd(45,45) = pd(45,45) - rrt(598) * density(36) 
  pd(14,17) = pd(14,17) + rrt(599) * density(37) 
  pd(14,37) = pd(14,37) + rrt(599) * density(17) 
  pd(17,17) = pd(17,17) - rrt(599) * density(37) 
  pd(17,37) = pd(17,37) - rrt(599) * density(17) 
  pd(21,17) = pd(21,17) + rrt(599) * density(37) 
  pd(21,37) = pd(21,37) + rrt(599) * density(17) 
  pd(37,17) = pd(37,17) - rrt(599) * density(37) 
  pd(37,37) = pd(37,37) - rrt(599) * density(17) 
  pd(01,18) = pd(01,18) + rrt(600) * density(37) 
  pd(01,37) = pd(01,37) + rrt(600) * density(18) 
  pd(18,18) = pd(18,18) - rrt(600) * density(37) 
  pd(18,37) = pd(18,37) - rrt(600) * density(18) 
  pd(21,18) = pd(21,18) + rrt(600) * density(37) 
  pd(21,37) = pd(21,37) + rrt(600) * density(18) 
  pd(37,18) = pd(37,18) - rrt(600) * density(37) 
  pd(37,37) = pd(37,37) - rrt(600) * density(18) 
  pd(21,33) = pd(21,33) + rrt(601) * density(37) 
  pd(21,37) = pd(21,37) + rrt(601) * density(33) 
  pd(29,33) = pd(29,33) + rrt(601) * density(37) 
  pd(29,37) = pd(29,37) + rrt(601) * density(33) 
  pd(33,33) = pd(33,33) - rrt(601) * density(37) 
  pd(33,37) = pd(33,37) - rrt(601) * density(33) 
  pd(37,33) = pd(37,33) - rrt(601) * density(37) 
  pd(37,37) = pd(37,37) - rrt(601) * density(33) 
  pd(21,34) = pd(21,34) + rrt(602) * density(37) * 2.0d0
  pd(21,37) = pd(21,37) + rrt(602) * density(34) * 2.0d0
  pd(34,34) = pd(34,34) - rrt(602) * density(37) 
  pd(34,37) = pd(34,37) - rrt(602) * density(34) 
  pd(37,34) = pd(37,34) - rrt(602) * density(37) 
  pd(37,37) = pd(37,37) - rrt(602) * density(34) 
  pd(21,37) = pd(21,37) + rrt(603) * density(45) 
  pd(21,45) = pd(21,45) + rrt(603) * density(37) 
  pd(37,37) = pd(37,37) - rrt(603) * density(45) 
  pd(37,45) = pd(37,45) - rrt(603) * density(37) 
  pd(40,37) = pd(40,37) + rrt(603) * density(45) 
  pd(40,45) = pd(40,45) + rrt(603) * density(37) 
  pd(45,37) = pd(45,37) - rrt(603) * density(45) 
  pd(45,45) = pd(45,45) - rrt(603) * density(37) 
  pd(17,17) = pd(17,17) - rrt(604) * density(36) 
  pd(17,36) = pd(17,36) - rrt(604) * density(17) 
  pd(36,17) = pd(36,17) - rrt(604) * density(36) 
  pd(36,36) = pd(36,36) - rrt(604) * density(17) 
  pd(40,17) = pd(40,17) + rrt(604) * density(36) 
  pd(40,36) = pd(40,36) + rrt(604) * density(17) 
  pd(18,18) = pd(18,18) - rrt(605) * density(36) 
  pd(18,36) = pd(18,36) - rrt(605) * density(18) 
  pd(36,18) = pd(36,18) - rrt(605) * density(36) 
  pd(36,36) = pd(36,36) - rrt(605) * density(18) 
  pd(41,18) = pd(41,18) + rrt(605) * density(36) 
  pd(41,36) = pd(41,36) + rrt(605) * density(18) 
  pd(21,33) = pd(21,33) + rrt(606) * density(36) 
  pd(21,36) = pd(21,36) + rrt(606) * density(33) 
  pd(33,33) = pd(33,33) - rrt(606) * density(36) 
  pd(33,36) = pd(33,36) - rrt(606) * density(33) 
  pd(36,33) = pd(36,33) - rrt(606) * density(36) 
  pd(36,36) = pd(36,36) - rrt(606) * density(33) 
  pd(32,34) = pd(32,34) + rrt(607) * density(36) 
  pd(32,36) = pd(32,36) + rrt(607) * density(34) 
  pd(34,34) = pd(34,34) - rrt(607) * density(36) 
  pd(34,36) = pd(34,36) - rrt(607) * density(34) 
  pd(36,34) = pd(36,34) - rrt(607) * density(36) 
  pd(36,36) = pd(36,36) - rrt(607) * density(34) 
  pd(36,36) = pd(36,36) - rrt(608) * density(45) 
  pd(36,45) = pd(36,45) - rrt(608) * density(36) 
  pd(42,36) = pd(42,36) + rrt(608) * density(45) 
  pd(42,45) = pd(42,45) + rrt(608) * density(36) 
  pd(45,36) = pd(45,36) - rrt(608) * density(45) 
  pd(45,45) = pd(45,45) - rrt(608) * density(36) 
  pd(17,17) = pd(17,17) - rrt(609) * density(37) 
  pd(17,37) = pd(17,37) - rrt(609) * density(17) 
  pd(37,17) = pd(37,17) - rrt(609) * density(37) 
  pd(37,37) = pd(37,37) - rrt(609) * density(17) 
  pd(42,17) = pd(42,17) + rrt(609) * density(37) 
  pd(42,37) = pd(42,37) + rrt(609) * density(17) 
  pd(32,33) = pd(32,33) + rrt(610) * density(37) 
  pd(32,37) = pd(32,37) + rrt(610) * density(33) 
  pd(33,33) = pd(33,33) - rrt(610) * density(37) 
  pd(33,37) = pd(33,37) - rrt(610) * density(33) 
  pd(37,33) = pd(37,33) - rrt(610) * density(37) 
  pd(37,37) = pd(37,37) - rrt(610) * density(33) 
  pd(37,37) = pd(37,37) - rrt(611) * density(45) 
  pd(37,45) = pd(37,45) - rrt(611) * density(37) 
  pd(43,37) = pd(43,37) + rrt(611) * density(45) 
  pd(43,45) = pd(43,45) + rrt(611) * density(37) 
  pd(45,37) = pd(45,37) - rrt(611) * density(45) 
  pd(45,45) = pd(45,45) - rrt(611) * density(37) 
  pd(14,17) = pd(14,17) + rrt(612) * density(38) 
  pd(14,38) = pd(14,38) + rrt(612) * density(17) 
  pd(17,17) = pd(17,17) - rrt(612) * density(38) 
  pd(17,38) = pd(17,38) - rrt(612) * density(17) 
  pd(32,17) = pd(32,17) + rrt(612) * density(38) 
  pd(32,38) = pd(32,38) + rrt(612) * density(17) 
  pd(38,17) = pd(38,17) - rrt(612) * density(38) 
  pd(38,38) = pd(38,38) - rrt(612) * density(17) 
  pd(01,18) = pd(01,18) + rrt(613) * density(38) 
  pd(01,38) = pd(01,38) + rrt(613) * density(18) 
  pd(18,18) = pd(18,18) - rrt(613) * density(38) 
  pd(18,38) = pd(18,38) - rrt(613) * density(18) 
  pd(32,18) = pd(32,18) + rrt(613) * density(38) 
  pd(32,38) = pd(32,38) + rrt(613) * density(18) 
  pd(38,18) = pd(38,18) - rrt(613) * density(38) 
  pd(38,38) = pd(38,38) - rrt(613) * density(18) 
  pd(29,33) = pd(29,33) + rrt(614) * density(38) 
  pd(29,38) = pd(29,38) + rrt(614) * density(33) 
  pd(32,33) = pd(32,33) + rrt(614) * density(38) 
  pd(32,38) = pd(32,38) + rrt(614) * density(33) 
  pd(33,33) = pd(33,33) - rrt(614) * density(38) 
  pd(33,38) = pd(33,38) - rrt(614) * density(33) 
  pd(38,33) = pd(38,33) - rrt(614) * density(38) 
  pd(38,38) = pd(38,38) - rrt(614) * density(33) 
  pd(21,34) = pd(21,34) + rrt(615) * density(38) 
  pd(21,38) = pd(21,38) + rrt(615) * density(34) 
  pd(32,34) = pd(32,34) + rrt(615) * density(38) 
  pd(32,38) = pd(32,38) + rrt(615) * density(34) 
  pd(34,34) = pd(34,34) - rrt(615) * density(38) 
  pd(34,38) = pd(34,38) - rrt(615) * density(34) 
  pd(38,34) = pd(38,34) - rrt(615) * density(38) 
  pd(38,38) = pd(38,38) - rrt(615) * density(34) 
  pd(32,38) = pd(32,38) + rrt(616) * density(45) 
  pd(32,45) = pd(32,45) + rrt(616) * density(38) 
  pd(38,38) = pd(38,38) - rrt(616) * density(45) 
  pd(38,45) = pd(38,45) - rrt(616) * density(38) 
  pd(40,38) = pd(40,38) + rrt(616) * density(45) 
  pd(40,45) = pd(40,45) + rrt(616) * density(38) 
  pd(45,38) = pd(45,38) - rrt(616) * density(45) 
  pd(45,45) = pd(45,45) - rrt(616) * density(38) 
  pd(32,38) = pd(32,38) + rrt(617) * density(46) 
  pd(32,46) = pd(32,46) + rrt(617) * density(38) 
  pd(38,38) = pd(38,38) - rrt(617) * density(46) 
  pd(38,46) = pd(38,46) - rrt(617) * density(38) 
  pd(41,38) = pd(41,38) + rrt(617) * density(46) 
  pd(41,46) = pd(41,46) + rrt(617) * density(38) 
  pd(46,38) = pd(46,38) - rrt(617) * density(46) 
  pd(46,46) = pd(46,46) - rrt(617) * density(38) 
  pd(32,38) = pd(32,38) + rrt(618) * density(47) 
  pd(32,47) = pd(32,47) + rrt(618) * density(38) 
  pd(38,38) = pd(38,38) - rrt(618) * density(47) 
  pd(38,47) = pd(38,47) - rrt(618) * density(38) 
  pd(42,38) = pd(42,38) + rrt(618) * density(47) 
  pd(42,47) = pd(42,47) + rrt(618) * density(38) 
  pd(47,38) = pd(47,38) - rrt(618) * density(47) 
  pd(47,47) = pd(47,47) - rrt(618) * density(38) 
  pd(14,17) = pd(14,17) + rrt(619) * density(48) 
  pd(14,48) = pd(14,48) + rrt(619) * density(17) 
  pd(17,17) = pd(17,17) - rrt(619) * density(48) 
  pd(17,48) = pd(17,48) - rrt(619) * density(17) 
  pd(40,17) = pd(40,17) + rrt(619) * density(48) 
  pd(40,48) = pd(40,48) + rrt(619) * density(17) 
  pd(48,17) = pd(48,17) - rrt(619) * density(48) 
  pd(48,48) = pd(48,48) - rrt(619) * density(17) 
  pd(01,18) = pd(01,18) + rrt(620) * density(48) 
  pd(01,48) = pd(01,48) + rrt(620) * density(18) 
  pd(18,18) = pd(18,18) - rrt(620) * density(48) 
  pd(18,48) = pd(18,48) - rrt(620) * density(18) 
  pd(40,18) = pd(40,18) + rrt(620) * density(48) 
  pd(40,48) = pd(40,48) + rrt(620) * density(18) 
  pd(48,18) = pd(48,18) - rrt(620) * density(48) 
  pd(48,48) = pd(48,48) - rrt(620) * density(18) 
  pd(29,33) = pd(29,33) + rrt(621) * density(48) 
  pd(29,48) = pd(29,48) + rrt(621) * density(33) 
  pd(33,33) = pd(33,33) - rrt(621) * density(48) 
  pd(33,48) = pd(33,48) - rrt(621) * density(33) 
  pd(40,33) = pd(40,33) + rrt(621) * density(48) 
  pd(40,48) = pd(40,48) + rrt(621) * density(33) 
  pd(48,33) = pd(48,33) - rrt(621) * density(48) 
  pd(48,48) = pd(48,48) - rrt(621) * density(33) 
  pd(21,34) = pd(21,34) + rrt(622) * density(48) 
  pd(21,48) = pd(21,48) + rrt(622) * density(34) 
  pd(34,34) = pd(34,34) - rrt(622) * density(48) 
  pd(34,48) = pd(34,48) - rrt(622) * density(34) 
  pd(40,34) = pd(40,34) + rrt(622) * density(48) 
  pd(40,48) = pd(40,48) + rrt(622) * density(34) 
  pd(48,34) = pd(48,34) - rrt(622) * density(48) 
  pd(48,48) = pd(48,48) - rrt(622) * density(34) 
  pd(40,45) = pd(40,45) + rrt(623) * density(48) * 2.0d0
  pd(40,48) = pd(40,48) + rrt(623) * density(45) * 2.0d0
  pd(45,45) = pd(45,45) - rrt(623) * density(48) 
  pd(45,48) = pd(45,48) - rrt(623) * density(45) 
  pd(48,45) = pd(48,45) - rrt(623) * density(48) 
  pd(48,48) = pd(48,48) - rrt(623) * density(45) 
  pd(40,46) = pd(40,46) + rrt(624) * density(48) 
  pd(40,48) = pd(40,48) + rrt(624) * density(46) 
  pd(41,46) = pd(41,46) + rrt(624) * density(48) 
  pd(41,48) = pd(41,48) + rrt(624) * density(46) 
  pd(46,46) = pd(46,46) - rrt(624) * density(48) 
  pd(46,48) = pd(46,48) - rrt(624) * density(46) 
  pd(48,46) = pd(48,46) - rrt(624) * density(48) 
  pd(48,48) = pd(48,48) - rrt(624) * density(46) 
  pd(40,47) = pd(40,47) + rrt(625) * density(48) 
  pd(40,48) = pd(40,48) + rrt(625) * density(47) 
  pd(42,47) = pd(42,47) + rrt(625) * density(48) 
  pd(42,48) = pd(42,48) + rrt(625) * density(47) 
  pd(47,47) = pd(47,47) - rrt(625) * density(48) 
  pd(47,48) = pd(47,48) - rrt(625) * density(47) 
  pd(48,47) = pd(48,47) - rrt(625) * density(48) 
  pd(48,48) = pd(48,48) - rrt(625) * density(47) 
  pd(14,17) = pd(14,17) + rrt(626) * density(49) 
  pd(14,49) = pd(14,49) + rrt(626) * density(17) 
  pd(17,17) = pd(17,17) - rrt(626) * density(49) 
  pd(17,49) = pd(17,49) - rrt(626) * density(17) 
  pd(41,17) = pd(41,17) + rrt(626) * density(49) 
  pd(41,49) = pd(41,49) + rrt(626) * density(17) 
  pd(49,17) = pd(49,17) - rrt(626) * density(49) 
  pd(49,49) = pd(49,49) - rrt(626) * density(17) 
  pd(01,18) = pd(01,18) + rrt(627) * density(49) 
  pd(01,49) = pd(01,49) + rrt(627) * density(18) 
  pd(18,18) = pd(18,18) - rrt(627) * density(49) 
  pd(18,49) = pd(18,49) - rrt(627) * density(18) 
  pd(41,18) = pd(41,18) + rrt(627) * density(49) 
  pd(41,49) = pd(41,49) + rrt(627) * density(18) 
  pd(49,18) = pd(49,18) - rrt(627) * density(49) 
  pd(49,49) = pd(49,49) - rrt(627) * density(18) 
  pd(29,33) = pd(29,33) + rrt(628) * density(49) 
  pd(29,49) = pd(29,49) + rrt(628) * density(33) 
  pd(33,33) = pd(33,33) - rrt(628) * density(49) 
  pd(33,49) = pd(33,49) - rrt(628) * density(33) 
  pd(41,33) = pd(41,33) + rrt(628) * density(49) 
  pd(41,49) = pd(41,49) + rrt(628) * density(33) 
  pd(49,33) = pd(49,33) - rrt(628) * density(49) 
  pd(49,49) = pd(49,49) - rrt(628) * density(33) 
  pd(21,34) = pd(21,34) + rrt(629) * density(49) 
  pd(21,49) = pd(21,49) + rrt(629) * density(34) 
  pd(34,34) = pd(34,34) - rrt(629) * density(49) 
  pd(34,49) = pd(34,49) - rrt(629) * density(34) 
  pd(41,34) = pd(41,34) + rrt(629) * density(49) 
  pd(41,49) = pd(41,49) + rrt(629) * density(34) 
  pd(49,34) = pd(49,34) - rrt(629) * density(49) 
  pd(49,49) = pd(49,49) - rrt(629) * density(34) 
  pd(40,45) = pd(40,45) + rrt(630) * density(49) 
  pd(40,49) = pd(40,49) + rrt(630) * density(45) 
  pd(41,45) = pd(41,45) + rrt(630) * density(49) 
  pd(41,49) = pd(41,49) + rrt(630) * density(45) 
  pd(45,45) = pd(45,45) - rrt(630) * density(49) 
  pd(45,49) = pd(45,49) - rrt(630) * density(45) 
  pd(49,45) = pd(49,45) - rrt(630) * density(49) 
  pd(49,49) = pd(49,49) - rrt(630) * density(45) 
  pd(41,46) = pd(41,46) + rrt(631) * density(49) * 2.0d0
  pd(41,49) = pd(41,49) + rrt(631) * density(46) * 2.0d0
  pd(46,46) = pd(46,46) - rrt(631) * density(49) 
  pd(46,49) = pd(46,49) - rrt(631) * density(46) 
  pd(49,46) = pd(49,46) - rrt(631) * density(49) 
  pd(49,49) = pd(49,49) - rrt(631) * density(46) 
  pd(41,47) = pd(41,47) + rrt(632) * density(49) 
  pd(41,49) = pd(41,49) + rrt(632) * density(47) 
  pd(42,47) = pd(42,47) + rrt(632) * density(49) 
  pd(42,49) = pd(42,49) + rrt(632) * density(47) 
  pd(47,47) = pd(47,47) - rrt(632) * density(49) 
  pd(47,49) = pd(47,49) - rrt(632) * density(47) 
  pd(49,47) = pd(49,47) - rrt(632) * density(49) 
  pd(49,49) = pd(49,49) - rrt(632) * density(47) 
  pd(14,17) = pd(14,17) + rrt(633) * density(50) 
  pd(14,50) = pd(14,50) + rrt(633) * density(17) 
  pd(17,17) = pd(17,17) - rrt(633) * density(50) 
  pd(17,50) = pd(17,50) - rrt(633) * density(17) 
  pd(42,17) = pd(42,17) + rrt(633) * density(50) 
  pd(42,50) = pd(42,50) + rrt(633) * density(17) 
  pd(50,17) = pd(50,17) - rrt(633) * density(50) 
  pd(50,50) = pd(50,50) - rrt(633) * density(17) 
  pd(01,18) = pd(01,18) + rrt(634) * density(50) 
  pd(01,50) = pd(01,50) + rrt(634) * density(18) 
  pd(18,18) = pd(18,18) - rrt(634) * density(50) 
  pd(18,50) = pd(18,50) - rrt(634) * density(18) 
  pd(42,18) = pd(42,18) + rrt(634) * density(50) 
  pd(42,50) = pd(42,50) + rrt(634) * density(18) 
  pd(50,18) = pd(50,18) - rrt(634) * density(50) 
  pd(50,50) = pd(50,50) - rrt(634) * density(18) 
  pd(29,33) = pd(29,33) + rrt(635) * density(50) 
  pd(29,50) = pd(29,50) + rrt(635) * density(33) 
  pd(33,33) = pd(33,33) - rrt(635) * density(50) 
  pd(33,50) = pd(33,50) - rrt(635) * density(33) 
  pd(42,33) = pd(42,33) + rrt(635) * density(50) 
  pd(42,50) = pd(42,50) + rrt(635) * density(33) 
  pd(50,33) = pd(50,33) - rrt(635) * density(50) 
  pd(50,50) = pd(50,50) - rrt(635) * density(33) 
  pd(21,34) = pd(21,34) + rrt(636) * density(50) 
  pd(21,50) = pd(21,50) + rrt(636) * density(34) 
  pd(34,34) = pd(34,34) - rrt(636) * density(50) 
  pd(34,50) = pd(34,50) - rrt(636) * density(34) 
  pd(42,34) = pd(42,34) + rrt(636) * density(50) 
  pd(42,50) = pd(42,50) + rrt(636) * density(34) 
  pd(50,34) = pd(50,34) - rrt(636) * density(50) 
  pd(50,50) = pd(50,50) - rrt(636) * density(34) 
  pd(40,45) = pd(40,45) + rrt(637) * density(50) 
  pd(40,50) = pd(40,50) + rrt(637) * density(45) 
  pd(42,45) = pd(42,45) + rrt(637) * density(50) 
  pd(42,50) = pd(42,50) + rrt(637) * density(45) 
  pd(45,45) = pd(45,45) - rrt(637) * density(50) 
  pd(45,50) = pd(45,50) - rrt(637) * density(45) 
  pd(50,45) = pd(50,45) - rrt(637) * density(50) 
  pd(50,50) = pd(50,50) - rrt(637) * density(45) 
  pd(41,46) = pd(41,46) + rrt(638) * density(50) 
  pd(41,50) = pd(41,50) + rrt(638) * density(46) 
  pd(42,46) = pd(42,46) + rrt(638) * density(50) 
  pd(42,50) = pd(42,50) + rrt(638) * density(46) 
  pd(46,46) = pd(46,46) - rrt(638) * density(50) 
  pd(46,50) = pd(46,50) - rrt(638) * density(46) 
  pd(50,46) = pd(50,46) - rrt(638) * density(50) 
  pd(50,50) = pd(50,50) - rrt(638) * density(46) 
  pd(42,47) = pd(42,47) + rrt(639) * density(50) * 2.0d0
  pd(42,50) = pd(42,50) + rrt(639) * density(47) * 2.0d0
  pd(47,47) = pd(47,47) - rrt(639) * density(50) 
  pd(47,50) = pd(47,50) - rrt(639) * density(47) 
  pd(50,47) = pd(50,47) - rrt(639) * density(50) 
  pd(50,50) = pd(50,50) - rrt(639) * density(47) 
  pd(14,17) = pd(14,17) + rrt(640) * density(51) 
  pd(14,51) = pd(14,51) + rrt(640) * density(17) 
  pd(17,17) = pd(17,17) - rrt(640) * density(51) 
  pd(17,51) = pd(17,51) - rrt(640) * density(17) 
  pd(43,17) = pd(43,17) + rrt(640) * density(51) 
  pd(43,51) = pd(43,51) + rrt(640) * density(17) 
  pd(51,17) = pd(51,17) - rrt(640) * density(51) 
  pd(51,51) = pd(51,51) - rrt(640) * density(17) 
  pd(01,18) = pd(01,18) + rrt(641) * density(51) 
  pd(01,51) = pd(01,51) + rrt(641) * density(18) 
  pd(18,18) = pd(18,18) - rrt(641) * density(51) 
  pd(18,51) = pd(18,51) - rrt(641) * density(18) 
  pd(43,18) = pd(43,18) + rrt(641) * density(51) 
  pd(43,51) = pd(43,51) + rrt(641) * density(18) 
  pd(51,18) = pd(51,18) - rrt(641) * density(51) 
  pd(51,51) = pd(51,51) - rrt(641) * density(18) 
  pd(29,33) = pd(29,33) + rrt(642) * density(51) 
  pd(29,51) = pd(29,51) + rrt(642) * density(33) 
  pd(33,33) = pd(33,33) - rrt(642) * density(51) 
  pd(33,51) = pd(33,51) - rrt(642) * density(33) 
  pd(43,33) = pd(43,33) + rrt(642) * density(51) 
  pd(43,51) = pd(43,51) + rrt(642) * density(33) 
  pd(51,33) = pd(51,33) - rrt(642) * density(51) 
  pd(51,51) = pd(51,51) - rrt(642) * density(33) 
  pd(21,34) = pd(21,34) + rrt(643) * density(51) 
  pd(21,51) = pd(21,51) + rrt(643) * density(34) 
  pd(34,34) = pd(34,34) - rrt(643) * density(51) 
  pd(34,51) = pd(34,51) - rrt(643) * density(34) 
  pd(43,34) = pd(43,34) + rrt(643) * density(51) 
  pd(43,51) = pd(43,51) + rrt(643) * density(34) 
  pd(51,34) = pd(51,34) - rrt(643) * density(51) 
  pd(51,51) = pd(51,51) - rrt(643) * density(34) 
  pd(40,45) = pd(40,45) + rrt(644) * density(51) 
  pd(40,51) = pd(40,51) + rrt(644) * density(45) 
  pd(43,45) = pd(43,45) + rrt(644) * density(51) 
  pd(43,51) = pd(43,51) + rrt(644) * density(45) 
  pd(45,45) = pd(45,45) - rrt(644) * density(51) 
  pd(45,51) = pd(45,51) - rrt(644) * density(45) 
  pd(51,45) = pd(51,45) - rrt(644) * density(51) 
  pd(51,51) = pd(51,51) - rrt(644) * density(45) 
  pd(41,46) = pd(41,46) + rrt(645) * density(51) 
  pd(41,51) = pd(41,51) + rrt(645) * density(46) 
  pd(43,46) = pd(43,46) + rrt(645) * density(51) 
  pd(43,51) = pd(43,51) + rrt(645) * density(46) 
  pd(46,46) = pd(46,46) - rrt(645) * density(51) 
  pd(46,51) = pd(46,51) - rrt(645) * density(46) 
  pd(51,46) = pd(51,46) - rrt(645) * density(51) 
  pd(51,51) = pd(51,51) - rrt(645) * density(46) 
  pd(42,47) = pd(42,47) + rrt(646) * density(51) 
  pd(42,51) = pd(42,51) + rrt(646) * density(47) 
  pd(43,47) = pd(43,47) + rrt(646) * density(51) 
  pd(43,51) = pd(43,51) + rrt(646) * density(47) 
  pd(47,47) = pd(47,47) - rrt(646) * density(51) 
  pd(47,51) = pd(47,51) - rrt(646) * density(47) 
  pd(51,47) = pd(51,47) - rrt(646) * density(51) 
  pd(51,51) = pd(51,51) - rrt(646) * density(47) 
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
  rrt(020) = KVT10_N2N2*1.0D0
  rrt(021) = KVT10_N2N2*2.0D0
  rrt(022) = KVT10_N2N2*3.0D0
  rrt(023) = KVT10_N2N2*4.0D0
  rrt(024) = KVT10_N2N2*5.0D0
  rrt(025) = KVT10_N2N2*6.0D0
  rrt(026) = KVT10_N2N2*7.0D0
  rrt(027) = KVT10_N2N2*8.0D0
  rrt(028) = KVT01_N2N2*1.0D0
  rrt(029) = KVT01_N2N2*2.0D0
  rrt(030) = KVT01_N2N2*3.0D0
  rrt(031) = KVT01_N2N2*4.0D0
  rrt(032) = KVT01_N2N2*5.0D0
  rrt(033) = KVT01_N2N2*6.0D0
  rrt(034) = KVT01_N2N2*7.0D0
  rrt(035) = KVT01_N2N2*8.0D0
  rrt(036) = KVT10_N2N*1.0D0
  rrt(037) = KVT10_N2N*2.0D0
  rrt(038) = KVT10_N2N*3.0D0
  rrt(039) = KVT10_N2N*4.0D0
  rrt(040) = KVT10_N2N*5.0D0
  rrt(041) = KVT10_N2N*6.0D0
  rrt(042) = KVT10_N2N*7.0D0
  rrt(043) = KVT10_N2N*8.0D0
  rrt(044) = KVT01_N2N*1.0D0
  rrt(045) = KVT01_N2N*2.0D0
  rrt(046) = KVT01_N2N*3.0D0
  rrt(047) = KVT01_N2N*4.0D0
  rrt(048) = KVT01_N2N*5.0D0
  rrt(049) = KVT01_N2N*6.0D0
  rrt(050) = KVT01_N2N*7.0D0
  rrt(051) = KVT01_N2N*8.0D0
  rrt(052) = KVT10_N2O*1.0D0
  rrt(053) = KVT10_N2O*2.0D0
  rrt(054) = KVT10_N2O*3.0D0
  rrt(055) = KVT10_N2O*4.0D0
  rrt(056) = KVT10_N2O*5.0D0
  rrt(057) = KVT10_N2O*6.0D0
  rrt(058) = KVT10_N2O*7.0D0
  rrt(059) = KVT10_N2O*8.0D0
  rrt(060) = KVT01_N2O*1.0D0
  rrt(061) = KVT01_N2O*2.0D0
  rrt(062) = KVT01_N2O*3.0D0
  rrt(063) = KVT01_N2O*4.0D0
  rrt(064) = KVT01_N2O*5.0D0
  rrt(065) = KVT01_N2O*6.0D0
  rrt(066) = KVT01_N2O*7.0D0
  rrt(067) = KVT01_N2O*8.0D0
  rrt(068) = KVT10_O2O2*1.0D0
  rrt(069) = KVT10_O2O2*2.0D0
  rrt(070) = KVT10_O2O2*3.0D0
  rrt(071) = KVT10_O2O2*4.0D0
  rrt(072) = KVT01_O2O2*1.0D0
  rrt(073) = KVT01_O2O2*2.0D0
  rrt(074) = KVT01_O2O2*3.0D0
  rrt(075) = KVT01_O2O2*4.0D0
  rrt(076) = KVT10_O2O*1.0D0
  rrt(077) = KVT10_O2O*2.0D0
  rrt(078) = KVT10_O2O*3.0D0
  rrt(079) = KVT10_O2O*4.0D0
  rrt(080) = KVT01_O2O*1.0D0
  rrt(081) = KVT01_O2O*2.0D0
  rrt(082) = KVT01_O2O*3.0D0
  rrt(083) = KVT01_O2O*4.0D0
  rrt(084) = bolsig_rates(bolsig_pointer(20))
  rrt(085) = bolsig_rates(bolsig_pointer(21))
  rrt(086) = bolsig_rates(bolsig_pointer(22))
  rrt(087) = bolsig_rates(bolsig_pointer(23))
  rrt(088) = bolsig_rates(bolsig_pointer(24))
  rrt(089) = bolsig_rates(bolsig_pointer(25))
  rrt(090) = bolsig_rates(bolsig_pointer(26))
  rrt(091) = bolsig_rates(bolsig_pointer(27))
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
  rrt(113) = 1.8D-7*(300.0D0/TE)**0.39*0.50D0
  rrt(114) = 1.8D-7*(300.0D0/TE)**0.39*0.45D0
  rrt(115) = 1.8D-7*(300.0D0/TE)**0.39*0.05D0
  rrt(116) = 2.7D-7*(300.0D0/TE)**0.7*0.55D0
  rrt(117) = 2.7D-7*(300.0D0/TE)**0.7*0.40D0
  rrt(118) = 2.7D-7*(300.0D0/TE)**0.7*0.05D0
  rrt(119) = 4.2D-7*(300.0D0/TE)**0.85*0.20D0
  rrt(120) = 4.2D-7*(300.0D0/TE)**0.85*0.80D0
  rrt(121) = 2.0D-7*(300.0D0/TE)**0.5
  rrt(122) = 2.3D-6*(300.0D0/TE)**0.53
  rrt(123) = rrt(121)
  rrt(124) = rrt(121)
  rrt(125) = 1.4D-6*(300.0D0/TE)**0.5
  rrt(126) = 1.3D-6*(300.0D0/TE)**0.5
  rrt(127) = 7.0D-20*(300.0D0/TE)**4.5
  rrt(128) = rrt(127)
  rrt(129) = 6.0D-27*(300.0D0/TE)**1.5*ANY_NEUTRAL
  rrt(130) = rrt(129)
  rrt(131) = bolsig_rates(bolsig_pointer(49))
  rrt(132) = bolsig_rates(bolsig_pointer(50))
  rrt(133) = bolsig_rates(bolsig_pointer(51))
  rrt(134) = bolsig_rates(bolsig_pointer(52))
  rrt(135) = bolsig_rates(bolsig_pointer(53))
  rrt(136) = bolsig_rates(bolsig_pointer(54))
  rrt(137) = 1.0D-11
  rrt(138) = 1.0D-31
  rrt(139) = 1.0D-31
  rrt(140) = 1.0D-31*ANY_NEUTRAL
  rrt(141) = 8.0D-31*ANY_NEUTRAL
  rrt(142) = 6.0D-33*ANY_NEUTRAL
  rrt(143) = 1.1D-31*(300.0D0/TE)**2*EXP(-70.0D0/TGAS)*EXP(1500.0D0*(TE-TGAS)/(TE*TGAS))
  rrt(144) = 1.4D-10
  rrt(145) = 2.6D-10
  rrt(146) = 2.6D-10
  rrt(147) = 5.0D-13
  rrt(148) = 5.0D-15
  rrt(149) = 3.0D-10
  rrt(150) = 6.9D-10
  rrt(151) = 2.2D-9
  rrt(152) = 1.9D-9
  rrt(153) = 3.0D-10
  rrt(154) = 1.5D-10
  rrt(155) = 5.0D-10
  rrt(156) = 2.7D-10*(TEFFN2/300.0D0)**0.5*EXP(-5590.0D0/TEFFN2)
  rrt(157) = 2.0D-10
  rrt(158) = 3.6D-10
  rrt(159) = 1.9D-12*(TEFFN2/300.0D0)**0.5*EXP(-4990.0D0/TEFFN2)
  rrt(160) = 2.1D-9
  rrt(161) = 2.5D-9
  rrt(162) = 3.0D-10
  rrt(163) = 5.0D-10
  rrt(164) = 5.0D-10
  rrt(165) = 5.0D-10
  rrt(166) = 5.0D-10
  rrt(167) = 5.0D-10
  rrt(168) = 1.5D-10
  rrt(169) = 1.5D-10
  rrt(170) = 1.5D-10
  rrt(171) = 1.5D-10
  rrt(172) = 2.1D-9
  rrt(173) = 2.1D-9
  rrt(174) = 2.1D-9
  rrt(175) = 2.1D-9
  rrt(176) = 2.1D-9
  rrt(177) = 2.5D-9
  rrt(178) = 2.5D-9
  rrt(179) = 2.5D-9
  rrt(180) = 2.5D-9
  rrt(181) = 2.5D-9
  rrt(182) = 0.50D0
  rrt(183) = 1.34D5
  rrt(184) = 1.0D2
  rrt(185) = 2.45D7
  rrt(186) = 2.6D-4
  rrt(187) = 1.5D-3
  rrt(188) = 8.5D-2
  rrt(189) = 11.0D0
  rrt(190) = 7.0D-12
  rrt(191) = 2.1D-11
  rrt(192) = 2.0D-12
  rrt(193) = 4.0D-11*(300.0D0/TGAS)**0.667
  rrt(194) = 2.1D-12*(TGAS/300.0D0)**0.55
  rrt(195) = 2.0D-13*(TGAS/300.0D0)**0.55
  rrt(196) = rrt(195)
  rrt(197) = 2.0D-14*(TGAS/300.0D0)**0.55
  rrt(198) = 3.0D-16
  rrt(199) = 6.9D-11
  rrt(200) = 1.0D-11
  rrt(201) = 1.0D-12
  rrt(202) = 3.0D-10
  rrt(203) = 1.5D-10
  rrt(204) = 3.0D-11
  rrt(205) = 2.0D-12
  rrt(206) = 3.0D-10
  rrt(207) = 2.4D-10
  rrt(208) = 1.0D-11
  rrt(209) = 3.0D-10
  rrt(210) = 1.9D-13
  rrt(211) = 2.8D-11
  rrt(212) = 3.6D-10
  rrt(213) = 4.0D-12
  rrt(214) = 1.0D-11
  rrt(215) = 1.7D-33
  rrt(216) = 1.7D-33
  rrt(217) = 1.7D-33
  rrt(218) = 1.0D-32
  rrt(219) = 1.0D-32
  rrt(220) = 2.4D-33
  rrt(221) = 2.4D-33
  rrt(222) = 2.4D-33
  rrt(223) = 1.4D-32
  rrt(224) = 1.4D-32
  rrt(225) = 4.0D-13
  rrt(226) = 5.2D-12
  rrt(227) = 1.8D-10
  rrt(228) = 3.5D-12
  rrt(229) = 1.0D-13*EXP(-510.0D0/TGAS)
  rrt(230) = 1.8D-12
  rrt(231) = 1.0D-12
  rrt(232) = 6.0D-13
  rrt(233) = 6.0D-14
  rrt(234) = 1.0D-13
  rrt(235) = 2.6D-12
  rrt(236) = 3.0D-11
  rrt(237) = 7.0D-16
  rrt(238) = 2.0D-14*EXP(-600.0D0/TGAS)
  rrt(239) = 3.8D-18*EXP(-205.0D0/TGAS)
  rrt(240) = 3.0D-21
  rrt(241) = 2.5D-11
  rrt(242) = 5.2D-11*EXP(-2840.0D0/TGAS)
  rrt(243) = 7.0D-28*TGAS**3.8*EXP(700.0D0/TGAS)
  rrt(244) = 1.0D-11*EXP(-2300.0D0/TGAS)
  rrt(245) = 8.1D-14
  rrt(246) = 3.4D-11*(300.0D0/TGAS)**0.1*EXP(-4200.0D0/TGAS)
  rrt(247) = 4.3D-22*TGAS**2.4*EXP(-281.0D0/TGAS)
  rrt(248) = 1.7D-15*(TGAS/300.0D0)
  rrt(249) = 6.0D-14
  rrt(250) = 2.2D-11
  rrt(251) = 9.0D-12
  rrt(252) = 3.0D-13
  rrt(253) = 9.0D-15
  rrt(254) = 8.0D-12
  rrt(255) = 6.4D-12*EXP(67.0D0/TGAS)
  rrt(256) = 1.0D-12
  rrt(257) = 2.6D-11*EXP(67.0D0/TGAS)
  rrt(258) = 2.3D-11
  rrt(259) = 1.2D-10
  rrt(260) = 1.2D-10
  rrt(261) = 1.7D-10
  rrt(262) = 7.2D-11
  rrt(263) = 4.4D-11
  rrt(264) = 5.0D-11*EXP(-300.0D0/TGAS)
  rrt(265) = 1.0D-12
  rrt(266) = 1.3D-12*EXP(-850.0D0/TGAS)
  rrt(267) = 3.0D-12*EXP(-850.0D0/TGAS)
  rrt(268) = 1.0D-17
  rrt(269) = 1.1D-10
  rrt(270) = 2.9D-11
  rrt(271) = 3.2D-11
  rrt(272) = 2.9D-10
  rrt(273) = 5.1D-10
  rrt(274) = 2.9D-10
  rrt(275) = 2.9D-10
  rrt(276) = 6.3D-12
  rrt(277) = 3.1D-12
  rrt(278) = 1.8D-11*(TGAS/300.0)**0.5
  rrt(279) = 3.2D-12*(TGAS/300.0)*EXP(-3150.0D0/TGAS)
  rrt(280) = 9.1D-13
  rrt(281) = 3.0D-12
  rrt(282) = 7.0D-13
  rrt(283) = 2.3D-12
  rrt(284) = 3.0D-10*EXP(-38370.0D0/TGAS)
  rrt(285) = 7.5D-12*(TGAS/300.0)*EXP(-19500.0D0/TGAS)
  rrt(286) = 4.2D-18
  rrt(287) = 8.3D-12*EXP(-14000.0D0/TGAS)
  rrt(288) = 1.5D-10*EXP(-14090.0D0/TGAS)
  rrt(289) = 9.1D-12*(TGAS/300.0D0)**0.18
  rrt(290) = 1.0D-11
  rrt(291) = 2.5D-10*EXP(-50390.0D0/TGAS)
  rrt(292) = 3.3D-16*(300.0D0/TGAS)**0.5*EXP(-39200.0D0/TGAS)
  rrt(293) = 2.2D-12*EXP(-32100.0D0/TGAS)
  rrt(294) = 5.1D-13*EXP(-33660.0D0/TGAS)
  rrt(295) = 2.8D-12*EXP(-23400.0D0/TGAS)
  rrt(296) = 2.5D-13*EXP(-765.0D0/TGAS)
  rrt(297) = 4.6D-10*EXP(-25170.0D0/TGAS)
  rrt(298) = 1.7D-11
  rrt(299) = 2.0D-11*EXP(-49800.0D0/TGAS)
  rrt(300) = 2.8D-12*EXP(-25400.0D0/TGAS)
  rrt(301) = 3.3D-12*EXP(-13500.0D0/TGAS)
  rrt(302) = 4.5D-10*EXP(-18500.0D0/TGAS)
  rrt(303) = 1.2D-13*EXP(-2450.0D0/TGAS)
  rrt(304) = 2.3D-13*EXP(-1600.0D0/TGAS)
  rrt(305) = 1.5D-12*EXP(-15020.0D0/TGAS)
  rrt(306) = 4.3D-12*EXP(-3850.0D0/TGAS)
  rrt(307) = 2.7D-11*EXP(-6.74D4/TGAS)
  rrt(308) = 1.6D-12*(TGAS/300.0D0)**0.5*(0.19D0+8.6D0*TGAS)*EXP(-32000.0D0/TGAS)
  rrt(309) = 5.4D-8*(1.0D0-EXP(-3354.0D0/TGAS))*EXP(-113200.0D0/TGAS)*1.0D0
  rrt(310) = rrt(309)
  rrt(311) = rrt(309)
  rrt(312) = 5.4D-8*(1.0D0-EXP(-3354.0D0/TGAS))*EXP(-113200.0D0/TGAS)*6.6D0
  rrt(313) = rrt(312)
  rrt(314) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*1.0D0
  rrt(315) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*5.9D0
  rrt(316) = 6.1D-9*(1.0D0-EXP(-2240.0D0/TGAS))*EXP(-59380.0D0/TGAS)*21.D0
  rrt(317) = rrt(314)
  rrt(318) = rrt(314)
  rrt(319) = 8.7D-9*EXP(-75994.0D0/TGAS)*1.0D0
  rrt(320) = rrt(319)
  rrt(321) = 8.7D-9*EXP(-75994.0D0/TGAS)*20.D0
  rrt(322) = rrt(321)
  rrt(323) = rrt(321)
  rrt(324) = 6.6D-10*EXP(-11600.0D0/TGAS)*1.0D0
  rrt(325) = 6.6D-10*EXP(-11600.0D0/TGAS)*0.38D0
  rrt(326) = 6.6D-10*EXP(-11600.0D0/TGAS)*6.3D0*EXP(170.0D0/TGAS)
  rrt(327) = rrt(326)
  rrt(328) = 1.2D-8*(300.0D0/TGAS)*EXP(-29000.0D0/TGAS)*1.0D0
  rrt(329) = rrt(328)
  rrt(330) = 1.2D-8*(300.0D0/TGAS)*EXP(-29000.0D0/TGAS)*2.0D0
  rrt(331) = 1.2D-8*(300.0D0/TGAS)*EXP(-29000.0D0/TGAS)*4.0D0
  rrt(332) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*1.0D0
  rrt(333) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*0.78D0
  rrt(334) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*7.8D0
  rrt(335) = 6.8D-6*(300.0D0/TGAS)**2*EXP(-36180.0D0/TGAS)*5.9D0
  rrt(336) = 3.1D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*1.0D0
  rrt(337) = rrt(336)
  rrt(338) = rrt(336)
  rrt(339) = 3.1D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*10.D0
  rrt(340) = rrt(339)
  rrt(341) = 6.2D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*1.0D0
  rrt(342) = rrt(341)
  rrt(343) = rrt(341)
  rrt(344) = 6.2D-5*(300.0D0/TGAS)**2*EXP(-25000.0D0/TGAS)*12.D0
  rrt(345) = rrt(344)
  rrt(346) = 2.1D-11*(300.0D0/TGAS)**4.4*EXP(-11080.0D0/TGAS)*ANY_NEUTRAL
  rrt(347) = MAX(8.3D-34*EXP(500.0D0/TGAS),1.91D-33)
  rrt(348) = 1.8D-33*EXP(435.0D0/TGAS)*1.0D0
  rrt(349) = rrt(348)
  rrt(350) = 1.8D-33*EXP(435.0D0/TGAS)*3.0D0
  rrt(351) = rrt(350)
  rrt(352) = MAX(2.8D-34*EXP(720.0D0/TGAS),1.0D-33*(300.0D0/TGAS)**0.41)
  rrt(353) = 4.0D-33*(300.0D0/TGAS)**0.41*1.0D0
  rrt(354) = 4.0D-33*(300.0D0/TGAS)**0.41*0.8D0
  rrt(355) = 4.0D-33*(300.0D0/TGAS)**0.41*3.6D0
  rrt(356) = 4.0D-33*(300.0D0/TGAS)**0.41*0.17D0
  rrt(357) = 1.0D-32*(300.0D0/TGAS)**0.5
  rrt(358) = rrt(357)
  rrt(359) = 1.8D-31*(300.0D0/TGAS)
  rrt(360) = rrt(359)
  rrt(361) = rrt(359)
  rrt(362) = MAX(5.8D-34*(300.0D0/TGAS)**2.8,5.4D-34*(300.0D0/TGAS)**1.9)
  rrt(363) = 7.6D-34*(300.0D0/TGAS)**1.9
  rrt(364) = rrt(363)
  rrt(365) = MIN(3.9D-33*(300.0D0/TGAS)**1.9,1.1D-34*EXP(1060.0D0/TGAS))
  rrt(366) = rrt(365)
  rrt(367) = 3.9D-35*EXP(-10400.0D0/TGAS)*ANY_NEUTRAL
  rrt(368) = 1.2D-31*(300.0D0/TGAS)**1.8*1.0D0
  rrt(369) = 1.2D-31*(300.0D0/TGAS)**1.8*0.78D0
  rrt(370) = rrt(369)
  rrt(371) = 8.9D-32*(300.0D0/TGAS)**2*1.0D0
  rrt(372) = rrt(371)
  rrt(373) = 8.9D-32*(300.0D0/TGAS)**2*13.D0
  rrt(374) = rrt(373)
  rrt(375) = 8.9D-32*(300.0D0/TGAS)**2*2.4D0
  rrt(376) = 3.7D-30*(300.0D0/TGAS)**4.1*ANY_NEUTRAL
  rrt(377) = 1.0D-12
  rrt(378) = 2.8D-10
  rrt(379) = 2.5D-10
  rrt(380) = 2.8D-11
  rrt(381) = 5.0D-10
  rrt(382) = 8.0D-10
  rrt(383) = 3.0D-12
  rrt(384) = 1.0D-12
  rrt(385) = 5.5D-10
  rrt(386) = (1.5D0-2.0D-3*TEFFN+9.6D-7*TEFFN**2)*1.0D-12
  rrt(387) = 2.0D-11*(300.0D0/TEFFN)**0.5
  rrt(388) = 1.0D-10
  rrt(389) = 2.4D-11
  rrt(390) = 3.0D-12
  rrt(391) = 1.3D-10
  rrt(392) = 2.3D-10
  rrt(393) = 2.2D-10
  rrt(394) = 2.0D-11
  rrt(395) = 1.6D-9
  rrt(396) = 6.0D-11*(300.0D0/TEFFN2)**0.5
  rrt(397) = 1.3D-10*(300.0D0/TEFFN2)**0.5
  rrt(398) = 1.0D-10
  rrt(399) = 7.2D-13*(TEFFN2/300.0D0)
  rrt(400) = 3.3D-10
  rrt(401) = 5.0D-10
  rrt(402) = 4.0D-10
  rrt(403) = 1.0D-17
  rrt(404) = 1.2D-10
  rrt(405) = 6.3D-10
  rrt(406) = 1.0D-11
  rrt(407) = 6.6D-10
  rrt(408) = 2.3D-11
  rrt(409) = 4.4D-11
  rrt(410) = 6.6D-11
  rrt(411) = 7.0D-11
  rrt(412) = 7.0D-11
  rrt(413) = 2.9D-10
  rrt(414) = 2.9D-10
  rrt(415) = MIN(2.1D-16*EXP(TEFFN4/121.0D0),1.0D-10)
  rrt(416) = 2.5D-10
  rrt(417) = 2.5D-10
  rrt(418) = 1.0D-11
  rrt(419) = 4.0D-10
  rrt(420) = 4.6D-12*(TEFFN4/300.0D0)**2.5*EXP(-2650.0D0/TEFFN4)
  rrt(421) = 3.3D-6*(300.0D0/TEFFN4)**4*EXP(-5030.0D0/TEFFN4)
  rrt(422) = 1.0D-10
  rrt(423) = 1.0D-10
  rrt(424) = 3.0D-10
  rrt(425) = 1.0D-10
  rrt(426) = 1.1D-6*(300.0D0/TEFFN4)**5.3*EXP(-2360.0D0/TEFFN4)
  rrt(427) = 1.0D-9
  rrt(428) = 1.7D-29*(300.0D0/TEFFN)**2.1
  rrt(429) = 1.0D-29*ANY_NEUTRAL
  rrt(430) = rrt(429)
  rrt(431) = 6.0D-29*(300.0D0/TEFFN)**2*ANY_NEUTRAL
  rrt(432) = rrt(429)
  rrt(433) = rrt(429)
  rrt(434) = 5.2D-29*(300.0D0/TEFFN2)**2.2
  rrt(435) = 9.0D-30*EXP(400.0D0/TEFFN2)
  rrt(436) = 2.4D-30*(300.0D0/TEFFN2)**3.2
  rrt(437) = 9.0D-31*(300.0D0/TEFFN2)**2
  rrt(438) = 1.0D-10
  rrt(439) = 8.0D-10
  rrt(440) = 1.2D-9
  rrt(441) = 2.0D-10
  rrt(442) = 2.0D-12
  rrt(443) = 3.3D-10
  rrt(444) = 3.5D-10
  rrt(445) = 7.0D-10
  rrt(446) = 5.0D-10
  rrt(447) = 1.0D-11
  rrt(448) = 1.0D-11
  rrt(449) = 2.6D-12
  rrt(450) = 7.0D-11
  rrt(451) = 2.0D-11
  rrt(452) = 5.0D-10
  rrt(453) = 5.0D-10
  rrt(454) = 7.4D-10
  rrt(455) = 2.8D-14
  rrt(456) = 1.8D-11
  rrt(457) = 4.0D-12
  rrt(458) = 5.0D-10
  rrt(459) = 7.0D-10
  rrt(460) = 3.0D-15
  rrt(461) = 1.0D-10*EXP(-1044.0D0/TEFFN4)
  rrt(462) = rrt(461)
  rrt(463) = 4.0D-10
  rrt(464) = 3.0D-10
  rrt(465) = 1.0D-10
  rrt(466) = 1.0D-10
  rrt(467) = 2.5D-10
  rrt(468) = 1.1D-30*(300.0D0/TEFFN)*ANY_NEUTRAL
  rrt(469) = rrt(429)
  rrt(470) = 3.5D-31*(300.0D0/TEFFN2)*ANY_NEUTRAL
  rrt(471) = 2.0D-7*(300.0D0/TIONN)**0.5
  rrt(472) = rrt(471)
  rrt(473) = rrt(471)
  rrt(474) = rrt(471)
  rrt(475) = rrt(471)
  rrt(476) = rrt(471)
  rrt(477) = rrt(471)
  rrt(478) = rrt(471)
  rrt(479) = rrt(471)
  rrt(480) = rrt(471)
  rrt(481) = rrt(471)
  rrt(482) = rrt(471)
  rrt(483) = rrt(471)
  rrt(484) = rrt(471)
  rrt(485) = rrt(471)
  rrt(486) = rrt(471)
  rrt(487) = rrt(471)
  rrt(488) = rrt(471)
  rrt(489) = rrt(471)
  rrt(490) = rrt(471)
  rrt(491) = rrt(471)
  rrt(492) = rrt(471)
  rrt(493) = rrt(471)
  rrt(494) = rrt(471)
  rrt(495) = rrt(471)
  rrt(496) = rrt(471)
  rrt(497) = rrt(471)
  rrt(498) = rrt(471)
  rrt(499) = rrt(471)
  rrt(500) = rrt(471)
  rrt(501) = rrt(471)
  rrt(502) = rrt(471)
  rrt(503) = rrt(471)
  rrt(504) = rrt(471)
  rrt(505) = rrt(471)
  rrt(506) = rrt(471)
  rrt(507) = rrt(471)
  rrt(508) = rrt(471)
  rrt(509) = rrt(471)
  rrt(510) = rrt(471)
  rrt(511) = rrt(471)
  rrt(512) = rrt(471)
  rrt(513) = rrt(471)
  rrt(514) = rrt(471)
  rrt(515) = rrt(471)
  rrt(516) = rrt(471)
  rrt(517) = rrt(471)
  rrt(518) = rrt(471)
  rrt(519) = rrt(471)
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
  rrt(573) = 1.0D-7
  rrt(574) = 1.0D-7
  rrt(575) = 1.0D-7
  rrt(576) = 1.0D-7
  rrt(577) = 1.0D-7
  rrt(578) = 1.0D-7
  rrt(579) = 1.0D-7
  rrt(580) = 1.0D-7
  rrt(581) = 1.0D-7
  rrt(582) = 1.0D-7
  rrt(583) = 1.0D-7
  rrt(584) = 1.0D-7
  rrt(585) = 1.0D-7
  rrt(586) = 1.0D-7
  rrt(587) = 1.0D-7
  rrt(588) = 1.0D-7
  rrt(589) = 1.0D-7
  rrt(590) = 1.0D-7
  rrt(591) = 1.0D-7
  rrt(592) = 1.0D-7
  rrt(593) = 1.0D-7
  rrt(594) = 2.0D-25*(300.0D0/TIONN)**2.5*ANY_NEUTRAL
  rrt(595) = rrt(594)
  rrt(596) = rrt(594)
  rrt(597) = rrt(594)
  rrt(598) = rrt(594)
  rrt(599) = rrt(594)
  rrt(600) = rrt(594)
  rrt(601) = rrt(594)
  rrt(602) = rrt(594)
  rrt(603) = rrt(594)
  rrt(604) = rrt(594)
  rrt(605) = rrt(594)
  rrt(606) = rrt(594)
  rrt(607) = rrt(594)
  rrt(608) = rrt(594)
  rrt(609) = rrt(594)
  rrt(610) = rrt(594)
  rrt(611) = rrt(594)
  rrt(612) = 2.0D-25*(300.0D0/TIONN2)**2.5*ANY_NEUTRAL
  rrt(613) = rrt(612)
  rrt(614) = rrt(612)
  rrt(615) = rrt(612)
  rrt(616) = rrt(612)
  rrt(617) = rrt(612)
  rrt(618) = rrt(612)
  rrt(619) = rrt(612)
  rrt(620) = rrt(612)
  rrt(621) = rrt(612)
  rrt(622) = rrt(612)
  rrt(623) = rrt(612)
  rrt(624) = rrt(612)
  rrt(625) = rrt(612)
  rrt(626) = rrt(612)
  rrt(627) = rrt(612)
  rrt(628) = rrt(612)
  rrt(629) = rrt(612)
  rrt(630) = rrt(612)
  rrt(631) = rrt(612)
  rrt(632) = rrt(612)
  rrt(633) = rrt(612)
  rrt(634) = rrt(612)
  rrt(635) = rrt(612)
  rrt(636) = rrt(612)
  rrt(637) = rrt(612)
  rrt(638) = rrt(612)
  rrt(639) = rrt(612)
  rrt(640) = rrt(612)
  rrt(641) = rrt(612)
  rrt(642) = rrt(612)
  rrt(643) = rrt(612)
  rrt(644) = rrt(612)
  rrt(645) = rrt(612)
  rrt(646) = rrt(612)
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
