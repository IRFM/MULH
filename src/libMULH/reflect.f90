module reflect

contains
function reflect_electron(pcol,Ep,pvf1,tipo,seec,atype)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!			REFLECT ELECTRON
!
! Calculate the outgoing velocities of an electron that has been reflected off
! the wall, either elastically or inelastically.
!
! INPUTS:
!   - pcol: integer 1D array indicating which wall the electron collided with
!   - Ep: real, energy of the primary electron
!   - pvf1: 1D real array with velocity of primary electron at time of impact
!   - tipo: type of reflection, 1= elastic, 2= inelastic
!   - seec : integer, emission model
!   - atype: integer, type of analysis being performed
!
! OUTPUTS:
!   - reflect_electron: real 1D array with outgoing velocity
!
!

use precision_def
use constants
use ecuyer_taus_rng
implicit none
! Function arguments
integer, intent(in) :: pcol(4), tipo, seec
integer, intent(in), optional :: atype
real(long), intent(in) :: Ep, pvf1(3)
! Dummy variables
integer :: i, j
real(long) :: thetaf, vf, phi, rt1, rt2
! Function output
real(long) :: reflect_electron(1,3)

! Determine which is the normal component of the outgoing velocity
if (pcol(1)==1 .OR. pcol(3)==1) then
  i = 1          ! Normal component
  j = 2          ! Parallel component (not z)
elseif (pcol(2)==1 .OR. pcol(4)==1) then
  i = 2          ! Normal component
  j = 1          ! Parallel component (not z)
endif

! Obtain velocities for each case
if (tipo == 1) then      ! Electron is elastically reflected (specular reflection)
    
  reflect_electron(1,1) = ((-1)**i)*pvf1(1)    ! x-component of reflected velocity
  reflect_electron(1,2) = ((-1)**j)*pvf1(2)    ! y-component of reflected velocity
  reflect_electron(1,3) = pvf1(3)    ! z-component of reflected velocity


elseif (tipo == 2) then ! Electron is inelastically reflected
  if (seec == 3) then ! Cheng model with de Lara's apporach to calculating Re and Rr
    ! Generate a random angle to the normal from -pi/2 to pi/2
    if (atype == 8) then
      read(21,*) rt1
      ! Original code :
      !thetaf = asin(2*rt1 - 1)   ! Calculate emission angle of secondaries with respect to the normal
      ! Modification A Placais 2017.09.26:
      thetaf = asin(sqrt(rt1))
      ! Firstly, thetaf should be in [0,pi/2] rather than in [-pi/2,pi/2] (except if phis belongs to [0, pi])
      ! Secondly, asin(rt) will not give an uniform distribution. For furher details, see:
      ! John Greenwood, The correct and incorrect generation of a cosine distribution of scattered particles for Monte-Carlo
      ! modelling of vacuum systems, In Vacuum, Volume 67, Issue 2, 2002, Pages 217-222, ISSN 0042-207X,
      ! https://doi.org/10.1016/S0042-207X(02)00173-2.
      ! (http://www.sciencedirect.com/science/article/pii/S0042207X02001732)
      ! Keywords: Monte-Carlo; Cosine distribution; Gas scattering
      ! Generate a random velocity magnitude from 0 to pvf1 (random energy from 0 to Ep)
      read(21,*) rt2
      vf = sqrt(2.*rt2*Ep/me)
    elseif (atype == 9) then
      read(19,'(F17.15,1X,F17.15)') rt1, rt2
      ! Original code :
      !thetaf = asin(2*rt1 - 1)   ! Calculate emission angle of secondaries with respect to the normal
      ! Modification A Placais 2017.09.26:
      thetaf = asin(sqrt(rt1))
      ! Firstly, thetaf should be in [0,pi/2] rather than in [-pi/2,pi/2] (except if phis belongs to [0, pi])
      ! Secondly, asin(rt) will not give an uniform distribution. For furher details, see:
      ! John Greenwood, The correct and incorrect generation of a cosine distribution of scattered particles for Monte-Carlo
      ! modelling of vacuum systems, In Vacuum, Volume 67, Issue 2, 2002, Pages 217-222, ISSN 0042-207X,
      ! https://doi.org/10.1016/S0042-207X(02)00173-2.
      ! (http://www.sciencedirect.com/science/article/pii/S0042207X02001732)
      ! Keywords: Monte-Carlo; Cosine distribution; Gas scattering
      ! Generate a random velocity magnitude from 0 to pvf1 (random energy from 0 to Ep)
      ! Generate a random velocity magnitude from 0 to pvf1 (random energy from 0 to Ep)
      vf = sqrt(2.*rt2*Ep/me)
    else
      ! Original code :
      !thetaf = asin(2*taus88() - 1)   ! Calculate emission angle of secondaries with respect to the normal
      ! Modification A Placais 2017.09.26:
      thetaf = asin(sqrt(taus88()))
      ! Firstly, thetaf should be in [0,pi/2] rather than in [-pi/2,pi/2] (except if phis belongs to [0, pi])
      ! Secondly, asin(rt) will not give an uniform distribution. For furher details, see:
      ! John Greenwood, The correct and incorrect generation of a cosine distribution of scattered particles for Monte-Carlo
      ! modelling of vacuum systems, In Vacuum, Volume 67, Issue 2, 2002, Pages 217-222, ISSN 0042-207X,
      ! https://doi.org/10.1016/S0042-207X(02)00173-2.
      ! (http://www.sciencedirect.com/science/article/pii/S0042207X02001732)
      ! Keywords: Monte-Carlo; Cosine distribution; Gas scattering
      ! Generate a random velocity magnitude from 0 to pvf1 (random energy from 0 to Ep)
      ! Generate a random velocity magnitude from 0 to pvf1 (random energy from 0 to Ep)
      vf = sqrt(2.*taus88()*Ep/me)
    endif
    
    reflect_electron(1,i) = vf * cos(thetaf)      ! Normal component of reflected velocity
    
    if (pvf1(3) == 0.) then
        
      reflect_electron(1,j) = vf * sin(thetaf)	! If incident electron has no velocity along z,
      reflect_electron(1,3) = 0.		! the reflected will not have one neither.
        
    else
      phi = atan(abs(pvf1(j)/pvf1(3)))    ! tan(phi) = ratio of velocities parallel to the surface
    
      reflect_electron(1,j) = sin(phi) * vf * sin(thetaf) ! Reflected and incident electron have
      reflect_electron(1,3) = cos(phi) * vf * sin(thetaf) ! same azimuthal angle. Justified ??
       
    endif
    
    ! Make sure parallel components have the same signs as before collision
    if (pvf1(j) < 0) then
      reflect_electron(1,j) = - abs(reflect_electron(1,j))
    else
      reflect_electron(1,j) = abs(reflect_electron(1,j))
    endif
    
    if (pvf1(3) < 0) then
      reflect_electron(1,3) = - abs(reflect_electron(1,3))
    else
      reflect_electron(1,3) = abs(reflect_electron(1,3))
    endif
  
    ! Make sure normal velocity points into the inside of the waveguide
    if (pcol(1) == 1 .OR. pcol(2) == 1) then
      reflect_electron(1,i) = -abs(reflect_electron(1,i))
    elseif (pcol(3) == 1 .OR. pcol(4) == 1) then
      reflect_electron(1,i) = abs(reflect_electron(1,i))
    endif
  

  else ! Other models. Here reflected electrons' direction is completely random.
    vf = sqrt(2.*taus88()*Ep/me)

    phi = taus88() * 2 * pi        ! Calculate emission azimuthal angle of secondaries

    ! Original code :
    ! thetaf = asin(2*taus88() - 1)   ! Calculate emission angle of secondaries with respect to the normal
    ! Modification A Placais 2017.09.26:
    thetaf = asin(sqrt(taus88()))
    ! Firstly, thetaf should be in [0,pi/2] rather than in [-pi/2,pi/2] (except if phis belongs to [0, pi])
    ! Secondly, asin(rt) will not give an uniform distribution. For furher details, see:
    ! John Greenwood, The correct and incorrect generation of a cosine distribution of scattered particles for Monte-Carlo
    ! modelling of vacuum systems, In Vacuum, Volume 67, Issue 2, 2002, Pages 217-222, ISSN 0042-207X,
    ! https://doi.org/10.1016/S0042-207X(02)00173-2.
    ! (http://www.sciencedirect.com/science/article/pii/S0042207X02001732)
    ! Keywords: Monte-Carlo; Cosine distribution; Gas scattering
    ! Generate a random velocity magnitude from 0 to pvf1 (random energy from 0 to Ep)

    !open(unit=666, file='../data/reflected_electrons.txt', status='unknown', position='append')
    !  write(666,*) thetaf, phi, -Ep/e, -vf*vf/(2*eme)
    !close(unit=666)
    
    reflect_electron(1,i) = vf * cos(thetaf)   ! Normal component of secondary velocity

    reflect_electron(1,j) = sin(phi) * vf * sin(thetaf)
    reflect_electron(1,3) = cos(phi) * vf * sin(thetaf)

    if (pcol(1) == 1 .OR. pcol(2) == 1) then
      reflect_electron(1,i) = - abs(reflect_electron(1,i))   ! Normal velocity pointing into the inside of the waveguide
    elseif (pcol(3) == 1 .OR. pcol(4) == 1) then
      reflect_electron(1,i) = abs(reflect_electron(1,i))   ! Normal velocity pointing into the inside of the waveguide
    endif

  endif
    
endif

return
end function reflect_electron

end module reflect
