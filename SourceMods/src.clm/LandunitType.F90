module LandunitType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Landunit data type allocation 
  ! -------------------------------------------------------- 
  ! landunits types can have values of (see landunit_varcon.F90)
  ! -------------------------------------------------------- 
  !   1  => (istsoil)    soil (vegetated or bare soil landunit)
  !   2  => (istcrop)    crop (only for crop configuration)
  !   3  => (UNUSED)     (formerly non-multiple elevation class land ice; currently unused)
  !   4  => (istice_mec) land ice (multiple elevation classes) 
  !   5  => (istdlak)    deep lake
  !   6  => (istwet)     wetland
  !   7  => (isturb_tbd) urban tbd
  !   8  => (isturb_hd)  urban hd
  !   9  => (isturb_md)  urban md
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use clm_varcon     , only : ispval
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  type, public :: landunit_type
     ! g/l/c/p hierarchy, local g/l/c/p cells only
     integer , pointer :: gridcell     (:) ! index into gridcell level quantities
     real(r8), pointer :: wtgcell      (:) ! weight (relative to gridcell)
     integer , pointer :: coli         (:) ! beginning column index per landunit
     integer , pointer :: colf         (:) ! ending column index for each landunit
     integer , pointer :: ncolumns     (:) ! number of columns for each landunit
     integer , pointer :: patchi       (:) ! beginning patch index for each landunit
     integer , pointer :: patchf       (:) ! ending patch index for each landunit
     integer , pointer :: npatches     (:) ! number of patches for each landunit

     ! topological mapping functionality
     integer , pointer :: itype        (:) ! landunit type
     logical , pointer :: ifspecial    (:) ! true=>landunit is not vegetated
     logical , pointer :: lakpoi       (:) ! true=>lake point
     logical , pointer :: urbpoi       (:) ! true=>urban point
     logical , pointer :: glcmecpoi    (:) ! true=>glacier_mec point
     logical , pointer :: active       (:) ! true=>do computations on this landunit 

     ! urban properties
     real(r8), pointer :: canyon_hwr   (:) ! urban landunit canyon height to width ratio (-)   
     real(r8), pointer :: wtroad_perv  (:) ! urban landunit weight of pervious road column to total road (-)
     real(r8), pointer :: wtlunit_roof (:) ! weight of roof with respect to urban landunit (-)
     real(r8), pointer :: ht_roof      (:) ! height of urban roof (m)
     real(r8), pointer :: z_0_town     (:) ! urban landunit momentum roughness length (m)
     real(r8), pointer :: z_d_town     (:) ! urban landunit displacement height (m)
!KO
     real(r8), pointer :: wall_to_plan_area_ratio (:) ! urban landunit ratio of wall area to plan area (-)
!KO

   contains

     procedure, public :: Init    ! Allocate and initialize
     procedure, public :: Clean   ! Clean up memory
     
  end type landunit_type
  ! Singleton instance of the landunitType
  type(landunit_type), public, target :: lun  !geomorphological landunits
  !------------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, begl, endl)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Allocate memory and initialize to signalling NaN to require
    ! data be properly initialized somewhere else.
    !
    ! !ARGUMENTS:
    class(landunit_type) :: this
    integer, intent(in) :: begl,endl
    !------------------------------------------------------------------------

    ! The following is set in InitGridCellsMod
    allocate(this%gridcell     (begl:endl)); this%gridcell  (:) = ispval
    allocate(this%wtgcell      (begl:endl)); this%wtgcell   (:) = nan
    allocate(this%coli         (begl:endl)); this%coli      (:) = ispval
    allocate(this%colf         (begl:endl)); this%colf      (:) = ispval
    allocate(this%ncolumns     (begl:endl)); this%ncolumns  (:) = ispval
    allocate(this%patchi       (begl:endl)); this%patchi    (:) = ispval
    allocate(this%patchf       (begl:endl)); this%patchf    (:) = ispval
    allocate(this%npatches     (begl:endl)); this%npatches  (:) = ispval
    allocate(this%itype        (begl:endl)); this%itype     (:) = ispval 
    allocate(this%ifspecial    (begl:endl)); this%ifspecial (:) = .false.
    allocate(this%lakpoi       (begl:endl)); this%lakpoi    (:) = .false.
    allocate(this%urbpoi       (begl:endl)); this%urbpoi    (:) = .false.
    allocate(this%glcmecpoi    (begl:endl)); this%glcmecpoi (:) = .false.

    ! The following is initialized in routine setActive in module reweightMod
    allocate(this%active       (begl:endl))

    ! The following is set in routine urbanparams_inst%Init in module UrbanParamsType
    allocate(this%canyon_hwr   (begl:endl)); this%canyon_hwr   (:) = nan
    allocate(this%wtroad_perv  (begl:endl)); this%wtroad_perv  (:) = nan
    allocate(this%ht_roof      (begl:endl)); this%ht_roof      (:) = nan
    allocate(this%wtlunit_roof (begl:endl)); this%wtlunit_roof (:) = nan
    allocate(this%z_0_town     (begl:endl)); this%z_0_town     (:) = nan
    allocate(this%z_d_town     (begl:endl)); this%z_d_town     (:) = nan
!KO
    allocate(this%wall_to_plan_area_ratio (begl:endl)); this%wall_to_plan_area_ratio (:) = nan
!KO

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine Clean(this)
    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Clean up memory use
    !
    ! !ARGUMENTS:
    class(landunit_type) :: this
    !------------------------------------------------------------------------

    deallocate(this%gridcell     )
    deallocate(this%wtgcell      )
    deallocate(this%coli         )
    deallocate(this%colf         )
    deallocate(this%ncolumns     )
    deallocate(this%patchi       )
    deallocate(this%patchf       )
    deallocate(this%npatches     )
    deallocate(this%itype        )
    deallocate(this%ifspecial    )
    deallocate(this%lakpoi       )
    deallocate(this%urbpoi       )
    deallocate(this%glcmecpoi    )
    deallocate(this%active       )
    deallocate(this%canyon_hwr   )
    deallocate(this%wtroad_perv  )
    deallocate(this%ht_roof      )
    deallocate(this%wtlunit_roof )
    deallocate(this%z_0_town     )
    deallocate(this%z_d_town     )
!KO
    deallocate(this%wall_to_plan_area_ratio)
!KO

  end subroutine Clean

end module LandunitType
