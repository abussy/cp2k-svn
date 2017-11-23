!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
 
! **************************************************************************************************
!> \brief Methods for X-Ray absorption spectroscopy (XAS) using TDDFPT  
!> \author AB (11.2017)
! **************************************************************************************************

MODULE xas_tdp_methods

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              deallocate_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_constants,                 ONLY: xas_tdp_by_index,&
                                              xas_tdp_by_kind
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_get_subs_vals
   USE qs_environment_types,            ONLY: qs_environment_type,&
                                              get_qs_env
   USE xas_tdp_types,                   ONLY: xas_tdp_control_type,&
                                              read_xas_tdp_control,&
                                              xas_tdp_control_create,&
                                              xas_tdp_env_type,&
                                              xas_tdp_env_create,&
                                              xas_tdp_env_release,&
                                              set_xas_tdp_env
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xas_tdp_methods'

   PUBLIC :: xas_tdp

CONTAINS
                
! **************************************************************************************************
!> \brief Driver for XAS TDDFT calculations.
!> \param qs_env the inherited qs_environment
!> \param dft_control the set of input DFT keywords and values
!> \author AB
!> \note Empty for now...
! **************************************************************************************************
   SUBROUTINE xas_tdp(qs_env, dft_control)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'xas_tdp', routineP = moduleN//':'//routineN

      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: handle, output_unit
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control

      CALL timeset(routineN, handle)

      NULLIFY(xas_tdp_env, xas_tdp_control)
      
      NULLIFY (logger)
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      PRINT*, "Saying 'hi' from the ", routineP, " routine!"

!   Initialization
      CALL xas_tdp_init(xas_tdp_env,xas_tdp_control,qs_env,dft_control,logger)

!   Clean-up
      CALL xas_tdp_env_release(xas_tdp_env)

      CALL timestop(handle)
   END SUBROUTINE xas_tdp

!> *************************************************************************************************
!> \brief Create and initialize relevant types for TDP XAS calculations
!> \param xas_tdp_env the environment type to create/initialize
!> \param xas_tdp_control the control typer to create/initialize
!> \param qs_env the inherited qs_env
!> \param dft_control the control type where xas_tdp_control lives
!> \param logger ...
!> *************************************************************************************************
   SUBROUTINE xas_tdp_init(xas_tdp_env,xas_tdp_control,qs_env,dft_control,logger)
   
      TYPE(xas_tdp_env_type), POINTER                    :: xas_tdp_env
      TYPE(xas_tdp_control_type), POINTER                :: xas_tdp_control
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_logger_type), POINTER                      :: logger

      CHARACTER(len=*), PARAMETER :: routineN = 'xas_tdp_init', routineP = moduleN//':'//routineN

      TYPE(section_vals_type), POINTER                   :: dft_section, xas_tdp_section
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set, ex_atomic_kind_set
      INTEGER                                            :: i,nex_kinds,nex_atoms, natoms
      INTEGER, DIMENSION(:), POINTER                     :: indices
      CHARACTER(len=2)                                   :: symbol

      NULLIFY(dft_section,xas_tdp_section,atomic_kind_set, ex_atomic_kind_set)

!   Create and initialize the XAS TDP control type
      dft_section => section_vals_get_subs_vals(qs_env%input, "DFT")                                 
      xas_tdp_section => section_vals_get_subs_vals(dft_section, "XAS_TDP")                          
                                                                                                                 
      CALL xas_tdp_control_create(dft_control%xas_tdp_control)                                       
      CALL read_xas_tdp_control(dft_control%xas_tdp_control,xas_tdp_section)
      
      xas_tdp_control => dft_control%xas_tdp_control

!   Create and initialize the XAS TDP environment type
      CALL xas_tdp_env_create(xas_tdp_env)

!     Retrieving the excited atomic kind set
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)

      IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_kind) THEN
         ALLOCATE(ex_atomic_kind_set(SIZE(xas_tdp_control%ex_kinds_symbols)))
         ALLOCATE(indices(SIZE(xas_tdp_control%ex_kinds_symbols)))
         nex_kinds = 0
         nex_atoms = 0
         DO i = 1,SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind=atomic_kind_set(i), element_symbol = symbol,&
                                 natom=natoms)
            IF (ANY(xas_tdp_control%ex_kinds_symbols .EQ. symbol)) THEN
               nex_kinds = nex_kinds + 1
               nex_atoms = nex_atoms + natoms
               indices(nex_kinds) = i
            END IF
         END DO
!         ex_atomic_kind_set = atomic_kind_set
         PRINT*, "Test :", indices
!         CALL set_xas_tdp_env(xas_tdp_env=xas_tdp_env, nex_atoms = nex_atoms,  nex_kinds=nex_kinds,&
!                              ex_atomic_kind_set=ex_atomic_kind_set)
         DEALLOCATE(indices)
         DEALLOCATE(ex_atomic_kind_set)

      ELSE IF (xas_tdp_control%define_excited .EQ. xas_tdp_by_index) THEN
         PRINT*, "asdf"
      END IF

!     for testing purposes
      CPASSERT(ASSOCIATED(logger))

   END SUBROUTINE xas_tdp_init
END MODULE xas_tdp_methods
