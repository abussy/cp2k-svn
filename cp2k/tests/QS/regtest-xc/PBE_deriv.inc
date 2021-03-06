&DEBUG
&END DEBUG
&FORCE_EVAL
  &DFT
    &MGRID
      CUTOFF 200
    &END MGRID
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER WAVELET
    &END POISSON
    &QS
      EXTRAPOLATION USE_PREV_WF
    &END QS
    &SCF
      &DIAGONALIZATION
      &END DIAGONALIZATION
    &END SCF
    &XC
      &XC_FUNCTIONAL
        &PBE

@IF $FUNC == PBE
          PARAMETRIZATION ORIG
@ENDIF

@IF $FUNC == PBEsol
          PARAMETRIZATION PBESOL
@ENDIF

        &END PBE
      &END XC_FUNCTIONAL
    &END XC
    BASIS_SET_FILE_NAME GTH_BASIS_SETS
    MULTIPLICITY $MULT
    POTENTIAL_FILE_NAME GTH_POTENTIALS

@IF $MULT == 1
    UKS FALSE
@ENDIF

@IF $MULT == 3
    UKS TRUE
@ENDIF

  &END DFT
  &SUBSYS
    &CELL
      ABC 6.0 6.0 6.0
      PERIODIC NONE
    &END CELL
    &COORD
      H 3.0 3.0 2.6
      H 3.0 3.0 3.4
    &END COORD
    &KIND H
      BASIS_SET DZVP-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
&GLOBAL
  PRINT_LEVEL LOW
  PROJECT_NAME ${MULT}H2_${FUNC}
  RUN_TYPE DEBUG
&END GLOBAL

