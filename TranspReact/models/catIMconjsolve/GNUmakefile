TR_DIR ?= ../../
AMREX_HOME ?= ../../../../amrex
SUNDIALS_DIR?=/Users/hsitaram/projects/ccpc/multiscale-biomass-pyrolysis/TranspReact/sdbld

PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = TRUE
DEBUG      = FALSE

DIM       = 2

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE 
USE_CUDA   = FALSE
USE_SUNDIALS = FALSE
USE_HYPRE=TRUE
HYPRE_DIR=/Users/hsitaram/projects/ldrd_battery/echemAMR/hypre/src/hypre/

Bpack   := ./Make.package 
Blocs   := . 

ifeq ($(USE_SUNDIALS), TRUE)
    INCLUDE_LOCATIONS += ${SUNDIALS_DIR}/include
    LIBRARIES += -L$(SUNDIALS_DIR)/lib64 -lsundials_cvode -lsundials_arkode
    LIBRARIES += -L$(SUNDIALS_DIR)/lib64 -lsundials_core -lsundials_nvecserial -lsundials_nvecmanyvector
endif

include ${TR_DIR}/Make.tr
