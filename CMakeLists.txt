##
#  CMake script
##

MESSAGE("_________________________________________________________________________________")
MESSAGE("                                                                                 ")
MESSAGE("                                   ExaDG                                         ")
MESSAGE("            High-order Discontinuous Galerkin for the Exa-Scale                  ")
MESSAGE("_________________________________________________________________________________")


# this is the standard deal.II search mechanism, including check for Trilinos and p4est

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.4 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

IF(NOT DEAL_II_WITH_P4EST)
  MESSAGE(FATAL_ERROR
    "\nExaDG requires a deal.II installation built with support for p4est but it appears to be missing!\n"
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

include_directories(applications/macros)

# Set the source files to be compiled
SET( TARGET_SRC
     include/time_integration/bdf_time_integration.cpp
     include/time_integration/extrapolation_scheme.cpp
     include/time_integration/time_int_base.cpp
     include/time_integration/time_int_bdf_base.cpp
     include/time_integration/time_int_explicit_runge_kutta_base.cpp
     include/time_integration/enum_types.cpp
     include/matrix_free/resolve_templates_float_2d.cpp
     include/matrix_free/resolve_templates_float_3d.cpp
     include/matrix_free/resolve_templates_double_2d.cpp
     include/matrix_free/resolve_templates_double_3d.cpp
     include/solvers_and_preconditioners/preconditioner/enum_types.cpp
     include/solvers_and_preconditioners/solvers/enum_types.cpp
     include/solvers_and_preconditioners/multigrid/multigrid_preconditioner_base.cpp
     include/solvers_and_preconditioners/multigrid/multigrid_input_parameters.cpp     
     include/solvers_and_preconditioners/multigrid/transfer/mg_transfer_mf_p.cpp
     include/solvers_and_preconditioners/multigrid/transfer/mg_transfer_mf_c.cpp
     include/postprocessor/error_calculation.cpp
     include/postprocessor/lift_and_drag_calculation.cpp
     include/postprocessor/pressure_difference_calculation.cpp
     include/postprocessor/kinetic_energy_spectrum.cpp
     include/postprocessor/kinetic_energy_calculation.cpp
     include/postprocessor/statistics_manager.cpp 
     include/poisson/user_interface/enum_types.cpp
     include/poisson/user_interface/input_parameters.cpp
     include/poisson/spatial_discretization/operator.cpp
     include/poisson/spatial_discretization/laplace_operator.cpp     
     include/convection_diffusion/user_interface/enum_types.cpp
     include/convection_diffusion/user_interface/input_parameters.cpp
     include/convection_diffusion/spatial_discretization/operators/convective_operator.cpp
     include/convection_diffusion/spatial_discretization/operators/diffusive_operator.cpp
     include/convection_diffusion/spatial_discretization/operators/mass_matrix_operator.cpp
     include/convection_diffusion/spatial_discretization/operators/rhs_operator.cpp
     include/convection_diffusion/spatial_discretization/operators/combined_operator.cpp
     include/convection_diffusion/spatial_discretization/dg_operator.cpp
     include/convection_diffusion/spatial_discretization/interface.cpp
     include/convection_diffusion/time_integration/time_int_bdf.cpp
     include/convection_diffusion/time_integration/time_int_explicit_runge_kutta.cpp
     include/convection_diffusion/time_integration/driver_steady_problems.cpp
     include/convection_diffusion/postprocessor/postprocessor.cpp
     include/convection_diffusion/postprocessor/output_generator.cpp
     include/incompressible_navier_stokes/user_interface/input_parameters.cpp
     include/incompressible_navier_stokes/user_interface/enum_types.cpp
     include/incompressible_navier_stokes/spatial_discretization/operators/mass_matrix_operator.cpp     
     include/incompressible_navier_stokes/spatial_discretization/operators/convective_operator.cpp 
     include/incompressible_navier_stokes/spatial_discretization/operators/viscous_operator.cpp 
     include/incompressible_navier_stokes/spatial_discretization/operators/gradient_operator.cpp 
     include/incompressible_navier_stokes/spatial_discretization/operators/divergence_operator.cpp
     include/incompressible_navier_stokes/spatial_discretization/operators/rhs_operator.cpp 
     include/incompressible_navier_stokes/spatial_discretization/operators/divergence_penalty_operator.cpp 
     include/incompressible_navier_stokes/spatial_discretization/operators/continuity_penalty_operator.cpp 
     include/incompressible_navier_stokes/spatial_discretization/operators/momentum_operator.cpp 
     include/incompressible_navier_stokes/spatial_discretization/operators/projection_operator.cpp 
     include/incompressible_navier_stokes/spatial_discretization/turbulence_model.cpp 
     include/incompressible_navier_stokes/spatial_discretization/dg_navier_stokes_base.cpp
     include/incompressible_navier_stokes/spatial_discretization/dg_projection_methods.cpp
     include/incompressible_navier_stokes/spatial_discretization/dg_dual_splitting.cpp
     include/incompressible_navier_stokes/spatial_discretization/dg_pressure_correction.cpp
     include/incompressible_navier_stokes/spatial_discretization/dg_coupled_solver.cpp
     include/incompressible_navier_stokes/spatial_discretization/interface.cpp 
     include/incompressible_navier_stokes/spatial_discretization/calculators/divergence_calculator.cpp 
     include/incompressible_navier_stokes/spatial_discretization/calculators/vorticity_calculator.cpp 
     include/incompressible_navier_stokes/spatial_discretization/calculators/velocity_magnitude_calculator.cpp 
     include/incompressible_navier_stokes/spatial_discretization/calculators/q_criterion_calculator.cpp 
     include/incompressible_navier_stokes/spatial_discretization/calculators/streamfunction_calculator_rhs_operator.cpp 
     include/incompressible_navier_stokes/preconditioners/compatible_laplace_operator.cpp
     include/incompressible_navier_stokes/time_integration/time_int_bdf_navier_stokes.cpp
     include/incompressible_navier_stokes/time_integration/time_int_bdf_coupled_solver.cpp
     include/incompressible_navier_stokes/time_integration/time_int_bdf_dual_splitting.cpp
     include/incompressible_navier_stokes/time_integration/time_int_bdf_pressure_correction.cpp
     include/incompressible_navier_stokes/time_integration/driver_steady_problems.cpp
     include/incompressible_navier_stokes/postprocessor/output_generator.cpp 
     include/incompressible_navier_stokes/postprocessor/divergence_and_mass_error.cpp 
     include/incompressible_navier_stokes/postprocessor/inflow_data_calculator.cpp 
     include/incompressible_navier_stokes/postprocessor/kinetic_energy_dissipation_detailed.cpp 
     include/incompressible_navier_stokes/postprocessor/line_plot_calculation.cpp 
     include/incompressible_navier_stokes/postprocessor/line_plot_calculation_statistics.cpp 
     include/incompressible_navier_stokes/postprocessor/perturbation_energy_orr_sommerfeld.cpp 
     include/incompressible_navier_stokes/postprocessor/mean_velocity_calculator.cpp 
     include/incompressible_navier_stokes/postprocessor/flow_rate_calculator.cpp 
     include/incompressible_navier_stokes/postprocessor/postprocessor.cpp 
     include/compressible_navier_stokes/user_interface/enum_types.cpp
     include/compressible_navier_stokes/user_interface/input_parameters.cpp
     include/compressible_navier_stokes/time_integration/time_int_explicit_runge_kutta.cpp
     include/compressible_navier_stokes/spatial_discretization/interface.cpp
     include/compressible_navier_stokes/spatial_discretization/dg_operator.cpp
     include/compressible_navier_stokes/postprocessor/write_output.cpp
     include/compressible_navier_stokes/postprocessor/postprocessor.cpp
  )

# Set the include directory and the name of the project
INCLUDE_DIRECTORIES(include)

PROJECT(exadg)

ADD_LIBRARY(exadg  ${TARGET_SRC})

# Define custom targets to easily switch the build type:
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug\nPlease call 'make' to build the debug project"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMENT "Switch CMAKE_BUILD_TYPE to Release\nPlease call 'make' to build the release project"
  )

DEAL_II_INITIALIZE_CACHED_VARIABLES()

DEAL_II_SETUP_TARGET(exadg)

ADD_SUBDIRECTORY(applications)
ADD_SUBDIRECTORY(tests)


