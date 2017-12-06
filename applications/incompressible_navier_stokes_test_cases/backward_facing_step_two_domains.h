/*
 * TurbulentChannel.h
 *
 *  Created on: Oct 14, 2016
 *      Author: fehn
 */

#ifndef APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_TURBULENT_CHANNEL_H_
#define APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_TURBULENT_CHANNEL_H_

#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>

/**************************************************************************************/
/*                                                                                    */
/*                                 INPUT PARAMETERS                                   */
/*                                                                                    */
/**************************************************************************************/

// single or double precision?
//typedef float VALUE_TYPE;
typedef double VALUE_TYPE;

// set the number of space dimensions: dimension = 2, 3
unsigned int const DIMENSION = 3;

// set the polynomial degree of the shape functions for velocity and pressure.
// currently, one has to use the same polynomial degree for both domains.
unsigned int const FE_DEGREE_VELOCITY = 5;
unsigned int const FE_DEGREE_PRESSURE = FE_DEGREE_VELOCITY-1;

// set xwall specific parameters
unsigned int const FE_DEGREE_XWALL = 1;
unsigned int const N_Q_POINTS_1D_XWALL = 1;

// set the number of refine levels for DOMAIN 1
unsigned int const REFINE_STEPS_SPACE_DOMAIN1 = 3;

// set the number of refine levels for DOMAIN 2
unsigned int const REFINE_STEPS_SPACE_DOMAIN2 = 3;

// set the number of refine levels for temporal convergence tests
unsigned int const REFINE_STEPS_TIME_MIN = 0;
unsigned int const REFINE_STEPS_TIME_MAX = REFINE_STEPS_TIME_MIN;

// set problem specific parameters like physical dimensions, etc.
double const PI = numbers::PI;

// Height H
double const H = 0.041;

// channel
double const LENGTH_CHANNEL = 2.0*PI*H;
double const HEIGHT_CHANNEL = 2.0*H;
double const WIDTH_CHANNEL = 4.0*H;

// use a gap between both geometries for visualization purposes
double const GAP_CHANNEL_BFS = 2.0*H;

// backward facing step geometry
double const LENGTH_BFS_DOWN = 20.0*H;
double const LENGTH_BFS_UP = 2.0*H;
double const HEIGHT_BFS_STEP = H;
double const HEIGHT_BFS_INFLOW = HEIGHT_CHANNEL;
double const WIDTH_BFS = WIDTH_CHANNEL;

double const X1_COORDINATE_INFLOW = - LENGTH_BFS_UP;
double const X1_COORDINATE_OUTFLOW = LENGTH_BFS_DOWN;
double const X1_COORDINATE_OUTFLOW_CHANNEL = - LENGTH_BFS_UP - GAP_CHANNEL_BFS;

// consider a friction Reynolds number of Re_tau = 290 = u_tau * H / nu
//and body force f = tau_w/H with tau_w = u_tau^2.
double const VISCOSITY = 1.5268e-5;

// estimate the maximum velocity
double const MAX_VELOCITY = 2.0;

double const START_TIME = 0.0;
double const SAMPLE_START_TIME = 2.0;
double const END_TIME = 6.0;

QuantityStatistics QUANTITY_VELOCITY;
QuantityStatisticsSkinFriction<3> QUANTITY_SKIN_FRICTION;
QuantityStatistics QUANTITY_REYNOLDS;
QuantityStatistics QUANTITY_PRESSURE;
QuantityStatisticsPressureCoefficient<3> QUANTITY_PRESSURE_COEFF;
const unsigned int N_POINTS_LINE = 101;

// use a negative GRID_STRETCH_FAC to deactivate grid stretching
const double GRID_STRETCH_FAC = 1.8;

std::string OUTPUT_FOLDER = "output/bfs/test/";
std::string OUTPUT_FOLDER_VTU = OUTPUT_FOLDER + "vtu/";
std::string OUTPUT_NAME_1 = "precursor";
std::string OUTPUT_NAME_2 = "bfs";

// DOMAIN 1: turbulent channel problem: used to generate inflow data for the BFS
// DOMAIN 2: backward facing step (using results of the turbulent channel flow as velocity inflow profile)

// data structures that we need to apply the velocity inflow profile
// we currently use global variables for this purpose
unsigned int N_POINTS_Y = 101;
unsigned int N_POINTS_Z = N_POINTS_Y;
std::vector<double> Y_VALUES(N_POINTS_Y);
std::vector<double> Z_VALUES(N_POINTS_Z);
std::vector<Tensor<1,DIMENSION,double> > VELOCITY_VALUES(N_POINTS_Y*N_POINTS_Z);

// initial vectors
void initialize_y_and_z_values()
{
  AssertThrow(N_POINTS_Y >= 2, ExcMessage("Variable N_POINTS_Y is invalid"));
  AssertThrow(N_POINTS_Z >= 2, ExcMessage("Variable N_POINTS_Z is invalid"));

  for(unsigned int i=0; i<N_POINTS_Y; ++i)
    Y_VALUES[i] = double(i)/double(N_POINTS_Y-1)*HEIGHT_CHANNEL;

  for(unsigned int i=0; i<N_POINTS_Z; ++i)
    Z_VALUES[i] = -WIDTH_CHANNEL/2.0 + double(i)/double(N_POINTS_Z-1)*WIDTH_CHANNEL;
}

void initialize_velocity_values()
{
  AssertThrow(N_POINTS_Y >= 2, ExcMessage("Variable N_POINTS_Y is invalid"));
  AssertThrow(N_POINTS_Z >= 2, ExcMessage("Variable N_POINTS_Z is invalid"));

  for(unsigned int iy=0; iy<N_POINTS_Y; ++iy)
  {
    for(unsigned int iz=0; iz<N_POINTS_Z; ++iz)
    {
      Tensor<1,DIMENSION,double> velocity;
      VELOCITY_VALUES[iy*N_POINTS_Y + iz] = velocity;
    }
  }
}

// we do not need this function here (but have to implement it)
template<int dim>
void InputParametersNavierStokes<dim>::set_input_parameters()
{

}

/*
 *  To set input parameters for DOMAIN 1 and DOMAIN 2, use
 *
 *  if(domain_id == 1){}
 *  else if(domain_id == 2){}
 *
 *  Most of the input parameters are the same for both domains!
 */
template<int dim>
void InputParametersNavierStokes<dim>::set_input_parameters(unsigned int const domain_id)
{
  // MATHEMATICAL MODEL
  problem_type = ProblemType::Unsteady;
  equation_type = EquationType::NavierStokes;
  use_outflow_bc_convective_term = true;
  formulation_viscous_term = FormulationViscousTerm::DivergenceFormulation;
  right_hand_side = true;


  // PHYSICAL QUANTITIES
  start_time = START_TIME;
  end_time = END_TIME;
  viscosity = VISCOSITY;


  // TEMPORAL DISCRETIZATION
  solver_type = SolverType::Unsteady;
  temporal_discretization = TemporalDiscretization::BDFDualSplittingScheme; // BDFDualSplittingScheme; //BDFPressureCorrection; //BDFCoupledSolution;
  treatment_of_convective_term = TreatmentOfConvectiveTerm::Explicit; //Explicit;
  calculation_of_time_step_size = TimeStepCalculation::ConstTimeStepCFL; // AdaptiveTimeStepCFL
  max_velocity = MAX_VELOCITY;
  cfl = 0.15;
  cfl_exponent_fe_degree_velocity = 1.5;
  time_step_size = 1.0e-1;
  max_number_of_time_steps = 1e8;
  order_time_integrator = 2; // 1; // 2; // 3;
  start_with_low_order = true;


  // SPATIAL DISCRETIZATION

  // spatial discretization method
  spatial_discretization = SpatialDiscretization::DG;

  // convective term - currently no parameters

  // viscous term
  IP_formulation_viscous = InteriorPenaltyFormulation::SIPG;
  IP_factor_viscous = 1.0;
  penalty_term_div_formulation = PenaltyTermDivergenceFormulation::Symmetrized;

  // gradient term
  gradp_integrated_by_parts = true;
  gradp_use_boundary_data = true;

  // divergence term
  divu_integrated_by_parts = true;
  divu_use_boundary_data = true;

  // special case: pure DBC's
  if(domain_id == 1)
    pure_dirichlet_bc = true;
  else if(domain_id == 2)
    pure_dirichlet_bc = false;

  // div-div and continuity penalty
  use_divergence_penalty = true;
  divergence_penalty_factor = 1.0e0;
  use_continuity_penalty = true;
  continuity_penalty_components = ContinuityPenaltyComponents::All;
  continuity_penalty_use_boundary_data = false;
  continuity_penalty_factor = divergence_penalty_factor;

  // TURBULENCE
  use_turbulence_model = false;
  turbulence_model = TurbulenceEddyViscosityModel::Sigma;
  // Smagorinsky: 0.165
  // Vreman: 0.28
  // WALE: 0.50
  // Sigma: 1.35
  turbulence_model_constant = 1.35;

  // PROJECTION METHODS

  // pressure Poisson equation
  IP_factor_pressure = 1.0;
  solver_pressure_poisson = SolverPressurePoisson::PCG;
  preconditioner_pressure_poisson = PreconditionerPressurePoisson::GeometricMultigrid;
  multigrid_data_pressure_poisson.smoother = MultigridSmoother::Chebyshev; //Chebyshev; //Jacobi; //GMRES;
  //Chebyshev
  multigrid_data_pressure_poisson.coarse_solver = MultigridCoarseGridSolver::Chebyshev;

  abs_tol_pressure = 1.e-12;
  rel_tol_pressure = 1.e-6;

  // stability in the limit of small time steps
  use_approach_of_ferrer = false;
  deltat_ref = 1.e0;

  // projection step
  solver_projection = SolverProjection::PCG;
  preconditioner_projection = PreconditionerProjection::InverseMassMatrix; //BlockJacobi; //PointJacobi; //InverseMassMatrix;
  update_preconditioner_projection = true;
  abs_tol_projection = 1.e-12;
  rel_tol_projection = 1.e-6;



  // HIGH-ORDER DUAL SPLITTING SCHEME

  // formulations
  order_extrapolation_pressure_nbc = order_time_integrator <=2 ? order_time_integrator : 2;

  // convective step

  // nonlinear solver
  newton_solver_data_convective.abs_tol = 1.e-12;
  newton_solver_data_convective.rel_tol = 1.e-6;
  newton_solver_data_convective.max_iter = 100;
  // linear solver
  abs_tol_linear_convective = 1.e-12;
  rel_tol_linear_convective = 1.e-6;
  max_iter_linear_convective = 1e4;
  use_right_preconditioning_convective = true;
  max_n_tmp_vectors_convective = 100;

  // stability in the limit of small time steps and projection step
  small_time_steps_stability = false;

  // viscous step
  solver_viscous = SolverViscous::PCG;
  preconditioner_viscous = PreconditionerViscous::InverseMassMatrix; //GeometricMultigrid;
  abs_tol_viscous = 1.e-12;
  rel_tol_viscous = 1.e-6;


  // PRESSURE-CORRECTION SCHEME

  // formulation
  order_pressure_extrapolation = 1; // use 0 for non-incremental formulation
  rotational_formulation = true;

  // momentum step

  // Newton solver
  newton_solver_data_momentum.abs_tol = 1.e-12;
  newton_solver_data_momentum.rel_tol = 1.e-6;
  newton_solver_data_momentum.max_iter = 100;

  // linear solver
  abs_tol_momentum_linear = 1.e-12;
  rel_tol_momentum_linear = 1.e-6;
  max_iter_momentum_linear = 1e4;
  use_right_preconditioning_momentum = true;
  max_n_tmp_vectors_momentum = 100;
  update_preconditioner_momentum = false;

  solver_momentum = SolverMomentum::GMRES;
  preconditioner_momentum = MomentumPreconditioner::InverseMassMatrix;


  // COUPLED NAVIER-STOKES SOLVER
  use_scaling_continuity = false;
  scaling_factor_continuity = 1.0;

  // nonlinear solver (Newton solver)
  newton_solver_data_coupled.abs_tol = 1.e-12;
  newton_solver_data_coupled.rel_tol = 1.e-6;
  newton_solver_data_coupled.max_iter = 1e2;

  // linear solver
  solver_linearized_navier_stokes = SolverLinearizedNavierStokes::GMRES; //GMRES; //FGMRES;
  abs_tol_linear = 1.e-12;
  rel_tol_linear = 1.e-6;
  max_iter_linear = 1e3;
  max_n_tmp_vectors = 100;

  // preconditioning linear solver
  preconditioner_linearized_navier_stokes = PreconditionerLinearizedNavierStokes::BlockTriangular;
  update_preconditioner = false;

  // preconditioner velocity/momentum block
  momentum_preconditioner = MomentumPreconditioner::InverseMassMatrix;

  // preconditioner Schur-complement block
  schur_complement_preconditioner = SchurComplementPreconditioner::CahouetChabard; //PressureConvectionDiffusion;
  discretization_of_laplacian =  DiscretizationOfLaplacian::Classical;

  // Chebyshev moother
  multigrid_data_schur_complement_preconditioner.smoother = MultigridSmoother::Chebyshev;
  multigrid_data_schur_complement_preconditioner.coarse_solver = MultigridCoarseGridSolver::Chebyshev;


  if(domain_id == 1)
  {
    // OUTPUT AND POSTPROCESSING
    print_input_parameters = true;

    // write output for visualization of results
    output_data.write_output = true;
    output_data.output_folder = OUTPUT_FOLDER_VTU;
    output_data.output_name = OUTPUT_NAME_1;
    output_data.output_start_time = start_time;
    output_data.output_interval_time = (end_time-start_time)/20;
    output_data.write_divergence = true;
    output_data.number_of_patches = FE_DEGREE_VELOCITY;

    // output of solver information
    output_solver_info_every_timesteps = 1e3;

    // turbulent channel statistics
    turb_ch_data.calculate_statistics = false;
    turb_ch_data.sample_start_time = SAMPLE_START_TIME;
    turb_ch_data.sample_end_time = END_TIME;
    turb_ch_data.sample_every_timesteps = 10;
    turb_ch_data.viscosity = VISCOSITY;
    turb_ch_data.filename_prefix = OUTPUT_FOLDER + OUTPUT_NAME_1;

    inflow_data.write_inflow_data = true;
    inflow_data.x_coordinate = X1_COORDINATE_OUTFLOW_CHANNEL;
    inflow_data.n_points_y = N_POINTS_Y;
    inflow_data.n_points_z = N_POINTS_Z;
    inflow_data.y_values = &Y_VALUES;
    inflow_data.z_values = &Z_VALUES;
    inflow_data.array = &VELOCITY_VALUES;
  }
  else if(domain_id == 2)
  {
    // OUTPUT AND POSTPROCESSING
    print_input_parameters = true;

    // write output for visualization of results
    output_data.write_output = true;
    output_data.output_folder = OUTPUT_FOLDER_VTU;
    output_data.output_name = OUTPUT_NAME_2;
    output_data.output_start_time = start_time;
    output_data.output_interval_time = (end_time-start_time)/20;
    output_data.write_divergence = true;
    output_data.number_of_patches = FE_DEGREE_VELOCITY;

    // output of solver information
    output_solver_info_every_timesteps = 1e3;

    // turbulent channel statistics
    turb_ch_data.calculate_statistics = false;
    turb_ch_data.sample_start_time = SAMPLE_START_TIME;
    turb_ch_data.sample_end_time = END_TIME;
    turb_ch_data.sample_every_timesteps = 10;
    turb_ch_data.viscosity = VISCOSITY;
    turb_ch_data.filename_prefix = OUTPUT_FOLDER + OUTPUT_NAME_2;

    bfs_statistics.calculate_statistics = true;
    bfs_statistics.sample_start_time = SAMPLE_START_TIME;
    bfs_statistics.sample_end_time = end_time;
    bfs_statistics.sample_every_timesteps = 10;
    bfs_statistics.filename_prefix = OUTPUT_FOLDER + OUTPUT_NAME_2;

    QUANTITY_VELOCITY.type = QuantityType::Velocity;
    QUANTITY_VELOCITY.averaging_direction = 2;

    QUANTITY_REYNOLDS.type = QuantityType::ReynoldsStresses;
    QUANTITY_REYNOLDS.averaging_direction = 2;

    QUANTITY_PRESSURE.type = QuantityType::Pressure;
    QUANTITY_PRESSURE.averaging_direction = 2;
    
    QUANTITY_PRESSURE_COEFF.type = QuantityType::PressureCoefficient;
    QUANTITY_PRESSURE_COEFF.averaging_direction = 2;
    QUANTITY_PRESSURE_COEFF.reference_velocity = 1.0;
    QUANTITY_PRESSURE_COEFF.reference_point = Point<DIMENSION>(X1_COORDINATE_INFLOW,0,0); //Jovic&Driver p.17

    Tensor<1,dim,double> normal; normal[1] = 1.0;
    Tensor<1,dim,double> tangent; tangent[0] = 1.0;
    QUANTITY_SKIN_FRICTION.type = QuantityType::SkinFriction;
    QUANTITY_SKIN_FRICTION.averaging_direction = 2;
    QUANTITY_SKIN_FRICTION.reference_velocity = 1.0;
    QUANTITY_SKIN_FRICTION.normal_vector = normal;
    QUANTITY_SKIN_FRICTION.tangent_vector = tangent;
    QUANTITY_SKIN_FRICTION.viscosity = VISCOSITY;

    Line<dim> vel_1, vel_2, vel_3, vel_4, vel_5, vel_6, vel_7, vel_8, vel_9, vel_10, vel_11, Cp_1, Cp_2, Cf;

    vel_1.begin = Point<dim> (0*H,   0,0);
    vel_1.end =   Point<dim> (0*H, 2*H,0);
    vel_2.begin = Point<dim> (1*H,-1*H,0);
    vel_2.end =   Point<dim> (1*H, 2*H,0);
    vel_3.begin = Point<dim> (2*H,-1*H,0);
    vel_3.end =   Point<dim> (2*H, 2*H,0);
    vel_4.begin = Point<dim> (3*H,-1*H,0);
    vel_4.end =   Point<dim> (3*H, 2*H,0);
    vel_5.begin = Point<dim> (4*H,-1*H,0);
    vel_5.end =   Point<dim> (4*H, 2*H,0);
    vel_6.begin = Point<dim> (5*H,-1*H,0);
    vel_6.end =   Point<dim> (5*H, 2*H,0);
    vel_7.begin = Point<dim> (6*H,-1*H,0);
    vel_7.end =   Point<dim> (6*H, 2*H,0);
    vel_8.begin = Point<dim> (7*H,-1*H,0);
    vel_8.end =   Point<dim> (7*H, 2*H,0);
    vel_9.begin = Point<dim> (8*H,-1*H,0);
    vel_9.end =   Point<dim> (8*H, 2*H,0);
    vel_10.begin = Point<dim> (9*H,-1*H,0);
    vel_10.end =   Point<dim> (9*H, 2*H,0);
    vel_11.begin = Point<dim> (10*H,-1*H,0);
    vel_11.end =   Point<dim> (10*H, 2*H,0);
    
    Cp_1.begin = Point<dim> (X1_COORDINATE_INFLOW,0,0);
    Cp_1.end =   Point<dim> (0,0,0);
    Cp_2.begin = Point<dim> (0,-H,0);
    Cp_2.end =   Point<dim> (X1_COORDINATE_OUTFLOW,-H,0);
    Cf.begin = Point<dim> (0,-H,0);
    Cf.end =   Point<dim> (X1_COORDINATE_OUTFLOW,-H,0);


    vel_1.n_points = N_POINTS_LINE;
    vel_2.n_points = N_POINTS_LINE;
    vel_3.n_points = N_POINTS_LINE;
    vel_4.n_points = N_POINTS_LINE;
    vel_5.n_points = N_POINTS_LINE;
    vel_6.n_points = N_POINTS_LINE;
    vel_7.n_points = N_POINTS_LINE;
    vel_8.n_points = N_POINTS_LINE;
    vel_9.n_points = N_POINTS_LINE;
    vel_10.n_points = N_POINTS_LINE;
    vel_11.n_points = N_POINTS_LINE;
    Cp_1.n_points = N_POINTS_LINE;
    Cp_2.n_points = N_POINTS_LINE;
    Cf.n_points = N_POINTS_LINE;

    vel_1.quantities.push_back(&QUANTITY_VELOCITY);
    vel_1.quantities.push_back(&QUANTITY_REYNOLDS);
    vel_2.quantities.push_back(&QUANTITY_VELOCITY);
    vel_2.quantities.push_back(&QUANTITY_REYNOLDS);
    vel_3.quantities.push_back(&QUANTITY_VELOCITY);
    vel_3.quantities.push_back(&QUANTITY_REYNOLDS);
    vel_4.quantities.push_back(&QUANTITY_VELOCITY);
    vel_4.quantities.push_back(&QUANTITY_REYNOLDS);
    vel_5.quantities.push_back(&QUANTITY_VELOCITY);
    vel_5.quantities.push_back(&QUANTITY_REYNOLDS);
    vel_6.quantities.push_back(&QUANTITY_VELOCITY);
    vel_6.quantities.push_back(&QUANTITY_REYNOLDS);
    vel_7.quantities.push_back(&QUANTITY_VELOCITY);
    vel_7.quantities.push_back(&QUANTITY_REYNOLDS);
    vel_8.quantities.push_back(&QUANTITY_VELOCITY);
    vel_8.quantities.push_back(&QUANTITY_REYNOLDS);
    vel_9.quantities.push_back(&QUANTITY_VELOCITY);
    vel_9.quantities.push_back(&QUANTITY_REYNOLDS);
    vel_10.quantities.push_back(&QUANTITY_VELOCITY);
    vel_10.quantities.push_back(&QUANTITY_REYNOLDS);
    vel_11.quantities.push_back(&QUANTITY_VELOCITY);
    vel_11.quantities.push_back(&QUANTITY_REYNOLDS);
    Cp_1.quantities.push_back(&QUANTITY_PRESSURE);
    Cp_1.quantities.push_back(&QUANTITY_PRESSURE_COEFF);
    Cp_2.quantities.push_back(&QUANTITY_PRESSURE);
    Cp_2.quantities.push_back(&QUANTITY_PRESSURE_COEFF);
    Cf.quantities.push_back(&QUANTITY_SKIN_FRICTION);

    vel_1.name = "vel_1";
    vel_2.name = "vel_2";
    vel_3.name = "vel_3";
    vel_4.name = "vel_4";
    vel_5.name = "vel_5";
    vel_6.name = "vel_6";
    vel_7.name = "vel_7";
    vel_8.name = "vel_8";
    vel_9.name = "vel_9";
    vel_10.name = "vel_10";
    vel_11.name = "vel_11";
    Cp_1.name = "Cp_1";
    Cp_2.name = "Cp_2";
    Cf.name = "Cf";

    line_plot_data.lines.push_back(vel_1);
    line_plot_data.lines.push_back(vel_2);
    line_plot_data.lines.push_back(vel_3);
    line_plot_data.lines.push_back(vel_4);
    line_plot_data.lines.push_back(vel_5);
    line_plot_data.lines.push_back(vel_6);
    line_plot_data.lines.push_back(vel_7);
    line_plot_data.lines.push_back(vel_8);
    line_plot_data.lines.push_back(vel_9);
    line_plot_data.lines.push_back(vel_10);
    line_plot_data.lines.push_back(vel_11);
    line_plot_data.lines.push_back(Cp_1);
    line_plot_data.lines.push_back(Cp_2);
    line_plot_data.lines.push_back(Cf);

    line_plot_data.write_output = true;

    mean_velocity_data.calculate_statistics = true;
    mean_velocity_data.sample_start_time = bfs_statistics.sample_start_time;
    mean_velocity_data.sample_end_time = bfs_statistics.sample_end_time;
    mean_velocity_data.sample_every_timesteps = bfs_statistics.sample_every_timesteps;
    mean_velocity_data.filename_prefix = OUTPUT_FOLDER + OUTPUT_NAME_2;
    Tensor<1,dim> normal_vector;
    normal_vector[0] = 1;
    mean_velocity_data.normal_vector = normal_vector;
    mean_velocity_data.boundary_IDs.insert(2);
    mean_velocity_data.area = HEIGHT_CHANNEL * WIDTH_CHANNEL;
  }
}

/**************************************************************************************/
/*                                                                                    */
/*    FUNCTIONS (ANALYTICAL SOLUTION, BOUNDARY CONDITIONS, VELOCITY FIELD, etc.)      */
/*                                                                                    */
/**************************************************************************************/

template<int dim>
class InitialSolutionVelocity : public Function<dim>
{
public:
  InitialSolutionVelocity (const unsigned int  n_components = dim,
                           const double        time = 0.)
    :
    Function<dim>(n_components, time)
  {
  }

  virtual ~InitialSolutionVelocity(){};

  virtual double value (const Point<dim>    &p,
                        const unsigned int  component = 0) const;
};

template<int dim>
double InitialSolutionVelocity<dim>::value(const Point<dim>   &p,
                                           const unsigned int component) const
{
  double result = 0.0;
  double const y = (p[1] - HEIGHT_CHANNEL/2.0)/(HEIGHT_CHANNEL/2.0);

  if(dim==3)
  {
    if(component == 0)
    {
      double factor = 0.5;

      if(std::abs(y)<1.0)
        result = -MAX_VELOCITY*(pow(y,6.0)-1.0)*(1.0+((double)rand()/RAND_MAX-0.5)*factor);
      else
        result = 0.0;
    }

//    double factor = 0.5;
//    double const x = p[0]/LENGTH_CHANNEL;
//    double const z = p[2]/WIDTH_CHANNEL;
//
//    if(component == 0)
//    {
//      if(std::abs(y)<1.0)
//        result = -MAX_VELOCITY*(pow(y,6.0)-1.0)*(1.0 + (((double)rand()/RAND_MAX-1.0) + std::sin(z*8.)*0.5)*factor);
//    }
//    if(component == 2)
//    {
//      if(std::abs(y)<1.0)
//        result = -MAX_VELOCITY*(pow(y,6.0)-1.0)*std::sin(x*8.)*0.5*factor;
//    }
  }
  else
  {
    AssertThrow(false, ExcMessage("Dimension has to be dim==3."));
  }

  return result;
}

#include "../../include/incompressible_navier_stokes/postprocessor/inflow_data_calculator.h"

template<int dim>
class InflowProfile : public Function<dim>
{
public:
  InflowProfile (const unsigned int  n_components = dim,
                 const double        time = 0.)
    :
    Function<dim>(n_components, time)
  {
    initialize_y_and_z_values();
    initialize_velocity_values();
  }

  virtual ~InflowProfile(){};

  virtual double value (const Point<dim>    &p,
                        const unsigned int  component = 0) const;
};

template<int dim>
double InflowProfile<dim>::value(const Point<dim>   &p,
                                 const unsigned int component) const
{
  double result = linear_interpolation_2d_cartesian(p,Y_VALUES,Z_VALUES,VELOCITY_VALUES,component);

  return result;
}


/*
 *  Neumann boundary conditions for velocity
 *
 *  - Laplace formulation of viscous term
 *    -> prescribe velocity gradient (grad U)*n on Gamma_N
 *
 *  - Divergence formulation of viscous term
 *    -> prescribe (grad U + (grad U)^T)*n on Gamma_N
 */
template<int dim>
class NeumannBoundaryVelocity : public Function<dim>
{
public:
  NeumannBoundaryVelocity (const double time = 0.)
    :
    Function<dim>(dim, time)
  {}

  virtual ~NeumannBoundaryVelocity(){};

  virtual double value (const Point<dim> &p,const unsigned int component = 0) const;
};

template<int dim>
double NeumannBoundaryVelocity<dim>::value(const Point<dim> &/* p */,const unsigned int /* component */) const
{
  double result = 0.0;
  return result;
}

/*
 *  PressureBC_dudt:
 *
 *  This functions is only used when applying the high-order dual splitting scheme and
 *  is evaluated on Dirichlet boundaries (where the velocity is prescribed).
 *  Hence, this is the function that is set in the dirichlet_bc map of boundary_descriptor_pressure.
 *
 *  Note:
 *    When using a couples solution approach we do not have to evaluate something like
 *    pressure Neumann BC's on Dirichlet boundaries (we only have p⁺ = p⁻ on Dirichlet boundaries,
 *    i.e., no boundary data used). So it doesn't matter when writing this function into the
 *    dirichlet_bc map of boundary_descriptor_pressure because this function will never be evaluated
 *    in case of a coupled solution approach.
 *
 */

template<int dim>
class PressureBC_dudt : public Function<dim>
{
public:
  PressureBC_dudt (const double time = 0.)
    :
    Function<dim>(dim, time)
  {}

  virtual ~PressureBC_dudt(){};

  virtual double value (const Point<dim>    &p,
                        const unsigned int  component = 0) const;
};

template<int dim>
double PressureBC_dudt<dim>::value(const Point<dim>   &/* p */,
                                   const unsigned int /* component */) const
{
  double result = 0.0;

  return result;
}

/*
 *  Right-hand side function: Implements the body force vector occuring on the
 *  right-hand side of the momentum equation of the Navier-Stokes equations
 */
template<int dim>
 class RightHandSide : public Function<dim>
 {
 public:
   RightHandSide (const double time = 0.)
     :
     Function<dim>(dim, time)
   {}

   virtual ~RightHandSide(){};

   virtual double value (const Point<dim>    &p,
                         const unsigned int  component = 0) const;
 };

 template<int dim>
 double RightHandSide<dim>::value(const Point<dim>   &/* p */,
                                  const unsigned int component) const
 {
   double result = 0.0;

   //channel flow with periodic bc
   if(component==0)
     return 0.2844518;
   else
     return 0.0;

   return result;
 }


/**************************************************************************************/
/*                                                                                    */
/*         GENERATE GRID, SET BOUNDARY INDICATORS AND FILL BOUNDARY DESCRIPTOR        */
/*                                                                                    */
/**************************************************************************************/

// TODO
double grid_transform_y(const double &eta)
{
 double y = eta;

 return y;
}

// /*
//  *  maps eta in [0,1] --> y in [-1,1]*length_y/2.0 (using a hyperbolic mesh stretching)
//  */
//double grid_transform_y(const double &eta)
//{
//  double y = 0.0;
//
//  if(GRID_STRETCH_FAC >= 0)
//    y = DIMENSIONS_X2/2.0*std::tanh(GRID_STRETCH_FAC*(2.*eta-1.))/std::tanh(GRID_STRETCH_FAC);
//  else // use a negative GRID_STRETCH_FACto deactivate grid stretching
//    y = DIMENSIONS_X2/2.0*(2.*eta-1.);
//
//  return y;
//}
//
///*
// * inverse mapping:
// *
// *  maps y in [-1,1]*length_y/2.0 --> eta in [0,1]
// */
//double inverse_grid_transform_y(const double &y)
//{
//  double eta = 0.0;
//
//  if(GRID_STRETCH_FAC >= 0)
//    eta = (std::atanh(y*std::tanh(GRID_STRETCH_FAC)*2.0/DIMENSIONS_X2)/GRID_STRETCH_FAC+1.0)/2.0;
//  else // use a negative GRID_STRETCH_FACto deactivate grid stretching
//    eta = (2.*y/DIMENSIONS_X2+1.)/2.0;
//
//  return eta;
//}
//
//template <int dim>
//Point<dim> grid_transform (const Point<dim> &in)
//{
//  Point<dim> out = in;
//
//  out[0] = in(0)-numbers::PI;
//  out[1] = grid_transform_y(in[1]);
//
//  if(dim==3)
//    out[2] = in(2)-0.5*numbers::PI;
//  return out;
//}
//

//#include <deal.II/grid/manifold_lib.h>
//
//template <int dim>
//class ManifoldTurbulentChannel : public ChartManifold<dim,dim,dim>
//{
//public:
//  ManifoldTurbulentChannel(Tensor<1,dim> &dimensions_in)
//  {
//    dimensions = dimensions_in;
//  }
//
//  /*
//   *  push_forward operation that maps point xi in reference coordinates [0,1]^d to
//   *  point x in physical coordinates
//   */
//  Point<dim> push_forward(const Point<dim> &xi) const
//  {
//    Point<dim> x;
//
//    x[0] = xi[0]*dimensions[0]-dimensions[0]/2.0;
//    x[1] = grid_transform_y(xi[1]);
//
//    if(dim==3)
//      x[2] = xi[2]*dimensions[2]-dimensions[2]/2.0;
//
//    return x;
//  }
//
//  /*
//   *  pull_back operation that maps point x in physical coordinates
//   *  to point xi in reference coordinates [0,1]^d
//   */
//  Point<dim> pull_back(const Point<dim> &x) const
//  {
//    Point<dim> xi;
//
//    xi[0] = x[0]/dimensions[0]+0.5;
//    xi[1] = inverse_grid_transform_y(x[1]);
//
//    if(dim==3)
//      xi[2] = x[2]/dimensions[2]+0.5;
//
//    return xi;
//  }
//
//private:
// Tensor<1,dim> dimensions;
//};


template<int dim>
void create_grid_and_set_boundary_conditions_1(
    parallel::distributed::Triangulation<dim>              &triangulation,
    unsigned int const                                     n_refine_space,
    std::shared_ptr<BoundaryDescriptorNavierStokesU<dim> > boundary_descriptor_velocity,
    std::shared_ptr<BoundaryDescriptorNavierStokesP<dim> > boundary_descriptor_pressure,
    std::vector<GridTools::PeriodicFacePair<typename
      Triangulation<dim>::cell_iterator> >                 &periodic_faces)
{
  /* --------------- Generate grid ------------------- */
  if(dim==2)
  {
    AssertThrow(false, ExcMessage("NotImplemented"));
  }
  else if(dim==3)
  {
    Tensor<1,dim> dimensions;
    dimensions[0] = LENGTH_CHANNEL;
    dimensions[1] = HEIGHT_CHANNEL;
    dimensions[2] = WIDTH_CHANNEL;

    Tensor<1,dim> center;
    center[0] = - (LENGTH_BFS_UP + GAP_CHANNEL_BFS + LENGTH_CHANNEL/2.0);
    center[1] = HEIGHT_CHANNEL/2.0;

    GridGenerator::subdivided_hyper_rectangle (triangulation,
                                               std::vector<unsigned int>({1,1,1}), //refinements
                                               Point<dim>(center-dimensions/2.0),
                                               Point<dim>(center+dimensions/2.0));
  }

  // TODO
//  // manifold
//  unsigned int manifold_id = 1;
//  for (typename Triangulation<dim>::cell_iterator cell = triangulation.begin(); cell != triangulation.end(); ++cell)
//  {
//    cell->set_all_manifold_ids(manifold_id);
//  }
//
//  // apply mesh stretching towards no-slip boundaries in y-direction
//  static const ManifoldTurbulentChannel<dim> manifold(dimensions);
//  triangulation.set_manifold(manifold_id, manifold);

  //periodicity in x- and z-direction (add 10 to avoid conflicts with dirichlet boundary, which is 0)
  triangulation.begin()->face(0)->set_all_boundary_ids(0+10);
  triangulation.begin()->face(1)->set_all_boundary_ids(1+10);
  //periodicity in z-direction
  if (dim == 3)
  {
    triangulation.begin()->face(4)->set_all_boundary_ids(2+10);
    triangulation.begin()->face(5)->set_all_boundary_ids(3+10);
  }

  GridTools::collect_periodic_faces(triangulation, 0+10, 1+10, 0, periodic_faces);
  if (dim == 3)
    GridTools::collect_periodic_faces(triangulation, 2+10, 3+10, 2, periodic_faces);

  triangulation.add_periodicity(periodic_faces);

  // perform global refinements: use one level finer for the channel
  triangulation.refine_global(n_refine_space);

  // fill boundary descriptor velocity
  // no slip boundaries at lower and upper wall with ID=0
  std::shared_ptr<Function<dim> > zero_function_velocity;
  zero_function_velocity.reset(new ZeroFunction<dim>(dim));
  boundary_descriptor_velocity->dirichlet_bc.insert(std::pair<types::boundary_id,std::shared_ptr<Function<dim> > >
                                                   (0,zero_function_velocity));

  // fill boundary descriptor pressure
  // no slip boundaries at lower and upper wall with ID=0
  std::shared_ptr<Function<dim> > pressure_bc_dudt;
  pressure_bc_dudt.reset(new PressureBC_dudt<dim>());
  boundary_descriptor_pressure->neumann_bc.insert(std::pair<types::boundary_id,std::shared_ptr<Function<dim> > >
                                                   (0,pressure_bc_dudt));
}

template<int dim>
void create_grid_and_set_boundary_conditions_2(
    parallel::distributed::Triangulation<dim>              &triangulation,
    unsigned int const                                     n_refine_space,
    std::shared_ptr<BoundaryDescriptorNavierStokesU<dim> > boundary_descriptor_velocity,
    std::shared_ptr<BoundaryDescriptorNavierStokesP<dim> > boundary_descriptor_pressure,
    std::vector<GridTools::PeriodicFacePair<typename
      Triangulation<dim>::cell_iterator> >                 &periodic_faces)
{
  /* --------------- Generate grid ------------------- */
  if(dim==2)
  {
    AssertThrow(false, ExcMessage("NotImplemented"));
  }
  else if(dim==3)
  {
    Triangulation<dim> tria_1, tria_2, tria_3;

    // inflow part of BFS
    GridGenerator::subdivided_hyper_rectangle(tria_1,
                                              std::vector<unsigned int>({1,1,1}),
                                              Point<dim>(-LENGTH_BFS_UP,0.0,-WIDTH_BFS/2.0),
                                              Point<dim>(0.0,HEIGHT_BFS_INFLOW,WIDTH_BFS/2.0));

    // downstream part of BFS (upper)
    GridGenerator::subdivided_hyper_rectangle(tria_2,
                                              std::vector<unsigned int>({10,1,1}),
                                              Point<dim>(0.0,0.0,-WIDTH_BFS/2.0),
                                              Point<dim>(LENGTH_BFS_DOWN,HEIGHT_BFS_INFLOW,WIDTH_BFS/2.0));

    // downstream part of BFS (lower = step)
    GridGenerator::subdivided_hyper_rectangle(tria_3,
                                              std::vector<unsigned int>({10,1,1}),
                                              Point<dim>(0.0,0.0,-WIDTH_BFS/2.0),
                                              Point<dim>(LENGTH_BFS_DOWN,-HEIGHT_BFS_STEP,WIDTH_BFS/2.0));

    Triangulation<dim> tmp1;
    GridGenerator::merge_triangulations (tria_1, tria_2, tmp1);
    GridGenerator::merge_triangulations (tmp1, tria_3, triangulation);
  }

  // set boundary ID's
  typename Triangulation<dim>::cell_iterator cell = triangulation.begin(), endc = triangulation.end();
  for(;cell!=endc;++cell)
  {
    for(unsigned int face_number=0; face_number < GeometryInfo<dim>::faces_per_cell; ++face_number)
    {
      // outflow boundary on the right has ID = 1
      if ((std::fabs(cell->face(face_number)->center()(0) - X1_COORDINATE_OUTFLOW)< 1e-12))
        cell->face(face_number)->set_boundary_id (1);
      // inflow boundary on the left has ID = 2
      if ((std::fabs(cell->face(face_number)->center()(0) - X1_COORDINATE_INFLOW)< 1e-12))
        cell->face(face_number)->set_boundary_id (2);

      // periodicity in z-direction (add 10 to avoid conflicts with other boundaries)
      if((std::fabs(cell->face(face_number)->center()(2) - WIDTH_BFS/2.0)< 1e-12))
        cell->face(face_number)->set_all_boundary_ids (2+10);
      // periodicity in z-direction (add 10 to avoid conflicts with other boundaries)
      if((std::fabs(cell->face(face_number)->center()(2) + WIDTH_BFS/2.0)< 1e-12))
        cell->face(face_number)->set_all_boundary_ids (3+10);
    }
  }

  // TODO
//  // manifold
//  unsigned int manifold_id = 1;
//  for (typename Triangulation<dim>::cell_iterator cell = triangulation.begin(); cell != triangulation.end(); ++cell)
//  {
//    cell->set_all_manifold_ids(manifold_id);
//  }

  // periodicity in z-direction
  GridTools::collect_periodic_faces(triangulation, 2+10, 3+10, 2, periodic_faces);
  triangulation.add_periodicity(periodic_faces);

  // perform global refinements
  triangulation.refine_global(n_refine_space);

  // fill boundary descriptor velocity
  // no slip boundaries at the upper and lower wall with ID=0
  std::shared_ptr<Function<dim> > zero_function_velocity;
  zero_function_velocity.reset(new ZeroFunction<dim>(dim));
  boundary_descriptor_velocity->dirichlet_bc.insert(std::pair<types::boundary_id,std::shared_ptr<Function<dim> > >
                                                   (0,zero_function_velocity));

  // inflow boundary condition at left boundary with ID=2: prescribe velocity profile which
  // is obtained as the results of the simulation on DOMAIN 1
  std::shared_ptr<Function<dim> > inflow_profile;
  inflow_profile.reset(new InflowProfile<dim>(dim));
  boundary_descriptor_velocity->dirichlet_bc.insert(std::pair<types::boundary_id,std::shared_ptr<Function<dim> > >
                                                   (2,inflow_profile));

  // outflow boundary condition at right boundary with ID=1
  boundary_descriptor_velocity->neumann_bc.insert(std::pair<types::boundary_id,std::shared_ptr<Function<dim> > >
                                                   (1,zero_function_velocity));

  // fill boundary descriptor pressure
  // no slip boundaries at the upper and lower wall with ID=0
  std::shared_ptr<Function<dim> > pressure_bc_dudt;
  pressure_bc_dudt.reset(new PressureBC_dudt<dim>());
  boundary_descriptor_pressure->neumann_bc.insert(std::pair<types::boundary_id,std::shared_ptr<Function<dim> > >
                                                   (0,pressure_bc_dudt));

  // inflow boundary condition at left boundary with ID=2
  // the inflow boundary condition is time dependent (du/dt != 0) but, for simplicity,
  // we assume that this is negligible when using the dual splitting scheme
  boundary_descriptor_pressure->neumann_bc.insert(std::pair<types::boundary_id,std::shared_ptr<Function<dim> > >
                                                   (2,pressure_bc_dudt));

  // outflow boundary condition at right boundary with ID=1: set pressure to zero
  std::shared_ptr<Function<dim> > zero_function_pressure;
  zero_function_pressure.reset(new ZeroFunction<dim>(1));
  boundary_descriptor_pressure->dirichlet_bc.insert(std::pair<types::boundary_id,std::shared_ptr<Function<dim> > >
                                                   (1,zero_function_pressure));
}


template<int dim>
void set_field_functions_1(std::shared_ptr<FieldFunctionsNavierStokes<dim> > field_functions)
{
  // initialize functions (analytical solution, rhs, boundary conditions)
  std::shared_ptr<Function<dim> > initial_solution_velocity;
  initial_solution_velocity.reset(new InitialSolutionVelocity<dim>());
  std::shared_ptr<Function<dim> > initial_solution_pressure;
  initial_solution_pressure.reset(new ZeroFunction<dim>(1));

  // use a constant body force for the turbulent channel (DOMAIN 1)
  std::shared_ptr<Function<dim> > right_hand_side;
  right_hand_side.reset(new RightHandSide<dim>());

  field_functions->initial_solution_velocity = initial_solution_velocity;
  field_functions->initial_solution_pressure = initial_solution_pressure;
  // This function will not be used since no analytical solution is available for this flow problem
  field_functions->analytical_solution_pressure = initial_solution_pressure;
  field_functions->right_hand_side = right_hand_side;
}

template<int dim>
void set_field_functions_2(std::shared_ptr<FieldFunctionsNavierStokes<dim> > field_functions)
{
  // initialize functions (analytical solution, rhs, boundary conditions)
  std::shared_ptr<Function<dim> > initial_solution_velocity;
  initial_solution_velocity.reset(new InitialSolutionVelocity<dim>());
//  initial_solution_velocity.reset(new ZeroFunction<dim>(dim));
  std::shared_ptr<Function<dim> > initial_solution_pressure;
  initial_solution_pressure.reset(new ZeroFunction<dim>(1));

  // no body forces for the second domain
  std::shared_ptr<Function<dim> > right_hand_side;
  right_hand_side.reset(new ZeroFunction<dim>(dim));

  field_functions->initial_solution_velocity = initial_solution_velocity;
  field_functions->initial_solution_pressure = initial_solution_pressure;
  // This function will not be used since no analytical solution is available for this flow problem
  field_functions->analytical_solution_pressure = initial_solution_pressure;
  field_functions->right_hand_side = right_hand_side;
}

template<int dim>
void set_analytical_solution(std::shared_ptr<AnalyticalSolutionNavierStokes<dim> > analytical_solution)
{
  analytical_solution->velocity.reset(new ZeroFunction<dim>(dim));
  analytical_solution->pressure.reset(new ZeroFunction<dim>(1));
}

// Postprocessor

#include "../../include/incompressible_navier_stokes/postprocessor/postprocessor.h"
#include "../../include/incompressible_navier_stokes/postprocessor/line_plot_calculation_statistics.h"
#include "../../include/incompressible_navier_stokes/postprocessor/mean_velocity_calculator.h"
#include "../../include/incompressible_navier_stokes/postprocessor/statistics_manager.h"
#include "../../include/incompressible_navier_stokes/postprocessor/inflow_data_calculator.h"

template<int dim>
struct PostProcessorDataBFS
{
  PostProcessorData<dim> pp_data;
  TurbulentChannelData turb_ch_data;
  InflowData<dim> inflow_data;
  BFSStatistics bfs_data;
  LinePlotData<dim> line_data;
  MeanVelocityCalculatorData<dim> mean_velocity_data;
};

template<int dim, int fe_degree_u, int fe_degree_p, typename Number>
class PostProcessorBFS : public PostProcessor<dim, fe_degree_u, fe_degree_p, Number>
{
public:
  PostProcessorBFS(PostProcessorDataBFS<dim> const & pp_data_bfs_in)
    :
    PostProcessor<dim,fe_degree_u,fe_degree_p, Number>(pp_data_bfs_in.pp_data),
    write_final_output(true),
    write_final_output_lines(true),
    write_final_output_mean_velocity(true),
    pp_data_bfs(pp_data_bfs_in)
  {
    inflow_data_calculator.reset(new InflowDataCalculator<dim,Number>(pp_data_bfs_in.inflow_data));
  }

  void setup(DoFHandler<dim> const                                  &dof_handler_velocity_in,
             DoFHandler<dim> const                                  &dof_handler_pressure_in,
             Mapping<dim> const                                     &mapping_in,
             MatrixFree<dim,Number> const                           &matrix_free_data_in,
             DofQuadIndexData const                                 &dof_quad_index_data_in,
             std::shared_ptr<AnalyticalSolutionNavierStokes<dim> >  analytical_solution_in)
  {
    // call setup function of base class
    PostProcessor<dim,fe_degree_u,fe_degree_p,Number>::setup(
        dof_handler_velocity_in,
        dof_handler_pressure_in,
        mapping_in,
        matrix_free_data_in,
        dof_quad_index_data_in,
        analytical_solution_in);

    // perform setup of turbulent channel related things
    statistics_turb_ch.reset(new StatisticsManager<dim>(dof_handler_velocity_in,mapping_in));

    bool individual_cells_are_stretched = false;

    // TODO
//    individual_cells_are_stretched = true;

    statistics_turb_ch->setup(&grid_transform_y,individual_cells_are_stretched);

    // inflow data
    if(pp_data_bfs.inflow_data.write_inflow_data == true)
    {
      inflow_data_calculator->setup(dof_handler_velocity_in,mapping_in);
    }

    if(pp_data_bfs.bfs_data.calculate_statistics == true)
    {
      statistics_bfs.reset(new LineStatisticsCalculator<dim>(dof_handler_velocity_in, dof_handler_pressure_in, mapping_in));
      statistics_bfs->setup(pp_data_bfs.line_data);
    }

    if(pp_data_bfs.mean_velocity_data.calculate_statistics == true)
    {
      centerline_velocity.reset(new MeanVelocityCalculator<dim,fe_degree_u,Number>(
          matrix_free_data_in, dof_quad_index_data_in,pp_data_bfs.mean_velocity_data));
    }
  }

  void do_postprocessing(parallel::distributed::Vector<Number> const &velocity,
                         parallel::distributed::Vector<Number> const &intermediate_velocity,
                         parallel::distributed::Vector<Number> const &pressure,
                         parallel::distributed::Vector<Number> const &vorticity,
                         std::vector<SolutionField<dim,Number> > const &additional_fields,
                         double const                                time,
                         int const                                   time_step_number)
  {
    Timer timer;
    timer.restart();

    PostProcessor<dim,fe_degree_u,fe_degree_p,Number>::do_postprocessing(
        velocity,
        intermediate_velocity,
        pressure,
        vorticity,
        additional_fields,
        time,
        time_step_number);

    // EPSILON: small number which is much smaller than the time step size
    const double EPSILON = 1.0e-10;
    if(pp_data_bfs.turb_ch_data.calculate_statistics == true)
    {
      if((time > pp_data_bfs.turb_ch_data.sample_start_time-EPSILON) &&
         (time < pp_data_bfs.turb_ch_data.sample_end_time+EPSILON) &&
         (time_step_number % pp_data_bfs.turb_ch_data.sample_every_timesteps == 0))
      {
        // evaluate statistics
        statistics_turb_ch->evaluate(velocity);

        // write intermediate output
        if(time_step_number % (pp_data_bfs.turb_ch_data.sample_every_timesteps * 100) == 0)
        {
          statistics_turb_ch->write_output(pp_data_bfs.turb_ch_data.filename_prefix,
                                           pp_data_bfs.turb_ch_data.viscosity);
        }
      }
      // write final output
      if((time > pp_data_bfs.turb_ch_data.sample_end_time-EPSILON) && write_final_output)
      {
        statistics_turb_ch->write_output(pp_data_bfs.turb_ch_data.filename_prefix,
                                         pp_data_bfs.turb_ch_data.viscosity);
        write_final_output = false;
      }
    }

    // inflow data
    if(pp_data_bfs.inflow_data.write_inflow_data == true)
    {
      inflow_data_calculator->calculate(velocity);
    }

   // calculate statistics for a set of lines
   if(pp_data_bfs.bfs_data.calculate_statistics == true)
   {
     if((time > pp_data_bfs.bfs_data.sample_start_time-EPSILON) &&
        (time < pp_data_bfs.bfs_data.sample_end_time+EPSILON) &&
        (time_step_number % pp_data_bfs.bfs_data.sample_every_timesteps == 0))
     {
       // evaluate statistics
       statistics_bfs->evaluate(velocity, pressure);

       // write intermediate output
       if(time_step_number % (pp_data_bfs.bfs_data.sample_every_timesteps) == 0)
       {
         statistics_bfs->write_output(pp_data_bfs.bfs_data.filename_prefix);
       }
     }
     // write final output
     if((time > pp_data_bfs.bfs_data.sample_end_time-EPSILON) && write_final_output_lines)
     {
       statistics_bfs->write_output(pp_data_bfs.bfs_data.filename_prefix);
       write_final_output_lines = false;
     }
   }

   // calculate mean centerline velocity
   if(pp_data_bfs.mean_velocity_data.calculate_statistics == true)
   {
     if((time > pp_data_bfs.mean_velocity_data.sample_start_time-EPSILON) &&
        (time < pp_data_bfs.mean_velocity_data.sample_end_time+EPSILON) &&
        (time_step_number % pp_data_bfs.mean_velocity_data.sample_every_timesteps == 0))
     {
       // evaluate statistics
       centerline_velocity->evaluate(velocity);

       // write intermediate output
       if(time_step_number % (pp_data_bfs.mean_velocity_data.sample_every_timesteps) == 0)
       {
       centerline_velocity->write_output(pp_data_bfs.mean_velocity_data.filename_prefix);
       }
     }
     // write final output
     if((time > pp_data_bfs.mean_velocity_data.sample_end_time-EPSILON) && write_final_output_mean_velocity)
     {
       centerline_velocity->write_output(pp_data_bfs.mean_velocity_data.filename_prefix);
       write_final_output_mean_velocity = false;
     }
   }
  }

  bool write_final_output;
  bool write_final_output_lines;
  bool write_final_output_mean_velocity;
  PostProcessorDataBFS<dim> pp_data_bfs;
  std::shared_ptr<StatisticsManager<dim> > statistics_turb_ch;
  std::shared_ptr<InflowDataCalculator<dim, Number> > inflow_data_calculator;
  std::shared_ptr<LineStatisticsCalculator<dim> > statistics_bfs;
  std::shared_ptr<MeanVelocityCalculator<dim,fe_degree_u,Number> > centerline_velocity;
};

#include "../../include/incompressible_navier_stokes/postprocessor/postprocessor.h"

template<int dim, typename Number>
std::shared_ptr<PostProcessorBase<dim,Number> >
construct_postprocessor(InputParametersNavierStokes<dim> const &param)
{
  PostProcessorData<dim> pp_data;
  pp_data.output_data = param.output_data;
  pp_data.error_data = param.error_data;
  pp_data.lift_and_drag_data = param.lift_and_drag_data;
  pp_data.pressure_difference_data = param.pressure_difference_data;
  pp_data.mass_data = param.mass_data;

  PostProcessorDataBFS<dim> pp_data_bfs;
  pp_data_bfs.pp_data = pp_data;
  pp_data_bfs.turb_ch_data = param.turb_ch_data;
  pp_data_bfs.inflow_data = param.inflow_data;
  pp_data_bfs.bfs_data = param.bfs_statistics;
  pp_data_bfs.line_data = param.line_plot_data;
  pp_data_bfs.mean_velocity_data = param.mean_velocity_data;

  std::shared_ptr<PostProcessorBase<dim,Number> > pp;
  pp.reset(new PostProcessorBFS<dim,FE_DEGREE_VELOCITY,FE_DEGREE_PRESSURE,Number>(pp_data_bfs));

  return pp;
}


#endif /* APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_TURBULENT_CHANNEL_H_ */