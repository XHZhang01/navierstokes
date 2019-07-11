/*
*3D_tube_precursor_poiseuille.h
*
* Created on: 02.Mai.2019
* Author: Xuhui Zhang
*/

#ifndef APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_3D_tube_precursor_poiseuille_H_
#define APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_3D_tube_precursor_poiseuille_H_

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

// set the polynomial degree of the shape functions for velocity and pressure
unsigned int const FE_DEGREE_VELOCITY = 2;
unsigned int const FE_DEGREE_PRESSURE = FE_DEGREE_VELOCITY-1;

// set the number of refine levels for spatial convergence tests
unsigned int const REFINE_STEPS_SPACE_MIN = 1;
unsigned int const REFINE_STEPS_SPACE_MAX = REFINE_STEPS_SPACE_MIN;

// set the number of refine levels for temporal convergence tests
unsigned int const REFINE_STEPS_TIME_MIN = 0;
unsigned int const REFINE_STEPS_TIME_MAX = 0; //REFINE_STEPS_TIME_MIN;

// set problem specific parameters like physical dimensions, etc.
ProblemType PROBLEM_TYPE = ProblemType::Unsteady;
//const unsigned int TEST_CASE = 3; // 1, 2 or 3
//const double Um = (DIMENSION == 2 ? (TEST_CASE==1 ? 0.3 : 1.5) : (TEST_CASE==1 ? 0.45 : 2.25));

// output folders
std::string OUTPUT_FOLDER = "output/3dtpp/Re3500/";
std::string OUTPUT_FOLDER_VTU = OUTPUT_FOLDER + "vtu/";
std::string OUTPUT_NAME_1 = "precursor";
std::string FILENAME_FLOWRATE = "precursor_mean_velocity";

// radius
double const R = 0.002;
double const R_INNER = R;
double const R_OUTER = 3.0*R;
double const D = 2.0*R_OUTER;

// lengths (dimensions in flow direction z)
double const LENGTH_PRECURSOR = 8.0*R_OUTER;
double const LENGTH_INFLOW = 8.0*R_OUTER;
double const LENGTH_CONE = (R_OUTER-R_INNER)/std::tan(20.0/2.0*numbers::PI/180.0);
double const LENGTH_THROAT = 0.04;
double const LENGTH_OUTFLOW = 20.0*R_OUTER;
double const OFFSET = 2.0*R_OUTER;

// mesh parameters
unsigned int const N_CELLS_AXIAL = 2;
unsigned int const N_CELLS_AXIAL_PRECURSOR = 4*N_CELLS_AXIAL;
//unsigned int const N_CELLS_AXIAL_INFLOW = 4*N_CELLS_AXIAL;
//unsigned int const N_CELLS_AXIAL_CONE = 2*N_CELLS_AXIAL;
//unsigned int const N_CELLS_AXIAL_THROAT = 4*N_CELLS_AXIAL;
//unsigned int const N_CELLS_AXIAL_OUTFLOW = 10*N_CELLS_AXIAL;

//unsigned int const MANIFOLD_ID_CYLINDER = 1234;
//unsigned int const MANIFOLD_ID_OFFSET_CONE = 7890;

// z-coordinates
//double const Z2_OUTFLOW = LENGTH_OUTFLOW;
//double const Z1_OUTFLOW = 0.0;

//double const Z2_THROAT = 0.0;
//double const Z1_THROAT = - LENGTH_THROAT;

//double const Z2_CONE = - LENGTH_THROAT;
//double const Z1_CONE = - LENGTH_THROAT - LENGTH_CONE;

//double const Z2_INFLOW = - LENGTH_THROAT - LENGTH_CONE;
//double const Z1_INFLOW = - LENGTH_THROAT - LENGTH_CONE - LENGTH_INFLOW;

double const Z2_PRECURSOR = - LENGTH_THROAT - LENGTH_CONE - LENGTH_INFLOW - OFFSET;
double const Z1_PRECURSOR = - LENGTH_THROAT - LENGTH_CONE - LENGTH_INFLOW - OFFSET - LENGTH_PRECURSOR;

double const AREA_INFLOW = R_OUTER*R_OUTER*numbers::PI;
double const AREA_THROAT = R_INNER*R_INNER*numbers::PI;

// kinematic viscosity (same viscosity for all Reynolds numbers)
double const VISCOSITY = 3.31e-6;

// set the throat Reynolds number Re_throat = U_{mean,throat} * (2 R_throat) / nu
double const RE = 3500; //500; //2000; //3500; //5000; //6500; //8000;

double const MEAN_VELOCITY_THROAT = RE * VISCOSITY / (2.0*R_INNER);
double const TARGET_FLOW_RATE = MEAN_VELOCITY_THROAT*AREA_THROAT;
double const MEAN_VELOCITY_INFLOW = TARGET_FLOW_RATE/AREA_INFLOW;

double const MAX_VELOCITY = 2.0*TARGET_FLOW_RATE/AREA_INFLOW;
double const MAX_VELOCITY_CFL = 2.0*TARGET_FLOW_RATE/AREA_THROAT;

// start and end time

// estimation of flow-through time T_0 (through nozzle section)
// based on the mean velocity through throat
double const T_0 = LENGTH_THROAT/MEAN_VELOCITY_THROAT;
double const START_TIME_PRECURSOR = -500.0*T_0; // let the flow develop
//double const START_TIME_NOZZLE = 0.0*T_0;
double const END_TIME = 0.0*T_0; //150.0*T_0;

// output
bool const WRITE_OUTPUT = true;
double const OUTPUT_START_TIME_PRECURSOR = START_TIME_PRECURSOR;
//double const OUTPUT_START_TIME_NOZZLE = START_TIME_NOZZLE;
double const OUTPUT_INTERVAL_TIME = 5.0*T_0;  //10.0*T_0;

// sampling

// sampling interval should last over (100-200) * T_0 according to preliminary results.
//double const SAMPLE_START_TIME = 50.0*T_0; // let the flow develop
//double const SAMPLE_END_TIME = END_TIME; // that's the only reasonable choice
//unsigned int SAMPLE_EVERY_TIMESTEPS = 1;
//unsigned int WRITE_OUTPUT_EVERY_TIMESTEPS = SAMPLE_EVERY_TIMESTEPS*100;

// data structures that we need to control the mass flow rate:
// NOTA BENE: these variables will be modified by the postprocessor!
double FLOW_RATE = 0.0;
// the flow rate controller also needs the time step size as parameter
double TIME_STEP_FLOW_RATE_CONTROLLER = 1.0;

class FlowRateController
{
public:
  FlowRateController()
    :
    // initialize the body force such that the desired flow rate is obtained
    // under the assumption of a parabolic velocity profile in radial direction
    f(4.0*VISCOSITY*MAX_VELOCITY/std::pow(R_OUTER,2.0)) // f(t=t_0) = f_0
  {}

  double get_body_force()
  {
    return f;
  }

  void update_body_force()
  {
    // use an I-controller to asymptotically reach the desired target flow rate

    // dimensional analysis: [k] = 1/(m^2 s^2) -> k = const * U_{mean,inflow}^2 / D^4
    // constant: choose a default value of 1
    double const k = 1.0e0*std::pow(MEAN_VELOCITY_INFLOW,2.0)/std::pow(D,4.0);
    f += k*(TARGET_FLOW_RATE - FLOW_RATE)*TIME_STEP_FLOW_RATE_CONTROLLER;
  }

private:
  double f;
};

// use a global variable which will be called by the postprocessor
// in order to update the body force.
FlowRateController FLOW_RATE_CONTROLLER;

// - choose a large number of points to ensure a smooth inflow profile
unsigned int N_POINTS_R = 10 * (FE_DEGREE_VELOCITY+1) * std::pow(2.0, 1);
unsigned int N_POINTS_PHI = N_POINTS_R;
std::vector<double> R_VALUES(N_POINTS_R);
std::vector<double> PHI_VALUES(N_POINTS_PHI);
std::vector<Tensor<1,DIMENSION,double> > VELOCITY_VALUES(N_POINTS_R*N_POINTS_PHI);

// initialize vectors
void initialize_r_and_phi_values()
{
  AssertThrow(N_POINTS_R >= 2, ExcMessage("Variable N_POINTS_R is invalid"));
  AssertThrow(N_POINTS_PHI >= 2, ExcMessage("Variable N_POINTS_PHI is invalid"));

  // 0 <= radius <= R_OUTER
  for(unsigned int i=0; i<N_POINTS_R; ++i)
    R_VALUES[i] = double(i)/double(N_POINTS_R-1)*R_OUTER;

  // - pi <= phi <= pi
  for(unsigned int i=0; i<N_POINTS_PHI; ++i)
    PHI_VALUES[i] = -numbers::PI + double(i)/double(N_POINTS_PHI-1)*2.0*numbers::PI;
}

void initialize_velocity_values()
{
  AssertThrow(N_POINTS_R >= 2, ExcMessage("Variable N_POINTS_R is invalid"));
  AssertThrow(N_POINTS_PHI >= 2, ExcMessage("Variable N_POINTS_PHI is invalid"));

  for(unsigned int iy=0; iy<N_POINTS_R; ++iy)
  {
    for(unsigned int iz=0; iz<N_POINTS_PHI; ++iz)
    {
      Tensor<1,DIMENSION,double> velocity;
      // flow in z-direction
      velocity[2] = MAX_VELOCITY*(1.0-std::pow(R_VALUES[iy]/R_OUTER,2.0));

      VELOCITY_VALUES[iy*N_POINTS_PHI + iz] = velocity;
    }
  }
}
/*
 *  This function returns the radius of the cross-section at a
 *  specified location z in streamwise direction.
 */
/*
double radius_function(double const z)
{
  double radius = R_OUTER;

  if(z >= Z1_INFLOW && z <= Z2_INFLOW)
    radius = R_OUTER;
  else if(z >= Z1_CONE && z <= Z2_CONE)
    radius = R_OUTER * (1.0 - (z-Z1_CONE)/(Z2_CONE-Z1_CONE)*(R_OUTER-R_INNER)/R_OUTER);
  else if(z >= Z1_THROAT && z <= Z2_THROAT)
    radius = R_INNER;
  else if(z > Z1_OUTFLOW && z <= Z2_OUTFLOW)
    radius = R_OUTER;

  return radius;
}
*/
template<int dim>
void InputParameters<dim>::set_input_parameters()
{
  // MATHEMATICAL MODEL
  problem_type = ProblemType::Unsteady;
  equation_type = EquationType::NavierStokes;
  use_outflow_bc_convective_term = true;
  formulation_viscous_term = FormulationViscousTerm::LaplaceFormulation;
  formulation_convective_term = FormulationConvectiveTerm::DivergenceFormulation;
  right_hand_side = false;

  // PHYSICAL QUANTITIES
  start_time = START_TIME_PRECURSOR;
  end_time = END_TIME;
  viscosity = VISCOSITY;

  // TEMPORAL DISCRETIZATION
    solver_type = SolverType::Unsteady;

  temporal_discretization = TemporalDiscretization::BDFDualSplittingScheme;
  treatment_of_convective_term = TreatmentOfConvectiveTerm::Explicit;
  calculation_of_time_step_size = TimeStepCalculation::CFL;
  adaptive_time_stepping = true;
  //temporal_discretization = TemporalDiscretization::BDFPressureCorrection;
  //treatment_of_convective_term = TreatmentOfConvectiveTerm::Implicit;
  //calculation_of_time_step_size = TimeStepCalculation::CFL;
  adaptive_time_stepping_limiting_factor = 3.0;
  max_velocity = MAX_VELOCITY_CFL;
  cfl = 0.3;//0.4
  cfl_exponent_fe_degree_velocity = 1.5;
  time_step_size = 1.0e-1;
  order_time_integrator = 2;
  start_with_low_order = true;

  // SPATIAL DISCRETIZATION

  // triangulation
  triangulation_type = TriangulationType::Distributed;

  // mapping
  degree_mapping = FE_DEGREE_VELOCITY;

  // convective term
  if(formulation_convective_term == FormulationConvectiveTerm::DivergenceFormulation)
      upwind_factor = 0.5; // allows using larger CFL values for explicit formulations

  // viscous term
  IP_formulation_viscous = InteriorPenaltyFormulation::SIPG;
  IP_factor_viscous = 1.0;

  // special case: pure DBC's
  //if(domain_id == 1)
    pure_dirichlet_bc = false;
  //else if(domain_id == 2)
    //pure_dirichlet_bc = false;

  // div-div and continuity penalty terms
  use_divergence_penalty = true;
  divergence_penalty_factor = 1.0e0;
  use_continuity_penalty = true;
  continuity_penalty_factor = divergence_penalty_factor;
  add_penalty_terms_to_monolithic_system = false;

  // TURBULENCE
  use_turbulence_model = false;
  turbulence_model = TurbulenceEddyViscosityModel::Sigma;
  // Smagorinsky: 0.165, Vreman: 0.28, WALE: 0.50, Sigma: 1.35
  turbulence_model_constant = 1.35;

  // PROJECTION METHODS

  // pressure Poisson equation
  IP_factor_pressure = 1.0;
  solver_data_pressure_poisson = SolverData(1000,1.e-12,1.e-3,100);
  solver_pressure_poisson = SolverPressurePoisson::CG; //FGMRES;
  preconditioner_pressure_poisson = PreconditionerPressurePoisson::Multigrid;
  multigrid_data_pressure_poisson.type = MultigridType::phMG;
  multigrid_data_pressure_poisson.p_sequence = PSequenceType::Bisect;
  multigrid_data_pressure_poisson.dg_to_cg_transfer = DG_To_CG_Transfer::Fine;

  //Variant 1:
  //multigrid_data_pressure_poisson.smoother_data.smoother = MultigridSmoother::Chebyshev;
  //multigrid_data_pressure_poisson.coarse_problem.preconditioner = MultigridCoarseGridPreconditioner::PointJacobi;

  //Variant 2 (Standardfall):
  multigrid_data_pressure_poisson.coarse_problem.solver = MultigridCoarseGridSolver::CG;
  multigrid_data_pressure_poisson.coarse_problem.preconditioner = MultigridCoarseGridPreconditioner::AMG;

  // projection step
  solver_projection = SolverProjection::CG;
  solver_data_projection = SolverData(1000, 1.e-12, 1.e-3);
  preconditioner_projection = PreconditionerProjection::InverseMassMatrix;
  update_preconditioner_projection = true;

  // HIGH-ORDER DUAL SPLITTING SCHEME

  // formulations
  order_extrapolation_pressure_nbc = order_time_integrator <=2 ? order_time_integrator : 2;

  // viscous step
  solver_viscous = SolverViscous::CG;
  solver_data_viscous = SolverData(1000,1.e-12,1.e-3);
  preconditioner_viscous = PreconditionerViscous::InverseMassMatrix;

  // PRESSURE-CORRECTION SCHEME

  // formulation
  order_pressure_extrapolation = 1; // use 0 for non-incremental formulation
  rotational_formulation = true; // use false for standard formulation

  // momentum step

  // Newton solver
  newton_solver_data_momentum = NewtonSolverData(100,1.e-12,1.e-3);

  // linear solver
  if(treatment_of_convective_term == TreatmentOfConvectiveTerm::Implicit)
    solver_data_momentum = SolverData(1e4, 1.e-12, 1.e-1, 100);
  else
    solver_data_momentum = SolverData(1e4, 1.e-12, 1.e-3, 100);

  solver_momentum = SolverMomentum::GMRES;
  preconditioner_momentum = MomentumPreconditioner::InverseMassMatrix;
  update_preconditioner_momentum = false;

  // COUPLED NAVIER-STOKES SOLVER
  use_scaling_continuity = false;

  // nonlinear solver (Newton solver)
  newton_solver_data_coupled = NewtonSolverData(100,1.e-20,1.e-3);

  // linear solver
  solver_coupled = SolverCoupled::GMRES; //GMRES; //FGMRES;
  if(treatment_of_convective_term == TreatmentOfConvectiveTerm::Implicit)
    solver_data_coupled = SolverData(1e4, 1.e-12, 1.e-1, 100);
  else
    solver_data_coupled = SolverData(1e4, 1.e-12, 1.e-3, 100);

  // preconditioning linear solver
  preconditioner_coupled = PreconditionerCoupled::BlockTriangular;
  update_preconditioner_coupled = false;

  // preconditioner velocity/momentum block
  preconditioner_velocity_block = MomentumPreconditioner::InverseMassMatrix;

  // preconditioner Schur-complement block
  preconditioner_pressure_block = SchurComplementPreconditioner::CahouetChabard; //PressureConvectionDiffusion;
  discretization_of_laplacian =  DiscretizationOfLaplacian::Classical;

  // Chebyshev moother
  multigrid_data_pressure_block.smoother_data.smoother = MultigridSmoother::Chebyshev;
  multigrid_data_pressure_block.coarse_problem.solver = MultigridCoarseGridSolver::Chebyshev;

  // OUTPUT AND POSTPROCESSING

  // output of solver information
  solver_info_data.print_to_screen = true;
  solver_info_data.interval_time = T_0;

  // write output for visualization of results
  output_data.write_output = WRITE_OUTPUT;
  output_data.output_folder = OUTPUT_FOLDER_VTU;
  output_data.output_name = OUTPUT_NAME_1;
  output_data.output_start_time = OUTPUT_START_TIME_PRECURSOR;
  output_data.output_interval_time = OUTPUT_INTERVAL_TIME;
  output_data.write_divergence = true;
  output_data.write_processor_id = true;
  output_data.mean_velocity.calculate = true;
  //output_data.mean_velocity.sample_start_time = SAMPLE_START_TIME;
  //output_data.mean_velocity.sample_end_time = SAMPLE_END_TIME;
  //output_data.mean_velocity.sample_every_timesteps = 1;
  output_data.degree = FE_DEGREE_VELOCITY;

  // calculation of flow rate (use volume-based computation)
  mean_velocity_data.calculate = true;
  mean_velocity_data.filename_prefix = OUTPUT_FOLDER + FILENAME_FLOWRATE;
  Tensor<1,dim,double> direction; direction[2] = 1.0;
  mean_velocity_data.direction = direction;
  mean_velocity_data.write_to_file = true;
}


/**************************************************************************************/
/*                                                                                    */
/*                        GENERATE GRID AND SET BOUNDARY INDICATORS                   */
/*                                                                                    */
/**************************************************************************************/

#include "../../include/functionalities/one_sided_cylindrical_manifold.h"

template<int dim>
void create_grid_and_set_boundary_ids(
    std::shared_ptr<parallel::Triangulation<dim>>     triangulation,
    unsigned int const                                n_refine_space,
    std::vector<GridTools::PeriodicFacePair<typename
      Triangulation<dim>::cell_iterator> >            &/*periodic_faces*/)
{
  /*
   *   PRECURSOR
   */
  Triangulation<2> tria_2d;
  GridGenerator::hyper_ball(tria_2d, Point<2>(), R_OUTER);
  GridGenerator::extrude_triangulation(tria_2d,N_CELLS_AXIAL_PRECURSOR+1,LENGTH_PRECURSOR,*triangulation);
  Tensor<1,dim> offset = Tensor<1,dim>();
  offset[2] = Z1_PRECURSOR;
  GridTools::shift(offset,*triangulation);

  /*
   *  MANIFOLDS
   */
  triangulation->set_all_manifold_ids(0);

  // first fill vectors of manifold_ids and face_ids
  std::vector<unsigned int> manifold_ids;
  std::vector<unsigned int> face_ids;

  for (typename Triangulation<dim>::cell_iterator cell = triangulation->begin();cell != triangulation->end(); ++cell)
  {
    for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    {
      bool face_at_sphere_boundary = true;
      for (unsigned int v=0; v<GeometryInfo<dim-1>::vertices_per_cell; ++v)
      {
        Point<dim> point = Point<dim>(0,0,cell->face(f)->vertex(v)[2]);

        if (std::abs((cell->face(f)->vertex(v)-point).norm()-R_OUTER) > 1e-12)
          face_at_sphere_boundary = false;
      }
      if (face_at_sphere_boundary)
      {
        face_ids.push_back(f);
        unsigned int manifold_id = manifold_ids.size() + 1;
        cell->set_all_manifold_ids(manifold_id);
        manifold_ids.push_back(manifold_id);
      }
    }
  }

  // generate vector of manifolds and apply manifold to all cells that have been marked
  static std::vector<std::shared_ptr<Manifold<dim> > > manifold_vec;
  manifold_vec.resize(manifold_ids.size());

  for(unsigned int i=0;i<manifold_ids.size();++i)
  {
    for (typename Triangulation<dim>::cell_iterator cell = triangulation->begin(); cell != triangulation->end(); ++cell)
    {
      if(cell->manifold_id() == manifold_ids[i])
      {
        manifold_vec[i] = std::shared_ptr<Manifold<dim> >(
            static_cast<Manifold<dim>*>(new OneSidedCylindricalManifold<dim>(cell,face_ids[i],Point<dim>())));
        triangulation->set_manifold(manifold_ids[i],*(manifold_vec[i]));
      }
    }
  }


  std::cout<< "\nError 2 " << std::endl;

  /*
   *  BOUNDARY ID's
   */
  typename Triangulation<dim>::cell_iterator cell = triangulation->begin(), endc = triangulation->end();
  for(;cell!=endc;++cell)
  {
    for(unsigned int face_number=0; face_number < GeometryInfo<dim>::faces_per_cell; ++face_number)
    {
      // left boundary
      if ((std::fabs(cell->face(face_number)->center()[2] - Z1_PRECURSOR) < 1e-12))
      {
        cell->face(face_number)->set_boundary_id (1);
      }

      // right boundary
      if ((std::fabs(cell->face(face_number)->center()[2] - Z2_PRECURSOR) < 1e-12))
      {
        cell->face(face_number)->set_boundary_id (2);
      }
    }
  }
  std::cout<< "\nn_refine_space = " << n_refine_space << std::endl;

  // perform global refinements
  triangulation->refine_global(n_refine_space);

  std::cout<< "\nError 4 " << std::endl;
}


/**************************************************************************************/
/*                                                                                    */
/*    FUNCTIONS (ANALYTICAL SOLUTION, BOUNDARY CONDITIONS, VELOCITY FIELD, etc.)      */
/*                                                                                    */
/**************************************************************************************/


namespace IncNS
{

template<int dim>
class PressureInflowBC : public Function<dim>
{
public:
  PressureInflowBC (const double time = 0.)
    :
    Function<dim>(1 /*n_components*/, time)
  {}

  double value (const Point<dim>   &p,
                const unsigned int /*component = 0*/) const
  {
    double t = this->get_time();

    // assuming a parabolic velocity profile
    double const pressure_gradient = -2.*VISCOSITY*MAX_VELOCITY;
    double const pi = numbers::PI;
    double const T = 5;
    std::cout << "\nPressure: " << (p[0]-Z2_PRECURSOR) * pressure_gradient * std::sin(pi*t/T);

    return (p[0]-Z2_PRECURSOR) * pressure_gradient * std::sin(pi*t/T);
  }
};

template<int dim>
class RightHandSide : public Function<dim>
{
public:
 RightHandSide (const double time = 0.)
   :
   Function<dim>(dim, time)
 {}

 double value (const Point<dim>    & /*p*/,
               const unsigned int  component = 0) const
 {
   double result = 0.0;

   // Channel flow with periodic bc in z-direction:
   // The flow is driven by body force in z-direction
   if(component==2)
   {
     result = FLOW_RATE_CONTROLLER.get_body_force();
   }

   return result;
 }
};

template<int dim>
void set_boundary_conditions(
    std::shared_ptr<BoundaryDescriptorU<dim> > boundary_descriptor_velocity,
    std::shared_ptr<BoundaryDescriptorP<dim> > boundary_descriptor_pressure)
{
  typedef typename std::pair<types::boundary_id,std::shared_ptr<Function<dim> > > pair;

  // fill boundary descriptor velocity
  boundary_descriptor_velocity->dirichlet_bc.insert(pair(0,new Functions::ZeroFunction<dim>(dim)));
  boundary_descriptor_velocity->neumann_bc.insert(pair(1,new Functions::ZeroFunction<dim>(dim)));
  boundary_descriptor_velocity->neumann_bc.insert(pair(2,new Functions::ZeroFunction<dim>(dim)));

  // fill boundary descriptor pressure
  boundary_descriptor_pressure->neumann_bc.insert(pair(0,new Functions::ZeroFunction<dim>(1)));
  boundary_descriptor_pressure->dirichlet_bc.insert(pair(1,new PressureInflowBC<dim>()));
  boundary_descriptor_pressure->dirichlet_bc.insert(pair(2,new Functions::ZeroFunction<dim>(1)));
}

template<int dim>
void set_field_functions(std::shared_ptr<FieldFunctions<dim> > field_functions)
{
  field_functions->initial_solution_velocity.reset(new Functions::ZeroFunction<dim>(dim));
  field_functions->initial_solution_pressure.reset(new Functions::ZeroFunction<dim>(1));
  field_functions->analytical_solution_pressure.reset(new Functions::ZeroFunction<dim>(1));
  field_functions->right_hand_side.reset(new Functions::ZeroFunction<dim>(dim));
}

template<int dim>
void set_analytical_solution(std::shared_ptr<AnalyticalSolution<dim> > analytical_solution)
{
  analytical_solution->velocity.reset(new Functions::ZeroFunction<dim>(dim));
  analytical_solution->pressure.reset(new Functions::ZeroFunction<dim>(1));
}

#include "../../include/incompressible_navier_stokes/postprocessor/postprocessor.h"

template<int dim, int degree_u, int degree_p, typename Number>
std::shared_ptr<PostProcessorBase<dim, degree_u, degree_p, Number> >
construct_postprocessor(InputParameters<dim> const &param)
{
  PostProcessorData<dim> pp_data;

  pp_data.output_data = param.output_data;
  pp_data.error_data = param.error_data;
  pp_data.lift_and_drag_data = param.lift_and_drag_data;
  pp_data.pressure_difference_data = param.pressure_difference_data;
  pp_data.mass_data = param.mass_data;
  pp_data.line_plot_data = param.line_plot_data;

  std::shared_ptr<PostProcessor<dim,degree_u,degree_p,Number> > pp;
  pp.reset(new PostProcessor<dim,degree_u,degree_p,Number>(pp_data));

  return pp;
}

}

#endif /* APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_3D_tube_precursor_poiseuille_H_ */
