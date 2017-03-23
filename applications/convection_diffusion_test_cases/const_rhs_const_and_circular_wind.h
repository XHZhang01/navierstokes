/*
 * ConstWindConstRHS.h
 *
 *  Created on: Nov 23, 2016
 *      Author: fehn
 */

#ifndef APPLICATIONS_CONVECTION_DIFFUSION_TEST_CASES_CONST_RHS_CONST_AND_CIRCULAR_WIND_H_
#define APPLICATIONS_CONVECTION_DIFFUSION_TEST_CASES_CONST_RHS_CONST_AND_CIRCULAR_WIND_H_

#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>


/**************************************************************************************/
/*                                                                                    */
/*                                 INPUT PARAMETERS                                   */
/*                                                                                    */
/**************************************************************************************/

// constant source term inside rectangular domain
// pure dirichlet boundary conditions (homogeneous)
// use constant or circular advection velocity

// set the number of space dimensions: DIMENSION = 2, 3
const unsigned int DIMENSION = 2;

// set the polynomial degree of the shape functions
const unsigned int FE_DEGREE = 8;

// set the number of refine levels for spatial convergence tests
const unsigned int REFINE_STEPS_SPACE_MIN = 1;
const unsigned int REFINE_STEPS_SPACE_MAX = 8;

// set the number of refine levels for temporal convergence tests
const unsigned int REFINE_STEPS_TIME_MIN = 0;
const unsigned int REFINE_STEPS_TIME_MAX = 0;

// problem specific parameters
const double START_TIME = 0.0;
const double DIFFUSIVITY = 1.0e0;

enum class TypeVelocityField { Constant, Circular, CircularZeroAtBoundary };
TypeVelocityField const TYPE_VELOCITY_FIELD = TypeVelocityField::Circular; //CircularZeroAtBoundary; //Circular; //Constant;

void InputParametersConvDiff::set_input_parameters()
{
  // MATHEMATICAL MODEL
  problem_type = ProblemType::Steady;
  equation_type = EquationType::ConvectionDiffusion;
  right_hand_side = true;

  // PHYSICAL QUANTITIES
  start_time = START_TIME;
  end_time = 1.0;
  diffusivity = DIFFUSIVITY;

  // TEMPORAL DISCRETIZATION
  temporal_discretization = TemporalDiscretization::BDF;
  treatment_of_convective_term = TreatmentOfConvectiveTerm::Implicit;
  order_time_integrator = 2;
  start_with_low_order = true;
  calculation_of_time_step_size = TimeStepCalculation::ConstTimeStepUserSpecified;
  time_step_size = 1.0e-1;
  cfl_number = 0.2;
  diffusion_number = 0.01;

  // SPATIAL DISCRETIZATION
  // convective term
  numerical_flux_convective_operator = NumericalFluxConvectiveOperator::LaxFriedrichsFlux;

  // viscous term
  IP_factor = 1.0;

  // SOLVER
  solver = Solver::GMRES;
  abs_tol = 1.e-20;
  rel_tol = 1.e-8;
  max_iter = 1e3;
  max_n_tmp_vectors = 100;
  preconditioner = Preconditioner::MultigridConvectionDiffusion; //PointJacobi; //BlockJacobi; //MultigridDiffusion; //MultigridConvectionDiffusion;
  // MG smoother
  multigrid_data.smoother = MultigridSmoother::GMRES; //Jacobi; //GMRES; //Chebyshev; //ChebyshevNonsymmetricOperator;

  // MG smoother data: GMRES smoother
  multigrid_data.gmres_smoother_data.preconditioner = PreconditionerGMRESSmoother::BlockJacobi; //None; //PointJacobi; //BlockJacobi;
  multigrid_data.gmres_smoother_data.number_of_iterations = 10;

  // MG smoother data: Jacobi smoother
  multigrid_data.jacobi_smoother_data.preconditioner = PreconditionerJacobiSmoother::BlockJacobi; //None; //PointJacobi; //BlockJacobi;
  multigrid_data.jacobi_smoother_data.number_of_smoothing_steps = 5;
  multigrid_data.jacobi_smoother_data.damping_factor = 0.8;

  // MG coarse grid solver
  multigrid_data.coarse_solver = MultigridCoarseGridSolver::GMRES_NoPreconditioner; //Chebyshev; //GMRES_Jacobi;

  update_preconditioner = false;

  // NUMERICAL PARAMETERS
  runtime_optimization = false;

  // OUTPUT AND POSTPROCESSING
  print_input_parameters = false;
  output_data.write_output = false; //true;
  output_data.output_prefix = "const_wind_const_rhs";
  output_data.output_start_time = start_time;
  output_data.output_interval_time = (end_time-start_time);// /20;
  output_data.number_of_patches = FE_DEGREE;

  error_data.analytical_solution_available = false;
  error_data.error_calc_start_time = start_time;
  error_data.error_calc_interval_time = output_data.output_interval_time;

  output_solver_info_every_timesteps = 1e2;
}


/**************************************************************************************/
/*                                                                                    */
/*    FUNCTIONS (ANALYTICAL SOLUTION, BOUNDARY CONDITIONS, VELOCITY FIELD, etc.)      */
/*                                                                                    */
/**************************************************************************************/

/*
 *  Analytical solution
 */

template<int dim>
class AnalyticalSolution : public Function<dim>
{
public:
  AnalyticalSolution (const unsigned int  n_components = 1,
                      const double        time = 0.)
    :
    Function<dim>(n_components, time)
  {}

  virtual ~AnalyticalSolution(){};

  virtual double value (const Point<dim>   &p,
                        const unsigned int component = 0) const;
};

template<int dim>
double AnalyticalSolution<dim>::value(const Point<dim>    &/* p */,
                                      const unsigned int  /* component */) const
{
  double t = this->get_time();
  double result = 0.0;

  return result;
}

/*
 *  Right-hand side
 */

template<int dim>
class RightHandSide : public Function<dim>
{
public:
  RightHandSide (const unsigned int   n_components = 1,
                 const double         time = 0.)
    :
    Function<dim>(n_components, time)
  {}

  virtual ~RightHandSide(){};

  virtual double value (const Point<dim>    &p,
                       const unsigned int  component = 0) const;
};

template<int dim>
double RightHandSide<dim>::value(const Point<dim>     & /* p */,
                                const unsigned int   /* component */) const
{
  double result = 1.0;
  return result;
}

/*
 *  Neumann boundary condition
 */

template<int dim>
class NeumannBoundary : public Function<dim>
{
public:
  NeumannBoundary (const unsigned int n_components = 1,
                   const double       time = 0.)
    :
    Function<dim>(n_components, time)
  {}

  virtual ~NeumannBoundary(){};

  virtual double value (const Point<dim>    &p,
                        const unsigned int  component = 0) const;
};

template<int dim>
double NeumannBoundary<dim>::value(const Point<dim>   &/* p */,
                                   const unsigned int /* component */) const
{
  double result = 0.0;
  return result;
}

/*
 *  Velocity field
 */

template<int dim>
class VelocityField : public Function<dim>
{
public:
  VelocityField (const unsigned int n_components = dim,
                 const double       time = 0.)
    :
    Function<dim>(n_components, time)
  {}

  virtual ~VelocityField(){};

  virtual double value (const Point<dim>    &p,
                        const unsigned int  component = 0) const;
};

template<int dim>
double VelocityField<dim>::value(const Point<dim>   &point,
                                 const unsigned int component) const
{
  double value = 0.0;

  if(TYPE_VELOCITY_FIELD == TypeVelocityField::Constant)
  {
    // constant velocity field (u,v) = (1,1)
    value = 1.0;
  }
  else if(TYPE_VELOCITY_FIELD == TypeVelocityField::Circular)
  {
    // circular velocity field (u,v) = (-y,x)
    if(component == 0)
      value = - point[1];
    else if(component == 1)
      value = point[0];
    else
      AssertThrow(component <= 1, ExcMessage("Velocity field for 3-dimensional problem is not implemented!"));
  }
  else if(TYPE_VELOCITY_FIELD == TypeVelocityField::CircularZeroAtBoundary)
  {
    const double pi = numbers::PI;
    double sinx = std::sin(pi*point[0]);
    double siny = std::sin(pi*point[1]);
    double sin2x = std::sin(2.*pi*point[0]);
    double sin2y = std::sin(2.*pi*point[1]);
    if (component == 0)
      value = pi*sin2y*std::pow(sinx,2.);
    else if (component == 1)
      value = -pi*sin2x*std::pow(siny,2.);
  }
  else
  {
    AssertThrow(false, ExcMessage("Invalid type of velocity field prescribed for this problem."));
  }

  return value;
}

/**************************************************************************************/
/*                                                                                    */
/*         GENERATE GRID, SET BOUNDARY INDICATORS AND FILL BOUNDARY DESCRIPTOR        */
/*                                                                                    */
/**************************************************************************************/

template<int dim>
void create_grid_and_set_boundary_conditions(
    parallel::distributed::Triangulation<dim>               &triangulation,
    unsigned int const                                      n_refine_space,
    std_cxx11::shared_ptr<BoundaryDescriptorConvDiff<dim> > boundary_descriptor)
{
  // hypercube: line in 1D, square in 2D, etc., hypercube volume is [left,right]^dim
  const double left = -1.0, right = 1.0;
  GridGenerator::hyper_cube(triangulation,left,right);

  triangulation.refine_global(n_refine_space);

  std_cxx11::shared_ptr<Function<dim> > analytical_solution;
  analytical_solution.reset(new AnalyticalSolution<dim>());
  boundary_descriptor->dirichlet_bc.insert(std::pair<types::boundary_id,std_cxx11::shared_ptr<Function<dim> > >(0,analytical_solution));
}

template<int dim>
void set_field_functions(std_cxx11::shared_ptr<FieldFunctionsConvDiff<dim> > field_functions)
{
  // initialize functions (analytical solution, rhs, boundary conditions)
  std_cxx11::shared_ptr<Function<dim> > analytical_solution;
  analytical_solution.reset(new AnalyticalSolution<dim>());

  std_cxx11::shared_ptr<Function<dim> > right_hand_side;
  right_hand_side.reset(new RightHandSide<dim>());

  std_cxx11::shared_ptr<Function<dim> > velocity;
  velocity.reset(new VelocityField<dim>());

  field_functions->analytical_solution = analytical_solution;
  field_functions->right_hand_side = right_hand_side;
  field_functions->velocity = velocity;
}

template<int dim>
void set_analytical_solution(std_cxx11::shared_ptr<AnalyticalSolutionConvDiff<dim> > analytical_solution)
{
  analytical_solution->solution.reset(new AnalyticalSolution<dim>(1));
}


#endif /* APPLICATIONS_CONVECTION_DIFFUSION_TEST_CASES_CONST_RHS_CONST_AND_CIRCULAR_WIND_H_ */