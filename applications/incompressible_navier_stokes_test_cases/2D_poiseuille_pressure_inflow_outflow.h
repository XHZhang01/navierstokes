/*
 * Poiseuille.h
 *
 *  Created on: May 20, 2016
 *      Author: Zhang
 */

#ifndef APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_POISEUILLE_H_
#define APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_POISEUILLE_H_

#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdio.h>

using namespace std;

// This function is used to get measured values from .asc file and generate a vector with measured values.
// There are 4 options, Atemwegsdruck (airway pressure), Fluss (volumetric flow) and Oesophagusdruck (Esophagus pressure).


/**************************************************************************************/
/*                                                                                    */
/*                                 GET DATA FROM .ASC FILE                            */
/*                                                                                    */
/**************************************************************************************/

float* getdata ()
{

// Declare a file variable and open the .asc file.
ifstream file;
file.open ("/scratch/zhang/workspace/navierstokes/applications/Patientendaten/cossmoothAtemwegsdruck.asc",ios::in);
/*
if(!file)
{
  cout << "failed";
  return false;
}
*/
string temp,temp2data;

float value;
vector <float> datavector;

// vector <float> Atemwegsdruck;
// vector <float> Fluss;
// vector <float> Oesophagusdruck;

int lengthofdatavector = 0;

// Vector datavector contains all the measured values.
// Variable lengthofdatavector is the total number of measured values.
while (getline(file,temp))
{
  istringstream LineBand(temp);

  while (LineBand>>temp2data)
  {
    value = atof(temp2data.data());
    datavector.push_back(value);
    lengthofdatavector++;
//    printf ("%.2f\n",value);
  }
}

file.close();

float *p = new float [lengthofdatavector];

for (int i=0; i<lengthofdatavector; i++)
  p[i] = datavector[i];


//for (int i = 0;i <= (lengthofdatavector-1);i++)
//{
  //int remainder = (i-2)%5; // Get measured values of Atemwegsdruck (airway pressure) with unit cmH2O.
  //int remainder = (i-3)%5; // Get measured values of Fluss (volumetric flow) with unit L/min.
  //int remainder = (i-4)%5; // Get measured values of Oesophagusdruck (Esophagus pressure) with unit mmHg.

  //if (remainder == 0)
  //{
    //Atemwegsdruck.push_back(datavector[i]);
    //Fluss.push_back(datavector[i]);
    //Oesophagusdruck.push_back(datavector[i]);

    //printf ("%.2f\n" , datavector[i]);
  //}
//}

//float current_value  = datavector[position];
return p;
}

/**************************************************************************************/
/*                                                                                    */
/*                                 INPUT PARAMETERS                                   */
/*                                                                                    */
/**************************************************************************************/

// single or double precision?
//typedef float VALUE_TYPE;
typedef double VALUE_TYPE;

// set the number of space dimensions: dimension = 2, 3
unsigned int const DIMENSION = 2;

// set the polynomial degree of the shape functions for velocity and pressure
unsigned int const FE_DEGREE_VELOCITY = 3;
unsigned int const FE_DEGREE_PRESSURE = FE_DEGREE_VELOCITY-1;

// set the number of refine levels for spatial convergence tests
unsigned int const REFINE_STEPS_SPACE_MIN = 1;
unsigned int const REFINE_STEPS_SPACE_MAX = REFINE_STEPS_SPACE_MIN;

// set the number of refine levels for temporal convergence tests
unsigned int const REFINE_STEPS_TIME_MIN = 0;
unsigned int const REFINE_STEPS_TIME_MAX = REFINE_STEPS_TIME_MIN;

// set problem specific parameters like physical dimensions, etc.
const ProblemType PROBLEM_TYPE = ProblemType::Unsteady;
const double MAX_VELOCITY = 15.09;
double const VISCOSITY = 1.7e-5;  // m^2/s
double const DENSITY = 1;       // kg/m^3 (@ 20??C)

const double H = 0.0075;
const double L = 0.312;

double const PERIOD = 3; // in period lasts 3 s
unsigned int const N_PERIODS = 1;
double const START_TIME = 49.50;
double const END_TIME = START_TIME+PERIOD*N_PERIODS;
double const PEEP_KINEMATIC = 1475; //14.8 * 98.0665 / DENSITY;      // 8 cmH20, 1 cmH20 = 98.0665 Pa, transform to kinematic pressure
//double const TIDAL_VOLUME = 500.0e-6;                       // 500 ml = 500 * 10^{-6} m^3
double const C_RS_KINEMATIC = DENSITY * 0.04464e-5;   // total respiratory compliance C_rs = 100 ml/cm H20
//double const DELTA_P_INITIAL = TIDAL_VOLUME/C_RS_KINEMATIC; // initialize pressure difference in order to obtain desired tidal volume

// outlet boundary IDs
types::boundary_id const OUTLET_ID_FIRST = 2;
types::boundary_id OUTLET_ID_LAST = 3;

// output
bool const WRITE_OUTPUT = true;
bool const HIGH_ORDER_OUTPUT = true;
double const OUTPUT_START_TIME = START_TIME;
double const OUTPUT_INTERVAL_TIME = PERIOD/300;

std::string const OUTPUT_FOLDER = "output/poiseuille/Re8000/";
std::string const OUTPUT_FOLDER_VTU = OUTPUT_FOLDER + "vtu/";
std::string const OUTPUT_NAME = "2D_poiseuille_pressure_inflow_outflow 2";

// restart
//bool const WRITE_RESTART = false;
//double const RESTART_INTERVAL_TIME = PERIOD;

class OutflowBoundary
{
public:
  OutflowBoundary(types::boundary_id const id)
    :
      boundary_id(id),
      resistance(6.18e5), // in preliminary tests with 5 generations we used a constant value of 1.0e7
      compliance(C_RS_KINEMATIC), // note that one could use a statistical distribution as in Roth et al. (2018)
      //volume(compliance * PEEP_KINEMATIC), // p = 1/C * V -> V = C * p (initialize volume so that p(t=0) = PEEP_KINEMATIC)
      volume(0.0),
      flow_rate(0.0),
      time_old(START_TIME)
  {}

  void
  set_flow_rate(double const flow_rate_)
  {
    flow_rate = flow_rate_;
  }

  void
  integrate_volume(double const time)
  {
    // currently use BDF1 time integration // TODO one could use a higher order time integrator
    volume += flow_rate*(time-time_old);
    time_old = time;
  }

  double
  get_pressure() const
  {
    return resistance*flow_rate + volume/compliance + PEEP_KINEMATIC;
  }

  double
  get_volume() const
  {
    return volume;
  }

  types::boundary_id get_boundary_id() const
  {
    return boundary_id;
  }

private:
  types::boundary_id const boundary_id;
  double resistance;
  double compliance;
  double volume;
  double flow_rate;
  double time_old;
};

// we need individual outflow boundary conditions for each outlet
std::vector<std::shared_ptr<OutflowBoundary>> OUTFLOW_BOUNDARIES;

// we need to compute the flow rate for each outlet
std::map<types::boundary_id, double> FLOW_RATES;

template<int dim>
class PressureOutlet : public Function<dim>
{
public:
  PressureOutlet (std::shared_ptr<OutflowBoundary> outflow_boundary_,
                  double const time = 0.)
    :
    Function<dim>(1 /*n_components*/, time),
    outflow_boundary(outflow_boundary_)
  {}

  double value (const Point<dim>   &/*p*/,
                const unsigned int /*component*/) const
  {
    return outflow_boundary->get_pressure();
  }

private:
  std::shared_ptr<OutflowBoundary> outflow_boundary;
};

template<int dim>
void InputParameters<dim>::set_input_parameters()
{
  // MATHEMATICAL MODEL
  problem_type = PROBLEM_TYPE; // PROBLEM_TYPE is also needed somewhere else
  equation_type = EquationType::NavierStokes;
  formulation_viscous_term = FormulationViscousTerm::LaplaceFormulation;
  formulation_convective_term = FormulationConvectiveTerm::DivergenceFormulation;
  use_outflow_bc_convective_term = true;
  right_hand_side = false;


  // PHYSICAL QUANTITIES
  start_time = START_TIME;
  end_time = END_TIME;
  viscosity = VISCOSITY; // VISCOSITY is also needed somewhere else


  // TEMPORAL DISCRETIZATION
  solver_type = SolverType::Unsteady;
  temporal_discretization = TemporalDiscretization::BDFDualSplittingScheme;
  treatment_of_convective_term = TreatmentOfConvectiveTerm::Explicit;
  calculation_of_time_step_size = TimeStepCalculation::CFL;
  adaptive_time_stepping = true;
  max_velocity = MAX_VELOCITY; // MAX_VELOCITY is also needed somewhere else
  cfl = 3.0e-1;
  time_step_size = 1.0e-1;
  //max_number_of_time_steps = 1e8;
  cfl_exponent_fe_degree_velocity = 1.5;
  order_time_integrator = 2; // 1; // 2; // 3;
  start_with_low_order = true; // true; // false;

  //convergence_criterion_steady_problem = ConvergenceCriterionSteadyProblem::SolutionIncrement;
  //abs_tol_steady = 1.e-12;
  //rel_tol_steady = 1.e-3;

  // SPATIAL DISCRETIZATION

  // triangulation
  triangulation_type = TriangulationType::Distributed;

  // mapping
  degree_mapping = FE_DEGREE_VELOCITY;

  // convective term
  if(formulation_convective_term == FormulationConvectiveTerm::DivergenceFormulation)
    upwind_factor = 0.5;

  // viscous term
  IP_formulation_viscous = InteriorPenaltyFormulation::SIPG;

  // special case: pure DBC's
  pure_dirichlet_bc = false;

  // div-div and continuity penalty
  use_divergence_penalty = true;
  divergence_penalty_factor = 1.0e0;
  use_continuity_penalty = true;
  continuity_penalty_factor = divergence_penalty_factor;
  add_penalty_terms_to_monolithic_system = false;

  // PROJECTION METHODS
  /*
  // pressure Poisson equation
  solver_pressure_poisson = SolverPressurePoisson::CG;
  solver_data_pressure_poisson = SolverData(1000,1.e-12,1.e-3,100);
  preconditioner_pressure_poisson = PreconditionerPressurePoisson::Multigrid;

  // projection step
  solver_projection = SolverProjection::CG;
  solver_data_projection = SolverData(1000, 1.e-12, 1.e-3);
  preconditioner_projection = PreconditionerProjection::InverseMassMatrix;
  */
  // pressure Poisson equation
  solver_data_pressure_poisson = SolverData(1000,1.e-12,1.e-3,100);
  preconditioner_pressure_poisson = PreconditionerPressurePoisson::Multigrid;
  multigrid_data_pressure_poisson.type = MultigridType::phMG;
  multigrid_data_pressure_poisson.p_sequence = PSequenceType::Bisect;
  multigrid_data_pressure_poisson.dg_to_cg_transfer = DG_To_CG_Transfer::Fine;

  // Variant 1:

  //multigrid_data_pressure_poisson.coarse_problem.solver = MultigridCoarseGridSolver::Chebyshev;
  //multigrid_data_pressure_poisson.coarse_problem.preconditioner = MultigridCoarseGridPreconditioner::PointJacobi;

  // Variant 2 (Standardfall):

  multigrid_data_pressure_poisson.coarse_problem.solver = MultigridCoarseGridSolver::CG;
  multigrid_data_pressure_poisson.coarse_problem.preconditioner = MultigridCoarseGridPreconditioner::AMG;

  // projection step
  solver_projection = SolverProjection::CG;
  solver_data_projection = SolverData(1000, 1.e-12, 1.e-3);
  preconditioner_projection = PreconditionerProjection::InverseMassMatrix;

  // HIGH-ORDER DUAL SPLITTING SCHEME

  // formulations
  order_extrapolation_pressure_nbc = order_time_integrator <=2 ? order_time_integrator : 2;

  // viscous step
  solver_viscous = SolverViscous::CG;
  solver_data_viscous = SolverData(1000,1.e-12,1.e-3);
  preconditioner_viscous = PreconditionerViscous::InverseMassMatrix;

  // PRESSURE-CORRECTION SCHEME

  // momentum step

  // Newton solver
  newton_solver_data_momentum = NewtonSolverData(100,1.e-14,1.e-6);

  // linear solver
  solver_momentum = SolverMomentum::GMRES;
  solver_data_momentum = SolverData(1e4, 1.e-12, 1.e-3, 100);
  preconditioner_momentum = MomentumPreconditioner::InverseMassMatrix;
  update_preconditioner_momentum = false;

  // formulation
  order_pressure_extrapolation = 1;
  rotational_formulation = true;

  // OUTPUT AND POSTPROCESSING

  // write output for visualization of results
  output_data.write_output = true;
  output_data.output_folder = OUTPUT_FOLDER;
  output_data.output_name = OUTPUT_NAME;
  output_data.output_start_time = OUTPUT_START_TIME;
  output_data.output_interval_time = OUTPUT_INTERVAL_TIME;
  output_data.write_velocity_magnitude = true;
  output_data.degree = FE_DEGREE_VELOCITY;

  // output of solver information
  solver_info_data.print_to_screen = true;
  solver_info_data.interval_time = OUTPUT_INTERVAL_TIME;
}

/**************************************************************************************/
/*                                                                                    */
/*                        GENERATE GRID AND SET BOUNDARY INDICATORS                   */
/*                                                                                    */
/**************************************************************************************/

template<int dim>
void create_grid_and_set_boundary_ids(
    std::shared_ptr<parallel::Triangulation<dim>>     triangulation,
    unsigned int const                                n_refine_space,
    std::vector<GridTools::PeriodicFacePair<typename
      Triangulation<dim>::cell_iterator> >            &/*periodic_faces*/)
{
  std::vector<unsigned int> repetitions({1,1});
  Point<dim> point1(0.0,-H/2.), point2(L,H/2.);
  GridGenerator::subdivided_hyper_rectangle(*triangulation,repetitions,point1,point2);

  // left boundary has ID=1, right boundary has ID=2, all other boundaries have ID=0
  typename Triangulation<dim>::cell_iterator cell = triangulation->begin(), endc = triangulation->end();
  for(;cell!=endc;++cell)
  {
    for(unsigned int face_number=0;face_number < GeometryInfo<dim>::faces_per_cell;++face_number)
    {
     if ((std::fabs(cell->face(face_number)->center()(0) - 0.0)< 1e-12))
         cell->face(face_number)->set_boundary_id (1);
     if ((std::fabs(cell->face(face_number)->center()(0) - L)< 1e-12))
        cell->face(face_number)->set_boundary_id (2);
    }
  }

  triangulation->refine_global(n_refine_space);
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
  float *pointer = getdata();
  PressureInflowBC (const double time = 0.)
    :
    Function<dim>(1 /*n_components*/, time)
  {}

  double value (const Point<dim>   &p,
                const unsigned int /*component = 0*/) const
  {
    /*
    double t = this->get_time();

    // assuming a parabolic velocity profile
    //double const pressure_gradient = -2.*VISCOSITY*MAX_VELOCITY;
    double const pi = numbers::PI;
    double const T = 3.0;
    //std::cout << "\nTime: " << t;

    //return (p[0]-L) * pressure_gradient * std::sin(pi*t/T);
    //return (21 + 28/pi * (std::sin(2.1*t) + std::sin(3*2.1*t)/3 + std::sin(5*2.1*t)/5 + std::sin(7+2.1*t)/7 + std::sin(9*2.1*t)/9 + std::sin(11*2.1*t)/11) + std::sin(13*2.1*t)/13 + std::sin(15*2.1*t)/15 + std::sin(17*2.1*t)/17 + std::sin(19*2.1*t)/19) * 98.0665;
    return (21 + 7 * std::sin(2*pi*t/T)) + 98.0665;
    */

    double t = this->get_time();
    // using measured values of Atemwegsdruck (airway pressure)
    double times = t*100;
    // std::cout << "\nTime: " << t;

    int floor = (int)times;
    // cout << "\nEntry: " << floor;
    int ceil = floor+1;
    // cout << ceil;
    double valueoffloor = pointer[floor];
    // cout << valueoffloor;
    double valueofceil = pointer[ceil];
    // cout << valueofceil;
    // cout << "\nPressure: " << (valueofceil*(times-floor)-valueoffloor*(times-ceil))*98.0665;
    return (valueofceil*(times-floor)-valueoffloor*(times-ceil))*98.0665;
    // cout << "Error 3";

  }
};


template<int dim>
void set_boundary_conditions(
    std::shared_ptr<BoundaryDescriptorU<dim> > boundary_descriptor_velocity,
    std::shared_ptr<BoundaryDescriptorP<dim> > boundary_descriptor_pressure)
{
  /*
  typedef typename std::pair<types::boundary_id,std::shared_ptr<Function<dim> > > pair;

  // fill boundary descriptor velocity
  boundary_descriptor_velocity->dirichlet_bc.insert(pair(0,new Functions::ZeroFunction<dim>(dim)));
  boundary_descriptor_velocity->neumann_bc.insert(pair(1,new Functions::ZeroFunction<dim>(dim)));
  boundary_descriptor_velocity->neumann_bc.insert(pair(2,new Functions::ZeroFunction<dim>(dim)));

  // fill boundary descriptor pressure
  boundary_descriptor_pressure->neumann_bc.insert(pair(0,new Functions::ZeroFunction<dim>(1)));
  boundary_descriptor_pressure->dirichlet_bc.insert(pair(1,new PressureInflowBC<dim>()));
  boundary_descriptor_pressure->dirichlet_bc.insert(pair(2,new Functions::ZeroFunction<dim>(1)));
  */

  // set boundary conditions
  typedef typename std::pair<types::boundary_id,std::shared_ptr<Function<dim> > > pair;

  // 0 = walls
  boundary_descriptor_velocity->dirichlet_bc.insert(pair(0, new Functions::ZeroFunction<dim>(dim)));
  boundary_descriptor_pressure->neumann_bc.insert(pair(0, new Functions::ZeroFunction<dim>(dim)));

  // 1 = inlet
  boundary_descriptor_velocity->neumann_bc.insert(pair(1, new Functions::ZeroFunction<dim>(dim)));
  boundary_descriptor_pressure->dirichlet_bc.insert(pair(1, new PressureInflowBC<dim>()));

  // outlets
  for(types::boundary_id id = OUTLET_ID_FIRST; id < OUTLET_ID_LAST; ++id)
  {
    std::shared_ptr<OutflowBoundary> outflow_boundary;
    outflow_boundary.reset(new OutflowBoundary(id));
    OUTFLOW_BOUNDARIES.push_back(outflow_boundary);

    boundary_descriptor_velocity->neumann_bc.insert(pair(id, new Functions::ZeroFunction<dim>(dim)));
    boundary_descriptor_pressure->dirichlet_bc.insert(pair(id, new PressureOutlet<dim>(outflow_boundary)));
  }

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

template<int dim>
struct PostProcessorDataLung
{
  PostProcessorData<dim> pp_data;
  FlowRateCalculatorData<dim> flow_rate_data;
};

template<int dim, int degree_u, int degree_p, typename Number>
class PostProcessorLung : public PostProcessor<dim, degree_u, degree_p, Number>
{
public:
  typedef PostProcessor<dim, degree_u, degree_p, Number> Base;

  typedef LinearAlgebra::distributed::Vector<Number> VectorType;

  typedef typename Base::NavierStokesOperator NavierStokesOperator;

  PostProcessorLung(PostProcessorDataLung<dim> const & pp_data_in)
    :
    Base(pp_data_in.pp_data),
    pp_data_lung(pp_data_in),
    time_last(START_TIME)
  {
  }

  void setup(NavierStokesOperator const                &navier_stokes_operator_in,
             DoFHandler<dim> const                     &dof_handler_velocity_in,
             DoFHandler<dim> const                     &dof_handler_pressure_in,
             Mapping<dim> const                        &mapping_in,
             MatrixFree<dim,Number> const              &matrix_free_data_in,
             DofQuadIndexData const                    &dof_quad_index_data_in,
             std::shared_ptr<AnalyticalSolution<dim> > analytical_solution_in)
  {
    // call setup function of base class
    Base::setup(
        navier_stokes_operator_in,
        dof_handler_velocity_in,
        dof_handler_pressure_in,
        mapping_in,
        matrix_free_data_in,
        dof_quad_index_data_in,
        analytical_solution_in);

    // fill FLOW_RATES map
    for(auto iterator = OUTFLOW_BOUNDARIES.begin(); iterator != OUTFLOW_BOUNDARIES.end(); ++iterator)
    {
      FLOW_RATES.insert(std::pair<types::boundary_id, double>((*iterator)->get_boundary_id(),0.0));
    }

    // flow rates:
    // In a first step, fill the set with boundary IDs. Note that we have to do that here and cannot do this step
    // in the function set_input_parameters() since the list of outflow boundaries is not known at that time.
    // The list of outflow boundaries is known once the grid has been created and the outflow boundary IDs have been set.
    for(auto iterator = OUTFLOW_BOUNDARIES.begin(); iterator != OUTFLOW_BOUNDARIES.end(); ++iterator)
    {
      pp_data_lung.flow_rate_data.boundary_IDs.insert((*iterator)->get_boundary_id());
    }

    flow_rate_calculator.reset(new FlowRateCalculator<dim,degree_u,Number>(
        matrix_free_data_in, dof_handler_velocity_in, dof_quad_index_data_in, pp_data_lung.flow_rate_data));
  }

  void do_postprocessing(VectorType const &velocity,
                         VectorType const &intermediate_velocity,
                         VectorType const &pressure,
                         double const     time,
                         int const        time_step_number)
  {
    Base::do_postprocessing(
        velocity,
        intermediate_velocity,
        pressure,
        time,
        time_step_number);

    // calculate flow rates for all outflow boundaries
    AssertThrow(pp_data_lung.flow_rate_data.calculate == true, ExcMessage("Activate flow rate computation."));
    flow_rate_calculator->calculate_flow_rates(velocity, time, FLOW_RATES);

    // set flow rate for all outflow boundaries and update volume (i.e., integrate flow rate over time)
    Number volume = 0.0;
    for(auto iterator = OUTFLOW_BOUNDARIES.begin(); iterator != OUTFLOW_BOUNDARIES.end(); ++iterator)
    {
      (*iterator)->set_flow_rate(FLOW_RATES.at((*iterator)->get_boundary_id()));
      (*iterator)->integrate_volume(time);
      volume += (*iterator)->get_volume();
    }

    // write volume to file
    if(pp_data_lung.flow_rate_data.write_to_file)
    {
      std::ostringstream filename;
      filename << OUTPUT_FOLDER + OUTPUT_NAME + "_volume";
      write_output(volume, time, "Volume in [m^3]", time_step_number, filename);

      // write time step size
      std::ostringstream filename_dt;
      filename_dt << OUTPUT_FOLDER + OUTPUT_NAME + "_time_step_size";
      write_output(time-time_last, time, "Time step size in [s]", time_step_number, filename_dt);
      time_last = time;
    }

  }

private:
  void
  write_output(double const &             value,
               double const &             time,
               std::string const &        name,
               unsigned int const         time_step_number,
               std::ostringstream const & filename)
  {
    // write output file
    if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      std::ofstream f;
      if(time_step_number == 1)
      {
        f.open(filename.str().c_str(), std::ios::trunc);
        f << std::endl << "  Time                " + name << std::endl;
      }
      else
      {
        f.open(filename.str().c_str(), std::ios::app);
      }

      unsigned int precision = 12;
      f << std::scientific << std::setprecision(precision) << std::setw(precision + 8) << time
        << std::setw(precision + 8) << value << std::endl;
    }
  }

  // postprocessor data supplemented with data required for lung
  PostProcessorDataLung<dim> pp_data_lung;

  // calculate flow rates for all outflow boundaries
  std::shared_ptr<FlowRateCalculator<dim,degree_u,Number> > flow_rate_calculator;

  double time_last;
};

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
  pp_data.kinetic_energy_data = param.kinetic_energy_data;
  pp_data.kinetic_energy_spectrum_data = param.kinetic_energy_spectrum_data;

  // Lung specific modules
  PostProcessorDataLung<dim> pp_data_lung;
  pp_data_lung.pp_data = pp_data;
  pp_data_lung.flow_rate_data = param.flow_rate_data;

  std::shared_ptr<PostProcessorLung<dim,degree_u,degree_p,Number> > pp;
  pp.reset(new PostProcessorLung<dim,degree_u,degree_p,Number>(pp_data_lung));

  return pp;
}

#endif /* APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_POISEUILLE_H_ */
