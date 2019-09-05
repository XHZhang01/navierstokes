/*
 * 3D_Poiseuille_Pressure_inflow_outflow.h
 *
 *  Created on: July 1, 2019
 *      Author: Zhang
 */

#ifndef APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_POISEUILLE_H_
#define APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_POISEUILLE_H_

#include "../../include/incompressible_navier_stokes/postprocessor/postprocessor.h"
#include "../../include/incompressible_navier_stokes/postprocessor/inflow_data_calculator.h"
#include "../../include/functionalities/one_sided_cylindrical_manifold.h"
#include "../grid_tools/dealii_extensions.h"

#include "../../include/incompressible_navier_stokes/postprocessor/flow_rate_calculator.h"
#include "convection_diffusion/user_interface/boundary_descriptor.h"
#include "convection_diffusion/user_interface/field_functions.h"
#include "convection_diffusion/user_interface/input_parameters.h"


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

double *getdata ()
{

// Declare a file variable and open the .asc file.
ifstream file;
file.open ("/home/zhang/workspace/navierstokes/applications/Patientendaten/cossmoothAtemwegsdruck.asc",ios::in);
/*
if(!file)
{
  cout << "failed";
  return false;
}
*/
string temp,temp2data;

double value;
vector <double> datavector;

// vector <double> Atemwegsdruck;
// vector <double> Fluss;
// vector <double> Oesophagusdruck;

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

double *p = new double [lengthofdatavector];

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

/************************************************************************************************************/
/*                                                                                                          */
/*                                              INPUT PARAMETERS                                            */
/*                                                                                                          */
/************************************************************************************************************/

// convergence studies in space or time
unsigned int const DEGREE_MIN = 3;
unsigned int const DEGREE_MAX = DEGREE_MIN;

unsigned int const REFINE_SPACE_MIN = 2;
unsigned int const REFINE_SPACE_MAX = REFINE_SPACE_MIN;

unsigned int const REFINE_TIME_MIN = 0;
unsigned int const REFINE_TIME_MAX = 0;

// problem specific parameters

// space dimensions
unsigned int const DIMENSION = 3;

// polynomial degree (velocity)
unsigned int const DEGREE_U = DEGREE_MIN;

// set the number of refine levels
unsigned int const REFINE_STEPS_SPACE = REFINE_SPACE_MIN + 1;

// set the throat Reynolds number Re_throat = U_{mean,throat} * (2 R_throat) / nu
double const RE = 8000; //500; //2000; //3500; //5000; //6500; //8000;

// output folders
std::string const OUTPUT_FOLDER = "output/poiseuille/Re8000/";
std::string const OUTPUT_FOLDER_VTU = OUTPUT_FOLDER + "vtu/";
std::string const OUTPUT_NAME = "3D_poiseuille_pressure_inflow_outflow_tube_trachea_d3r2n2";

// set problem specific parameters like physical dimensions, etc.
const double MAX_VELOCITY = 15.09;
double const VISCOSITY = 1.7e-5;  // m^2/s
double const DENSITY = 1.0;       // kg/m^3 (@ 20°C)

// radius
double const R_TUBE1 = 0.00375;
double const R_TUBE2 = 3*R_TUBE1;

// lengths (dimensions in flow direction z)
double const LENGTH_TUBE1 = 0.312;
double const LENGTH_TUBE2 = 0.06;

// z-coordinates
double const Z2_TUBE1 = LENGTH_TUBE1;
double const Z1_TUBE1 = 0.0;

double const Z1_TUBE2 = LENGTH_TUBE1;
double const Z2_TUBE2 = LENGTH_TUBE1 + LENGTH_TUBE2;

// mesh parameters
unsigned int const N_CELLS_AXIAL = 2;
unsigned int const N_CELLS_AXIAL_TUBE1 = 4*N_CELLS_AXIAL;
unsigned int const N_CELLS_AXIAL_TUBE2 = 10*N_CELLS_AXIAL;

unsigned int const MANIFOLD_ID_CYLINDER = 1234;

double const PERIOD = 3.0; // in period lasts 3 s
unsigned int const N_PERIODS = 1;
double const START_TIME = 49.50;
double const END_TIME = START_TIME+PERIOD*N_PERIODS;
double const PEEP_KINEMATIC = 1400; //14.8 * 98.0665 / DENSITY;      // 8 cmH20, 1 cmH20 = 98.0665 Pa, transform to kinematic pressure
//double const TIDAL_VOLUME = 500.0e-6;                       // 500 ml = 500 * 10^{-6} m^3
double const C_RS_KINEMATIC = DENSITY * 0.04464e-5;   // total respiratory compliance C_rs = 100 ml/cm H20
//double const DELTA_P_INITIAL = TIDAL_VOLUME/C_RS_KINEMATIC; // initialize pressure difference in order to obtain desired tidal volume

// inlet boundary IDs
types::boundary_id const INLET_ID_FIRST = 1;
types::boundary_id INLET_ID_LAST = 2;

// outlet boundary IDs
types::boundary_id const OUTLET_ID_FIRST = 2;
types::boundary_id OUTLET_ID_LAST = 3;


// start and end time
bool const WRITE_OUTPUT = true;
bool const HIGH_ORDER_OUTPUT = true;
double const OUTPUT_START_TIME = START_TIME;
double const OUTPUT_INTERVAL_TIME = PERIOD/300;


class OutflowBoundary
{
public:
  OutflowBoundary(types::boundary_id const id)
    :
      boundary_id(id),
      resistance(3.5e5), // in preliminary tests with 5 generations we used a constant value of 1.0e7
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
//std::map<types::boundary_id, double> FLOW_RATES;

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
    return (outflow_boundary->get_pressure());
  }

private:
  std::shared_ptr<OutflowBoundary> outflow_boundary;
};


namespace IncNS
{
void set_input_parameters(InputParameters &param)
{
  // MATHEMATICAL MODEL
  param.dim = DIMENSION;
  param.problem_type = ProblemType::Unsteady;
  param.equation_type = EquationType::NavierStokes;
  param.use_outflow_bc_convective_term = true;
  param.formulation_viscous_term = FormulationViscousTerm::LaplaceFormulation;
  param.formulation_convective_term = FormulationConvectiveTerm::DivergenceFormulation;
  param.right_hand_side = false;


  // PHYSICAL QUANTITIES
  param.start_time = START_TIME;
  param.end_time = END_TIME;
  param.viscosity = VISCOSITY;


  // TEMPORAL DISCRETIZATION
  param.solver_type = SolverType::Unsteady;

  param.temporal_discretization = TemporalDiscretization::BDFDualSplittingScheme;
  param.treatment_of_convective_term = TreatmentOfConvectiveTerm::Explicit;
  param.calculation_of_time_step_size = TimeStepCalculation::CFL;
  param.adaptive_time_stepping = true;
//  param.temporal_discretization = TemporalDiscretization::BDFPressureCorrection;
//  param.treatment_of_convective_term = TreatmentOfConvectiveTerm::Implicit;
//  param.calculation_of_time_step_size = TimeStepCalculation::CFL;
//  param.adaptive_time_stepping_limiting_factor = 3.0;
  param.max_velocity = MAX_VELOCITY;
  param.cfl = 0.4;
  param.cfl_exponent_fe_degree_velocity = 1.5;
  param.time_step_size = 1.0e-1;
  param.order_time_integrator = 2;
  param.start_with_low_order = true;
  param.dt_refinements = 0;

  // output of solver information
  param.solver_info_data.print_to_screen = true;
  param.solver_info_data.interval_time = OUTPUT_INTERVAL_TIME;


  // SPATIAL DISCRETIZATION
  param.triangulation_type = TriangulationType::Distributed;
  param.degree_u = DEGREE_U;
  param.degree_p = DegreePressure::MixedOrder;
  param.mapping = MappingType::Isoparametric;

  param.h_refinements = REFINE_STEPS_SPACE;

  // convective term
  param.upwind_factor = 0.5;

  // viscous term
  param.IP_formulation_viscous = InteriorPenaltyFormulation::SIPG;
  param.IP_factor_viscous = 1.0;

  // special case: pure DBC's
  param.pure_dirichlet_bc = false;

  // div-div and continuity penalty terms
  param.use_divergence_penalty = true;
  param.divergence_penalty_factor = 1.0e0;
  param.use_continuity_penalty = true;
  param.continuity_penalty_factor = param.divergence_penalty_factor;
  param.add_penalty_terms_to_monolithic_system = false;

  // TURBULENCE
  param.use_turbulence_model = false;
  param.turbulence_model = TurbulenceEddyViscosityModel::Sigma;
  // Smagorinsky: 0.165, Vreman: 0.28, WALE: 0.50, Sigma: 1.35
  param.turbulence_model_constant = 1.35;

  // PROJECTION METHODS

  // pressure Poisson equation
  param.solver_data_pressure_poisson = SolverData(1000,1.e-12,1.e-3,100);
  param.preconditioner_pressure_poisson = PreconditionerPressurePoisson::Multigrid;
  param.multigrid_data_pressure_poisson.type = MultigridType::phMG;
  param.multigrid_data_pressure_poisson.p_sequence = PSequenceType::Bisect;
  param.multigrid_data_pressure_poisson.dg_to_cg_transfer = DG_To_CG_Transfer::Fine;

  // Variant 1:

//  param.multigrid_data_pressure_poisson.coarse_problem.solver = MultigridCoarseGridSolver::Chebyshev;
//  param.multigrid_data_pressure_poisson.coarse_problem.preconditioner = MultigridCoarseGridPreconditioner::PointJacobi;

  // Variant 2 (Standardfall):

  param.multigrid_data_pressure_poisson.coarse_problem.solver = MultigridCoarseGridSolver::CG;
  param.multigrid_data_pressure_poisson.coarse_problem.preconditioner = MultigridCoarseGridPreconditioner::AMG;

  // projection step
  param.solver_projection = SolverProjection::CG;
  param.solver_data_projection = SolverData(1000, 1.e-12, 1.e-3);
  param.preconditioner_projection = PreconditionerProjection::InverseMassMatrix;


  // HIGH-ORDER DUAL SPLITTING SCHEME

  // formulations
  param.order_extrapolation_pressure_nbc = param.order_time_integrator <=2 ? param.order_time_integrator : 2;

  // viscous step
  param.solver_viscous = SolverViscous::CG;
  param.solver_data_viscous = SolverData(1000,1.e-12,1.e-3);
  param.preconditioner_viscous = PreconditionerViscous::InverseMassMatrix;


  // PRESSURE-CORRECTION SCHEME

  // formulation
  param.order_pressure_extrapolation = 1; // use 0 for non-incremental formulation
  param.rotational_formulation = true; // use false for standard formulation

  // momentum step

  // Newton solver
  param.newton_solver_data_momentum = NewtonSolverData(100,1.e-14,1.e-6);

  // linear solver
  param.solver_momentum = SolverMomentum::GMRES;
  param.solver_data_momentum = SolverData(1e4, 1.e-12, 1.e-3, 100);
  param.preconditioner_momentum = MomentumPreconditioner::InverseMassMatrix;
  param.update_preconditioner_momentum = false;

  // formulation
  param.order_pressure_extrapolation = 1;
  param.rotational_formulation = true;

}

}

/************************************************************************************************************/
/*                                                                                                          */
/*                                       CREATE GRID AND SET BOUNDARY IDs                                   */
/*                                                                                                          */
/************************************************************************************************************/

template<int dim>
void create_grid_and_set_boundary_ids(
    std::shared_ptr<parallel::Triangulation<dim>>     triangulation,
    unsigned int const                                n_refine_space,
    std::vector<GridTools::PeriodicFacePair<typename
      Triangulation<dim>::cell_iterator> >            &/*periodic_faces*/)
{
  /*
   *   TUBE 1
   */
  Triangulation<2> tria_2d_tube1;
  Triangulation<dim> tria_tube1;
  GridGenerator::hyper_ball(tria_2d_tube1, Point<2>(), R_TUBE1);
  GridGenerator::extrude_triangulation(tria_2d_tube1,N_CELLS_AXIAL_TUBE1+1,LENGTH_TUBE1,tria_tube1);
  Tensor<1,dim> offset_tube1;
  offset_tube1[2] = Z1_TUBE1;
  GridTools::shift(offset_tube1,tria_tube1);

  Triangulation<dim> * current_tria = &tria_tube1;


  /*
   *   TUBE 2
   */
  const unsigned int n_cells_circle = 4;
  double const R_1 = R_TUBE1 + 1.0/3.0*(R_TUBE2-R_TUBE1);
  double const R_2 = R_TUBE1 + 2.0/3.0*(R_TUBE2-R_TUBE1);

  Triangulation<2> tria_2d_tube2_inner, circle_1, circle_2, circle_3, tria_tmp_2d_1, tria_tmp_2d_2, tria_2d_tube2;
  GridGenerator::hyper_ball(tria_2d_tube2_inner, Point<2>(), R_TUBE1);

  GridGenerator::hyper_shell(circle_1, Point<2>(), R_TUBE1, R_1, n_cells_circle, true);
  GridTools::rotate(numbers::PI/4, circle_1);
  GridGenerator::hyper_shell(circle_2, Point<2>(), R_1, R_2, n_cells_circle, true);
  GridTools::rotate(numbers::PI/4, circle_2);
  GridGenerator::hyper_shell(circle_3, Point<2>(), R_2, R_TUBE2, n_cells_circle, true);
  GridTools::rotate(numbers::PI/4, circle_3);

  // merge 2d triangulations
  GridGenerator::merge_triangulations (tria_2d_tube2_inner, circle_1, tria_tmp_2d_1);
  GridGenerator::merge_triangulations (circle_2, circle_3, tria_tmp_2d_2);
  GridGenerator::merge_triangulations (tria_tmp_2d_1, tria_tmp_2d_2, tria_2d_tube2);

  // extrude in z-direction
  Triangulation<dim> tria_tube2;
  GridGenerator::extrude_triangulation(tria_2d_tube2,N_CELLS_AXIAL_TUBE2+1,LENGTH_TUBE2,tria_tube2);
  Tensor<1,dim> offset_tube2; offset_tube2[2] = Z1_TUBE2;
  GridTools::shift(offset_tube2,tria_tube2);

  /*
   *  MERGE TRIANGULATIONS
   */
  GridGenerator::merge_triangulations (tria_tube1, tria_tube2, *triangulation);

  /*
   *  MANIFOLDS
   */
  current_tria = &(*triangulation);
  current_tria->set_all_manifold_ids(0);

  // first fill vectors of manifold_ids and face_ids
  std::vector<unsigned int> manifold_ids;
  std::vector<unsigned int> face_ids;

  for (typename Triangulation<dim>::cell_iterator cell = current_tria->begin();cell != current_tria->end(); ++cell)
  {
    // TUBE 1
    if(cell->center()[2] < Z2_TUBE1)
    {
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      {
        bool face_at_sphere_boundary = true;
        for (unsigned int v=0; v<GeometryInfo<dim-1>::vertices_per_cell; ++v)
        {
          Point<dim> point = Point<dim>(0,0,cell->face(f)->vertex(v)[2]);
          if (std::abs((cell->face(f)->vertex(v)-point).norm()-R_TUBE1) > 1e-12)
          {
            face_at_sphere_boundary = false;
          }
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

    // TUBE 2
    else if(cell->center()[2] > Z1_TUBE2 && cell->center()[2] < Z2_TUBE2)
    {
      Point<dim> point2 = Point<dim>(0,0,cell->center()[2]);

      // cylindrical manifold for outer cell layers
      if((cell->center()-point2).norm() > R_TUBE1/std::sqrt(2.0))
        cell->set_all_manifold_ids(MANIFOLD_ID_CYLINDER);

      // one-sided cylindrical manifold for core region
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      {
        bool face_at_sphere_boundary = true;
        for (unsigned int v=0; v<GeometryInfo<dim-1>::vertices_per_cell; ++v)
        {
          Point<dim> point = Point<dim>(0,0,cell->face(f)->vertex(v)[2]);
          if (std::abs((cell->face(f)->vertex(v)-point).norm()-R_TUBE1) > 1e-12 ||
              (cell->center()-point2).norm() > R_TUBE1/std::sqrt(2.0))
          {
            face_at_sphere_boundary = false;
          }
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
    else
    {
      AssertThrow(false, ExcMessage("Should not arrive here."));
    }
  }


  // one-sided spherical manifold
  // generate vector of manifolds and apply manifold to all cells that have been marked
  static std::vector<std::shared_ptr<Manifold<dim> > > manifold_vec;
  manifold_vec.resize(manifold_ids.size());

  for(unsigned int i=0;i<manifold_ids.size();++i)
  {
    for (typename Triangulation<dim>::cell_iterator cell = current_tria->begin(); cell != current_tria->end(); ++cell)
    {
      if(cell->manifold_id() == manifold_ids[i])
      {
        Point<dim> center = Point<dim>();
        manifold_vec[i] = std::shared_ptr<Manifold<dim> >(
            static_cast<Manifold<dim>*>(new OneSidedCylindricalManifold<dim>(cell,face_ids[i],center)));
        current_tria->set_manifold(manifold_ids[i],*(manifold_vec[i]));
      }
    }
  }

  // set cylindrical manifold
  static std::shared_ptr<Manifold<dim> > cylinder_manifold;
  cylinder_manifold = std::shared_ptr<Manifold<dim> >(static_cast<Manifold<dim>*>(new MyCylindricalManifold<dim>(Point<dim>())));
  current_tria->set_manifold(MANIFOLD_ID_CYLINDER, *cylinder_manifold);



  /*
   *  BOUNDARY ID's
   */

  typename Triangulation<dim>::cell_iterator cell = current_tria->begin(), endc = current_tria->end();
  for(;cell!=endc;++cell)
  {
    for(unsigned int face_number=0; face_number < GeometryInfo<dim>::faces_per_cell; ++face_number)
    {
      // left boundary
      if ((std::fabs(cell->face(face_number)->center()[2] - Z1_TUBE1) < 1e-12))
      {
        cell->face(face_number)->set_boundary_id (1);
      }

      // right boundary
      if ((std::fabs(cell->face(face_number)->center()[2] - Z2_TUBE2) < 1e-12))
      {
        cell->face(face_number)->set_boundary_id (2);
      }
    }
  }

  // perform global refinements
  triangulation->refine_global(n_refine_space);
}


/************************************************************************************************************/
/*                                                                                                          */
/*                         FUNCTIONS (INITIAL/BOUNDARY CONDITIONS, RIGHT-HAND SIDE, etc.)                   */
/*                                                                                                          */
/************************************************************************************************************/

namespace IncNS
{

template<int dim>
class PressureInflowBC : public Function<dim>
{
public:
  double *pointer;

  PressureInflowBC (const double time = 0.)
    :
    Function<dim>(1 /*n_components*/, time)
  {pointer = getdata();}

  double value (const Point<dim>   &p,
                  const unsigned int /*component = 0*/) const
  {
    /*
    double t = this->get_time();

    // assuming a parabolic velocity profile
    //double const pressure_gradient = -2.*VISCOSITY*MAX_VELOCITY;
    double const pi = numbers::PI;
    double const T = 3.0;

    //return (p[0]-L) * pressure_gradient * std::sin(pi*t/T);
    //return (21 + 28/pi * (std::sin(2.31*t) + std::sin(3*2.31*t)/3 + std::sin(5*2.31*t)/5 + std::sin(7+2.31*t)/7 + std::sin(9*2.31*t)/9 + std::sin(11*2.31*t)/11) + std::sin(13*2.31*t)/13 + std::sin(15*2.31*t)/15 + std::sin(17*2.31*t)/17 + std::sin(19*2.31*t)/19) * 98.0665;
    return (21 + 7 * std::sin(2*pi*t/T)) * 98.0665;
    */



    double t = this->get_time();
    // using measured values of Atemwegsdruck (airway pressure)
    double times = t*100;
    // std::cout << "\nTime: " << t;

    // float *pointer = getdata();

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
  for(types::boundary_id id = INLET_ID_FIRST; id < INLET_ID_LAST; ++id)
  {

    //std::shared_ptr<OutflowBoundary> inflow_boundary;
    //inflow_boundary.reset(new OutflowBoundary(id));
    //INFLOW_BOUNDARIES.push_back(inflow_boundary);


    boundary_descriptor_velocity->neumann_bc.insert(pair(id, new Functions::ZeroFunction<dim>(dim)));
    boundary_descriptor_pressure->dirichlet_bc.insert(pair(id, new PressureInflowBC<dim>()));
  }

  // 2 = outlet
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

/************************************************************************************************************/
/*                                                                                                          */
/*                                              POSTPROCESSOR                                               */
/*                                                                                                          */
/************************************************************************************************************/

template<int dim>
struct PostProcessorDataLung
{
  PostProcessorData<dim> pp_data;
  FlowRateCalculatorData<dim> flow_rate_data;
};

template<int dim, typename Number>
class PostProcessorLung : public PostProcessor<dim, Number>
{
public:
  typedef PostProcessor<dim, Number> Base;

  typedef LinearAlgebra::distributed::Vector<Number> VectorType;

  typedef typename Base::Operator Operator;

  PostProcessorLung(PostProcessorDataLung<dim> const & pp_data_in)
    :
    Base(pp_data_in.pp_data),
    pp_data_lung(pp_data_in),
    time_last(START_TIME)
  {
  }

  unsigned int OUTFLOW_ID = 2;
  unsigned int INFLOW_ID = 1;

  void setup(Operator const & pde_operator)
  {
    // call setup function of base class
    Base::setup(pde_operator);

    // fill flow_rates map

    flow_rates.insert(std::pair<types::boundary_id, double>(OUTFLOW_ID,0.0));

    flow_rates.insert(std::pair<types::boundary_id, double>(INFLOW_ID,0.0));

    /*
    for(auto iterator_out = OUTFLOW_BOUNDARIES.begin(); iterator_out != OUTFLOW_BOUNDARIES.end(); ++iterator_out)
    {
      flow_rates.insert(std::pair<types::boundary_id, double>((*iterator_out)->get_boundary_id(),0.0));
    }
    */

    /*
    for(auto iterator_in = INFLOW_BOUNDARIES.begin(); iterator_in != INFLOW_BOUNDARIES.end(); ++iterator_in)
    {
      flow_rates.insert(std::pair<types::boundary_id, double>((*iterator_in)->get_boundary_id(),0.0));
    }
    */

    flow_rate_calculator.reset(new FlowRateCalculator<dim,Number>(
        pde_operator.get_matrix_free(),
        pde_operator.get_dof_handler_u(),
        pde_operator.get_dof_index_velocity(),
        pde_operator.get_quad_index_velocity_linear(),
        pp_data_lung.flow_rate_data));
  }

  void do_postprocessing(VectorType const &velocity,
                         VectorType const &pressure,
                         double const     time,
                         int const        time_step_number)
  {
    Base::do_postprocessing(
        velocity,
        pressure,
        time,
        time_step_number);

    // calculate flow rates for all outflow boundaries
    AssertThrow(pp_data_lung.flow_rate_data.calculate == true,
        ExcMessage("Activate flow rate computation."));

    flow_rate_calculator->calculate_flow_rates(velocity, time, flow_rates);

    // set flow rate for all outflow boundaries and update volume (i.e., integrate flow rate over time)
    Number volume = 0.0;

    for(auto iterator_out = OUTFLOW_BOUNDARIES.begin(); iterator_out != OUTFLOW_BOUNDARIES.end(); ++iterator_out)
    {
      (*iterator_out)->set_flow_rate(flow_rates.at((*iterator_out)->get_boundary_id()));
      (*iterator_out)->integrate_volume(time);
      volume += (*iterator_out)->get_volume();
    }

    /*
    for(auto iterator_in = INFLOW_BOUNDARIES.begin(); iterator_in != INFLOW_BOUNDARIES.end(); ++iterator_in)
    {
       (*iterator_in)->set_flow_rate(flow_rates.at((*iterator_in)->get_boundary_id()));
       (*iterator_in)->integrate_volume(time);
       volume += (*iterator_in)->get_volume();
    }
    */

    // write flow rates to file

    std::ostringstream filename_in;
    filename_in << OUTPUT_FOLDER + OUTPUT_NAME + "_flow_rate_inflow";
    write_output(flow_rates.at(INFLOW_ID), time, "Flow rate in [m^3/sec]", time_step_number, filename_in);

    std::ostringstream filename_out;
    filename_out << OUTPUT_FOLDER + OUTPUT_NAME + "_flow_rate_outflow";
    write_output(flow_rates.at(OUTFLOW_ID), time, "Flow rate in [m^3/sec]", time_step_number, filename_out);

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

  // we need to compute the flow rate for each outlet
  std::map<types::boundary_id, double> flow_rates;

  // calculate flow rates for all outflow boundaries
  std::shared_ptr<FlowRateCalculator<dim,Number> > flow_rate_calculator;

  double time_last;
};

template<int dim, typename Number>
std::shared_ptr<PostProcessorBase<dim, Number> >
construct_postprocessor(InputParameters const &param)
{
  (void)param;

  PostProcessorData<dim> pp_data;

  // write output for visualization of results
  pp_data.output_data.write_output = WRITE_OUTPUT;
  pp_data.output_data.output_folder = OUTPUT_FOLDER_VTU;
  pp_data.output_data.output_name = OUTPUT_NAME;
  pp_data.output_data.output_start_time = OUTPUT_START_TIME;
  pp_data.output_data.output_interval_time = OUTPUT_INTERVAL_TIME;
  pp_data.output_data.write_vorticity = false;
  pp_data.output_data.write_divergence = false;
  pp_data.output_data.write_velocity_magnitude = true;
  pp_data.output_data.write_vorticity_magnitude = false;
  pp_data.output_data.write_q_criterion = false;
  pp_data.output_data.write_processor_id = false;
  pp_data.output_data.degree = param.degree_u;
  pp_data.output_data.write_higher_order = HIGH_ORDER_OUTPUT;

  // Lung specific modules
  PostProcessorDataLung<dim> pp_data_lung;
  pp_data_lung.pp_data = pp_data;

  // calculation of flow rate
  pp_data_lung.flow_rate_data.calculate = true;
  pp_data_lung.flow_rate_data.write_to_file = true;
  pp_data_lung.flow_rate_data.filename_prefix = OUTPUT_FOLDER + OUTPUT_NAME + "_flow_rate";

  std::shared_ptr<PostProcessorBase<dim,Number> > pp;
  pp.reset(new PostProcessorLung<dim,Number>(pp_data_lung));

  return pp;
}

}

#endif /* APPLICATIONS_INCOMPRESSIBLE_NAVIER_STOKES_TEST_CASES_POISEUILLE_H_ */
