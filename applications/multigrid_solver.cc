//---------------------------------------------------------------------------
//    $Id: program.cc 56 2015-02-06 13:05:10Z kronbichler $
//    Version: $Name$
//
//    Copyright (C) 2013 - 2015 by Katharina Kormann and Martin Kronbichler
//
//---------------------------------------------------------------------------

// This file acts as a test for the correctness of the Poisson solver

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/integrators/laplace.h>

#include <fstream>
#include <sstream>

#include "../include/poisson/laplace_operator.h"

#include "../include/poisson/multigrid_preconditioner_laplace.h"
#include "solvers_and_preconditioners/iterative_solvers.h"

namespace Step37
{
using namespace dealii;

// dimension
const unsigned int DIM = 2;
// polynomial degree of shape functions
const unsigned int FE_DEGREE = 2;

// penalty factor interior penalty method
const double PENALTY_FACTOR = 0.25;

// mesh refinement strategy
const bool         RUN_VARIABLE_SIZES = false; // true;
const unsigned int N_REFINEMENTS      = 5;     // only relevant if RUN_VARIABLE_SIZES == false
const bool         DO_ADAPTIVE        = false;

// number of repetitions used to calculate minimum wall time
const unsigned int N_REPETITIONS = 3;


// Definition of analytic solution for testing the Poisson solver
template<int dim>
class SolutionBase
{
protected:
  static const unsigned int n_source_centers = 3;
  static const Point<dim>   source_centers[n_source_centers];
  static const double       width;
};


template<>
const Point<1> SolutionBase<1>::source_centers[SolutionBase<1>::n_source_centers] =
  {Point<1>(-1.0 / 3.0), Point<1>(0.0), Point<1>(+1.0 / 3.0)};

template<>
const Point<2> SolutionBase<2>::source_centers[SolutionBase<2>::n_source_centers] =
  {Point<2>(-0.5, +0.5), Point<2>(-0.5, -0.5), Point<2>(+0.5, -0.5)};

template<>
const Point<3> SolutionBase<3>::source_centers[SolutionBase<3>::n_source_centers] =
  {Point<3>(-0.5, +0.5, 0.25), Point<3>(-0.6, -0.5, -0.125), Point<3>(+0.5, -0.5, 0.5)};


template<int dim>
const double SolutionBase<dim>::width = 1. / 3.;



template<int dim>
class Solution : public Function<dim>, protected SolutionBase<dim>
{
public:
  Solution() : Function<dim>()
  {
  }

  virtual ~Solution()
  {
  }

  virtual double
  value(const Point<dim> & p, const unsigned int component = 0) const;

  virtual Tensor<1, dim>
  gradient(const Point<dim> & p, const unsigned int component = 0) const;
};



template<int dim>
double
Solution<dim>::value(const Point<dim> & p, const unsigned int) const
{
  const double pi           = numbers::PI;
  double       return_value = 0;
  for(unsigned int i = 0; i < this->n_source_centers; ++i)
  {
    const Tensor<1, dim> x_minus_xi = p - this->source_centers[i];
    return_value += std::exp(-x_minus_xi.norm_square() / (this->width * this->width));
  }

  return return_value / Utilities::fixed_power<dim>(std::sqrt(2 * pi) * this->width);
}



template<int dim>
Tensor<1, dim>
Solution<dim>::gradient(const Point<dim> & p, const unsigned int) const
{
  const double   pi = numbers::PI;
  Tensor<1, dim> return_value;

  for(unsigned int i = 0; i < this->n_source_centers; ++i)
  {
    const Tensor<1, dim> x_minus_xi = p - this->source_centers[i];

    return_value +=
      (-2 / (this->width * this->width) *
       std::exp(-x_minus_xi.norm_square() / (this->width * this->width)) * x_minus_xi);
  }

  return return_value / Utilities::fixed_power<dim>(std::sqrt(2 * pi) * this->width);
}



template<int dim>
class RightHandSide : public Function<dim>, protected SolutionBase<dim>
{
public:
  RightHandSide() : Function<dim>()
  {
  }

  virtual ~RightHandSide()
  {
  }

  virtual double
  value(const Point<dim> & p, const unsigned int component = 0) const;
};


template<int dim>
double
RightHandSide<dim>::value(const Point<dim> & p, const unsigned int) const
{
  const double pi           = numbers::PI;
  double       return_value = 0;
  for(unsigned int i = 0; i < this->n_source_centers; ++i)
  {
    const Tensor<1, dim> x_minus_xi = p - this->source_centers[i];

    // The first contribution is the Laplacian:
    return_value += ((2 * dim - 4 * x_minus_xi.norm_square() / (this->width * this->width)) /
                     (this->width * this->width) *
                     std::exp(-x_minus_xi.norm_square() / (this->width * this->width)));
  }

  return return_value / Utilities::fixed_power<dim>(std::sqrt(2 * pi) * this->width);
}



template<int dim>
class SolutionNegGrad : public Function<dim>, protected SolutionBase<dim>
{
public:
  SolutionNegGrad() : Function<dim>(dim)
  {
  }

  virtual ~SolutionNegGrad()
  {
  }

  virtual void
  vector_value(const Point<dim> & p, Vector<double> & v) const
  {
    AssertDimension(v.size(), dim);
    Tensor<1, dim> grad = solution.gradient(p);
    for(unsigned int d = 0; d < dim; ++d)
      v[d] = -grad[d];
  }

private:
  Solution<dim> solution;
};



template<int dim>
class MatrixIntegrator : public MeshWorker::LocalIntegrator<dim>
{
public:
  void
  cell(MeshWorker::DoFInfo<dim> & dinfo, typename MeshWorker::IntegrationInfo<dim> & info) const;
  void
  boundary(MeshWorker::DoFInfo<dim> &                  dinfo,
           typename MeshWorker::IntegrationInfo<dim> & info) const;
  void
  face(MeshWorker::DoFInfo<dim> &                  dinfo1,
       MeshWorker::DoFInfo<dim> &                  dinfo2,
       typename MeshWorker::IntegrationInfo<dim> & info1,
       typename MeshWorker::IntegrationInfo<dim> & info2) const;
};

template<int dim>
void
MatrixIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &                  dinfo,
                            typename MeshWorker::IntegrationInfo<dim> & info) const
{
  // return;
  LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0, false).matrix, info.fe_values());
}


// Interior faces use the interior penalty method
template<int dim>
void
MatrixIntegrator<dim>::face(MeshWorker::DoFInfo<dim> &                  dinfo1,
                            MeshWorker::DoFInfo<dim> &                  dinfo2,
                            typename MeshWorker::IntegrationInfo<dim> & info1,
                            typename MeshWorker::IntegrationInfo<dim> & info2) const
{
  // return;
  const unsigned int deg = info1.fe_values(0).get_fe().tensor_degree();
  LocalIntegrators::Laplace::ip_matrix(
    dinfo1.matrix(0, false).matrix,
    dinfo1.matrix(0, true).matrix,
    dinfo2.matrix(0, true).matrix,
    dinfo2.matrix(0, false).matrix,
    info1.fe_values(0),
    info2.fe_values(0),
    LocalIntegrators::Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));
}


template<int dim>
void
MatrixIntegrator<dim>::boundary(MeshWorker::DoFInfo<dim> &                  dinfo,
                                typename MeshWorker::IntegrationInfo<dim> & info) const
{
  // return;
  const unsigned int deg = info.fe_values(0).get_fe().tensor_degree();
  LocalIntegrators::Laplace::nitsche_matrix(
    dinfo.matrix(0, false).matrix,
    info.fe_values(0),
    LocalIntegrators::Laplace::compute_penalty(dinfo, dinfo, deg, deg));
}


template<int dim, int fe_degree>
class LaplaceProblem
{
public:
  typedef float Number;

  /*
   *  Constructor
   */
  LaplaceProblem(const bool use_dg);

  void
  run();

private:
  void
  setup_system();
  void
  assemble_system(const LaplaceOperator<dim, fe_degree, double> & laplace_operator);
  void
  solve();
  void
  output_results(const unsigned int cycle) const;

  parallel::distributed::Triangulation<dim> triangulation;
  std::shared_ptr<FiniteElement<dim>>       fe;
  MappingQGeneric<dim>                      mapping;
  DoFHandler<dim>                           dof_handler;
  ConstraintMatrix                          constraints;
  MatrixFree<dim, double>                   matrix_free;

  parallel::distributed::Vector<double> solution;
  parallel::distributed::Vector<double> solution_update;
  parallel::distributed::Vector<double> system_rhs;

  double                     setup_time;
  mutable ConditionalOStream pcout;
  ConditionalOStream         time_details;

  ConvergenceTable convergence_table;
};



template<int dim, int fe_degree>
LaplaceProblem<dim, fe_degree>::LaplaceProblem(const bool use_dg)
  : triangulation(MPI_COMM_WORLD,
                  Triangulation<dim>::none,
                  parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy),
    mapping(fe_degree),
    dof_handler(triangulation),
    setup_time(0.0),
    pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0),
    time_details(std::cout, true && Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
{
  if(!use_dg)
    fe.reset(new FE_Q<dim>(QGaussLobatto<1>(fe_degree + 1)));
  else
    fe.reset(new FE_DGQArbitraryNodes<dim>(QGaussLobatto<1>(fe_degree + 1)));
}


template<int dim, int fe_degree>
void
LaplaceProblem<dim, fe_degree>::setup_system()
{
  Timer time;
  time.start();
  setup_time = 0;

  dof_handler.distribute_dofs(*fe);
  dof_handler.distribute_mg_dofs(*fe);

  pcout << "Number of active cells:       " << triangulation.n_global_active_cells() << std::endl;
  pcout << "Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;

  IndexSet relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, relevant_dofs);
  constraints.clear();
  constraints.reinit(relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  VectorTools::interpolate_boundary_values(dof_handler,
                                           1,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  // constraints.add_line(0);
  constraints.close();

  setup_time += time.wall_time();
  time_details << "Distribute DoFs & B.C.     (CPU/wall) " << time() << "s/" << time.wall_time()
               << "s" << std::endl;
  time.restart();

  const QGauss<1>                                  quad(dof_handler.get_fe().degree + 1);
  typename MatrixFree<dim, double>::AdditionalData addit_data;
  addit_data.tasks_parallel_scheme = MatrixFree<dim, double>::AdditionalData::none;
  addit_data.build_face_info       = true;
  matrix_free.reinit(mapping, dof_handler, constraints, quad, addit_data);
  matrix_free.initialize_dof_vector(system_rhs);
  solution.reinit(system_rhs);
  solution_update.reinit(system_rhs);
  for(unsigned int i = 0; i < system_rhs.local_size(); ++i)
    system_rhs.local_element(i) = (double)rand() / RAND_MAX;

  setup_time += time.wall_time();
  time_details << "Setup matrix-free system   (CPU/wall) " << time() << "s/" << time.wall_time()
               << "s" << std::endl;
}



template<int dim, int fe_degree>
void
LaplaceProblem<dim, fe_degree>::assemble_system(
  const LaplaceOperator<dim, fe_degree, double> & laplace_operator)
{
  Timer                                     time;
  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(dof_handler, 1, Solution<dim>(), boundary_values);
  for(typename std::map<types::global_dof_index, double>::const_iterator it =
        boundary_values.begin();
      it != boundary_values.end();
      ++it)
    if(dof_handler.locally_owned_dofs().is_element(it->first))
      solution(it->first) = it->second;
  solution.update_ghost_values();

  system_rhs = 0;
  QGauss<dim>       quadrature_formula(fe->degree + 1);
  QGauss<dim - 1>   quadrature_face(fe->degree + 1);
  FEValues<dim>     fe_values(*fe,
                          quadrature_formula,
                          update_values | update_gradients | update_JxW_values |
                            update_quadrature_points);
  FEFaceValues<dim> fe_face_values(*fe,
                                   quadrature_face,
                                   update_values | update_gradients | update_quadrature_points |
                                     update_normal_vectors);

  const bool use_feq = fe->dofs_per_vertex > 0;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q_points    = quadrature_formula.size();

  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  RightHandSide<dim>                   right_hand_side;
  std::vector<double>                  rhs_values(n_q_points);
  Solution<dim>                        solution_function;
  std::vector<Tensor<1, dim>>          solution_gradients(n_q_points);
  Vector<double>                       local_rhs(dofs_per_cell);

  for(unsigned int c = 0; c < matrix_free.n_macro_cells(); ++c)
    for(unsigned int v = 0; v < matrix_free.n_components_filled(c); ++v)
    {
      typename DoFHandler<dim>::cell_iterator cell = matrix_free.get_cell_iterator(c, v);
      fe_values.reinit(cell);
      if(use_feq)
        fe_values.get_function_gradients(solution, solution_gradients);

      right_hand_side.value_list(fe_values.get_quadrature_points(), rhs_values);
      for(unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        double rhs_val = 0;
        if(use_feq)
          for(unsigned int q = 0; q < n_q_points; ++q)
            rhs_val += ((fe_values.shape_value(i, q) * rhs_values[q] -
                         fe_values.shape_grad(i, q) * solution_gradients[q]) *
                        fe_values.JxW(q));
        else
          for(unsigned int q = 0; q < n_q_points; ++q)
            rhs_val += ((fe_values.shape_value(i, q) * rhs_values[q]) * fe_values.JxW(q));
        local_rhs(i) = rhs_val;
      }
      if(!use_feq)
        for(unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
          if(cell->at_boundary(face))
          {
            fe_face_values.reinit(cell, face);
            const double sigmaF = laplace_operator.get_penalty_factor() *
                                  laplace_operator.get_array_penalty_parameter()[c][v] * 2.;
            for(unsigned int q = 0; q < quadrature_face.size(); ++q)
            {
              const double solution_value =
                solution_function.value(fe_face_values.quadrature_point(q));
              for(unsigned int i = 0; i < dofs_per_cell; ++i)
                local_rhs(i) +=
                  (sigmaF * fe_face_values.shape_value(i, q) -
                   fe_face_values.shape_grad(i, q) * fe_face_values.normal_vector(q)) *
                  solution_value * fe_face_values.JxW(q);
            }
          }

      cell->get_dof_indices(local_dof_indices);
      constraints.distribute_local_to_global(local_rhs, local_dof_indices, system_rhs);
    }
  system_rhs.compress(VectorOperation::add);

  setup_time += time.wall_time();
  time_details << "Assemble right hand side   (CPU/wall) " << time() << "s/" << time.wall_time()
               << "s" << std::endl;
}



template<int dim, int fe_degree>
void
LaplaceProblem<dim, fe_degree>::solve()
{
  Timer time;

  std::shared_ptr<BoundaryDescriptorLaplace<dim>> boundary_descriptor;
  boundary_descriptor.reset(new BoundaryDescriptorLaplace<dim>());

  std::shared_ptr<Function<dim>> zero_function;
  zero_function.reset(new Functions::ZeroFunction<dim>());

  boundary_descriptor->dirichlet.insert(
    std::pair<types::boundary_id, std::shared_ptr<Function<dim>>>(1, zero_function));
  boundary_descriptor->neumann.insert(
    std::pair<types::boundary_id, std::shared_ptr<Function<dim>>>(0, zero_function));

  LaplaceOperatorData<dim> laplace_operator_data;
  laplace_operator_data.bc             = boundary_descriptor;
  laplace_operator_data.penalty_factor = PENALTY_FACTOR;

  LaplaceOperator<dim, fe_degree, double> laplace_operator;
  laplace_operator.reinit(matrix_free, mapping, laplace_operator_data);

  MultigridData mg_data;
  mg_data.chebyshev_smoother_data.smoother_smoothing_range = 15;
  mg_data.chebyshev_smoother_data.smoother_poly_degree     = 4;
  mg_data.coarse_solver                                    = MultigridCoarseGridSolver::Chebyshev;

  typedef float                               Number;
  std::shared_ptr<PreconditionerBase<double>> preconditioner;

  typedef MyMultigridPreconditionerLaplace<dim,
                                           double,
                                           LaplaceOperator<dim, fe_degree, Number>,
                                           LaplaceOperatorData<dim>>
    MULTIGRID;
  preconditioner.reset(new MULTIGRID());

  std::shared_ptr<MULTIGRID> mg_preconditioner =
    std::dynamic_pointer_cast<MULTIGRID>(preconditioner);
  mg_preconditioner->initialize(
    mg_data, dof_handler, mapping, laplace_operator_data, laplace_operator_data.bc->dirichlet);

  // setup solver data
  CGSolverData solver_data;
  solver_data.solver_tolerance_rel = 1e-8;
  solver_data.use_preconditioner   = true;

  // setup solver
  CGSolver<LaplaceOperator<dim, fe_degree, double>,
           PreconditionerBase<double>,
           parallel::distributed::Vector<double>>
    solver(laplace_operator, *preconditioner, solver_data);

  setup_time += time.wall_time();
  pcout << "Initialize multigrid solver(CPU/wall) " << time() << "s/" << time.wall_time() << "s\n";

  assemble_system(laplace_operator);
  laplace_operator.apply_nullspace_projection(system_rhs);

  pcout << "Total setup time               (wall) " << setup_time << "s\n";
  pcout << "Number of multigrid levels: " << triangulation.n_global_levels() << std::endl;

  // measure time for applying the multigrid V-cycle
  for(unsigned int i = 0; i < N_REPETITIONS; ++i)
  {
    solution_update = 0;
    time.restart();

    preconditioner->vmult(solution_update, system_rhs);

    pcout << "Time V-cycle precondition  (CPU/wall) " << time() << "s/" << time.wall_time()
          << "s\n";
  }

  // measure time for applying the matrix-vector product
  for(unsigned int i = 0; i < N_REPETITIONS; ++i)
  {
    solution_update = 0;
    time.restart();

    laplace_operator.vmult(solution_update, system_rhs);

    pcout << "Time matrix-vector         (CPU/wall) " << time() << "s/" << time.wall_time()
          << "s\n";
  }

  // solve problem
  double sol_time = 1e10;
  for(unsigned int i = 0; i < N_REPETITIONS; ++i)
  {
    solution_update = 0;
    time.restart();

    unsigned int n_iter = solver.solve(solution_update, system_rhs);

    sol_time = std::min(sol_time, time.wall_time());
    pcout << "Time solve (" << n_iter << " iterations)  (CPU/wall) " << time() << "s/"
          << time.wall_time() << "s\n";
  }

  if(fe->dofs_per_vertex > 0)
    constraints.distribute(solution_update);
  solution += solution_update;

  solution.update_ghost_values();

  {
    Utilities::System::MemoryStats stats;
    Utilities::System::get_memory_stats(stats);
    Utilities::MPI::MinMaxAvg memory =
      Utilities::MPI::min_max_avg(stats.VmRSS / 1024., triangulation.get_communicator());
    pcout << "Memory stats [MB]: " << memory.min << " [p" << memory.min_index << "] " << memory.avg
          << " " << memory.max << " [p" << memory.max_index << "]" << std::endl;
  }

  // calculate L2 error
  Vector<float> difference_per_cell(triangulation.n_active_cells());
  VectorTools::integrate_difference(dof_handler,
                                    solution,
                                    Solution<dim>(),
                                    difference_per_cell,
                                    QGauss<dim>(fe->degree + 2),
                                    VectorTools::L2_norm);
  const double L2_error =
    std::sqrt(Utilities::MPI::sum(difference_per_cell.norm_sqr(), MPI_COMM_WORLD));

  pcout << "L2 error: " << std::scientific << std::setprecision(6) << L2_error << std::endl;

  // calculate H1 error
  VectorTools::integrate_difference(dof_handler,
                                    solution,
                                    Solution<dim>(),
                                    difference_per_cell,
                                    QGauss<dim>(fe->degree + 2),
                                    VectorTools::H1_seminorm);
  const double H1_error =
    std::sqrt(Utilities::MPI::sum(difference_per_cell.norm_sqr(), MPI_COMM_WORLD));

  pcout << "H1 error: " << std::scientific << std::setprecision(6) << H1_error << std::endl;

  const unsigned int level          = triangulation.n_global_levels() - 1;
  const unsigned int n_active_cells = triangulation.n_global_active_cells();
  const unsigned int n_dofs         = dof_handler.n_dofs();

  convergence_table.add_value("cycle", level);
  convergence_table.add_value("cells", n_active_cells);
  convergence_table.add_value("dofs", n_dofs);
  convergence_table.add_value("L2", L2_error);
  convergence_table.add_value("H1", H1_error);
  convergence_table.add_value("Sol time", sol_time);
}

template<int dim, int fe_degree>
void
LaplaceProblem<dim, fe_degree>::output_results(const unsigned int cycle) const
{
  return;
  DataOut<dim> data_out;

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches();

  std::ostringstream filename;
  filename << "solution-" << cycle << ".vtu";

  std::ofstream output(filename.str().c_str());
  data_out.write_vtu(output);
}

template<int dim, int fe_degree>
void
LaplaceProblem<dim, fe_degree>::run()
{
  pcout << "Number of MPI processes: " << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD)
        << std::endl;
  pcout << "Testing " << fe->get_name() << std::endl << std::endl;

  unsigned int n_cycles = 1;
  unsigned int sizes[]  = {1,  2,  3,  4,  5,  6,  7,  8,  10, 12, 14,  16,
                          20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128};
  if(RUN_VARIABLE_SIZES)
    n_cycles = sizeof(sizes) / sizeof(unsigned int);
  else
    n_cycles = N_REFINEMENTS;

  for(unsigned int cycle = 0; cycle < n_cycles; ++cycle)
  {
    pcout << "Cycle " << cycle << std::endl;

    if(RUN_VARIABLE_SIZES)
    {
      triangulation.clear();
      unsigned int n_refinements = 0;
      unsigned int n_subdiv      = sizes[cycle];
      if(n_subdiv > 1)
        while(n_subdiv % 2 == 0)
        {
          n_refinements += 1;
          n_subdiv /= 2;
        }
      if(dim == 2)
        n_refinements += 3;
      unsigned int njobs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
      while(njobs > 0)
      {
        njobs >>= dim;
        n_refinements++;
      }
      GridGenerator::subdivided_hyper_cube(triangulation, n_subdiv, -1, 1);
      triangulation.refine_global(n_refinements);
      if(DO_ADAPTIVE)
      {
        for(typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
            cell != triangulation.end();
            ++cell)
          if(cell->is_locally_owned() && cell->center().norm() < 0.55)
            cell->set_refine_flag();
        triangulation.execute_coarsening_and_refinement();
        for(typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
            cell != triangulation.end();
            ++cell)
          if(cell->is_locally_owned() && cell->center().norm() > 0.3 &&
             cell->center().norm() < 0.42)
            cell->set_refine_flag();
        triangulation.execute_coarsening_and_refinement();
        for(typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
            cell != triangulation.end();
            ++cell)
          if(cell->is_locally_owned() && cell->center().norm() > 0.335 &&
             cell->center().norm() < 0.39)
            cell->set_refine_flag();
        triangulation.execute_coarsening_and_refinement();
      }
    }
    else
    {
      if(cycle == 0)
      {
        GridGenerator::hyper_cube(triangulation, -1, 1.);
        // triangulation.refine_global (3-dim);
      }
      triangulation.refine_global(1);
      // triangulation.begin_active()->set_refine_flag();
      // triangulation.execute_coarsening_and_refinement();
    }

    // set boundary ID's
    for(typename Triangulation<dim>::cell_iterator cell = triangulation.begin();
        cell != triangulation.end();
        ++cell)
    {
      for(unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {
        if(cell->face(f)->at_boundary())
          cell->face(f)->set_all_boundary_ids(1);
      }
    }

    setup_system();
    solve();
    output_results(cycle);

    pcout << std::endl;
    if(dof_handler.n_dofs() >
       types::global_dof_index(4000000) * Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
      break;
  }

  convergence_table.set_precision("L2", 3);
  convergence_table.set_precision("H1", 3);

  convergence_table.set_scientific("L2", true);
  convergence_table.set_scientific("H1", true);

  convergence_table.evaluate_convergence_rates("L2",
                                               "cells",
                                               ConvergenceTable::reduction_rate_log2,
                                               dim);
  convergence_table.evaluate_convergence_rates("H1",
                                               "cells",
                                               ConvergenceTable::reduction_rate_log2,
                                               dim);

  if(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    convergence_table.write_text(std::cout);
}
} // namespace Step37

int
main(int argc, char ** argv)
{
  try
  {
    using namespace Step37;
    Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

    deallog.depth_console(0);
    bool do_dg = true;
    if(argc > 1)
      do_dg = atoi(argv[1]);

    LaplaceProblem<DIM, FE_DEGREE> laplace_problem(do_dg);
    laplace_problem.run();
  }
  catch(std::exception & exc)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------" << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------" << std::endl;
    return 1;
  }
  catch(...)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------" << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------" << std::endl;
    return 1;
  }

  return 0;
}
