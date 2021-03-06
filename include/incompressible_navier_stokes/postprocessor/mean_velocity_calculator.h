/*
 * mean_velocity_calculator.h
 *
 *  Created on: Nov 20, 2017
 *      Author: fehn
 */

#ifndef INCLUDE_INCOMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_MEAN_VELOCITY_CALCULATOR_H_
#define INCLUDE_INCOMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_MEAN_VELOCITY_CALCULATOR_H_

#include <deal.II/matrix_free/fe_evaluation_notemplate.h>


#include "../../functionalities/print_functions.h"

namespace IncNS
{
template<int dim>
struct MeanVelocityCalculatorData
{
  MeanVelocityCalculatorData()
    : calculate(false),
      write_to_file(false),
      direction(Tensor<1, dim, double>()),
      filename_prefix("mean_velocity")
  {
  }

  void
  print(ConditionalOStream & pcout)
  {
    if(calculate == true)
    {
      pcout << "  Mean velocity/flow rate calculator:" << std::endl;

      print_parameter(pcout, "Calculate mean velocity/flow rate", calculate);
      print_parameter(pcout, "Write results to file", write_to_file);
      if(write_to_file == true)
        print_parameter(pcout, "Filename", filename_prefix);
    }
  }

  // calculate mean velocity?
  bool calculate;

  // Set containing boundary ID's of the surface area
  // for which we want to calculate the mean velocity.
  // This parameter is only relevant for area-based computation.
  std::set<types::boundary_id> boundary_IDs;

  // write results to file?
  bool write_to_file;

  // Direction in which we want to compute the flow rate
  // This parameter is only relevant for volume-based computation.
  Tensor<1, dim, double> direction;

  // filename
  std::string filename_prefix;
};

template<int dim, typename Number>
class MeanVelocityCalculator
{
public:
  typedef LinearAlgebra::distributed::Vector<Number> VectorType;

  typedef CellIntegrator<dim, dim, Number> CellIntegratorU;
  typedef FaceIntegrator<dim, dim, Number> FaceIntegratorU;

  typedef MeanVelocityCalculator<dim, Number> This;

  typedef VectorizedArray<Number> scalar;

  MeanVelocityCalculator(MatrixFree<dim, Number> const &         matrix_free_in,
                         unsigned int const                      dof_index_in,
                         unsigned int const                      quad_index_in,
                         MeanVelocityCalculatorData<dim> const & data_in);

  /*
   * Calculates the mean velocity through a given cross section of the domain by dividing
   * the flow rate through the cross section area. This function is more general than
   * calculate_mean_velocity_volume() and can be used for domains with varying cross-section
   * area in streamwise direction.
   */
  Number
  calculate_mean_velocity_area(VectorType const & velocity, double const & time);

  /*
   * Calculate mean velocity (only makes sense if the domain has a constant cross-section area in
   * streamwise direction.
   */
  Number
  calculate_mean_velocity_volume(VectorType const & velocity, double const & time);

  /*
   * Calculate flow rate in m^3/s, for example for problems with non-constant cross-section area. To
   * obtain the flow rate, the length of the domain in streamwise direction has to be specified.
   */
  Number
  calculate_flow_rate_volume(VectorType const & velocity,
                             double const &     time,
                             double const &     length) const;

  /*
   * Calculates the flow rate through a given cross section of the domain.
   */
  Number
  calculate_flow_rate_area(VectorType const & velocity, double const & time) const;


private:
  void
  write_output(Number const & value, double const & time, std::string const & name) const;

  Number
  calculate_area() const;

  Number
  calculate_volume() const;

  void
  local_calculate_volume(MatrixFree<dim, Number> const & data,
                         std::vector<Number> &           dst,
                         VectorType const &,
                         std::pair<unsigned int, unsigned int> const & cell_range) const;

  Number
  do_calculate_flow_rate_area(VectorType const & velocity) const;

  Number
  do_calculate_mean_velocity_volume(VectorType const & velocity) const;

  Number
  do_calculate_flow_rate_volume(VectorType const & velocity) const;

  void
  local_calculate_flow_rate_volume(MatrixFree<dim, Number> const &               data,
                                   std::vector<Number> &                         dst,
                                   VectorType const &                            src,
                                   std::pair<unsigned int, unsigned int> const & cell_range) const;

  MeanVelocityCalculatorData<dim> const & data;
  MatrixFree<dim, Number> const &         matrix_free;
  unsigned int                            dof_index, quad_index;
  bool                                    area_has_been_initialized, volume_has_been_initialized;
  double                                  area, volume;
  mutable bool                            clear_files;
};

} // namespace IncNS

#endif /* INCLUDE_INCOMPRESSIBLE_NAVIER_STOKES_POSTPROCESSOR_MEAN_VELOCITY_CALCULATOR_H_ */
