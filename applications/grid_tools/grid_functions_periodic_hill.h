/*
 * grid_functions_periodic_hill.h
 *
 *  Created on: Mar 29, 2019
 *      Author: fehn
 */

#ifndef APPLICATIONS_GRID_TOOLS_GRID_FUNCTIONS_PERIODIC_HILL_H_
#define APPLICATIONS_GRID_TOOLS_GRID_FUNCTIONS_PERIODIC_HILL_H_

double
m_to_mm(double coordinate)
{
  return 1000.0 * coordinate;
}

double
mm_to_m(double coordinate)
{
  return 0.001 * coordinate;
}

/*
 * This function returns the distance by which points are shifted in y-direction due to the hill,
 * i.e., a value of 0 is returned at x=0H, 9H, and a value of -H at x=4.5H.
 */
double
f(double x_m)
{
  if(x_m > LENGTH / 2.0)
    x_m = LENGTH - x_m;

  double x = m_to_mm(x_m);
  double y = 0.0;

  AssertThrow(x_m <= LENGTH / 2.0 + 1.e-12, ExcMessage("Parameter out of bounds."));

  if(x <= 9.0)
    y =
      -m_to_mm(H) +
      std::min(m_to_mm(H), m_to_mm(H) + 6.775070969851e-3 * x * x - 2.124527775800e-3 * x * x * x);
  else if(x > 9.0 && x <= 14.0)
    y = -m_to_mm(H) + 2.507355893131e1 + 9.754803562315e-1 * x - 1.016116352781e-1 * x * x +
        1.889794677828e-3 * x * x * x;
  else if(x > 14.0 && x <= 20.0)
    y = -m_to_mm(H) + 2.579601052357e1 + 8.206693007457e-1 * x - 9.055370274339e-2 * x * x +
        1.626510569859e-3 * x * x * x;
  else if(x > 20.0 && x <= 30.0)
    y = -m_to_mm(H) + 4.046435022819e1 - 1.379581654948 * x + 1.945884504128e-2 * x * x -
        2.070318932190e-4 * x * x * x;
  else if(x > 30.0 && x <= 40.0)
    y = -m_to_mm(H) + 1.792461334664e1 + 8.743920332081e-1 * x - 5.567361123058e-2 * x * x +
        6.277731764683e-4 * x * x * x;
  else if(x > 40.0 && x <= 54.0)
    y = -m_to_mm(H) + std::max(0.0,
                               5.639011190988e1 - 2.010520359035 * x + 1.644919857549e-2 * x * x +
                                 2.674976141766e-5 * x * x * x);
  else if(x > 54.0)
    y = -m_to_mm(H);
  else
    AssertThrow(false, ExcMessage("Not implemented."));

  return mm_to_m(y);
}

/*
 * Grid stretch factor gamma as a function of x:
 *  - manipulates (increases) the grid stretch factor where the friction coefficient reaches a
 * maximum (between 8.3H < x <= 9H), and also between 6.3H < x <= 8.3H
 */
double
gamma_x(double x_m)
{
  (void)x_m;

  double gamma = 1.0;

  // TODO
  //  double xmod = x_m/H;
  //  if (xmod > 6.3 && xmod <= 8.3)
  //  {
  //    xmod = 8.3 - xmod;
  //    gamma = (-0.02*std::cos(xmod*numbers::PI)+1.02);
  //  }
  //  else if (xmod > 8.3) //move mesh closer to the wall to get a lower y+ value at the peak
  //  {
  //    xmod = 9.0 - xmod;
  //    gamma = (-0.05*std::cos(xmod*numbers::PI*2./0.7)+1.05);
  //  }

  return GRID_STRETCH_FAC * gamma;
}

template<int dim>
class PeriodicHillManifold : public ChartManifold<dim>
{
public:
  PeriodicHillManifold() : ChartManifold<dim>()
  {
  }

  Point<dim>
  push_forward(Point<dim> const & xi) const
  {
    Point<dim> x = xi;

    // transform y-coordinate only
    double const gamma           = gamma_x(xi[0]);
    double const xi_1_normalized = (xi[1] - H) / HEIGHT;
    double       xi_1_hat = std::tanh(gamma * (2.0 * xi_1_normalized - 1.0)) / std::tanh(gamma);
    x[1]                  = H + (xi_1_hat + 1.0) / 2.0 * HEIGHT + (1.0 - xi_1_hat) / 2.0 * f(xi[0]);

    return x;
  }

  Point<dim>
  pull_back(Point<dim> const & x) const
  {
    Point<dim> xi = x;

    // transform y-coordinate only
    double const f_x      = f(x[0]);
    double const xi_1_hat = (2.0 * x[1] - 2.0 * H - HEIGHT - f_x) / (HEIGHT - f_x);
    double const gamma    = gamma_x(x[0]);
    xi[1] = HEIGHT / 2.0 * (std::atanh(xi_1_hat * std::tanh(gamma)) / gamma + 1) + H;

    return xi;
  }

  std::unique_ptr<Manifold<dim>>
  clone() const override
  {
    return std_cxx14::make_unique<PeriodicHillManifold<dim>>();
  }
};

#endif /* APPLICATIONS_GRID_TOOLS_GRID_FUNCTIONS_PERIODIC_HILL_H_ */
