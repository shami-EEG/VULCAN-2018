
#include "fastchem.h"

#include <vector>
#include <cmath>



namespace fastchem {


template <class double_type>
bool FastChem<double_type>::calcTotalHydrogenDensity(const double temperature_gas, const double pressure, const unsigned int grid_index,
                                                     double_type& h_density, double_type& density_iteration_lambda, double_type& density_iteration_error)
{
  //this value will be fixed
  double_type total_density = pressure/(CONST_K * temperature_gas);

  double_type total_density_calc = 0;

  for (size_t i=0; i<nb_species; ++i)
    total_density_calc += species[i]->number_density[grid_index];


  double_type current_total_density_error = (total_density - total_density_calc)/total_density;

  if (density_iteration_error * current_total_density_error < 0)
    density_iteration_lambda = 0.1 * density_iteration_lambda + 0.9;


  bool is_converged = false;

  if (std::fabs(current_total_density_error) < accuracy_delta)
     is_converged = true;
  else
  {
    if (total_density_calc > total_density)
      h_density = density_iteration_lambda * h_density;
    else
      h_density = 1./density_iteration_lambda * h_density;

    is_converged = false;
  }


  density_iteration_error = current_total_density_error;


  if (std::isnan(total_density_calc)) h_density = total_density_calc;


  return is_converged;
}


template class FastChem<double>;
template class FastChem<long double>;


}
