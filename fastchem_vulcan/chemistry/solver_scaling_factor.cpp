
#include "fastchem.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>


namespace fastchem {


//computes the scaling factor psi to avoid numerical overflow
template <class double_type>
double_type FastChem<double_type>::solverScalingFactor(Element<double_type>& species, const double_type number_density_min, const double_type h_density, const unsigned int grid_index)
{
  double_type scaling_factor = 0.0;


  for (size_t j=0; j<species.molecule_list.size(); ++j)
  {
    unsigned int i = species.molecule_list[j];


    if (molecules[i].stoichometric_vector[species.index] < 1 || molecules[i].stoichometric_vector[species.index] > species.solver_order)
      continue;


    if (molecules[i].abundance == species.abundance)
    {
      molecules[i].sum[grid_index] = 0.0;


      for (size_t k=0; k<molecules[i].element_indices.size(); ++k)
      {
        unsigned int l = molecules[i].element_indices[k];

        if (l != species.index)
          molecules[i].sum[grid_index] += molecules[i].stoichometric_vector[l] * std::log(elements[l].number_density[grid_index]);

      }


      molecules[i].sum[grid_index] += molecules[i].mass_action_constant[grid_index];
    }


    if (molecules[i].sum[grid_index] > scaling_factor)
      scaling_factor = molecules[i].sum[grid_index];

  }



  double_type xi = (number_density_min + species.number_density[grid_index]) * std::exp(-scaling_factor);

  for (size_t i=0; i<species.molecule_list.size(); ++i)
  {
    unsigned int j = species.molecule_list[i];

    if (species.abundance == molecules[j].abundance)
      xi += molecules[j].stoichometric_vector[species.index] * std::pow(molecules[j].number_density[grid_index], molecules[j].stoichometric_vector[species.index]);
  }


  if (xi == 0)
    xi = std::numeric_limits<double_type>::max_exponent / 6.0;
  else
    xi = std::numeric_limits<double_type>::max_exponent - std::log(xi);

  xi = std::sqrt(xi);


  scaling_factor -= xi;


  return scaling_factor;
}




template class FastChem<double>;
template class FastChem<long double>;

}



