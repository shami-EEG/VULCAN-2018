
#include "fastchem.h"

#include <vector>
#include <cmath>



namespace fastchem {


template <class double_type>
double FastChem<double_type>::calcMeanMolecularWeight(const double total_density, const unsigned int grid_index)
{
   double mean_molecular_weight = 0.0;

   for(size_t i=0; i<nb_species; ++i)
     mean_molecular_weight += species[i]->molecular_weight * species[i]->number_density[grid_index];

   mean_molecular_weight /= total_density;


   return mean_molecular_weight;
}



template class FastChem<double>;
template class FastChem<long double>;


}
