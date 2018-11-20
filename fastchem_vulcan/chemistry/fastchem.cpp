
#include "fastchem.h"

#include <string>
#include <vector>



namespace fastchem {


template <class double_type>
FastChem<double_type>::FastChem(const std::string& model_parameter_file, const unsigned int verbose_level_init)
{
  verbose_level = verbose_level_init;

  bool parameter_file_loaded = false;

  if (model_parameter_file != "")
    parameter_file_loaded = readParameterFile(model_parameter_file);


  if (!parameter_file_loaded)
  {
    std::cout << "Error reading parameters\n";
    is_initialized = false;
  }


  if (parameter_file_loaded) init();
}



//Copy constructor
//Could be made more pretty, but this one does the job as well...
template <class double_type>
FastChem<double_type>::FastChem(const FastChem &obj)
{
  nb_chemical_elements = obj.nb_chemical_elements;
  nb_species = obj.nb_species;
  nb_molecules = obj.nb_molecules;
  nb_elements = obj.nb_elements;

  nb_max_fastchem_iter = obj.nb_max_fastchem_iter;
  nb_max_pressure_iter = obj.nb_max_pressure_iter;
  nb_max_bisection_iter = obj.nb_max_bisection_iter;
  nb_max_neldermead_iter = obj.nb_max_neldermead_iter;
  nb_max_newton_iter = obj.nb_max_newton_iter;

  element_density_minlimit = obj.element_density_minlimit;
  molecule_density_minlimit = obj.molecule_density_minlimit;

  accuracy = obj.accuracy;
  accuracy_delta = obj.accuracy_delta;
  newton_err = obj.newton_err;

  verbose_level = obj.verbose_level;
  use_scaling_factor = obj.use_scaling_factor;
  is_initialized = obj.is_initialized;


  chemical_element_file = obj.chemical_element_file;
  species_data_file = obj.species_data_file;
  element_abundances_file = obj.element_abundances_file;


  chemical_elements = obj.chemical_elements;
  elements = obj.elements;
  molecules = obj.molecules;

  element_calculation_order = obj.element_calculation_order;

  for (size_t i=0; i<nb_elements; ++i)
    species.push_back(&elements[i]);

  for (size_t i=0; i<nb_molecules; ++i)
    species.push_back(&molecules[i]);
}



template class FastChem<double>;
template class FastChem<long double>;


}
