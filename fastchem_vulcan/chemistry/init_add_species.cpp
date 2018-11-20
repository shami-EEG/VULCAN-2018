
#include "fastchem.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include <cmath>



namespace fastchem {



template <class double_type>
void FastChem<double_type>::addAtom(std::string symbol)
{
  Element<double_type> species;

  species.symbol = symbol;
  species.element_index = getChemicalElementIndex(symbol);


  if (species.element_index == FASTCHEM_UNKNOWN_SPECIES)
    std::cout << "Element " << symbol << " from element abundance file not found in elements.dat. Neglected!\n";
  else
  {
    species.name = chemical_elements[species.element_index].name;
    species.molecular_weight = chemical_elements[species.element_index].atomic_weight;
    species.abundance = chemical_elements[species.element_index].abundance;
    elements.push_back(species);

    elements.back().index = elements.size()-1;
  }


}



template <class double_type>
void FastChem<double_type>::setElementAbundance(const std::string symbol, const double abundance)
{
  unsigned int index = getChemicalElementIndex(symbol);

  if (index == FASTCHEM_UNKNOWN_SPECIES)
    std::cout << "Element " << symbol << " for setting abundances not found. Neglected!\n";
  else
    chemical_elements[index].abundance = abundance;
}



template <class double_type>
void FastChem<double_type>::addMolecule(const std::string name, const std::string symbol,
                                        const std::vector<std::string> species_elements, const std::vector<int> stoichometric_coeff,
                                        const std::vector<double_type> mass_action_coeff, const int charge)
{
  Molecule<double_type> species;

  species.name = name;
  species.symbol = symbol;

  species.mass_action_coeff = mass_action_coeff;


  species.stoichometric_vector.assign(nb_elements, 0);


  bool is_stoichometry_complete = true;
  unsigned int nb_species_elements = 0;

  for (size_t i=0; i<species_elements.size(); ++i)
  {
    unsigned int index = getElementIndex(species_elements[i]);

    if (index == FASTCHEM_UNKNOWN_SPECIES)
      is_stoichometry_complete = false;
    else
    {
      species.stoichometric_vector[index] = stoichometric_coeff[i];
      species.element_indices.push_back(index);
    }

    nb_species_elements += stoichometric_coeff[i];
  }



  if (!is_stoichometry_complete)
    std::cout << "Stoichometry of species " << symbol << " incomplete. Neglected!\n";
  else
  {
    for(size_t j=0; j<nb_elements; ++j)
     species.sigma += species.stoichometric_vector[j];

    species.sigma = 1 - species.sigma;

    species.charge = charge;

    //definition of the molecular abundances
    species.abundance = 10.;

    for (size_t j=0; j<species.element_indices.size(); ++j)
      if (species.abundance > elements[species.element_indices[j]].abundance && elements[species.element_indices[j]].symbol != "e-")
        species.abundance = elements[species.element_indices[j]].abundance;

    //scaled abundances
    species.abundance_scaled = 10.;

    for (size_t j=0; j<species.element_indices.size(); ++j)
    {
      unsigned element_index = species.element_indices[j];

      if (species.abundance_scaled > elements[element_index].abundance/species.stoichometric_vector[element_index] && elements[species.element_indices[j]].symbol != "e-")
        species.abundance_scaled = elements[element_index].abundance/species.stoichometric_vector[element_index];
    }



    for (size_t j=0; j<species.element_indices.size(); ++j)
      species.molecular_weight += elements[species.element_indices[j]].molecular_weight * std::fabs(species.stoichometric_vector[species.element_indices[j]]);

    molecules.push_back(species);
  }


  //add the current molecule index to their respective elements
  if (is_stoichometry_complete)
    for (size_t i=0; i<species.element_indices.size(); ++i)
      elements[species.element_indices[i]].molecule_list.push_back(molecules.size()-1);
}




template class FastChem<double>;
template class FastChem<long double>;


}
