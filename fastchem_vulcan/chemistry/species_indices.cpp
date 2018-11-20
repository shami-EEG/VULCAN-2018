
#include "fastchem.h"

#include <string>
#include <vector>



namespace fastchem {


template <class double_type>
unsigned int FastChem<double_type>::getChemicalElementIndex(const std::string symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<nb_chemical_elements; ++i)
    if (symbol == chemical_elements[i].symbol)
    {
      index = i;
      break;
    }

  if (index == FASTCHEM_UNKNOWN_SPECIES && verbose_level >= 2) std::cout << "Warning, element " << symbol << " not found!\n";

  return index;
}



template <class double_type>
unsigned int FastChem<double_type>::getMoleculeIndex(const std::string symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<molecules.size(); ++i)
    if (symbol == molecules[i].symbol)
    {
      index = i;
      break;
    }

  if (index == FASTCHEM_UNKNOWN_SPECIES && verbose_level >= 2) std::cout << "Warning, molecule " << symbol << " not found!\n";

  return index;
}



template <class double_type>
unsigned int FastChem<double_type>::getElementIndex(const std::string symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<elements.size(); ++i)
    if (symbol == elements[i].symbol)
    {
      index = i;
      break;
    }

  if (index == FASTCHEM_UNKNOWN_SPECIES && verbose_level >= 2) std::cout << "Warning, element " << symbol << " not found!\n";

  return index;
}



template <class double_type>
unsigned int FastChem<double_type>::getSpeciesIndex(const std::string symbol)
{
  unsigned int index = FASTCHEM_UNKNOWN_SPECIES;

  for (size_t i=0; i<nb_species; ++i)
    if (symbol == species[i]->symbol)
    {
      index = i;
      break;
    }


  if (index == FASTCHEM_UNKNOWN_SPECIES && verbose_level >= 2) std::cout << "Warning, species " << symbol << " not found!\n";

  return index;
}



template <class double_type>
std::string FastChem<double_type>::getSpeciesName(const unsigned int species_index)
{
  if (species_index < nb_species)
    return species[species_index]->name;
  else
    return "";
}


template <class double_type>
std::string FastChem<double_type>::getSpeciesSymbol(const unsigned int species_index)
{

  if (species_index < nb_species)
    return species[species_index]->symbol;
  else
    return "";

}



template <class double_type>
std::string FastChem<double_type>::getElementName(const unsigned int species_index)
{
  if (species_index < nb_elements)
    return elements[species_index].name;
  else
    return "";
}



template <class double_type>
std::string FastChem<double_type>::getElementSymbol(const unsigned int species_index)
{

  if (species_index < nb_elements)
    return elements[species_index].symbol;
  else
    return "";

}



template class FastChem<double>;
template class FastChem<long double>;

}
