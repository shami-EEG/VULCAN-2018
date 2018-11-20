
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
void FastChem<double_type>::determineSolverOrder()
{

  for (unsigned int i=0; i<nb_elements; i++)
    elements[i].solver_order = determineSolverOrder(elements[i]);


}




template <class double_type>
unsigned int FastChem<double_type>::determineSolverOrder(const Element<double_type>& species)
{
  unsigned int solver_order = 0;

  for (size_t i=0; i<nb_molecules; ++i)
    for (size_t j=0; j<molecules[i].element_indices.size(); ++j)
    {
      if (molecules[i].element_indices[j] == species.index && species.symbol != "e-")
        if ( unsigned (molecules[i].stoichometric_vector[ molecules[i].element_indices[j] ]) > solver_order && molecules[i].abundance == species.abundance)
        {
          solver_order = molecules[i].stoichometric_vector[ molecules[i].element_indices[j] ];

          break;
        }

    }


  return solver_order;
}




template <class double_type>
void FastChem<double_type>::determineElementCalculationOrder()
{
  element_calculation_order.reserve(nb_elements);

  element_calculation_order.push_back(elements[0].index);


  for (size_t i=1; i<elements.size(); i++)
  {
    if (elements[i].abundance <= elements[element_calculation_order.back()].abundance)
    {
      element_calculation_order.push_back(elements[i].index);

      continue;
    }


    if (elements[i].abundance >= elements[element_calculation_order.front()].abundance)
    {
      element_calculation_order.insert(element_calculation_order.begin(), elements[i].index);

      continue;
    }


    for (std::vector<unsigned int>::iterator it = element_calculation_order.begin()+1; it<element_calculation_order.end(); it++)
    {
       if (elements[*(it-1)].abundance > elements[i].abundance && elements[*it].abundance <= elements[i].abundance)
       {
         element_calculation_order.insert(it, elements[i].index);

         break;
       }


    }

  }


}



template class FastChem<double>;
template class FastChem<long double>;


}
