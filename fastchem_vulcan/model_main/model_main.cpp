
#include "../chemistry/fastchem.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <functional>

#include <cmath>


bool read_config_file(std::string file_path, std::string& fastchem_options_file, std::string& atmosphere_file,
                      std::string& chem_output_file, std::string& monitor_output_file)
{
  std::fstream file;

  file.open(file_path.c_str(), std::ios::in);

  if (file.fail())
  {
    std::cout << "Unable to read config file: " << file_path << "\n";
    return false;
  }


  std::string line;

  std::getline(file, line);

  std::getline(file, line);
  fastchem_options_file = line;

  std::getline(file, line);
  std::getline(file, line);

  std::getline(file, line);
  atmosphere_file = line;

  std::getline(file, line);
  std::getline(file, line);

  std::getline(file, line);
  chem_output_file = line;


  std::getline(file, line);
  std::getline(file, line);

  std::getline(file, line);
  monitor_output_file = line;

  /*std::cout << fastchem_options_file << "\n"
            << atmosphere_file << "\n"
            << chem_output_file << "\n"
            << monitor_output_file << "\n";*/

  return true;
}




int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    std::cout << "config file command line parameter missing!\n";

    return 1;
  }


  std::string config_file_name = argv[1];

  std::string fastchem_options_file;
  std::string atmosphere_file;
  std::string chem_output_file;
  std::string monitor_output_file;


  if (!read_config_file(config_file_name, fastchem_options_file,atmosphere_file,chem_output_file, monitor_output_file))
    return 1;



  fastchem::FastChem<long double> fastchem(fastchem_options_file, 1);


  std::vector<double> temperature;
  std::vector<double> pressure;


  //read the input data
  std::fstream file(atmosphere_file.c_str(), std::ios::in);

  std::string line;

  while(std::getline(file, line))
  {
    std::stringstream  line_stream(line);

    double temperature_in;
    double pressure_in;


    if (!(line_stream >> pressure_in >> temperature_in)) break;
    //if (!(line_stream >> temperature_in >> pressure_in)) break;


    pressure.push_back(pressure_in*1.e+6);  //pressure in dyn cm-2
    temperature.push_back(temperature_in);
  }

  file.close();




  unsigned int nb_grid_points = pressure.size();




  std::vector<double> total_density(pressure);

  for (unsigned int i=0; i<nb_grid_points; i++)
    total_density[i] /= fastchem::CONST_K * temperature[i];

  // for printing the t-p outout
  //for (unsigned int i=0; i<nb_grid_points; i++)
  //  std::cout << i << "\t" << pressure[i] << "\t" << temperature[i] << std::endl;


  std::vector< std::vector<double> > densities;
  densities.resize(nb_grid_points);

  std::vector<double> h_densities(nb_grid_points, 0.0);
  std::vector<double> mean_molecular_weights(nb_grid_points, 0.0);


  std::vector< std::vector<unsigned int> > element_conserved;
  element_conserved.resize(nb_grid_points);

  std::vector<unsigned int> nb_pressure_iterations (nb_grid_points, 0);
  std::vector<unsigned int> nb_chemistry_iterations (nb_grid_points, 0);
  std::vector<unsigned int> fastchem_flags (nb_grid_points, 0);


  //set terminal output of FastChem
  fastchem.setVerboseLevel(1);


  //call FastChem point-wise
  //for (unsigned int i=0; i<nb_grid_points; i++)
    //fastchem.calcDensities(temperature[i],pressure[i],densities[i],h_densities[i], mean_molecular_weights[i]);
    //fastchem_flags[i] = fastchem.calcDensities(temperature[i],pressure[i],densities[i],h_densities[i], mean_molecular_weights[i],
    //                                           element_conserved[i], nb_iterations[i]);


  //call FastChem for the entire temperature-pressure profile
  fastchem.calcDensities(temperature, pressure, densities, h_densities, mean_molecular_weights,
                         element_conserved, fastchem_flags, nb_pressure_iterations, nb_chemistry_iterations);


  //use OpenMP for acceleration
  /*unsigned int nb_omp_threads = 6;
  omp_set_num_threads(nb_omp_threads);

  std::vector< fastchem::FastChem<long double>  > fastchems(nb_omp_threads, fastchem);

  #pragma omp parallel for schedule(dynamic, 1)
  for (unsigned int i=0; i<nb_grid_points; i++)
  {
    std::cout << i << "\n";
    fastchem_flags[i] = fastchems[omp_get_thread_num()].calcDensities(temperature[i],pressure[i],densities[i],h_densities[i], mean_molecular_weights[i],
                                                                      element_conserved[i], nb_pressure_iterations[i], nb_chemistry_iterations[i]);

  }*/



  file.open(chem_output_file.c_str(), std::ios::out);

  unsigned int nb_species = fastchem.getSpeciesNumber();


  file << std::setw(16) << std::left << "P" << "\t"
       << std::setw(16) << std::left << "T" << "\t"
       << std::setw(16) << std::left << "m";
  for (unsigned int i=0; i<nb_species; i++)
    file << "\t" << std::setw(16) << std::left << fastchem.getSpeciesSymbol(i);

  file << "\n";


  for (unsigned int i=0; i<nb_grid_points; i++)
  {
    file << std::setprecision(10) << std::scientific
         << pressure[i]/1.e6 << "\t" << temperature[i] << "\t" << mean_molecular_weights[i];

    for (unsigned int j=0; j<nb_species; j++)
      //file << "\t" << densities[i][j];
		file << "\t" << densities[i][j] /total_density[i];

    file << "\n";
  }

  file.close();



  file.open(monitor_output_file.c_str(), std::ios::out);

  unsigned int nb_elements = fastchem.getElementNumber();


  file << std::setw(10) << std::left << "#grid point" << "\t"
       << std::setw(10) << std::left << "p_iterations" << "\t"
       << std::setw(10) << std::left << "c_iterations" << "\t"
       << std::setw(16) << std::left << "p_convergence" << "\t"
       << std::setw(16) << std::left << "c_convergence" << "\t"
       << std::setw(16) << std::left << "P (bar)" << "\t"
       << std::setw(16) << std::left << "T(k)" << "\t"
       << std::setw(16) << std::left << "n_<H> (cm-3)" << "\t"
       << std::setw(16) << std::left << "n_g (cm-3)" << "\t"
       << std::setw(16) << std::left << "m(u)";
  for (unsigned int i=0; i<nb_elements; i++)
    file << "\t" << std::setw(5) << std::left << fastchem.getElementSymbol(i);

  file << "\n";


  std::vector<std::string> output_flags {"fail", "ok"};
  std::vector<std::string> convergence_flags {"yes", "p conv fail", "chem conv fail", "no conv", "init fail"};


  for (unsigned int i=0; i<nb_grid_points; i++)
  {
    std::string p_conv;

    if (fastchem_flags[i] == fastchem::FASTCHEM_NO_PRESSURE_CONVERGENCE || fastchem_flags[i] == fastchem::FASTCHEM_NO_CONVERGENCE)
      p_conv = output_flags[0];
    else
      p_conv = output_flags[1];


    std::string c_conv;

    if (fastchem_flags[i] == fastchem::FASTCHEM_NO_FASTCHEM_CONVERGENCE || fastchem_flags[i] == fastchem::FASTCHEM_NO_CONVERGENCE)
      c_conv = output_flags[0];
    else
      c_conv = output_flags[1];




    file << std::setw(10) << i << "\t"
         << std::setw(10) << nb_pressure_iterations[i] << "\t"
         << std::setw(10) << nb_chemistry_iterations[i] << "\t"
         << std::setw(16) << p_conv << "\t"
         << std::setw(16) << c_conv << "\t";


    file << std::setprecision(10) << std::scientific
         << pressure[i]/1.e6 << "\t" << temperature[i] << "\t" << h_densities[i] << "\t" << total_density[i] << "\t" << mean_molecular_weights[i];

    for (unsigned int j=0; j<nb_elements; j++)
      file << "\t" << std::setw(5) << output_flags[element_conserved[i][j]];

    file << "\n";
  }

  file.close();


  /*file.open("output/chem_species.dat", std::ios::out);


  for (unsigned int j=0; j<nb_species; j++)
    file << fastchem.getSpeciesSymbol(j) << "\t" << j << "\t" << densities[0][j] /total_density[0] << "\t" << densities[nb_grid_points-1][j] /total_density[nb_grid_points-1] << "\n";

  file.close();*/

  std::cout << "model finished: " << std::endl;

  return 0;
}
