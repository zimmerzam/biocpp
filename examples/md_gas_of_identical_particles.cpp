#include <BioCpp/molecular-dynamics/dynamics.h>
#include <BioCpp/molecular-dynamics/integrator.h>
#include <BioCpp/molecular-dynamics/phase_space.h>
#include <BioCpp/molecular-dynamics/potential.h>
#include <BioCpp/molecular-dynamics/thermostat.h>

#include <BioCpp.h>

#include <iostream>

class atom : public BioCpp::base_atom{
  public:
    double mass;
};
typedef BioCpp::pdb::model<atom>::type model;


class config_space : public phase_space<atom, 2, 1>{
  protected:
    typedef phase_space<atom, 2,1>::ext_vector_t ext_vector_t;  
  public:
    typedef phase_space<atom, 2,1>::vector_t vector_t;
    typedef phase_space<atom, 2,1>::matrix_t matrix_t;
    
    Eigen::Vector3d barycenter;
    
    config_space( model& sys ) : phase_space(sys) {
      barycenter = (system[1].mass*system[1].coordinate + system[2].mass*system[2].coordinate)/(system[1].mass+system[2].mass);
    };
    
    vector_t project( ext_vector_t& vector ){
      Eigen::Vector3d vec1 = vector.head(3);
      Eigen::Vector3d vec2 = vector.tail(3);
      
      // this vector stores the force applied to the atom-atom distance
      vector_t proj;
      proj << (vec2-vec1).dot(system[2].coordinate-system[1].coordinate);
      return proj;
    }
    
    void set( vector_t& coordinate ){
      double total_mass = system[1].mass+system[2].mass;
      Eigen::Vector3d dir = (system[2].coordinate - system[1].coordinate);
      dir.normalize();
      system[1].coordinate = barycenter - system[2].mass/total_mass*dir*coordinate(0);
      system[2].coordinate = barycenter + system[1].mass/total_mass*dir*coordinate(0);
      gen_coordinate = coordinate;
    }
};

struct params{
  double epsilon;
  double sigma;
};

int main(){
  atom at1;
  at1.serial=1;
  at1.id = BioCpp::atom::id::CA;
  at1.altLoc = ' ';
 	at1.resName = BioCpp::amino_acid::id::ALA;
 	at1.chainId = 'A';
 	at1.resSeq = 1;
 	at1.iCode = ' ';
 	at1.coordinate << 0.,0.,0.;
 	at1.occupancy = 0.;
 	at1.tempFactor = 0.;
 	at1.element = BioCpp::element::id::C;
 	at1.charge = 0.;
  at1.mass = 1.;
  
  atom at2;
  at2.serial=2;
  at2.id = BioCpp::atom::id::CA;
  at2.altLoc = ' ';
 	at2.resName = BioCpp::amino_acid::id::ALA;
 	at2.chainId = 'B';
 	at2.resSeq = 2;
 	at2.iCode = ' ';
 	at2.coordinate << 2.,4.,1.;
 	at2.occupancy = 0.;
 	at2.tempFactor = 0.;
 	at2.element = BioCpp::element::id::C;
 	at2.charge = 0.;
  at2.mass = 2.;
  
  model mdl;
  mdl.Append(1, at1);
  mdl.Append(2, at2);
  
  BioCpp::pdb::print_atom_line printer(std::cout);
  BioCpp::Iterate<atom>(mdl, printer);
 
  params param;
  param.epsilon = 10e-1;
  param.sigma = 1.;
  config_space space(mdl);
  dynamics dyn(0.,10e-2,10e-7);
  dyn.start( mdl, params,  );
  return 0;
}
