/*!
    \file statistics_o_position.cpp
*/

#include <iostream>
#include <BioCpp/Core>
#include <BioCpp/Standard>

typedef BioCpp::base::atom atom;
typedef BioCpp::standard::base::residue residue;
typedef BioCpp::standard::base::chain chain;
typedef BioCpp::standard::base::complex complex;
typedef BioCpp::standard::base::complex_constructor complex_constructor;
typedef BioCpp::io::model<atom>::type model;

struct positions{
  Eigen::Vector3d pca;
  Eigen::Vector3d pc;
  Eigen::Vector3d n;
  Eigen::Vector3d h;
	Eigen::Vector3d ca;
	Eigen::Vector3d ha;
	Eigen::Vector3d cb;
	Eigen::Vector3d c;
	Eigen::Vector3d o;
	Eigen::Vector3d sn;
	Eigen::Vector3d sh;
	Eigen::Vector3d sca;
};

int main( int argc, char* argv[] ){
  std::list<positions> position_list;
  
	for( int file_num = 1; file_num!=argc; ++file_num ){
	  const char* filename = argv[file_num];
  	BioCpp::io::pdb::file PDB(filename, 0); // read the pdb file. 
  	bool is_good_structure = true;
  	
  	// gets the structure, align the sequence from ATOM lines to the one from SEQRES and builds the topology
		// if SEQRES has not been found and the sequence in ATOM has at least one hole raise an error. if the
		// sequence from ATOM lines has no holes no error will be raised and the topology will be built.
		model mdl = PDB.readModel<atom>(1); 
		BioCpp::error err = PDB.error;
  	BioCpp::warning war = PDB.warning;
		if(PDB.TseqRes.size() != 0){
  		for( BioCpp::io::seqres_record::iterator seqres = PDB.TseqRes.begin(); seqres!=PDB.TseqRes.end(); ++seqres ){
 	    	if( PDB.RseqRes.find(seqres->first)!=PDB.RseqRes.end() and seqres->second!=PDB.RseqRes[seqres->first]){
      		PDB.RseqRes[seqres->first].insert( PDB.RseqRes[seqres->first].begin(), '-' );
      		PDB.RseqRes[seqres->first].insert( PDB.RseqRes[seqres->first].end(), '-' );
        	double score = BioCpp::fasta::StrictNeedlemanWunsch(seqres->second, PDB.RseqRes[seqres->first], BioCpp::fasta::ZIMM1, err, war);
        	if(score < 0){
          	is_good_structure = false;
        	}
      	}
    	}
  	}
  	else{
   		for( BioCpp::io::seqres_record::iterator seqres = PDB.RseqRes.begin(); seqres!=PDB.RseqRes.end(); ++seqres ){
      	for( std::string::iterator str = seqres->second.begin(); str != seqres->second.end(); ++str ){
        	if((*str)=='-'){
          	is_good_structure = false;
        	}
      	}
    	}
  	}
  	if( not is_good_structure ){
  		break;
  	}
  	// Order atoms in 'all_info' by using informations stored in the pdb: missing residues are 
  	//detected and stored as empty residues 
  	BioCpp::standard::base::complex_constructor constr(BioCpp::residue::dictionary);
  	BioCpp::standard::base::complex cmp = constr(mdl, PDB.RseqRes, PDB.TseqRes);
  	
  	for( complex::iterator ch = cmp.begin(); ch!=cmp.end(); ++ch ){
  		chain::iterator last = ch->end(); --last;
  		chain::iterator secn = ch->begin(); ++secn;
  		for( chain::iterator res = secn; res!=last; ++res ){
  			chain::iterator pres = res; --pres;
  			chain::iterator sres = res; ++sres;
  			if( pres->exists( BioCpp::atom::CA ) and pres->exists( BioCpp::atom::C )
  			    and res->exists( BioCpp::atom::N ) and res->exists( BioCpp::atom::CA ) and res->exists( BioCpp::atom::C )  
  			    and res->exists( BioCpp::atom::H ) and res->exists( BioCpp::atom::CB ) and res->exists( BioCpp::atom::O )
  			    and res->exists( BioCpp::atom::HA )
  			    and sres->exists( BioCpp::atom::N ) and sres->exists( BioCpp::atom::H ) and sres->exists( BioCpp::atom::CA ) ){
  				positions pos;
  				pos.pca = (*pres)[BioCpp::atom::CA].coordinate;
  				pos.pc =  (*pres)[BioCpp::atom::C].coordinate;
  				pos.n = (*res)[BioCpp::atom::N].coordinate;
  				pos.ca = (*res)[BioCpp::atom::CA].coordinate;
  				pos.ha = (*res)[BioCpp::atom::HA].coordinate;
  				pos.c =  (*res)[BioCpp::atom::C].coordinate;
  				pos.h =  (*res)[BioCpp::atom::H].coordinate;
  				pos.cb =  (*res)[BioCpp::atom::CB].coordinate;
  				pos.o =  (*res)[BioCpp::atom::O].coordinate;
  				pos.sn = (*sres)[BioCpp::atom::N].coordinate;
  				pos.sh = (*sres)[BioCpp::atom::H].coordinate;
  				pos.sca = (*sres)[BioCpp::atom::CA].coordinate;
  				position_list.push_back( pos );
  			} 
  		}
  	}
  }
  
  // First output: positions of H, CB, O are in local coordinate system
  // For instance the coordinate system for O is centered in C and is defined by
  // x=ca_c+c_n, y = ca_c-c_n, z = c_n.cross(ca_c).
  std::cout << "###############################################" << std::endl;
  std::cout << "Now using a local coordinate frame" << std::endl;
  std::cout << "###############################################" << std::endl;
  positions average1;
  std::list<positions> frame1_list;
  average1.h = Eigen::Vector3d::Zero();
  average1.cb = Eigen::Vector3d::Zero();
  average1.o = Eigen::Vector3d::Zero();
  average1.ha = Eigen::Vector3d::Zero();
  for(std::list<positions>::iterator p = position_list.begin(); p != position_list.end(); ++p){
    positions pos;
  	// H position
  	{
    	Eigen::Vector3d pc_n = p->n - p->pc;
    	Eigen::Vector3d n_ca = p->ca - p->n;
    	Eigen::Vector3d n_h = p->h - p->n;
    	pc_n.normalize();
    	n_ca.normalize();
    	Eigen::Vector3d x = pc_n+n_ca;
    	Eigen::Vector3d y = pc_n-n_ca;
    	Eigen::Vector3d z = x.cross(y);
    	x.normalize();
    	y.normalize();
    	z.normalize();
    	pos.h = Eigen::Vector3d(n_h.dot(x),n_h.dot(y),n_h.dot(z) );
      average1.h+=pos.h;
    }
    // HA position
  	{
    	Eigen::Vector3d n_ca = p->ca - p->n;
    	Eigen::Vector3d ca_c = p->c - p->ca;
    	Eigen::Vector3d ca_ha = p->ha - p->ca;
    	n_ca.normalize();
    	ca_c.normalize();
    	Eigen::Vector3d x = n_ca+ca_c;
    	Eigen::Vector3d y = n_ca-ca_c;
    	Eigen::Vector3d z = x.cross(y);
    	x.normalize();
    	y.normalize();
    	z.normalize();
    	pos.ha = Eigen::Vector3d(ca_ha.dot(x),ca_ha.dot(y),ca_ha.dot(z) );
      average1.ha+=pos.ha;
    }
  	// CB position
  	{
    	Eigen::Vector3d n_ca = p->ca - p->n;
    	Eigen::Vector3d ca_c = p->c - p->ca;
    	Eigen::Vector3d ca_cb = p->cb - p->ca;
    	n_ca.normalize();
    	ca_c.normalize();
    	Eigen::Vector3d x = n_ca+ca_c;
    	Eigen::Vector3d y = n_ca-ca_c;
    	Eigen::Vector3d z = x.cross(y);
    	x.normalize();
    	y.normalize();
    	z.normalize();
    	pos.cb = Eigen::Vector3d(ca_cb.dot(x),ca_cb.dot(y),ca_cb.dot(z) );
      average1.cb+=pos.cb;
    }
  	// O position
  	{
    	Eigen::Vector3d ca_c = p->c - p->ca;
    	Eigen::Vector3d c_n = p->n - p->c;
    	Eigen::Vector3d c_o = p->o - p->c;
    	ca_c.normalize();
    	c_n.normalize();
    	Eigen::Vector3d x = ca_c+c_n;
    	Eigen::Vector3d y = ca_c-c_n;
    	Eigen::Vector3d z = x.cross(y);
    	x.normalize();
    	y.normalize();
    	z.normalize();
    	pos.o = Eigen::Vector3d(c_o.dot(x),c_o.dot(y),c_o.dot(z) );
      average1.o += pos.o;
    }
    frame1_list.push_back(pos);
  }
  average1.h /= position_list.size();
  average1.cb /= position_list.size();
  average1.ha /= position_list.size();
  average1.o /= position_list.size();
  double var_h = 0., var_cb = 0., var_ha = 0., var_o = 0.;
  for(std::list<positions>::iterator p = frame1_list.begin(); p != frame1_list.end(); ++p){
    var_h += (p->h-average1.h).squaredNorm();
    var_cb += (p->cb-average1.cb).squaredNorm();
    var_ha += (p->ha-average1.ha).squaredNorm();
    var_o += (p->o-average1.o).squaredNorm();
  }
  var_h = sqrt(var_h/frame1_list.size());
  var_cb = sqrt(var_cb/frame1_list.size());
  var_ha = sqrt(var_ha/frame1_list.size());
  var_o = sqrt(var_o/frame1_list.size());
  std::cout << "BACKBONE H average position: " << std::endl
            << average1.h.transpose() << "  " << var_h << std::endl;
  std::cout << "BACKBONE HA average position: " << std::endl
            << average1.ha.transpose() << "  " << var_ha << std::endl;
  std::cout << "BACKBONE CB average position: " << std::endl
            << average1.cb.transpose() << "  " << var_cb << std::endl;
  std::cout << "BACKBONE O average position: " << std::endl
            << average1.o.transpose() << "  " << var_o << std::endl << std::endl;
  // Second output: positions of H, CB, C are in a local coordinate system based on CA coordinate

  std::cout << "###############################################" << std::endl;
  std::cout << "Now using a CA-based coordinate frame" << std::endl;
  std::cout << "###############################################" << std::endl;
  positions average2;
  std::list<positions> frame2_list;
  average2.h = Eigen::Vector3d::Zero();
  average2.cb = Eigen::Vector3d::Zero();
  average2.o = Eigen::Vector3d::Zero();
  average2.ha = Eigen::Vector3d::Zero();
  for(std::list<positions>::iterator p = position_list.begin(); p != position_list.end(); ++p){
    positions pos;
   	Eigen::Vector3d pca_ca = p->ca - p->pca;
   	Eigen::Vector3d ca_sca = p->sca - p->ca;
   	Eigen::Vector3d x = pca_ca+ca_sca;
   	Eigen::Vector3d y = pca_ca-ca_sca;
   	Eigen::Vector3d z = x.cross(y);
   	x.normalize();
   	y.normalize();
   	z.normalize();
  	// H position
  	Eigen::Vector3d ca_h = p->h - p->ca;
  	pos.h = Eigen::Vector3d(ca_h.dot(x),ca_h.dot(y),ca_h.dot(z) );
  	average2.h+=pos.h;
  	// HA position
  	Eigen::Vector3d ca_ha = p->ha - p->ca;
  	pos.ha = Eigen::Vector3d(ca_ha.dot(x),ca_ha.dot(y),ca_ha.dot(z) );
  	average2.ha+=pos.ha;
  	// CB position
  	Eigen::Vector3d ca_cb = p->cb - p->ca;
  	pos.cb = Eigen::Vector3d(ca_cb.dot(x),ca_cb.dot(y),ca_cb.dot(z) );
  	average2.cb+=pos.cb;
  	// CO position
  	Eigen::Vector3d ca_o = p->o - p->ca;
  	pos.o = Eigen::Vector3d(ca_o.dot(x),ca_o.dot(y),ca_o.dot(z) );
  	average2.o+=pos.o;
  	frame2_list.push_back(pos);
  }
  average2.h /= position_list.size();
  average2.cb /= position_list.size();
  average2.ha /= position_list.size();
  average2.o /= position_list.size();
  double var2_h = 0., var2_ha = 0., var2_cb = 0., var2_o = 0.;
  for(std::list<positions>::iterator p = frame2_list.begin(); p != frame2_list.end(); ++p){
    var2_h += (p->h-average2.h).squaredNorm();
    var2_cb += (p->cb-average2.cb).squaredNorm();
    var2_ha += (p->ha-average2.ha).squaredNorm();
    var2_o += (p->o-average2.o).squaredNorm();
  }
  var2_h = sqrt(var2_h/frame2_list.size());
  var2_cb = sqrt(var2_cb/frame2_list.size());
  var2_ha = sqrt(var2_ha/frame2_list.size());
  var2_o = sqrt(var2_o/frame2_list.size());
  std::cout << "BACKBONE H average position: " << std::endl
            << average2.h.transpose() << "  " << var2_h << std::endl;
  std::cout << "BACKBONE HA average position: " << std::endl
            << average2.ha.transpose() << "  " << var2_ha << std::endl;
  std::cout << "BACKBONE CB average position: " << std::endl
            << average2.cb.transpose() << "  " << var2_cb << std::endl;
  std::cout << "BACKBONE O average position: " << std::endl
            << average2.o.transpose() << "  " << var2_o << std::endl << std::endl;
  
  // Third output: angle between co and nh
  std::cout << "###############################################" << std::endl;
  std::cout << "Peptide bond planarity" << std::endl;
  std::cout << "###############################################" << std::endl;
  double average_angle = 0., dev_angle = 0.;
  double average_dist_o = 0., dev_dist_o = 0.;
  std::list<double> angles_list, dist_o_list;
  for(std::list<positions>::iterator p = position_list.begin(); p != position_list.end(); ++p){
    Eigen::Vector3d co = p->o-p->c;
    Eigen::Vector3d nh = p->sh-p->sn;
    co.normalize();
    nh.normalize();
    double angle = acos(co.dot(nh));
    angles_list.push_back(angle);
    average_angle+=angle;
    
    Eigen::Vector3d o = p->c - nh*(p->o-p->c).norm();
    double dist_o = (o-p->o).norm();
    dist_o_list.push_back(dist_o);
    average_dist_o += dist_o;
  }
  average_angle/=position_list.size();
  for(std::list<double>::iterator a = angles_list.begin(); a != angles_list.end(); ++a){
    dev_angle += ((*a)-average_angle)*((*a)-average_angle);
  }
  dev_angle = sqrt(dev_angle/angles_list.size());
  average_dist_o/=position_list.size();
  for(std::list<double>::iterator a = dist_o_list.begin(); a != dist_o_list.end(); ++a){
    dev_dist_o += ((*a)-average_dist_o)*((*a)-average_dist_o);
  }
  dev_dist_o = sqrt(dev_dist_o/dist_o_list.size());
  std::cout << "Average angle between C-O and N-H: " << std::endl;
  std::cout << average_angle << "  " << dev_angle << std::endl;
  std::cout << "Average distance between reconstructed O and original one: " << std::endl;
  std::cout << average_dist_o << "  " << dev_dist_o << std::endl;
  
  // Fourth output: bond length and angles
  std::cout << "###############################################" << std::endl;
  std::cout << "Bond length and angles" << std::endl;
  std::cout << "###############################################" << std::endl;
  double c_n = 0., n_ca = 0., ca_c = 0., ca_ha = 0., n_ca_c = 0., ca_c_n = 0., c_n_ca = 0.;
  std::list<double> c_n_list, n_ca_list, ca_c_list, ca_ha_list, n_ca_c_list, ca_c_n_list, c_n_ca_list;
  double dev_c_n = 0., dev_n_ca = 0., dev_ca_c = 0.,  dev_ca_ha = 0., dev_n_ca_c = 0., dev_ca_c_n = 0., dev_c_n_ca = 0.;
  for(std::list<positions>::iterator p = position_list.begin(); p != position_list.end(); ++p){
    Eigen::Vector3d cn = p->n-p->pc;
    c_n += cn.norm();
    c_n_list.push_back( cn.norm() );
    Eigen::Vector3d nca = p->ca-p->n;
    n_ca += nca.norm();
    n_ca_list.push_back( nca.norm() );
    Eigen::Vector3d cac = p->c-p->ca;
    ca_c += cac.norm();
    ca_c_list.push_back( cac.norm() );
    Eigen::Vector3d csn = p->sn-p->c;
    
    Eigen::Vector3d caha = p->ha-p->ca;
    ca_ha += caha.norm();
    ca_ha_list.push_back( caha.norm() );
    
    cn.normalize();
    nca.normalize();
    cac.normalize();
    csn.normalize();
    c_n_ca += acos( nca.dot(cn) );
    n_ca_c += acos( cac.dot(nca) );
    ca_c_n += acos( csn.dot(cac) );
    c_n_ca_list.push_back( acos( nca.dot(cn) ) );
    n_ca_c_list.push_back( acos( cac.dot(nca) ) );
    ca_c_n_list.push_back( acos( csn.dot(cac) ) );
  }
  c_n/=position_list.size();
  n_ca/=position_list.size();
  ca_c/=position_list.size();
  ca_ha/=position_list.size();
  n_ca_c/=position_list.size();
  ca_c_n/=position_list.size();
  c_n_ca/=position_list.size();
  
  std::cout << c_n << "  " << dev_c_n << std::endl;
  std::cout << n_ca << "  " << dev_n_ca << std::endl;
  std::cout << ca_c << "  " << dev_ca_c << std::endl;
  std::cout << ca_ha << "  " << dev_ca_ha << std::endl;
  std::cout << n_ca_c << "  " << dev_n_ca_c << std::endl;
  std::cout << ca_c_n << "  " << dev_ca_c_n << std::endl;
  std::cout << c_n_ca << "  " << dev_c_n_ca << std::endl;
  std::cout << "###############################################"  << std::endl;
  return 0;
}
