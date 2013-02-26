/* ************************************************************************** */
/*                                                                            */
/*    Copyright 2013 Stefano Zamuner                                          */
/*                                                                            */
/*    This file is part of BioCpp.                                            */
/*                                                                            */
/*    BioCpp is free software: you can redistribute it and/or modify          */
/*    it under the terms of the GNU General Public License as published by    */
/*    the Free Software Foundation, either version 3 of the License, or       */
/*    (at your option) any later version.                                     */
/*                                                                            */
/*    BioCpp is distributed in the hope that it will be useful,               */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*    GNU General Public License for more details.                            */
/*                                                                            */
/*    You should have received a copy of the GNU General Public License       */
/*    along with BioCpp.  If not, see <http://www.gnu.org/licenses/>.         */
/*                                                                            */
/* ************************************************************************** */

#ifndef STRICT_NEEDLEMAN_WUNSCH_ALIGNMENT
#define STRICT_NEEDLEMAN_WUNSCH_ALIGNMENT

#include <string>
#include <algorithm>

#include "../utils/errors_and_warnings.h"

namespace BioCpp{
namespace fasta{

/*! \brief strict version of Needleman-Wunsch alignment algorithm */
class StrictNeedlemanWunschAlignment{
  private:
    double Match(char a, char b, double same, double hole, double different);
    double HoleInsert(char a, double allowed, double forbidden);
  public:
    StrictNeedlemanWunschAlignment(){};
    error operator() (std::string& tseq, std::string& rseq);
};

double StrictNeedlemanWunschAlignment::Match(char a, char b, double same, double hole, double different){
  /* scoring parameters */
  char A = std::toupper(a);
  char B = std::toupper(b);
  if( A==B and A!='X' )
    return same;
  else if( A == 'X' or B == 'X' )
    return same;
  else if( A=='X' and B=='X' )
    return hole;
  return different;
}

double StrictNeedlemanWunschAlignment::HoleInsert(char a, double allowed, double forbidden){
  if(std::islower(a)){
    return allowed;
  }
  return forbidden;
}

error StrictNeedlemanWunschAlignment::operator()(std::string& tseq, std::string& rseq){
  double epsilon = 0.0001;
  double same = 1., all_hole = -0.011, for_hole=-1200.758, diff = -1000.531;

  /* input string have to be uppercase */
  std::transform(tseq.begin(), tseq.end(),tseq.begin(), ::toupper);
  std::transform(rseq.begin(), rseq.end(),rseq.begin(), ::toupper);
  
  
  unsigned int tsize = tseq.size(), rsize = rseq.size();
  std::string n_tseq, n_rseq;
  
  bool t_beg_hole = tseq[0]=='-' ? true : false;
  bool t_end_hole = tseq[tseq.size()-1]=='-' ? true : false;
  
  /* control if input strings are valid */
  for(unsigned int t = 0; t < tsize; ++t){
    char tres = tseq[t];
    unsigned int nsize = n_tseq.size();
    if( tres == '-' ){
      if(nsize!=0)
        n_tseq[nsize-1] = std::tolower(n_tseq[nsize-1]);
      continue;
    }
    else if( tres >='A' and tres <='Z' ){
      n_tseq += tres  ;
    }
    if( tres=='O' or tres=='J' or tres=='U' or tres<'A' or tres>'Z' ){
      return ALIGN_NOT_A_VALID_TSEQRES;
    }
  }
  for(unsigned int r = 0; r < rsize; ++r){
    char rres = rseq[r];
    unsigned int nsize = n_rseq.size();
    if( rres == '-' ){
      if(nsize!=0)
        n_rseq[nsize-1] = std::tolower(n_rseq[nsize-1]);
      continue;
    }
    else if( rres >='A' and rres <='Z' ){
      n_rseq += rres  ;
    }
    if( rres=='O' or rres=='J' or rres=='U' or rres<'A' or rres>'Z' ){
      return ALIGN_NOT_A_VALID_RSEQRES;
    }
  }
  
  /* fill the matrix */
  tsize = n_tseq.size(); rsize = n_rseq.size();
  Eigen::MatrixXd F = Eigen::MatrixXd::Zero(tsize+1,rsize+1);
  for(unsigned int t = 0; t < tsize+1; ++t){
    F(t,0) = t*all_hole;
  }
  for(unsigned int r = 0; r < rsize+1; ++r){
    F(0,r) = r*all_hole;
  }
  for(unsigned int t = 0; t < tsize; ++t){
    for(unsigned int r = 0; r < rsize; ++r){
      double w = Match( n_tseq[t], n_rseq[r], same, all_hole, diff );
      double match = F(t,r) + w;
      double deletion = F(t,r+1) + HoleInsert(n_rseq[r], all_hole, for_hole);
      double insertion = F(t+1,r) + HoleInsert(n_tseq[t], all_hole, for_hole);
      F(t+1,r+1) = std::max( match, std::max(deletion, insertion) );
    }
  }
  
  std::string a_tseq, a_rseq;
  /* check if alignment is possible */
  Eigen::MatrixXf::Index maxcol, maxrow;
  F.maxCoeff(&maxrow, &maxcol);
  if( maxrow!=(int)tsize and maxcol!=(int)rsize )
    return ALIGN_FAILED;
  /* align TAIL */
  else if( maxrow==(int)tsize and maxcol!=(int)rsize){
    while( (int)rsize > maxcol ){
      --rsize;
      if(std::islower(n_rseq[rsize]) ){
        a_tseq = '-' + a_tseq;
        a_rseq = '-' + a_rseq;
      }
      char tres = std::toupper(n_rseq[rsize]), rres = t_end_hole ? tres : '*';
      a_tseq = tres + a_tseq;
      a_rseq = rres + a_rseq;
    }
  }
  else if( maxcol==(int)rsize and maxrow!=(int)tsize){
    while( (int)tsize > maxrow ){
      --tsize;
      if(std::islower(n_tseq[tsize]) ){
        a_tseq = '-' + a_tseq;
        a_rseq = '-' + a_rseq;
      }
      char tres = std::toupper(n_tseq[tsize]), rres = '^';
      a_rseq += rres;
      a_tseq =  tres + a_tseq;
    }
  }
  /* align CENTER */
  while( tsize > 0 and rsize > 0 ){
    double w = Match( n_tseq[tsize-1], n_rseq[rsize-1], same, all_hole, diff );
    if(tsize > 0 and rsize > 0 and fabs( F(tsize, rsize)-( F(tsize-1, rsize-1)+w ) ) <= epsilon ){
      --tsize;
      --rsize;
      if( std::islower(n_tseq[tsize]) and std::islower(n_rseq[rsize]) ){
        a_tseq = '-' + a_tseq;
        a_rseq = '-' + a_rseq;
      }
      char t_res = std::toupper(n_tseq[tsize]) == 'X' ? std::toupper(n_rseq[rsize]) : std::toupper(n_tseq[tsize]);
      char r_res = std::toupper(n_rseq[rsize]) == 'X' ? t_res : std::toupper(n_rseq[rsize]);
      a_tseq = t_res + a_tseq;
      a_rseq = r_res + a_rseq;
    }
    else if( tsize > 0 and fabs( F(tsize,rsize) - ( F(tsize-1,rsize)+HoleInsert(n_rseq[rsize-1], all_hole, for_hole) ) ) < epsilon ){
      --tsize;
      if( std::islower(n_tseq[tsize]) ){
        a_tseq = '-' + a_tseq;
        a_rseq = '-' + a_rseq;
      }
      char t_res = std::toupper(n_tseq[tsize]), r_res = '^';
      a_tseq = t_res + a_tseq;
      a_rseq = r_res + a_rseq;
    }
    else if( rsize > 0 and fabs( F(tsize,rsize) - ( F(tsize,rsize-1)+HoleInsert(n_tseq[tsize-1], all_hole, for_hole) ) ) < epsilon ){
    	--rsize;
      if( std::islower(n_rseq[rsize]) ){
        a_tseq = '-' + a_tseq;
        a_rseq = '-' + a_rseq;
      }
    	char r_res = std::toupper(n_rseq[rsize]), t_res = r_res;
      a_tseq = t_res + a_tseq;
      a_rseq = r_res + a_rseq;
    }
    else{
      return UNSPECIFIED;
    }
    
    /* align HEAD */
    if( tsize==0 and rsize!=0 ){
    	while(rsize > 0){
    	  --rsize;
    	  char t_res = std::toupper( n_rseq[rsize] ), r_res = t_beg_hole ? t_res : '*';
    	  a_tseq = t_res + a_tseq;
    	  a_rseq = r_res + a_rseq;
    	}
    }
    else if ( rsize==0 and tsize!=0 ){
    	while( tsize>0 ){
    	  --tsize;
    	  char t_res = std::toupper( n_tseq[tsize] ), r_res = '^';
    	  a_tseq = t_res + a_tseq;
    	  a_rseq = r_res + a_rseq;
    	}
    }

  }

  tseq = a_tseq;
  rseq = a_rseq;
	return ERR_NONE;
}

} // end fasta
} // end biocpp
#endif
