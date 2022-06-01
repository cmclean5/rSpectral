//define guards, so headers are declare only once.
#ifndef SPECTRALMODULARITY_H
#define SPECTRALMODULARITY_H

/*
  Use of the stack for memory allocation. This 
 should be faster for large networks, but will need to reset NSIZE 
 large enough for your network size, and then re-Make. 

 */

#include "Headers.h"
#include "network.h"

class SpectralModularity {

 public:
   SpectralModularity();
  SpectralModularity(network *, edgelist *, double *, int, int, bool=false);
  ~SpectralModularity();
  int calculateSpectralModularity();
  void setMinCn( int );
  void settol  ( double );
  void setPrint( bool );
  void setEignOpts( double, int, int );
  
 private:
  void calculateB( double *, int );
  int  delta( int, int );
  void split( double *, int, int *, const char* );
  void deltaModularityMax( int, double & );
  void maxModularity( double & );
  void modifySplit( int );
  void deltaModularity( double & );
  void maximiseIndexVectors();
  void calculateEigenVectors();
  void setupMatrices();
  void assignSpace();
  void freeSpace();

  
  //global variable values  
  //static const bool PRINT = true;
  bool PRINT;
  static const int dummy  = -1000;
  //static const int NSIZE  = 5000;
  //static const int MSIZE  = NSIZE * NSIZE;

  double tol;//the tolerance value, 10^-5; eigenvalues below this threshold are not used
  int MINCn;//The minimum cluster size

  network *gg;
  double *A;
  double *Bgi;  //The Modularity matrix, Bgi
  int    NR_Bgi;//Number of rows of Bgi
  int    NC_Bgi;//Number of cols of Bgi
  int    M;//number of edges
  bool   usedBgi;
  
  double specQ;
  double NORM;
  double betai;

  int    MAXK;  //Counter storing the maximum community number so far
  
  double *u;
  double *Bgi_temp;  
  int    *SI;
  int    *si;
  int    *visited;
  int    *keys_p;
  int    *keys_n;

  bool usedu;  
  bool usedBgi_temp;
  bool usedSI;
  bool usedsi;
  bool usedvisited;
  bool usedkeys_p;
  bool usedkeys_n;
  
  arma::eigs_opts opts;
  
};

#endif
