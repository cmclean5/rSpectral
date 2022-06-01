#include "Headers.h"
#include "network.h"
#include "readfile.h"
#include "SpectralModularity.h"

// [[Rcpp::export]]
Rcpp::List spectral( Rcpp::DataFrame     DF,
                     Rcpp::IntegerVector CnMIN=1,
                     Rcpp::NumericVector TOL=0.00001,
                     Rcpp::IntegerVector names=1,
                     Rcpp::IntegerVector verbose=0){
  //Rcpp::List spectral( Rcpp::List params ){

  //For more information wrapping and packaging C/C++ in R see:
  //[1] https://www.gormanalysis.com/blog/exposing-a-cpp-student-class-with-rcpp/ (building R package and adding your C/C++ code)
  //[2] https://www.youtube.com/watch?v=DWkIbk_HE9o (setting up roxygen in R package)
  //[3] http://web.mit.edu/insong/www/pdf/rpackage_instructions.pdf
  //[4] http://r-pkgs.had.co.nz/src.html

  //Development steps
  //1) Open RStudio
  //2) edit/change code
  //3) Run pkgbuild::compile_dll() to compile your C++ code
  //4) Run devtools::document() to automatically build the NAMESPACE file package documentation.
  //5) Run devtools::load_all() to compile the C++ code and load our package

  //To build and install package into R
  //1)  First need to edit the NAMESPACE file, and add:
  //2)  export(spectral), i.e. the function names (from this file) we want to use
  //3)  cd /afs/inf.ed.ac.uk/user/c/cmclean5/WORK/STUDIES
  //4)  Run R CMD build CDMSuite
  //5   Run R CMD INSTALL CDMSuite_0.1.0.tar.gz
  //6)  Start R
  //7)  library(CDMSuite)
  //8)  CDMSuite::spectral(...)

  int i,j,k,KK;

  int Cn_min       = 1;
  double tol       = 0.00001;
  int N            = 0;
  int M            = 0;
  double *A        = 0;
  edgelist *el     = 0;
  bool useLoops    = false;
  bool checkM      = true;
  bool print       = false;
  int alphaNumeric = 1;

  //initialise network
  readfile *reader          = 0;
  network *gg               = new network();
  SpectralModularity *model = 0;


  int ncols = DF.length();
  int nrows = DF.nrows();

  if( (ncols > 0) && (nrows > 0) ){

    //set size for our DATASET
    KK              = nrows*ncols;
    string *DATASET = new string[KK];
    
    if( CnMIN.length() == 1 ){
      if( (CnMIN[0] > 0) ){
	Cn_min = CnMIN[0];
      }
    }

    if( (TOL.length() == 1) ){
      if( (TOL[0] > 0) ){
	  tol = TOL[0];
	}
    }

    if( (names.length() == 1) ){
      if( names[0] == 0 ){
      	alphaNumeric = 0;
      }
    }

    if( verbose.length() == 1 ){
      if( verbose[0] == 1 ){
        print = 1;
      }
    }

    if( ncols == 2 ){
      //unweighted networks
      
      Rcpp::StringVector       V1 = DF[0];
      Rcpp::StringVector       V2 = DF[1];
      
      for(k=0; k<KK; k++){
	i = floor(k/ncols);
	j = k % ncols;

	Rcpp::String v1(V1[i]);
	Rcpp::String v2(V2[i]);
	
	//string s1 = Rcpp::as< std::string >(V1(i)); 
	//string s2 = Rcpp::as< std::string >(V2(i));
	
	if( j == 0 ){ DATASET[(i*ncols)+j] = v1.get_cstring(); }
	if( j == 1 ){ DATASET[(i*ncols)+j] = v2.get_cstring(); }
	
      }    
      
    }

    if( ncols == 3 ){
      //for moment, lets not considering weighted.
      
      Rcpp::StringVector       V1 = DF[0];
      Rcpp::StringVector       V2 = DF[1];
      Rcpp::StringVector       V3 = DF[2];
      
      for(k=0; k<KK; k++){
	i = floor(k/ncols);
	j = k % ncols;

	Rcpp::String v1(V1[i]);
	Rcpp::String v2(V2[i]);
	Rcpp::String v3(V3[i]);
	
	if( j == 0 ){ DATASET[(i*ncols)+j] = v1.get_cstring(); }
	if( j == 1 ){ DATASET[(i*ncols)+j] = v2.get_cstring(); }
	if( j == 2 ){ DATASET[(i*ncols)+j] = v3.get_cstring(); }
	
      }    
      
    }

    
    
    //load edgelist into network
    reader = new readfile( gg, DATASET, ncols, nrows, alphaNumeric );

    //build Adjaceny Matrix
    gg->buildNetworkReps( useLoops, checkM );
    //---

    N = gg->getN();
    M = gg->getM2();
    A = gg->getA();

    //gg->setPrint(true);
    //gg->printVertices();
    
    if( N != 0 && M != 0 ){

      //set-up clustering alg.
      model = new SpectralModularity(gg,el,A,N,M,print);
      //model->setPrint(print);
      model->settol( tol );
      model->setMinCn( Cn_min );
      
      //--- run spectral clustering
      int cal = model->calculateSpectralModularity();

      //--- reorder community numbers in network
      gg->reorderK();

    }
        

  }

  
  if( gg->getN() > 0 ){
    
    N = gg->getN();

    //--- output the node label and its cluster
    Rcpp::StringVector  ID   (N);
    Rcpp::NumericVector Coms (N);

    for( i=0; i<N; i++ ){
      ID[i]   = gg->V[i].label;
      Coms[i] = gg->V[i].K;
    }

    // Create a named list with the above quantities
    return Rcpp::List::create(Rcpp::Named("ID") = ID,
			      Rcpp::Named("K")  = Coms);

  } else {
    
    //--- output the node label and its cluster
    Rcpp::StringVector  ID   (0);
    Rcpp::NumericVector Coms (0);

    // Create a named list with the above quantities
    return Rcpp::List::create(Rcpp::Named("ID") = ID,
			      Rcpp::Named("K")  = Coms);

    }


}//spectral

