#include "Headers.h"
#include "network.h"
#include "readfile.h"
#include "SpectralModularity.h"

// global
readfile *reader          = nullptr;
string *dataset           = nullptr;
network *gg               = nullptr;
SpectralModularity *model = nullptr;

//' Clean things up
//'
//' @export
// [[Rcpp::export]]
 void freeSpace(){

   // delete c++ objects
   if( gg      != nullptr ){ delete gg; }
   if( reader  != nullptr ){ delete reader; }
   if( dataset != nullptr ){ delete[] dataset; }
   
 }


//' Load edge list for analysis
//' 
//' This function read edge list and create the network for analysis. It does 
//' not return anything it just creates required structures in memory.
//'
//' @param df edge list
//' @param names are we dealing with alphaNumeric (1) or numeric (!1) ids
//'
//' @export
//'
//' @examples
//' library(igraph)
//' g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
//' g <- add_edges(g, c(1,6, 1,11, 6, 11))
//' el = as.data.frame(get.edgelist(g,names=TRUE))
//' rSpectral::load_data(df=el)
//' status = rSpectral::spectral(fix_neig=0)
//' spec   = rSpectral::membership(detach_graph=1)
// [[Rcpp::export]]
void load_data ( Rcpp::DataFrame     df,
                 Rcpp::IntegerVector names=1 ){

  int i,j,k,KK;
 
  bool useLoops    = false;
  bool checkM      = true;   
  int alphaNumeric = 1;

  int ncols = df.length();
  int nrows = df.nrows();

  if( (ncols > 0) && (nrows > 0) ){
    
    if( (names.length() == 1) ){
      if( names[0] == 0 ){
      	alphaNumeric = 0;
      }
    }    
    
    //set size for our dataset
    KK              = nrows*ncols;
    dataset         = new string[KK];
      
    
    if( ncols == 2 ){
      //unweighted networks
      
      Rcpp::StringVector       V1 = df[0];
      Rcpp::StringVector       V2 = df[1];
      
      for(k=0; k<KK; k++){
        i = floor(k/ncols);
        j = k % ncols;

        Rcpp::String v1(V1[i]);
        Rcpp::String v2(V2[i]);
	
        if( j == 0 ){ dataset[(i*ncols)+j] = v1.get_cstring(); }
        if( j == 1 ){ dataset[(i*ncols)+j] = v2.get_cstring(); }
	
      }    
      
    }

    if( ncols == 3 ){
      //for moment, lets not considering weighted.
      
      Rcpp::StringVector       V1 = df[0];
      Rcpp::StringVector       V2 = df[1];
      Rcpp::StringVector       V3 = df[2];
      
      for(k=0; k<KK; k++){
        i = floor(k/ncols);
        j = k % ncols;

        Rcpp::String v1(V1[i]);
        Rcpp::String v2(V2[i]);
        Rcpp::String v3(V3[i]);
        
        if( j == 0 ){ dataset[(i*ncols)+j] = v1.get_cstring(); }
        if( j == 1 ){ dataset[(i*ncols)+j] = v2.get_cstring(); }
        if( j == 2 ){ dataset[(i*ncols)+j] = v3.get_cstring(); }
	
      }    
      
    }
    
    
    //load edgelist into network
    gg     = new network();
    reader = new readfile( gg, dataset, ncols, nrows, alphaNumeric );

    //build Adjaceny Matrix    
    gg->buildNetworkReps( useLoops, checkM );
   
  
  }

}//load_data

//' Spectral modularity clustering
//'
//' @param names are we dealing with alphaNumeric (1) or numeric (!1) ids
//' @param verbose should be removed
//' @param summary should be removed
//' @param Cn_min minimum cluster size
//' @param tol tolerance
//' @param fix_neig wether fixing neighbouring nodes found in same community
//'
//' @return status (does it really return something?)
//' @export
//' @examples
//' library(igraph)
//' g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
//' g <- add_edges(g, c(1,6, 1,11, 6, 11))
//' el = as.data.frame(get.edgelist(g,names=TRUE))
//' rSpectral::load_data(df=el)
//' status = rSpectral::spectral(fix_neig=0)
//' spec   = rSpectral::membership(detach_graph=1)
//' V(g)$color[as.numeric(spec$ID)]<-RColorBrewer::brewer.pal(max(spec$K),'Set1')[spec$K]
//' plot(g)
// [[Rcpp::export]]
void spectral( Rcpp::IntegerVector Cn_min=1,
               Rcpp::NumericVector tol=0.00001,
               Rcpp::IntegerVector names=1,
               Rcpp::IntegerVector fix_neig=0,
               Rcpp::IntegerVector verbose=0,
               Rcpp::IntegerVector summary=0){

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

  int CnMIN        = 1;
  double TOL       = 0.00001;
  int N            = 0;
  int M            = 0;
  double *A        = nullptr;
  edgelist *el     = nullptr;
  bool print       = false;
  bool modelSummary= false;
  bool neigFix     = false;
  int alphaNumeric = 1;

  if( gg != nullptr ){    
    
    if( Cn_min.length() == 1 ){
      if( (Cn_min[0] > 0) ){
        CnMIN = Cn_min[0];
      }
    }

    if( (tol.length() == 1) ){
      if( (tol[0] > 0) ){
        TOL = tol[0];
      }
    }

    if( (names.length() == 1) ){
      if( names[0] == 0 ){
      	alphaNumeric = 0;
      }
    }

    if( fix_neig.length() == 1 ){
      if( fix_neig[0] == 1 ){
        neigFix = 1;
      }
    }
    
    if( verbose.length() == 1 ){
      if( verbose[0] == 1 ){
        print = 1;
      }
    }

    if( summary.length() == 1 ){
      if( summary[0] == 1 ){
         modelSummary = 1;
      }
    }
     
    N = gg->getN();
    M = gg->getM2();
    A = gg->getA();

    //gg->setPrint(true);
    //gg->printVertices();
    
    if( N != 0 && M != 0 ){

      //set-up clustering alg.
      model = new SpectralModularity(gg,el,A,N,M,neigFix,print,modelSummary);
      //model->setPrint(print);
      model->settol( TOL );
      model->setMinCn( CnMIN );
      
      //--- run spectral clustering
      int cal = model->calculateSpectralModularity();

      //--- reorder community numbers in network
      gg->reorderK();

    }
        

  }

  // delete model
  if( model != nullptr ){ delete model; }
  
}//spectral

//' Return membership data
//'
//' @param detach_graph whether you want to keep graph in memory
//'
//' @return membership vector
//' @export
//'
//' @examples
//' library(igraph)
//' g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
//' g <- add_edges(g, c(1,6, 1,11, 6, 11))
//' el = as.data.frame(get.edgelist(g,names=TRUE))
//' rSpectral::load_data(df=el)
//' status1 = rSpectral::spectral(fix_neig=0)
//' spec1   = rSpectral::membership(detach_graph=0)
//' status2 = rSpectral::spectral(fix_neig=1)
//' spec2   = rSpectral::membership(detach_graph=0)
// [[Rcpp::export]]
 Rcpp::List membership( Rcpp::IntegerVector detach_graph=1 ){

   int i,N,detach;

   detach = 1;
   
  if( gg->getN() > 0 ){

    if( detach_graph.length() == 1 ){
      if( detach_graph[0] == 0 ){
         detach = 0;
      }
    }     
    
    N = gg->getN();
    
    //--- output the node label and its cluster
    Rcpp::StringVector  ID   (N);
    Rcpp::NumericVector Coms (N);

    for( i=0; i<N; i++ ){
      ID[i]   = gg->V[i].label;
      Coms[i] = gg->V[i].K;
    }   


    //detach graph
    if( detach ){ freeSpace(); }
    
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


 }//membership
