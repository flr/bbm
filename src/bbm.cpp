#define TMB_LIB_INIT R_init_bbm
#include <TMB.hpp>
template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_SCALAR(grec);
  DATA_SCALAR(gadu);
  DATA_INTEGER(nindex);		// number of surveys
  DATA_MATRIX(Crec);		// matrix of recruits catch with dimensions nyr x nstep
  DATA_MATRIX(Cadu);		// matrix of total catch with dimensions nyr x nstep
  DATA_VECTOR(Bobs);		// vector of observed total biomass
  DATA_VECTOR(Pobs);		// vector of observed recruit proportion 
  DATA_IMATRIX(idxBobs);	// matrix 2 colums: survey and index of year
  DATA_IMATRIX(idxPobs);	// matrix 2 colums: survey and index of year
  DATA_IVECTOR(indexper);	// index of the time step of each survey	
  DATA_VECTOR(f);			// fraction of the year corresponding to each time step

//  DATA_IVECTOR(idxpobs);
//  DATA_IVECTOR(idxbobs);

  PARAMETER_VECTOR(logq);
  PARAMETER_VECTOR(logpsi);
  PARAMETER_VECTOR(xi);
  PARAMETER(logB0);
  PARAMETER_VECTOR(logR);
  PARAMETER(mur);
  PARAMETER(logpsir);

  int nbobs = Bobs.rows();
  int npobs = Pobs.rows();
  int nyr = Cadu.rows();
  int nstep = Cadu.cols();
 
  vector<Type> q(nindex);
  q = exp(logq);
  vector<Type> psi(nindex);
  psi = exp(logpsi);
  Type B0 = exp(logB0);
  Type psir = exp(logpsir);
  vector<Type> R(nyr);
  R = exp(logR);
  
  matrix<Type> Brec(nyr, nstep+1);
  matrix<Type> Badu(nyr, nstep+1);
  matrix<Type> Btot(nyr, nstep+1);
  matrix<Type> P(nyr, nstep+1);

  vector<Type> Bpred(nbobs);
  vector<Type> Ppred(npobs);

  ADREPORT(q);
  ADREPORT(psi);
  ADREPORT(B0);
  ADREPORT(R);
  ADREPORT(psir);  
  
  Type nll;
  nll= 0.0;

  Type eps;
  eps=1.0e-3;
  Type pen;

  // First year
  
  Brec(0,0) = R(0);
  Badu(0,0) = B0;
  Btot(0,0) = Brec(0,0) + Badu(0,0);
  P(0,0) = Brec(0,0)/Btot(0,0);
  
  for (int j=1; j<=nstep; j++){
	pen=0; 
//	Brec(0,j)=100000*grec;
	Brec(0,j) = Brec(0,j-1)*exp(-grec*f(j-1)) - Crec(0, j-1)*exp(-grec*f(j-1)/2.0);
	Brec(0,j) = posfun(Brec(0,j), eps, pen);
	nll += pen;
	
	pen=0; 
//    Badu(0,j) = 100000*gadu;
	Badu(0,j) = Badu(0,j-1)*exp(-gadu*f(j-1)) - Cadu(0, j-1)*exp(-gadu*f(j-1)/2.0);
	Badu(0,j) = posfun(Badu(0,j), eps, pen);
	nll += pen;

	Btot(0,j) = Brec(0,j) + Badu(0,j);
	P(0,j) = Brec(0,j)/Btot(0,j);
  }
  
  // rest of years 
  
  for (int i=1; i<nyr; i++){  
	Brec(i,0) = R(i);
	Badu(i,0) = Btot(i-1, nstep); 	
	Btot(i,0) = Brec(i,0) + Badu(i,0);
	P(i,0) = Brec(i,0)/Btot(i,0);
	for (int j=1; j<=nstep; j++){
		pen=0; 
		//Brec(i,j) = 100000;
		Brec(i,j) = Brec(i,j-1)*exp(-grec*f(j-1)) - Crec(i, j-1)*exp(-grec*f(j-1)/Type(2.0));
		Brec(i,j) = posfun(Brec(i,j), eps, pen);
		nll += pen;
		
		pen=0; 
		// Badu(i,j) = 100000;
		Badu(i,j) = Badu(i,j-1)*exp(-gadu*f(j-1)) - Cadu(i, j-1)*exp(-gadu*f(j-1)/Type(2.0));
		Badu(i,j) = posfun(Badu(i,j), eps, pen);
		nll += pen;

		Btot(i,j) = Brec(i,j) + Badu(i,j);
		P(i,j) = Brec(i,j)/Btot(i,j);
	}
  }

  // prediction for biomass observations
  
  for(int i=0; i<nbobs; ++i){
	int locsurv = idxBobs(i,0);
	int step = indexper(locsurv-1);
	int locyear = idxBobs(i,1);	
	Type val = Bobs(i);
	nll -= dnorm(log(val), log(q(locsurv-1))+log(Btot(locyear-1, step-1)), Type(1.0)/sqrt(psi(locsurv-1)), true);	 
 }
  
  // prediction for biomass proportion observations

  for(int i=0; i<npobs; ++i){
	int locsurv = idxPobs(i,0);
	int step = indexper(locsurv-1);
	int locyear = idxPobs(i,1);	
	Type val = Pobs(i);
	nll -= dbeta(val, exp(xi(locsurv-1))*P(locyear-1, step-1), exp(xi(locsurv-1))*(Type(1.0)-P(locyear-1, step-1)), true);
 }
  // recruitment deviations
  nll -= sum(dnorm(log(R), mur, Type(1.0)/sqrt(psir), true));	 
  
  return nll;
}
