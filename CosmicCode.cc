//-----------------------------------------------------------------------------
// File: CosmicCode
// Description: Fit 3-parameter models to the latest compilation of Type Ia 
//              supernova data.
//
//              a    - universal scale factor
//
//              OM   - Omega_M
//              OL   - Omega_L
//              H    - roughly related to Hubble's constant
//
//              For some models OM and OL may have different meanings.
//
// Created: June 2008 HBP
// Updated for ESHEP 2012, La Pommeraye, Anjou, France
//-----------------------------------------------------------------------------
#include <cmath>
#include <cassert>
#include <iostream>
#include <string>
using namespace std;
//-----------------------------------------------------------------------------
struct CosmicCode
{
  struct CModel
  {
    int model;

    CModel(int _model=0) : model(_model)
    {
      switch (model)
	{
	case 0: // LCDM
	default:
	  cout << endl << "\tLCDM model" << endl << endl;
	  break;
	case 1: // phantom
	  cout << endl << "\tphantom model" << endl << endl;
	  break;
	}    
    }
    ~CModel() {}
  
    double operator()(double a, double* p)
    {
      // p[0] = OL
      // p[1] = OM
      // p[2] = H0
      // p[3] = n
      assert(p);
      double y = 0;
      switch (model)
	{
	case 0: // LCDM
	default:
	  {
	    // a^3 * [Omega_M/a^3 + (1-Omega_M-Omega_L)/a^2 + Omega_L]
	    double OM = p[0];
	    double OL = p[1];
	    y = OM + (1 - OM - OL)*a + OL*a*a*a;
	    if ( y < 0 ) y = 1.e20;
	  }
	  break;
	case 1: // phantom
	  {
	    // a^3 * [ exp(a^n-1)/a^3 ]
	    double n = p[3];
	    y = exp(pow(a, n)-1);
	  }
	  break;
	}
      return y;
    }
  };
    
    
  CModel model;
  double offset;
  int N;
  
  CosmicCode(int _ID, int _N=200)
    : model(CModel(_ID)),
      N(_N),
      offset(5*log10(2.99*pow(10.0, 5.0)) + 25)
  {}
  ~CosmicCode() {}

  // compute distance modulus
  double distanceModulus(double z, double* p)
  {
    double OM = p[0];
    double OL = p[1];
    double H0 = p[2];
    
    double a = 1.0/(1+z);
    double h = (1-a) / N;
    double F = 0;
    for(int i=0; i < N; i++)
      {
	double x = a + (i+0.5)*h;
	F = F + 1.0/sqrt(x * model(x, p));
      }
    F = F*h;

    double OK = 1 - OM - OL;
    if ( OK != 0 )
      {
	double rootOK = sqrt(fabs(OK));
	double theta  = rootOK * F;
	if ( OK > 0 )
	  F = sinh(theta);  // sinh
	else
	  F = sin(theta);   // sin
	F /= rootOK;
      }

    double y = 5*log10( (1+z) * F / H0) + offset;
    return y;
  }
  

  // compute lifetime vs a
  void scaleFactor(double amax, double* p, double* t, double* a)
  {
    double F = 0;
    double h = amax / N;
    for(int i=0; i < N; i++)
      {
	double x = (i+0.5)*h;
	F = F + sqrt(x / model(x, p));
	a[i] = x + 0.5*h;
	t[i] = F*h;
      }
  }


  // compute comoving distance vs. a
  void comovingDistance(double amax, double* p, double* chi, double* a)
  {
    double F = 0;
    double h = amax / N;
    for(int i=0; i < N; i++)
      {
	double x = (i+0.5)*h;
	F = F + 1.0 / sqrt(x * model(x, p));
	a[i] = x + 0.5*h;
	chi[i] = F*h;
      }
  }

  // compute Omega(a)
  void Omega(double amax, double* p, double* a, double* O)
  {
    double h = amax / N;
    for(int i=0; i < N; i++)
      {
	double x = (i+0.5)*h;
	a[i] = x + 0.5*h;
	O[i] = model(a[i], p) / pow(a[i], 3);
      }
  }

};
