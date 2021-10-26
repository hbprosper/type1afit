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
    int ID;

    CModel(int _id=0) : ID(_id)
    {
      switch (ID)
        {
        case 0: // LCDM
        default:
          cout << endl << "\tLCDM model" << endl << endl;
          break;
        case 1: // phantom
          cout << endl << "\tphantom model" << endl << endl;
          break;
        case 2: // CM
          cout << endl << "\tCM model" << endl << endl;
          break;
        }    
    }
    ~CModel() {}
  
    double operator()(double a, double* p)
    {
      assert(p);
      double y = 0;
      switch (ID)
	{
	case 0: // LCDM
	default:
	  {
	    // p[0]: OM
	    // p[1]: OL
	    // p[2]: H0
	   
	    // a^3 * [Omega_M/a^3 + (1-Omega_M-Omega_L)/a^2 + Omega_L]
	    double OM = p[0];
	    double OL = p[1];
	    y = OM + (1 - OM - OL)*a + OL*a*a*a;
	    if ( y < 0 ) y = 1.e20;
	  }
	  break;
	case 1: // phantom
	  {
	    // p[0]: n
	    // p[1]: H0
	    //
	    // a^3 * [ exp(a^n-1)/a^3 ]
	    double n = p[0];
	    y = exp(pow(a, n)-1);
	  }
	  break;
	case 2: // CM
	  {
	    // p[0]: OM
        // p[1]: H0
	    //
	    // a^3 * [Omega_M/a^3 + (1-Omega_M)/a^2]
	    double OM = p[0];
	    y = OM + (1 - OM)*a;
        if ( y < 0 ) y = 1.e20;
	  }
	  break;
	}
      return y;
    }
  };
    
    
  CModel model;
  double offset;
  int N;
  
  void* address()
  {
      return (void*)this;
  }
    
  CosmicCode(int _id, int _N=200)
    : model(CModel(_id)),
      N(_N),
      offset(5*log10(2.99*pow(10.0, 5.0)) + 25)
  {}
  ~CosmicCode() {}

  // compute distance modulus
  double distanceModulus(double z, double* p)
  {
    double OM;
    double OL;
    double H0;

    switch (model.ID)
      {
      case 0: // LCDM
      default:
	{
	  OM = p[0];
	  OL = p[1];
	  H0 = p[2];
	}
	break;
      case 1: // phantom
	{
	  OM = 1;
	  OL = 0;
	  H0 = p[1];
	}
	break;
      case 2: // CM
	{
	  OM = p[0];
      OL = 0;
	  H0 = p[1];
	}
	break;
      }
    
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
