
#ifndef DFTMETHOD_H
#define DFTMETHOD_H

#include "iostream"
#include "iomanip"
#include <math.h>

#include "TObject.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TGraph.h"

#include "SPEResponse.h"

#include "fftw3.h"

#include "TMinuit.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Minuit2/FunctionMinimum.h" 
#include "Minuit2/MnMigrad.h" 
#include "Minuit2/MnUserParameters.h" 
#include "Minuit2/MnPrint.h" 
#include "Minuit2/FCNBase.h" 
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

class DFTmethod : public TObject
{
 private:

  SPEResponse spef;

  TH1D* hdist;
      
  Double_t wbin;
  Double_t lo_edge;
  Double_t hi_edge;
  
  unsigned int N;
  unsigned int M;
      
  ROOT::Minuit2::Minuit2Minimizer *mFFT;
  
 public:
  
  DFTmethod();
  
  virtual ~DFTmethod();
  
  DFTmethod( SPEResponse _spef );

  void SetHisto( TH1D* _h );
    
  void CalculateValues();
  TGraph* GetGraph();

  Double_t Norm;

  Double_t Q0;
  Double_t s0;

  Double_t mu;

  Double_t fftPhase( Double_t vy, Double_t vz );
  
  void Fit();

  Double_t vals[20];
  Double_t errs[20];
  
  ClassDef( DFTmethod, 1 )
        
};

#endif
