
#ifndef PMT_H
#define PMT_H

#include "iostream"
#include "iomanip"

#include "TObject.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom3.h"

#include "SPEResponse.h"

class PMT : public TObject
{
 private:
  
  TH1D *spectrum;
  
  Int_t nbins;
  Double_t min;
  Double_t max;
  
  //Pedestal ped;
  SPEResponse res;
  
  
 public:
  
  PMT();
  
  virtual ~PMT();
  
  PMT( Int_t _nbins, Double_t _min, Double_t _max, SPEResponse _res );
  // Pedestal _ped ); , PMTModel _mod );
    
  void GenSpectrum( Int_t ntot );
  
  TH1D* GetSpectrum(){ return spectrum; };
  
  void DrawSpectrum()
  {
    spectrum->SetMaximum( 2.5*spectrum->GetBinContent( spectrum->GetMaximumBin() ) );
    spectrum->Draw( "PEZ" );
      
  };
    
  ClassDef( PMT, 1 )
    
};

#endif
