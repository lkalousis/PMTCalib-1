
#ifndef PMTMODEL_H
#define PMTMODEL_H

#include "TObject.h"
#include "TF1.h"

#include "PMType.h"

class PMTModel : public TObject
{
 private:

  Double_t params[20]={-1.0};
        
 public:
  
  PMTModel();
  
  virtual ~PMTModel();

  PMTModel( PMType::Model _modtype );

  PMType::Model modtype;
  
  Int_t nparams;
   
  void SetParams( Double_t _params[] );
  Double_t Value( Double_t xx );
  
  Double_t F1( Double_t xx ); // SIMPLE GAUSS
    
  
  ClassDef( PMTModel, 1 )
    
};

#endif
