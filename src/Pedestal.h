
#ifndef PEDESTAL_H
#define PEDESTAL_H

#include "TObject.h"
#include "TMath.h"
#include "TF1.h"

class Pedestal : public TObject
{
 private:
    
  TF1 *pedfunc;
      
 public:
  
  Pedestal();
  
  virtual ~Pedestal();

  Double_t Q0;
  Double_t s0;
  
  Pedestal( Double_t _Q0, Double_t _s0 );
  
  Double_t GenQ();

  
  ClassDef( Pedestal, 1 )
        
};

#endif
