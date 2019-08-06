
#include "Pedestal.h"

using namespace std;

ClassImp( Pedestal )

Double_t _pedfunc( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];
    
  Double_t Q_0 = par[0];
  Double_t s_0 = par[1];
        
  Double_t arg = 0.0; 
  if ( s_0!=0.0 ) arg = ( xx - Q_0 )/s_0;    
  else cout << "Error: The code tries to divide by zero." << endl;
  
  Double_t result = 1.0/( sqrt( 2 * TMath::Pi() ) * s_0 ) * TMath::Exp( -0.5*arg*arg );
    
  return result;
  
}

Pedestal::Pedestal()
{}

Pedestal::~Pedestal()
{}

Pedestal::Pedestal( Double_t _Q0, Double_t _s0 )
{
  Q0 = _Q0;
  s0 = _s0;
  
  pedfunc = new TF1( "pedfunc", _pedfunc, Q0-25.0*s0, Q0+25.0*s0, 2 );
  pedfunc->SetLineColor( kRed );
  pedfunc->SetLineWidth( 2.0 );
  pedfunc->SetNpx( 10000 );
  pedfunc->SetParameters( Q0, s0 );
  
}

Double_t Pedestal::GenQ()
{
  Double_t _x = pedfunc->GetRandom();

  return _x;
  
}
