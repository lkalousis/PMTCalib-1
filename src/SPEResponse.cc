
#include "SPEResponse.h"

using namespace std;

ClassImp( SPEResponse )


Double_t _gausexpfunc( Double_t *x, Double_t *par )
{
  Double_t xx = x[0];
    
  Double_t Q = par[0];
  Double_t s = par[1];

  Double_t lambda = par[2];
  
  Double_t w = par[3];
  
  Double_t arg = 0.0; 
  
  if ( s!=0.0 ) arg = ( xx - Q )/s;    
  else cout << "Error: The code tries to divide by zero." << endl;

  Double_t gn = 0.5*( 1.0+TMath::Erf( Q/( sqrt(2.0)*s ) ) );
  
  Double_t result = 0.0;
  if ( xx>=0.0 ) result = w*1.0/lambda*TMath::Exp( -xx/lambda ) + ( 1.0-w )/( sqrt( 2.0*TMath::Pi() )*s*gn )*TMath::Exp( -0.5*arg*arg );
    
  return result;
  
}

SPEResponse::SPEResponse()
{}

SPEResponse::~SPEResponse()
{}

SPEResponse::SPEResponse( PMType::Response _spetype, Double_t _params[] )
{
  spetype = _spetype;
    
  
  if ( spetype==PMType::GAUSS )
    {
      params[0] = _params[0];
      params[1] = _params[1];
      
      params[2] = _params[2];
      
      params[3] = _params[3];
      
      spefunc = new TF1( "spefunc", _gausexpfunc, params[0]-80.0*params[1], params[0]+80.0*params[1], 4 );
      spefunc->SetParameters( params[0], params[1], params[2], params[3] );

      nparams = 4;
      
    }
   
   spefunc->SetLineColor( kBlue );
   spefunc->SetLineWidth( 2.0 );
   spefunc->SetNpx( 10000 );
    
}

void SPEResponse::SetParams( Double_t _params[] )
{
  for ( Int_t i=0; i<nparams; i++ )
    {
      params[i] = _params[i];
      spefunc->SetParameter( i, params[i] );
      
    }
      
}

Double_t SPEResponse::GetValue( Double_t xx )
{
  Double_t result = spefunc->Eval( xx );

  return result;
  
}


Double_t SPEResponse::GenQ()
{
  Double_t _x = spefunc->GetRandom();

  return _x;
  
}
