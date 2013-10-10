#ifndef ZHUtils_H
#define ZHUtils_H

/** \class ZHUtils
 *  No description available.
 *
 *  \author R.-J. Wang 
 */

#include <iostream>
#include <TString.h>
#include <map>

#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <utility>
#include <vector>

#include "Math/LorentzVector.h"
#include "TVector2.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;



namespace ZHUtils{

  double Collins_Soper(const LorentzVector& lepton1, const LorentzVector& lepton2); 
  float weightNLOEWKsignal(float pt);
  float weightNLOEWKzz(float pt);
  float weightNLOEWKwz(float pt);
  //float weightNLOEWKwplusz(float pt);
  //float weightNLOEWKwminuz(float pt); 

}
#endif

