/** \class ZHUtils
 *
 *  \author R.-J. Wang
 */


#include "CMGTools/HtoZZ2l2nu/interface/ZHUtils.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRandom.h"

using namespace std;

namespace ZHUtils {

double Collins_Soper(const LorentzVector& lepton1, const LorentzVector& lepton2)
{

    TLorentzVector LVlep1(lepton1.Px(), lepton1.Py(), lepton1.Pz(), lepton1.E());
    TLorentzVector LVlep2(lepton2.Px(), lepton2.Py(), lepton2.Pz(), lepton2.E());
    TLorentzVector LVZ = LVlep1+LVlep2;

    //do transformation to Collins-Soper frame (1 rotation, 2 boosts)

    // 1st transormation - rotate to ptZ direction
    double zrot = -LVZ.Phi();
    LVlep1.RotateZ(zrot);
    LVlep2.RotateZ(zrot);

    // 2nd transformation - boost in z direction
    double beta_boostz = -LVZ.Pz()/LVZ.E();
    LVlep1.Boost(0.,0.,beta_boostz);
    LVlep2.Boost(0.,0.,beta_boostz);


    // 3rd transformation: boost in transverse direction (x-prime)
    double beta_boostx = -(LVlep1.Px()+LVlep2.Px())/(LVlep1.E()+LVlep2.E());
    LVlep1.Boost(beta_boostx,0.,0.);
    LVlep2.Boost(beta_boostx,0.,0.);

    // compute cos(theta*) in Colin-Soper frame
    double cos_theta_CS = LVlep1.CosTheta();
    if (LVZ.Pz() < 0) {
        cos_theta_CS *= -1.;
    }

    return cos_theta_CS;

}


float weightNLOEWKsignal(float pt)
{
    if(pt < 50) return 1;
    return 0.94-(0.2-0.068)/400.*(pt-50.);
}


float weightNLOEWKzz(float pt)
{
    if(pt < 50) return 1;
    return 1.+(-0.071*pt+0.55)/100.; 
}

float weightNLOEWKwz(float pt)
{
    if(pt < 50) return 1;
    return 1.+(-0.037*pt+1.9)/100.;
}


/*
float weightNLOEWKzz(float pt)
{
    if(pt < 50) return 1;
    double p_0 = -1.28637;
    double p_1 = -0.0967502;
    double p_2 = 5.38621e-05;
    double p_3 = -1.54371e-08;
    double a = p_0+p_1*pt+p_2*pt*pt+p_3*pow(pt,3);
    return 1.+a/100.; 
}


float weightNLOEWKwplusz(float pt)
{
    if(pt < 50) return 1;
    double p_0 = 1.31664;
    double p_1 = -0.0517186;
    double p_2 = -7.28883e-06;
    double p_3 = 2.5116e-08;
    double a = p_0+p_1*pt+p_2*pt*pt+p_3*pow(pt,3);
    return 1.+a/100.;
}

float weightNLOEWKwminuz(float pt)
{
    if(pt < 50) return 1;
    double p_0 = 1.54202;
    double p_1 = -0.0497632;
    double p_2 = -1.26063e-05;
    double p_3 = 2.81725e-08;
    double a = p_0+p_1*pt+p_2*pt*pt+p_3*pow(pt,3);
    return 1.+a/100.;
}
*/
}
