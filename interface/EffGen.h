#ifndef effgen_h
#define effgen_h

// Standard libraries
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

#include <ROOT/RVec.hxx>

#include "Base/Modules/interface/correctionWrapper.h"

typedef ROOT::VecOps::RVec<float> fRVec;
typedef ROOT::VecOps::RVec<int> iRVec;
typedef ROOT::VecOps::RVec<bool> bRVec;

class EffGen {
  public:
    EffGen (std::string filename, std::string tag);
    ~EffGen ();
    float get_weight(
        int displaced_pdgid, float input_ctau, float output_ctau,
		iRVec GenPart_pdgId, iRVec GenPart_genPartIdxMother, fRVec GenPart_mass,
		fRVec GenPart_pt, fRVec GenPart_vertex_x, fRVec GenPart_vertex_y, fRVec GenPart_vertex_z,
		fRVec GenPart_eta, fRVec GenPart_phi
    );
    float get_weight_noeff(
        int displaced_pdgid, float input_ctau, float output_ctau,
		iRVec GenPart_pdgId, iRVec GenPart_genPartIdxMother, fRVec GenPart_mass,
		fRVec GenPart_pt, fRVec GenPart_vertex_x, fRVec GenPart_vertex_y, fRVec GenPart_vertex_z,
		fRVec GenPart_eta, fRVec GenPart_phi
    );

  private:
    MyCorrections corr;

};


#endif
