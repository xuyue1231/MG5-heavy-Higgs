/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class GenericPdgFilter
 *
 *  Removes particles with specific PDG codes
 *
 *  \author M. Selvaggi
 *
 */

#include "modules/GenericPdgFilter.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

GenericPdgFilter::GenericPdgFilter() :
  fItInputArray(0)
{
}

//------------------------------------------------------------------------------

GenericPdgFilter::~GenericPdgFilter()
{
}

//------------------------------------------------------------------------------

void GenericPdgFilter::Init()
{

  ExRootConfParam param;
  Size_t i, size;

  // PT threshold
  fPTMin = GetDouble("PTMin", 0.0);
  fEtaMax = GetDouble("EtaMax", 999.0);

  // import input array
  fInputAllArray = ImportArray(GetString("InputAllArray", "Delphes/allParticles"));

  fInputArray = ImportArray(GetString("InputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  param = GetParam("PdgCode");
  size = param.GetSize();

  fMother = GetBool("RequireMotherPdg", false);

  // read PdgCodes to be filtered out from the data card

  fPdgCodes.clear();
  for(i = 0; i < size; ++i)
  {
    fPdgCodes.push_back(param[i].GetInt());
  }

  // create output array
  fOutputArray = ExportArray(GetString("OutputArray", "filteredParticles"));
}

//------------------------------------------------------------------------------

void GenericPdgFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

Bool_t GenericPdgFilter::hasMother(Candidate *candidate, Int_t pdg) {
  if(candidate->PID==pdg) return true;
  if(candidate->M1>=0) {
    for(Int_t i=candidate->M1; i<=TMath::Max(candidate->M2,0); i++) {
      Candidate *cand = static_cast<Candidate*>(fInputAllArray->At(i));
      if(cand->PID==pdg) return true;
    }
    Candidate *cand = static_cast<Candidate*>(fInputAllArray->At(candidate->M1));
    return hasMother(cand,pdg);
  }
  return false;
}

void GenericPdgFilter::Process()
{
  Candidate *candidate;
  Int_t pdgCode;
  Bool_t pass;
  Double_t pt;
  Double_t eta;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    pdgCode = candidate->PID;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    pt = candidateMomentum.Pt();
    eta = candidateMomentum.Eta();

    pass = kTRUE;
    if(fMother) {
      for(std::vector<Int_t>::size_type j=0; j<fPdgCodes.size(); j++) {
	if(hasMother(candidate,fPdgCodes[j])) {
	  pass = kFALSE; break;
	}
      }
    }
    else {
      if(find(fPdgCodes.begin(), fPdgCodes.end(), pdgCode) != fPdgCodes.end()) pass = kFALSE;
    }

    if(pt < fPTMin || fabs(eta) >= fEtaMax) pass = kFALSE;

    if(pass) fOutputArray->Add(candidate);
  }
}
