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


/** \class GenericFilter
 *
 *
 */

#include "modules/GenericFilter.h"

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

GenericFilter::GenericFilter()
{
}

//------------------------------------------------------------------------------

GenericFilter::~GenericFilter()
{
}

//------------------------------------------------------------------------------

void GenericFilter::Init()
{
  // PT threshold
  fElectronPTMin = GetDouble("ElectronPTMin", 5.0);
  fMuonPTMin = GetDouble("MuonPTMin", 5.0);
  fTauPTMin = GetDouble("TauPTMin", 5.0);
  
  // import input array(s)

  fInputArray = ImportArray(GetString("InputArray", "EFlowMerger/eflow"));
  fItInputArray = fInputArray->MakeIterator();

  fElectronInputArray = ImportArray(GetString("ElectronInputArray", "ElectronFilter/electrons"));
  fItElectronInputArray = fElectronInputArray->MakeIterator();

  fMuonInputArray = ImportArray(GetString("MuonInputArray", "MuonMomentumSmearing/muons"));
  fItMuonInputArray = fMuonInputArray->MakeIterator();

  fTauInputArray = ImportArray(GetString("TauInputArray", "TauFilter/taus"));
  fItTauInputArray = fTauInputArray->MakeIterator();

  fOutputArray = ExportArray(GetString("OutputArray", "eflow"));
}

//------------------------------------------------------------------------------

void GenericFilter::Finish()
{
  if(fItInputArray) delete fItInputArray;
  if(fItElectronInputArray) delete fItElectronInputArray;
  if(fItMuonInputArray) delete fItMuonInputArray;
  if(fItTauInputArray) delete fItTauInputArray;
}

//------------------------------------------------------------------------------

void GenericFilter::Process()
{
  Candidate *candidate, *electron, *muon, *jet;
  TLorentzVector v_clus, v_ele, v_muo, v_jet;

  // loop over input objects
  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    v_clus = candidate->Momentum;

    bool overlap = false;
    fItElectronInputArray->Reset();
    while((electron = static_cast<Candidate*>(fItElectronInputArray->Next()))) {
      v_ele = electron->Momentum;
      if(v_ele.Pt()>fElectronPTMin && v_clus.DeltaR(v_ele)<0.2) {
	overlap = true;
	break;
      }
    }

    if(!overlap) {
      fItMuonInputArray->Reset();
      while((muon = static_cast<Candidate*>(fItMuonInputArray->Next()))) {
	v_muo = muon->Momentum;
	if(v_muo.Pt()>fMuonPTMin && v_clus.DeltaR(v_muo)<0.2) {
	  overlap = true;
	  break;
	}
      }

      if(!overlap) {
	fItTauInputArray->Reset();
	while((jet = static_cast<Candidate*>(fItTauInputArray->Next()))) {
	  if(jet->TauTag==1) {
	    v_jet = jet->Momentum;
	    if(v_jet.Pt()>fTauPTMin && v_clus.DeltaR(v_jet)<0.4) {
	      overlap = true;
	      break;
	    }
	  }
	}
      }
    }

    if(!overlap) {
      fOutputArray->Add(candidate);
    }
  }
}

//------------------------------------------------------------------------------
