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


/** \class AngImpSmearing
 *
 *  Performs transverse angular resolution smearing.
 *
 *  \author M. Selvaggi - UCL, Louvain-la-Neuve
 *
 */

#include "modules/AngImpSmearing.h"

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

AngImpSmearing::AngImpSmearing() :
  fFormulaEta(0), fFormulaPhi(0), fFormula(0), fItInputArray(0)
{
  fFormulaEta = new DelphesFormula;
  fFormulaPhi = new DelphesFormula;
  fFormula = new DelphesFormula;
}

//------------------------------------------------------------------------------

AngImpSmearing::~AngImpSmearing()
{
  if(fFormulaEta) delete fFormulaEta;
  if(fFormulaPhi) delete fFormulaPhi;
  if(fFormula) delete fFormula;
}

//------------------------------------------------------------------------------

void AngImpSmearing::Init()
{
  // read resolution formula

  fFormulaEta->Compile(GetString("EtaResolutionFormula", "0.0"));
  fFormulaPhi->Compile(GetString("PhiResolutionFormula", "0.0"));
  fFormula->Compile(GetString("ImpactParResolutionFormula", "0.0"));

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "TrackMerger/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void AngImpSmearing::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void AngImpSmearing::Process()
{
  Candidate *candidate, *particle, *mother;
  Double_t pt, px, py, eta, phi, e;
  Double_t xd, yd, zd, d0, sx, sy, sz, dd0;

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->Position;
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    eta = candidatePosition.Eta();
    phi = candidatePosition.Phi();
    pt = candidateMomentum.Pt();
    e = candidateMomentum.E();
    
    px = candidateMomentum.Px();
    py = candidateMomentum.Py();

    if(pt <= 0.0) continue;

    // apply smearing formula for eta,phi

    Double_t deta = gRandom->Gaus(0, fFormulaEta->Eval(pt, eta, phi, e));
    Double_t dphi = gRandom->Gaus(0, fFormulaPhi->Eval(pt, eta, phi, e));
    
    // calculate coordinates of closest approach to track circle in transverse plane xd, yd, zd
    xd =  candidate->Xd;
    yd =  candidate->Yd;
    zd =  candidate->Zd;

    // calculate smeared values
    sx = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));
    sy = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));
    sz = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));

    xd += sx;
    yd += sy;
    zd += sz;

    // calculate impact parameter (after-smearing)
    d0 = (xd*py - yd*px)/pt;

    dd0 = gRandom->Gaus(0.0, fFormula->Eval(pt, eta, phi, e));

    // fill smeared values in candidate
    mother = candidate;
    
    candidate = static_cast<Candidate*>(candidate->Clone());
    eta = candidateMomentum.Eta();
    phi = candidateMomentum.Phi();
    candidate->Momentum.SetPtEtaPhiE(pt, eta+deta, phi+dphi, pt*TMath::CosH(eta+deta));
    candidate->Xd = xd;
    candidate->Yd = yd;
    candidate->Zd = zd;

    candidate->D0 = d0;
    candidate->ErrorD0 = dd0;
    
    candidate->AddCandidate(mother);    
    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
