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

//------------------------------------------------------------------------------

#ifndef GenericPdgFilter_h
#define GenericPdgFilter_h

/** \class GenericPdgFilter
 *
 *  Removes particles with specific PDG codes
 *
 *  \author M. Selvaggi
 *
 */

#include "classes/DelphesModule.h"
#include "classes/DelphesClasses.h"
#include <vector>

class TIterator;
class TObjArray;

class GenericPdgFilter: public DelphesModule
{
public:

  GenericPdgFilter();
  ~GenericPdgFilter();

  void Init();
  void Process();
  void Finish();

  Bool_t hasMother(Candidate *candidate, Int_t pdg);
  
private:

  Double_t fPTMin; //!
  Double_t fEtaMax; //!
  Bool_t fMother; //!

  std::vector<Int_t> fPdgCodes;

  const TObjArray *fInputAllArray; //!

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(GenericPdgFilter, 1)
};

#endif
