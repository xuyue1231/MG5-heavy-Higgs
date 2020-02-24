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

#ifndef GenericFilter_h
#define GenericFilter_h

/** \class GenericFilter
 *
 *
 */

#include "classes/DelphesModule.h"

#include <vector>
#include <map>

class TIterator;
class TObjArray;

class GenericFilter: public DelphesModule
{
public:

  GenericFilter();
  ~GenericFilter();

  void Init();
  void Process();
  void Finish();

private:

  Double_t fElectronPTMin; //!
  Double_t fMuonPTMin; //!
  Double_t fTauPTMin; //!
  
  const TObjArray *fInputArray; //!
  TIterator *fItInputArray; //!
  
  const TObjArray *fElectronInputArray; //!
  TIterator *fItElectronInputArray; //!
  
  const TObjArray *fMuonInputArray; //!
  TIterator *fItMuonInputArray; //!
  
  const TObjArray *fTauInputArray; //!
  TIterator *fItTauInputArray; //!

  TObjArray *fOutputArray; //!

  ClassDef(GenericFilter, 1)
};

#endif
