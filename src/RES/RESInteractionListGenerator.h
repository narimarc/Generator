//____________________________________________________________________________
/*!

\class    genie::RESInteractionListGenerator

\brief    Creates a list of all the interactions that can be generated by the
          RES thread (generates semi-inclusive resonance reactions).
          Concrete implementations of the InteractionListGeneratorI interface.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  May 13, 2005

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _RES_INTERACTION_LIST_GENERATOR_H_
#define _RES_INTERACTION_LIST_GENERATOR_H_

#include "BaryonResonance/BaryonResList.h"
#include "EVGCore/InteractionListGeneratorI.h"

namespace genie {

class RESInteractionListGenerator : public InteractionListGeneratorI {

public :
  RESInteractionListGenerator();
  RESInteractionListGenerator(string config);
 ~RESInteractionListGenerator();

  // implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfigData(void);

  bool          fIsCC;
  bool          fIsNC;
  bool          fIsEM;
  BaryonResList fResList;
};

}      // genie namespace
#endif // _RES_INTERACTION_LIST_GENERATOR_H_

