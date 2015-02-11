//____________________________________________________________________________
/*!

\class    genie::DFRInteractionListGenerator

\brief    Concrete implementations of the InteractionListGeneratorI interface.
          Generates a list of all the interactions that can be generated by the 
          DFR EventGenerator.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

\created  Feb 15, 2008

\cpright  Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _DIFFRACTIVE_INTERACTION_LIST_GENERATOR_H_
#define _DIFFRACTIVE_INTERACTION_LIST_GENERATOR_H_

#include "EVGCore/InteractionListGeneratorI.h"

namespace genie {

class DFRInteractionListGenerator : public InteractionListGeneratorI {

public :
  DFRInteractionListGenerator();
  DFRInteractionListGenerator(string config);
 ~DFRInteractionListGenerator();

  //-- implement the InteractionListGeneratorI interface
  InteractionList * CreateInteractionList(const InitialState & init) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void LoadConfigData(void);

  bool fIsCC;
  bool fIsNC;
};

}      // genie namespace
#endif // _DIFFRACTIVE_INTERACTION_LIST_GENERATOR_H_
