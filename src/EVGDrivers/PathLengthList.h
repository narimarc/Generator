//____________________________________________________________________________
/*!

\class   genie::PathLengthList

\brief   Object to be filled with the neutrino path-length, for all detector
         geometry materials, when starting from a position x and travelling
         along the direction of the neutrino 4-momentum.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 24, 2005

*/
//____________________________________________________________________________

#ifndef _PATH_LENGTH_LIST_H_
#define _PATH_LENGTH_LIST_H_

#include <map>
#include <ostream>

class TLorentzVector;

using std::map;
using std::ostream;

namespace genie {

class PDGCodeList;

class PathLengthList : public map<int, double> {

public :

  PathLengthList(const PDGCodeList & pdglist);
  ~PathLengthList();

  void   AddPathLength  (int pdgc, double pl); // path-legth(pdgc) += pl
  void   SetPathLength  (int pdgc, double pl); // path-legth(pdgc)  = pl
  void   SetAllToZero   (void);
  double PathLength     (int pdgc) const;
  void   Print          (ostream & stream) const;
  
  friend ostream & operator << (ostream & stream, const PathLengthList & list);
};

}      // genie namespace

#endif // _PATH_LENGTH_LIST_H_
