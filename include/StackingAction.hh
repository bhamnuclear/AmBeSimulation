#ifndef STACKINGACTION_HH
#define STACKINGACTION_HH

#include <G4UserStackingAction.hh>
#include <G4Track.hh>

#include "EventAction.hh"

class RunAction;
class EventAction;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class StackingAction : public G4UserStackingAction
{
public:
  //! constructor
  StackingAction(RunAction*,EventAction*);

  //! destructor
  ~StackingAction(){;};

  //! Main interface
  G4ClassificationOfNewTrack   ClassifyNewTrack (const G4Track*);

private:
  RunAction* fRunAction;
  EventAction* fEventAction;

};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
