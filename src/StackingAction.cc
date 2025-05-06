#include "StackingAction.hh"
#include "RunAction.hh"

#include <G4SystemOfUnits.hh>

StackingAction::StackingAction(RunAction* runAction, EventAction* eventAction) :
  G4UserStackingAction(),
  fRunAction(runAction),
  fEventAction(eventAction)
{
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack (const G4Track*
 aTrack)
{
  /* -- from extended/electromagnetic/TestEm1
  //keep primary particle
  if (track->GetParentID() == 0) return fUrgent;

    // Task 4a.1: If the track has energy < 1 MeV, return fKill

  //energy spectrum of secondaries
  G4double energy = aTrack->GetKineticEnergy();
  G4double charge = aTrack->GetDefinition()->GetPDGCharge();

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if (charge != 0.) analysisManager->FillH1(5,energy);
  else              analysisManager->FillH1(6,energy);
  return fUrgent;
  //*/

  // Register only secondaries, i.e. tracks having ParentID > 0
  /*
  if (aTrack->GetParentID()>0)
  {
    fEventAction->AddSecondary(aTrack->GetParticleDefinition(), aTrack->GetKineticEnergy());
  }
  */
  // Do not affect track classification. Just return what would have
  // been returned by the base class
  return G4UserStackingAction::ClassifyNewTrack(aTrack);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
