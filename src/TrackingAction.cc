#include "TrackingAction.hh"

#include <G4StepStatus.hh>
#include <G4IonTable.hh>
#include <G4TrackingManager.hh>

TrackingAction::TrackingAction(EventAction* EventAct, bool FissFragments):
    fEventAction(EventAct),
    fFissFragments(FissFragments)
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)//when track is created
{
  if (aTrack->GetParentID()>0)
  {
    if (fFissFragments)
    {
      if(aTrack->GetParticleDefinition()->GetParticleName()=="neutron" and aTrack->GetCreatorProcess()->GetProcessName()=="nFissionHP")
        fEventAction->AddFissNeut(aTrack->GetKineticEnergy()/MeV);

      if(aTrack->GetParticleDefinition()->IsGeneralIon() and aTrack->GetCreatorProcess()->GetProcessName()=="nFissionHP")
        fEventAction->AddFissIon(aTrack->GetParticleDefinition()->GetAtomicMass());
    }

    if(aTrack->GetParticleDefinition()->GetParticleName()=="neutron" and aTrack->GetCreatorProcess()->GetProcessName()!="nFissionHP")
    {
      fEventAction->AddSecondary(aTrack->GetParticleDefinition(), aTrack->GetKineticEnergy());
    }

    if(aTrack->GetParticleDefinition()->GetParticleName()!="neutron")
    {
      fEventAction->AddSecondary(aTrack->GetParticleDefinition(), aTrack->GetKineticEnergy());
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)//when track is killed
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
