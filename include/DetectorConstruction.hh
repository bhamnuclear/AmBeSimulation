#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include <G4VUserDetectorConstruction.hh>
#include <G4GenericMessenger.hh>
#include <G4Tubs.hh>

class G4LogicalVolume;
class G4GenericMessenger;
class G4Material;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction(G4String IsotopeString);
  virtual ~DetectorConstruction();

  void ConstructSDandField() override;

  virtual G4VPhysicalVolume* Construct() override;
  G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

  void DetectorMessenger();
  G4Tubs* GetAmBeSolid(){return AbsorberSolid;};

protected:
    G4LogicalVolume*  fScoringVolume = nullptr;

private:
  G4GenericMessenger* fMessenger = nullptr;

  G4bool fDebug;
  G4String fRadioIsotope = "241Am";

  G4Tubs* ContainSolid = nullptr;
  G4LogicalVolume* ContainLog = nullptr;
  G4Tubs* AbsorberSolid = nullptr;
  G4LogicalVolume* AbsorberLog = nullptr;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
