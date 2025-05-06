#include "DetectorConstruction.hh"

#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4VisAttributes.hh>

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Orb.hh>
#include <G4UnionSolid.hh>
#include <G4SubtractionSolid.hh>

#include <G4SDManager.hh>
#include <G4SDParticleFilter.hh>

#include <G4MultiFunctionalDetector.hh>
#include <G4PSEnergyDeposit.hh>
#include <G4PSCellFlux.hh>
#include <G4PSDoseDeposit.hh>

#include <sstream>

//using namespace std;


DetectorConstruction::DetectorConstruction(G4String IsotopeString)
: G4VUserDetectorConstruction(),
  fScoringVolume(0),
  fRadioIsotope(IsotopeString)
{
  fDebug = false;

  DetectorMessenger();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4int ncomponents, natoms;
  G4bool checkOverlaps = true;
  G4NistManager* nist = G4NistManager::Instance();

  //--Elements
  G4Element* C = nist->FindOrBuildElement("C");
  G4Element* O = nist->FindOrBuildElement("O");
  G4Element* Si = nist->FindOrBuildElement("Si");
  G4Element* P = nist->FindOrBuildElement("P");
  G4Element* S = nist->FindOrBuildElement("S");
  G4Element* Cr = nist->FindOrBuildElement("Cr");
  G4Element* Mn = nist->FindOrBuildElement("Mn");
  G4Element* Fe = nist->FindOrBuildElement("Fe");
  G4Element* Ni = nist->FindOrBuildElement("Ni");
  G4Element* Mo = nist->FindOrBuildElement("Mo");


  //--Materials
  G4Material* matBe = nist->FindOrBuildMaterial("G4_Be");
  G4Material* matAir = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* matVacuum = nist->FindOrBuildMaterial("G4_Galactic");

  G4Material* Steel304 = new G4Material("Steel304", 8*g/cm3, ncomponents=3);
  Steel304->AddElement(Fe, 72*perCent);
  Steel304->AddElement(Cr, 18*perCent);
  Steel304->AddElement(Ni,  10*perCent);

  //Steel: 02X17H14M2
  //from https://www.sciencedirect.com/science/article/pii/S0022311504002521?via=ihub
  G4Material* matInox = new G4Material("Stainless-Steel", 7.917*g/cm3, ncomponents=9);
  matInox->AddElement(Fe, 64.653*perCent);
  matInox->AddElement(C, 0.02*perCent);
  matInox->AddElement(Si, 0.4*perCent);
  matInox->AddElement(Mn, 1.5*perCent);
  matInox->AddElement(P, 0.017*perCent);
  matInox->AddElement(S, 0.01*perCent);
  matInox->AddElement(Cr, 17.0*perCent);
  matInox->AddElement(Ni, 14.0*perCent);
  matInox->AddElement(Mo, 2.4*perCent);
  //*/


  //Am241
  G4Isotope* isoAm = new G4Isotope("isoAm", 95, 241, 241. * g/mole);
  G4Element* elAm = new G4Element("Am241", "Am241", 1);
  elAm->AddIsotope(isoAm, 100.*perCent);

  //Pu239
  G4Isotope* isoPu = new G4Isotope("isoAm", 94, 239, 239. * g/mole);
  G4Element* elPu = new G4Element("Pu239", "Pu239", 1);
  elPu->AddIsotope(isoPu, 100.*perCent);

  //AmO2
  G4Material* matAmO2 = new G4Material("AmO2", 11.68*g/cm3, ncomponents=2);
  matAmO2->AddElement(elAm,natoms=1);
  matAmO2->AddElement(O,natoms=2);

  //PuO2
  G4Material* matPuO2 = new G4Material("PuO2", 11.5*g/cm3, ncomponents=2);
  matPuO2->AddElement(elPu,natoms=1);
  matPuO2->AddElement(O,natoms=2);

  //AmBe
  G4double PerBe = .890326;
  G4double PerAmO2 = 1-PerBe;
  G4double densityAmBe = (1.848*PerBe+11.68*PerAmO2);//2.9265704g/cm3
  G4double densityPuBe = (1.848*PerBe+11.5*PerAmO2);//2.9265704g/cm3

  G4Material* matAmBe = new G4Material("AmBe", densityAmBe*g/cm3, ncomponents=2);
  matAmBe->AddMaterial(matBe, PerBe*100*perCent);
  if (fRadioIsotope=="241Am")
    matAmBe->AddMaterial(matAmO2, PerAmO2*100*perCent);
  else if (fRadioIsotope=="239Pu")
    matAmBe->AddMaterial(matPuO2, PerAmO2*100*perCent);

  //matAmBe->AddMaterial(matVacuum, PerAmO2*100*perCent);//TEST ONLY


  //--Colours
  //red
  G4VisAttributes* red = new G4VisAttributes(G4Colour::Red());
  red->SetVisibility(true);
  //red->SetForceSolid(true);
  //yellow
  G4VisAttributes* yellow = new G4VisAttributes(G4Colour::Yellow());
  yellow->SetVisibility(true);
  //yellow->SetForceSolid(true);

  //--Measurements
  //absorber
  G4double AbsorberRadius = 17.6*mm/2.;
  //G4double AbsorberRadius = 17.6*mm/2.;
  G4double AbsorberLength = 17.5*mm;
  G4double ContainThickness = 2.4*mm;

  /*
  //absorber -- RUS/6137/S-96 IBN-241-15-4
  G4double AbsorberRadius = 17.4*mm/2.;
  G4double AbsorberLength = 17.6*mm;
  G4double ContainThicknessR = 2.54*mm;
  G4double ContainThicknessL = 2.4*mm;
  //*/


  //--World
  G4double WorldSizeXY = 50.*cm *3.*10;
  G4double WorldSizeZ  = 51.*cm *3.*10;
  // compute dimensions
  G4double ContainRadius = AbsorberRadius + ContainThickness;
  G4double ContainLength = AbsorberLength + 2*ContainThickness;



  //-- GEOMETRY
  G4VisAttributes* visAttr = new G4VisAttributes();
  visAttr->SetVisibility(true);

  //--World
  G4Tubs* WorldSolid = new G4Tubs("WorldSolid", 0., WorldSizeXY, WorldSizeZ, 0., 360.*deg);
  //G4LogicalVolume* WorldLog = new G4LogicalVolume(WorldSolid, matAir, "WorldLogical");
  G4LogicalVolume* WorldLog = new G4LogicalVolume(WorldSolid, matVacuum, "WorldLogical");//TestOnly
  WorldLog->SetVisAttributes(visAttr);
  G4VPhysicalVolume* worldPhys = new G4PVPlacement(0, G4ThreeVector(), WorldLog, "World", 0, false, 0);

  //Place source
  //Container (Steel)
  ContainSolid = new G4Tubs("ContainerSolid", 0., ContainRadius, 0.5*ContainLength, 0., 360.*deg);
  ContainLog = new G4LogicalVolume(ContainSolid, matInox, "ContainerLogical");
  ContainLog->SetVisAttributes(yellow);
  new G4PVPlacement(0, G4ThreeVector(), ContainLog, "Container", WorldLog, false, 0, checkOverlaps);
  //Absorber (BeO)
  AbsorberSolid = new G4Tubs("AmBeSolid", 0., AbsorberRadius, 0.5*AbsorberLength, 0., 360.*deg);
  AbsorberLog = new G4LogicalVolume(AbsorberSolid, matAmBe, "AmBeLogical");
  AbsorberLog->SetVisAttributes(red);
  new G4PVPlacement(0, G4ThreeVector(), AbsorberLog, "AmBe", ContainLog, false, 0, checkOverlaps);


  //Show the material table
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  // The Construct() method has to return the final (physical) world volume:
  return worldPhys;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fMessenger, fScoringVolume;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DetectorMessenger()
{
  fMessenger = new G4GenericMessenger(this, "/AmBe/detector/", "Geometry control");
  fMessenger->DeclareProperty("Debug",fDebug);
  fMessenger->DeclareProperty("RadioIsotope",fRadioIsotope);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
