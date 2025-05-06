#ifndef NNDCLoader_HH
#define NNDCLoader_HH

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <random>

#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>

class NNDCLoader
{
public:
  NNDCLoader(G4String filename);

  G4double GetEnergy(G4double prob);
  G4double GetEnergy();

private:
  G4int fDBG;

  std::vector<G4double> energy, probability;
  std::vector<G4double> energyUnsorted, probabilityUnsorted;

  void loadData(G4String filename);
  void SortData();
};

#endif
