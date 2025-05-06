#include "NNDCLoader.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
NNDCLoader::NNDCLoader(G4String filename)
{
  fDBG = 0;
  loadData(filename);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double NNDCLoader::GetEnergy(G4double prob)
{
  //Get probability (prob) between 0 and 1
  for (unsigned int i=0; i<=probability.size(); i++)
  {
    if (prob < probability[i])
      return energy[i];
  }
  //G4cout << "NNDCLoader\tprobability:"<<prob<<G4endl;//DBG
  G4cout << "NNDCLoader::GetEnergy\tNo valid energy found"<<G4endl;
  return -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double NNDCLoader::GetEnergy()
{
  return GetEnergy(G4UniformRand());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NNDCLoader::loadData(G4String filename)
{
  //just to be safe
  energy.clear();
  probability.clear();

  //open file
  std::ifstream file(filename);
  if (!file)
  {
    G4cerr << "NNDCLoader: Error opening file: " << filename << G4endl;
    return;
  }

  G4double energy_, intensity_;//, energyUnc_, intensityUnc_;
  G4double totalIntensity = 0.0;
  std::vector<G4double> intensity;

  //while (file >> energy_ >> energyUnc_ >> intensity_ >> intensityUnc_)
  while (file >> energy_ >> intensity_)
  {
      energyUnsorted.push_back(energy_*keV);
      intensity.push_back(intensity_);
      totalIntensity += intensity_;
  }
  file.close();

  // Normalize intensities to sum to 100%, intensities are already in percent
  //normalise to 1
  for (G4double intensity_: intensity) {
      probabilityUnsorted.push_back((intensity_/ totalIntensity));
  }
  SortData();

  if (energy.size()!=probability.size())
  {
    G4cerr << "NNDCLoader: energy and probability have different size" << G4endl;
  }
  else if (energy.empty() || probability.empty())
  {
    G4cerr << "NNDCLoader: No data available for sampling." << G4endl;
    return;
  }

  if (fDBG)//Debug
  {
    G4cout << "NNDCLoader: Dumping energy and probabilities\n";
    for (unsigned int i=0; i<energy.size(); i++)
    {
      G4cout << energy[i]<<"\t"<<probability[i]<<"\n";
    }
    G4cout << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NNDCLoader::SortData()
{
  //sort by ascending intensity, from lowest to highest
  std::vector<size_t> indices(energyUnsorted.size());
  for (size_t i = 0; i < indices.size(); ++i)
  {
    indices[i] = i;
  }
  std::sort(indices.begin(), indices.end(), [this](size_t a, size_t b)
  {
    return probabilityUnsorted[a] < probabilityUnsorted[b];
  });

  for (size_t i : indices) {
    energy.push_back(energyUnsorted[i]);
    probability.push_back(probabilityUnsorted[i]);
  }

  for (unsigned int i=0; i<probability.size(); i++)
  {
    if (i>0) probability[i]+=probability[i-1];
  }
}
