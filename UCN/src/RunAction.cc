#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
using   namespace       std;

#define	OUTPUT_FILE	"FinalSim_EnergyOutput.txt"

RunAction::RunAction()
: G4UserRunAction()
{ }


RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run* run)
{
  ofstream outfile;
  outfile.open(OUTPUT_FILE, ios::app);
  outfile << "Particle species \t Momentum Direction: x \t y \t z \t Initial placement: x \t y \t z \t Energy Deposited (keV): East Scint \t East MWPC \t West Scint \t West MWPC \n";
  outfile.close();

  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}


void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }

  ofstream outfile;
  outfile.open(OUTPUT_FILE, ios::app);
  outfile << "Total number of simulated events during this run: " << run -> GetNumberOfEvent() << "\n";
  outfile.close();

}
