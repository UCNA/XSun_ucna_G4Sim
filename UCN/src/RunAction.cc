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

#define	OUTPUT_FILE	"TrackerTestOutput.txt"

RunAction::RunAction()
: G4UserRunAction()
{ }


RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run* run)
{
  fKillCount = 0;
  ofstream outfile;
  outfile.open(OUTPUT_FILE, ios::app);
/*  outfile << "Event ID \t Species \t Init KE (keV) \t xPos (m) \t yPos (m) \t zPos (m) \t xMomentum \t yMo \t zMo \t time (ns) \t weight \t"
	<< "Trapped? \t Comp Time (s) \t Energy Deposited (keV): East Scint \t East MWPC \t West Scint \t West MWPC \n";
*/
  outfile << "Event ID \t Init KE (keV) \t"
	<< "Trapped? \t Edep (keV): East Scint \t East MWPC \t West Scint \t West MWPC \n";
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
     << "--------------------End of Global Run----------------------- \n";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------ \n";
  }

  G4cout << "Number of trapped events killed: " << fKillCount << G4endl;

  ofstream outfile;
  outfile.open(OUTPUT_FILE, ios::app);
  outfile << "Total number events: " << run -> GetNumberOfEvent()
	  << " || Final kill count: " << fKillCount << "\n";
  outfile.close();

}
