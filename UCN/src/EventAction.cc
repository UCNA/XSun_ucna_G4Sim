#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
using   namespace       std;

#define	OUTPUT_FILE	"EnergyOutput.txt"

EventAction::EventAction()
: G4UserEventAction()
{}


EventAction::~EventAction()
{}


void EventAction::BeginOfEventAction(const G4Event* evt)
{
  fEdep_East_Scint = 0;	// Ensuring these values are reset.
  fEdep_West_Scint = 0;
  fEdep_East_MWPC = 0;
  fEdep_West_MWPC = 0;

  if((evt->GetEventID())%100 == 0)
  {
    G4cout << "/n -------------- Begin of event: " << evt->GetEventID() << G4endl;
  }
}


void EventAction::EndOfEventAction(const G4Event* evt)
{
  ofstream outfile;
  outfile.open(OUTPUT_FILE, ios::app);
  outfile << fEdep_East_Scint/keV << "\t \t" << fEdep_East_MWPC/keV << "\t \t"
	  << fEdep_West_Scint/keV << "\t \t" << fEdep_West_MWPC/keV << "\n";
  outfile.close();
}

// typeFlag = 0 -> Scint
//	      1 -> MWPC
void EventAction::AddEdep(G4double edep, int typeFlag, int locFlag)
{
  if(locFlag == 0)
  {
    if(typeFlag == 0)
    {
      fEdep_East_Scint = fEdep_East_Scint + edep;
    }
    else if(typeFlag == 1)
    {
      fEdep_East_MWPC = fEdep_East_MWPC + edep;
    }
  }
  else if(locFlag == 1)
  {
    if(typeFlag == 0)
    {
      fEdep_West_Scint = fEdep_West_Scint + edep;
    }
    else if(typeFlag == 1)
    {
      fEdep_West_MWPC = fEdep_West_MWPC + edep;
    }
  }
}
