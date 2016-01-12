#include "EventAction.hh"
#include "RunAction.hh"
#include "TrackerHit.hh"
#include "TrackerSD.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
using   namespace       std;

#define	OUTPUT_FILE	"TrackerTestOutput.txt"

EventAction::EventAction()
: G4UserEventAction(), fStartTime(0)
{
  fMyDetectorConstruction = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  // initialize all values to -1 since that doesn't correspond to any actual SD

  for(int i = 0; i < fNbSDs; i++)	// see comment in EventAction.hh on fNbSDs
  {
    fHitsCollectionIDs[i] = -1;
  }

  fScintEast_index = -1;
  fScintWest_index = -1;
  fActiveWireVolEast_index = -1;
  fActiveWireVolWest_index = -1;
}


EventAction::~EventAction()
{}

TrackerHitsCollection* EventAction::GetHitsCollection(int hcID, const G4Event* event) const
{
  TrackerHitsCollection* hitsCollection = static_cast<TrackerHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

  return hitsCollection;
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  fTrapped = false;
  fEdep_East_Scint = 0;	// Ensuring these values are reset.
  fEdep_West_Scint = 0;
  fEdep_East_MWPC = 0;
  fEdep_West_MWPC = 0;
  fStartTime = 0;

  // Sets the start of the C++ 'clock' used for tracking trapped ptcl's
  fStartTime = clock();

  if((evt->GetEventID())%1000 == 0)
  {
    G4cout << "\n -------------- Begin of event: " << evt->GetEventID() << "\n" << G4endl;
  }
}


void EventAction::EndOfEventAction(const G4Event* evt)
{
  if(fTrapped == true)
  {
    G4cout << "Event " << evt->GetEventID() << " trapped. Incrementing kill count." << G4endl;
    // some code that M. Mendenhall came up with that allows me to access my own RunAction class
    ((RunAction*)G4RunManager::GetRunManager()->GetUserRunAction()) -> IncrementKillCount();
  }

  // store how much time the event ran for
  clock_t timeOfEvent = clock() - fStartTime;
  double compTime = ((double)timeOfEvent)/CLOCKS_PER_SEC;

  // get hits collection IDs only once in an overall run though.
  if(fHitsCollectionIDs[0] == -1)
  {
    for(int i = 0; i < fNbSDs; i++)
    {
      fHitsCollectionIDs[i] = G4SDManager::GetSDMpointer()->GetCollectionID((*fMyDetectorConstruction).fHCNamesArray[i]);
//      G4cout << "HC ID: " << fHitsCollectionIDs[i] << ", belongs to HC name: " << (*fMyDetectorConstruction).fHCNamesArray[i] << G4endl;
      // these are the special ones we want to record extra info from.
      if((*fMyDetectorConstruction).fHCNamesArray[i] == "HC_scint_EAST")
      {
        fScintEast_index = i;
      }
      else if((*fMyDetectorConstruction).fHCNamesArray[i] == "HC_scint_WEST")
      {
        fScintWest_index = i;
      }
      else if((*fMyDetectorConstruction).fHCNamesArray[i] == "HC_wireVol_EAST")
      {
        fActiveWireVolEast_index = i;
      }
      else if((*fMyDetectorConstruction).fHCNamesArray[i] == "HC_wireVol_WEST")
      {
        fActiveWireVolWest_index = i;
      }
    }
  }

  TrackerHitsCollection* SD_totalHC[fNbSDs];
  TrackerHit* SD_hits[fNbSDs];
  G4double SD_edep[fNbSDs];
  for(int i = 0; i < fNbSDs; i++)
  {
    SD_totalHC[i] = GetHitsCollection(fHitsCollectionIDs[i], evt);
    SD_hits[i] = (*SD_totalHC[i])[0];	// 0th entry is the hit item in HC we are using to record everything
    SD_edep[i] = SD_hits[i] -> GetEdep();
  }


  // print out of class member variables due to stepping action
  ofstream outfile;
  outfile.open(OUTPUT_FILE, ios::app);
  outfile << fTrapped << "\t"
	  << compTime << "\t"
	  << fEdep_East_Scint/keV << "\t" << SD_edep[fScintEast_index]/keV << "\t"
	  << fEdep_East_MWPC/keV << "\t" << SD_edep[fActiveWireVolEast_index]/keV << "\t"
	  << fEdep_West_Scint/keV << "\t" << SD_edep[fScintWest_index]/keV << "\t"
	  << fEdep_West_MWPC/keV << "\t" << SD_edep[fActiveWireVolWest_index]/keV << "\n";
  outfile.close();
}

// locFlag = 0 -> EAST
// 	     1 -> WEST
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
