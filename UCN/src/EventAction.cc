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

/*TrackerHitsCollection* EventAction::GetHitsCollection(int hcID, const G4Event* event) const
{
  TrackerHitsCollection* hitsCollection = static_cast<TrackerHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

  return hitsCollection;
}*/

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  fTrapped = false;
  fStartTime = 0;

  // Sets the start of the C++ 'clock' used for tracking trapped ptcl's
  fStartTime = clock();

  if((evt->GetEventID())%1000 == 0)
  {
    G4cout << "\n -------------- Begin of event: " << evt->GetEventID() << "\n" << G4endl;
  }

  // get hits collection IDs only once in an overall run though.
  if(fHitsCollectionIDs[0] == -1)
  {
    for(int i = 0; i < fNbSDs; i++)
    {
      fHitsCollectionIDs[i] = G4SDManager::GetSDMpointer()->GetCollectionID((*fMyDetectorConstruction).fHCNamesArray[i]);
      G4cout << "HC ID: " << fHitsCollectionIDs[i] << ", belongs to HC name: " << (*fMyDetectorConstruction).fHCNamesArray[i] << G4endl;
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

  G4HCofThisEvent* hce = evt->GetHCofThisEvent();
  if(!hce)
  {
    G4ExceptionDescription msg;
    msg << "No hits collection of this event found." << G4endl;
    G4Exception("EventAction::EndOfEventAction()", "Error", JustWarning, msg);
    return;
  }


  TrackerHitsCollection* SD_totalHC[fNbSDs];
  TrackerHit* SD_hitsFullStep[fNbSDs];
  TrackerHit* SD_hitsTrack[fNbSDs];
  G4double SD_edepStepTrack[fNbSDs];
  G4double SD_edepSimpleTrack[fNbSDs];
  G4double SD_edepFullTrack[fNbSDs];

  for(int j = 0; j < fNbSDs; j++)
  {
    SD_edepStepTrack[j] = 0;
    SD_edepSimpleTrack[j] = 0;
    SD_edepFullTrack[j] = 0;
  }

  G4cout << "No seg fault" << G4endl;

  for(int i = 0; i < fNbSDs; i++)
  {
    SD_totalHC[i] = static_cast<TrackerHitsCollection*>(hce->GetHC(fHitsCollectionIDs[i]));

    if(SD_totalHC[i]->GetSize() == 0)
    {
      continue;
    }
    SD_hitsFullStep[i] = (*SD_totalHC[i])[0];
    SD_edepStepTrack[i] = SD_hitsFullStep[i] -> GetStepEdep();

    for(int t = 1; t < (SD_totalHC[i] -> GetSize()); t++)
    {
    SD_hitsTrack[i] = (*SD_totalHC[i])[t];
    SD_edepSimpleTrack[i] += SD_hitsTrack[i] -> GetTestEdep();
    SD_edepFullTrack[i] += SD_hitsTrack[i] -> GetEdep();
    }
    G4cout << "Index " << i << " has HC of size = " << SD_totalHC[i] -> GetSize() << G4endl;

  }

  // print out of class member variables due to stepping action
  ofstream outfile;
  outfile.open(OUTPUT_FILE, ios::app);
  outfile << fTrapped << "\t"
	  << compTime << "\n"
	  << "\t" << SD_edepStepTrack[fScintEast_index]/keV << "\t"
	  << "\t" << SD_edepStepTrack[fActiveWireVolEast_index]/keV << "\t"
	  << "\t" << SD_edepStepTrack[fScintWest_index]/keV << "\t"
	  << "\t" << SD_edepStepTrack[fActiveWireVolWest_index]/keV << "\n"
          << "\t" << SD_edepSimpleTrack[fScintEast_index]/keV << "\t"
          << "\t" << SD_edepSimpleTrack[fActiveWireVolEast_index]/keV << "\t"
          << "\t" << SD_edepSimpleTrack[fScintWest_index]/keV << "\t"
          << "\t" << SD_edepSimpleTrack[fActiveWireVolWest_index]/keV << "\n"
          << "\t" << SD_edepFullTrack[fScintEast_index]/keV << "\t"
          << "\t" << SD_edepFullTrack[fActiveWireVolEast_index]/keV << "\t"
          << "\t" << SD_edepFullTrack[fScintWest_index]/keV << "\t"
          << "\t" << SD_edepFullTrack[fActiveWireVolWest_index]/keV << "\n";

  outfile.close();
}

// locFlag = 0 -> EAST
// 	     1 -> WEST
// typeFlag = 0 -> Scint
//	      1 -> MWPC

