#include "Run.hh"

Run::Run()
: G4Run()
{}


Run::~Run()
{}


void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  G4cout << "Call to merge" << G4endl;

  G4Run::Merge(run);
}


