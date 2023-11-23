#pragma once

//#include "n4-volumes.hh"
#include <n4-all.hh>
#include "G4Step.hh"
#include <G4ThreeVector.hh>
#include "colina.hh"

G4PVPlacement* colina_detector(G4ThreeVector& init_pos, params& params);
n4::sensitive_detector* sensitive(G4ThreeVector& init_pos, unsigned int sensorID, G4String sensor_name);
