#pragma once

#include <n4-all.hh>
#include <G4ThreeVector.hh>

//void el_generator    (G4Event* event, G4ThreeVector& init_pos);
void el_generator(G4Event* event, G4ThreeVector& init_pos, G4double rad_i, G4double z_i);
void sphere_generator(G4Event* event, G4ThreeVector& init_pos);