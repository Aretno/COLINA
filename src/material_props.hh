#pragma once

#include "G4MaterialPropertiesTable.hh"
#include <G4SystemOfUnits.hh>   // physical units such as `m` for metre

#include "materials.hh"
#include <n4-all.hh>

const G4double optPhotMinE_   =  0.2  * eV;
const G4double optPhotMaxE_   = 11.5  * eV;
const G4double noAbsLength_   = 1.e8  * m;

G4MaterialPropertiesTable* PTFE_props();
G4MaterialPropertiesTable* PTFE_props(double refl);
G4MaterialPropertiesTable* Al2O3_mirror_props(double refl);
G4MaterialPropertiesTable* Sensor_mirror_props(double refl, double trans);
G4MaterialPropertiesTable* LXe_props();
G4MaterialPropertiesTable* quartz_props();
G4MaterialPropertiesTable* steel_props();
G4MaterialPropertiesTable* tpb_props();
G4MaterialPropertiesTable* fakegrid_props(G4double transparency,
                                          G4double thickness);
                                          