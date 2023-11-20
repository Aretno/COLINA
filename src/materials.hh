#pragma once

#include <n4-all.hh>
//#include "nain4.hh"
//#include "g4-mandatory.hh"

#include <G4SystemOfUnits.hh>   // physical units such as `m` for metre
#include "G4Material.hh"
#include "cmath"
#include <CLHEP/Units/PhysicalConstants.h>

double XenonRefractiveIndex(double energy, double density);
double LXeScintillation(double energy);

G4Material* quartz_with_props();
G4Material* steel_with_props();
G4Material* fakegrid_with_props(G4Material* model_mat, double grid_transparency, double grid_thickness);
G4Material* lxe_with_props();
G4Material* CopyMaterial(G4Material* original, const G4String& newname);   