#include "generator.hh"

#include <G4RandomDirection.hh> // for launching particles in random directions
#include <G4Event.hh>           // needed to inject primary particles into an event
#include <G4ThreeVector.hh>

//#include "nain4.hh"
//#include "g4-mandatory.hh"
//#include "n4_ui.hh"
//#include "n4-utils.hh"
//#include "n4-volumes.hh"

void el_generator(G4Event* event, G4ThreeVector& init_pos, G4double rad_i, G4double z_i) {
  auto particle_type = n4::find_particle("opticalphoton");
  auto particle = new G4PrimaryParticle(particle_type);
  auto p = G4RandomDirection();
  auto energy = 7.085 * eV; // 7.085 eV for 175 nm
  
  double el_rad  = rad_i;//40.5*mm; // same as in geometry
  double sphere_rad = 20.*cm;
  double anode_z = 0. *mm;
  double r      = el_rad * sqrt(G4UniformRand());
  double phi    = G4UniformRand() * 2*M_PI;

  init_pos = G4ThreeVector(r * cos(phi),r * sin(phi), z_i);
  particle->SetMomentumDirection(p);
  particle->SetKineticEnergy(energy);

  auto vertex   = new G4PrimaryVertex(init_pos, 0);
  vertex -> SetPrimary(particle);
  event  -> AddPrimaryVertex(vertex);
}

void sphere_generator(G4Event* event, G4ThreeVector& init_pos) {
  auto particle_type = n4::find_particle("opticalphoton");
  auto particle = new G4PrimaryParticle(particle_type);
  auto p = G4RandomDirection();
  auto energy = 7.085 * eV; // 7.085 eV for 175 nm
  
  double target_rad  = 20.*cm;
  double r           = target_rad * sqrt(G4UniformRand());

  double u = G4UniformRand();
  double v = G4UniformRand();

  double theta = - u * M_PI;
  double phi   = acos(2.0 * v - 1.0);

  double x = r * sin(phi) * cos(theta);
  double y = r * sin(phi) * sin(theta);
  double z = r * cos(phi);

  init_pos = G4ThreeVector(x, y, z);
  particle->SetMomentumDirection(p);
  particle->SetKineticEnergy(energy);

  auto vertex   = new G4PrimaryVertex(init_pos, 0);
  vertex -> SetPrimary(particle);
  event  -> AddPrimaryVertex(vertex);
}