#include "materials.hh"
#include "material_props.hh"

double XenonRefractiveIndex(double energy, double density)
{
  // Formula for the refractive index taken from
  // A. Baldini et al., "Liquid Xe scintillation calorimetry
  // and Xe optical properties", arXiv:physics/0401072v1 [physics.ins-det]

  // The Lorentz-Lorenz equation (also known as Clausius-Mossotti equation)
  // relates the refractive index of a fluid with its density:
  // (n^2 - 1) / (n^2 + 2) = - A · d_M,     (1)
  // where n is the refractive index, d_M is the molar density and
  // A is the first refractivity viral coefficient:
  // A(E) = \sum_i^3 P_i / (E^2 - E_i^2),   (2)
  // with:
  double P[3] = {71.23, 77.75, 1384.89}; // [eV^3 cm3 / mole]
  double E[3] = {8.4, 8.81, 13.2};       // [eV]

  // Note.- Equation (1) has, actually, a sign difference with respect
  // to the one appearing in the reference. Otherwise, it yields values
  // for the refractive index below 1.

  // Let's calculate the virial coefficient.
  // We won't use the implicit system of units of Geant4 because
  // it results in loss of numerical precision.

  energy = energy / eV;
  double virial = 0.;

  for (unsigned int i=0; i<3; i++)
  virial = virial + P[i] / (energy*energy - E[i]*E[i]);

  // Need to use g/cm3
  density = density / g * cm3;

  double mol_density = density / 131.29;
  double alpha = virial * mol_density;

  // Isolating now the n2 from equation (1) and taking the square root
  double n2 = (1. - 2*alpha) / (1. + alpha);
  if (n2 < 1.) {
    //      G4String msg = "Non-physical refractive index for energy "
    // + bhep::to_string(energy) + " eV. Use n=1 instead.";
    //      G4Exception("[XenonProperties]", "RefractiveIndex()",
    // 	  JustWarning, msg);
    n2 = 1.;
  }

  return sqrt(n2);
}

double LXeScintillation(double energy)
  {
    // K. Fuji et al., "High accuracy measurement of the emission spectrum of liquid xenon
    // in the vacuum ultraviolet region",
    // Nuclear Instruments and Methods in Physics Research A 795 (2015) 293–297
    // http://ac.els-cdn.com/S016890021500724X/1-s2.0-S016890021500724X-main.pdf?_tid=83d56f0a-3aff-11e7-bf7d-00000aacb361&acdnat=1495025656_407067006589f99ae136ef18b8b35a04
    double Wavelength_peak = 174.8*nm;
    double Wavelength_FWHM = 10.2*nm;
    double Wavelength_sigma = Wavelength_FWHM/2.35;

    double Energy_peak = (CLHEP::h_Planck*CLHEP::c_light/Wavelength_peak);
    double Energy_sigma = (CLHEP::h_Planck*CLHEP::c_light*Wavelength_sigma/pow(Wavelength_peak,2));
    // double bin = 6*Energy_sigma/500;

    double intensity =
	  exp(-pow(Energy_peak/eV-energy/eV,2)/(2*pow(Energy_sigma/eV, 2)))/(Energy_sigma/eV*sqrt(CLHEP::pi*2.));

    return intensity;
  }

G4Material* CopyMaterial(G4Material* original, const G4String& newname) {
  auto newmat = n4::material(newname);  
    if (newmat == 0) {

      G4double density     = original->GetDensity();
      G4double temperature = original->GetTemperature();
      G4double pressure    = original->GetPressure();
      G4State  state       = original->GetState();
      G4int    n_elem      = original->GetNumberOfElements();

      if (n_elem == 1) {
        G4double z = original->GetZ();
        G4double a = original->GetA();
        newmat = new G4Material(newname, z, a, density, state, temperature, pressure);
      }
      else {
        const G4double* fractions = original->GetFractionVector();
        newmat = new G4Material(newname, density, n_elem, state, temperature, pressure);
        for (G4int i = 0; i < n_elem; ++i)
          newmat->AddElement(new G4Element(original->GetElement(i)->GetName(),
                                          original->GetElement(i)->GetSymbol(),
                                          original->GetElement(i)->GetZ(),
                                          original->GetElement(i)->GetA()),
                                          fractions[i]);
        }
    }

    return newmat;   
}

G4Material* quartz_with_props() {
  auto quartz = n4::material("G4_SILICON_DIOXIDE");  
  quartz -> SetMaterialPropertiesTable(quartz_props());
  return quartz;
}

G4Material* steel_with_props() {
  auto steel = n4::material("G4_STAINLESS-STEEL");  
  steel -> SetMaterialPropertiesTable(steel_props());
  return steel;
}

G4Material* lxe_with_props() {
  auto lxe = n4::material("G4_lXe");  
  lxe -> SetMaterialPropertiesTable(LXe_props());
  return lxe;
}

G4Material* fakegrid_with_props(G4Material* model_mat, double grid_transparency, double grid_thickness) {
  auto fake = CopyMaterial(model_mat, "Grid");//n4::material("G4_lXe");  
  fake -> SetMaterialPropertiesTable(fakegrid_props(grid_transparency, grid_thickness));
  return fake;
}

G4Material* TPB()
{
  G4String name = "TPB"; // Tetraphenyl butadiene
  G4Material* mat = G4Material::GetMaterial(name, false);
  if (mat == 0) {
    G4NistManager* nist = G4NistManager::Instance();
    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");
    mat = new G4Material(name, 1*g/cm3, 2, kStateSolid);
    mat->AddElement(H, 22);
    mat->AddElement(C, 28);
  }
  return mat;
}
