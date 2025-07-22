import astropy.units as u
from typing import Optional, Union, List
import os
os.environ['GALFIND_CONFIG_NAME'] = 'galfind_config_Griley.ini'

from galfind.Data import morgan_version_to_dir
from galfind import Catalogue, EAZY, SED_code
from galfind import galfind_logger, Redshift_Bin_Selector, Redshift_Extractor, Multiple_Mask_Selector, Multiple_SED_fit_Selector, Min_Instrument_Unmasked_Band_Selector, Unmasked_Band_Selector, Bluewards_LyLim_Non_Detect_Selector, Bluewards_Lya_Non_Detect_Selector, Redwards_Lya_Detect_Selector, Chi_Sq_Lim_Selector, Chi_Sq_Diff_Selector, Robust_zPDF_Selector, Sextractor_Bands_Radius_Selector    
from galfind import Bagpipes

# Define config (same as used in your main2 call)
survey = "JADES-DR3-GS-East"
version = "v13"
instrument_names = ["ACS_WFC", "NIRCam"]
aper_diams = [0.32] * u.arcsec
forced_phot_band = ["F277W", "F356W", "F444W"]

# Define the SED fitter used
eazy = EAZY({"templates": "fsps_larson", "lowz_zmax": None})
bagpipes = Bagpipes({
    "fix_z": eazy,
    "sfh": "continuity_bursty",
    "z_calculator": None,
    "sps_model": "BPASS",
})

# Load the catalogue (with SEDs already saved from previous main run)
cat = Catalogue.pipeline(
    survey,
    version,
    instrument_names=instrument_names,
    aper_diams=aper_diams,
    forced_phot_band=forced_phot_band,
    version_to_dir_dict=morgan_version_to_dir
)

# --- Extract SED for galaxy 13553 ---
galaxy_id = 13553
aper_diam = aper_diams[0]

# Get the stored SED dictionary for this galaxy and this fitter
sed = cat.get_SED(galaxy_id, aper_diam, bagpipes)

# Extract SED components
wavelength = sed["wavelength"]           # in Angstroms
spectrum = sed["spectrum_total"]         # in erg/s/cm^2/Å
phot_wave = sed["phot_wavelength"]       # filter effective wavelengths in Angstroms
phot_obs = sed["photometry_obs"]         # observed fluxes in μJy
phot_err = sed["photometry_err"]         # observed errors in μJy
phot_model = sed["photometry_total"]     # model photometry in μJy

# --- Plot ---
plt.figure(figsize=(10, 6), facecolor='white')

# Plot model SED (convert Å → μm)
plt.plot(wavelength / 1e4, spectrum, label='Model SED (Bagpipes)', color='navy')

# Observed photometry with error bars
plt.errorbar(
    phot_wave / 1e4, phot_obs, yerr=phot_err,
    fmt='o', color='black', label='Observed photometry', capsize=3
)

# Model photometry points
plt.scatter(phot_wave / 1e4, phot_model, color='crimson', label='Model photometry', zorder=3)

# Final formatting
plt.xlabel("Wavelength [μm]")
plt.ylabel("Flux [μJy]")
plt.title(f"Galaxy {galaxy_id}: SED from Bagpipes fit")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.4)
plt.tight_layout()
plt.show()
