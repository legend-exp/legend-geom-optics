"""
VM2000 reflective film inside water tank.

.. [Geis2017] Ch. Geis et al. 2017 “Optical response of highly reflective film used in the water Cherenkov muon veto of
    the XENON1T dark matter experiment” In: Journal of Instrumentation, Vol. 12, 2017, https://dx.doi.org/10.1088/1748-0221/12/06/P06017

WLS properties taken from MaGe.
"""

from __future__ import annotations

import logging
import math

import numpy as np
import pint
from pint import Quantity

from legendoptics import store
from legendoptics.utils import InterpolatingGraph, readdatafile

log = logging.getLogger(__name__)
u = pint.get_application_registry()


def vm2000_calculate_wls_mfp(yield_value):
    # Set total path length (currently hardcoded to 1 mm)
    total_path = 1.0e-3 * u.m

    # Handle edge cases
    if yield_value == 0:
        return 10.0  # 10 m,  Large mean free path, no absorption
    if yield_value == 1:
        return 0.01e-3  # 0.01 mm, Very small mean free path, 100% absorption

    # Calculate mean free path for valid yield values
    help_value = math.log(1.0 - yield_value)

    return -total_path / help_value


@store.register_pluggable
def vm2000_refractive_index() -> float:
    """Refractive index.

    .. optics-const::
    """
    return 1.15


@store.register_pluggable
def vm2000_absorption_length() -> float:
    """Absorption length.

    .. optics-const::
    """
    return 50.0


@store.register_pluggable
def vm2000_particle_scintillationyields() -> (
    tuple[Quantity, Quantity, Quantity, Quantity, Quantity, Quantity]
):
    """Set all particle Scintillations yields to 0 (except of electron scintillation yield).

    .. optics-plot::
    """

    scintillation_params_energy = [1.0, 10000000.0]
    vm2000_protonscintillation_yield = [0.0, 0.0]
    vm2000_alphascintillation_yield = [0.0, 0.0]
    vm2000_deuteronscintillation_yield = [0.0, 0.0]
    vm2000_ionscintillation_yield = [0.0, 0.0]
    vm2000_tritonscintillation_yield = [0.0, 0.0]

    return (
        scintillation_params_energy,
        vm2000_protonscintillation_yield,
        vm2000_alphascintillation_yield,
        vm2000_deuteronscintillation_yield,
        vm2000_ionscintillation_yield,
        vm2000_tritonscintillation_yield,
    )


@store.register_pluggable
@u.with_context("sp")
def vm2000_parameters() -> tuple[Quantity, Quantity, Quantity, Quantity, Quantity]:
    """Wavelength-shifting parameters for the reflective foil VM2000."""
    # Constants
    wls_yield = 0.075  # 0.6 MaGe, 0.075 XENON paper

    # Populate VM2000_energy_range array with energy values
    ppsci_high_e = (115 * u.nm).to("eV")
    ppsc_low_e = (650 * u.nm).to("eV")

    num1 = 251
    dee = (ppsci_high_e - ppsc_low_e) / (num1 - 2)
    vm2000_energy_range = np.zeros(num1) * u.eV
    for ji in range(1, num1):
        vm2000_energy_range[ji] = ppsc_low_e + ji * dee
        vm2000_energy_range[0] = 1.8 * u.eV

    # Create arrays for energy and optical properties
    vm2000_reflectivity = np.zeros(num1)
    vm2000_efficiency = np.zeros(num1)
    wls_absorption = np.zeros(num1) * u.m

    # Set reflectivity, absorption, and emission
    for ji in range(num1):
        if vm2000_energy_range[ji] < (370 * u.nm).to(
            "eV"
        ):  # 370 nm < (related to energy)
            vm2000_reflectivity[ji] = 0.95  # Visible light 0.95, 0.99
        else:
            vm2000_reflectivity[ji] = 0.12  # UV light 0.15, 0.3 (paper)

        if vm2000_energy_range[ji] > 3.35 * u.eV:  # 5 eV 3.35
            # depending on path length in foil --> angle
            wls_absorption[ji] = vm2000_calculate_wls_mfp(wls_yield)  # Absorbs UV
        else:
            wls_absorption[ji] = (
                1.0 * u.m
            )  # Imperturbed, no absorption of visible light

    g = InterpolatingGraph(
        *readdatafile("vm2000_em_spec.dat", pkg="l200geom.materials"), zero_outside=True
    )
    wls_emission = g(vm2000_energy_range.to("nm")).to("dimensionless")

    # Copy the first element to 0th position
    wls_absorption[0] = wls_absorption[1]  # depending on path length in foil --> angle
    wls_emission[0] = wls_emission[1]

    return (
        vm2000_energy_range,
        vm2000_reflectivity,
        vm2000_efficiency,
        wls_absorption,
        wls_emission,
    )


@store.register_pluggable
def vm2000_scint_timeconstant() -> Quantity:
    """Time constant, from MaGe.

    .. optics-const::
    """
    return 0.5 * u.ns


def pyg4_vm2000_attach_rindex(mat, reg) -> None:
    """Attach the refractive index to the given VM2000 material instance.

    See Also
    --------
    .vm2000_refractive_index
    """
    λ = np.array([100, 600]) * u.nm
    r = [vm2000_refractive_index()] * 2

    with u.context("sp"):
        mat.addVecPropertyPint("RINDEX", λ.to("eV"), r)


def pyg4_vm2000_attach_absorption_length(mat, reg) -> None:
    """Attach the refractive index to the given VM2000 material instance.

    See Also
    --------
    .vm2000_absorption_length
    """
    λ = np.array([100, 600]) * u.nm
    r = [vm2000_absorption_length()] * 2 * u.m

    with u.context("sp"):
        mat.addVecPropertyPint("ABSLENGTH", λ.to("eV"), r)


def pyg4_vm2000_attach_particle_scintillationyields(mat, reg) -> None:
    """Attach the scintillation yiels (except of electron yield) to the given VM2000 material instance.

    See Also
    --------
    .vm2000_particle_scintillationyields
    """
    energy = vm2000_particle_scintillationyields()[0] * u.eV
    proton_yield = vm2000_particle_scintillationyields()[1]
    alpha_yield = vm2000_particle_scintillationyields()[2]
    deuteron_yield = vm2000_particle_scintillationyields()[3]
    ion_yield = vm2000_particle_scintillationyields()[4]
    triton_yield = vm2000_particle_scintillationyields()[5]

    mat.addVecPropertyPint("PROTONSCINTILLATIONYIELD", energy, proton_yield)
    mat.addVecPropertyPint("ALPHASCINTILLATIONYIELD", energy, alpha_yield)
    mat.addVecPropertyPint("DEUTERONSCINTILLATIONYIELD", energy, deuteron_yield)
    mat.addVecPropertyPint("IONSCINTILLATIONYIELD", energy, ion_yield)
    mat.addVecPropertyPint("TRITONSCINTILLATIONYIELD", energy, triton_yield)


def pyg4_vm2000_attach_reflectivity(mat, reg) -> None:
    """Attach the reflectivity to the given VM2000 material instance.

    See Also
    --------
    .vm2000_parameters
    """
    energy = vm2000_parameters()[0]
    r = vm2000_parameters()[1]

    mat.addVecPropertyPint("REFLECTIVITY", energy, r)


def pyg4_vm2000_attach_efficiency(mat, reg) -> None:
    """Attach the efficiency to the given VM2000 material instance.

    See Also
    --------
    .vm2000_parameters
    """
    vm2000_energy_range, _, vm2000_efficiency, _, _ = vm2000_parameters()

    mat.addVecPropertyPint("EFFICIENCY", vm2000_energy_range, vm2000_efficiency)


def pyg4_vm2000_attach_wls(mat, reg) -> None:
    """Attach wavelength shifting properties to the given VM2000 material instance.

    See Also
    --------
    .vm2000_parameters
    .vm2000_scint_timeconstant
    """

    vm2000_energy_range, _, _, wls_absorption, wls_emission = vm2000_parameters()

    mat.addVecPropertyPint("WLSABSLENGTH", vm2000_energy_range, wls_absorption)
    mat.addVecPropertyPint("WLSCOMPONENT", vm2000_energy_range, wls_emission)
    mat.addConstPropertyPint("WLSTIMECONSTANT", vm2000_scint_timeconstant())


def pyg4_vm2000_attach_border_params(mat, reg) -> None:
    """Attach border parameters between water and the given VM2000 material instance.

    See Also
    --------
    .vm2000_parameters
    """
    vm2000_energy_range, vm2000_reflectivity, vm2000_efficiency, _, _ = (
        vm2000_parameters()
    )

    reflectivity_front = vm2000_reflectivity * 0
    efficiency_border = vm2000_efficiency * 0
    transmittance_border = [1.0] * len(vm2000_energy_range)

    mat.addVecPropertyPint("REFLECTIVITY", vm2000_energy_range, reflectivity_front)
    mat.addVecPropertyPint("EFFICIENCY", vm2000_energy_range, efficiency_border)
    mat.addVecPropertyPint("TRANSMITTANCE", vm2000_energy_range, transmittance_border)
