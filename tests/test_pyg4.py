"""Test the attaching of all properties to Geant4 materials."""

from __future__ import annotations

import pint
import pyg4ometry.geant4 as g4

u = pint.get_application_registry()


def _create_dummy_mat() -> tuple[g4.Registry, g4.Material]:
    reg = g4.Registry()
    mat = g4.MaterialSingleElement("dummy", 1, 1, 1, reg)
    return reg, mat


def test_pyg4_attach_lar() -> None:
    import legendoptics.lar

    u = pint.get_application_registry()
    reg, mat = _create_dummy_mat()
    legendoptics.lar.pyg4_lar_attach_rindex(mat, reg)
    legendoptics.lar.pyg4_lar_attach_attenuation(mat, reg, 90 * u.K)
    legendoptics.lar.pyg4_lar_attach_scintillation(mat, reg)


def test_pyg4_attach_tpb() -> None:
    import legendoptics.tpb

    reg, mat = _create_dummy_mat()
    legendoptics.tpb.pyg4_tpb_attach_rindex(mat, reg)
    legendoptics.tpb.pyg4_tpb_attach_wls(mat, reg)


def test_pyg4_attach_fibers() -> None:
    import legendoptics.fibers

    reg, mat = _create_dummy_mat()
    legendoptics.fibers.pyg4_fiber_cladding1_attach_rindex(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.fibers.pyg4_fiber_cladding1_attach_rindex(mat, reg)
    reg, mat = _create_dummy_mat()
    legendoptics.fibers.pyg4_fiber_core_attach_rindex(mat, reg)
    legendoptics.fibers.pyg4_fiber_core_attach_wls(mat, reg)
    legendoptics.fibers.pyg4_fiber_core_attach_absorption(mat, reg)
    legendoptics.fibers.pyg4_fiber_core_attach_scintillation(mat, reg)


def test_pyg4_attach_tetratex() -> None:
    import legendoptics.tetratex

    reg, mat = _create_dummy_mat()
    legendoptics.tetratex.pyg4_tetratex_attach_reflectivity(mat, reg)


def test_pyg4_attach_germanium() -> None:
    import legendoptics.germanium

    reg, mat = _create_dummy_mat()
    legendoptics.germanium.pyg4_germanium_attach_reflectivity(mat, reg)


def test_pyg4_attach_copper() -> None:
    import legendoptics.copper

    reg, mat = _create_dummy_mat()
    legendoptics.copper.pyg4_copper_attach_reflectivity(mat, reg)


def test_pyg4_attach_silicon() -> None:
    import legendoptics.silicon

    reg, mat = _create_dummy_mat()
    legendoptics.silicon.pyg4_silicon_attach_complex_rindex(mat, reg)


def test_pyg4_attach_nylon() -> None:
    import legendoptics.nylon

    reg, mat = _create_dummy_mat()
    legendoptics.nylon.pyg4_nylon_attach_rindex(mat, reg)
    legendoptics.nylon.pyg4_nylon_attach_absorption(mat, reg)


def test_pyg4_attach_pen() -> None:
    import legendoptics.pen

    reg, mat = _create_dummy_mat()
    legendoptics.pen.pyg4_pen_attach_rindex(mat, reg)
    legendoptics.pen.pyg4_pen_attach_attenuation(mat, reg)
    legendoptics.pen.pyg4_pen_attach_wls(mat, reg)
    legendoptics.pen.pyg4_pen_attach_scintillation(mat, reg)
