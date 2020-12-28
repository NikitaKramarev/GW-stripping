import gw_stripping
import pytest
import numpy
import pandas
import astropy.constants
import os.path

@pytest.mark.parametrize('q, res', [
    # Next tests from Eggleton Paper
    (9.9901e-4, 0.0484),
    (9.901e-3, 0.1020),
    (0.0909, 0.2068),
    (0.2857, 0.3031),
    (0.5000, 0.3789),
    (0.7143, 0.4599),
    (0.9091, 0.5782),
    (0.9901, 0.7203),
    ])
def test_roche_nonrelativistic(q, res):
    f, _ = gw_stripping.roche_nonrelativistic(q)
    numpy.testing.assert_allclose(f, res, atol=0.0001)


def get_relativistic_asymtotic_data():
    M = [1, 1.5, 2.5, 3, 4]
    a = [50, 60, 70, 80, 90]
    q = [0.0909, 0.2857, 0.5000, 0.7143, 0.9091]
    return zip(q, M, a)


def get_relativistic_asymtotic_data_2():
    au = astropy.constants.au.to_value("km")
    return map(lambda x: (x[0], x[1], x[2] * au), get_relativistic_asymtotic_data())

@pytest.mark.parametrize('q, M, a', get_relativistic_asymtotic_data())
def test_roche_relativistic_asymtotic(q, M, a):
    res, _ = gw_stripping.roche_nonrelativistic(q)
    f, _, _ = gw_stripping.roche_relativistic(q, M, a)
    numpy.testing.assert_allclose(f, res, atol=0.10)

@pytest.mark.parametrize('q, M, a', get_relativistic_asymtotic_data_2())
def test_roche_relativistic_asymtotic_2(q, M, a):
    res, _ = gw_stripping.roche_nonrelativistic(q)
    f, _, _ = gw_stripping.roche_relativistic(q, M, a)
    numpy.testing.assert_allclose(f, res, atol=0.00001)


def get_mass_data(path, mode):
    df = pandas.read_table(path, sep=" ")
    eta = df["eta"].tolist()
    r = df["radius"].tolist()
    m = df["mass"].tolist()

    if (mode == "r2"):
        return map(lambda eta, r, m: (r, "r2", eta, m), eta, r, m)

    if (mode == "m2"):
        return map(lambda eta, r, m: (m, "m2", eta, r), eta, r, m)

    raise Exception("unknown mode: ", mode)


def local_file(path):
    return os.path.join(os.path.dirname(__file__), path)


@pytest.mark.parametrize('inp, m2_or_r2, eta, res', get_mass_data(local_file("mass.txt"), "r2"))
def test_mass_r2(inp, m2_or_r2, eta, res):
    f, _, _ = gw_stripping.radius_to_mass(inp, eta)
    numpy.testing.assert_allclose(f, res, atol=0.07)


@pytest.mark.parametrize('inp, m2_or_r2, eta, res', get_mass_data(local_file("mass.txt"), "m2"))
def test_mass_m2(inp, m2_or_r2, eta, res):
    f, _, _ = gw_stripping.mass_to_radius(inp, eta)
    numpy.testing.assert_allclose(f, res, rtol=0.1)
