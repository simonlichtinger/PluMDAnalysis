from filecmp import cmp
import os
import sys

import pytest
import MDAnalysis

sys.path.append(os.path.abspath('../PluMDAnalysis'))
import PluMDAnalysis
from PluMDAnalysis import PluMDAnalysis
from PluMDAnalysis.aux_classes import PLUMED_Distance, PLUMED_Group, PLUMED_Restraint

# load test data just once
u = MDAnalysis.Universe("test/equil.gro", "test/test.xtc")

def test_consistency_checks():
    # Check for correct handling of various consistency errors
    with pytest.raises(TypeError):
        pmda = PluMDAnalysis(1)
    pmda = PluMDAnalysis(u)
    helices = [
        (u.select_atoms('resid 68:73 and name CA'),
        u.select_atoms('resid 46:51 and name CA')),
        (u.select_atoms('resid 79:84 and name CA'),
        u.select_atoms('resid 98:103 and name CA'))
    ]
    with pytest.raises(TypeError):
        pmda.add_atom_group(1)
    with pytest.raises(ValueError):
        pmda.add_atom_group(helices[0][0], group_type='ABC')
    dist = PLUMED_Distance(helices[0][0], helices[1][0], "test")
    rest = PLUMED_Restraint([dist])
    with pytest.raises(RuntimeError): rest.get_plumed_str()
    with pytest.raises(RuntimeError): rest.add_distance_values([1])
    with pytest.raises(ValueError): rest.add_time_series([1,2],[[1],[1],[1]])
    with pytest.raises(ValueError): rest.add_time_series([1,2],[[1,2],[1]])
    with pytest.raises(ValueError): rest.add_time_series([1,2],[1,2])
    with pytest.raises(RuntimeError): rest.determine_distance_values(u,[(0,100), (4900,5000)])
    rest.add_time_series([1,2],[[1],[1]])
    with pytest.raises(ValueError): rest.add_distance_values([[1],[2],[3]])
    with pytest.raises(ValueError): rest.add_distance_values([1,1])
    with pytest.raises(ValueError): rest.add_distance_values([[1,2],[2]])
    with pytest.raises(ValueError): rest.determine_distance_values(u,[(0,100), (200,300), (4900,5000)])
    with pytest.raises(ValueError): rest.determine_distance_values(u,[(0,100), 4900])
    rest.add_distance_values([[1],[2]])
    
def test_simple_steer_manual():
    helices = [
        (u.select_atoms('resid 68:73 and name CA'),
        u.select_atoms('resid 46:51 and name CA')),
        (u.select_atoms('resid 79:84 and name CA'),
        u.select_atoms('resid 98:103 and name CA'))
    ]

    pmda = PluMDAnalysis(u)
    pmda.add_distance(helices[0][0], helices[1][0])
    pmda.add_distance(helices[0][1], helices[1][1])


    # This adds an atom group again, so it shouldn't have any effect! The assert later controls for that
    pmda.add_atom_group(helices[0][0])

    pmda.add_restraint_manual([0, 1000], [[0,0],[1000,1000]], [[5,2],[3,0.5]])
    pmda.generate_PLUMED_input()
    assert cmp("plumed.dat", "test/expected_simple_steer.dat")
    os.remove("plumed.dat")

def test_simple_steer_automatic():
    helices = [
        (u.select_atoms('resid 68:73 and name CA'),
        u.select_atoms('resid 46:51 and name CA')),
        (u.select_atoms('resid 79:84 and name CA'),
        u.select_atoms('resid 98:103 and name CA')),
        (u.select_atoms('resid 325:330 and name CA'),
        u.select_atoms('resid 290:295 and name CA')),
        (u.select_atoms('resid 341:346 and name CA'),
        u.select_atoms('resid 364:369 and name CA'))
    ]

    pmda = PluMDAnalysis(u)
    pmda.add_distance(helices[0][0], helices[2][0]) # corresponding to 1->7 and 2->8 tip distances
    pmda.add_distance(helices[1][0], helices[3][0])
    pmda.add_restraint_from_trajectory([0, 1000], [[0,0],[1000,1000]], [(0,100), (4900,5000)])
    pmda.generate_PLUMED_input()
    assert cmp("plumed.dat", "test/expected_simple_steer_automatic.dat")
    os.remove("plumed.dat")

#test_simple_steer_automatic()
