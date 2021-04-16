[![Run unit tests with OS coverage](https://github.com/simonlichtinger/PluMDAnalysis/actions/workflows/unit_tests_os.yml/badge.svg?branch=main)](https://github.com/simonlichtinger/PluMDAnalysis/actions/workflows/unit_tests_os.yml)
[![codecov](https://codecov.io/gh/simonlichtinger/PluMDAnalysis/branch/main/graph/badge.svg?token=OS2KKHONFI)](https://codecov.io/gh/simonlichtinger/PluMDAnalysis)

# PluMDAnalysis

This library automates the generation of PLUMED input files by providing an
interface with python trajectory analysis via MDAnalysis. The current focus is
on steered MD via distances between atom groups, but I will try to keep it 
extensible and will do so once the occasion arises.

So far, minimal functional version for making a simple steered MD. 
