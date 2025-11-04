#!/usr/bin/env python3
from json.encoder import JSONEncoder
from json.decoder import JSONDecoder

from qiskit_serverless import get_arguments, save_result

import pyscf
from pyscf import gto, scf
from pyscf.solvent import pcm
from pyscf.mcscf import avas 

### Argument retrieval
args                      = get_arguments() 
data                      = args["data"]      # Chemistry Data

i_data = JSONDecoder().decode(data) # What should we decode here?

[mol_geo, eps, ao_labels] = i_data

print(f">>>>> DEFINING MOLECULE")
# You must complete this coding section

print(f">>>>> BUILDING MOLECULE")
# You must complete this coding section


print(f">>>>> DEFINING PCM")
# You must complete this coding section


print(f">>>>> BUILDING RESTRICTED HARTREE FOCK")
# You must complete this coding section


print(f">>>>> RUNNING AVAS")
# You must complete this coding section


print(f">>>>> STARTING CASCI")
# You must complete this coding section

CASCI_E = None
print(f">>>>> CASCI_E: {CASCI_E}")

o_data = JSONEncoder().encode([CASCI_E])

# JSON-safe package
save_result({"outputs": o_data})  # single JSON blob returned to client

