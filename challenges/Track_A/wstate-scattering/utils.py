import cmath
import networkx as nx
import numpy as np
from math import sqrt, pi
from typing import Tuple, List
from functools import lru_cache
from qiskit import QuantumCircuit
from sympy import exp, simplify
from qiskit.transpiler import CouplingMap

def __truncate_amplitudes__(cs, phis, xzero, d, L):
    assert d % 2 == 1, "d must be odd"
    assert d >= 3, "d must be >=3"
    eta = int((d-1)/2)
    assert xzero - eta >= 0
    assert xzero + eta < L
    
    cs_truncated = cs[xzero - eta:xzero + eta + 1]
    truncation_norm = (cs_truncated**2).sum()
    cs_truncated /= np.sqrt(truncation_norm)
    
    return cs_truncated, phis[xzero - eta:xzero + eta + 1]

def __calc_thetas__(cs, d):
    assert d % 2 == 1, "d must be odd"
    assert d >= 3, "d must be >=3"
    assert len(cs) == d
    assert np.abs((cs**2).sum() - 1.) < 1e-9

    eta = int((d-1)/2)
    thetas = np.zeros((d))

    # central theta (index eta)
    thetas[eta] = 2. * np.arcsin(np.sqrt(((cs[eta:])**2).sum()))
    
    # theta at eta - 1 not needed (and not defined by recursion relation)
    thetas[eta-1] = None

    # thetas above eta:
    for j in range(eta):
        sine_prod = 1.
        for i in range(j+1):
            sine_prod *= np.sin(thetas[eta+i]/2.)
        thetas[eta + j + 1] = 2. * np.arccos(cs[eta+j]/sine_prod)

    # thetas below eta:
    coseta = np.cos(thetas[eta]/2.)
    for j in range(1, eta):
        sine_prod = 1.
        for i in range(2, j+1):
            sine_prod *= np.sin(thetas[eta-i]/2.)
        thetas[eta - j - 1] = 2. * np.arccos(cs[eta - j]/sine_prod/coseta)

    return thetas

def __prep_w_state__(qc, d, phis, thetas):
    assert d % 2 == 1, "d must be odd"
    assert d >= 3, "d must be >=3"
    assert qc.num_qubits == d
    
    eta = int((d-1)/2)
    qc.ry(thetas[eta], eta)
    qc.cry(thetas[eta+1], eta, eta+1)
    qc.cx(eta, eta-1, ctrl_state=0)
    qc.cx(eta+1, eta)
    for j in range(eta+2, d):
        qc.cry(thetas[j], j-1, j)
        qc.cx(j, j-1)

    for j in range(eta-2, -1, -1):
        qc.cry(thetas[j], j+1, j)
        qc.cx(j, j+1)

    for q in range(d):
        qc.rz(phis[q], q)

    return qc

def __prep_zero_w_state__(qc, d, phis, thetas):
    assert d % 2 == 1, "d must be odd"
    assert d >= 3, "d must be >=3"
    assert qc.num_qubits == d
    
    eta = int((d-1)/2)
    # "...with the angle of the first RY-gate set to zero..."
    qc.ry(0, eta)
    qc.cry(thetas[eta+1], eta, eta+1)
    # "...and the second CNOT changed from a control on |0⟩ to a control on |1⟩."
    qc.cx(eta, eta-1, ctrl_state=1)
    qc.cx(eta+1, eta)
    for j in range(eta+2, d):
        qc.cry(thetas[j], j-1, j)
        qc.cx(j, j-1)

    for j in range(eta-2, -1, -1):
        qc.cry(thetas[j], j+1, j)
        qc.cx(j, j+1)

    for q in range(d):
        qc.rz(phis[q], q)

    return qc

def __compute_normalization__(mprime:int, kzero:float, sigma_sq_k:float, L:int)-> float:
    m = mprime - L/2
    k = 2*pi*m/L
    norm_exp = ((kzero-k)**2)/(2*sigma_sq_k)
    return exp(-norm_exp)

def __construct_wavepacket__(n:int, mprime:int, xzero:float, sigma_sq_k: float, kzero:float, L:int, norm:float)-> complex:
    m = mprime - L/2
    k =2*pi*m/L
    exponent = ((kzero - k)**2)/(4*sigma_sq_k) 
    phase = (1j*k)*(n - xzero)    
    summand = norm/sqrt(L) * exp(phase) * exp(-exponent) # positive time evolution
    #summand = norm/sqrt(L) * exp(-phase) * exp(-exponent) # bakes-in the reversed momentum for negative time evolution
    return summand

def __return_params__(n:int, xzero:int, sigma_sq_k:float, kzero:float, L:int, norm:float)-> Tuple[float, float]:
    amplitude = 0.
    for mprime in range(1, L+1):
        summand = simplify(__construct_wavepacket__(n, mprime, xzero, sigma_sq_k, kzero, L, norm))
        amplitude += summand
    mag, phase = cmath.polar(amplitude)
    return mag, phase

def __W_Y__(circ, L, theta): 
    q = 0 
    for q in range(L): 
        circ.ry(theta, q)
        
    return circ

def __Ryz__(theta): 
    qc = QuantumCircuit(2, name='Ryz')
    qc.s(1)
    qc.cz(0, 1)
    qc.ry(theta, 0)
    qc.rx(-theta, 1)
    qc.cz(0, 1)
    qc.sdg(1)
    gate = qc.to_gate(label=fr'$e^{{-i \theta (YZ + ZY) / 2}}$')
    return gate
    

def __W_YZ__(circ, L, theta): 
    q = 0 
    while q+1 < L: 
        circ.append(__Ryz__(theta), [q, q+1])
        q += 2

    q = 1
    while q+1 < L: 
        circ.append(__Ryz__(theta), [q, q+1])
        q += 2
        
    return circ

def __Ryx__(theta): 
    qc = QuantumCircuit(2, name='Ryx')
    qc.s(0)
    qc.h(0)
    qc.z(1)
    qc.h(1)
    qc.s(1)
    qc.cx(0, 1)
    qc.ry(theta, 0)
    qc.rz(theta, 1)
    qc.cx(0, 1)
    qc.h(0)
    qc.sdg(0)
    qc.sdg(1)
    qc.h(1)
    qc.z(1)
    
    gate = qc.to_gate(label=fr'$e^{{-i \theta (YX + XY) / 2}}$')
    return gate

def __W_YX__(circ, L, theta): 
    q = 0 
    while q+1 < L: 
        circ.append(__Ryx__(theta), [q, q+1])
        q += 2

    q = 1
    while q+1 < L: 
        circ.append(__Ryx__(theta), [q, q+1])
        q += 2
        
    return circ

def __cz_layer__(circ, L):
    q = 0 
    while q+1 < L: 
        circ.cz(q, q+1)
        q += 2

    q = 1
    while q+1 < L: 
        circ.cz(q, q+1)
        q += 2
    return circ

def __W_ZYZ__(circ, L, theta):
    circ.compose(__cz_layer__(circ, L))
    for q in range(L-2): 
        circ.ry(theta, q+1)
    circ.compose(__cz_layer__(circ, L))
    return circ

def __H_tilde__(): 
    qc = QuantumCircuit(1, name='Htilde')
    qc.sdg(0)
    qc.h(0)
    qc.s(0)
    gate = qc.to_gate(label=r'$\tilde{H}$')
    return gate

def __W_ZXY__(circ, L, theta): 
    ### O1 
    for n in range(L): 
        if np.mod(n, 4) == 2 or np.mod(n, 4) == 3: 
            circ.append(__H_tilde__(), [n])
    
    circ.compose(__cz_layer__(circ, L))
        
    for n in range(L-2): 
        n += 1
        if np.mod(n, 4) == 2 or np.mod(n, 4) == 3: 
            circ.rx(-theta, n)
        else: 
            circ.rx(theta, n)
            
    circ.compose(__cz_layer__(circ, L))
    
    for n in range(L): 
        if np.mod(n, 4) == 2 or np.mod(n, 4) == 3: 
            circ.append(__H_tilde__(), [n])
    
    ### O2
    for n in range(L): 
        if np.mod(n, 4) == 0 or np.mod(n, 4) == 1: 
            circ.append(__H_tilde__(), [n])
    
    circ.compose(__cz_layer__(circ, L))
        
    for n in range(L-2): 
        n += 1
        if np.mod(n, 4) == 0 or np.mod(n, 4) == 1: 
            circ.rx(-theta, n)
        else: 
            circ.rx(theta, n)
            
    circ.compose(__cz_layer__(circ, L))
    
    for n in range(L): 
        if np.mod(n, 4) == 0 or np.mod(n, 4) == 1: 
            circ.append(__H_tilde__(), [n])
             
    return circ

def __circ_ADAPT_VQE__(L, theta_y1, theta_yz1, theta_y2, theta_zxy, 
                   theta_yz2, theta_yz3,theta_y3, theta_zyz):
    qc_ADAPT_VQE=QuantumCircuit(L)
    #qc_ADAPT_VQE.barrier()
    qc_ADAPT_VQE = __W_Y__(qc_ADAPT_VQE, L, -2*theta_y1)
    #qc_ADAPT_VQE.barrier()
    qc_ADAPT_VQE = __W_YZ__(qc_ADAPT_VQE, L, -2*theta_yz1)
    #qc_ADAPT_VQE.barrier()
    qc_ADAPT_VQE = __W_Y__(qc_ADAPT_VQE, L, -2*theta_y2)
    #qc_ADAPT_VQE.barrier()
    qc_ADAPT_VQE = __W_ZXY__(qc_ADAPT_VQE, L, -2*theta_zxy)
    #qc_ADAPT_VQE.barrier()
    qc_ADAPT_VQE = __W_YZ__(qc_ADAPT_VQE, L, -2*theta_yz2)
    #qc_ADAPT_VQE.barrier()
    qc_ADAPT_VQE = __W_YZ__(qc_ADAPT_VQE, L, -2*theta_yz3)
    #qc_ADAPT_VQE.barrier()
    qc_ADAPT_VQE = __W_Y__(qc_ADAPT_VQE, L, -2*theta_y3)
    #qc_ADAPT_VQE.barrier()
    qc_ADAPT_VQE = __W_ZYZ__(qc_ADAPT_VQE, L, -2*theta_zyz)

    return qc_ADAPT_VQE

@lru_cache()
def get_longest_path(coupling_map: CouplingMap, num_qubits: int = None) -> List:
    graph = nx.Graph(list(coupling_map))
    path = max(nx.all_simple_paths(graph, 14, 114), key=len)

    if num_qubits:
        if num_qubits % 2 != 0:
            raise ValueError("Must be even")
        path = (path[: len(path) - (len(path) % 2)])[:num_qubits]

    return path

#taken from the tables on Pg.42 of paper.
state_thetas_28={"theta_y1": 0.1212,"theta_yz1" : 0.0185, 'theta_y2' : -0.5452, 'theta_zxy' : 0.0397,
'theta_yz2' : 0.0599, 'theta_yz3' : 0.0556, 'theta_y3' : -0.2637, 'theta_zyz' : 0.0566}

state_thetas_100={"theta_y1": 0.1701,"theta_yz1" : 0.0173, 'theta_y2' : -0.5932, 'theta_zxy' : 0.04,
'theta_yz2' : 0.0631, 'theta_yz3' : 0.0554, 'theta_y3' : -0.2647, 'theta_zyz' : 0.0567}


def build_adapt_vqe_state(L:int) -> QuantumCircuit:
    if L <= 28: 
        thetas = state_thetas_28
    else: 
        thetas = state_thetas_100
    adapt_vqe_state = __circ_ADAPT_VQE__(L=L, **thetas)
    return adapt_vqe_state

def prep_wavepacket(kzero: float, L:int, d:int, var_k:float, prep_zero:bool=False)  -> QuantumCircuit:
    # removing x0 to just center at L/2
    norm = 0
    for mprime in range(1, L+1):
        norm += __compute_normalization__(mprime, kzero, var_k, L)
    norm = 1/sqrt(norm)

    results = np.array([__return_params__(n, int(np.floor(L/2)), var_k, kzero, L, norm) for n in range(L)])
    cs, phis = zip(*results)
    cs = np.array(cs)
    phis = np.array(phis)

    cs_truncated, phis_truncated = __truncate_amplitudes__(cs, phis, int(np.floor(L/2)), d, L)
    thetas = __calc_thetas__(cs_truncated, d)
    qc_wstate = QuantumCircuit(d)
    if prep_zero: 
        qc_wstate = __prep_zero_w_state__(qc_wstate, d, phis_truncated, thetas)
    else: 
        qc_wstate = __prep_w_state__(qc_wstate, d, phis_truncated, thetas)
    
    return qc_wstate, cs_truncated
