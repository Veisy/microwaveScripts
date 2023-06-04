from numpy import *


# Convert from phasor to rectangular
def phasor_to_rectangular(magnitude, angle):
    '''
    Converts a phasor to rectangular form

    Parameters
    ----------
    :param magnitude: Magnitude of phasor (float)
    :param angle: Angle of phasor (float)
    :return: Rectangular form of phasor (complex)

    Example
    -------
    magnitude = 1
    angle = 45
    rectangular = phasor_to_rectangular(magnitude, angle)
    rectangular
    (0.7071067811865476+0.7071067811865475j)
    '''
    rectangular = magnitude * exp(1j * angle * pi / 180)
    return rectangular


def convert_gamma_to_z(gamma, z0=50):
    '''
    Converts reflection coefficient Γ to impedance Z

    Parameters
    ----------
    :param z0:
    :param gamma: Γ (complex)
    :return: Z (complex)

    Example
    -------
    Γ = 0.3-0.5j
    Z = convert_gamma_to_z(Γ)
    Z
    (0.6000000000000001-0.5j)
    '''
    z = z0 * (1 + gamma) / (1 - gamma)
    return z


def transducer_gain(S, Gams, Gaml):
    '''
    Calculates the transducer gain for a given S-Matrix
    and Γ_s and Γ_l

    Parameters
    ----------
    :param S: (2x2 Matrix complex) or array of matrices
    :param Gams: Γ_s (complex)
    :param Gaml: Γ_l (complex)
    :return: Transducer Gain Gt in dB (float or array float)

    Example
    -------
     S11=S[0,0];S12=S[0,1];S21=S[1,0];S22=S[1,1]
     Gam_s, Gam_l = 0.3, 0.1-0.5j
     GtdB = transducerGain(S, Gam_s, Gam_l)
    GtdB
    16.51983047097562
    '''
    S11 = S[0, 0]
    S12 = S[0, 1]
    S21 = S[1, 0]
    S22 = S[1, 1]
    Gamin = S11 + (S12 * S21 * Gaml) / (1 - S22 * Gaml)
    Gamout = S22 + (S12 * S21 * Gams) / (1 - S11 * Gams)
    Gs = (1 - abs(Gams) ** 2) / abs(1 - Gamin * Gams) ** 2
    G0 = abs(S21) ** 2
    Gl = (1 - abs(Gaml) ** 2) / abs(1 - S22 * Gaml) ** 2
    Gt = Gs * G0 * Gl

    print("Gs: ", Gs)
    print("G0: ", G0)
    print("Gl: ", Gl)
    print("Gt: ", Gt)

    GsdB = 10 * log10(Gs)
    G0dB = 10 * log10(G0)
    GldB = 10 * log10(Gl)
    GtdB = GsdB + G0dB + GldB

    print("GsdB: ", GsdB)
    print("G0dB: ", G0dB)
    print("GldB: ", GldB)
    print("GtdB: ", GtdB)
    return GtdB


S11 = phasor_to_rectangular(0.5105, 105.7)
S12 = phasor_to_rectangular(0.1314, 19.9)
S21 = phasor_to_rectangular(3.454, 19.3)
S22 = phasor_to_rectangular(0.0930, 107.5)
print("S Parameters S11:", S11, ", S12:", S12, ", S21:", S21, ", S22:", S22)
print("S22 conjugate: ", conj(S22))

delta = (S11 * S22) - (S12 * S21)
C1 = S11 - delta * conj(S22)
C2 = S22 - delta * conj(S11)
B1 = 1 + abs(S11) ** 2 - abs(S22) ** 2 - abs(delta) ** 2
B2 = 1 + abs(S22) ** 2 - abs(S11) ** 2 - abs(delta) ** 2

Gam_s_pozitif = (B1 + sqrt(B1 ** 2 - 4 * abs(C1) ** 2)) / (2 * C1)
Gam_s_pozitif_impedance = convert_gamma_to_z(Gam_s_pozitif)

Gam_s_negatif = (B1 - sqrt(B1 ** 2 - 4 * abs(C1) ** 2)) / (2 * C1)
Gam_s_negatif_impedance = convert_gamma_to_z(Gam_s_negatif)

Gam_l_pozitif = (B2 + sqrt(B2 ** 2 - 4 * abs(C2) ** 2)) / (2 * C2)
Gam_l_pozitif_impedance = convert_gamma_to_z(Gam_l_pozitif)

Gam_l_negatif = (B2 - sqrt(B2 ** 2 - 4 * abs(C2) ** 2)) / (2 * C2)
Gam_l_negative_impedance = convert_gamma_to_z(Gam_l_negatif)

k = (1 - abs(S11) ** 2 - abs(S22) ** 2 + abs(delta) ** 2) / (2 * abs(S12 * S21))
print("K: ", k)
print("Gam_s positive: ", Gam_s_pozitif, "| Impedance: ", Gam_s_pozitif_impedance)
print("Gam_s negative: ", Gam_s_negatif, "| Impedance: ", Gam_s_negatif_impedance)
print("Gam_l positive: ", Gam_l_pozitif, "| Impedance: ", Gam_l_pozitif_impedance)
print("Gam_l negative: ", Gam_l_negatif, "| Impedance: ", Gam_l_negative_impedance)
print("-------------")

S = matrix([[S11, S12], [S21, S22]])

transducer_gain(S, Gam_s_negatif, Gam_l_negatif)
