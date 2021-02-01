from __future__ import division
import numpy as np
import math
import lal

G = lal.G_SI
c = lal.C_SI

"""
    BBH final-state fitting formulas from https://arxiv.org/abs/1611.00332
    (c) 2016-2017 Xisco Jimenez-Forteza, David Keitel, Sascha Husa, Mark Hannam, Sebastian Khan, Michael Puerrer
    also included in LALInference under GPL at
    https://versions.ligo.org/cgit/lalsuite/tree/lalinference/python/lalinference/imrtgr/nrutils.py
"""

def fits_setup(m1, m2, chi1, chi2):
    """
    common setup function for UIB final-state and luminosity fit functions
    """

    # Vectorize the function if arrays are provided as input
    m1   = np.vectorize(float)(np.array(m1))
    m2   = np.vectorize(float)(np.array(m2))
    chi1 = np.vectorize(float)(np.array(chi1))
    chi2 = np.vectorize(float)(np.array(chi2))

    if np.any(m1<0):
      raise ValueError("m1 must not be negative")
    if np.any(m2<0):
      raise ValueError("m2 must not be negative")

    if np.any(abs(chi1)>1):
      raise ValueError("chi1 has to be in [-1, 1]")
    if np.any(abs(chi2)>1):
      raise ValueError("chi2 has to be in [-1, 1]")

    # binary masses
    m    = m1+m2
    if np.any(m<=0):
      raise ValueError("m1+m2 must be positive")
    msq  = m*m
    m1sq = m1*m1
    m2sq = m2*m2

    # symmetric mass ratio
    eta  = m1*m2/msq
    if np.any(eta>0.25):
      print("Truncating eta from above to 0.25. This should only be necessary in some rounding corner cases, but better check your m1 and m2 inputs...")
      eta = np.minimum(eta,0.25)
    if np.any(eta<0.0):
      print("Truncating negative eta to 0.0. This should only be necessary in some rounding corner cases, but better check your m1 and m2 inputs...")
      eta = np.maximum(eta,0.0)
    eta2 = eta*eta
    eta3 = eta2*eta
    eta4 = eta2*eta2

    # spin variables (in m = 1 units)
    S1    = chi1*m1sq/msq # spin angular momentum 1
    S2    = chi2*m2sq/msq # spin angular momentum 2
    Stot  = S1+S2         # total spin
    Shat  = (chi1*m1sq+chi2*m2sq)/(m1sq+m2sq) # effective spin, = msq*Stot/(m1sq+m2sq)
    Shat2 = Shat*Shat
    Shat3 = Shat2*Shat
    Shat4 = Shat2*Shat2

    # spin difference, assuming m1>m2
    chidiff  = chi1 - chi2
    if np.any(m2>m1): # fit assumes m1>m2
      chidiff = np.sign(m1-m2)*chidiff
    chidiff2 = chidiff*chidiff

    # typical squareroots and functions of eta
    sqrt2 = 2.**0.5
    sqrt3 = 3.**0.5
    sqrt1m4eta = (1. - 4.*eta)**0.5

    return m, eta, eta2, eta3, eta4, Stot, Shat, Shat2, Shat3, Shat4, chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta

def final_mass(m1, m2, chi1, chi2, version="v2"):
    """
    Calculate the final mass with the aligned-spin NR fit
    by Xisco Jimenez Forteza, David Keitel, Sascha Husa et al.
    [LIGO-P1600270] [https://arxiv.org/abs/1611.00332]
    versions v1 and v2 use the same ansatz,
    with v2 calibrated to additional SXS and RIT data

    m1, m2: component masses
    chi1, chi2: dimensionless spins of two BHs
    Note: Here it is assumed that m1>m2.
    """

    m, eta, eta2, eta3, eta4, Stot, Shat, Shat2, Shat3, Shat4, chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta = fits_setup(m1, m2, chi1, chi2)

    if version == "v1":
        # rational-function Pade coefficients (exact) from Eq. (22) of 1611.00332v1
        b10 = 0.487
        b20 = 0.295
        b30 = 0.17
        b50 = -0.0717
        # fit coefficients from Tables VII-X of 1611.00332v1
        # values at increased numerical precision copied from
        # https://git.ligo.org/uib-papers/finalstate2016/blob/master/LALInference/EradUIB2016_pyform_coeffs.txt
        # git commit 7b47e0f35a8f960b99b24caf3ffea2ddefdc4e29
        a2 = 0.5635376058169299
        a3 = -0.8661680065959881
        a4 = 3.181941595301782
        b1 = -0.15800074104558132
        b2 = -0.15815904609933157
        b3 = -0.14299315232521553
        b5 = 8.908772171776285
        f20 = 3.8071100104582234
        f30 = 25.99956516423936
        f50 = 1.552929335555098
        f10 = 1.7004558922558886
        f21 = 0.
        d10 = -0.12282040108157262
        d11 = -3.499874245551208
        d20 = 0.014200035799803777
        d30 = -0.01873720734635449
        d31 = -5.1830734185518725
        f11 = 14.39323998088354
        f31 = -232.25752840151296
        f51 = -0.8427987782523847

    elif version == "v2":
        # rational-function Pade coefficients (exact) from Eq. (22) of 1611.00332v2
        b10 = 0.346
        b20 = 0.211
        b30 = 0.128
        b50 = -0.212
        # fit coefficients from Tables VII-X of 1611.00332v2
        # values at increased numerical precision copied from
        # https://git.ligo.org/uib-papers/finalstate2016/blob/master/LALInference/EradUIB2016v2_pyform_coeffs.txt
        # git commit f490774d3593adff5bb09ae26b7efc6deab76a42
        a2 = 0.5609904135313374
        a3 = -0.84667563764404
        a4 = 3.145145224278187
        b1 = -0.2091189048177395
        b2 = -0.19709136361080587
        b3 = -0.1588185739358418
        b5 = 2.9852925538232014
        f20 = 4.271313308472851
        f30 = 31.08987570280556
        f50 = 1.5673498395263061
        f10 = 1.8083565298668276
        f21 = 0.
        d10 = -0.09803730445895877
        d11 = -3.2283713377939134
        d20 = 0.01118530335431078
        d30 = -0.01978238971523653
        d31 = -4.91667749015812
        f11 = 15.738082204419655
        f31 = -243.6299258830685
        f51 = -0.5808669012986468

    else:
        raise ValueError('Unknown version -- should be either "v1" or "v2".')

    # Calculate the radiated-energy fit from Eq. (27) of 1611.00332
    Erad = (((1. + -2.0/3.0*sqrt2)*eta + a2*eta2 + a3*eta3 + a4*eta4)*(1. + b10*b1*Shat*(f10 + f11*eta + (16. - 16.*f10 - 4.*f11)*eta2) + b20*b2*Shat2*(f20 + f21*eta + (16. - 16.*f20 - 4.*f21)*eta2) + b30*b3*Shat3*(f30 + f31*eta + (16. - 16.*f30 - 4.*f31)*eta2)))/(1. + b50*b5*Shat*(f50 + f51*eta + (16. - 16.*f50 - 4.*f51)*eta2)) + d10*sqrt1m4eta*eta2*(1. + d11*eta)*chidiff + d30*Shat*sqrt1m4eta*eta*(1. + d31*eta)*chidiff + d20*eta3*chidiff2

    # Convert to actual final mass
    Mf = m*(1.-Erad)

    return Mf

def final_spin(m1, m2, chi1, chi2, version="v2"):
    """
    Calculate the final spin with the aligned-spin NR fit
    by Xisco Jimenez Forteza, David Keitel, Sascha Husa et al.
    [LIGO-P1600270] [https://arxiv.org/abs/1611.00332]
    versions v1 and v2 use the same ansatz,
    with v2 calibrated to additional SXS and RIT data

    m1, m2: component masses
    chi1, chi2: dimensionless spins of two BHs
    Note: Here it is assumed that m1>m2.
    """

    m, eta, eta2, eta3, eta4, Stot, Shat, Shat2, Shat3, Shat4, chidiff, chidiff2, sqrt2, sqrt3, sqrt1m4eta = fits_setup(m1, m2, chi1, chi2)

    if version == "v1":
        # rational-function Pade coefficients (exact) from Eqs. (7) and (8) of 1611.00332v1
        a20 = 5.28
        a30 = 1.27
        a50 = 2.89
        b10 = -0.194
        b20 = 0.075
        b30 = 0.00782
        b50 = -0.527
        # fit coefficients from Tables I-IV of 1611.00332v1
        # evalues at increased numerical precision copied from
        # https://git.ligo.org/uib-papers/finalstate2016/blob/master/LALInference/FinalSpinUIB2016_pyform_coeffs.txt
        # git commit 7b47e0f35a8f960b99b24caf3ffea2ddefdc4e29
        a2 = 3.772362507208651
        a3 = -9.627812453422376
        a5 = 2.487406038123681
        b1 = 1.0005294518146604
        b2 = 0.8823439288807416
        b3 = 0.7612809461506448
        b5 = 0.9139185906568779
        f21 = 8.887933111404559
        f31 = 23.927104476660883
        f50 = 1.8981657997557002
        f11 = 4.411041530972546
        f52 = 0.
        d10 = 0.2762804043166152
        d11 = 11.56198469592321
        d20 = -0.05975750218477118
        d30 = 2.7296903488918436
        d31 = -3.388285154747212
        f12 = 0.3642180211450878
        f22 = -40.35359764942015
        f32 = -178.7813942566548
        f51 = -5.556957394513334

    elif version == "v2":
        # rational-function Pade coefficients (exact) from Eqs. (7) and (8) of 1611.00332v2
        a20 = 5.24
        a30 = 1.3
        a50 = 2.88
        b10 = -0.194
        b20 = 0.0851
        b30 = 0.00954
        b50 = -0.579
        # fit coefficients from Tables I-IV of 1611.00332v2
        # values at increased numerical precision copied from
        # https://git.ligo.org/uib-papers/finalstate2016/blob/master/LALInference/FinalSpinUIB2016v2_pyform_coeffs.txt
        # git commit f490774d3593adff5bb09ae26b7efc6deab76a42
        a2 = 3.8326341618708577
        a3 = -9.487364155598392
        a5 = 2.5134875145648374
        b1 = 1.0009563702914628
        b2 = 0.7877509372255369
        b3 = 0.6540138407185817
        b5 = 0.8396665722805308
        f21 = 8.77367320110712
        f31 = 22.830033250479833
        f50 = 1.8804718791591157
        f11 = 4.409160174224525
        f52 = 0.
        d10 = 0.3223660562764661
        d11 = 9.332575956437443
        d20 = -0.059808322561702126
        d30 = 2.3170397514509933
        d31 = -3.2624649875884852
        f12 = 0.5118334706832706
        f22 = -32.060648277652994
        f32 = -153.83722669033995
        f51 = -4.770246856212403

    else:
        raise ValueError('Unknown version -- should be either "v1" or "v2".')

    # Calculate the fit for the Lorb' quantity from Eq. (16) of 1611.00332
    Lorb = (2.*sqrt3*eta + a20*a2*eta2 + a30*a3*eta3)/(1. + a50*a5*eta) + (b10*b1*Shat*(f11*eta + f12*eta2 + (64. - 16.*f11 - 4.*f12)*eta3) + b20*b2*Shat2*(f21*eta + f22*eta2 + (64. - 16.*f21 - 4.*f22)*eta3) + b30*b3*Shat3*(f31*eta + f32*eta2 + (64. - 16.*f31 - 4.*f32)*eta3))/(1. + b50*b5*Shat*(f50 + f51*eta + f52*eta2 + (64. - 64.*f50 - 16.*f51 - 4.*f52)*eta3)) + d10*sqrt1m4eta*eta2*(1. + d11*eta)*chidiff + d30*Shat*sqrt1m4eta*eta3*(1. + d31*eta)*chidiff + d20*eta3*chidiff2

    # Convert to actual final spin, Stot and Lorb were adimensional, so it is correct, simply M_in_tot=1
    chif = Lorb + Stot

    return chif


def frequency22_merger(m1, m2, a1, a2): #Merger frequency defined with respect to the 22 mode of the inspiral (2=quadrupole and 2=...) (Ref: Nagar)
    
    q = m1/m2                   #Husa conventions, m1>m2 [https://arxiv.org/abs/1611.00332]
    eta = q/(1+q)**2
    M_tot = m1+m2
    chi_1 = 0.5*(1.0+np.sqrt(1.0-4.0*eta))
    chi_2 = 1.0-chi_1
    chi_eff = chi_1*a1 + chi_2*a2
    
    A = -0.28562363*eta + 0.090355762
    B = -0.18527394*eta + 0.12596953
    C =  0.40527397*eta + 0.25864318
    
    res = (A*chi_eff**2 + B*chi_eff + C)*((2*math.pi*M_tot)*lal.G_SI*lal.C_SI**(-3))**(-1)
    return res



def omega_22(M,a):
    
    return (1.5251-1.1568*((1-a)**(0.1292)))/(M*G*(c**(-3)))

def Q_22(a):
    
    return 0.7+1.4187*(1.0-a)**(-0.4990)

def tau_22(M,a):
    
    return (2*Q_22(a)/omega_22(M,a))

def chi_eff(q, a1, a2):                       #eta=0.25 -> chi_1=chi_2=0.5 -> chi_eff=aritmetic mean of (a1,a2)

    eta = q/(1+q)**2
    chi_1 = 0.5*(1.0+np.sqrt(1.0-4.0*eta))
    chi_2 = 1.0-chi_1
    
    return chi_1*a1 + chi_2*a2

def chi_p(m1,m2,s1x,s1y,s2x,s2y):
 
    # Magnitude of the spin projections in the orbital plane
    S1_perp = m1*m1*np.sqrt(s1x*s1x + s1y*s1y)
    S2_perp = m2*m2*np.sqrt(s2x*s2x + s2y*s2y)
    A1 = 2. + (3.*m2) / (2.*m1)
    A2 = 2. + (3.*m1) / (2.*m2)
    ASp1 = A1*S1_perp
    ASp2 = A2*S2_perp
    if ASp2 > ASp1:
        return ASp2/(A1*m1**2)
    else:
        return ASp1/(A1*m1**2)

def McQ2Masses(mc, q):
    """
    Simple utility to convert between mc an d q to component masses
    """
    factor = mc * np.power(1. + q, 1.0/5.0);
    m1 = factor * np.power(q, -3.0/5.0);
    m2 = factor * np.power(q, +2.0/5.0);
    return m1, m2

def Masses2McQ(m1, m2):
    """
    Simple utility to convert between component masses to mc q
    """
    q = m2/m1
    eta = m1*m2/(m1+m2)
    mc = (m1*m2)**(3./5.)/(m1+m2)**(1./5.)
    return mc, q

def PolarToCartesian(a, th, ph):
    """
    Simple utility function to convert between polar and cartesian representations
    for the spin quantities
    """
    return a*np.cos(th)*np.cos(ph), a*np.cos(th)*np.sin(ph), a*np.sin(th)

def RedshiftedMass(m,z):
    return m*(1+z)
