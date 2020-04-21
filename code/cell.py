def reaction(z, t, k1, k_1, k2, k_2, E_0):
    '''Reaction model'''

    #Calculate constants
    Vs = k2*E_0
    Vp = k_1*E_0
    Ks = (k_1 + k2)/k1
    Kp = (k_1 + k2)/k_2

    #Variables
    S, P = z 

    #Evaluate model
    dPdt = (Vs/Ks*S - Vp/Kp*P)/(1 + S/Ks + P/Kp) 
    dsdt = -1*dpdt 
    dzdt = [dSdt, dPdt]

    return dzdt
