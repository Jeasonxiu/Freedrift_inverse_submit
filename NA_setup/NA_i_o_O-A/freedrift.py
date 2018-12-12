# functions for playing with free drift calculations.
# dependancies
import numpy as np

# constants
rhoa = 1.25 
rhoi = 917.0 
rhoo = 1026.0 
fc   = 1.46e-4 

# characteristic values
def nansen(Ca,Co):
    return np.sqrt((rhoa*Ca)/(rhoo*Co))

def nansen2_Ua(Ca,Co,Ua):
    return ((rhoa*Ca)/(rhoo*Co))*Ua

def rossby_p(Co,hifc):
    return rhoi*hifc/(rhoo*Co)

def rossby_pp(Co):
    return rhoi/(rhoo*Co)

def nansen_f(Ca,Co,rhoa,rhoo):
    return np.sqrt((rhoa*Ca)/(rhoo*Co))

def nansen2_Ua_f(Ca,Co,Ua,rhoa,rhoo):
    return ((rhoa*Ca)/(rhoo*Co))*Ua

def rossby_p_f(Co,hifc,rhoi,rhoo):
    return rhoi*hifc/(rhoo*Co)

def ice_mag(uwnd,vwnd,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):
    wnd_mag = np.hypot(uwnd,vwnd)
    ocn_mag = np.hypot(uocn,vocn)
    Eq11 = np.tan(ThetaO)
    Eq12 = Ro/wnd_mag
    Eq13 = np.cos(ThetaO)
	
    Eq21 = 2*np.sin(ThetaO)*Ro/wnd_mag
    Eq22 = (Ro**2)/(wnd_mag**2)
    Eq23 = Na2**2
	
    y1 = np.sqrt(Na2)
    yi = y1
    eps = 1.0e-8
    conv = 0
    while conv == 0:
        fi =   yi**4 +   Eq21*(yi**3) +   Eq22*(yi**2) - Eq23
        dfi= 4*yi**3 + 3*Eq21*(yi**2) + 2*Eq22*yi
        yi = yi - fi/dfi
        if (np.abs(fi/dfi) < eps):
            conv = 1
    # x = ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
    x = - ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
    alpha = yi
    ThetaI = x
	
    uice = alpha*(np.cos(-ThetaI)*uwnd - np.sin(-ThetaI)*vwnd) + uocn
    vice = alpha*(np.cos(-ThetaI)*vwnd + np.sin(-ThetaI)*uwnd) + vocn

    return uice, vice

def ocn_mag(uwnd,vwnd,uice,vice,
	    Na2,Ro,ThetaA,ThetaO):
    wnd_mag = np.hypot(uwnd,vwnd)
    ice_mag = np.hypot(uice,vice)
    Eq11 = np.tan(ThetaO)
    Eq12 = Ro/wnd_mag
    Eq13 = np.cos(ThetaO)
	
    Eq21 = 2*np.sin(ThetaO)*Ro/wnd_mag
    Eq22 = (Ro**2)/(wnd_mag**2)
    Eq23 = Na2**2
	
    y1 = np.sqrt(Na2)
    yi = y1
    eps = 1.0e-9
    conv = 0
    while conv == 0:
        fi =   yi**4 +   Eq21*(yi**3) +   Eq22*(yi**2) - Eq23
        dfi= 4*yi**3 + 3*Eq21*(yi**2) + 2*Eq22*yi
        yi = yi - fi/dfi
        if (np.abs(fi/dfi) < eps):
            conv = 1
    # x = ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
    x = - ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
    alpha = yi
    ThetaI = x
	
    uocn = uice - alpha*(np.cos(-ThetaI)*uwnd - np.sin(-ThetaI)*vwnd) 
    vocn = vice - alpha*(np.cos(-ThetaI)*vwnd + np.sin(-ThetaI)*uwnd) 

    return uocn, vocn

def wnd_mag(uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):
    ice_mag = np.hypot(uice,vice)
    ocn_mag = np.hypot(uocn,vocn)
    Eq11 = np.tan(ThetaO)
    Eq13 = np.cos(ThetaO)
    
    Eq23 = Na2**2
    y1 = np.sqrt(Na2)
    yi = y1
    eps = 1.0E-9
    conv = 0
    it = 0
    wind_mag_calc = 10.0
    while (conv == 0):
    # constants that change with iteration
        Eq12 = Ro/wind_mag_calc
    
        Eq21 = 2*np.sin(ThetaO)*Ro/wind_mag_calc
        Eq22 = (Ro**2)/(wind_mag_calc**2)
    
        fi =   yi**4 +   Eq21*(yi**3) +   Eq22*(yi**2) - Eq23
        dfi= 4*yi**3 + 3*Eq21*(yi**2) + 2*Eq22*(yi)
        yi = yi - fi/dfi
    
    #     c need to put the solution through to the wind mag to properly converge
        x = - ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
        uwind_calc = (np.cos(x)*(uice - uocn) - np.sin(x)*(vice - vocn))/yi
        vwind_calc = (np.cos(x)*(vice - vocn) + np.sin(x)*(uice - uocn))/yi
        wind_mag_calc = np.hypot(uwind_calc, vwind_calc)
    
        it = it + 1
        if(np.abs(fi/dfi) < eps):
            conv = 1
        if(it < 10):
            eps = 1.0E-9
        elif(it < 50):
            eps = 1.0E-8
        elif(it < 200):
            eps = 1.0E-7

    return uwind_calc, vwind_calc

def misfit_fb(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):
    # now try back in the mtm balance
    wnd_mag = np.hypot(uwnd,vwnd)

    u_ice_ocn = uocn - uice
    v_ice_ocn = vocn - vice
    ice_ocn_mag = np.hypot(u_ice_ocn,v_ice_ocn)
    
    atmo_stressx = Na2*wnd_mag*( uwnd*np.cos(ThetaA) - vwnd*np.sin(ThetaA) )
    atmo_stressy = Na2*wnd_mag*( vwnd*np.cos(ThetaA) + uwnd*np.sin(ThetaA) )
    
    ocn_stressx  = ice_ocn_mag*( u_ice_ocn*np.cos(ThetaO) - v_ice_ocn*np.sin(ThetaO) )
    ocn_stressy  = ice_ocn_mag*( v_ice_ocn*np.cos(ThetaO) + u_ice_ocn*np.sin(ThetaO) )
    
    cor_tiltx = -v_ice_ocn*Ro
    cor_tilty =  u_ice_ocn*Ro
    
    tot_mtm_x = atmo_stressx + ocn_stressx + cor_tiltx
    tot_mtm_y = atmo_stressy + ocn_stressy + cor_tilty

    return np.hypot(tot_mtm_x, tot_mtm_y)

def misfit_fb_scaled(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):

    u_ice_ocn = uocn - uice
    v_ice_ocn = vocn - vice
    ice_ocn_mag = np.hypot(u_ice_ocn,v_ice_ocn)
    
    wnd_mag = np.hypot(uwnd,vwnd)

    misfit = misfit_fb(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO)
    
    return misfit / (ice_ocn_mag**2 * wnd_mag**2)

def misfit_a(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):
    uwnd_calc, vwnd_calc = wnd_mag(uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO)

    norm_u = (uwnd - uwnd_calc)**2
    norm_v = (vwnd - vwnd_calc)**2

    return norm_u + norm_v

def misfit_a_scaled(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):

    u_ice_ocn = uocn - uice
    v_ice_ocn = vocn - vice
    ice_ocn_mag = np.hypot(u_ice_ocn,v_ice_ocn)
    
    wnd_mag = np.hypot(uwnd,vwnd)

    misfit = misfit_a(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO)
    
    return misfit * ice_ocn_mag**2 / wnd_mag**2

def misfit_i(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):
    uice_calc, vice_calc = ice_mag(uwnd,vwnd,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO)

    norm_u = (uice - uice_calc)**2
    norm_v = (vice - vice_calc)**2

    return norm_u + norm_v

def misfit_i_scaled(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):

    wnd_mag = np.hypot(uwnd,vwnd)

    misfit = misfit_i(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO)
    
    return misfit / wnd_mag

def misfit_o(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):
    uocn_calc, vocn_calc = ocn_mag(uwnd,vwnd,uice,vice,
	    Na2,Ro,ThetaA,ThetaO)

    norm_u = (uocn - uocn_calc)**2
    norm_v = (vocn - vocn_calc)**2

    return norm_u + norm_v

def misfit_o_scaled(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):

    wnd_mag = np.hypot(uwnd,vwnd)

    misfit = misfit_o(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO)
    
    return misfit / wnd_mag

