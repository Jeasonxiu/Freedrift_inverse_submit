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
    its = 0
    while conv == 0:
        fi =   yi**4 +   Eq21*(yi**3) +   Eq22*(yi**2) - Eq23
        dfi= 4*yi**3 + 3*Eq21*(yi**2) + 2*Eq22*yi
        yi = yi - fi/dfi
        its = its + 1
        if (np.abs(fi/dfi) < eps):
            conv = 1
        if(its < 10):
            eps = 1.0E-9
        elif (its < 20):
            eps = 1.0e-8
        elif (its < 50):
            eps = 1.0e-7
        elif (its < 200):
            eps = 1.0e-6
        elif (its < 1000):
            eps = 1.0e-2
    # x = ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
    # x = - ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
    x = - ThetaA + np.arctan2((yi*Eq13),Eq11*(yi*Eq13) + Eq12)
    alpha = yi
    ThetaI = x
	
    uice = alpha*(np.cos(-ThetaI)*uwnd - np.sin(-ThetaI)*vwnd) + uocn
    vice = alpha*(np.cos(-ThetaI)*vwnd + np.sin(-ThetaI)*uwnd) + vocn

    return uice, vice

def o_i_mag(uwnd,vwnd,
	    Na2,Ro,ThetaA,ThetaO):
    wnd_mag = np.hypot(uwnd,vwnd)
    Eq11 = np.sin(ThetaO)
    Eq12 = Ro/wnd_mag
    Eq13 = np.cos(ThetaO)
	
    Eq21 = 2*np.sin(ThetaO)*Ro/wnd_mag
    Eq22 = (Ro**2)/(wnd_mag**2)
    Eq23 = Na2**2
	
    # y1 = np.sqrt(Na2)
    y0 = 0.04
    y1 = 0.035
    yi1 = y0
    yi  = y1
    # prep for new method
    fi =   yi1**4 +   Eq21*(yi1**3) +   Eq22*(yi1**2) - Eq23
    mi  =0.0
    eps = 1.0E-9
    conv = 0
    its = 0
    while (conv == 0):
        fi1 = fi
        fi  =   yi**4 +   Eq21*(yi**3) +   Eq22*(yi**2) - Eq23
        dfi = 4*yi**3 + 3*Eq21*(yi**2) + 2*Eq22*(yi) #+1E-10
        mi  = (fi - fi1)/(yi - yi1 + 1e-10)
        mi  = mi / dfi
        chng= (yi1 - yi)/(1.0 - mi*fi/fi1)
        yin = yi1 - chng
        yi1 = yi
        yi  = yin
        its = its + 1
        if (np.abs(chng) < eps):
            conv = 1
    	    # x = ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
            # x = - ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
            #x = - ThetaA + np.arctan2((yi*Eq13),Eq11*(yi*Eq13) + Eq12)
            x = - ThetaA + np.arctan2(Eq11*yi + Eq12,yi*Eq13)
            #x =   ThetaA - np.arctan2(Eq11*yi + Eq12,yi*Eq13)
            alpha = yi
            ThetaI = x

            uo_i = -( alpha*(np.cos(-ThetaI)*uwnd - np.sin(-ThetaI)*vwnd) )
            vo_i = -( alpha*(np.cos(-ThetaI)*vwnd + np.sin(-ThetaI)*uwnd) )
	
        if(its < 10):
            eps = 1.0E-9
        elif (its < 20):
            eps = 1.0e-8
        elif (its < 50):
            eps = 1.0e-7
        elif (its < 200):
            eps = 1.0e-6
        elif (its < 1000):
            eps = 1.0e-5
        elif (its < 2000):
            conv = 1
            uo_i = 0.1
            vo_i = 0.1
            print(ThetaO)
    return uo_i, vo_i

def ice_alpha(hifc_ua,
	    Na2,Ro,ThetaA,ThetaO):
    Eq11 = np.tan(ThetaO)
    Eq12 = Ro*hifc_ua
    Eq13 = np.cos(ThetaO)
	
    Eq21 = 2*np.sin(ThetaO)*Ro*hifc_ua
    Eq22 = (Ro**2)*(hifc_ua**2)
    Eq23 = Na2**2
	
    #y1 = np.sqrt(Na2)
    y1 = 0.02
    yi = y1
    eps = 1.0e-8
    conv = 0
    its = 0
    while conv == 0:
        fi =   yi**4 +   Eq21*(yi**3) +   Eq22*(yi**2) - Eq23
        dfi= 4*yi**3 + 3*Eq21*(yi**2) + 2*Eq22*yi
        yi = yi - fi/dfi
        its = its + 1
        if (np.abs(fi/dfi) < eps):
            conv = 1
        if(its < 10):
            eps = 1.0E-9
        elif (its < 20):
            eps = 1.0e-8
        elif (its < 50):
            eps = 1.0e-7
        elif (its < 200):
            eps = 1.0e-6
        elif (its < 1000):
            eps = 1.0e-2
    # x = ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
    # x =  ThetaA - np.arctan(Eq11 + Eq12/(yi*Eq13))
    x = - ThetaA + np.arctan2(Eq11*(yi*Eq13) + Eq12,yi*Eq13)
    alpha = yi
    ThetaI = x
    
    return alpha, ThetaI

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
    # x = - ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
    x = - ThetaA + np.arctan2((yi*Eq13),Eq11*(yi*Eq13) + Eq12)
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
    its = 0
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
        # x = - ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
        x = - ThetaA + np.arctan2((yi*Eq13),Eq11*(yi*Eq13) + Eq12)
        uwind_calc = (np.cos(x)*(uice - uocn) - np.sin(x)*(vice - vocn))/yi
        vwind_calc = (np.cos(x)*(vice - vocn) + np.sin(x)*(uice - uocn))/yi
        wind_mag_calc = np.hypot(uwind_calc, vwind_calc)
    
        its = its + 1
        if(np.abs(fi/dfi) < eps):
            conv = 1
        if(its < 10):
            eps = 1.0E-9
        elif (its < 20):
            eps = 1.0e-8
        elif (its < 50):
            eps = 1.0e-7
        elif (its < 200):
            eps = 1.0e-6
        elif (its < 1000):
            eps = 1.0e-2
        elif (its < 2000):
            eps = 1.2e-0

    return uwind_calc, vwind_calc

def wnd_o_i_mag(uo_i,vo_i,
	    Na2,Ro,ThetaA,ThetaO):
    o_i_mag = np.hypot(uo_i,vo_i)
    
    #Eq11 = np.tan(ThetaO)
    Eq11 = np.sin(ThetaO)
    Eq13 = np.cos(ThetaO)
    
    Eq23 = Na2**2

# prep need for new two point method
    wind_mag_calc = 1.0
    # y1 = np.sqrt(Na2)
    y0 = 0.04
    y1 = 0.035
    yi1 = y0
    yi  = y1
    # prep for new method
    Eq12 = Ro/wind_mag_calc
       
    Eq21 = 2*np.sin(ThetaO)*Ro/wind_mag_calc
    Eq22 = (Ro**2)/(wind_mag_calc**2)
    fi =   yi1**4 +   Eq21*(yi1**3) +   Eq22*(yi1**2) - Eq23
    mi  =0.0
    eps = 1.0E-9
    conv = 0
    its = 0
    while (conv == 0):
        # constants that change with iteration
        Eq12 = Ro/wind_mag_calc
       
        Eq21 = 2*np.sin(ThetaO)*Ro/wind_mag_calc
        Eq22 = (Ro**2)/(wind_mag_calc**2)
       
        fi1 = fi
        fi  =   yi**4 +   Eq21*(yi**3) +   Eq22*(yi**2) - Eq23
        dfi = 4*yi**3 + 3*Eq21*(yi**2) + 2*Eq22*(yi) #+1E-10
        mi  = (fi - fi1)/(yi - yi1 + 1e-10)
        mi  = mi / dfi
        chng= (yi1 - yi)/(1.0 - mi*fi/fi1)
        yin = yi1 - chng
        yi1 = yi
        yi  = yin

    #     c need to put the solution through to the wind mag to properly converge
        # x = - ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))
        #tan_top = (yi+ Eq12*np.sin(ThetaO))*np.sin(ThetaA-ThetaO) - Eq12*np.cos(ThetaO)*np.cos(ThetaA-ThetaO)
        #tan_bot = (yi+ Eq12*np.sin(ThetaO))*np.cos(ThetaA-ThetaO) + Eq12*np.cos(ThetaO)*np.sin(ThetaA-ThetaO)
        #x = np.arctan2(tan_top,tan_bot)
        #x = - ThetaA + np.arctan2((yi*Eq13),Eq11*(yi*Eq13) + Eq12)
        x = - ThetaA + np.arctan2(Eq11*(yi*Eq13) + Eq12,(yi*Eq13))
        #x =   ThetaA - np.arctan2(Eq11*yi + Eq12,yi*Eq13)
        uwind_calc = -(np.cos(x)*(uo_i) - np.sin(x)*(vo_i))/yi
        vwind_calc = -(np.cos(x)*(vo_i) + np.sin(x)*(uo_i))/yi
        wind_mag_calc = np.hypot(uwind_calc, vwind_calc)
    
        its = its + 1
        if(np.abs(chng) < eps):
            conv = 1
            # print('its = ',its)
        if(its < 10):
            eps = 1.0E-9
        elif (its < 20):
            eps = 1.0e-8
        elif (its < 50):
            eps = 1.0e-7
        elif (its < 200):
            eps = 1.0e-6
        elif (its < 1000):
            eps = 1.0e-5
        elif (its < 2000):
            conv = 1
            uwind_calc = 0.1
            vwind_calc = 0.1
            wind_mag_calc = 1.0


    return uwind_calc, vwind_calc

def misfit_fb(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):
    # now try back in the mtm balance
    wnd_mag = np.hypot(uwnd,vwnd)

    u_ice_ocn = uocn - uice
    v_ice_ocn = vocn - vice
    ice_ocn_mag = np.hypot(u_ice_ocn,v_ice_ocn)

    ThetaAO = ThetaA - ThetaO
    
    atmo_stressx = Na2*wnd_mag*( uwnd*np.cos(ThetaAO) - vwnd*np.sin(ThetaAO) )
    atmo_stressy = Na2*wnd_mag*( vwnd*np.cos(ThetaAO) + uwnd*np.sin(ThetaAO) )
    
    ocn_stressx  = ice_ocn_mag*( u_ice_ocn )
    ocn_stressy  = ice_ocn_mag*( v_ice_ocn )
    
    cor_tiltx = -v_ice_ocn*Ro
    cor_tilty =  u_ice_ocn*Ro
    cor_tiltx = -(v_ice_ocn*np.cos(-ThetaO) - u_ice_ocn*np.sin(-ThetaO))*Ro
    cor_tilty =  (u_ice_ocn*np.cos(-ThetaO) + v_ice_ocn*np.sin(-ThetaO))*Ro
    
    tot_mtm_x = atmo_stressx + ocn_stressx + cor_tiltx
    tot_mtm_y = atmo_stressy + ocn_stressy + cor_tilty

    # return np.hypot(tot_mtm_x, tot_mtm_y)**2
    return tot_mtm_x**2 +  tot_mtm_y**2

def misfit_fb2(uwnd,vwnd,uo_i,vo_i,
	    Capa,Copo,rhohifc,ThetaA,ThetaO):
    # now try back in the mtm balance
    wnd_mag = np.hypot(uwnd,vwnd)

    ice_ocn_mag = np.hypot(uo_i,vo_i)

    atmo_stressx = Capa*wnd_mag*( uwnd*np.cos(ThetaA) - vwnd*np.sin(ThetaA) )
    atmo_stressy = Capa*wnd_mag*( vwnd*np.cos(ThetaA) + uwnd*np.sin(ThetaA) )
    
    ocn_stressx = Copo*ice_ocn_mag*( uo_i*np.cos(ThetaO) - vo_i*np.sin(ThetaO) )
    ocn_stressy = Copo*ice_ocn_mag*( vo_i*np.cos(ThetaO) + uo_i*np.sin(ThetaO) )
    
    cor_tiltx = -vo_i*rhohifc
    cor_tilty =  uo_i*rhohifc
    
    tot_mtm_x = atmo_stressx + ocn_stressx + cor_tiltx
    tot_mtm_y = atmo_stressy + ocn_stressy + cor_tilty

    # return np.hypot(tot_mtm_x, tot_mtm_y)**2
    return tot_mtm_x**2 +  tot_mtm_y**2

def misfit_fb_scaled(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):

    u_ice_ocn = uocn - uice
    v_ice_ocn = vocn - vice
    ice_ocn_mag = np.hypot(u_ice_ocn,v_ice_ocn)
    
    wnd_mag = np.hypot(uwnd,vwnd)

    misfit = misfit_fb(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO)
    
    #return misfit / (ice_ocn_mag**1 * wnd_mag**4)
    return misfit / (wnd_mag**2)

def misfit_a(uwnd,vwnd,uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO):
    uwnd_calc, vwnd_calc = wnd_mag(uice,vice,uocn,vocn,
	    Na2,Ro,ThetaA,ThetaO)

    norm_u = (uwnd - uwnd_calc)**2
    norm_v = (vwnd - vwnd_calc)**2

    return norm_u + norm_v

def misfit_a_io(uwnd,vwnd,uo_i,vo_i,
	    Na2,Ro,ThetaA,ThetaO):
    uwnd_calc, vwnd_calc = wnd_o_i_mag(uo_i,vo_i,
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
    #return misfit * ice_ocn_mag**2 

def misfit_o_i(uwnd,vwnd,uo_i,vo_i,
	    Na2,Ro,ThetaA,ThetaO):
    uo_i_calc, vo_i_calc = o_i_mag(uwnd,vwnd,
	    Na2,Ro,ThetaA,ThetaO)

    norm_u = (uo_i - uo_i_calc)**2
    norm_v = (vo_i - vo_i_calc)**2

    return norm_u + norm_v

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

def alpha(uwnd,vwnd,u_ice_ocn,v_ice_ocn):
    wnd_mag = np.hypot(uwnd,vwnd)
    ice_ocn_mag = np.hypot( u_ice_ocn, v_ice_ocn )

    ln_alpha = np.log( ice_ocn_mag ) - np.log( wnd_mag )
    alpha    = np.exp( ln_alpha )
    return alpha

#       theta = atan(v_ice_ocn/u_ice_ocn) - atan(vwnd/uwnd)


def theta(uwnd,vwnd,u_ice_ocn,v_ice_ocn):
    # tan_a = vwnd / uwnd
    # tan_io= v_ice_ocn / u_ice_ocn

    # theta = np.arctan2(vwnd,uwnd) - np.arctan2(v_ice_ocn,u_ice_ocn)
    theta = np.arctan2(v_ice_ocn*uwnd - u_ice_ocn*vwnd,u_ice_ocn*uwnd + v_ice_ocn*vwnd )
    return theta

def t_theta(uwnd,vwnd,u_ice_ocn,v_ice_ocn):
    #tan_a = vwnd / uwnd
    #tan_io= v_ice_ocn / u_ice_ocn

    top = v_ice_ocn*uwnd - u_ice_ocn*vwnd
    bot = u_ice_ocn*uwnd + v_ice_ocn*vwnd

    #t_theta = (tan_io - tan_a) / ( 1 + tan_io*tan_a )
    t_theta = top/bot
    return t_theta

def alpha_sig(uwnd,vwnd,u_ice_ocn,v_ice_ocn,
		wnd_sig,i_o_sig):
    alpha_use = alpha(uwnd,vwnd,u_ice_ocn,v_ice_ocn)
 
    wnd_mag = np.hypot(uwnd,vwnd)
    ice_ocn_mag = np.hypot( u_ice_ocn, v_ice_ocn )

    i_o_sig2 = i_o_sig**2

    i_o_sig2_o = i_o_sig2 / ice_ocn_mag**2
    wnd_sig2_o = wnd_sig**2 / wnd_mag**2

    alpha_sig2 = alpha_use**2 * (i_o_sig2_o + wnd_sig2_o)
    alpha_sig = np.sqrt(alpha_sig2)
    return alpha_sig

#       tricky angle error propogation
def t_theta_sig(uwnd,vwnd,u_ice_ocn,v_ice_ocn,
		wnd_sig,i_o_sig):


    wnd_mag = np.hypot(uwnd,vwnd)
    ice_ocn_mag = np.hypot( u_ice_ocn, v_ice_ocn )

    t_theta_use = t_theta(uwnd,vwnd,u_ice_ocn,v_ice_ocn)
    theta_use = theta(uwnd,vwnd,u_ice_ocn,v_ice_ocn)

    i_o_sig2 = i_o_sig**2

    sig_mag = np.sqrt(ice_ocn_mag**2*wnd_sig**2 + wnd_mag**2*i_o_sig2)

    over_1 = (v_ice_ocn*uwnd - u_ice_ocn*vwnd)**-2
    over_2 = (u_ice_ocn*uwnd + v_ice_ocn*vwnd)**-2

    #tan_a_sig = wnd_sig/uwnd * np.sqrt( 1 + tan_a**2 )
    #tan_io_sig= ice_ocn_sig/u_ice_ocn * np.sqrt( 1 + tan_io**2 )

    #atan_a_sig = tan_a_sig / (1 + tan_a**2)
    #atan_io_sig= tan_io_sig/ (1 + tan_io**2)

#     tan_sig_1 = (tan_io_sig**2 + tan_a_sig**2) / (tan_io - tan_a)**2
#     tan_sig_2 = tan_a**2*tan_io_sig**2 + tan_io**2*tan_a_sig**2
#     tan_sig_3 = tan_sig_2 / (1 + tan_a*tan_io)**2

    t_theta_sig = np.sqrt(over_1 + over_2)*sig_mag*np.abs(t_theta_use)

    return t_theta_sig 


def theta_sig(uwnd,vwnd,u_ice_ocn,v_ice_ocn,
		wnd_sig,i_o_sig):

    t_theta_use = t_theta(uwnd,vwnd,u_ice_ocn,v_ice_ocn)
    t_theta_sig_use = t_theta_sig(uwnd,vwnd,u_ice_ocn,v_ice_ocn,
				wnd_sig,i_o_sig)

    theta_sig = t_theta_sig_use / (1 + t_theta_use**2)
    return theta_sig

##def t_theta_sig(uwnd,vwnd,uice,vice,uocn,vocn,
#		wnd_sig,ice_sig,ocn_sig):
#
#    u_ice_ocn = uocn - uice
#    v_ice_ocn = vocn - vice
#
#    tan_a = vwnd / uwnd
#    tan_io= v_ice_ocn / u_ice_ocn
#
#    theta_use = theta(uwnd,vwnd,uice,vice,uocn,vocn)
#
#    ice_ocn_sig2 = ice_sig**2 + ocn_sig**2
#    ice_ocn_sig  = np.sqrt(ice_ocn_sig2)
#
#    tan_a_sig = wnd_sig/uwnd * np.sqrt( 1 + tan_a**2 )
#    tan_io_sig= ice_ocn_sig/u_ice_ocn * np.sqrt( 1 + tan_io**2 )
#
#    atan_a_sig = tan_a_sig / (1 + tan_a**2)
#    atan_io_sig= tan_io_sig/ (1 + tan_io**2)
#
#    tan_sig_1 = (tan_io_sig**2 + tan_a_sig**2) / (tan_io - tan_a)**2
#    tan_sig_2 = tan_a**2*tan_io_sig**2 + tan_io**2*tan_a_sig**2
#    tan_sig_3 = tan_sig_2 / (1 + tan_a*tan_io)**2
#
#    t_theta_sig = theta_use* np.sqrt(tan_sig_1 + tan_sig_3)
#
#    return t_theta_sig 
#

def misfit_alpha_theta(alpha, theta, hifc_ua,
                Na2,Ro,ThetaA,ThetaO):
     alpha_m, theta_m = ice_alpha(hifc_ua,
            			     Na2,Ro,ThetaA,ThetaO)	
     misfit_alpha = ( alpha - alpha_m )**2
     misfit_theta = ( theta - theta_m )**2

     return misfit_alpha, misfit_theta

def force_bal(uwnd,vwnd,ui_o,vi_o,
	    Capa,Copo,rhohifc,ThetaA,ThetaO):
    # now try back in the mtm balance
    wnd_mag = np.hypot(uwnd,vwnd)

    ice_ocn_mag = np.hypot(ui_o,vi_o)

    ThetaAO = ThetaA - ThetaO

    atmo_stressx = Capa*wnd_mag*( uwnd*np.cos(ThetaA) - vwnd*np.sin(ThetaA) )
    atmo_stressy = Capa*wnd_mag*( vwnd*np.cos(ThetaA) + uwnd*np.sin(ThetaA) )
    
    ocn_stressx = Copo*ice_ocn_mag*( ui_o*np.cos(ThetaO) - vi_o*np.sin(ThetaO) )
    ocn_stressy = Copo*ice_ocn_mag*( vi_o*np.cos(ThetaO) + ui_o*np.sin(ThetaO) )
    
    cor_tiltx = -vi_o*rhohifc
    cor_tilty =  ui_o*rhohifc
    

    return atmo_stressx, atmo_stressy, ocn_stressx, ocn_stressy, cor_tiltx, cor_tilty


def wnd_mag_geo(uo_i,vo_i,ugeo,vgeo,
	    Na2,Ro,ThetaA,ThetaO):

    o_i_mag = np.hypot(uo_i,vo_i)

    y11 = -(o_i_mag/Na2) * (np.cos(ThetaO-ThetaA)*uo_i - np.sin(ThetaO-ThetaA)*vo_i)
    y12 = -(o_i_mag/Na2) * (np.cos(ThetaO-ThetaA)*vo_i + np.sin(ThetaO-ThetaA)*uo_i)

    y21 = -(Ro/Na2) * (np.sin(ThetaA)*ugeo - np.cos(ThetaA)*vgeo)
    y22 = -(Ro/Na2) * (np.sin(ThetaA)*vgeo + np.cos(ThetaA)*ugeo)

    y1 = y11 + y21
    y2 = y12 + y22
    y_mag = np.sqrt(np.hypot(y1,y2))
    y_arg = np.arctan2(y2,y1)

    uwnd = y_mag*np.cos(y_arg)
    vwnd = y_mag*np.sin(y_arg)

    return uwnd, vwnd


def o_i_mag_geo(uwnd,vwnd,ugeo,vgeo,
	    Na2,Ro,ThetaA,ThetaO):

    wnd_mag = np.hypot(uwnd,vwnd)

    y11 = -wnd_mag*Na2 * (np.cos(ThetaA-ThetaO)*uwnd - np.sin(ThetaA-ThetaO)*vwnd)
    y12 = -wnd_mag*Na2 * (np.cos(ThetaA-ThetaO)*vwnd + np.sin(ThetaA-ThetaO)*uwnd)

    y21 = -Ro * (-np.sin(-ThetaO)*ugeo - np.cos(-ThetaO)*vgeo)
    y22 = -Ro * (-np.sin(-ThetaO)*vgeo + np.cos(-ThetaO)*ugeo)

    y1 = y11 + y21
    y2 = y12 + y22
    y_mag = np.sqrt(np.hypot(y1,y2))
    y_arg = np.arctan2(y2,y1)

    uo_i = y_mag*np.cos(y_arg)
    vo_i = y_mag*np.sin(y_arg)

    return uo_i, vo_i

def geo_mag_geo(uwnd,vwnd,uo_i,vo_i,
	    Na2,Ro,ThetaA,ThetaO):

    wnd_mag = np.hypot(uwnd,vwnd)
    o_i_mag = np.hypot(uo_i,vo_i)

    y11 = (wnd_mag*Na2) * (-np.sin(ThetaA)*uwnd - np.cos(ThetaA)*vwnd)/Ro
    y12 = (wnd_mag*Na2) * (-np.sin(ThetaA)*vwnd + np.cos(ThetaA)*uwnd)/Ro
    y11 = (wnd_mag*Na2) * (-np.sin(ThetaA)*uwnd - np.cos(ThetaA)*vwnd)
    y12 = (wnd_mag*Na2) * (-np.sin(ThetaA)*vwnd + np.cos(ThetaA)*uwnd)

    y21 = (o_i_mag) * (-np.sin(ThetaO)*uo_i - np.cos(ThetaO)*vo_i)/Ro
    y22 = (o_i_mag) * (-np.sin(ThetaO)*vo_i + np.cos(ThetaO)*uo_i)/Ro
    y21 = (o_i_mag) * (-np.sin(ThetaO)*uo_i - np.cos(ThetaO)*vo_i)
    y22 = (o_i_mag) * (-np.sin(ThetaO)*vo_i + np.cos(ThetaO)*uo_i)

    y1 = y11 + y21
    y2 = y12 + y22
    y_mag = np.hypot(y1,y2)
    y_arg = np.arctan2(y2,y1)

    # return y1, y2
    return y1/Ro, y2/Ro
    # return y_mag*np.cos(y_arg), y_mag*np.sin(y_arg)

def misfit_wnd_geo(uwnd,vwnd,uo_i,vo_i,ugeo,vgeo,
	    Na2,Ro,ThetaA,ThetaO):
    uwnd_calc, vwnd_calc = wnd_mag_geo(uo_i,vo_i,ugeo,vgeo,
	    Na2,Ro,ThetaA,ThetaO)

    norm_u = (uwnd - uwnd_calc)**2
    norm_v = (vwnd - vwnd_calc)**2

    return norm_u + norm_v

def misfit_o_i_geo(uwnd,vwnd,uo_i,vo_i,ugeo,vgeo,
	    Na2,Ro,ThetaA,ThetaO):
    uo_i_calc, vo_i_calc = o_i_mag_geo(uwnd,vwnd,ugeo,vgeo,
	    Na2,Ro,ThetaA,ThetaO)

    norm_u = (uo_i - uo_i_calc)**2
    norm_v = (vo_i - vo_i_calc)**2

    return norm_u + norm_v

def misfit_geo_geo(uwnd,vwnd,uo_i,vo_i,ugeo,vgeo,
	    Na2,Ro,ThetaA,ThetaO):
    ugeo_calc, vgeo_calc = geo_mag_geo(uwnd,vwnd,uo_i,vo_i,
	    Na2,Ro,ThetaA,ThetaO)

    norm_u = (ugeo - ugeo_calc)**2
    norm_v = (vgeo - vgeo_calc)**2

    return norm_u + norm_v


def force_bal_geo(uwnd,vwnd,uo_i,vo_i,ugeo,vgeo,
	    Capa,Copo,rhohifc,ThetaA,ThetaO):

    wnd_mag = np.hypot(uwnd,vwnd)
    o_i_mag = np.hypot(uo_i,vo_i)

    y11 = wnd_mag*Capa * (np.cos(ThetaA)*uwnd - np.sin(ThetaA)*vwnd)
    y12 = wnd_mag*Capa * (np.cos(ThetaA)*vwnd + np.sin(ThetaA)*uwnd)

    y21 = o_i_mag*Copo * (np.cos(ThetaO)*uo_i - np.sin(ThetaO)*vo_i)
    y22 = o_i_mag*Copo * (np.cos(ThetaO)*vo_i + np.sin(ThetaO)*uo_i)

    y31 = -vgeo*rhohifc
    y32 =  ugeo*rhohifc


    return y11, y12, y21, y22, y31, y32 


def misfit_fb_geo(uwnd,vwnd,uo_i,vo_i,ugeo,vgeo,
	    Capa,Copo,rhohifc,ThetaA,ThetaO):

    wnd_mag = np.hypot(uwnd,vwnd)
    o_i_mag = np.hypot(uo_i,vo_i)

    y11 = wnd_mag*Capa * (np.cos(ThetaA)*uwnd - np.sin(ThetaA)*vwnd)
    y12 = wnd_mag*Capa * (np.cos(ThetaA)*vwnd + np.sin(ThetaA)*uwnd)

    y21 = o_i_mag*Copo * (np.cos(ThetaO)*uo_i - np.sin(ThetaO)*vo_i)
    y22 = o_i_mag*Copo * (np.cos(ThetaO)*vo_i + np.sin(ThetaO)*uo_i)

    y31 = -vgeo*rhohifc
    y32 =  ugeo*rhohifc

    mtm_x = y11 + y21 + y31
    mtm_y = y12 + y22 + y32

    return mtm_x**2 + mtm_y**2


