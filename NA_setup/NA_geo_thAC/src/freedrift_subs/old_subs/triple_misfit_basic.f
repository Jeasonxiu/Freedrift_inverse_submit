c-------------------------------------------------------------------------
c
c	Subroutine triple_misfit - calculates a misfit value between 
c				observed data and predicted data
c				three times, calculating vels for
c				wind, ice, ocn from the other two
c				Also adds the FB residual
c
c	Note: Calls no other routines
c
c------------------------------------------------------------------------
c 
        subroutine triple_misfit(
     &          model, 
     & 		uwind, vwind, wind_mag, wind_sig,
     & 		uice, vice, ice_mag, ice_sig,
     & 		uocn, vocn, ocn_mag, ocn_sig,
     & 		hifc_scale, wind_mag_scale,
     &          fb_ratio, misfitval )
c
c
c
        real            model(4)
     	real   		uwind, vwind, wind_mag, wind_sig
     	real   		uice, vice, ice_mag, ice_sig
     	real   		uocn, vocn, ocn_mag, ocn_sig
     	real   		u_ice_ocn, v_ice_ocn, ice_ocn_mag
     	real   		uwind_calc, vwind_calc, wind_mag_calc
     	real   		uice_calc, vice_calc, ice_mag_calc
     	real   		uocn_calc, vocn_calc, ocn_mag_calc
     	real   		hifc_scale, wind_mag_scale
     	real            fb_ratio, misfitval 
     	real            misfit_i, misfit_o
     	real            misfit_a, misfit_fb
        real    	Eq11, Eq12, Eq13
        real    	Eq21, Eq22, Eq23
        real    	thA, thO
        real    	Ro, Na2
        real    	fi, dfi, yi, y1, x, eps, pi
        real            atmo_stressx, ocn_stressx, cor_tiltx, tot_mtm_x
        real            atmo_stressy, ocn_stressy, cor_tilty, tot_mtm_y

        parameter       ( pi = 3.141259E0 )
c 
        integer        conv
c
c						misfit between observed and
c						predicted
c model contains:
c 1: Na2Ua - scaled by wind_mag to get the correct Ua
c 2: Rossby radius, Na independant - scaled by hifc
c 3: ThetaO
c 4: ThetaA - ThetaO

c 	write(*,*)' model = ',model

c note scaling for model(1,2)
c when removing Ua from Na2Ua, Na2 = Na2Ua/(wind_mag*wind_mag_scale)
c When using Na2Ua, Na2Ua = Na2Ua/wind_mag_scale
c (wind_mag_scale = wind_mag/wind_mag_* )
c When using Ro, Ro = Ro/hifc_scale
c (hifc_scale = hifc/hifc_*)

	thA = model(3) + model(4)
	thO = model(3)
	Ro  = model(2)/hifc_scale
	Na2 = model(1)/(wind_mag*wind_mag_scale)
c 	write(*,*)' Na2 = ',Na2
c 	write(*,*)' Ro  = ',Ro

c constants used to solve free drift balance
        Eq11 = TAN(thO)
        Eq12 = Ro / wind_mag
        Eq13 = COS(thO)

        Eq21 = 2*SIN(thO)*Ro/wind_mag
        Eq22 = (Ro**2)/(wind_mag**2)
        Eq23 = Na2**2

c now to try newtons method
c easy to start thetaI = 0, alpha = Na, alpha strictly < Na
c f  =   y^4 +   Eq21 y^3 +   Eq22 y^2 - Eq23
c f' = 4 y^3 + 3 Eq21 y^2 + 2 Eq22 y

        y1 = SQRT(Na2)
        yi = y1
        eps = 1.0E-7
        conv = 0
        do while (conv == 0)
            fi =   yi**4 +   Eq21*(yi**3) +   Eq22*(yi**2) - Eq23
            dfi= 4*yi**3 + 3*Eq21*(yi**2) + 2*Eq22*(yi)
            yi = yi - fi/dfi
c             write(*,*)' yi = ',yi,' fi/dfi = ',fi/dfi
            if(ABS(fi/dfi) < eps)then
                conv = 1
            end if
        end do
	
c	c x = ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))

        x = - thA + ATAN(Eq11 + Eq12/(yi*Eq13))


        uice_calc = yi*(COS(-x)*uwind - SIN(-x)*vwind) + uocn
        vice_calc = yi*(COS(-x)*vwind + SIN(-x)*uwind) + vocn

	misfit_i = (uice - uice_calc)**2 + (vice - vice_calc)**2

c        ice_mag_calc = SQRT(uice_calc**2 + vice_calc**2)
c
c         write(*,*)' Ice original   = ',uice,vice
c         write(*,*)' Ice calculated = ',uice_calc,vice_calc

        uocn_calc = uice - yi*(COS(-x)*uwind - SIN(-x)*vwind)
        vocn_calc = vice - yi*(COS(-x)*vwind + SIN(-x)*uwind)

	misfit_o = (uocn - uocn_calc)**2 + (vocn - vocn_calc)**2

c        ocn_mag_calc = SQRT(uocn_calc**2 + vocn_calc**2)
c
c         write(*,*)' Ocn original   = ',uocn,vocn
c         write(*,*)' Ocn calculated = ',uocn_calc,vocn_calc

c Atmo needs a special solver Na2Ua constant, Eq** change per interation

c easy to start thetaI = 0, alpha = Na, alpha strictly < Na
c f  =   y^4 +   Eq21 y^3 +   Eq22 y^2 - Eq23
c f' = 4 y^3 + 3 Eq21 y^2 + 2 Eq22 y

c constants that remain constant
            Eq11 = TAN(thO)
            Eq13 = COS(thO)

            Eq23 = Na2**2
        yi = y1
        eps = 1.0E-7
        conv = 0
        wind_mag_calc = wind_mag
        do while (conv == 0)
c constants that change with iteration
            Eq12 = Ro/wind_mag_calc

            Eq21 = 2*SIN(thO)*Ro/wind_mag_calc
            Eq22 = (Ro**2)/(wind_mag_calc**2)

            fi =   yi**4 +   Eq21*(yi**3) +   Eq22*(yi**2) - Eq23
            dfi= 4*yi**3 + 3*Eq21*(yi**2) + 2*Eq22*(yi)
            yi = yi - fi/dfi

c need to put the solution through to the wind mag to properly converge
            x = - thA + ATAN(Eq11 + Eq12/(yi*Eq13))
c             write(*,*)' yi = ',yi,' fi/dfi = ',fi/dfi
c             write(*,*)' x  = ',x
            uwind_calc = (COS(x)*(uice - uocn) 
     &			- SIN(x)*(vice - vocn))/yi
            vwind_calc = (COS(x)*(vice - vocn) 
     &			+ SIN(x)*(uice - uocn))/yi
            wind_mag_calc = SQRT(uwind_calc**2 + vwind_calc**2)

            if(ABS(fi/dfi) < eps)then
                conv = 1
            end if
        end do

c         write(*,*)' Wind original   = ',uwind,vwind
c         write(*,*)' Wind calculated = ',uwind_calc,vwind_calc

	misfit_a = (uwind - uwind_calc)**2 + (vwind - vwind_calc)**2

c fb residual misfit

        u_ice_ocn = uocn - uice
        v_ice_ocn = vocn - vice
        ice_ocn_mag = SQRT(u_ice_ocn**2.0E0 + v_ice_ocn**2.0E0)

c components indiidually
c atmo stress component

        atmo_stressx = model(1)/wind_mag_scale
     &*( uwind*COS(thA) - vwind*SIN(thA) )
        atmo_stressy = model(1)/wind_mag_scale
     &*( vwind*COS(thA) + uwind*SIN(thA) )

c ocean stress component

        ocn_stressx  = ice_ocn_mag
     &*( u_ice_ocn*COS(thO) - v_ice_ocn*SIN(thO) )
        ocn_stressy  = ice_ocn_mag
     &*( v_ice_ocn*COS(thO) + u_ice_ocn*SIN(thO) )

c coriolis acceleration + ocean tilt (geostrophic) component

        cor_tiltx = -v_ice_ocn*Ro
        cor_tilty =  u_ice_ocn*Ro

c add it all up !!

        tot_mtm_x = atmo_stressx + ocn_stressx + cor_tiltx
        tot_mtm_y = atmo_stressy + ocn_stressy + cor_tilty

c simple misfit  - just the total momentum, that should be 0

        misfit_fb = tot_mtm_x**2.0E0 + tot_mtm_y**2.0E0

c 	write(*,*)' misfit_a = ',misfit_a
c 	write(*,*)' misfit_i = ',misfit_i
c 	write(*,*)' misfit_o = ',misfit_o
c 	write(*,*)' misfit_fb= ',misfit_fb

	misfitval = misfit_i/ice_sig  + misfit_o/ocn_sig +
     &              misfit_a/wind_sig + fb_ratio*misfit_fb

 
        return
        end
