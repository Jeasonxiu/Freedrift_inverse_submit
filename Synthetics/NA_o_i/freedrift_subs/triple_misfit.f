c-------------------------------------------------------------------------
c
c	Subroutine triple_misfit - calculates a misfit value between 
c				observed data and predicted data
c				three times, calculating vels for
c				wind, ice, i_o from the other two
c				Also adds the FB residual
c
c	Note: Calls no other routines
c
c------------------------------------------------------------------------
c 
        subroutine triple_misfit(
     &          model, 
     & 		uwind, vwind, wind_mag, wind_sig,
     & 		ui_o, vi_o, i_o_mag, i_o_sig,
     & 		hifc, 
     & 		beta_a, beta_i, beta_fb,
     &          fb_ratio, misfitval )
c
c
c
        real            model(4)
     	real   		uwind, vwind, wind_mag, wind_sig
     	real   		ui_o, vi_o, i_o_mag, i_o_sig
     	real   		uwind_calc, vwind_calc, wind_mag_calc
     	real   		ui_o_calc, vi_o_calc, i_o_mag_calc
     	real   		hifc
        real 		beta_a, beta_i, beta_fb
     	real            fb_ratio, misfitval 
     	real            misfit_o
     	real            misfit_a, misfit_fb
        real    	Eq11, Eq12, Eq13
        real    	Eq21, Eq22, Eq23
        real    	thA, thO, thAO
        real    	Ro, Na2
        real    	fi, dfi, yi, y1, x, eps, pi, Na2_scale
        real            atmo_stressx, ocn_stressx, cor_tiltx, tot_mtm_x
        real            atmo_stressy, ocn_stressy, cor_tilty, tot_mtm_y

        parameter       ( pi = 3.141259E0 )
        parameter       ( Na2_scale = 0.01 )
c 
        integer        conv, it
c
c						misfit between observed and
c						predicted
c model contains:
c 1: Na2- no longer scaled by wind_mag to get the correct Ua
c 2: Rossby radius, Na independant - independent of hifc
c 3: ThetaO
c 4: ThetaA - ThetaO

c	write(*,*)' model = ',model

c When using Ro, Ro = Ro*hifc

	thA = model(3) + model(4)
	thO = model(3)
	thAO= model(4)
	Ro  = model(2)*hifc
	Na2 = model(1)
c	write(*,*)' Na2 = ',Na2
c	write(*,*)' Ro  = ',Ro

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
        eps = 1.0E-9
        conv = 0
	it = 0
        do while (conv == 0)
            fi =   yi**4 +   Eq21*(yi**3) +   Eq22*(yi**2) - Eq23
            dfi= 4*yi**3 + 3*Eq21*(yi**2) + 2*Eq22*(yi)
            yi = yi - fi/dfi
c            write(*,*)' yi = ',yi,' fi/dfi = ',fi/dfi
	    it = it + 1
c            write(*,*)' it = ',it,' eps = ',eps
            if(ABS(fi/dfi) < eps)then
                conv = 1
            end if
            if(it < 10) then
                eps = 1.0E-9
	    else if(it < 50) then
                eps = 1.0E-8
	    else if(it < 200) then
                eps = 1.0E-7
	    else if(it < 500) then
                eps = 1.0E-3
	    else if(it < 1000) then
                eps = 1.0E-2
            end if
        end do
	
c	c x = ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))

        x = - thA + ATAN(Eq11 + Eq12/(yi*Eq13))


        ui_o_calc = yi*(COS(-x)*uwind - SIN(-x)*vwind) 
        vi_o_calc = yi*(COS(-x)*vwind + SIN(-x)*uwind)

	misfit_o = (ui_o - ui_o_calc)**2 + (vi_o - vi_o_calc)**2
c       scaling misfit to get it order 1
	misfit_o = misfit_o * beta_i
c	misfit_o = misfit_o 

c        i_o_mag_calc = SQRT(ui_o_calc**2 + vi_o_calc**2)
c
c         write(*,*)' Ocn original   = ',ui_o,vi_o
c         write(*,*)' Ocn calculated = ',ui_o_calc,vi_o_calc

c Atmo needs a special solver Na2Ua constant, Eq** change per interation

c easy to start thetaI = 0, alpha = Na, alpha strictly < Na
c f  =   y^4 +   Eq21 y^3 +   Eq22 y^2 - Eq23
c f' = 4 y^3 + 3 Eq21 y^2 + 2 Eq22 y

c constants that remain constant
            Eq11 = TAN(thO)
            Eq13 = COS(thO)

            Eq23 = Na2**2
        yi = y1
        eps = 1.0E-9
        conv = 0
	it = 0
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
c           write(*,*)' yi = ',yi,' fi/dfi = ',fi/dfi
c            write(*,*)' x  = ',x
            uwind_calc = (COS(x)*( ui_o) 
     &			- SIN(x)*( vi_o))/yi
            vwind_calc = (COS(x)*( vi_o) 
     &			+ SIN(x)*( ui_o))/yi
            wind_mag_calc = SQRT(uwind_calc**2 + vwind_calc**2)

	    it = it + 1
c           write(*,*)' it = ',it,' eps = ',eps
            if(ABS(fi/dfi) < eps)then
                conv = 1
            end if
            if(it < 10) then
                eps = 1.0E-9
	    else if(it < 50) then
                eps = 1.0E-8
	    else if(it < 200) then
                eps = 1.0E-7
	    else if(it < 500) then
                eps = 1.0E-3
	    else if(it < 1000) then
                eps = 1.0E-2
            end if
        end do

c        write(*,*)' Wind original   = ',uwind,vwind
c        write(*,*)' Wind calculated = ',uwind_calc,vwind_calc

	misfit_a = (uwind - uwind_calc)**2 + (vwind - vwind_calc)**2
c       scaling misfit to get it order 1
c	misfit_a = misfit_a * ice_i_o_mag**2 / wind_mag**2
	misfit_a = misfit_a * beta_a

c fb residual misfit

c components indiidually
c atmo stress component

        atmo_stressx = Na2*wind_mag
     &*( uwind*COS(thAO) - vwind*SIN(thAO) )
        atmo_stressy = Na2*wind_mag
     &*( vwind*COS(thAO) + uwind*SIN(thAO) )

c ocean stress component

        ocn_stressx  = i_o_mag
     &*( u_ice_i_o )
        ocn_stressy  = i_o_mag
     &*( v_ice_i_o )

c coriolis acceleration + ocean tilt (geostrophic) component

        cor_tiltx = -(vi_o*COS(-thO) - ui_o*SIN(-thO))*Ro
        cor_tilty =  (ui_o*COS(-thO) + vi_o*SIN(-thO))*Ro

c add it all up !!

        tot_mtm_x = atmo_stressx + ocn_stressx + cor_tiltx
        tot_mtm_y = atmo_stressy + ocn_stressy + cor_tilty

c simple misfit  - just the total momentum, that should be 0

        misfit_fb = tot_mtm_x**2.0E0 + tot_mtm_y**2.0E0
  	misfit_fb = misfit_fb * beta_fb

c 	write(*,*)' misfit_a = ',misfit_a
c	 write(*,*)' misfit_i = ',misfit_i
c 	write(*,*)' misfit_o = ',misfit_o
c 	write(*,*)' misfit_fb= ',misfit_fb

c 	write(*,*)' Error scaled......'
c 	write(*,*)' misfit_a = ',misfit_a/wind_sig
c 	write(*,*)' misfit_i = ',misfit_i/ice_sig/ice_mag
c 	write(*,*)' misfit_o = ',misfit_o/i_o_sig
c	write(*,*)' misfit_fb= ',misfit_fb*fb_ratio
c 
	misfitval = misfit_o/i_o_sig**2 +
     &              misfit_a/wind_sig**2 + misfit_fb/fb_ratio**2

c note that the errors are dimensional now
 
        return
        end
