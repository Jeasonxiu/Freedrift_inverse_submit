C-------------------------------------------------------------------------
c
c	Subroutine triple_misfit - calculates a misfit value between 
c				observed data and predicted data
c				three times, calculating vels for
c				wind, ice, o_i from the other two
c				Also adds the FB residual
c				This new version allows for full circle angles
c				Misift set within the loop, 
c				default value = 10*vector
c
c	Note: Calls no other routines
c
c------------------------------------------------------------------------
c 
        subroutine triple_misfit(
     &          model, 
     & 		uwind, vwind, wind_mag, wind_sig,
     & 		uo_i, vo_i, o_i_mag, o_i_sig,
     & 		hifc, 
     & 		beta_a, beta_i, beta_fb,
     &          fb_ratio, misfitval, bflag )
c
c
c
        real            model(4)
     	real   		uwind, vwind, wind_mag, wind_sig
     	real   		uo_i, vo_i, o_i_mag, o_i_sig
     	real   		uwind_calc, vwind_calc, wind_mag_calc
     	real   		uo_i_calc, vo_i_calc, o_i_mag_calc
     	real   		hifc
        real 		beta_a, beta_i, beta_fb
     	real            fb_ratio, misfitval 
     	real            misfit_o
     	real            misfit_a, misfit_fb
        real    	Eq11, Eq12, Eq13
        real    	Eq21, Eq22, Eq23
        real    	thA, thO, thAO
        real    	Ro, Na2, Capa, Copo, rhoi, bigp
        real(8)    	fi, dfi, yi, y1, x, eps, pi, Na2_scale
        real(8)    	y0, yi1, yin, mi, chng, fi1
        real            atmo_stressx, ocn_stressx, cor_tiltx, tot_mtm_x
        real            atmo_stressy, ocn_stressy, cor_tilty, tot_mtm_y

        parameter       ( pi = 3.141259E0 )
        parameter       ( rhoi = 917.0E0 )
        parameter       ( bigp = 1.0e7 )
        parameter       ( Na2_scale = 0.01 )
c 
        integer        conv, it, bflag

c Set up no convergence flag
	bflag = 0

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

	Copo = rhoi / model(2)
	Capa = model(1) * Copo

c constants used to solve free drift balance
        Eq11 = SIN(thO)
c       Eq11 = TAN(thO)
        Eq12 = Ro / wind_mag
        Eq13 = COS(thO)

        Eq21 = 2*SIN(thO)*Ro/wind_mag
        Eq22 = (Ro**2)/(wind_mag**2)
        Eq23 = Na2**2


c now to try newtons method
c easy to start thetaI = 0, alpha = Na, alpha strictly < Na
c f  =   y^4 +   Eq21 y^3 +   Eq22 y^2 - Eq23
c f' = 4 y^3 + 3 Eq21 y^2 + 2 Eq22 y
c Now with a two point method to remove oscillation - check paper citations

	y0 = 0.04E0
    	y1 = 0.035E0
    	yi1 = y0
    	yi  = y1
c   	# prep for new method
    	fi =   yi1**4 +   Eq21*(yi1**3) +   Eq22*(yi1**2) - Eq23
    	mi  =0.0
        chng = 0.0
    	eps = 1.0D-10
    	conv = 0
    	its = 0
    	do while (conv == 0)

    	    fi1 = fi
    	    fi  =   yi**4 +   Eq21*(yi**3) +   Eq22*(yi**2) - Eq23
    	    dfi = 4*yi**3 + 3*Eq21*(yi**2) + 2*Eq22*(yi) 
    	    mi  = (fi - fi1)/(yi - yi1 + 1.0D-9)
    	    mi  = mi / (dfi + 1.0D-9)
    	    chng= (yi1 - yi)/(1.0D0 - mi*fi/fi1 + 1.0D-9)
    	    yin = yi1 - chng
    	    yi1 = yi
    	    yi  = yin

c            write(*,*)' yi = ',yi,' fi/dfi = ',fi/dfi
c           write(*,*)' yi = ',yi,' yi1  = ',yi1 ,fi,fi1,dfi,mi,ABS(chng),eps
	    it = it + 1
c            write(*,*)' it = ',it,' eps = ',eps
            if(ABS(chng) < eps)then
                conv = 1
c       	x = - thA + ATAN(Eq11 + Eq12/(yi*Eq13))
    		x =   thA - ATAN2(Eq11*yi + Eq12,yi*Eq13)

        	uo_i_calc = -( yi*(COS(-x)*uwind - SIN(-x)*vwind) )
        	vo_i_calc = -( yi*(COS(-x)*vwind + SIN(-x)*uwind) )
            end if
            if(it < 10) then
                eps = 1.0D-9
	    else if(it < 50) then
                eps = 1.0D-8
	    else if(it < 200) then
                eps = 1.0D-7
	    else if(it < 500) then
                eps = 1.0D-6
	    else if(it < 1000) then
                eps = 1.0D-4
	    else if(it < 2000) then
                eps = 1.0D-2
	    else 
                conv = 1
		bflag = bflag + 1
        	uo_i_calc = bigp * uo_i
        	vo_i_calc = bigp * vo_i
c          write(*,*)'No convergence -- ocn',mi,(fi - fi1)/(yi - yi1 + 0.0E-2)
c         write(*,*)'No convergence -- ocn',fi,fi1,dfi
            end if
        end do
	
c	c x = ThetaA + np.arctan(Eq11 + Eq12/(yi*Eq13))


	misfit_o = (uo_i - uo_i_calc)**2 + (vo_i - vo_i_calc)**2
c       scaling misfit to get it order 1
	misfit_o = misfit_o * beta_i
c	misfit_o = misfit_o 

c        o_i_mag_calc = SQRT(uo_i_calc**2 + vo_i_calc**2)
c
c         write(*,*)' Ocn original   = ',uo_i,vo_i
c         write(*,*)' Ocn calculated = ',uo_i_calc,vo_i_calc

c Atmo needs a special solver Na2Ua constant, Eq** change per interation

c easy to start thetaI = 0, alpha = Na, alpha strictly < Na
c f  =   y^4 +   Eq21 y^3 +   Eq22 y^2 - Eq23
c f' = 4 y^3 + 3 Eq21 y^2 + 2 Eq22 y

c constants that remain constant
c           Eq11 = TAN(thO)
            Eq11 = SIN(thO)
            Eq13 = COS(thO)

            Eq23 = Na2**2
c prepping the non-constant constants...
	wind_mag_calc = 1.0E0
        Eq12 = Ro/wind_mag_calc

        Eq21 = 2*SIN(thO)*Ro/wind_mag_calc
        Eq22 = (Ro**2)/(wind_mag_calc**2)

	y0 = 0.04
    	y1 = 0.035
    	yi1 = y0
    	yi  = y1
c   	# prep for new method
    	fi =   yi1**4 +   Eq21*(yi1**3) +   Eq22*(yi1**2) - Eq23
    	mi  =0.0
        eps = 1.0D-10
        conv = 0
	it = 0
        wind_mag_calc = wind_mag
        do while (conv == 0)
c constants that change with iteration
            Eq12 = Ro/wind_mag_calc

            Eq21 = 2*SIN(thO)*Ro/wind_mag_calc
            Eq22 = (Ro**2)/(wind_mag_calc**2)

    	    fi1 = fi
    	    fi  =   yi**4 +   Eq21*(yi**3) +   Eq22*(yi**2) - Eq23
    	    dfi = 4*yi**3 + 3*Eq21*(yi**2) + 2*Eq22*(yi) 
    	    mi  = (fi - fi1)/(yi - yi1 + 1.0D-12)
    	    mi  = mi / (dfi + 1.0D-12)
    	    chng= (yi1 - yi)/(1.0D0 - mi*fi/fi1 + 1.0D-12)
    	    yin = yi1 - chng
    	    yi1 = yi
    	    yi  = yin

c need to put the solution through to the wind mag to properly converge
c           x = - thA + ATAN(Eq11 + Eq12/(yi*Eq13))
c   	    x = - thA + ATAN2(Eq11*(yi*Eq13) + Eq12,(yi*Eq13))
    	    x =   thA - ATAN2(Eq11*yi + Eq12,yi*Eq13)
c           write(*,*)' yi = ',yi,' fi/dfi = ',fi/dfi
c            write(*,*)' x  = ',x
            uwind_calc = -(COS(x)*( uo_i) 
     &			- SIN(x)*( vo_i))/yi
            vwind_calc = -(COS(x)*( vo_i) 
     &			+ SIN(x)*( uo_i))/yi
            wind_mag_calc = SQRT(uwind_calc**2 + vwind_calc**2)

	    it = it + 1
c           write(*,*)' it = ',it,' eps = ',eps
            if(ABS(chng) < eps)then
                conv = 1
            end if
            if(it < 10) then
                eps = 1.0D-9
	    else if(it < 50) then
                eps = 1.0D-8
	    else if(it < 200) then
                eps = 1.0D-7
	    else if(it < 500) then
                eps = 1.0D-6
	    else if(it < 1000) then
                eps = 1.0D-5
	    else if(it < 2000) then
                eps = 1.0D-4
	    else if(it < 3000) then
                eps = 1.0D-3
	    else
                conv = 1
		bflag = bflag + 1
        	uwind_calc = bigp * uwind
        	vwind_calc = bigp * vwind
c         write(*,*)'No convergence -- wnd',mi,(fi - fi1)/(yi - yi1 + 0.0E-2)
            end if
        end do

c        write(*,*)' Wind original   = ',uwind,vwind
c        write(*,*)' Wind calculated = ',uwind_calc,vwind_calc

	misfit_a = (uwind - uwind_calc)**2 + (vwind - vwind_calc)**2
c       scaling misfit to get it order 1
c	misfit_a = misfit_a * ice_o_i_mag**2 / wind_mag**2
	misfit_a = misfit_a * beta_a

c fb residual misfit

c components indiidually
 	atmo_stressx = Capa*wind_mag*( uwind*COS(thA) - vwind*SIN(thA) )
        atmo_stressy = Capa*wind_mag*( vwind*COS(thA) + uwind*SIN(thA) )

c ocean stress component

        ocn_stressx  = Copo*o_i_mag *( uo_i*COS(thO) - vo_i*SIN(thO) )
        ocn_stressy  = Copo*o_i_mag *( vo_i*COS(thO) + uo_i*SIN(thO) )

c coriolis acceleration + ocean tilt (geostrophic) component

        cor_tiltx = -vo_i*rhoi*hifc
        cor_tilty =  uo_i*rhoi*hifc

cc atmo stress component
c
c        atmo_stressx = Capa*wind_mag
c     &*( uwind*COS(thAO) - vwind*SIN(thAO) )
c        atmo_stressy = Capa*wind_mag
c     &*( vwind*COS(thAO) + uwind*SIN(thAO) )
c
cc ocean stress component
c
c        ocn_stressx  = Copo*o_i_mag
c     &*( u_ice_o_i )
c        ocn_stressy  = Copo*o_i_mag
c     &*( v_ice_o_i )
c
cc coriolis acceleration + ocean tilt (geostrophic) component
c
c        cor_tiltx = -(vo_i*COS(-thO) - uo_i*SIN(-thO))*Ro
c        cor_tilty =  (uo_i*COS(-thO) + vo_i*SIN(-thO))*Ro
c
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
c 	write(*,*)' misfit_o = ',misfit_o/o_i_sig
c	write(*,*)' misfit_fb= ',misfit_fb*fb_ratio
c 
	misfitval = misfit_o/o_i_sig**2 +
     &              misfit_a/wind_sig**2 + misfit_fb/fb_ratio**2

c note that the errors are dimensional now
 
        return
        end
