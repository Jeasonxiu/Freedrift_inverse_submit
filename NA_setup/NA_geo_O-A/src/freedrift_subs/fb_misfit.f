c-------------------------------------------------------------------------
c
c	Subroutine fb_misfit - calculates a misfit value 
c				from the 
c				 FB residual
c
c	Note: Calls no other routines
c	Note: this case for atmo errs too
c
c------------------------------------------------------------------------
c 
        subroutine fb_misfit(
     &          model, 
     & 		uwind, vwind, wind_mag, 
     & 		uo_i, vo_i, o_i_mag, 
     & 		ugeo, vgeo, geo_mag, 
     & 		hifc, rhoi,
     &          beta_fb, fb_ratio, misfitval )
c
c
c
        real            model(4)
     	real   		uwind, vwind, wind_mag 
     	real   		uo_i, vo_i, o_i_mag 
     	real   		ugeo, vgeo, geo_mag 
     	real   		uwind_use, vwind_use, wind_mag_use
     	real   		hifc
     	real            beta_fb, fb_ratio, misfitval 
     	real            misfit_fb
        real    	thA, thO
        real    	Capa, Copo, rhoi
        real    	Ro, Na2
        real            atmo_stressx, ocn_stressx, cor_tiltx, tot_mtm_x
        real            atmo_stressy, ocn_stressy, cor_tilty, tot_mtm_y

        parameter       ( pi = 3.141259E0 )
c 
        integer        conv, it
c
c						misfit between observed and
c						predicted
c model contains:
c 1: Na2Ua - scaled by wind_mag to get the correct Ua
c 2: Rossby radius, Na independant - scaled by hifc
c 3: ThetaO
c 4: ThetaA - ThetaO

c	write(*,*)' model = ',model

c note scaling for model(1,2)
c when removing Ua from Na2Ua, Na2 = Na2Ua/(wind_mag*wind_mag_scale)
c When using Na2Ua, Na2Ua = Na2Ua/wind_mag_scale
c (wind_mag_scale = wind_mag/wind_mag_* )
c When using Ro, Ro = Ro/hifc_scale
c (hifc_scale = hifc/hifc_*)

	thA = model(3) 
	thO = model(3) + model(4)
	Copo= rhoi/model(2) 
	Capa= model(1)*Copo
c	write(*,*)' Na2 = ',Na2
c	write(*,*)' Ro  = ',Ro

c extra wind errors too
	
	uwind_use = uwind 
	vwind_use = vwind

	wind_mag_use = HYPOT(uwind_use,vwind_use)

c components indiidually
c atmo stress component

        atmo_stressx = Capa*wind_mag_use*( uwind_use*COS(thA) - vwind_use*SIN(thA) )
        atmo_stressy = Capa*wind_mag_use*( vwind_use*COS(thA) + uwind_use*SIN(thA) )

c ocean stress component

        ocn_stressx  = Copo*o_i_mag *( uo_i*COS(thO) - vo_i*SIN(thO) )
        ocn_stressy  = Copo*o_i_mag *( vo_i*COS(thO) + uo_i*SIN(thO) )

c coriolis acceleration + ocean tilt (geostrophic) component

        cor_tiltx = -vgeo*rhoi*hifc
        cor_tilty =  ugeo*rhoi*hifc

c add it all up !!

        tot_mtm_x = atmo_stressx + ocn_stressx + cor_tiltx
        tot_mtm_y = atmo_stressy + ocn_stressy + cor_tilty

c simple misfit  - just the total momentum, that should be 0

        misfit_fb = tot_mtm_x**2.0E0 + tot_mtm_y**2.0E0


 	misfitval = (misfit_fb * beta_fb )/ fb_ratio**2
c	misfitval = (misfit_fb * beta_fb )

 
        return
        end
