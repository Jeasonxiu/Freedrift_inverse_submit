C-------------------------------------------------------------------------
c
c	Subroutine geo_misfit - calculates a misfit value between 
c				observed data and predicted data
c				three times, calculating vels for
c				wind,  o_i geo,from the other two
c				This new version allows for full circle angles
c				Misift set within the loop, 
c				default value = 10*vector
c
c	Note: Calls no other routines
c
c------------------------------------------------------------------------
c 
        subroutine geo_misfit(
     &          model, 
     & 		uwind, vwind, wind_mag, wind_sig,
     & 		uo_i, vo_i, o_i_mag, o_i_sig,
     & 		ugeo, vgeo, geo_mag, geo_sig,
     & 		hifc, rhoi,
     & 		beta_a, beta_i, beta_g, 
     &          misfitval )
c
c
c
        real            model(4)
     	real		uwind, vwind, wind_mag, wind_sig
     	real		uo_i, vo_i, o_i_mag, o_i_sig
     	real    	ugeo, vgeo, geo_mag, geo_sig
     	real		uwind_calc, vwind_calc
     	real		uo_i_calc, vo_i_calc
     	real		ugeo_calc, vgeo_calc
     	real		rhoi, hifc
        real		beta_a, beta_i, beta_g
     	real   	        fb_ratio, misfitval 
     	real*8         misfit_o, misfit_g
     	real*8         misfit_a
        real*8    	y1, y11, y12
        real*8    	y2, y21, y22
        real*8    	y_mag, y_arg
        real*8    	thA, thO, thAO
        real*8    	Ro, Na2, bigp

        parameter       ( pi = 3.141259D0 )
c       parameter       ( rhoi = 917.0D0 )
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
c 3: thO
c 4: thA - thO

c	write(*,*)' model = ',model

c When using Ro, Ro = Ro*hifc

	thA = model(3) + model(4)
	thO = model(3)
	thAO= model(4)
	Ro  = model(2)*hifc
	Na2 = model(1)
c	write(*,*)' Na2 = ',Na2
c	write(*,*)' Ro  = ',Ro

c Start with wind mag from uo_i and ugeo (ugeo is uo_i^geo)
	
	y11 = 0.0D0
	y12 = 0.0D0
	y21 = 0.0D0
	y22 = 0.0D0
	y1 = 0.0D0
	y2 = 0.0D0

        y11 = -(o_i_mag/Na2) * (COS(thO-thA)*uo_i - SIN(thO-thA)*vo_i)
        y12 = -(o_i_mag/Na2) * (COS(thO-thA)*vo_i + SIN(thO-thA)*uo_i)
    
c        y21 = -(Ro/Na2) * (-SIN(-thA)*ugeo - COS(-thA)*vgeo)
c        y22 = -(Ro/Na2) * (-SIN(-thA)*vgeo + COS(-thA)*ugeo)
    
        y21 = (Ro/Na2) * (SIN(-thA)*ugeo + COS(-thA)*vgeo)
        y22 = (Ro/Na2) * (SIN(-thA)*vgeo - COS(-thA)*ugeo)
    
c	write(*,*)y11,y12,y21,y22

        y1 = MAX(ABS(y11 + y21),1.0E-10)*SIGN(1.0D0,y11 + y21)
c       y1 = y11 + y21
        y2 = y12 + y22
c       y_mag = sqrt(hypot(y1,y2))
        y_mag = SQRT(HYPOT(y1,y2))
        y_arg = ATAN2(y2,y1)

c	write(*,*)y1,y2,y_mag,y_arg
    
        uwind_calc = y_mag*COS(y_arg)
        vwind_calc = y_mag*SIN(y_arg)

c	write(*,*)uwind_calc,vwind_calc
    
	misfit_a = (uwind - uwind_calc)**2 + (vwind - vwind_calc)**2
c	write(*,*)' misfit_a = ',misfit_a
	misfit_a = misfit_a * beta_a
	
c	write(*,*)' misfit_a = ',misfit_a

c Now for i_o

	y11 = 0.0D0
	y12 = 0.0D0
	y21 = 0.0D0
	y22 = 0.0D0
	y1 = 0.0D0
	y2 = 0.0D0

        y11 = -(wnd_mag*Na2) * (COS(thA-thO)*uwnd - SIN(thA-thO)*vwnd)
        y12 = -(wnd_mag*Na2) * (COS(thA-thO)*vwnd + SIN(thA-thO)*uwnd)
    
        y21 = -Ro * (-SIN(-thO)*ugeo - COS(-thO)*vgeo)
        y22 = -Ro * (-SIN(-thO)*vgeo + COS(-thO)*ugeo)
    
        y1 = MAX(ABS(y11 + y21),1.0E-10)*SIGN(1.0D0,y11 + y21)
c       y1 = y11 + y21
        y2 = y12 + y22
c       y_mag = sqrt(hypot(y1,y2))
        y_mag = SQRT(HYPOT(y1,y2))
        y_arg = ATAN2(y2,y1)
    
        uo_i_calc = y_mag*COS(y_arg)
        vo_i_calc = y_mag*SIN(y_arg)
    
	misfit_o = (uo_i - uo_i_calc)**2 + (vo_i - vo_i_calc)**2
	misfit_o = misfit_o * beta_i

c	write(*,*)' misfit_o = ',misfit_o
c Finally geo

	y11 = 0.0D0
	y12 = 0.0D0
	y21 = 0.0D0
	y22 = 0.0D0
	y1 = 0.0D0
	y2 = 0.0D0

        y11 = (wnd_mag*Na2/Ro) * (-SIN(thA)*uwnd - COS(thA)*vwnd)
        y12 = (wnd_mag*Na2/Ro) * (-SIN(thA)*vwnd + COS(thA)*uwnd)
    
c       y21 = (o_i_mag/Ro) * (-SIN(thO)*uo_i - COS(thO)*vo_i)
c       y22 = (o_i_mag/Ro) * (-SIN(thO)*vo_i + COS(thO)*uo_i)
        y21 = -1.0D0*(o_i_mag/Ro) * (SIN(thO)*uo_i + COS(thO)*vo_i)
        y22 = -1.0D0*(o_i_mag/Ro) * (SIN(thO)*vo_i - COS(thO)*uo_i)
    
        y1 = MAX(ABS(y11 + y21),1.0E-10)*SIGN(1.0D0,y11 + y21)
c       y1 = y11 + y21
        y2 = y12 + y22
    
        ugeo_calc = y1
        vgeo_calc = y2
c        y_mag = hypot(y1,y2)
c        y_arg = ATAN2(y2,y1)
c    
c        ugeo_calc = y_mag*COS(y_arg)
c        vgeo_calc = y_mag*SIN(y_arg)
    
c now just misfit the magnitude
c	misfit_g = (ugeo - ugeo_calc)**2 + (vgeo - vgeo_calc)**2
	misfit_g = SQRT(ugeo**2 + vgeo**2) - SQRT(ugeo_calc**2 + vgeo_calc**2)
c	misfit_g = misfit_g * beta_g
	misfit_g = misfit_g**2 * beta_g

c	write(*,*)' misfit_g = ',misfit_g

	misfitval = misfit_o/o_i_sig**2 +
     &              misfit_a/wind_sig**2 +
     &              misfit_g/geo_sig**2 
c	misfitval = misfit_o+
c     &              misfit_a+
c     &              misfit_g

c note that the errors are dimensional now
 
        return
        end
