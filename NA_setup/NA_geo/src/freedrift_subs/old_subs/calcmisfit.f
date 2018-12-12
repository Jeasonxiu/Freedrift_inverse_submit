c-------------------------------------------------------------------------
c
c	Subroutine calcmisfit - calculates a misfit value between 
c				observed data and predicted data
c				plus model roughness
c
c	Note: Calls no other routines
c
c------------------------------------------------------------------------
c 
        subroutine calcmisfit(
     &          predicted_data, observed_data,sd, ndata, 
     &          weight, nw, misfitval )
c
c
        include 'freedrift_param.inc'
c
c
        real*4          predicted_data(maxdata,maxwave),
     &                  observed_data(maxdata,maxwave),
     &                  sd(maxdata, maxwave),
     &                  weight(maxwave)
	real*4          misfitval
c 
        integer         nw, ndata(maxwave)
c
c						misfit between observed and
c						predicted
        misfitval=0.d0

	
        do iw=1,nw
          do i=1,ndata(iw)

            aval=(observed_data(i,iw)-predicted_data(i,iw))**2 / sd(i,iw)**2

            misfitval=misfitval+aval

          enddo
        enddo

	misfitval = misfitval * 10000.0	

 
        return
        end
