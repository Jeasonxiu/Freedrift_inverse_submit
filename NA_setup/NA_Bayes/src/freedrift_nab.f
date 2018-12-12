c
c ----------------------------------------------------------------------------
c
c       Driver program that calls the Neighbourhood algorithm 
c	ensemble inference routines to calculate Bayesian 
c	information measures from an input ensemble of models.
c
c
c               	MEMORY AND ARRAY SIZES
c
c                                       The NAB routines use the
c                                       include file 'nab_param.inc' 
c					to define all parameters controlling
c                                       memory required by the arrays.
c
c                                       A description of each parameter,
c                                       indicating which ones should be
c                                       changed to fit your application
c                                       can be found in 'nab_param.inc'
c
c
c               	INPUT AND OUTPUT
c
c       Input files:
c                   nab.in              Contains options for Neighbourhood
c                                       algorithm integration (see HTML manual)
c
c		    ensemble.nad	Direct access NAD file contain
c					input models and misfit values
c
c       Output files:
c                   nab.sum             summary of results
c
c                   nab.out             output file containing results of
c					all integrals and marginals
c					(This file is read by plot programs
c					 `naplot-x, naplot-p') 
c
c                   sobol.coeff         initializing data used for
c                                       quasi-random sequences
c                                       (This is output for reference only)
c
c       Comments:
c
c                Logical units 30-40 are reserved for use by NAB subroutines
c                The user specific routines should not use these values
c                for logical units. The NAB routines also write to 
c		 LU 6 as standard out.
c
c                                               M. Sambridge, RSES, ANU.
c                                               Last updated Sept. 1999.
c
c ----------------------------------------------------------------------------
c					
        Program freedrift_nab
c                                       Call NAB routine to do the work
	call nab

	stop
	end
c
c==========================================================================
c 
c		 USER SUPPLIED ROUTINES BELOW
c
c==========================================================================
c
c---------------------------------------------------------------------------
c
c     logPPD_NA - calculates log-PPD of posterior probability 
c	          density function for user supplied problem using the
c                 Neighbourhood approximation to the PPD.
c
c	       THIS IS A PROBLEM SPECIFIC USER SUPPLIED SUBROUTINE
c
c       Input:
c	      node		: index of model from NAD input file
c	      misfit(*)		: array of input data (misfit values)
c
c       Output:
c	      logPPD_NA	: log-posterior probability density
c			  function of model
c
c     Comments:
c	       This routine converts the input data array to a
c	       log-posteriori probability density function. 
c	       It is a user supplied routine.
c
c	       For example if the input data for each model is
c	       a simple sum of squares of residuals weighted
c	       by a priori data covariances (i.e. standard least squares)
c	       then the Log-PPD is just a factor of -0.5 times
c	       this. This routine allows the user to use other
c	       Likelihood functions (or forms of posteriori probability
c	       density) and if necessary rescale them, or include 
c	       a priori PDFs. 
c
c	       WARNING: All Bayesian integrals will be affected
c	       by the choice and scale of this function ! The weight
c	       of each model in the integrals depends on the 
c	       exponential of logPPD. It is very important that
c	       the data covariances are chosen suitably for the user problem. 
c
c
c                                       M. Sambridge, RSES, July 1998.
c
c---------------------------------------------------------------------------
c
        Function logPPD_NA(node,misfit)

	real*4		misfit(*)
	real*4		logPPD_NA

 	data		degofreedom/16/

c	logPPD = 1.0
c
c			Note that in this case the input misfit is assumed
c			to be the chi-square of the model fit and 
c			has a factor of 1/degofreedom in it which is cancelled 
c			out here.
c
c			In this example there are 116 degrees of freedom
c			(see Sambridge 1999b for details)
c
c			In this case 10 momentum balances * 2 data each - 4 param
c			16 degoff

 	logPPD_NA = -0.5*degofreedom*misfit(node)

	return
	end
c
c---------------------------------------------------------------------------
c
c     logPPD - calculates log of posterior probability 
c	       density function for user supplied problem
c	       using full forward modelling.
c
c	       THIS IS A PROBLEM SPECIFIC USER SUPPLIED SUBROUTINE
c
c       Input:
c	      x		: model at which logPPD is to be evaluated
c
c       Output:
c	      logPPD	: log-posterior probability density
c			  function of input model
c
c     Comments:
c              This routine was written to allow the program to 
c              be used for evaluation of Bayesian integrals where the
c              forward problem IS to be called. i.e. here we do not use the
c              neighbourhood approximation to the forward problem
c              but instead solve the actual forward problem
c
c              Definition of the logPPD:
c
c	       An example: If the input data for each model is
c	       a simple sum of squares of residuals weighted
c	       by a priori data covariances (i.e. standard least squares)
c	       then the Log-PPD is just a factor of -0.5 times
c	       this. This routine allows the user to use other
c	       Likelihood functions (or forms of posteriori probability
c	       density) and if necessary rescale them, or include 
c	       a priori PDFs. 
c
c	       WARNING: All Bayesian integrals will be affected
c	       by the choice and scale of this function ! The weight
c	       of each model in the integrals depends on the 
c	       exponential of logPPD. It is very important that
c	       the data covariances are chosen suitably for the user problem. 
c
c                                       M. Sambridge, RSES, July 1998.
c
c---------------------------------------------------------------------------
c
        Function logPPD(x)

	real*4		x(*)
	real*4		chi
	real*4		logPPD
c
	data		degoffreedom/16/

c			Note that in this case the input misfit is assumed
c			to be the chi-square of the model fit and 
c			has a factor of 1/degofreedom in it which is cancelled 
c			out here.
c
c			In this example there are 116 degrees of freedom
c			(see Sambridge 1999b for details)
c
	nd = 24
c			calculate synthetics for this model and compare to
c			data to calculate chi-square measure.
c
	call forward(nd,x,chi)

        logPPD = -0.5*degoffreedom*chi

	return
	end
c---------------------------------------------------------------------------
c
c     user_init - user supplied routine for any initialization 
c	          tasks needed for user supplied problem.
c
c		  In this case receiver receiver function inversion
c		  initilization is performed. Sac files and inverse
c		  data covariance matrix are read in
c		  for the observed receiver function and stored in
c		  a common block so that routine forward can access them.
c
c	       THIS IS A PROBLEM SPECIFIC USER SUPPLIED SUBROUTINE
c
c	Input: - none
c
c	Output: - none
c
c     Comments:
c	       This routine is called once at the start of the
c	       nab program to perform any initialization tasks
c	       that may be needed. For example, calculate a priori
c	       probability density functions, read in data etc.
c	       If no tasks are required then a simple dummy routine
c
c	       Receiver function setup:
c
c       	 This code is based on work by T. Shibutani 
c		 (RCEP, DPRI, KYOTO UNIV.) who wrote the receiver 
c		  function routines.  
c
c				M. Sambridge, RSES (ANU), April 2005. 
c
c-------------------------------------------------------------------------
c
	subroutine user_init()

c						initialize receiver
c						function forward modelling
c
	return
	end
c
c-------------------------------------------------------------------------
c
c	forward - performs forward modelling for user supplied problem.
c		  In this case it calculates predicted receiver function
c		  for a single model and calculates the misfit measure
c		  between observation and prediction using a data
c		  covariance matrix passed to it in a common block.
c
c	Input: 
c	      nd		:Number of dimensions in parameter space
c	      model(nd)		:input velocity model
c
c	Output:
c	      logppd		:negative log(ppd) 
c
c	Comments:
c		 This routine allows the user to perform the forward 
c		 modelling and define an a posterior probability density
c		 function using the resulting mismatch with the observed
c		 data. 
c
c		 The forward routine used here differs from that 
c		 in supplied to the NA search routine. Here it
c		 calculates the same synthetics but reads in a 
c		 full non-diagonal data covariance matrix from file
c	         `icov.in' and uses this to calculate
c		 the chi-square misfit statistic, or the -log PPD
c		 for this problem.  
c
c		 Note that in this example the chi-square statistic 
c	 	 calculated here is transformed into the actual logppd 
c		 by the user supplied routine `logPPD' called in the
c		 main NAB program.
c
c				M. Sambridge, RSES (ANU), Jan. 1999. 
c
c-------------------------------------------------------------------------
c
	subroutine forward(nd,model,chi)

c						initialize receiver
c						function forward modelling

	real*4		chi
c	real*4		misfitval
	real*4		model(nd)
c
c	integer		ndata(maxwave)

	return
	end
c
c---------------------------------------------------------------------------
c
c     evaluate_sp - evaluates the i-th user supplied special function of 
c		    the model for numerical intergation.
c
c	       THIS IS A PROBLEM SPECIFIC USER SUPPLIED SUBROUTINE
c
c	Input:
c	      i			:index of special function to evaluate
c	      x(nd)		:model values in raw  (input) units
c	      nd		:number of dimensions
c	      node		:node of Voronoi cell containing x.
c
c	Output:
c	      spf		:value of i-th special function
c
c       Comments:
c		 This routine is for the receiver function 
c		 parameterisation used by Shibutani et al. (1997).
c
c		 Parameter 1  = velocity jump across Moho
c		 Parameter 2  = depth to Moho
c		 Parameter 3  = velocity gradient in first layer
c		 Parameter 4  = velocity gradient in second layer
c		 Parameter 5  = velocity gradient in third layer
c		 Parameter 6  = velocity gradient in fourth layer
c		 Parameter 7  = velocity gradient in fifth layer
c		 Parameter 8  = depth of first interface 
c		 Parameter 9  = depth of second interface 
c		 Parameter 10 = depth of third interface 
c		 Parameter 11 = depth of fourth interface 
c		 Parameter 12 = depth of fifth interface 
c
c                                               M. Sambridge, March 1998.
c
c---------------------------------------------------------------------------
c
        Subroutine evaluate_sp (i,x,nd,node,range,scales,spf)

        real*4		range(2,*)
        integer		scales(*)
c
	real*4		x(nd)

	if(i.eq.1)then
c						evaluate DVmoho 
	   spf = x(12)-x(17)

	else if(i.eq.2)then
c						evaluate depth to Moho 

	   spf = x(1)+x(2)+x(3)+x(4)+x(5)

	else if(i.eq.3)then
	   t = x(1)
	   if(t.eq.0.0)t = 0.1
	   spf = (x(13)-x(7))/t
	else if(i.eq.4)then
	   t = x(2)
	   if(t.eq.0.0)t = 0.1
	   spf = (x(14)-x(8))/t
	else if(i.eq.5)then
	   t = x(3)
	   if(t.eq.0.0)t = 0.1
	   spf = (x(15)-x(9))/t
	else if(i.eq.6)then
	   t = x(4)
	   if(t.eq.0.0)t = 0.1
	   spf = (x(16)-x(10))/t
	else if(i.eq.7)then
	   t = x(5)
	   if(t.eq.0.0)t = 0.1
	   spf = (x(17)-x(11))/t
	else if(i.eq.8)then
	   spf = x(1)
	else if(i.eq.9)then
	   spf = x(1)+x(2)
	else if(i.eq.10)then
	   spf = x(1)+x(2)+x(3)
	else if(i.eq.11)then
	   spf = x(1)+x(2)+x(3)+x(4)
	else if(i.eq.12)then
	   spf = x(1)+x(2)+x(3)+x(4)+x(5)
	else
	end if

 100    format(/'  *****  Error in subroutine evaluate_sp  *****'/
     &          '  Number of special functions supplied =',
     &          i4,/'  Number requested =',i4//
     &          '  Check number of special functions in input file'/)

	return
	end
c
c---------------------------------------------------------------------------
c
c     evaluate_sp_range - returns the range of the i-th user supplied
c			  special function 
c
c	       THIS IS A PROBLEM SPECIFIC USER SUPPLIED SUBROUTINE
c
c	Input:
c	      i			:index of special function to evaluate
c             range             :ranges of all parameter values.
c             scales            :scale factors for each parameter.
c
c	Output:
c	      spmax		:maximum value of i-th special function
c	      spmin		:minimum value of i-th special function
c
c       Comments:
c		 This routine is for the receiver function 
c		 parameterisation used by Shibutani et al. (1997).
c
c		 Parameter 1 = velocity jump across Moho
c		 Parameter 2 = depth to Moho
c		 Parameter 3  = velocity gradient in first layer
c		 Parameter 4  = velocity gradient in second layer
c		 Parameter 5  = velocity gradient in third layer
c		 Parameter 6  = velocity gradient in fourth layer
c		 Parameter 7  = velocity gradient in fifth layer
c		 Parameter 8  = depth of first interface 
c		 Parameter 9  = depth of second interface 
c		 Parameter 10 = depth of third interface 
c		 Parameter 11 = depth of fourth interface 
c		 Parameter 12 = depth of fifth interface 
c
c                                               M. Sambridge, March 1998.
c
c---------------------------------------------------------------------------
c
        Subroutine evaluate_sp_range(i,range,scales,spmax,spmin)
c
	real*4		range(2,*)
        integer		scales(*)

	if(i.eq.1)then
c						evaluate DVmoho range
	   spmin = range(1,12)-range(2,17)
	   spmax = range(2,12)-range(1,17)

	else if(i.eq.2)then
c						evaluate Moho range

	   spmin = range(1,1)+range(1,2)+range(1,3)
     &             +range(1,4)+range(1,5)
	   spmax = range(2,1)+range(2,2)+range(2,3)
     &             +range(2,4)+range(2,5)

	else if(i.eq.3)then
	   t = range(1,1)
	   if(t.eq.0.0)t = 0.1
	   spmax = (range(2,7)-range(1,13))/t
	   if(range(1,13).lt.range(2,7))then
	      spmin = (range(1,13)-range(2,7))/t
	   else
	      spmin = (range(1,13)-range(2,7))/range(2,1)
	   end if
	else if(i.eq.4)then
	   t = range(1,2)
	   if(t.eq.0.0)t = 0.1
	   spmax = (range(2,14)-range(1,8))/t
	   if(range(1,14).lt.range(2,8))then
	      spmin = (range(1,14)-range(2,8))/t
	   else
	      spmin = (range(1,14)-range(2,8))/range(2,2)
	   end if
	else if(i.eq.5)then
	   spmax = (range(2,15)-range(1,9))/range(1,3)
	   if(range(1,15).lt.range(2,9))then
	      t = range(1,3)
	      if(t.eq.0.0)t = 0.1
	      spmin = (range(1,15)-range(2,9))/t
	   else
	      spmin = (range(1,15)-range(2,9))/range(2,3)
	   end if
	else if(i.eq.6)then
	   spmax = (range(2,16)-range(1,10))/range(1,4)
	   if(range(1,16).lt.range(2,10))then
	      t = range(1,4)
	      if(t.eq.0.0)t = 0.1
	      spmin = (range(1,16)-range(2,10))/t
	   else
	      spmin = (range(1,16)-range(2,10))/range(2,4)
	   end if
	else if(i.eq.7)then
	   spmax = (range(2,17)-range(1,11))/range(1,5)
	   if(range(1,17).lt.range(2,11))then
	      t = range(1,5)
	      if(t.eq.0.0)t = 0.1
	      spmin = (range(1,17)-range(2,11))/t
	   else
	      spmin = (range(1,17)-range(2,11))/range(2,5)
	   end if
	else if(i.eq.8)then
	   spmin = range(1,1)
	   spmax = range(2,1)
	else if(i.eq.9)then
	   spmin = range(1,1)+range(1,2)
	   spmax = range(2,1)+range(2,2)
	else if(i.eq.10)then
	   spmin = range(1,1)+range(1,2)+range(1,3)
	   spmax = range(2,1)+range(2,2)+range(2,3)
	else if(i.eq.11)then
	   spmin = range(1,1)+range(1,2)+range(1,3)
     &             +range(1,4)
	   spmax = range(2,1)+range(2,2)+range(2,3)
     &             +range(2,4)
	else if(i.eq.12)then
	   spmin = range(1,1)+range(1,2)+range(1,3)
     &             +range(1,4)+range(1,5)
	   spmax = range(2,1)+range(2,2)+range(2,3)
     &             +range(2,4)+range(2,5)

	end if

	return
	end
