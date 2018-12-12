c------------------------------------------------------------------------
c
c	Subroutine forward_modelling - calculates the array "predicted_
c				       data" from the model "rmodel"
c
c	Note_1: Calls are made to : theo
c
c       Note_2: rmodel(1:nlayer): thickness
c               rmodel(nlayer+1:2*nlayer): velocity at upper interface
c               rmodel(2*nlayer+1:3*nlayer): velocity at lower interface
c               rmodel(3*nlayer+1:4*nlayer): Vp/Vs ratio
c
c-------------------------------------------------------------------------
c 
        subroutine forward_modelling(
     &          rmodel, moddim, time_shift, ndata, fs,
     &          predicted_data,
     &          n_iter, iproc, model_path, moho,
     &          station, net, code )

c 
        include 'freedrift_param.inc'
 
 
        real*8          rmodel(maxmoddim)
c 
        real*4          incident_angle(maxwave),
     &                  constant_a(maxwave),
     &                  constant_c(maxwave),
     &                  time_shift(maxwave),
     &                  predicted_data(maxdata,maxwave),
     &			T0_synt(maxdata),
     &                  q_alpha(maxlayer),
     &                  q_beta(maxlayer)
 
        integer         ndata(maxwave),
     &                  nn(maxdata), ll(maxdata),
     &                  micio(maxdata)
c my variables
        integer         iproc, ipr1, ipr2, n_iter
	real		moho
        character*10	moho_string, n_iter_string
	character*200	model_path
	character*7     node_name, clat, clon, lat, lon
        character*1     typeo, q
        character*3     mot
        character*100   str, model, input, mdir, hdir, prem, cdir
 	character*100	path
	character*400	command
	character*100	outmodel
        character*400   ifile, ofile, cmd
        integer         ifanis, ifdeck, nm, nic, noc, im
        parameter       ( maxmod=450 )
        parameter       ( nmin = 0)
        parameter       ( lmin = 0)
        parameter       ( lmax = 300 )
        parameter       ( fmin = 0.0)
        parameter       ( fmax = 50.0 )
        integer*8       norder, lorder, eigid, nrow
        integer*4       ncol, mmax
        real*8          tref, phvel, grvel, attn
        real*8          dmi1, dmi2, mrho, mvpv, mvsv
        real*8          ray, U, dU, V, dV, P, dP, A2
        real*8          A(lmax),T(lmax),y2(lmax),nA(lmax),nT(lmax),
     &                  A0(maxdata),T0(maxdata),nA0(maxdata),nT0(maxdata)
 

        integer*4       icru, imoh, nc, ic, nmod
        real*8          crad(70), crho(70), cvpv(70), cvsv(70), cerr(70),
     &                  cqmu(70), cqka(70), cvph(70), cvsh(70), ceta(70)
        real*8          llat, llon, ldep
        character*15    clay
        real*8          t1(9), t2(9), t3(9), t4(9), x, y, z
	character*6	station, code, net
	
        character*2     str1
        character*20    str2
        character*40    model_name
        character*10    str3
        character*20    str4
        character*120   str5
        character*200   str6
        character*10    str7
	character*100	apriori_model_path

        real*4          ap_H(30),
     &                  ap_VP(30),
     &                  ap_VS(30),
     &                  ap_RHO(30),
     &			e,f,g,h, ii, l,
     &			thick(maxlayer),
     &                  vs(maxlayer),
     &                  vp(maxlayer),
     &                  rho(maxlayer)


	
	lw_sta = lofw(station)
	lw_net = lofw(net)
	lw_code = lofw(code)

	
	mdir='/scratch/scratch/zcahbe4/NA_tutorial/Inversion/NA_'// code(1:lw_code)// '/models'
	

c	create folder 
        call hostnm(node_name)
        write(model, '(a,".",i2.2)') trim(node_name), iproc
        write(path,'(a,a,a)') trim(mdir),"/model_", trim(model)
	

	write(command, '(a,a)') "rm -r ", path
	call system(command)

	write(command,'(a,a)') "mkdir ", path
	call system(command)	

	write(model_path, '(a,a)') trim(path),"/modello.d"

c-----	load apriori model

	write(apriori_model_path,'(a,a,a)' ) "/home/zcahbe4/NA_tutorial/CRUST1.0/", station(1:lw_sta), "_CRUST1.0.d"
	open(111, file=apriori_model_path, status="old", access="sequential")
        do i=1,12
                read(111, *)
        end do
        j = 1
        do i=12, 100
                read(111,*,end=99) ap_H(j), ap_VP(j), ap_VS(j), ap_RHO(j), e,f,g,h,ii,l
		j = j + 1
        end do
  99    close(111)


c------ define parameters

        nlayer=moddim/4
        k=0
        qa=1450.0
        qb=600.0
        etap=1.0
        etas=1.0
        frefp=1.0
        frefs=1.0


        thick(1) = 3.0 
        thick(2) = 8.0
        thick(3) = (moho - 11.0) / 2
        thick(4) = thick(3)

	
        do i=1,nlayer
                vs(i) = rmodel(nlayer+i)
        end do

c-------------------------------------------------------------------------------
        open(2,file=model_path, form='formatted', access='sequential', status='unknown')
        write(2,fmt="(a)")'MODEL.01'
        write(2,fmt="(a)")'MODELLO DEMO'
        write(2,fmt="(a)")'ISOTROPIC'
        write(2,fmt="(a)")'KGS'
        write(2,fmt="(a)")'SPHERICAL EARTH'
        write(2,fmt="(a)")"1-D"
        write(2,fmt="(a)")'CONSTANT VELOCITY'
        write(2,fmt="(a)")'LINE08'
        write(2,fmt="(a)")'LINE09'
        write(2,fmt="(a)")'LINE10'
        write(2,fmt="(a)")'LINE11'
        write(2,fmt="(a)")'H    VP      VS      RHO     QP      QS      ETAP    ETAS    FREFP   FREFS'


	do i=1, nlayer
        	call brocher_Vp_Vs(vs(i), vvp)
                vp_brocher = vvp
                call brocher_rho_Vs(vs(i), rr)
                rho_brocher=rr

                write(2,fmt="(f0.2,a,f0.4,a,f0.4,a,f0.4,a,f7.1,a,f5.1,a,f3.1,a,f3.1,a,f0.1,a,f0.1)")thick(i), ' ',
     &           vp_brocher,'   ',vs(i),'       ',rho_brocher,' ',qa,'  ',qb,'  ',etap,'        ',etas,'        ',frefp,'       ',frefs

                crustal_thick = crustal_thick + thick(i)

	end do

c----- 	Complete the model with the apriori model
        do i=5,measurement_number
                write(2,fmt="(f0.2,a,f0.4,a,f0.4,a,f0.4,a,f7.1,a,f5.1,a,f3.1,a,f3.1,a,f0.1,a,f0.1)") 
     &		ap_H(i), ' ', ap_VP(i),'       ',ap_VS(i),'    ',ap_RHO(i),'   ',qa,'  ',qb,'  ',etap,
     &          '        ',etas,'        ',frefp,'       ',frefs

        end do
        close(2)


c#-------------------------------------------------------------------------------------------------------------------


        str5 = 'python /home/zcahbe4/NA_tutorial/NA/src/rfi_subs/forward_problem.py  '
        command = str5//path
	
        call system(command)

        open(777, file=trim(path)//'/predicted_data.txt', status='unknown')
        do i=1, 13
                read(777,*) T0_synt(i), predicted_data(i,1)
	end do
        close(777)



	return
	end





c=========================================================================================================
	subroutine brocher_rho_Vs(Vs, rho)
                call brocher_Vp_Vs(Vs, Vp)
                call brocher_rho_Vp(Vp, rho)
        return
        end
c---------------------------
        subroutine brocher_rho_Vp(Vp, rho)
                rho = 1.6612*Vp - 0.4721*(Vp**2) + 0.0671*(Vp**3) - 0.0043*(Vp**4) + 0.000106*(Vp**5)
        return
        end
c---------------------------
	subroutine brocher_Vp_Vs(Vs, Vp)
                Vp = 0.9409 + 2.0947*Vs - 0.8206*(Vs**2) + 0.2683*(Vs**3) - 0.0251*(Vs**4)
        return
        end

