!****************************************************************
!**                                                            **
!**      cgod  for                                           **
!**                                                            **
!****************************************************************
	
module cgod

	contains
	
      !@cuf attributes(host, device) & 
subroutine cgod3d(KOBL,i,j,k,kdir,idgod, &
                                 al,be,ge,el, &
                                 w,qqq1,qqq2, &
                                 dsl,dsp,dsc,ythll,ga, &
                                 qqq, kontact_)
!     include  'cyl.fl'
      implicit real*8 (a-h,o-z)

      dimension qqq(8),qqq1(8),qqq2(8)
      LOGICAL, intent(in), OPTIONAL :: kontact_
      LOGICAL :: kontact

      data x0,x1,x2,x3,x4,x5,x6,x7,x8,x9/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./

      if(.not. present(kontact_)) then
            kontact = .False. 
      else 
            kontact = kontact_ 
	end if

      x10=10.d0

        idgod=0
        kvak=0
!       if(kdir.eq.1) krsp=irsp
!       if(kdir.eq.2) krsp=jrsp

        eps=x10**(-x10)
!       eps=x10**(-20)
       eps8=x10**(-8)
         g1=ga-1.d0
         g2=ga+1.d0
        gg1=ga-1.d0 !g1
        gga=ga
        gg2=ga+1.d0 !g2

                  el2 = dsqrt(al**2+be**2)
                  el3 = dsqrt(al**2+ge**2)
                  if(el2.ne.x0)then
                  al2 = -be/el2
                  be2 =  al/el2
                  ge2 =  x0
                  al3 = -be2*ge
                  be3 =  al2*ge
                  ge3 =  el2
                  elseif(el3.ne.x0)then
                  al2 = -ge/el3
                  be2 =  x0
                  ge2 =  al/el3
                  al3 = -ge2*be
                  be3 =  el3
                  ge3 =  al2*be
                  else
                  if(mode.gt.0)print*,'cgod: el,el2,el3:',el,el2,el3
                  stop
                  endif

                 enI  = al * qqq1(2) + be * qqq1(3) + ge * qqq1(4)
                 teI2 = al2* qqq1(2) + be2* qqq1(3) + ge2* qqq1(4)
                 teI3 = al3* qqq1(2) + be3* qqq1(3) + ge3* qqq1(4)
                 enII = al * qqq2(2) + be * qqq2(3) + ge * qqq2(4)
                 teII2= al2* qqq2(2) + be2* qqq2(3) + ge2* qqq2(4)
                 teII3= al3* qqq2(2) + be3* qqq2(3) + ge3* qqq2(4)

                 pI   = qqq1(5)
                 pII  = qqq2(5)
                 rI   = qqq1(1)
                 rII  = qqq2(1)



        ipiz = 0
        if(pI.gt.pII) then

           eno2 = enII
           teo22= teII2
           teo23= teII3
           p2   = pII
           r2   = rII

           eno1 = enI
           teo12= teI2
           teo13= teI3
           p1   = pI
           r1   = rI

           enI  =-eno2
           teI2 = teo22
           teI3 = teo23
           pI   = p2
           rI   = r2

           enII =-eno1
           teII2= teo12
           teII3= teo13
           pII  = p1
           rII  = r1
           w    = -w
        ipiz = 1
        endif

         gp=g2/ga/x2
         gm=g1/ga/x2
         cI =x0
         cII =x0
         if(pI.ne.x0) cI =dsqrt(ga*pI /rI )
         if(pII.ne.x0) cII =dsqrt(ga*pII /rII )
         
        a   = dsqrt( rI * (  g2 * pII + g1 * pI ) / x2 )
        Uud = ( pII - pI ) / a
        Urz = - x2 * cII / g1 * ( x1 - ( pI / pII )**gm )
        Uvk = - x2 * ( cII + cI ) / g1
        Udf = enI - enII


        
        if (Udf.lt.Uvk) then
              il=-1
              ip=-1
        endif

        if (Udf.ge.Uvk.and.Udf.le.Urz) then
             p = pI * ( ( Udf - Uvk ) / ( Urz - Uvk ) )**(x1/gm)
             il = 0
             ip = 0
        endif

        if (Udf.gt.Urz.and.Udf.le.Uud) then
             call  devtwo(enI ,pI ,rI,&
                         enII,pII,rII,i,j,k,KOBL,kdir,idgod,ga,&
                         w   ,p)
             if(idgod.eq.2)return
             il=1
             ip=0
         endif

        if (Udf.gt.Uud) then
             call  newton(enI ,pI ,rI, & 
                         enII,pII,rII,i,j,k,KOBL,kdir,idgod,ga, &
                         w   ,p)
             if(idgod.eq.2)return
             il=1
             ip=1
        endif


!*********TWO SHOCKS**********************************************
         if(il.eq.1.and.ip.eq.1) then

         aI =dsqrt(rI *( g2/x2*p+g1/x2*pI  ))
         aII=dsqrt(rII*( g2/x2*p+g1/x2*pII ))

         u=(aI*enI+aII*enII+pI-pII)/(aI+aII)
         dI =enI -aI /rI 
         dII=enII+aII/rII
         dsl=dI
         dsp=dII
         dsc = u
               if(w.le.dI ) then
                           en=enI
                           p =pI
                           r =rI
                          te2=teI2
                          te3=teI3
               endif                           
               if(w.gt.dI.and.w.le.u) then
                           en=u
                           p =p
                           r =rI*aI/(aI-rI*(enI-u))
                          te2=teI2
                          te3=teI3
               endif                           
               if(w.gt.u.and.w.lt.dII) then
                           en=u
                           p =p
                           r =rII*aII/(aII+rII*(enII-u))
                          te2=teII2
                          te3=teII3
               endif                           
               if(w.ge.dII) then
                           en=enII
                           p =pII
                           r =rII
                          te2=teII2
                          te3=teII3
               endif                           
          endif


!*********LEFT - SHOCK, RIGHT - EXPANSION FAN*******************
         if(il.eq.1.and.ip.eq.0) then

         aI=dsqrt(rI*( g2/x2*p+g1/x2*pI ))

         if(dabs(p-pII).lt.eps) then
         aII =rII*cII
         else
         aII=gm*rII*cII*(x1-p/pII)/(x1-(p/pII)**gm)
         endif

         u=(aI*enI+aII*enII+pI-pII)/(aI+aII)
         dI  =enI-aI/rI
         dII =enII +cII 
         ddII=u   +cII -g1*(enII-u)/x2
         dsl=dI
         dsp=dII
         dsc = u

               if(w.le.dI ) then
                           en=enI
                           p =pI
                           r =rI
                          te2=teI2
                          te3=teI3
               endif                           
               if(w.gt.dI.and.w.le.u) then
                           en=u
                           p =p
                           r =rI*aI/(aI-rI*(enI-u))
                          te2=teI2
                          te3=teI3
               endif                           
               if(w.gt.u.and.w.le.ddII) then
                           ce=cII-g1/x2*(enII-u)
                           en=u
                           p =p
                           r =ga*p/ce/ce
                          te2=teII2
                          te3=teII3
               endif                           
               if(w.gt.ddII.and.w.lt.dII) then
                           ce=-g1/g2*(enII-w)+x2/g2*cII
                           en=w-ce
                           p =pII*(ce/cII)**(x1/gm)
                           r =ga*p/ce/ce
                          te2=teII2
                          te3=teII3
               endif                           
               if(w.ge.dII) then
                           en=enII
                           p =pII
                           r =rII
                          te2=teII2
                          te3=teII3
               endif                           
          endif

!*********TWO EXPANSION FANS**************************************
         if(il.eq.0.and.ip.eq.0) then

         if(dabs(p-pI).lt.eps) then
         aI =rI*cI
         else
         aI=gm*rI*cI*(x1-p/pI)/(x1-(p/pI)**gm)
         endif

         if(dabs(p-pII).lt.eps) then
         aII =rII*cII
         else
         aII=gm*rII*cII*(x1-p/pII)/(x1-(p/pII)**gm)
         endif

         u=(aI*enI+aII*enII+pI-pII)/(aI+aII)
         dI  =enI  -cI  
         ddI =u   -cI  -g1*(enI -u)/x2
         dII =enII +cII 
         ddII=u   +cII -g1*(enII-u)/x2
         dsl=dI
         dsp=dII
         dsc = u


               if(w.le.dI ) then
                           en=enI
                           p =pI
                           r =rI
                          te2=teI2
                          te3=teI3
               endif                           
               if(w.gt.dI.and.w.lt.ddI) then
                           ce=g1/g2*(enI-w)+x2/g2*cI
                           en=w+ce
                           p =pI*(ce/cI)**(x1/gm)
                           r =ga*p/ce/ce
                          te2=teI2
                          te3=teI3
               endif                           
               if(w.ge.ddI.and.w.le.u) then
                           ce=cI+g1/x2*(enI-u)
                           en=u
                           p =p
                           r =ga*p/ce/ce
                          te2=teI2
                          te3=teI3
               endif                           
               if(w.gt.u.and.w.le.ddII) then
                           ce=cII-g1/x2*(enII-u)
                           en=u
                           p =p
                           r =ga*p/ce/ce
                          te2=teII2
                          te3=teII3
               endif                           
               if(w.gt.ddII.and.w.lt.dII) then
                           ce=-g1/g2*(enII-w)+x2/g2*cII
                           en=w-ce
                           p =pII*(ce/cII)**(x1/gm)
                           r =ga*p/ce/ce
                          te2=teII2
                          te3=teII3
               endif                           
               if(w.ge.dII) then
                           en=enII
                           p =pII
                           r =rII
                          te2=teII2
                          te3=teII3
               endif                           
          endif

!*********VAKUUM ************************************************
         if(il.eq.-1.and.ip.eq.-1) then
         dI  =enI  -cI  
         ddI =enI  +x2/gg1*cI  
         dII =enII +cII 
         ddII=enII -x2/gg1*cII 

         dsl=dI
         dsp=dII
         dsc = (dI+dII)/x2


               if(w.le.dI ) then
                           en=enI
                           p =pI
                           r =rI
                          te2=teI2
                          te3=teI3
               endif                           
               if(w.gt.dI.and.w.lt.ddI) then
                           ce=gg1/gg2*(enI-w)+x2/gg2*cI
                           en=w+ce
                           p =pI*(ce/cI)**(x1/gm)
                           r =gga*p/ce/ce
                          te2=teI2
                          te3=teI3
               endif                           
               if(w.ge.ddI.and.w.le.ddII) then
                           en=w
                           p =x0
                           r =x0
                          te2=x0
                          te3=x0
!             if(mode.ne.0) print *,'nastupil vakuum ',i,j,KOBL,kdir
               endif                           
               if(w.gt.ddII.and.w.lt.dII) then
                           ce=-gg1/gg2*(enII-w)+x2/gg2*cII
                           en=w-ce
                           p =pII*(ce/cII)**(x1/gm)
                           r =gga*p/ce/ce
                          te2=teII2
                          te3=teII3
               endif                           
               if(w.ge.dII) then
                           en=enII
                           p =pII
                           r =rII
                          te2=teII2
                          te3=teII3
               endif                           
          endif


        if(ipiz.eq.1) en = - en
        if(ipiz.eq.1) then
             dsl1 = dsl
             dsp1 = dsp
             dsl =   -dsp1
             dsp =   -dsl1
             dsc = -dsc
             w = -w
        endif
                       uo = al * en  + al2 * te2 + al3 * te3
                       vo = be * en  + be2 * te2 + be3 * te3
                       wo = ge * en  + ge2 * te2 + ge3 * te3

        if(KOBL.eq.3)then
          if(katom.eq.1.and.i0_bound.eq.3.and.kdir.eq.1)then
            if(j.le.n0)then
              if(i.eq.mf)then
!                              r = x0
              endif
            endif
          endif
        endif

                       eo = p / g1 + .5 * r *(uo**2+vo**2+wo**2)
                       en=al*uo+be*vo+ge*wo
                  if(kontact == .False.) then
                       qqq(1)=  (r*(en-w))
                       qqq(2)=  (r*(en-w)*uo+al*p)
                       qqq(3)=  (r*(en-w)*vo+be*p)
                       qqq(4)=  (r*(en-w)*wo+ge*p)
                       qqq(5)=  (  (en-w)*eo+en*p)
                       qqq(6)=x0
                       qqq(7)=x0
                       qqq(8)=x0
                  else
                       qqq(1)=  0.0_8
                       qqq(2)=  (al*p)
                       qqq(3)=  (be*p)
                       qqq(4)=  (ge*p)
                       qqq(5)=  (en*p)
                       qqq(6)=x0
                       qqq(7)=x0
                       qqq(8)=x0
                  end if


         return
         end



!****************************************************************
!**                                                            **
!**      newton  for                                           **
!**                                                            **
!****************************************************************
!@cuf attributes(host, device) & 
        subroutine newton(enI ,pI ,rI, &
                         enII,pII,rII,i,j,k,KOBL,kdir,idgod,ga, &
                         w   ,p)
!c       include 'cyl.fl'
      implicit real*8 (a-h,o-z)

      data x0,x1,x2,x3,x4,x5,x6,x7,x8,x9/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./
      x10=10.d0

!        ga=ga
         g1=ga-x1

         eps=x10**(-x10)
!        eps=1.D-12

         gp=(ga+x1)/ga/x2
         gm=(ga-x1)/ga/x2

         cI =dsqrt(ga*pI /rI )
         cII=dsqrt(ga*pII/rII)
         
         pn=pI*rII*cII+pII*rI*cI+(enI-enII)*rI*cI*rII*cII
         pn=pn/(rI*cI+rII*cII)
!        pn=(pI+pII)/x2

         pee = pn

          kiter=0
   1      p=pn
          if(p.le.x0) then
!            if(mode.gt.0) print *, 'negative pressure, newton',i
             stop
          endif
          kiter=kiter+1

         
         fI=(p-pI)/(rI*cI*sqrt(gp*p/pI+gm))
         fIs=(ga+x1)*p/pI+(x3*ga-x1)
         fIs=fIs/(x4*ga*rI*cI*(gp*p/pI+gm)**(x3/x2))

         fII=(p-pII)/(rII*cII*dsqrt(gp*p/pII+gm))
         fIIs=(ga+x1)*p/pII+(x3*ga-x1)
         fIIs=fIIs/(x4*ga*rII*cII*(gp*p/pII+gm)**(x3/x2))

         if(kiter.eq.1100) then
!        if(mode.ne.0) then
!            print *, 'zaciklilsya v raspade,i,j,k,KOBL,kdir:'
!            print *, i,j,k,KOBL,kdir
!            print *, pn,pee,p,dabs(pn/pee-p/pee),eps
!         endif
!        stop
           idgod=2
           return
         endif
         pn=p-(fI+fII-(enI-enII))/(fIs+fIIs)

         if(dabs(pn/pee-p/pee).ge.eps) goto 1
         p=pn


         return
         end


!****************************************************************
!**                                                            **
!**      devtwo  for                                           **
!**                                                            **
!****************************************************************
!@cuf attributes(host, device) & 
        subroutine devtwo(enI ,pI ,rI, &
                         enII,pII,rII,i,j,k,KOBL,kdir,idgod,ga,&
                         w   ,p)
!       include 'cyl.fl'
      implicit real*8 (a-h,o-z)
                real*8  kl,kp,kc,ksi,ksir,ksit

      data x0,x1,x2,x3,x4,x5,x6,x7,x8,x9/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./
      x10=10.d0
                epsil=1.D-9
!               epsil=1.D-11


                kl =  pI
                kp =  pII


	      	call  lev  ( enI,pI,rI,enII,pII,rII,kl,ksi,ga)
      		call  lev  ( enI,pI,rI,enII,pII,rII,kp,ksir,ga)


	        if (dabs(ksi) .le.epsil) then
                                    um=kl
				    goto 1
	        end if

	        if (dabs(ksir).le.epsil)  then
                                    um=kp
                                    goto 1
	        end if

          kpizd=0
    2 continue
          kpizd=kpizd+1
          if(kpizd.eq.1100) then
!         if(mode.gt.0) then
!            print *, 'zaciklilsya, devtwo.f i,j,k,KOBL,kdir:'
!            print *, i,j,k,KOBL,kdir
!            print *, ksi,ksir,kc,epsil
!         endif
!         stop
           idgod=2
           return
          endif


                kc=(kl+kp)/x2

 		call  lev  ( enI,pI,rI,enII,pII,rII,kc,ksit,ga)

                if (abs(ksit) .le. epsil)         goto 3
	        if ((ksi*ksit) .le. 0.d0)   then
                                     kp  =kc
				     ksir=ksit
				else
                                     kl  =kc
				     ksi =ksit
	        end if

	        goto 2

    3           um=kc
    1 continue


		p=um

      return
	end

!*******************************************************************
	!@cuf attributes(host, device) & 
      subroutine lev (enI,pI,rI,enII,pII,rII,uuu,fee,ga)
!	include    'cyl.fl'
      implicit real*8 (a-h,o-z)

      data x0,x1,x2,x3,x4,x5,x6,x7,x8,x9/0.,1.,2.,3.,4.,5.,6.,7.,8.,9./
      x10=10.d0

         g1=ga-1.d0
         g2=ga+1.d0

         gp=(ga+x1)/ga/x2
         gm=(ga-x1)/ga/x2
         cI =dsqrt(ga*pI /rI )
         cII=dsqrt(ga*pII/rII)

         fI=(uuu-pI)/(rI*cI*dsqrt(gp*uuu/pI+gm))

         fII=x2/g1*cII*((uuu/pII)**gm-x1)

		f1= fI+fII
      		f2= enI-enII
		fee=f1-f2

        return
	end subroutine lev






end module cgod


