      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ, cutFactor, 
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7)
c*************************************************************************
c
c      Author: Kamil Banas    
c      Supervisor: Janusz Badur      
c      The Szewalski Institute of Fluid Flow Machinery POLISH ACADEMY OF SCIENCE
c      
ccc
c      
c       The following demonstrates a USERMAT subroutine for
c       a plasticity model of 3D solid elements
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ,   cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c***************** User defined part *************************************
c
c --- parameters
c
      INTEGER          mcomp, NEWTON
      DOUBLE PRECISION HALF, THIRD, ONE, TWO, SMALL, ONEHALF,
     &                 ZERO, TWOTHIRD, ONEDM02, ONEDM05, sqTiny,
     &				  TOL, k, THREE, FOUR, SIX
      PARAMETER       (ZERO       = 0.d0,
     &                 HALF       = 0.5d0,
     &                 THIRD      = 1.d0/3.d0,
     &                 ONE        = 1.d0,
     &                 TWO        = 2.d0,
     &                 THREE      = 3.d0,
     &                 FOUR       = 4.d0,
     &                 SIX        = 6.d0,
     &                 SMALL      = 1.d-08,
     &                 sqTiny     = 1.d-20,
     &                 ONEDM02    = 1.d-02,
     &                 ONEDM05    = 1.d-05,
     &                 ONEHALF    = 1.5d0,
     &                 TWOTHIRD   = 2.0d0/3.0d0,
     &                 mcomp      = 6,
     &		         TOL        = 1.001d0,
     &				 NEWTON     = 50,
     &                 k          = 1.13d0
     &                 )
c
c --- local variables
c
      EXTERNAL         vzero, vmove, get_ElmData
      DOUBLE PRECISION sigElp(mcomp), dsdeEl(mcomp,mcomp), G(mcomp),
     &                 sigDev(mcomp), JM    (mcomp,mcomp), dfds(mcomp),
     &                 sigi  (mcomp), strainEl(mcomp),
     &		         depsPL(mcomp)

      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7

      DATA G/1.0D0,1.0D0,1.0D0,0.0D0,0.0D0,0.0D0/
c
      INTEGER          i, j, size, interp
      DOUBLE PRECISION pEl,   qEl,     pleq_t,  sigy_t , sigy,
     &                 dpleq, pleq, 
     &                 young, posn,    sigy0,   dsigdep, 
     &                 elast1,elast2,
     &                 twoG,  threeG,  oneOv3G, qElOv3G, threeOv2qEl, 
     &                 fratio,  con1,    con2, dperr(3),
cccccccccccc own variable
     &		       S_m, S_e, S_mtr, S_etr, S_B_eq,
     &		       dF_dSe, dF_dSm, d2F_dSmdSe, d2F_d2Sm, d2F_d2Se,
     &		       ddeltaEpl_eq_ddelta_Ep, ddeltaEpl_eq_ddelta_Eq,
     &		       A11, A12, A21, A22, b1, b2, 
     &		       deltaEp, deltaEq, d_deltaEp, d_deltaEq, 
     &		       b_1, b_2,
     &		       BULK, mi, H, pierw, temp2,
     &               Temp_end
      
       DOUBLE PRECISION   Actual_Temp_plastic_curve (nProp)
ccccccccccccccccccccc
      DOUBLE PRECISION d0, d1, d2, d3, d4, C_11, C_12, C_21, C_22,
     &                A_11, A_12, A_21, A_22, B_11, B_12, B_21, B_22,
     &                delta_C
      
      DOUBLE PRECISION NN(6)
c*************************************************************************
c************************************************************************* 
cccccccccccc
        

cccccccccccc
      keycut   = 0
      cutFactor = 0.d0
      H  = ZERO 
      pleq_t   = ustatev(1)
      pleq     = pleq_t
c *** get Young's modulus and Poisson's ratio, initial yield stress and others
      interp = prop(3)
      Temp_end = dTemp + Temp
      
      call material_properties( nProp, prop, Temp_end,
     &                                young, posn)
c *** calculate the plastic slope

      call Plastic_Curve( Temp_end, nProp, prop,
     &                            Actual_Temp_plastic_curve, size)
ccc
      sigy0 = Actual_Temp_plastic_curve(2)
      
      call HARD_SLOPE(nProp, Actual_Temp_plastic_curve,size,pleq,interp,
     &                                    H, sigy)
cc    

      twoG     = young / (ONE+posn) 
      threeG   = ONEHALF * twoG
      mi = young/(TWO*(1.0D0+posn)) 
      BULK = young / (THREE*(ONE-TWO*posn))
c *** plastic strain tensor
      call vmove(ustatev(2), epsPl(1), ncomp)

      elast1=young*posn/((1.0D0+posn)*(1.0D0-TWO*posn))
      elast2=HALF*twoG
c *** define tsstif(1) since it is used for calculation of hourglass stiffness
      tsstif(1) = elast2
c
c *** calculate elastic stiffness matrix (3d)
c
      dsdeEl(1,1)=(elast1+TWO*elast2)*G(1)*G(1)
      dsdeEl(1,2)=elast1*G(1)*G(2)+elast2*TWO*G(4)*G(4)
      dsdeEl(1,3)=elast1*G(1)*G(3)+elast2*TWO*G(5)*G(5)
      dsdeEl(1,4)=elast1*G(1)*G(4)+elast2*TWO*G(1)*G(4)
      dsdeEl(1,5)=elast1*G(1)*G(5)+elast2*TWO*G(1)*G(5)
      dsdeEl(1,6)=elast1*G(1)*G(6)+elast2*TWO*G(4)*G(5)
      dsdeEl(2,2)=(elast1+TWO*elast2)*G(2)*G(2)
      dsdeEl(2,3)=elast1*G(2)*G(3)+elast2*TWO*G(6)*G(6)
      dsdeEl(2,4)=elast1*G(2)*G(4)+elast2*TWO*G(1)*G(4)
      dsdeEl(2,5)=elast1*G(2)*G(5)+elast2*TWO*G(1)*G(5)
      dsdeEl(2,6)=elast1*G(2)*G(6)+elast2*TWO*G(2)*G(6)
      dsdeEl(3,3)=(elast1+TWO*elast2)*G(3)*G(3)
      dsdeEl(3,4)=elast1*G(3)*G(4)+elast2*TWO*G(5)*G(6)
      dsdeEl(3,5)=elast1*G(3)*G(5)+elast2*TWO*G(5)*G(3)
      dsdeEl(3,6)=elast1*G(3)*G(6)+elast2*TWO*G(6)*G(3)
      dsdeEl(4,4)=elast1*G(4)*G(4)+elast2*(G(1)*G(2)+G(4)*G(4))
      dsdeEl(4,5)=elast1*G(4)*G(5)+elast2*(G(1)*G(6)+G(5)*G(4))
      dsdeEl(4,6)=elast1*G(4)*G(6)+elast2*(G(4)*G(6)+G(5)*G(2))
      dsdeEl(5,5)=elast1*G(5)*G(5)+elast2*(G(1)*G(3)+G(5)*G(5))
      dsdeEl(5,6)=elast1*G(5)*G(6)+elast2*(G(4)*G(3)+G(5)*G(6))
      dsdeEl(6,6)=elast1*G(6)*G(6)+elast2*(G(2)*G(3)+G(6)*G(6))
      do i=1,ncomp-1
        do j=i+1,ncomp
          dsdeEl(j,i)=dsdeEl(i,j)
        end do
      end do
c
c
c *** get initial stress
      call vzero(sigi(1),ncomp)
      i = ncomp
      call get_ElmData ('ISIG', elemId,kDomIntPt, i, sigi)
c
c *** calculate the trial stress and
c     copy elastic moduli dsdeEl to material Jacobian matrix
      do i=1,ncomp
         strainEl(i) = Strain(i) + dStrain(i) - epsPl(i)
      end do
      call vzero(sigElp, 6)
      do i=1,ncomp
         do j=1,ncomp
            dsdePl(j,i) = dsdeEl(j,i)
            sigElp(i) = sigElp(i)+dsdeEl(j,i)*strainEl(j)
         end do
         sigElp(i) = sigElp(i) + sigi(i)
      end do

c *** hydrostatic pressure stress 
      S_mtr = THIRD * (sigElp(1) + sigElp(2) + sigElp(3)) 
c *** compute the deviatoric stress tensor 
      sigDev(1) = sigElp(1) - S_mtr 
      sigDev(2) = sigElp(2) - S_mtr 
      sigDev(3) = sigElp(3) - S_mtr 
      sigDev(4) = sigElp(4) 
      sigDev(5) = sigElp(5) 
      sigDev(6) = sigElp(6) 
c *** compute von-mises stress trial stress
      S_etr =  
     &  sigDev(1) * sigDev(1)+sigDev(2) * sigDev(2)+ 
     &  sigDev(3) * sigDev(3)+ 
     &  TWO*(sigDev(4) * sigDev(4)+ sigDev(5) * sigDev(5)+  
     &  sigDev(6) * sigDev(6)) 
      S_etr = sqrt( ONEHALF * S_etr) 
c ********* compute Burzynski trial stress     
      S_B_eq = (HALF/k)* (3.0*(k-1.0)*S_mtr +
     &  sqrt( 9.0*((k-1.0)**2.0)*(S_mtr**2.0) +4.0*k*(S_etr**2.0) ))
c *********      
 
      fratio = S_B_eq / sigy - ONE

c *** check for yielding 
      IF (sigy .LE. ZERO.or.fratio .LE. -SMALL) GO TO 500 
c 
      sigy_t = sigy 
      threeOv2qEl = ONEHALF / S_etr 
c *** compute derivative of the yield function 
      DO i=1, ncomp 
         dfds(i) = threeOv2qEl * sigDev(i) 
      END DO 
      oneOv3G  = ONE / threeG 
      qElOv3G  = S_etr * oneOv3G
c   
c 
      S_e = S_etr
      S_m = S_mtr

      deltaEq = 0
      deltaEp = 0
      
c *** Newton-Raphosn procedure for return mapping iteration 
	DO i = 1,NEWTON
c         
          
	   d_deltaEp = ZERO
	   d_deltaEq = ZERO
	   
	   pierw = sqrt(9.0D0*((k-1.0D0)**TWO)*(S_m**2.0D0) + 4.0D0*k*(S_e**2))
		
	   dF_dSe = TWO*S_e/pierw 
	   dF_dSm = ONEHALF*((k-1.0D0)/k) *(1.0D0 + 3.0D0*(k-1.0D0)*S_m/pierw)
         
         temp2 = 18.0D0*((k-1.0D0)**TWO)/(pierw**THREE)
	   d2F_dSmdSe = -S_e*S_m*temp2
	   d2F_d2Sm = (S_e**TWO)*temp2
	   d2F_d2Se = (S_m**TWO)*temp2
         
	   ddeltaEpl_eq_ddelta_Ep = S_m / sigy
	   ddeltaEpl_eq_ddelta_Eq = S_e / sigy

	   A11 = dF_dSe + deltaEp * ( - BULK * d2F_dSmdSe) +
     &	   		deltaEq * (  BULK * d2F_d2Sm)
	   A12 = -dF_dSm + deltaEp * ( - threeG * d2F_d2Se) +
     &	   		deltaEq * (  threeG  * d2F_dSmdSe)
	   A21 = - BULK * dF_dSm - H * ddeltaEpl_eq_ddelta_Ep
	   A22 =  - threeG * dF_dSe - H * ddeltaEpl_eq_ddelta_Eq
	   
	   b1 = - deltaEp * dF_dSe + deltaEq * dF_dSm
	   b2 = -(HALF/k)* 
     &(THREE*(k-1.0D0)*S_m + sqrt( 9.0D0*((k-1.0D0)**2.0D0)*(S_m**2.0D0)
     &                + 4.0D0*k*(S_e**2.0D0) ))  + sigy
         
	   d_deltaEp = ( b1*A22 - b2*A12 ) / ( A11*A22 - A21*A12 )
	   d_deltaEq = ( A11*b2 - A21*b1 ) / ( A11*A22 - A21*A12 )

	   
	   deltaEp = deltaEp + d_deltaEp
	   deltaEq = deltaEq + d_deltaEq
	   
	   S_m = S_mtr - BULK * deltaEp
	   S_e = S_etr - threeG * deltaEq
       
         
	   dpleq =(S_m*deltaEp + S_e*deltaEq ) / sigy

	   	   
c *** Avoid negative equivalent plastic strain
         dpleq  = max (dpleq, sqTiny) 
         pleq   = pleq_t + dpleq 
      
      call HARD_SLOPE(nProp, Actual_Temp_plastic_curve,size,pleq,interp,
     &                                    H, sigy)
   
c *** Check Convergence
	   b_2 = -(HALF/k)* 
     &(THREE*(k-1.0D0)*S_m + sqrt( 9.0D0*((k-1.0D0)**2.0D0)*(S_m**2.0D0)
     &                + 4.0D0*k*(S_e**2.0D0) ))  + sigy
	   b_1 = - deltaEp * dF_dSe + deltaEq * dF_dSm
c *** Check Convergence
	   IF ( (abs(b_2) .LT. 0.00000001 ).AND.
     &		(abs(b_1) .LT. 0.00000001 )) GO TO 100
	END DO
 
c 
c *** Uncovergence, set keycut to 1 for bisect/cut 
      keycut   = 1 
      GO TO 990 
 100  CONTINUE 
c 
c ***  update stresses 
      con1 =  twoG * dpleq 
      DO i = 1 , nDirect 
         stress(i) =  S_m + S_e*sigDev(i)/S_etr 
      END DO 
      
      DO i = nDirect + 1 , ncomp
         stress(i) =  S_e*sigDev(i)/S_etr
      END DO 
c 
c ***  update plastic strains 
      DO i = 1 , nDirect 
         epsPl(i) = epsPl(i) + THIRD * deltaEp  +
     &                          deltaEq*ONEHALF*sigDev(i)/S_etr 
      END DO 
      DO i = nDirect + 1 , ncomp 
         epsPl(i) = epsPl(i)  + TWO*deltaEq*ONEHALF*sigDev(i)/S_etr 
      END DO 

c ***********************
c *** Update state variables
      ustatev(1) = pleq
      epseq     = pleq
      do i=1,ncomp
         ustatev(i+1) = epsPl(i)
      end do
c *** Update plastic work
      sedPl = sedPl + HALF * (sigy_t+sigy)*dpleq
c
c *** Material Jcobian matrix 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc      
	pierw = sqrt(9.0D0*((k-1.0D0)**TWO)*(S_m**2.0D0) + 4.0D0*k*(S_e**2))
		
	dF_dSe = TWO*S_e/pierw 
	dF_dSm = ONEHALF*((k-1.0D0)/k) *(1.0D0 + 3.0D0*(k-1.0D0)*S_m/pierw)
         
      temp2 = 18.0D0*((k-1.0D0)**TWO)/(pierw**THREE)
	d2F_dSmdSe = -S_e*S_m*temp2
	d2F_d2Sm = (S_e**TWO)*temp2
	d2F_d2Se = (S_m**TWO)*temp2      
cccccccc      
      
     
      A_11 = dF_dSe
      A_12 = -dF_dSm
      A_21 = -H * (S_m / sigy)
      A_22 = -H * S_e / sigy
      
      B_11 = -THIRD*deltaEp*d2F_dSmdSe + THIRD*deltaEq*d2F_d2Sm
      B_12 = - deltaEp*d2F_d2Se + deltaEq*d2F_dSmdSe
      B_21 = THIRD*(-dF_dSm + H* ( deltaEp / sigy) )
      B_22 = -(dF_dSe -H* (deltaEq / sigy)  )      
      

      
      delta_C = (A_11 + 3*BULK*B_11)*(A_22 + 3*mi*B_22) - 
     &                        (A_12 + 3*mi*B_11)*(A_21 + 3*BULK*B_21)
      C_11 = ((A_22 + THREE*mi*B_22)  *B_11 - 
     &                            (A_12 + THREE*mi*B_12)*B_21)/delta_C
      C_21 = ((A_11 + THREE*BULK*B_11)*B_21 - 
     &                        (A_21 + THREE*BULK*B_21)*B_11)/delta_C
      C_12 = ((A_22 + THREE*mi*B_22)  *B_12 - 
     &                        (A_12 + THREE*mi*B_12)*B_22)/delta_C
      C_22 = ((A_11 + THREE*BULK*B_11)*B_22 - 
     &                        (A_21 + THREE*BULK*B_21)*B_12)/delta_C 
      
      d0 = twoG*S_e/S_etr
      d1 = BULK - THIRD*twoG*S_e/S_etr - THREE*(BULK**2)*C_11
      d2 = FOUR*(mi**2)*deltaEq/S_etr - FOUR*(mi**2)*C_22
      d3 =  -twoG*BULK*C_12
      d4 = -SIX*mi*BULK*C_21

            
      do i=1,6
          NN(i) = ONEHALF* (sigDev(i)/S_etr )
      end do
ccccc dodaje d2       
      do i=1,6
          do j=1,6
              if( (i.ge.4) .or. (j.ge.4) ) then
                  dsdePl(i,j) = HALF *d2* NN(i)*NN(j)
                  GOTO 200
              end if
              
              dsdePl(i,j) = d2* NN(i)*NN(j)
              
  200 continue 
          end do
      end do
ccccc dodaje d0 i d1  dla 1 cwiartki   
      do i=1,3
          do j=1,3
                    dsdePl(i,j) = dsdePl(i,j) + d1
          end do
          
          dsdePl(i,i) = dsdePl(i,i) + d0
      end do
ccccc dodaje d0 dla 4 cwiartki      
      do i=4,6
          dsdePl(i,i) = dsdePl(i,i) + HALF*d0
      end do
ccccc dodaje d3 do gornej polowki    
      do i=1,3
          do j=1,6

              if(j .ge. 4) then
                      dsdePl(i,j) = dsdePl(i,j) + HALF*d3* NN(j)
                      GOTO 300
              end if

              dsdePl(i,j) = dsdePl(i,j) + d3*NN(j)
  300 continue 
              
           end do
      end do


ccccc dodaje d4 do lewej polowki   
      do i=1,6
          do j=1,3

              if(i .ge. 4) then
                      dsdePl(i,j) = dsdePl(i,j) + HALF*d4* NN(i)
                      GOTO 400
              end if

              dsdePl(i,j) = dsdePl(i,j) + d4*NN(i)
  400 continue 
              
           end do
      end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      goto 600 
  500 continue

c *** Update stress in case of elastic/unloading
      do i=1,ncomp
         stress(i) = sigElp(i)
      end do

  600 continue
      sedEl = ZERO
      DO i = 1 , ncomp
         sedEl = sedEl + stress(i)*(Strain(i)+dStrain(i)-epsPl(i))
      END DO
      sedEl    = sedEl * HALF
      ustatev(nStatev) = sigy
cccccccccc               
cccccccccc      
  990 continue
c
!      cutFactor = 0.5d0
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine Plastic_Curve( Temp, nProp, prop,
     &                            Actual_Temp_plastic_curve, size) 
c      
       INTEGER     nProp, size, size_i, size_01, size_02,
     &                i, j, m, licznik, 
     &                Point_Temp, Point_Temp1, Point_Temp2
      
       DOUBLE PRECISION   prop (nProp), Actual_Temp_plastic_curve(nProp)
       
       DOUBLE PRECISION  Temp,
     &                   T1, T2, 
     &                   sigT1, sigT2, sig, temporaryA_sig,
     &                   max_eps_pl, max_eps_pl_01, max_eps_pl_02, 
     &                   d_eps_pl    
 
       
      size = 0
cccccccccccc
cccccccccccc jezeli jeden punkt plastycznosc
        IF (prop(2) . eq. 1) THEN
              Point_Temp = 3 + prop(1) + 2*prop(2) + 3*prop(1) + 1
              size = 2*prop(3 + prop(1) + prop(2) + 1)

              do j=1,size
                  Actual_Temp_plastic_curve(j) = prop(Point_Temp + j)
              end do 
          GOTO 10
       END IF      
ccccccccccccccc       
ccccccccccccccc       
       
 
c "size" podaje  ilosc punktow plastycznosci dla aktualnej krzywej, tzn. (eps_pl + sigma)
       licznik = 0
c "licznik" zlicza liczbe punktow plastycznosci dla poprzedniego kroku "i"         
       DO i=1, ( prop(2)-1 )
c "i" petla po temperaturach dla zakresu plastycznego         
           
          IF (   Temp.le.prop(3 + prop(1) + 1)  ) THEN
c jezeli Temp jest mniejsza lub rowna od pierwszego punktu temperatutowego               
              Point_Temp = 3 + prop(1) + 2*prop(2) + 3*prop(1) + 1
              size = 2*prop(3 + prop(1) + prop(2) + 1)

              do j=1,size
                  Actual_Temp_plastic_curve(j) = prop(Point_Temp + j)
              end do
 
              GOTO 10
               
          ELSE IF (Temp .eq. prop(3 + prop(1) + i) ) THEN
c jezeli Temp jest rowna danemu punktowi temperatutowemu                
              Point_Temp = 3 + prop(1) + 2*prop(2) +
     &                                3*prop(1) + (i-1) + licznik + 1
              size = 2*prop(3 + prop(1) + prop(2) + i)
             
              do j=1,size
                Actual_Temp_plastic_curve(j) = prop(Point_Temp + j)
              end do
 
              GOTO 10
               
          ELSE IF ( (Temp .gt. prop(3 + prop(1) + i) ) .and. 
     &              (Temp .lt. prop(3 + prop(1) + i +1) )    ) THEN
c jezeli Temp jest w przedziale
cccccccccccccccccccccccccccccccccccc
              Point_Temp1 = 3 + prop(1) + 2*prop(2) +
     &                          3*prop(1) + (i-1) + licznik + 1
              size_01 = 2*prop(3 + prop(1) + prop(2) + i )
              
              Point_Temp2 = Point_Temp1 + size_01 + 1
              size_02 = 2*prop(3 + prop(1) + prop(2) + i + 1)
    
              call interpolate(nProp, prop, Temp, Point_Temp1,
     &                        Point_Temp2, size_01, size_02, 
     &                               Actual_Temp_plastic_curve, size)

c             
              GOTO 10
          ELSE IF ( ( i .eq. (prop(2) - 1) ) .and. 
     &                     (Temp .eq. prop(3 + prop(1)+ i + 1) ) ) THEN
c jezeli Temp jest rowna ostaniemu punktowi temperaturowemu
                 
              size_i = 2*prop(3 + prop(1) + prop(2) + i )
              Point_Temp = 3 + prop(1) + 2*prop(2) +
     &                          3*prop(1) + i + licznik + size_i + 1
             
              size = 2*prop(3 + prop(1) + prop(2) + i +1 ) 
             
              do j=1,size
                  Actual_Temp_plastic_curve(j) = prop(Point_Temp + j)
              end do
 
              GOTO 10
               
          ELSE IF (  ( i .eq. (prop(2)-1)  ) .and. 
     &                ( Temp .gt. prop(3+prop(1)+i+1) )   ) THEN    
c jezeli Temp jest poza danymi temperaturowmi
             
              if (prop(3) .eq. 0) then
                  
                   size_i = 2*prop(3 + prop(1) + prop(2) + i )
                   Point_Temp = 3 + prop(1) + 2*prop(2) +
     &                          3*prop(1) + i + licznik + size_i + 1
                   
                   size = 2*prop(3 + prop(1) + prop(2) + i +1)
             
                   do j=1,size
                      Actual_Temp_plastic_curve(j) =prop(Point_Temp + j)
                   end do

                   GOTO 10                 
              else if (prop(3) .eq. 1) then
ccccccccccccccc                  
                  Point_Temp1 = 3 + prop(1) + 2*prop(2) +
     &                          3*prop(1) + (i-1) + licznik + 1
                  size_01 = 2*prop(3 + prop(1) + prop(2) + i )
              
                  Point_Temp2 = Point_Temp1 + size_01 + 1
                  size_02 = 2*prop(3 + prop(1) + prop(2) + i + 1)
                            
                  call interpolate(nProp, prop, Temp, Point_Temp1,
     &                        Point_Temp2, size_01, size_02, 
     &                               Actual_Temp_plastic_curve, size)
 
                   GOTO 10
               end if
          END IF
      
          licznik = licznik + (  2*prop(3+ prop(1) +prop(2) + i)  )
c licznik zlicza liczbe punktow plastycznosci dla poprzedniego kroku "i"     
       END DO
       
 10    CONTINUE 
c
      
c       
      return 
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine HARD_SLOPE(nProp, Actual_Temp_plastic_curve, size,pleq,
     &                        interp,    H, sigy) 

       INTEGER     nProp, size, interp,
     &             Points_number, i
      
       DOUBLE PRECISION   Actual_Temp_plastic_curve (nProp)
       
        DOUBLE PRECISION  pleq, H,
     &                    pleq0, pleq1, delta_pleq,
     &                    syield0, syield1, delta_syield, sigy
c
       H = 0.00
 
       
       Points_number = size/2
       DO i=1, Points_number - 1
           pleq1 = Actual_Temp_plastic_curve(i*2 + 1)
c pleq1 -> wskazuje odksztalcenia plastyczne prawej strony aktualnego przedzialu           
           IF (pleq.le.pleq1) THEN
c pleq0 -> wskazuje odksztalcenia plastyczne lewej strony aktualnego przedzialu                
               pleq0 = Actual_Temp_plastic_curve(i*2 - 1)
               delta_pleq = pleq1 - pleq0
               
               syield0 = Actual_Temp_plastic_curve(i*2)
               syield1 = Actual_Temp_plastic_curve(i*2 + 2)
               delta_syield = syield1 - syield0
               H = delta_syield/delta_pleq
               sigy = syield0 + H*(pleq - pleq0)
               GOTO 10
           ELSE IF (i .eq. Points_number - 1) THEN
              if (interp .eq. 0) then
                  H = 0.0
                  syield1 = Actual_Temp_plastic_curve(i*2 + 2)
                  sigy = syield1
                  GOTO 10
              else if (interp .eq. 1) then
                  pleq0 = Actual_Temp_plastic_curve(i*2 - 1)
                  delta_pleq = pleq1 - pleq0
               
                  syield0 = Actual_Temp_plastic_curve(i*2)
                  syield1 = Actual_Temp_plastic_curve(i*2 + 2)
                  delta_syield = syield1 - syield0
                  H = delta_syield/delta_pleq
                  sigy = syield1 + H*(pleq - pleq1)
                  GOTO 10
              end if
           END IF
      END DO

 10   CONTINUE
c       
      return 
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine material_properties(nProp, prop, Temp,
     &                                                    E, mi)
      
      
      INTEGER     nProp,i, j

      DOUBLE PRECISION   prop (nProp)
       
      DOUBLE PRECISION    Temp, Point_Temp, Point_temp1, Point_temp2,
     &                    temporaryA_E, temporaryA_mi,
     &                    T1, T2, E1, E2, mi1, mi2,
     &                    E, mi

      
cccccccc

       E = 0.00
       mi  = 0.00
ccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccc jezeli jeden punkt sprezystych danych
       IF (prop(1) . eq. 1) THEN
          Point_Temp = 3 + prop(1) + 2*prop(2) + 1
          E = prop(Point_Temp + 1)
          mi  = prop(Point_Temp + 2) 
          GOTO 10
       END IF
cccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccc
       DO i=1, ( prop(1)-1 )
c "i" petla po temperaturach zakresu sprezystego
      
          IF (   Temp.le.prop(3 + 1)  ) THEN
c jezeli Temp jest mniejsza lub rowna od pierwszego punktu temperaturowego
                  Point_Temp = 3 + prop(1) + 2*prop(2) + 1
                  E = prop(Point_Temp + 1)
                  mi  = prop(Point_Temp + 2)
                  GOTO 10
              
          ELSE IF (Temp .eq. prop(3+i) ) THEN
c jezeli Temp jest rowna danemu punktowi temperaturowemu 
  
                  Point_Temp = 3 + prop(1) + 2*prop(2) + 3*(i-1) + 1

                  E = prop(Point_Temp + 1)
                  mi  = prop(Point_Temp + 2)
                  GOTO 10          

          ELSE IF ( (Temp .gt. prop(3+i) ) .and. 
     &                            (Temp .lt. prop(3+i+1)) ) THEN
c jezeli Temp jest w przedziale  
                  Point_Temp1 = 3 + prop(1) + 2*prop(2) + 3*(i-1) + 1
                  Point_Temp2 = 3 + prop(1) + 2*prop(2) + 3*i + 1
                  
                  T1 = prop(Point_temp1)
                  E1 = prop(Point_temp1 + 1)
                  mi1 = prop(Point_temp1 + 2)                  

                  T2 = prop(Point_temp2)
                  E2 = prop(Point_temp2 + 1)
                  mi2 = prop(Point_temp2 + 2)
                          
                  temporaryA_E = (E1 - E2)/(T1 - T2)
                  temporaryA_mi = (mi1 - mi2)/(T1 - T2)                 
                  
                  E = Temp*temporaryA_E + E1 - temporaryA_E*T1
                  mi  = Temp*temporaryA_mi + mi1 - temporaryA_mi*T1 
                  GOTO 10

           ELSE IF ( ( i .eq. (prop(1) - 1) ) .and. 
     &                            (Temp .eq. prop(3 + i + 1) ) ) THEN
c jezeli Temp jest rowna ostaniemu punktowi temperaturowemu
                      
                  Point_Temp = 3 + prop(1) + 2*prop(2) + 3*i + 1

                  E = prop(Point_Temp + 1)
                  mi  = prop(Point_Temp + 2)
                  GOTO 10            
           ELSE IF ( (Temp .gt. prop(3 + i + 1) ) .and. 
     &                            ( i .eq. (prop(1) - 1) ) ) THEN
c jezeli Temp jest poza danymi temperaturowmi 

				if (prop(3) .eq. 0) then
                      
                      Point_Temp = 3 + prop(1) + 2*prop(2) + 3*i + 1
	                E = prop(Point_Temp + 1)
	                mi = prop(Point_Temp + 2)
					GOTO 10
                        
				else if (prop(3) .eq. 1) then
                      Point_Temp1 = 3 + prop(1) + 2*prop(2) + 3*(i-1) +1
                      Point_Temp2 = 3 + prop(1) + 2*prop(2) + 3*(i) +1

                      T1 = prop(Point_temp1)
                      E1 = prop(Point_temp1 + 1)
                      mi1 = prop(Point_temp1 + 2)                 

                      T2 = prop(Point_temp2)
                      E2 = prop(Point_temp2 + 1)
                      mi2 = prop(Point_temp2 + 2)
                          
                      temporaryA_E = (E1 - E2)/(T1 - T2)
                      temporaryA_mi = (mi1 - mi2)/(T1 - T2)
                  
                      E = Temp*temporaryA_E + E1 - temporaryA_E*T1
                      mi  = Temp*temporaryA_mi + mi1 - temporaryA_mi*T1 
                      GOTO 10
                      
				end if     
                  
           END IF
         
       END DO
      

      
 10   CONTINUE
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine interpolate(nProp, prop, Temp, Point_Temp1,
     &                        Point_Temp2, size_01, size_02, 
     &                                          Atpc, size)
      
      INTEGER     nProp,i, j, interp,
     &            size, size_01, size_02,
     &            Point_Temp, Point_Temp1, Point_Temp2,
     &            licznik, start_T2

      DOUBLE PRECISION   prop (nProp), Atpc(nProp)
       
      DOUBLE PRECISION    Temp, T1, T2,
     &                    max_eps_pl, max_eps_pl_01, max_eps_pl_02,
     &                    d_eps_pl, pleqT, pleqT1, pleqT2, 
     &                    sigT1, sigT2, temporaryA_sig, sig
              
              
      interp = prop(3)
      
      T1 = prop(Point_Temp1)
      T2 = prop(Point_Temp2)
      
      max_eps_pl_01 = prop(Point_Temp1 + size_01 -1)
      max_eps_pl_02 = prop(Point_Temp2 + size_02 -1)  

              
c dane dla pierwszego punktu
      Atpc(1) = 0.0

      sigT1 = prop(Point_Temp1 + 2)
      sigT2 = prop(Point_Temp2 + 2)
      
      temporaryA_sig = (sigT1 - sigT2)/(T1 - T2)
      sig = temporaryA_sig*Temp + sigT1 - temporaryA_sig*T1
      
      Atpc(2) = sig   
      size = 1
      licznik = 1
      
      IF (max_eps_pl_01 .ge. max_eps_pl_02) THEN
          start_T2 = 3
          do i=3,(size_01-1),2
c i petla po pierwszej krzywej, dla temperatury T1          
              pleqT1 = prop(Point_Temp1 + i)
              do j=start_T2,(size_02-1),2
c j petla po drugiej krzywej, dla temperatury T2              
                  pleqT2 = prop(Point_Temp2 + j)
                  if (pleqT2 .le. pleqT1) then
                      size = size + 2
                      Atpc(size) = pleqT2
                      start_T2 = start_T2 + 2
                      licznik = licznik + 1
                      GOTO 20
                  else if (pleqT1 .lt. pleqT2) then
                      size = size + 2
                      Atpc(size) = pleqT1
                      GOTO 30
c wychodzi z petli po 1 kroku
                  end if
                  
 20               CONTINUE
              end do
              
              if (licznik .eq. (size_02/2) ) then
                  size = size + 2
                  Atpc(size) = pleqT1
              end if
              
 30           CONTINUE              
          end do
          
          size = size + 1    
          GOTO 10  
      ELSE IF (max_eps_pl_01 .lt. max_eps_pl_02) THEN
ccccccccccccccccccc
          start_T1 = 3
          do i=3,(size_02-1),2
c i petla po pierwszej krzywej, dla temperatury T2          
              pleqT2 = prop(Point_Temp2 + i)
              do j=start_T1,(size_01-1),2
c j petla po drugiej krzywej, dla temperatury T1              
                  pleqT1 = prop(Point_Temp1 + j)
                  if (pleqT1 .le. pleqT2) then
                      size = size + 2
                      Atpc(size) = pleqT1
                      start_T1 = start_T1 + 2
                      licznik = licznik + 1
                      GOTO 40
                  else if (pleqT2 .lt. pleqT1) then
                      size = size + 2
                      Atpc(size) = pleqT2
                      GOTO 50
c wychodzi z petli po 1 kroku
                  end if
                  
 40               CONTINUE
              end do
              
              if (licznik .eq. (size_01/2) ) then
                  size = size + 2
                  Atpc(size) = pleqT2
              end if
              
 50           CONTINUE              
          end do
          
          size = size + 1 
          
ccccccccccccccccccc          
          GOTO 10
      END IF
      
 10   CONTINUE   

      call Interp_SIG(prop, nProp, Temp, Point_Temp1, Point_Temp2, size,
     &                            size_01, size_02,        Atpc)
      
c      if(elemId . eq. 1) then
c          if(TEMP .gt. 855)then
c            OPEN(144,FILE='My_output.dat') 
c             write(144,*) Temp,  Point_Temp1,Point_Temp2,size,
c     &        max_eps_pl_01, max_eps_pl_02, size_01, size_02,
c     &       Atpc(1), Atpc(2), Atpc(3), Atpc(4), Atpc(5), 
c     &       Atpc(6), Atpc(7), Atpc(8), Atpc(9), Atpc(10),  
c     &       Atpc(11), Atpc(12), Atpc(13), Atpc(14), Atpc(15), 
c     &       Atpc(16), Atpc(17), Atpc(18), Atpc(19), Atpc(20), 
c     &       Atpc(21), Atpc(22), Atpc(23), Atpc(24), Atpc(25),  
c     &       Atpc(26), Atpc(27), Atpc(28), Atpc(29), Atpc(30), 
c     &       Atpc(31), Atpc(32), Atpc(33), Atpc(34), Atpc(35), 
c     &       Atpc(36), Atpc(37), Atpc(38), Atpc(39), Atpc(40),
c     &       Atpc(41), Atpc(42)           
c            CLOSE(144)
c          end if
c      end if

      return
          end
ccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccc
      subroutine    Interp_SIG(prop, nProp, Temp, Point_Temp1, 
     &                 Point_Temp2, size, size_01, size_02,        Atpc)
      
      INTEGER     nProp, Point_Temp, Point_Temp1, Point_Temp2,
     &              size, size_01, size_02, interp, i, j,
     &              start_T1, start_T2

      DOUBLE PRECISION   prop (nProp), Atpc(nProp)
       
      DOUBLE PRECISION    Temp, T1, T2,
     &                    pleqT, pleqT_0, pleqT_1,
     &                    pleqT1_0, pleqT1_1,
     &                    sigT1_0, sigT1_1, temporaryA_sigT1, sigT1,
     &                    temporaryA_sigT, sigT
      
      
      interp = prop(3)
      T1 = prop(Point_Temp1)
      T2 = prop(Point_Temp2)
      
      start_T1 = 3
      start_T2 = 3
      
      DO i=3,(size-1),2
          
          pleqT = Atpc(i)

c znajduje sigT1 dla pierwszej temperatury          
          do m = start_T1,(size_01-1),2

              IF ( pleqT .le. prop(Point_Temp1 + m) ) THEN
                  
                  if (pleqT .eq. prop(Point_Temp1 + m)) then
                      sigT1 = prop(Point_Temp1 + m + 1)
                      start_T1 = start_T1 + 2
                      GOTO 20
                  else
                      pleqT1_0 = prop(Point_Temp1 + m-2)
                      pleqT1_1 = prop(Point_Temp1 + m)
                                  
                      sigT1_0 = prop(Point_Temp1 + m -2 + 1)
                      sigT1_1 = prop(Point_Temp1 + m + 1)
                                  
                      temporaryA_sigT1 =
     &                        (sigT1_0 - sigT1_1)/(pleqT1_0 - pleqT1_1)
                      sigT1 = temporaryA_sigT1 * pleqT + sigT1_0 -
     &                                       temporaryA_sigT1*pleqT1_0
                      GOTO 20
                  end if
                  
              ELSE IF (m . eq. (size_01-1) ) THEN
                  if (interp .eq. 0) then

                      sigT1 = prop(Point_Temp1 + size_01)
                      GOTO 20
                  else if (interp .eq. 1) then
                      pleqT1_0 = prop(Point_Temp1 + m-2)
                      pleqT1_1 = prop(Point_Temp1 + m)
                                  
                      sigT1_0 = prop(Point_Temp1 + m -2 + 1)
                      sigT1_1 = prop(Point_Temp1 + m + 1)
                                  
                      temporaryA_sigT1 =
     &                        (sigT1_0 - sigT1_1)/(pleqT1_0 - pleqT1_1)
                      sigT1 = temporaryA_sigT1 * pleqT + sigT1_0 -
     &                                       temporaryA_sigT1*pleqT1_0
                      GOTO 20
                  end if                              
              END IF
          end do
      
 20       CONTINUE
 
          
c znajduje sigT1 dla pierwszej temperatury          
          do m = start_T2,(size_02-1),2
              IF ( pleqT .le. prop(Point_Temp2 + m) ) THEN
                  
                  if (pleqT .eq. prop(Point_Temp2 + m)) then
                      sigT2 = prop(Point_Temp2 + m + 1)
                      start_T2 = start_T2 + 2
                      GOTO 30
                  else
                      pleqT2_0 = prop(Point_Temp2 + m-2)
                      pleqT2_1 = prop(Point_Temp2 + m)
                                  
                      sigT2_0 = prop(Point_Temp2 + m -2 + 1)
                      sigT2_1 = prop(Point_Temp2 + m + 1)
                                  
                      temporaryA_sigT2 =
     &                        (sigT2_0 - sigT2_1)/(pleqT2_0 - pleqT2_1)
                      sigT2 = temporaryA_sigT2 * pleqT + sigT2_0 -
     &                                       temporaryA_sigT2*pleqT2_0
                      GOTO 30
                  end if
                  
              ELSE IF (m . eq. (size_02-1) ) THEN
                  if (interp .eq. 0) then

                      sigT2 = prop(Point_Temp2 + size_02)
                      GOTO 30
                  else if (interp .eq. 1) then
                      pleqT2_0 = prop(Point_Temp2 + m-2)
                      pleqT2_1 = prop(Point_Temp2 + m)
                                  
                      sigT2_0 = prop(Point_Temp2 + m -2 + 1)
                      sigT2_1 = prop(Point_Temp2 + m + 1)
                                  
                      temporaryA_sigT2 =
     &                        (sigT2_0 - sigT2_1)/(pleqT2_0 - pleqT2_1)
                      sigT2 = temporaryA_sigT2 * pleqT + sigT2_0 -
     &                                       temporaryA_sigT2*pleqT2_0
                      GOTO 30
                  end if                              
              END IF
          end do
          
 30       CONTINUE
          
          temporaryA_sigT = (sigT1 - sigT2)/(T1 - T2)
c
          sigT = temporaryA_sigT*Temp + sigT1 - temporaryA_sigT*T1
          Atpc(i+1) = sigT  
          
      END DO
      
      return
      end
