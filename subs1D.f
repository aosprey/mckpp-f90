C ******************************************************************
C
      REAL FUNCTION CPSW(S,T1,P0)
C
C ******************************************************************
C UNITS:      
C       PRESSURE        P0       DECIBARS
C       TEMPERATURE     T        DEG CELSIUS (IPTS-68)
C       SALINITY        S        (IPSS-78)
C       SPECIFIC HEAT   CPSW     J/(KG DEG C)
C ***
C REF: MILLERO ET AL,1973,JGR,78,4499-4507
C       MILLERO ET AL, UNESCO REPORT NO. 38 1981 PP. 99-188.
C PRESSURE VARIATION FROM LEAST SQUARES POLYNOMIAL
C DEVELOPED BY FOFONOFF 1980.
C ***
C CHECK VALUE: CPSW = 3849.500 J/(KG DEG. C) FOR S = 40 (IPSS-78),
C T = 40 DEG C, P0= 10000 DECIBARS
c
c   check that temperature is above -2
      T = T1
      if(T.lt.-2.) T = -2.
c
C   SCALE PRESSURE TO BARS
      P=P0/10.
C ***
C SQRT SALINITY FOR FRACTIONAL TERMS
      SR = SQRT(ABS(S))
C SPECIFIC HEAT CP0 FOR P=0 (MILLERO ET AL ,UNESCO 1981)
      A = (-1.38385E-3*T+0.1072763)*T-7.643575
      B = (5.148E-5*T-4.07718E-3)*T+0.1770383
      C = (((2.093236E-5*T-2.654387E-3)*T+0.1412855)*T
     X    -3.720283)*T+4217.4
      CP0 = (B*SR + A)*S + C
C CP1 PRESSURE AND TEMPERATURE TERMS FOR S = 0
      A = (((1.7168E-8*T+2.0357E-6)*T-3.13885E-4)*T+1.45747E-2)*T
     X   -0.49592
      B = (((2.2956E-11*T-4.0027E-9)*T+2.87533E-7)*T-1.08645E-5)*T
     X   +2.4931E-4
      C = ((6.136E-13*T-6.5637E-11)*T+2.6380E-9)*T-5.422E-8
      CP1 = ((C*P+B)*P+A)*P
C CP2 PRESSURE AND TEMPERATURE TERMS FOR S > 0
      A = (((-2.9179E-10*T+2.5941E-8)*T+9.802E-7)*T-1.28315E-4)*T
     X   +4.9247E-3
      B = (3.122E-8*T-1.517E-6)*T-1.2331E-4
      A = (A+B*SR)*S
      B = ((1.8448E-11*T-2.3905E-9)*T+1.17054E-7)*T-2.9558E-6
      B = (B+9.971E-8*SR)*S
      C = (3.513E-13*T-1.7682E-11)*T+5.540E-10
      C = (C-1.4300E-12*T*SR)*S
      CP2 = ((C*P+B)*P+A)*P
C SPECIFIC HEAT RETURN
      CPSW = CP0 + CP1 + CP2
      RETURN
      END
C **********************************************************************
c Subroutine to Coordinate computation of the expansion coefficients of 
c seawater with respect to temperature (alpha), salinity (beta), and
c pressure (kappa) using the 1980 equation of state for seawater. 
c 
c Subroutines Bet80,Alf80,Kap80 are called IN THAT ORDER to assure a
c minimum of repetitive work between the subroutines, as well as to 
c minimize storage in common (alpha overwrites coefficients that beta 
c obtains through common from Sig80). 
c 
c On Entry, the calling program should specify a nonzero value for
c any or all of: alpha,beta,kappa in the argument list.  This indicates 
c to THIS subroutine which values are desired.
c 
c If either alpha or beta (or both) are requested, then Sig80 is called 
c to compute the density at S,T,P.  Intermediate quantities in that 
c computation are passed in common for usage by Bet80 (FIRST!) for very 
c little extra work.  Alf80 doesn't gain too much from Sig80, and must
c redo many temperature power series.  To save space, the intermediate
c coefficients in common are overwritten as Alpha is computed since 
c beta will have already used them (if it was computed at all). 
c 
c Since Kappa doesn't require density values to be known, but does
c require the Bulk Modulus, an Entry point to compute only K in Sig80 
c has been set up.  However, if alpha or beta are computed then K and 
c other necessary terms to compute kappa will be in common and the
c call to BLKMOD in Sig80 won't be necessary.  There is, however, a 
c special case when we must still call BlkMod in Sig80: if P=0.  This 
c results because the Bulk Modulus won't be needed to compute density 
c (sigma-t) but they WILL be needed to compute kappa.  In this case 
c we again go to the entry point BlkMod in Sig80. 
c 
c Coordination of whether or not density has been computed (Sig80 
c called) is handled by logical flag KapFlg.  This is used to control 
c entry BlkMod in Sig80, as well as to minimize extra work in Kappa if
c density is already known, so that K and the other terms needed for 
c Kappa are already in common.
c 
c Units used in this subroutine and external subroutines:
c Sig80,Alf80,Bet80,Kap80 are:
c
c S(1978 practical salinity),T(degC),P(dbar),
c Alpha(degC**-1),Beta(nondim. {divided by 10**-3 for 1978 prac. sal.}),
c Kappa(bar**-1),Sig & Sig0(kg/m**3).
c
c Check Values:  S=35,T=15(degC),P=0(dbar)-->Alpha=2.14136e-4,
c                                            Beta =7.51638e-4,
c                                            Kappa=4.32576e-5.
c
c                S=40,T=0(degC),P=10,000(dbar)-->Alpha=2.69822e-4,
c                                                Beta =6.88317e-4,
c                                                Kappa=3.55271e-5.
c
c References: UNESCO(1981): Background papers and supporting data on
c                           the International Equation of State of
c                           Seawater, 1980.  UNESCO Tech. Pap. in Mar.
c                           Sci., No. 38, 192 pp.
c
c             Fofonoff, N.P., and R.C. Millard, Jr. (1983): Algorithms
c                           for computation of fundamental properties
c                           of seawater.  UNESCO Tech. Pap. in Mar.
c                           Sci., No. 44, 53 pp.
c
c             Lillibridge, J.L. (1988): Computing the seawater expansion
c                           coefficients directly from the 1980 equation
c                           of state.  Jour. Atm. Ocean. Tech., in press.
c
c John L. Lillibridge @ URI/GSO: 25 March 1987 
c 
      Subroutine ABK80(S,T1,P,Alpha,Beta,Kappa,Sig0,Sig) 
c    *,Exp. Coeff. of SeaWater 1980 <870330.1557> 
c 
      Implicit None 
c 
      Real P,P0,T,S,SR,Sig,Sig0,R1,R2,R3,R4,T1 
      Real PK,A,B,Alpha,Beta,Kappa  
      Real A1,B1,C,D,E,K
cc    Real A1,B1,C,D,E,K,DelK 
      Real Rho,Rho0,ABFac 
      Logical KapFlg,ABFlg
c
      Common /EOS/R1,R2,R3,R4,A,B,C,D,E,A1,B1,K,SR,P0,PK,Rho,Rho0 
     *           ,ABFac,ABFlg 
#ifdef OPENMP
!$OMP THREADPRIVATE(/EOS/)
      INTEGER tid,omp_get_thread_num
      tid=OMP_GET_THREAD_NUM()
#endif

c
c
c Check that temperature is above freezing
      T = T1
      if(T.lt.-2.) T = -2.
c 
c Set the KapFlg to .TRUE. to indicate Sig80 not called yet 
c and set the ABFlg to .TRUE. to be reset only by Beta to 
c save alpha a little work
c 
      KapFlg=.True. 
      ABFlg=.True.
c 
c Begin with Beta 
c Note that this MUST be called before Alpha due to destruction 
c of Intermediate sums used between Sig80 and Beta
c 
c      WRITE(6,*) tid,'ABK80',S,T1,P,Alpha,Beta,Kappa,Sig0,Sig
      If(Beta.ne.0)Then      
*!Main Program wants Beta 
c         WRITE(6,*) tid,'Calling Sig80',S,T,P,KapFlg,Sig0,Sig
         Call Sig80(S,T,P,KapFlg,Sig0,Sig) 
c         WRITE(6,*) tid,'Calling Bet80',S,T,P,Beta
         Call Bet80(S,T,P,Beta)
c         WRITE(6,*) tid,'Called Bet80',S,T,P,Beta
      EndIf 
c 
c Compute Alpha next overwriting many common variables used by
c Sig80 
c 
      If(Alpha.ne.0)Then     
*!Main Program wants Alpha
          If(KapFlg)Then
c             WRITE(6,*) tid,'Call Sig80 ',S,T,P,KapFlg,Sig0,Sig
              Call Sig80(S,T,P,KapFlg,Sig0,Sig) 
          EndIf 
c          WRITE(6,*) tid,'Call Alf80',S,T,P,Alpha
          Call Alf80(S,T,P,Alpha) 
      EndIf 
c 
c If we haven't called Sig80 yet, zero out Sig,Sig0 to be safe. 
c 
      If(KapFlg)Then
          Sig=0.
          Sig0=0. 
      EndIf 
c 
c Lastly compute Kappa using quantities in common from Sig80 if 
c possible, otherwise go to entry point BlkMod in Sig80 from
c within Kap80
c 
      If(Kappa.ne.0)Then
c         WRITE(6,*) 'Call Kap80',S,T,P,KapFlg,Kappa
          Call Kap80(S,T,P,KapFlg,Kappa)
      EndIf 
c 
c All Done, Return with appropriate desired values calculated 
c 
      Return
      End 
c 
c Coefficient of Haline Contraction using algebraic derivation of 
c the formulae for the 1980 Equation of State.
c 
c Units: P(db),T(DegC),S(1978 Practical Salinity),Sig & Sg80(kg/m**3)
c        Rho & Rho0(kg/m**3), Beta (nondim. {divided by 10**-3 for
c        1978 prac. sal.})
c 
c Common variables are polynomials of T evaluated by Sig80 in 
c computing density, which are then used again by Beta.  In addition
c the sigma-t value (sig) and secant bulk modulus (K) are stored for
c use in the overall expression for Beta. 
c 
c John L. Lillibridge at URI/GSO: 12/22/86 
c 
      Subroutine Bet80(S,T,P,Beta)  
c    *,Coeff. of Haline Contraction <870330.1552> 
c 
      Implicit None 
c 
      Real P,P0,T,S,SR,R1,R2,R3,R4 
cc    Real P,P0,T,S,SR,Sig,Sg80,R1,R2,R3,R4 
      Real Beta,PK,A,B,SR5
      Real A1,B1,C,D,E,K
      Real Rho,Rho0 
      Real DRho,DK,DK0,DA,DB
      Real ABFac
      Logical ABFlg 
c 
      Common /EOS/R1,R2,R3,R4,A,B,C,D,E,A1,B1,K,SR,P0,PK,Rho,Rho0 
     *           ,ABFac,ABFlg
#ifdef OPENMP
!$OMP THREADPRIVATE(/EOS/)
#endif
c 
c Compute the Sigma-t derivative term 
c Note that SR is the Sq.Root of Salinity from Sig80
c 
      SR5=SR*1.5
      DRho=R2+SR5*R3+(S+S)*R4 
      If(P.eq.0)Then
          Beta=DRho/Rho 
          Return
      EndIf 
c 
c Next Compute the Derivative Terms for the Bulk Modulus
c 
      DK0=A1+SR5*B1 
      DA=C+SR5*D
      DB=E
c 
c Derivative DK/DS
c 
      DK=(DB*P0+DA)*P0+DK0
c 
c Assemble Beta from all the terms (Rho0,K,PK from Sig80) 
c 
      ABFac=Rho0*P0/((K-P0)*(K-P0)) 
      ABFlg=.False. 
      Beta=DRho/(1.-PK) - ABFac*DK
      Beta=Beta/Rho 
      Return
      End 
C 
C Thermal Expansion Coefficient for Sea Water 
C 
C*******************************************************************
      Subroutine Alf80(S,T,P,Alpha) 
C    *,1980 Thermal Exp. Coeff. <870330.1549> 
C*******************************************************************
C 
C Derivative of the 1980 equation of state with respect to
C Temperature.  Constants in the EOS have been replaced simply
C with c(i)-->i*c(i), and the powers in T are reduced by one
C in the Horner's Rule power series sums. 
C 
C UNITS: P(DBAR), T(DEG C), S(1978 Practical Salinity),
C        P0,K(BAR), Alpha (DegC-1)
C 
C Updated for use with alpha,beta,kappa package by JLL March 24, 1987 
C Including use of extensive common to save intermediate summations.
C 
C Note that Rho,Rho0,P0(Bars),K,PK(=P0/K) & SR will be defined on 
C entry from a previous evaluation of Rho 
C 
C John L. Lillibridge @ URI/GSO: March 24, 1987
C
      Implicit None 
      Real P,P0,PK                             
*!Pressure Related Terms
      Real T,S,SR                              
*!Temp & Salinity Terms 
      Real R1,R2,R3,R4                         
*!Atm. Prs. Density Terms  
      REAL A,B,C,D,E,A1,B1,AW,BW,K,K0,KW       
*!Bulk Modulus Terms
      Real Rho,Rho0                            
*!Sigma + 1000
      Real Alph0,AlphaA,AlphB,AlphK,Alpha      
*!Alpha Terms 
      Real ABFac                               
*!Common to Alpha/Beta
      Logical ABFlg                            
*!  "     "     " 
C 
C Common Block for the Equation Of State subroutines Alf80,Bet80,Kap80
C Which can utilize Intermediate sums to minimize overhead once density 
C computed by subroutine Sig80. 
C 
      COMMON /EOS/R1,R2,R3,R4,A,B,C,D,E,A1,B1,K,SR,P0,PK,Rho,Rho0 
     *           ,ABFac,ABFlg 
#ifdef OPENMP
!$OMP THREADPRIVATE(/EOS/)
#endif
C 
C COMPUTE Alpha PURE WATER AT ATM PRESS 
C 
      R1=(((.3268166E-7*T-.4480332e-5)*T+.3005055e-3)*T 
     X-.1819058E-1)*T+6.793952E-2 
C 
C SEAWATER Alpha AT ATM PRESS 
C 
      R2=((.215500E-7*T-.247401E-5)*T+.152876E-3)*T-4.0899E-3 
      R3=-.33092E-5*T+1.0227E-4 
      Alph0=(R3*SR+R2)*S+R1 
C 
C Alpha AT ATM PRESS
C 
      IF(P.EQ.0.0)Then
          Alpha=-Alph0/Rho
          Return
      EndIf 
C 
C COMPUTE COMPRESSION TERMS 
C 
      B1=-.106018E-2*T+1.6483E-2
      A1=(-.18501E-3*T+.219974E-1)*T-0.603459 
      KW=((-.2062115E-3*T+.4081431E-1)*T-.4654210E+1)*T 
     X+148.4206 
      K0=(B1*SR+A1)*S+KW
C 
C Pressure Terms in Bulk Modulus
C 
      E=.183394E-8*T+2.0816E-8
      BW=.105574E-6*T-6.12293E-6
      AlphB=BW+E*S
C 
      C=-.32156E-5*T-1.0981E-5
      AW=(-.1733715E-5*T+.232184E-3)*T+1.43713E-3 
      AlphaA=C*S+AW 
C 
C EVALUATE PRESS POLYNOMIAL AND RETURN
C 
      AlphK=(AlphB*P0+AlphaA)*P0+K0 
C 
      If(ABFlg) Then
          ABFac=Rho0*P0/((K-P0)*(K-P0)) 
      EndIf 
      Alpha=Alph0/(1.-PK) - ABFac*AlphK 
      Alpha=-Alpha/Rho
      RETURN
      END 
c 
c Coefficient of Compressibility using algebraic derivation of
c the formulae for the 1980 Equation of State.
c 
c Units: P(db),T(degC),S(1978 Practical Salinity),Sig & Sg80(kg/m**3)
c        Rho & Rho0(kg/m**3), Kappa (bar-1) 
c 
c If density has already been computed, for alpha or beta, then 
c the coefficients A,B and the secant bulk modulus K will have
c already been computed so simply evaluate the formula for kappa. 
c However, if P=0, then neither A,B, nor K will be computed in
c Sig80, therefore we must go to that subroutine at entry point 
c BlkMod to obtain those coefficients.  Similarly, if density 
c has not been computed, we go to entry BlkMod in subroutine
c Sig80 to obtain A,B,K for kappa.
c 
c John L. Lillibridge @ URI/GSO: 25 March 1987 
c 
      Subroutine Kap80(S,T,P,KapFlg,Kappa)  
c    *,Compressibility of Sea Water <870330.1553> 
c 
      Implicit None 
c 
      Real P,P0,T,S,SR,R1,R2,R3,R4 
cc    Real P,P0,T,S,SR,Sig,Sg80,R1,R2,R3,R4 
      Real PK,A,B,Kappa 
      Real A1,B1,C,D,E,K,DelK 
      Real Rho,Rho0,ABFac
      Logical KapFlg, ABFlg                    
*           Signals we need entry to BlkMod 
c 
      Common /EOS/R1,R2,R3,R4,A,B,C,D,E,A1,B1,K,SR,P0,PK,Rho,Rho0 
     *           ,ABFac,ABFlg 
#ifdef OPENMP
!$OMP THREADPRIVATE(/EOS/)
#endif
c 
c If P=0 must go to Sig80 Subroutine at Entry BlkMod
c 
      If(P.eq.0) Then 
          KapFlg=.True. 
          Call BlkMod(S,T,P,KapFlg) 
          Kappa=1.0/K                   
*!K in Common (Bars), Kappa in bar-1 
          Return
      EndIf 
c 
c If Density hasn't been computed must also go to BlkMod
c 
      If(KapFlg) Then 
          Call BlkMod(S,T,P,KapFlg) 
      EndIf 
c 
c Full nonzero pressure formula for Kappa given A,B,K 
c 
      DelK=A+(P0+P0)*B
      Kappa=(1.-PK*DelK)/(K-P0) 
c 
      Return
      End 
C 
C EQUATION OF STATE FOR SEA WATER 
C 
C*******************************************************************
      Subroutine Sig80(S,T,P,KapFlg,Sig0,Sig) 
C    *,1980 EQ. OF STATE <870330.1544>  
C*******************************************************************
C 
C EQUATION OF STATE FOR SEAWATER PROPOSED BY JPOTS 1980 
C 
C REF: MILLERO, ET. AL. (1980) D.S.R.,27A,255-264.
C UNITS: P(DBAR), T(DEG C), S(P.S.U.), Sig,Sig0(KG/M**3)  
C        P0,K(BAR)
C 
C NPF,OCT 7 80; JLL,FEB 24 82.
C UPDATED TO USE ONLY SIGMA UNITS BY JLL MAR 10 82. 
C 
C Updated for use with alpha,beta,kappa package by JLL March 24, 1987 
C Including use of extensive common to save intermediate summations.
C 
C John L. Lillibridge @ URI/GSO: March 24 1987
C
      Implicit None 
      Real P,P0,PK                             
*!Pressure Related Terms
      Real T,S,SR                              
*!Temp & Salinity Terms 
      Real Sig0,R1,R2,R3,R4                    
*!Atm. Prs. Density Terms 
      REAL A,B,C,D,E,A1,B1,AW,BW,K,K0,KW       
*!Bulk Modulus Terms
*      Real BlkMod                              
*!Entry Point for Kappa 
      Real Sig,Rho,Rho0,ABFac                    
*!Sigma + 1000
      Logical KapFlg, ABFlg                           
*!True when Entry at BlkMod 
C 
C Common Block for the Equation Of State subroutines Alf80,Bet80,Kap80
C Which can utilize Intermediate sums to minimize overhead once density 
C computed by this routine. 
C 
      Common /EOS/R1,R2,R3,R4,A,B,C,D,E,A1,B1,K,SR,P0,PK,Rho,Rho0 
     *           ,ABFac,ABFlg 
#ifdef OPENMP
!$OMP THREADPRIVATE(/EOS/)
#endif
C 
C CONVERT PRESSURE TO BARS AND SQ ROOT SALINITY 
C And Set the Logical Flag used to control Entry Point Below for Kappa
C 
      P0=P/10.0 
      SR=SQRT(ABS(S)) 
      KapFlg=.False.
C 
C COMPUTE SIGMA PURE WATER AT ATM PRESS 
C 
      R1=((((6.536332E-9*T-1.120083E-6)*T+1.001685E-4)*T
     X-9.095290E-3)*T+6.793952E-2)*T-.157406 
*!Note Constant for Sigma vs. Rho!
C 
C SEAWATER SIGMA AT ATM PRESS 
C 
      R2=(((5.3875E-9*T-8.2467E-7)*T+7.6438E-5)*T-4.0899E-3)*T
     X+8.24493E-1 
      R3=(-1.6546E-6*T+1.0227E-4)*T-5.72466E-3
      R4=4.8314E-4
      Sig0=(R4*S+R3*SR+R2)*S+R1 
      Rho0=1000.0 + Sig0
C 
C SIGMA AT ATM PRESS
C 
      IF(P.EQ.0.0)Then
          Sig=Sig0  
          Rho=Rho0
          Return
      EndIf 
C 
C COMPUTE COMPRESSION TERMS 
C 
C This is the entry point when the Bulk Modulus K is desired, which 
C doesn't require knowledge of density. Entry from Kappa will have
C KapFlg set .True. so that we exit before computing Sig. 
C 
      Entry BlkMod(S,T,P,KapFlg)           
*!For Computation of Kappa  
C 
C If Entry without density computed need these two common terms 
C 
      If(KapFlg)Then
          P0=P/10.0 
          SR=SQRT(ABS(S)) 
      EndIf 
C 
      B1=(-5.3009E-4*T+1.6483E-2)*T+7.944E-2
      A1=((-6.1670E-5*T+1.09987E-2)*T-0.603459)*T+54.6746 
      KW=(((-5.155288E-5*T+1.360477E-2)*T-2.327105)*T 
     X+148.4206)*T+19652.21 
      K0=(B1*SR+A1)*S+KW
C 
C If computing Kappa and P=0, return here 
C 
      If(P.eq.0.0)Then
          K=K0
          Return
      EndIf 
C 
C Pressure Terms in Bulk Modulus
C 
      E=(9.1697E-10*T+2.0816E-8)*T-9.9348E-7
      BW=(5.2787E-8*T-6.12293E-6)*T+8.50935E-5
      B=BW+E*S
C 
      D=1.91075E-4
      C=(-1.6078E-6*T-1.0981E-5)*T+2.2838E-3
      AW=((-5.77905E-7*T+1.16092E-4)*T+1.43713E-3)*T+3.239908 
      A=(D*SR+C)*S+AW 
C 
C EVALUATE PRESS POLYNOMIAL AND RETURN
C 
      K=(B*P0+A)*P0+K0
C 
C If Entry at BlkMod Exit Now 
C 
      PK=P0/K 
      If (KapFlg) Return
      Sig=(1000.0*PK+Sig0)/(1.0-PK)   
      Rho=1000.0 + Sig  
      RETURN
      END 
********************************************************************
      subroutine acopy(a,b,ndim)
c                                 copy array a to b
      dimension a(ndim),b(ndim)

      do 5 n=1,ndim
      b(n) = a(n)
  5   continue
      return
      end
******************************************************
      subroutine aset(array,ndim,value)
c                         set array to value 
      dimension array(ndim)
      do 5 n=1,ndim
 5    array(n) = value
      return
      end
******************************************************
