PROGRAM ThermWire

  IMPLICIT NONE
  
  !**************VAR. DECLARATION**************
  
  !*****************LEGEND*********************
  !C: NUMBER OF QUANTUM DOTS, N: MATRIX DIM., M: PH. ORDER 
  !WMAX: MAX. OMEGA, WINTERVAL: OMEGA INTERVAL             
  !ETA:COMPLEX PART OF OMEGA DDH:DOTDOT HOPPING			
  !MU_L,R:CHEMICALPOTENTIAL							
  !    W: OMEGA, T_L,R: LEAD'S TEMPERATURE     			
  ! LAMBDA=PH. COUPLING PHENERGY= PH. ENERGY  			
  !GR, GA, GLT, GGT: GREEN'S FUNCTIONS(R,A,<,>)			
  !A: COEFF. MATRIX, AUX: AUXILIAR MATRIX,TERM: INDEP. TERM 
  !SE,SELT,SEGT:SELF-ENERGY(R/A,<,>),T: HOPPING			    
  !SPECTRAL:: SPECTRAL FUNCTION						
  ! TEMPERATURE,ENERGY:: TEMP & ENERGY PROFILE 			
  !
  !UPDATE JAN-29TH: L_M STANDS FOR THE L_m PARAMETERS, 
  !WHICH DEPENDS ON THE BESSEL FUNCTIONS
  !********************************************

  INTEGER,PARAMETER::C=2,N=C*C,M=1
  INTEGER,PARAMETER::output_1=31, output_2=32, output_3=33, output_4=34, output_5=35, output_6=36
  DOUBLE PRECISION,PARAMETER:: &
       WMAX=3.0D0,WINTERVAL=0.005D0,Tmax=1.5D0,dTT=0.05D0,Vmax=0.0D0,dV=0.01D0,DDH=0.0D0,&
       MU_L=0.0D0,T_R=0.01D0,LAMBDA=0.5D0,b=1.0D0,&
       PHENERGY=1.00D0
  COMPLEX*16,PARAMETER::ETA=(0,1)*1.0D-6
  COMPLEX*16::GRR(C,C),GAA(C,C),GR(C,C),GA(C,C),GRLT(C,C),GRGT(C,C),&
       GLT(C,C),GGT(C,C),SPECTRAL(C,C),TRANSCOEFF(C,C),TERM(N),X(N),&
       WX,AUX(C,C),AUX2(C,C),AUX3(C,C)
  COMPLEX*16::A(N,N),SE(C,C),SELT(C,C),SEGT(C,C),T(C,C),GAMMAL(C,C),GAMMAR(C,C)
  DOUBLE PRECISION::TEMPERATURE(C),ENERGY(N),W,TT,T_L,V,MU_R,L_M(C,C),AUXSPECTRAL,AUXCOEFF,AUXAUX,AUXAUX2,AUXAUX3,JE,JQ

  INTEGER::I,J,K,L,Q,P,R,IX,IY,IZ,auxrow,auxcolumn

  OPEN(unit=output_1,file="greater.txt",action="write",status="replace")
  OPEN(unit=output_2,file="lesser.txt",action="write",status="replace")
  OPEN(unit=output_3,file="spectral.txt",action="write",status="replace")
  OPEN(unit=output_4,file="transmision.txt",action="write",status="replace")
  OPEN(unit=output_5,file="current.txt",action="write",status="replace")
  OPEN(unit=output_6,file="heatcurrent.txt",action="write",status="replace")

  !!**********END VAR. DECLARATION**************

  PRINT*,"ThermWire (c) 2015, Diaz/Dominguez-Adame/Roemer"

 
  T_L=0.0+T_R
 
  !DO WHILE(T_L<=Tmax)
  
     

  !FILL THE MATRICES AND TEMPERATURE AND ENERGY PROFILES!
  !!******************TEMPERATURE PROFILE****************!
  DO I=1,C
     TEMPERATURE(I)=T_L+(I)*(T_R-T_L)/C!-((T_R-T_L)/C)      !QUITE SIMPLE LINEAR PROFILE 
     !write(8,*),TEMPERATURE(I)
  END DO
!   TEMPERATURE(5)=0.90000000000000013D0
!   TEMPERATURE(4)=0.70000000000000007D0
!   TEMPERATURE(3)=0.50000000000000011D0
!   TEMPERATURE(2)=0.30000000000000004D0
!   TEMPERATURE(1)=0.10000000000000009D0
  !**********************HOPPING************************!
  DO I=1,C
     DO J=1,C
        T(I,J)=0
        IF (J==I+1) T(I,J)= DDH*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))&
             *PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I+1))
        IF (J==I-1) T(I,J)= DDH*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))&
             *PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I-1))
        !write(9,*),I,J,T(I,J)
     END DO
  END DO

  !*********************ENERGY PROFILE******************!
!!$  DO I=1,N
!!$     ENERGY(I)=0.0D0 !RANDOM ENERGIES MUST BE INCLUDED
!!$  END DO
! ENERGY=0.0D0
 ENERGY(1)=0.25D0-LAMBDA*LAMBDA/PHENERGY
 ENERGY(2)=0.25D0-LAMBDA*LAMBDA/PHENERGY
  !*********************GAMMA MATRICES*************!
  ! 
!!$  DO I=1,Cy
!!$     DO J=1,C
!!$        GAMMAL(I,J)=0.0D0
!!$        GAMMAR(I,J)=0.0D0
!!$     END DO
!!$  END DO
  GAMMAL=0.0D0
  GAMMAR=0.0D0

  ! PLA RESULTS
  GAMMAL(1,1)=0.2D0 
  GAMMAL(1,2)=0.2D0*SQRT(b)
  GAMMAL(2,1)=0.2D0*SQRT(b)
  GAMMAL(C,C)=0.2D0*b

  GAMMAR(1,1)=0.2D0*b 
  GAMMAR(1,2)=0.2D0*SQRT(b)
  GAMMAR(2,1)=0.2D0*SQRT(b)
  GAMMAR(C,C)=0.2D0

!   ! LINEAR CHAIN
!  GAMMAL(1,1)=0.2 
!  GAMMAR(C,C)=0.2
  !*****************************************************!

  !*****************************************************!
  !MAIN LOOP
  !*****************************************************!

  V=-Vmax
  DO WHILE(V<=Vmax)
     W=-WMAX
     !W=0.98D0	 

     MU_R=-V

     !WE CALCULATE FIRST THE GGT FUNCTION AND THEN THE GLT
     !THE TWO PARTS ARE BASICALLY THE SAME
     !WITH MINOR DIFFERENCES
     JE=0.0D0
     JQ=0.0D0
     DO WHILE(W<=WMAX)
        !*********************************GREATER THAN FUNCTION***

        !PRINT*,"ThermWire(): W=", W

        ! ! *********************************CHECKING SOLVER******
        ! DO I=1,N
        !         TERM(I)=0.0+(0,1)*0.0
        !         X(I)=0.0
        ! END DO
        ! TERM(4)=1.0+(0,1)*0.0
        ! PRINT*,TERM
        ! 
        ! DO I=1,N
        ! 	DO J=1,N
        ! 		A(I,J)=0.0+(0,1)*0.0
        ! 		IF (I==J) A(I,J)=1.0+(0,1)*0.0
        ! 	END DO
        ! 
        ! END DO
        ! PRINT*,A
        ! 
        ! CALL SOLVER(A,TERM,X) 
        ! DO I=1,N
        !         write(9,*),X(I)
        ! END DO
        ! PRINT*,X
        ! 
        ! ! !*********************************END CHECKING SOLVER*

!!$        DO I=1,C
!!$           DO J=1,C
!!$              GLT(I,J)=0.0D0
!!$              GGT(I,J)=0.0D0
!!$              GRR(I,J)=0.0D0
!!$              GAA(I,J)=0.0D0
!!$           END DO
!!$        END DO
        GLT=0.0D0
        GGT=0.0D0
        GRR=0.0D0
        GAA=0.0D0
       
        !Loop for W-PHENERGY*R to get GGT
        DO R=-M,M

!!$           DO I=1,C
!!$              DO J=1,C
!!$                 GR(I,J)=0.0D0
!!$                 GA(I,J)=0.0D0
!!$                 GRLT(I,J)=0.0D0
!!$                 GRGT(I,J)=0.0D0
!!$              END DO
!!$           END DO
           GR=0.0D0
           GA=0.0D0
           GRLT=0.0D0
           GRGT=0.0D0
           SE=0.0D0
           SEGT=0.0D0
           SELT=0.0D0
           L_M=0.0D0

           !First loop to set self energies and coeff L_M
           DO I=1,C
              DO J=1,C
                 SE(I,J)=-0.5*(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))&
                      *PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                 SE(I,J)=SE(I,J)*(GAMMAL(I,J)+GAMMAR(I,J))
                 !write(8,*),W,R,I,J,SE(I,J)

                 SEGT(I,J)=-(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))&
                      *PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                 SEGT(I,J)=SEGT(I,J)*((1-(FERMI(W-PHENERGY*R,MU_L,T_L)))&
                      *GAMMAL(I,J)+(1-FERMI(W-PHENERGY*R,MU_R,T_R))*GAMMAR(I,J))   
                 !write(8,*),W,R,I,J,SEGT(I,J)

                 SELT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))&
                      *PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                 SELT(I,J)=SELT(I,J)*(FERMI(W-PHENERGY*R,MU_L,T_L)&
                      *GAMMAL(I,J)+FERMI(W-PHENERGY*R,MU_R,T_R)*GAMMAR(I,J))
                 !write(8,*),W,R,I,J,SELT(I,J)

                 IF (I==J) THEN

                    L_M(I,J)=2*(LAMBDA/PHENERGY)**2*&
                         SQRT(BOSE(PHENERGY,TEMPERATURE(I))&
                         *(BOSE(PHENERGY,TEMPERATURE(I))+1))

                 !write(8,*),W,R,L_M(I,J)
!!$                         PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))**2,&
!!$                         EXP(R*PHENERGY/(TEMPERATURE(I)*2))
!!$                    write(8,*),W,R,BOSE(PHENERGY,TEMPERATURE(I)),&
!!$                         SQRT(BOSE(PHENERGY,TEMPERATURE(I))*&
!!$                         (BOSE(PHENERGY,TEMPERATURE(I))+1))
                  !  write(8,*),R,L_M(I,J),BESSEL(R,L_M(I,J)),BESSEL(R,L_M(I,J))&
                  !       *EXP(R*PHENERGY/(TEMPERATURE(I)*2))

                    !There is some problem here cause for R=1 we should obtain 
                    !something finite but it seems to be related with 
                    !higher decimals
                    L_M(I,J)=PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))**2&
                         *EXP(R*PHENERGY/(TEMPERATURE(I)*2))*BESSEL(R,L_M(I,J))

                 !write(8,*),R,L_M(I,J)
                         
!                    IF (R==-3) L_M(I,J)=0.0D0
!                    IF (R==-1) L_M(I,J)=0.0D0
!                     IF (R==-2) L_M(I,J)=0.0D0
!                     IF (R==0) L_M(I,J)=0.778801D0
!                     IF (R==1) L_M(I,J)=0.1947D0
!                     IF (R==2) L_M(I,J)=0.0243375D0
!                     IF (R==3) L_M(I,J)=0.00202813D0
!!$                    write(8,*),W,R,L_M(I,J),&
!!$                         PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))**2,&
!!$                         EXP(R*PHENERGY/(TEMPERATURE(I)*2))

                 ELSE
                    L_M(I,J)=PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*&
                         PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                   

                 END IF
                !write(8,*),W,R,I,J,L_M(I,J)
                         
              END DO
           END DO

           !DRESSED RETARDED FUNCTION
           WX=W+ETA-PHENERGY*R

           !FILL THE COEFF. MATRIX 
!!$           DO Q=1,N
!!$              DO L=1,N
!!$                 A(L,Q)=0.0d0
!!$              END DO
!!$           END DO
           A=0.0D0

           DO IX=1,C 
              DO Q=1,N
                 DO L=1,N

                    IF (L==Q) THEN 
                       A(L,Q)=1.0 !FILL THE DIAGONAL
                       !write(8,*),W-PHENERGY*R,R,L,Q,A(L,Q)

                    ELSE IF (Q==L+IX*C) THEN !FILL UPPER DIAGONALS
                       auxrow=INT((L - 1)/C) + 1
                       auxcolumn=INT((Q - 1)/C) + 1
                       A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))&
                            /(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
                       !write(8,*),W-PHENERGY*R,R,L,Q,auxrow,auxcolumn,A(L,Q)

                    ELSE IF (Q==L-IX*C) THEN !FILL LOWER DIAGONALS
                       auxcolumn=INT((Q - 1)/C) + 1
                       auxrow=INT((L - 1)/C) + 1
                       A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))&
                            /(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
                       !write(8,*),W-PHENERGY*R,R,L,Q,auxrow,auxcolumn,A(L,Q)

                    END IF
                 END DO
              END DO
           END DO
           DO L=1,N
                 DO Q=1,N
                 !write(8,*),L,Q,Real(A(L,Q)),Real((0,1)*A(L,Q))
                 END DO
           END DO
           !FILL THE INDEPENDENT TERM 
           DO L=1,N
              TERM(L)=0
              !write(9,*),I,J,L,TERM(L)
              IF (MOD(L,C+1)==1) THEN
                 P= INT(L/(C+1))+1
                 TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
              
              END IF
              !write(8,*),Real(TERM(L)),Real((0,1)*TERM(L))
           END DO

!!$           DO L=1,N 
!!$              X(L)=0.0D0
!!$           END DO
           x=0.0D0

           !CALL THE LINEAR SYSTEM SOLVER TO GET DRESSED GR
           CALL SOLVER(A,TERM,X) 

           DO L=1,N
              auxrow=INT((L - 1)/C) + 1
              auxcolumn=MOD((L - 1),C) + 1
              GR(auxrow,auxcolumn)=X(L) 
              !write(8,*),L,auxrow,auxcolumn
           END DO
           !END IF

           !*****************************************************!
           !NOW OBTAIN THE RETARDED DRESSED GREEN FUNCTION
           ! DO I=1,C
           ! 	DO J=1,C
           !            IF(I==J) THEN
           !               GR(I,J)=GR(I,J)
           !               !write(9,*),I,J,GR(I,J)
           !            ELSE
           !               GR(I,J)=AUX(I,J)
           !               !write(9,*),I,J,GR(I,J)
           !            END IF
           !        END DO
           ! END DO
           !Only nondiagonal elements with R=0 are relevant

           !DRESSED ADVANCED FUNCTION
           DO I=1,C
              DO J=1,C
                 SE(I,J)=-SE(I,J)
              END DO
           END DO

           WX=W-ETA-PHENERGY*R

           !FILL THE COEFF. MATRIX 
!!$           DO Q=1,N
!!$              DO L=1,N
!!$                 A(L,Q)=0.0D0
!!$              END DO
!!$           END DO
           A=0.0D0
           
           DO IX=1,C 
              DO Q=1,N
                 DO L=1,N

                    IF (L==Q) THEN 
                       A(L,Q)=1.0 !FILL THE DIAGONAL
                       !write(8,*),I,J,L,Q,A(L,Q)

                    ELSE IF (Q==L+IX*C) THEN !FILL UPPER DIAGONALS
                       auxrow=INT((L - 1)/C) + 1
                       auxcolumn=INT((Q - 1)/C) + 1
                       A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))&
                            /(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
                       !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)

                    ELSE IF (Q==L-IX*C) THEN !FILL LOWER DIAGONALS
                       auxcolumn=INT((Q - 1)/C) + 1
                       auxrow=INT((L - 1)/C) + 1
                       A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))&
                            /(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
                       !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)

                    END IF

                 END DO
              END DO
           END DO

           !FILL THE INDEPENDENT TERM 
           DO L=1,N
              TERM(L)=0
              !write(9,*),I,J,L,TERM(L)
              IF (MOD(L,C+1)==1) THEN
                 P= INT(L/(C+1))+1
                 TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
                 !write(9,*),L,P,TERM(L)
              END IF

           END DO

!!$           DO L=1,N 
!!$              X(L)=0.0D0
!!$           END DO
           X=0.0D0

           !CALL THE LINEAR SYSTEM SOLVER TO GET DRESSED GR
           CALL SOLVER(A,TERM,X) 
           !                                 
           !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET DRESSED GR IN X
           !HERE, WE WANT JUST THE DIAGONAL TERMS (THOSE OF GR(1,1)...)
           !                                
           !*****THIS SAVES THE DESIRED TERMS IN THE DIAGONAL*****!
           ! 				

           DO L=1,N
              auxrow=INT((L - 1)/C) + 1
              auxcolumn=MOD((L - 1),C) + 1
              GA(auxrow,auxcolumn)=X(L) 
              !write(8,*),L,auxrow,auxcolumn
           END DO

           !*****************************************************!
           ! 			
           !NOW OBTAIN THE ADVANCED GREEN FUNCTION BUT ALREADY 
           !MULTIPLIED BY THE COEFF L_M
           ! DO I=1,C
           ! 	DO J=1,C
           !            IF(I==J) THEN
           !               GA(I,J)=GA(I,J)
           !               
           !            ELSE
           !               GA(I,J)=AUX(I,J)
           !               
           !            END IF
           !        END DO
           ! END DO
           !Only nondiagonal elements with R=0 are relevant
           !write(8,*),W-PHENERGY*R,R,GA(2,1),GR(2,1)
           !write(9,*),W-PHENERGY*R,R,SEGT(1,2),SELT(1,2)

           !Here we get GRGT at w-R*w0
           GRGT=MATMUL(GR,MATMUL(SEGT,GA))
           GRLT=MATMUL(GR,MATMUL(SELT,GA))
!           DO I=1,C
!               DO J=1,C
!                 write(81,*),R,I,J,Real(GRGT(I,J)),Real((0,1)*GRGT(I,J))
!                 write(91,*),R,I,J,Real(GRLT(I,J)),Real((0,1)*GRLT(I,J))
!                END DO
!           END DO
           DO I=1,C
              DO J=1,C
                 IF(I==J) THEN
                    GGT(I,J)=GGT(I,J)+L_M(I,J)*GRGT(I,J)
                    !write(8,*),W-PHENERGY*R,R,GRGT(1,1),GRGT(1,2)
                 ELSE 
                    IF(R==0) THEN
                       GGT(I,J)=L_M(I,J)*GRGT(I,J)
                       !write(8,*),W-PHENERGY*R,R,GRGT(1,1),GRGT(1,2)
                    END IF
                 END IF
                 !write(81,*),R,I,J,Real(GGT(I,J)),Real((0,1)*GGT(I,J))
              END DO
           END DO
           !write(88,*),W-PHENERGY*R,W-PHENERGY*R

           DO I=1,C
              DO J=1,C
                     
                  IF(I==J) THEN
                     GRR(I,J)=GRR(I,J)+L_M(I,J)*(GR(I,J)+0.5*GRLT(I,J))
                     !GAA(I,J)=GAA(I,J)+L_M(I,J)*(GR(I,J)+0.5*GRLT(I,J)-GRGT(I,J))
                  ELSE 
                     IF(R==0) THEN
                        GRR(I,J)=L_M(I,J)*(GR(I,J)+GRLT(I,J)-GRGT(I,J))
                        !GAA(I,J)=L_M(I,J)*(GA(I,J)
                        !write(81,*),W,I,J,Real(GRR(I,J)),Real((0,1)*GRR(I,J))
                        !write(91,*),W,I,J,Real(GAA(I,J)),Real((0,1)*GAA(I,J))
                     END IF
                  END IF
                  END DO
                 
           END DO

        END DO !Loop in R


        !Loop for W-PHENERGY*R to get GLT
        DO R=-M,M

!!$           DO I=1,C
!!$              DO J=1,C
!!$                 GR(I,J)=0.0D0
!!$                 GA(I,J)=0.0D0
!!$                 GRLT(I,J)=0.0D0
!!$                 GRGT(I,J)=0.0D0
!!$              END DO
!!$           END DO
           GR=0.0D0
           GA=0.0D0
           GRLT=0.0D0
           GRGT=0.0D0
           SE=0.0D0
           SEGT=0.0D0
           SELT=0.0D0
           L_M=0.0D0
           !First loop to set self energies and coeff L_M

           DO I=1,C
              DO J=1,C
                 SE(I,J)=-0.5*(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))&
                      *PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                 SE(I,J)=SE(I,J)*(GAMMAL(I,J)+GAMMAR(I,J))
                 !write(8,*),W,R,I,J,SE(I,J)

                 SEGT(I,J)=-(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))&
                      *PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                 SEGT(I,J)=SEGT(I,J)*((1-FERMI(W+PHENERGY*R,MU_L,T_L))&
                      *GAMMAL(I,J)+(1-FERMI(W+PHENERGY*R,MU_R,T_R))*GAMMAR(I,J))   
                 !write(8,*),W,R,I,J,SEGT(I,J)

                 SELT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))&
                      *PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                 SELT(I,J)=SELT(I,J)*(FERMI(W+PHENERGY*R,MU_L,T_L)&
                      *GAMMAL(I,J)+FERMI(W+PHENERGY*R,MU_R,T_R)*GAMMAR(I,J))
                 !write(8,*),W,R,I,J,SELT(I,J)

                 IF (I==J) THEN
                     L_M(I,J)=2*(LAMBDA/PHENERGY)**2*&
                         SQRT(BOSE(PHENERGY,TEMPERATURE(I))&
                         *(BOSE(PHENERGY,TEMPERATURE(I))+1))

                 !write(8,*),W,R,L_M(I,J)
!!$                         PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))**2,&
!!$                         EXP(R*PHENERGY/(TEMPERATURE(I)*2))
!!$                    write(8,*),W,R,BOSE(PHENERGY,TEMPERATURE(I)),&
!!$                         SQRT(BOSE(PHENERGY,TEMPERATURE(I))*&
!!$                         (BOSE(PHENERGY,TEMPERATURE(I))+1))
                  !  write(8,*),R,L_M(I,J),BESSEL(R,L_M(I,J)),BESSEL(R,L_M(I,J))&
                  !       *EXP(R*PHENERGY/(TEMPERATURE(I)*2))

                    !There is some problem here cause for R=1 we should obtain 
                    !something finite but it seems to be related with 
                    !higher decimals
                    L_M(I,J)=PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))**2&
                         *EXP(R*PHENERGY/(TEMPERATURE(I)*2))*BESSEL(R,L_M(I,J))
                    
          

                 ELSE

                    L_M(I,J)=PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))&
                         *PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                    !write(8,*),W,R,I,J,SEGT(I,J),SELT(I,J),L_M(I,J)
                   
                 END IF

              END DO
           END DO

           !DRESSED RETARDED FUNCTION
           WX=W+ETA+PHENERGY*R

           !FILL THE COEFF. MATRIX 
!!$           DO Q=1,N
!!$              DO L=1,N
!!$                 A(L,Q)=0.0D0
!!$              END DO
!!$           END DO
           A=0.0D0

           DO IX=1,C 
              DO Q=1,N
                 DO L=1,N

                    IF (L==Q) THEN 
                       A(L,Q)=1.0 !FILL THE DIAGONAL
                       !write(8,*),I,J,L,Q,A(L,Q)

                    ELSE IF (Q==L+IX*C) THEN !FILL UPPER DIAGONALS
                       auxrow=INT((L - 1)/C) + 1
                       auxcolumn=INT((Q - 1)/C) + 1
                       A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))&
                            /(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
                       !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)

                    ELSE IF (Q==L-IX*C) THEN !FILL LOWER DIAGONALS
                       auxcolumn=INT((Q - 1)/C) + 1
                       auxrow=INT((L - 1)/C) + 1
                       A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))&
                            /(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
                       !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)

                    END IF

                 END DO
              END DO
           END DO

           !FILL THE INDEPENDENT TERM 
           DO L=1,N
              TERM(L)=0
              !write(9,*),I,J,L,TERM(L)
              IF (MOD(L,C+1)==1) THEN
                 P= INT(L/(C+1))+1
                 TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
                 !write(9,*),L,P,TERM(L)
              END IF
           END DO

!!$           DO L=1,N 
!!$              X(L)=0.0D0
!!$           END DO
           X=0.0D0

           !CALL THE LINEAR SYSTEM SOLVER TO GET DRESSED GR
           CALL SOLVER(A,TERM,X) 

           !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET DRESSED GR IN X
           !HERE, WE WANT JUST THE DIAGONAL TERMS (THOSE OF GR(1,1)...)

           !*****THIS SAVES THE DESIRED TERMS IN THE DIAGONAL*****!
           DO L=1,N
              auxrow=INT((L - 1)/C) + 1
              auxcolumn=MOD((L - 1),C) + 1
              GR(auxrow,auxcolumn)=X(L) 
              !write(8,*),L,auxrow,auxcolumn
           END DO
           !END IF
           !*****************************************************!

           !NOW OBTAIN THE RETARDED DRESSED GREEN FUNCTION
           ! DO I=1,C
           ! 	DO J=1,C
           !            IF(I==J) THEN
           !               GR(I,J)=GR(I,J)
           !               !write(9,*),I,J,GR(I,J)
           !            ELSE
           !               GR(I,J)=AUX(I,J)
           !               !write(9,*),I,J,GR(I,J)
           !            END IF
           !        END DO
           ! END DO
           !write(8,*),W+PHENERGY*R,GR(1,1),GR(2,2)

           !DRESSED ADVANCED FUNCTION
           DO I=1,C
              DO J=1,C
                 SE(I,J)=-SE(I,J)
              END DO
           END DO
           ! 

           WX=W-ETA+PHENERGY*R

           !FILL THE COEFF. MATRIX 
!!$           DO Q=1,N
!!$              DO L=1,N
!!$                 A(L,Q)=0.0D0
!!$              END DO
!!$           END DO
           A=0.0D0

           DO IX=1,C 
              DO Q=1,N
                 DO L=1,N

                    IF (L==Q) THEN 
                       A(L,Q)=1.0 !FILL THE DIAGONAL
                       !write(8,*),I,J,L,Q,A(L,Q)

                    ELSE IF (Q==L+IX*C) THEN !FILL UPPER DIAGONALS
                       auxrow=INT((L - 1)/C) + 1
                       auxcolumn=INT((Q - 1)/C) + 1
                       A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))&
                            /(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
                       !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)

                    ELSE IF (Q==L-IX*C) THEN !FILL LOWER DIAGONALS
                       auxcolumn=INT((Q - 1)/C) + 1
                       auxrow=INT((L - 1)/C) + 1
                       A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))&
                            /(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
                       !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)

                    END IF


                 END DO
              END DO

           END DO

           !FILL THE INDEPENDENT TERM 
           DO L=1,N
              TERM(L)=0
              !write(9,*),I,J,L,TERM(L)
              IF (MOD(L,C+1)==1) THEN
                 P= INT(L/(C+1))+1
                 TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
                 !write(9,*),L,P,TERM(L)
              END IF

           END DO

!!$           DO L=1,N 
!!$              X(L)=0.0D0
!!$           END DO
           X=0.0D0

           !CALL THE LINEAR SYSTEM SOLVER TO GET DRESSED GR
           CALL SOLVER(A,TERM,X) 

           !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET DRESSED GR IN X
           !HERE, WE WANT JUST THE DIAGONAL TERMS (THOSE OF GR(1,1)...)

           !*****THIS SAVES THE DESIRED TERMS IN THE DIAGONAL*****!

           DO L=1,N
              auxrow=INT((L - 1)/C) + 1
              auxcolumn=MOD((L - 1),C) + 1
              GA(auxrow,auxcolumn)=X(L) 
              !write(8,*),L,auxrow,auxcolumn
           END DO

           !*****************************************************!

      
           GRGT=MATMUL(GR,MATMUL(SEGT,GA))
           GRLT=MATMUL(GR,MATMUL(SELT,GA))

           DO I=1,C
              DO J=1,C
                 IF(I==J) THEN
                    GLT(I,J)=GLT(I,J)+L_M(I,J)*GRLT(I,J)
                 ELSE 
                    IF(R==0) THEN
                       GLT(I,J)=L_M(I,J)*GRLT(I,J)
                    END IF
                 END IF 
                 !write(91,*),R,I,J,Real(GLT(I,J)),Real((0,1)*GLT(I,J))
              END DO
           END DO

           DO I=1,C
              DO J=1,C
                   IF(I==J) THEN
                    GRR(I,J)=GRR(I,J)-L_M(I,J)*(0.5*GRGT(I,J))
                    !GAA(I,J)=GAA(I,J)+L_M(I,J)*(0.5*GRLT(I,J))
                   !ELSE 
                   !IF(R==0) THEN
                      !GRR(I,J)=L_M(I,J)*GR(I,J)
                      !GAA(I,J)=L_M(I,J)*GA(I,J)
                      !write(81,*),W,I,J,Real(GRR(I,J)),Real((0,1)*GRR(I,J))
                      !write(91,*),W,I,J,Real(GAA(I,J)),Real((0,1)*GAA(I,J))
                   !END IF
                 END IF
              END DO
           END DO
           !write(81,*),W,I,J,Real(GRR(1,1)),Real((0,1)*GRR(1,1))
           !write(91,*),W,I,J,Real(GAA(1,1)),Real((0,1)*GAA(1,1))

        END DO !Loop in R

        GAA=GRR+GLT-GGT          
        !SPECTRAL=(0,1)*(GGT-GLT)
        !SPECTRAL=(0,1)*(GRR-GAA)
        TRANSCOEFF=MATMUL(MATMUL(GLT,GAMMAR),MATMUL(GGT,GAMMAL))
        SPECTRAL=MATMUL(GGT,GAMMAL)
        !TRANSCOEFF=GAA
        
        write(81,*),W,GLT(1,1),GLT(1,1)
        write(91,*),W,GAA(1,1),GAA(1,1)
     
        
        AUX=MATMUL(GAMMAL-GAMMAR,(0,1)*GLT)+MATMUL(GAMMAL*FERMI(W,MU_L,T_L)-GAMMAR*FERMI(W,MU_R,T_R),SPECTRAL)        
        AUX2=MATMUL(GAMMAL,(0,1)*GLT)+MATMUL(GAMMAL*FERMI(W,MU_L,T_L),SPECTRAL)
        AUX3=MATMUL(GAMMAR,(0,1)*GLT)+MATMUL(GAMMAR*FERMI(W,MU_R,T_L),SPECTRAL)
        
        AUXSPECTRAL=0.0
        AUXCOEFF=0.0
        AUXAUX=0.0
        AUXAUX2=0.0
        AUXAUX3=0.0
        DO J=1,C
                AUXSPECTRAL=AUXSPECTRAL+Real(SPECTRAL(J,J))
                AUXCOEFF= AUXCOEFF+Real(TRANSCOEFF(J,J))
                AUXAUX= AUXAUX+AUX(J,J)
                AUXAUX2= AUXAUX2+AUX2(J,J)
                AUXAUX3= AUXAUX3+AUX3(J,J)
        END DO
        JE=JE+WINTERVAL*AUXAUX
        JQ=JQ+WINTERVAL*( (W-MU_L) * AUXAUX2 - (W-MU_R)*AUXAUX3  )
         
       

        
       WRITE(output_3,*) W,AUXSPECTRAL/C
       WRITE(output_4,*) W,AUXCOEFF

        
        W=W+WINTERVAL
        
     END DO

     WRITE(output_5,*) V,JE 
     WRITE(output_6,*) T_L-T_R,JQ
     !WRITE(output_6,*) V,JQ
     PRINT*,"V", V
     PRINT*,"TL", T_L-T_R
     V=V+dV
     !T_L=T_L+dTT
  !END DO
  END DO
  CLOSE(output_1)
  CLOSE(output_2)
  CLOSE(output_3)
  CLOSE(output_4)  
  CLOSE(output_5)
  CLOSE(output_6)
  PRINT*,"--- finished!"   

  ! PRINT*,"x=1.5,n=-1,jn=", BESSEL(-1,1.5)
  ! PRINT*,"x=1.5,n=1,jn=", BESSEL(1,1.5)
  ! PRINT*,"x=1.5,n=-2,jn=", BESSEL(-2,1.5)
  ! PRINT*,"x=1.5,n=2,jn=", BESSEL(2,1.5)
  ! PRINT*,"x=1.5,n=-3,jn=", BESSEL(-3,1.5)

  STOP

CONTAINS

  DOUBLE PRECISION FUNCTION FERMI(FREC,MU,T)
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN)::FREC,MU,T
    DOUBLE PRECISION::OUTPUT
    IF(T==0) THEN
       IF(FREC>MU) OUTPUT=0
       IF(FREC<=MU) OUTPUT=1
    ELSE
       OUTPUT=1.D0/(1+EXP((FREC-MU)/T))
    END IF
    FERMI=OUTPUT
  END FUNCTION FERMI

  DOUBLE PRECISION FUNCTION BOSE(FREC,T)
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN)::FREC,T
    DOUBLE PRECISION::OUTPUT
    IF(T==0) THEN
       OUTPUT=0.D0
    ELSE
       OUTPUT=1.D0/(-1.D0+EXP(FREC/T))
    END IF

    BOSE=OUTPUT
  END FUNCTION BOSE

  DOUBLE PRECISION FUNCTION XHI(FREC,T)
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN)::FREC,T
    DOUBLE PRECISION::OUTPUT
    OUTPUT=((LAMBDA/PHENERGY)**2)*(2.D0*BOSE(FREC,T)+1.D0)

    XHI=OUTPUT
  END FUNCTION XHI

  DOUBLE PRECISION FUNCTION PHONON_PART(LAMBDA,PHENERGY,T)
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN)::LAMBDA,PHENERGY,T
    DOUBLE PRECISION::OUTPUT
    OUTPUT=((LAMBDA/PHENERGY)**2)*(2.D0*BOSE(PHENERGY,T)+1.D0)/2
    OUTPUT=EXP(-OUTPUT)
    PHONON_PART=OUTPUT
  END FUNCTION PHONON_PART

  INTEGER FUNCTION FACTORIAL(N)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I,OUTPUT
    OUTPUT=1
    DO I=1,N
       OUTPUT=OUTPUT*I
    END DO
    FACTORIAL=OUTPUT
  END FUNCTION FACTORIAL

  DOUBLE PRECISION FUNCTION TRACE(N,A)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    DOUBLE PRECISION::A(N,N)
    DOUBLE PRECISION::OUTPUT
    INTEGER::I,J
    DO I=1,N
       DO J=1,N
          IF (I==J) OUTPUT=OUTPUT+A(I,J)
       END DO
    END DO
    TRACE=OUTPUT
  END FUNCTION TRACE
  
  
  DOUBLE PRECISION FUNCTION BESSEL(N,X)
    IMPLICIT NONE
    INTEGER::N,N1,K,P,M
    DOUBLE PRECISION::X,OUTPUT,BK,BESSJ,BESSJ0,BESSJ1,BESSI,BESSI0,BESSI1
    N1=ABS(N)

    SELECT CASE(N1)
    CASE(0)
       BESSEL= BESSI0(X)
    CASE(1)
       BESSEL= BESSI1(X)
    CASE DEFAULT
       BESSEL= BESSI(N1,X)
    END SELECT

  END FUNCTION BESSEL

  ! DOUBLE PRECISION FUNCTION BESSEL(N,X)
  ! 	IMPLICIT NONE
  ! 	INTEGER::N,N1
  ! 	DOUBLE PRECISION::X
  ! 	DOUBLE PRECISION::OUTPUT,BK
  ! 	INTEGER::K,P,M
  ! 		N1=ABS(N)
  ! 		M=6
  ! 		OUTPUT=0
  ! 		P=1
  ! 		IF(N<0) P=(-1)**N
  ! 		DO K=0,M
  ! 			BK=1./((2**(2*K+N1))*FACTORIAL(K)*FACTORIAL(K+N1))
  ! 			OUTPUT=OUTPUT+BK*((X)**(2*K+N1))
  ! 		END DO
  ! 		OUTPUT=OUTPUT*P
  ! 	BESSEL=ABS(OUTPUT)
  ! END FUNCTION BESSEL

  SUBROUTINE SOLVER(MAT, RHS, VEC)
    COMPLEX*16::MAT(N,N), RHS(N), VEC(N)
    INTEGER::INFO
    INTEGER::IPIV(N)

    !PRINT*,"SOLVER()"

    !PRINT*,"before: mat=", MAT, "rhs=", RHS

    CALL ZGESV( N, 1, MAT, N, IPIV, RHS, N, INFO )

    IF(INFO.NE.0) THEN
       PRINT*,"SOLVER(): INFO=", INFO," --- WRNG!"
    ENDIF

    !PRINT*,"after: mat=", MAT, "rhs=", RHS

    VEC= RHS

  END SUBROUTINE SOLVER



END PROGRAM ThermWire
