PROGRAM linearchain
IMPLICIT NONE
!!!**************VAR. DECLARATION**************!!!
!!!********************************************!!!
!!!*****************LEGEND*********************!!!
!!!C: NUMBER OF QUANTUM DOTS, N: MATRIX DIM., M: PH. ORDER !!!
!!!WMAX: MAX. OMEGA, WINTERVAL: OMEGA INTERVAL             !!!
!!!ETA:COMPLEX PART OF OMEGA DDH:DOTDOT HOPPING			!!!
!!!MU_L,R:CHEMICALPOTENTIAL							!!!
!!!    W: OMEGA, T_L,R: LEAD'S TEMPERATURE     			!!!
!!! LAMBDA=PH. COUPLING PHENERGY= PH. ENERGY  			!!!
!!!GR, GA, GLT, GGT: GREEN'S FUNCTIONS(R,A,<,>)			!!!
!!!A: COEFF. MATRIX, AUX: AUXILIAR MATRIX,TERM: INDEP. TERM !!!
!!!SE,SELT,SEGT:SELF-ENERGY(R/A,<,>),T: HOPPING			!!!    
!!!SPECTRAL:: SPECTRAL FUNCTION						!!!
!!! TEMPERATURE,ENERGY:: TEMP & ENERGY PROFILE 			!!!
!!!UPDATE JAN-29TH: L_M STANDS FOR THE L_m PARAMETERS, WICH DEPENDS ON THE BESSEL FUNCTIONS!!!
INTEGER,PARAMETER::C=2,N=C*C,M=1
INTEGER,PARAMETER::output_1=1, output_2=2, output_3=3, output_4=4, output_5=5
REAL,PARAMETER:: WMAX=2.5,WINTERVAL=0.01,Vmax=2.0,dV=0.01,DDH=0.0,&
     MU_L=0.0,MU_R=0.0,T_L=0.01,T_R=0.01,LAMBDA=0.5,b=0.5,PHENERGY=1.00
COMPLEX*16,PARAMETER::ETA=(0,1)*1.E-6
COMPLEX*16::GRR(C,C),GAA(C,C),GR(C,C),GA(C,C),GRLT(C,C),GRGT(C,C),&
     GLT(C,C),GGT(C,C),SPECTRAL(C,C),TRANSCOEFF(C,C),TERM(N),X(N),&
     WX,AUX(C,C),AUX2(C,C)
COMPLEX*16::A(N,N),SE(C,C),SELT(C,C),SEGT(C,C),T(C,C),GAMMAL(C,C),GAMMAR(C,C)
REAL::TEMPERATURE(C),ENERGY(N),W,V,L_M(C,C),JE,JQ

INTEGER::I,J,K,L,Q,P,R,IX,IY,IZ,auxrow,auxcolumn

OPEN(unit=output_1,file="greater.txt",action="write",status="replace")
OPEN(unit=output_2,file="lesser.txt",action="write",status="replace")
OPEN(unit=output_3,file="spectral.txt",action="write",status="replace")
OPEN(unit=output_4,file="transmision.txt",action="write",status="replace")
OPEN(unit=output_5,file="current.txt",action="write",status="replace")

  !DOUBLE PRECISION jn
  !EXTERNAL jn

  PRINT*,"ThermWire"

!!**********END VAR. DECLARATION**************!!!

!FILL THE MATRICES AND TEMPERATURE AND ENERGY PROFILES!
!!******************TEMPERATURE PROFILE****************!
DO I=1,C
	TEMPERATURE(I)=T_L+I*(T_R-T_L)/C      !QUITE SIMPLE LINEAR PROFILE 
	!write(9,*),TEMPERATURE(I),T_L
END DO
!**********************HOPPING************************!
DO I=1,C
	DO J=1,C
		T(I,J)=0
		IF (J==I+1) T(I,J)= DDH*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I+1))
		IF (J==I-1) T(I,J)= DDH*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I-1))
		!write(9,*),I,J,T(I,J)
	END DO

END DO
!*********************ENERGY PROFILE******************!
DO I=1,N
        ENERGY(I)=0.0 !RANDOM ENERGIES MUST BE INCLUDED
END DO
ENERGY(1)=0.25-LAMBDA*LAMBDA/PHENERGY
ENERGY(2)=-0.25-LAMBDA*LAMBDA/PHENERGY
!*********************GAMMA MATRICES*************!
! 
DO I=1,C
	DO J=1,C
		GAMMAL(I,J)=0.0
		GAMMAR(I,J)=0.0
	END DO

END DO
! PLA RESULTS
 GAMMAL(1,1)=0.2 
 GAMMAL(1,2)=0.2*SQRT(b)
 GAMMAL(2,1)=0.2*SQRT(b)
 GAMMAL(C,C)=0.2*b

 GAMMAR(1,1)=0.2*b 
 GAMMAR(1,2)=0.2*SQRT(b)
 GAMMAR(2,1)=0.2*SQRT(b)
 GAMMAR(C,C)=0.2

! ! LINEAR CHAIN
! GAMMAL(1,1)=0.2 
! GAMMAR(C,C)=0.2
! !*****************************************************!
! 

V=-Vmax
DO WHILE(V<=-Vmax)
W=-WMAX	     !HERE WE START THE MAIN LOOP
			!WE CALCULATE FIRST THE GGT FUNCTION AND THEN THE GLT
			!THE TWO PARTS ARE BASICALLY THE SAME
			!WITH MINOR DIFFERENCES
JE=0.0
DO WHILE(W<=WMAX)
!*********************************GREATER THAN FUNCTION************************************************** 

!PRINT*,"ThermWire(): W=", W

! ! *********************************CHECKING SOLVER************************************************** 
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
! ! !*********************************END CHECKING SOLVER************************************************** 

DO I=1,C
	DO J=1,C
	GLT(I,J)=0.0
	GGT(I,J)=0.0
	GRR(I,J)=0.0
	GAA(I,J)=0.0
        END DO
END DO

!Loop for W-PHENERGY*R to get GGT
DO R=-M,M

DO I=1,C
	DO J=1,C
	GR(I,J)=0.0
	GA(I,J)=0.0
       	GRLT(I,J)=0.0
	GRGT(I,J)=0.0
	END DO
END DO

!First loop to set self energies and coeff L_M
DO I=1,C
	DO J=1,C
		SE(I,J)=-0.5*(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
  		SE(I,J)=SE(I,J)*(GAMMAL(I,J)+GAMMAR(I,J))
  		!write(8,*),W,R,I,J,SE(I,J)

                SEGT(I,J)=-(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
		SEGT(I,J)=SEGT(I,J)*((1-(FERMI(W-PHENERGY*R,MU_L,T_L)))*GAMMAL(I,J)+(1-FERMI(W-PHENERGY*R,MU_R,T_R))*GAMMAR(I,J))   
		!write(8,*),W,R,I,J,SEGT(I,J)
		
                SELT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                SELT(I,J)=SELT(I,J)*(FERMI(W-PHENERGY*R,MU_L,T_L)*GAMMAL(I,J)+FERMI(W-PHENERGY*R,MU_R,T_R)*GAMMAR(I,J))
                !write(8,*),W,R,I,J,SELT(I,J)
                
		IF (I==J) THEN
		    
                L_M(I,J)=2*(LAMBDA/PHENERGY)**2*SQRT(BOSE(PHENERGY,TEMPERATURE(I))*(BOSE(PHENERGY,TEMPERATURE(I))+1))
                 
                !write(8,*),W,R,L_M(I,J),PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))**2,EXP(R*PHENERGY/(TEMPERATURE(I)*2))
                write(8,*),W,R,BOSE(PHENERGY,TEMPERATURE(I)),SQRT(BOSE(PHENERGY,TEMPERATURE(I))*(BOSE(PHENERGY,TEMPERATURE(I))+1))
                write(8,*),BESSEL(R,L_M(I,J)),BESSEL(R,L_M(I,J))*EXP(R*PHENERGY/(TEMPERATURE(I)*2))
                
                !There is some problem here cause for R=1 we should obtain something finite but it seems to 
                !be related with higher decimals
                L_M(I,J)=PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))**2*EXP(R*PHENERGY/(TEMPERATURE(I)*2))*BESSEL(R,L_M(I,J))
                
                IF (R==-3) L_M(I,J)=0.0
                IF (R==-1) L_M(I,J)=0.0
                IF (R==-2) L_M(I,J)=0.0
                IF (R==0) L_M(I,J)=0.778801
                IF (R==1) L_M(I,J)=0.1947
                IF (R==2) L_M(I,J)=0.0243375
                IF (R==3) L_M(I,J)=0.00202813
                !write(8,*),W,R,L_M(I,J),PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))**2,EXP(R*PHENERGY/(TEMPERATURE(I)*2))
                                
                ELSE
                L_M(I,J)=PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                !write(8,*),W,R,I,J,SEGT(I,J)
                !write(8,*),W,R,I,J,SELT(I,J)
                
             END IF
             
	END DO

END DO
              
              !DRESSED RETARDED FUNCTION
              WX=W+ETA-PHENERGY*R
		
				
				!FILL THE COEFF. MATRIX 
 	                        DO Q=1,N
				            DO L=1,N
				            A(L,Q)=0.0
                                            END DO
					
				END DO
				
				DO IX=1,C 
	                 		DO Q=1,N
				                DO L=1,N
						
						IF (L==Q) THEN 
						A(L,Q)=1.0 !FILL THE DIAGONAL
						!write(8,*),W-PHENERGY*R,R,L,Q,A(L,Q)
						
						ELSE IF (Q==L+IX*C) THEN !FILL UPPER DIAGONALS
						auxrow=INT((L - 1)/C) + 1
						auxcolumn=INT((Q - 1)/C) + 1
						A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))/(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
					        !write(8,*),W-PHENERGY*R,R,L,Q,auxrow,auxcolumn,A(L,Q)
					        
						ELSE IF (Q==L-IX*C) THEN !FILL LOWER DIAGONALS
						auxcolumn=INT((Q - 1)/C) + 1
						auxrow=INT((L - 1)/C) + 1
						A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))/(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
					        !write(8,*),W-PHENERGY*R,R,L,Q,auxrow,auxcolumn,A(L,Q)
						
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
				
				DO L=1,N 
					X(L)=0.0
				END DO
				
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

! 				
! 				!FILL THE COEFF. MATRIX 
 	                        DO Q=1,N
				            DO L=1,N
				            A(L,Q)=0.0
                                            END DO
					
				END DO
				
				DO IX=1,C 
	                 		DO Q=1,N
				                DO L=1,N
						
						IF (L==Q) THEN 
						A(L,Q)=1.0 !FILL THE DIAGONAL
						!write(8,*),I,J,L,Q,A(L,Q)
						
						ELSE IF (Q==L+IX*C) THEN !FILL UPPER DIAGONALS
						auxrow=INT((L - 1)/C) + 1
						auxcolumn=INT((Q - 1)/C) + 1
						A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))/(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
					        !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)
					        
						ELSE IF (Q==L-IX*C) THEN !FILL LOWER DIAGONALS
						auxcolumn=INT((Q - 1)/C) + 1
						auxrow=INT((L - 1)/C) + 1
						A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))/(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
					        !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)
						
						END IF
						
						
						END DO
					END DO
					
				END DO
				
			
! 
!                                 !FILL THE INDEPENDENT TERM 
                                DO L=1,N
                                        TERM(L)=0
                                        !write(9,*),I,J,L,TERM(L)
                                        IF (MOD(L,C+1)==1) THEN
                                                           P= INT(L/(C+1))+1
                                                           TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
                                                           !write(9,*),L,P,TERM(L)
                                        END IF
                 
                                END DO	
                                 

! 				
				DO L=1,N 
					X(L)=0.0
				END DO
! 				
! 				!CALL THE LINEAR SYSTEM SOLVER TO GET DRESSED GR
! 				
                                CALL SOLVER(A,TERM,X) 
!                                 
!                                 !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET DRESSED GR IN X
! 				!HERE, WE WANT JUST THE DIAGONAL TERMS (THOSE OF GR(1,1)...)
!                                
! 				!*****THIS SAVES THE DESIRED TERMS IN THE DIAGONAL*****!
! 				
				DO L=1,N
                                        
                                        IF (MOD(L,C+1)==1) THEN
                                                           P= INT(L/(C+1))+1
                                                           GA(P,P)=X(L)
                                                           !write(9,*),L,P,GA(P,P)
                                        END IF
                 
                                END DO	


				DO L=1,N
					auxrow=INT((L - 1)/C) + 1
					auxcolumn=MOD((L - 1),C) + 1
                                        GA(auxrow,auxcolumn)=X(L) 
                                        !write(8,*),L,auxrow,auxcolumn
				END DO

! 				!*****************************************************!
! 			

! ! !NOW OBTAIN THE ADVANCED GREEN FUNCTION BUT ALREADY MULTIPLIED BY THE COEFF L_M
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
       END DO
END DO
write(88,*),W-PHENERGY*R,W-PHENERGY*R

DO I=1,C
	DO J=1,C
	IF(I==J) THEN
          GRR(I,J)=GRR(I,J)+L_M(I,J)*(GR(I,J)+0.5*GRLT(I,J))
          GAA(I,J)=GAA(I,J)+L_M(I,J)*(GR(I,J)+0.5*GRLT(I,J)-GRGT(I,J))
        ELSE 
              IF(R==0) THEN
              GRR(I,J)=L_M(I,J)*GR(I,J)
              GAA(I,J)=L_M(I,J)*GA(I,J)
        END IF
           END IF
        END DO
END DO

 
END DO !Loop in R


!Loop for W-PHENERGY*R to get GLT
DO R=-M,M

DO I=1,C
	DO J=1,C
	GR(I,J)=0.0
	GA(I,J)=0.0
	GRLT(I,J)=0.0
	GRGT(I,J)=0.0
	END DO
END DO

!First loop to set self energies and coeff L_M

DO I=1,C
	DO J=1,C
		SE(I,J)=-0.5*(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
  		SE(I,J)=SE(I,J)*(GAMMAL(I,J)+GAMMAR(I,J))
  		!write(8,*),W,R,I,J,SE(I,J)

		SEGT(I,J)=-(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
		SEGT(I,J)=SEGT(I,J)*((1-FERMI(W+PHENERGY*R,MU_L,T_L))*GAMMAL(I,J)+(1-FERMI(W+PHENERGY*R,MU_R,T_R))*GAMMAR(I,J))   
		!write(8,*),W,R,I,J,SEGT(I,J)
		
                SELT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                SELT(I,J)=SELT(I,J)*(FERMI(W+PHENERGY*R,MU_L,T_L)*GAMMAL(I,J)+FERMI(W+PHENERGY*R,MU_R,T_R)*GAMMAR(I,J))
                !write(8,*),W,R,I,J,SELT(I,J)
                
                IF (I==J) THEN
                
                L_M(I,J)=2*(LAMBDA/PHENERGY)**2*SQRT(BOSE(PHENERGY,TEMPERATURE(I))*(BOSE(PHENERGY,TEMPERATURE(I))+1))
                !There is some problem here cause for R=1 we should obtain something finite but it seems to 
                !be related with higher decimals
                L_M(I,J)=PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))**2*EXP(R*PHENERGY/(TEMPERATURE(I)*2))*BESSEL(R,L_M(I,J)) 
                IF (R==-3) L_M(I,J)=0.0
                IF (R==-1) L_M(I,J)=0.0
                IF (R==-2) L_M(I,J)=0.0
                IF (R==0) L_M(I,J)=0.778801
                IF (R==1) L_M(I,J)=0.1947
                IF (R==2) L_M(I,J)=0.0243375
                IF (R==3) L_M(I,J)=0.00202813
                !write(8,*),W,R,I,J,L_M(I,J)
                
                    
                ELSE
             
                L_M(I,J)=PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                !write(8,*),W,R,I,J,SEGT(I,J),SELT(I,J),L_M(I,J)
                
             END IF
             
	END DO

END DO

                                !DRESSED RETARDED FUNCTION
				WX=W+ETA+PHENERGY*R
			
				
				!FILL THE COEFF. MATRIX 
 	                        DO Q=1,N
				            DO L=1,N
				            A(L,Q)=0.0
                                            END DO
					
				END DO
				
				DO IX=1,C 
	                 		DO Q=1,N
				                DO L=1,N
						
						IF (L==Q) THEN 
						A(L,Q)=1.0 !FILL THE DIAGONAL
						!write(8,*),I,J,L,Q,A(L,Q)
						
						ELSE IF (Q==L+IX*C) THEN !FILL UPPER DIAGONALS
						auxrow=INT((L - 1)/C) + 1
						auxcolumn=INT((Q - 1)/C) + 1
						A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))/(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
					        !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)
					        
						ELSE IF (Q==L-IX*C) THEN !FILL LOWER DIAGONALS
						auxcolumn=INT((Q - 1)/C) + 1
						auxrow=INT((L - 1)/C) + 1
						A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))/(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
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
                       
				
				DO L=1,N 
					X(L)=0.0
				END DO
				
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

! 				
! 				!FILL THE COEFF. MATRIX 
 	                        DO Q=1,N
				            DO L=1,N
				            A(L,Q)=0.0
                                            END DO
					
				END DO
				
				DO IX=1,C 
	                 		DO Q=1,N
				                DO L=1,N
						
						IF (L==Q) THEN 
						A(L,Q)=1.0 !FILL THE DIAGONAL
						!write(8,*),I,J,L,Q,A(L,Q)
						
						ELSE IF (Q==L+IX*C) THEN !FILL UPPER DIAGONALS
						auxrow=INT((L - 1)/C) + 1
						auxcolumn=INT((Q - 1)/C) + 1
						A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))/(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
					        !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)
					        
						ELSE IF (Q==L-IX*C) THEN !FILL LOWER DIAGONALS
						auxcolumn=INT((Q - 1)/C) + 1
						auxrow=INT((L - 1)/C) + 1
						A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))/(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
					        !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)
						
						END IF
						
						
						END DO
					END DO
					
				END DO

!                                 !FILL THE INDEPENDENT TERM 
                                DO L=1,N
                                        TERM(L)=0
                                        !write(9,*),I,J,L,TERM(L)
                                        IF (MOD(L,C+1)==1) THEN
                                                           P= INT(L/(C+1))+1
                                                           TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
                                                           !write(9,*),L,P,TERM(L)
                                        END IF
                 
                                END DO	

! 				
				DO L=1,N 
					X(L)=0.0
				END DO
! 				
! 				!CALL THE LINEAR SYSTEM SOLVER TO GET DRESSED GR
! 				
                                CALL SOLVER(A,TERM,X) 
!                                 
!                                 !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET DRESSED GR IN X
! 				!HERE, WE WANT JUST THE DIAGONAL TERMS (THOSE OF GR(1,1)...)
!                                
! 				!*****THIS SAVES THE DESIRED TERMS IN THE DIAGONAL*****!

				DO L=1,N
					auxrow=INT((L - 1)/C) + 1
					auxcolumn=MOD((L - 1),C) + 1
                                        GA(auxrow,auxcolumn)=X(L) 
                                        !write(8,*),L,auxrow,auxcolumn
				END DO

! 				!*****************************************************!
! 			

! !NOW OBTAIN THE ADVANCED GREEN FUNCTION BUT ALREADY MULTIPLIED BY THE COEFF L_M
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
!write(8,*),W+PHENERGY*R,GA(1,1),GA(2,2)
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
       END DO
END DO


DO I=1,C
	DO J=1,C
	IF(I==J) THEN
          GRR(I,J)=GRR(I,J)-L_M(I,J)*(0.5*GRLT(I,J))
          GAA(I,J)=GAA(I,J)+L_M(I,J)*(0.5*GRLT(I,J))
        ELSE 
              IF(R==0) THEN
              GRR(I,J)=L_M(I,J)*GR(I,J)
              GAA(I,J)=L_M(I,J)*GA(I,J)
        END IF
           END IF
        END DO
END DO

END DO !Loop in R


!SPECTRAL=(0,1)*(GRR-GAA)
SPECTRAL=(0,1)*(GGT-GLT)
TRANSCOEFF=MATMUL(GAA,MATMUL(GAMMAR,MATMUL(GRR,GAMMAL)))

!AUX=MATMUL(GAMMAL-GAMMAR,(0,1)*GLT)+MATMUL(GAMMAL*FERMI(W,MU_L,T_L)+GAMMAR*FERMI(W,MU_R,T_R),SPECTRAL)

 
WRITE(output_1,*) W,REAL(TRANSCOEFF(1,1)),-REAL((0.1)*TRANSCOEFF(1,1))
WRITE(output_2,*) W,REAL(TRANSCOEFF(1,1)+TRANSCOEFF(2,2)),-REAL((0.1)*(TRANSCOEFF(1,1)+TRANSCOEFF(2,2)))

WRITE(output_3,*) W,0.5*(REAL(SPECTRAL(1,1))+REAL(SPECTRAL(2,2)))
!WRITE(output_3,*) w,SPECTRAL(1,1),SPECTRAL(2,2)
WRITE(output_4,*) W,REAL(TRANSCOEFF(1,1)+TRANSCOEFF(2,2))

W=W+WINTERVAL

END DO

WRITE(output_5,*) V,JE
V=V+dV

END DO
CLOSE(output_1)
CLOSE(output_2)
CLOSE(output_3)
CLOSE(output_4)  
CLOSE(output_5)

PRINT*,"--- finished!"   

! PRINT*,"x=1.5,n=-1,jn=", BESSEL(-1,1.5)
! PRINT*,"x=1.5,n=1,jn=", BESSEL(1,1.5)
! PRINT*,"x=1.5,n=-2,jn=", BESSEL(-2,1.5)
! PRINT*,"x=1.5,n=2,jn=", BESSEL(2,1.5)
! PRINT*,"x=1.5,n=-3,jn=", BESSEL(-3,1.5)

STOP

                         

CONTAINS


REAL FUNCTION FERMI(FREC,MU,T)
	IMPLICIT NONE
	REAL,INTENT(IN)::FREC,MU,T
	REAL::OUTPUT
	IF(T==0) THEN
		IF(FREC>MU) OUTPUT=0
		IF(FREC<=MU) OUTPUT=1
	ELSE
		OUTPUT=1./(1+EXP((FREC-MU)/T))
	END IF
	FERMI=OUTPUT
END FUNCTION FERMI


REAL FUNCTION BOSE(FREC,T)
	IMPLICIT NONE
	REAL,INTENT(IN)::FREC,T
	REAL::OUTPUT
	IF(T==0) THEN
		OUTPUT=0.
	ELSE
		OUTPUT=1./(-1+EXP(FREC/T))
	END IF
     
	BOSE=OUTPUT
END FUNCTION BOSE

REAL FUNCTION XHI(FREC,T)
	IMPLICIT NONE
	REAL,INTENT(IN)::FREC,T
	REAL::OUTPUT
	OUTPUT=((LAMBDA/PHENERGY)**2)*(2*BOSE(FREC,T)+1)
	
	XHI=OUTPUT
END FUNCTION XHI

REAL FUNCTION PHONON_PART(LAMBDA,PHENERGY,T)
	IMPLICIT NONE
	REAL,INTENT(IN)::LAMBDA,PHENERGY,T
	REAL::OUTPUT
	OUTPUT=((LAMBDA/PHENERGY)**2)*(2*BOSE(PHENERGY,T)+1)/2
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

COMPLEX FUNCTION TRACE(N,A)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	COMPLEX,INTENT(IN)::A(N,N)
	COMPLEX::OUTPUT
	INTEGER::I,J
	DO I=1,N
		DO J=1,N
			IF (I==J) OUTPUT=OUTPUT+A(I,J)
		END DO
	END DO
	TRACE=OUTPUT
END FUNCTION TRACE

REAL FUNCTION BESSEL(N,X)
    IMPLICIT NONE
    INTEGER::N,N1,K,P,M
    REAL::X,OUTPUT,BK,BESSJ,BESSJ0,BESSJ1,BESSI,BESSI0,BESSI1
    N1=ABS(N)

    SELECT CASE(N1)
    CASE(0)
       BESSEL= BESSI0(N1,X)
    CASE(1)
       BESSEL= BESSI1(N1,X)
    CASE DEFAULT
       BESSEL= BESSI(N1,X)
    END SELECT

END FUNCTION BESSEL
  

! REAL FUNCTION BESSEL(N,X)
! 	IMPLICIT NONE
! 	INTEGER::N,N1
! 	REAL::X
! 	REAL::OUTPUT,BK
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

COMPLEX(8) FUNCTION Tr(A)

	IMPLICIT NONE

	COMPLEX(8),INTENT(IN)::A(2,2)

	COMPLEX(8)::output

	INTEGER::i,j

	output=A(1,1)+A(2,2)

	Tr=output

END FUNCTION

END PROGRAM
