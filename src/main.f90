PROGRAM linearchain
IMPLICIT NONE
!!!**************VAR. DECLARATION**************!!!
!!!********************************************!!!
!!!*****************LEYEND*********************!!!
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
INTEGER,PARAMETER::C=3,N=C*C,M=2
INTEGER,PARAMETER::output_1=1, output_2=2, output_3=3
REAL,PARAMETER:: WMAX=3.0,WINTERVAL=0.1,DDH=0.2,MU_L=0.0,MU_R=0.0,T_L=1.0,T_R=1.0,LAMBDA=0.1,PHENERGY=1.00
COMPLEX*16,PARAMETER::ETA=0.1*(0,1)
COMPLEX*16::GR(C,C),GA(C,C),GLT(C,C),GGT(C,C),SPECTRAL(C,C),TERM(N),X(N),WX,AUX(C,C)
COMPLEX*16::A(N,N),SE(C,C),SELT(C,C),SEGT(C,C),T(C,C),GAMMAL(C,C),GAMMAR(C,C)
REAL::TEMPERATURE(C),ENERGY(N),W,L_M(C,C),AR(N,N),XR(N),TERMR(N)

INTEGER::I,J,K,L,Q,P,R,IX,IY,IZ,auxrow,auxcolumn

OPEN(unit=output_1,file="greater.txt",action="write",status="replace")
OPEN(unit=output_2,file="lesser.txt",action="write",status="replace")
OPEN(unit=output_3,file="spectral.txt",action="write",status="replace")

  !DOUBLE PRECISION jn
  !EXTERNAL jn

  PRINT*,"ThermWire"

!!**********END VAR. DECLARATION**************!!!

!FILL THE MATRICES AND TEMPERATURE AND ENERGY PROFILES!
!**********************HOPPING************************!
DO I=1,C
	DO J=1,C
		T(I,J)=0
		IF (I==J+1) T(I,J)=DDH
		IF (I==J-1) T(I,J)=DDH
	END DO

END DO
!*********************ENERGY PROFILE******************!
DO I=1,N
        ENERGY(I)=0.1 !RANDOM ENERGIES MUST BE INCLUDED
END DO
!*********************GAMMA MATRICES*************!
GAMMAL(1,1)=0.2
GAMMAL(C,C)=0.2
GAMMAR=GAMMAL

!******************TEMPERATURE PROFILE****************!
DO I=1,C
	TEMPERATURE(I)=T_L+I*(T_R-T_L)/C      !QUITE SIMPLE LINEAR PROFILE 
END DO
! !*****************************************************!
! 
W=-WMAX	     !HERE WE START THE MAIN LOOP
			!WE CALCULATE FIRST THE GGT FUNCTION AND THEN THE GLT
			!THE TWO PARTS ARE BASICALLY THE SAME
			!WITH MINOR DIFFERENCES
DO WHILE(W<=WMAX)
!*********************************GREATER THAN FUNCTION************************************************** 

PRINT*,"ThermWire(): W=", W

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
	END DO
END DO

!DRESSED RETARDED FUNCTION
DO R=-M,M

!First loop to set self energies and coeff L_M
DO I=1,C
	DO J=1,C
		
		IF (I==J) THEN
		SE(I,J)=0.5*(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
  		SE(I,J)=SE(I,J)*(GAMMAL(I,J)+GAMMAR(I,J))
  		!write(8,*),I,J,SE(I,J)
  		SEGT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
		SEGT(I,J)=SEGT(I,J)*((1-FERMI(W-PHENERGY*R,MU_L,T_L))*GAMMAL(I,J)+(1-FERMI(W-PHENERGY*R,MU_R,T_R))*GAMMAR(I,J))   
		SELT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                SELT(I,J)=SELT(I,J)*(FERMI(W-PHENERGY*R,MU_L,T_L)*GAMMAL(I,J)+FERMI(W-PHENERGY*R,MU_R,T_R)*GAMMAR(I,J))
                
                L_M(I,J)=(LAMBDA/PHENERGY)*SQRT((2*BOSE(PHENERGY,TEMPERATURE(I))+1)**2-1) 
		L_M(I,J)=PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*EXP(R*PHENERGY/(TEMPERATURE(I)*2))*BESSEL(R,L_M(I,J))
                !write(8,*),R,I,J,L_M(I,J)
                
                ELSE
		SEGT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
	        SEGT(I,J)=SEGT(I,J)*((1-FERMI(W,MU_L,T_L))*GAMMAL(I,J)+(1-FERMI(W,MU_R,T_R))*GAMMAR(I,J))
		SELT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                SELT(I,J)=SELT(I,J)*(FERMI(W,MU_L,T_L)*GAMMAL(I,J)+FERMI(W,MU_R,T_R)*GAMMAR(I,J))
                L_M(I,J)=PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                !write(8,*),R,I,J,L_M(I,J)
             END IF
	END DO

END DO


DO I=1,C
	DO J=1,C
		IF (I==J) THEN
				WX=W+ETA-PHENERGY*R
! 				
				
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
						auxrow=INT((Q - 1)/C) + 1
						auxcolumn=INT((L - 1)/C) + 1
						A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))/(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
					        !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)
						
						END IF
						
						
						END DO
					END DO
					
				END DO
				
! 				!Cheking the matrix 
! 	                        
	                         DO Q=1,N
				            DO L=1,N
				            !write(8,*),I,J,L,Q,A(L,Q)
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
					!write(8,*),TERM(L)
				END DO
				
				DO L=1,N 
					X(L)=0.0
				END DO
				
				!CALL THE LINEAR SYSTEM SOLVER TO GET DRESSED GR
				
                                CALL SOLVER(A,TERM,X) 
                                DO L=1,N 
					!write(8,*),I,J,X(L)
				END DO
                                
                                !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET DRESSED GR IN X
				!HERE, WE WANT JUST THE DIAGONAL TERMS (THOSE OF GR(1,1)...)
                               
				!*****THIS SAVES THE DESIRED TERMS IN THE DIAGONAL*****!
				
				DO L=1,N
                                        
                                        IF (MOD(L,C+1)==1) THEN
                                                           P= INT(L/(C+1))+1
                                                           GR(P,P)=X(L)
                                                           !write(9,*),L,P,GR(P,P)
                                        END IF
                 
                                END DO	
				
				
				!*****************************************************!
			
		ELSE
!			REPEAT FOR NON DIAGONAL TERMS SAME CODE BUT WX=W+ETA
				WX=W+ETA

				
!                               !FILL THE COEFF. MATRIX 
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
						auxrow=INT((Q - 1)/C) + 1
						auxcolumn=INT((L - 1)/C) + 1
						A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))/(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
					        !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)
						
						END IF
						
						
						END DO
					END DO
					
				END DO
				
! 				!Cheking the matrix 
	                        
! 	                        DO Q=1,N
! 				            DO L=1,N
! 				            write(8,*),I,J,L,Q,A(L,Q)
!                                             END DO
! 					
! 				END DO

                                !FILL THE INDEPENDENT TERM 
                                DO L=1,N
                                        TERM(L)=0
                                        !write(9,*),I,J,L,TERM(L)
                                        IF (MOD(L,C+1)==1) THEN
                                                           P= INT(L/(C+1))+1
                                                           TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
                                                           !write(9,*),I,J,L,TERM(L)
                                        END IF
                 
                                END DO	
                                
!                              DO L=1,N 
! 					write(9,*),TERM(L)
! 				END DO
				
				DO L=1,N 
					X(L)=0.0
				END DO
! 
! 				!RECALL THE LINEAR SOLVER AND GET AGAIN DRESSED GR. NOW WE SAVE THE TERMS IN AN AUXILIAR MATRIX
! 				!BEFORE THE SUSTITUTION
				CALL SOLVER(A,TERM,X)
                                DO L=1,N 
					!write(8,*),I,J,X(L)
				END DO
				
                             	DO L=1,N
					auxrow=INT((L - 1)/C) + 1
					auxcolumn=MOD((L - 1),C) + 1
                                        AUX(auxrow,auxcolumn)=X(L) 
                                        !write(9,*),L,auxrow,auxcolumn
				END DO
				
				

                               
 		END IF
 	END DO
END DO
!NOW OBTAIN THE RETARDED DRESSED GREEN FUNCTION
DO I=1,C
	DO J=1,C
           IF(I==J) THEN
              GR(I,J)=GR(I,J)
              !write(9,*),I,J,GR(I,J)
           ELSE
              GR(I,J)=AUX(I,J)
              !write(9,*),I,J,GR(I,J)
           END IF
       END DO
END DO

!DRESSED ADVANCED FUNCTION
DO I=1,C
	DO J=1,C
		
		IF (I==J) THEN
  		SE(I,J)=-SE(I,J)
  		!write(8,*),I,J,SE(I,J)
  		!SEGT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
		!SEGT(I,J)=SEGT(I,J)*((1-FERMI(W-PHENERGY*R,MU_L,T_L))*GAMMAL(I,J)+(1-FERMI(W-PHENERGY*R,MU_R,T_R))*GAMMAR(I,J))   
		ELSE
		!SEGT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
	        !SEGT(I,J)=SEGT(I,J)*((1-FERMI(W,MU_L,T_L))*GAMMAL(I,J)+(1-FERMI(W,MU_R,T_R))*GAMMAR(I,J))
		END IF
	END DO

END DO
! 
DO I=1,C
	DO J=1,C
		IF (I==J) THEN
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
						auxrow=INT((Q - 1)/C) + 1
						auxcolumn=INT((L - 1)/C) + 1
						A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))/(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
					        !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)
						
						END IF
						
						
						END DO
					END DO
					
				END DO
! 				
! ! 				!Cheking the matrix 
! 	                        DO Q=1,N
! 				            DO L=1,N
! 				            write(8,*),I,J,L,Q,A(L,Q)
!                                             END DO
! 					
! 				END DO
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
! !                              DO L=1,N 
! ! 					write(9,*),TERM(L)
! ! 				END DO
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
! 				
! 				
! 				!*****************************************************!
! 			
		ELSE
! !			REPEAT FOR NON DIAGONAL TERMS SAME CODE BUT WX=W+ETA
				WX=W-ETA

! ! 				
! !                            !FILL THE COEFF. MATRIX 
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
						auxrow=INT((Q - 1)/C) + 1
						auxcolumn=INT((L - 1)/C) + 1
						A(L,Q)=-(T(auxrow,auxcolumn)+SE(auxrow,auxcolumn))/(WX-ENERGY(auxrow)-SE(auxrow,auxrow))
					        !write(8,*),I,J,L,Q,auxrow,auxcolumn,A(L,Q)
						
						END IF
						
						
						END DO
					END DO
					
				END DO
! 				
! ! 				!Cheking the matrix 
! 	                        DO Q=1,N
! 				            DO L=1,N
! 				            write(8,*),I,J,L,Q,A(L,Q)
!                                             END DO
! 					
! 				END DO
! 
                                !FILL THE INDEPENDENT TERM 
                                DO L=1,N
                                        TERM(L)=0
                                        !write(9,*),I,J,L,TERM(L)
                                        IF (MOD(L,C+1)==1) THEN
                                                           P= INT(L/(C+1))+1
                                                           TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
                                                           !write(9,*),I,J,L,TERM(L)
                                        END IF
                 
                                END DO	
                                
! !                             DO L=1,N 
! ! 					write(9,*),TERM(L)
! ! 				END DO
! 				
				DO L=1,N 
					X(L)=0.0
				END DO
! ! 
! ! 				!RECALL THE LINEAR SOLVER AND GET AGAIN DRESSED GR. NOW WE SAVE THE TERMS IN AN AUXILIAR MATRIX
! ! 				!BEFORE THE SUSTITUTION
				CALL SOLVER(A,TERM,X)

                             	DO L=1,N
					auxrow=INT((L - 1)/C) + 1
					auxcolumn=MOD((L - 1),C) + 1
                                        AUX(auxrow,auxcolumn)=X(L) 
                                        !write(9,*),L,auxrow,auxcolumn
				END DO
				
				

                               
 		END IF
 	END DO
END DO
! !NOW OBTAIN THE ADVANCED GREEN FUNCTION BUT ALREADY MULTIPLIED BY THE COEFF L_M
DO I=1,C
	DO J=1,C
           IF(I==J) THEN
              GA(I,J)=L_M(I,J)*GA(I,J)
              !write(8,*),I,J,GA(I,J)
           ELSE
              GA(I,J)=L_M(I,J)*AUX(I,J)
              !write(8,*),I,J,GA(I,J)
           END IF
       END DO
END DO


GGT=GGT+MATMUL(GR,MATMUL(SEGT,GA))
GLT=GLT+MATMUL(GR,MATMUL(SELT,GA))
 
END DO !Loop in R

DO I=1,C
	DO J=1,C
           write(9,*),I,J,GGT(I,J),GLT(I,J)
          
       END DO
END DO
WRITE(output_1,*) w,GGT(1,1)
WRITE(output_2,*) w,GLT(1,1)

SPECTRAL=(0,1)*(GGT-GLT)
WRITE(output_3,*) w,REAL(SPECTRAL(1,1))     



W=W+WINTERVAL

END DO

CLOSE(output_1)
CLOSE(output_2)
CLOSE(output_3)
   
PRINT*,"--- finished!"
PRINT*,"x=1.0,n=1,jn=", BESSEL(-1,1.5)
STOP

                         

CONTAINS


REAL FUNCTION FERMI(W,MU,T)
	IMPLICIT NONE
	REAL,INTENT(IN)::W,MU,T
	REAL::OUTPUT
	IF(T==0) THEN
		IF(W>MU) OUTPUT=0
		IF(W<=MU) OUTPUT=1
	ELSE
		OUTPUT=1./(1+EXP((W-MU)/T))
	END IF
	FERMI=OUTPUT
END FUNCTION FERMI


REAL FUNCTION BOSE(W,T)
	IMPLICIT NONE
	REAL,INTENT(IN)::W,T
	REAL::OUTPUT
	IF(T==0) THEN
		OUTPUT=0.
	ELSE
		OUTPUT=1./(-1+EXP(W/T))
	END IF
	BOSE=OUTPUT
END FUNCTION BOSE

REAL FUNCTION XHI(W,T)
	IMPLICIT NONE
	REAL,INTENT(IN)::W,T
	REAL::OUTPUT
	OUTPUT=((LAMBDA/PHENERGY)**2)*(2*BOSE(W,T)+1)
	
	XHI=OUTPUT
END FUNCTION XHI

REAL FUNCTION PHONON_PART(LAMBDA,PHENERGY,T)
	IMPLICIT NONE
	REAL,INTENT(IN)::LAMBDA,PHENERGY,T
	REAL::OUTPUT
	OUTPUT=((LAMBDA/PHENERGY)**2)*(2*BOSE(W,T)+1)/2
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
    REAL::X,OUTPUT,BK,BESSJ,BESSJ0,BESSJ1
    N1=ABS(N)
    M=6
  
    IF (N1<2) THEN 
        OUTPUT=(X/2)**(N1)
        !PRINT*,"BESSEL()",N,N1,OUTPUT
	DO K=1,M
			BK=1./((2**(2*K+N1))*FACTORIAL(K)*FACTORIAL(K+N1))
  			OUTPUT=OUTPUT+BK*((X)**(2*K+N1)) 
                        !PRINT*,"BESSEL()",N,N1,K,BK,OUTPUT
	END DO
	BESSEL=OUTPUT
        !PRINT*,"BESSEL()",BESSEL

    ELSE IF (N1>=2) THEN 
		BESSEL= BESSJ(N1,X)
    END IF

END FUNCTION BESSEL



SUBROUTINE SOLVER(MAT, RHS, VEC)

    COMPLEX*16::MAT(N,N), RHS(N), VEC(N)
    
    INTEGER::INFO
    INTEGER::IPIV(N)
    
    PRINT*,"SOLVER()"
    
    !PRINT*,"before: mat=", MAT, "rhs=", RHS

    CALL ZGESV( N, 1, MAT, N, IPIV, RHS, N, INFO )
    
    IF(INFO.NE.0) THEN
       PRINT*,"SOLVER(): INFO=", INFO," --- WRNG!"
    ENDIF

    !PRINT*,"after: mat=", MAT, "rhs=", RHS

    VEC= RHS
    
END SUBROUTINE SOLVER

END PROGRAM
