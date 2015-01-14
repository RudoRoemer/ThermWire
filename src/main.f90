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

  INTEGER,PARAMETER::C=2,N=C*C,M=3
  REAL,PARAMETER:: WMAX=3,WINTERVAL=0.1,DDH=0.2,MU_L=0,MU_R=0,T_L=0,T_R=5,LAMBDA=0.1,PHENERGY=1
  COMPLEX,PARAMETER::ETA=0.001
  COMPLEX::GR(C,C),GA(C,C),GLT(C,C),GGT(C,C),SPECTRAL(C,C),TERM(N),X(N),WX,AUX(C,C)
  COMPLEX::A(N,N),SE(C,C),SELT(C,C),SEGT(C,C),T(C,C),GAMMAL(C,C),GAMMAR(C,C)
  REAL::TEMPERATURE(C),ENERGY(C),W
  INTEGER::I,J,K,L,Q,P,R

  PRINT*,"ThermWire"

  !**********END VAR. DECLARATION**************!!!

  !FILL THE MATRICES AND TEMPERATURE AND ENERGY PROFILES!
  !**********************HOPPING************************!
  DO I=1,C
     DO J=1,C
        T(I,J)=DDH !SAME HOPPING FOR ALL THE DOTS
        IF (I==J) T(I,J)=0
     END DO
  END DO
  !*********************ENERGY PROFILE******************!
  ENERGY=0.1 !RANDOM ENERGIES MUST BE INCLUDED
  !*********************GAMMA MATRICES*************!
  GAMMAL(1,1)=0.2
  GAMMAL(C,C)=0.2
  GAMMAR=GAMMAL

  !******************TEMPERATURE PROFILE****************!
  DO I=1,C
     TEMPERATURE(I)=I*(T_R-T_L)/C      !QUITE SIMPLE LINEAR PROFILE
  END DO
  !*****************************************************!

  W=-WMAX	     !HERE WE START THE MAIN LOOP
  !WE CALCULATE FIRST THE GGT FUNCTION AND THEN THE GLT
  !THE TWO PARTS ARE BASICALLY THE SAME
  !WITH MINOR DIFFERENCES
  DO WHILE(W<=WMAX)
     !*********************************GREATER THAN FUNCTION************************************************** 

     PRINT*,"ThermWire(): W=", W

     !DRESSED RETARDED FUNCTION
     DO I=1,C
	DO J=1,C
           !CALCULATE THE RETARDED SELF ENERGY AND GREATER THAN SELF ENERGY
           SE(I,J)=0.5*(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
           SE(I,J)=SE(I,J)*(GAMMAL(I,J)+GAMMAR(I,J))
           IF (I==J) THEN
              DO R=-M,M
                 SEGT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                 SEGT(I,J)=SEGT(I,J)*((1-FERMI(W-PHENERGY*M,MU_L,T_L))*GAMMAL(I,J)+(1-FERMI(W-PHENERGY*M,MU_R,T_R))*GAMMAR(I,J))
                 WX=W+ETA-PHENERGY*M
                 !FILL THE COEFF. MATRIX (BUT NOT THE DIAGONAL)
                 DO L=K,N
                    DO Q=K,N,C
                       IF (I==J) THEN 
                       ELSE
                          A(L,Q)=-(T(L,Q)+SE(L,Q))
                          A(L,:)=A(L,:)/(WX-ENERGY(L)-SE(L,L))
                       END IF
                    END DO
                    K=K+1
                 END DO
                 DO L=1,N
                    A(:,L)=A(L,:)
                    DO Q=1,N
                       IF (L==Q) A(L,Q)=1 !FILL THE DIAGONAL
                    END DO
                 END DO
                 !FILL THE INDEPENDENT TERM
                 DO L=1,N,C+1
                    !THIS IS QUITE TRICKY BECAUSE OF THE ORDERING OF THE TERMS IN THE SYSTEM
                    !TERM(1) IS RELATED WITH ENERGY(1), SE(1,1)
                    !TERM(5) IS RELATED WITH ENERGY(2), SE(2,2)
                    !TERM(9) IS RELATED WITH ENERGY(3),SE(3,3)
                    !AND SO ON
                    IF(L==1) P=1
                    IF(L==N) P=C
                    IF((L.NE.1).AND.(L.NE.N)) THEN
                       P=INT(SQRT(1.*(L*L-(C+2))))
                    END IF
                    TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
                 END DO
                 !CALL THE LINEAR SYSTEM SOLVER TO GET GR
                 CALL SOLVER(A,TERM,X) !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET GR IN X
                 !HERE, WE WANT JUST THE DIAGONAL TERMS (THOSE OF GR(1,1)...)
                 !*****THIS SAVES THE DESIRED TERMS IN THE DIAGONAL*****!
                 Q=1
                 DO L=1,C
                    GR(L,L)=X(Q)
                    Q=Q+1+C
                 END DO
                 !*****************************************************!
              END DO
           ELSE
              !REPEAT FOR NON DIAGONAL TERMS SAME CODE BUT WX=W+ETA
              WX=W+ETA
              SE(I,J)=0.5*(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
              SE(I,J)=SE(I,J)*(GAMMAL(I,J)+GAMMAR(I,J))
              SEGT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
              SEGT(I,J)=SEGT(I,J)*((1-FERMI(W,MU_L,T_L))*GAMMAL(I,J)+(1-FERMI(W,MU_R,T_R))*GAMMAR(I,J))
              !FILL THE COEFF. MATRIX
              DO L=K,N
                 DO Q=K,N,C
                    IF (I==J) THEN
                    ELSE
                       A(L,Q)=-(T(L,Q)+SE(L,Q))
                       A(L,:)=A(L,:)/(WX-ENERGY(L)-SE(L,L))
                    END IF
                 END DO
                 K=K+1
              END DO
              DO L=1,N
                 A(:,L)=A(L,:)
                 DO Q=1,N
                    IF (L==Q) A(L,Q)=1 !FILL THE DIAGONAL
                 END DO
              END DO
              !FILL THE INDEPENDENT TERM
              DO L=1,N,C+1
                 !THIS IS QUITE TRICKY BECAUSE OF THE ORDERING OF THE TERMS IN THE SYSTEM
                 !TERM(1) IS RELATED WITH ENERGY(1), SE(1,1)
                 !TERM(5) IS RELATED WITH ENERGY(2), SE(2,2)
                 !TERM(9) IS RELATED WITH ENERGY(3),SE(3,3)
                 !AND SO ON
                 IF(L==1) P=1
                 IF(L==N) P=C
                 IF((L.NE.1).AND.(L.NE.N)) THEN
                    P=INT(SQRT(1.*(L*L-(C+2))))
                 END IF
                 TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
              END DO
              !RECALL THE LINEAR SOLVER AND GET AGAIN GR. NOW WE SAVE THE TERMS IN AN AUXILIAR MATRIX
              !BEFORE THE SUSTITUTION
              CALL SOLVER(A,TERM,X) !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET GR IN X
              Q=1
              DO L=1,C
                 AUX(L,:)=X(Q:Q+(C-1))
                 Q=Q+C
              END DO
           END IF
	END DO
     END DO
     !NOW OBTAIN THE RETARDED GREEN FUNCTION
     DO I=1,C
	DO J=1,C
           IF(I==J) THEN
              GR(I,J)=GR(I,J)
           ELSE
              GR(I,J)=AUX(I,J)
           END IF
	END DO
     END DO
     !DRESSED ADVANCED FUNCTION
     SE=-SE  !NOW THE SELF-ENERGY IS CONJUGATED
     DO I=1,C
	DO J=1,C
        !CALCULATE THE LESSER THAN SELF ENERGY
           IF (I==J) THEN
              DO R=-M,M
                 WX=W-ETA-PHENERGY*M
                 !FILL THE COEFF. MATRIX (BUT NOT THE DIAGONAL)
                 DO L=K,N
                    DO Q=K,N,C
                       IF (I==J) THEN 
                       ELSE
                          A(L,Q)=-(T(L,Q)+SE(L,Q))
                          A(L,:)=A(L,:)/(WX-ENERGY(L)-SE(L,L))
                       END IF
                    END DO
                    K=K+1
                 END DO
                 DO L=1,N
                    A(:,L)=A(L,:)
                    DO Q=1,N
                       IF (L==Q) A(L,Q)=1 !FILL THE DIAGONAL
                    END DO
                 END DO
                 !FILL THE INDEPENDENT TERM
                 DO L=1,N,C+1
                    !THIS IS QUITE TRICKY BECAUSE OF THE ORDERING OF THE TERMS IN THE SYSTEM
                    !TERM(1) IS RELATED WITH ENERGY(1), SE(1,1)
                    !TERM(5) IS RELATED WITH ENERGY(2), SE(2,2)
                    !TERM(9) IS RELATED WITH ENERGY(3),SE(3,3)
                    !AND SO ON
                    IF(L==1) P=1
                    IF(L==N) P=C
                    IF((L.NE.1).AND.(L.NE.N)) THEN
                       P=INT(SQRT(1.*(L*L-(C+2))))
                    END IF
                    TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
                 END DO
                 !CALL THE LINEAR SYSTEM SOLVER TO GET GA
                 CALL SOLVER(A,TERM,X) !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET GA IN X
                 !HERE WE WANT JUST THE DIAGONAL TERMS GA(1,1)...

                 !*****THIS SAVES THE DESIRED TERMS IN THE DIAGONAL*****!
                 Q=1
                 DO L=1,C
                    GA(L,L)=X(Q)
                    Q=Q+1+C
                 END DO
                 !*****************************************************!
              END DO
           ELSE
              !REPEAT FOR NON DIAGONAL TERMS SAME CODE BUT WX=W-ETA
              WX=W-ETA

              !FILL THE COEFF. MATRIX
              DO L=K,N
                 DO Q=K,N,C
                    IF (I==J) THEN
                    ELSE
                       A(L,Q)=-(T(L,Q)+SE(L,Q))
                       A(L,:)=A(L,:)/(WX-ENERGY(L)-SE(L,L))
                    END IF
                 END DO
                 K=K+1
              END DO
              DO L=1,N
                 A(:,L)=A(L,:)
                 DO Q=1,N
                    IF (L==Q) A(L,Q)=1 !FILL THE DIAGONAL
                 END DO
              END DO
              !FILL THE INDEPENDENT TERM
              DO L=1,N,C+1
                 !THIS IS QUITE TRICKY BECAUSE OF THE ORDERING OF THE TERMS IN THE SYSTEM
                 !TERM(1) IS RELATED WITH ENERGY(1), SE(1,1)
                 !TERM(5) IS RELATED WITH ENERGY(2), SE(2,2)
                 !TERM(9) IS RELATED WITH ENERGY(3),SE(3,3)
                 !AND SO ON
                 IF(L==1) P=1
                 IF(L==N) P=C
                 IF((L.NE.1).AND.(L.NE.N)) THEN
                    P=INT(SQRT(1.*(L*L-(C+2))))
                 END IF
                 TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
              END DO
              !RECALL THE LINEAR SOLVER AND GET AGAIN GA. NOW WE SAVE THE TERMS IN AN AUXILIAR MATRIX
              !BEFORE THE SUSTITUTION
              CALL SOLVER(A,TERM,X) !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET GA IN X

              Q=1
              DO L=1,C
                 AUX(L,:)=X(Q:Q+(C-1))
                 Q=Q+C
              END DO
           END IF
	END DO
     END DO
     !NOW OBTAIN THE ADVANCED GREEN FUNCTION


     GGT=MATMUL(GR,MATMUL(SEGT,GA))

     !*********************************LESSER THAN FUNCTION************************************************** 

     !DRESSED RETARDED FUNCTION
     DO I=1,C
	DO J=1,C
           !CALCULATE THE RETARDED SELF ENERGY AND LESSER THAN SELF ENERGY
           SE(I,J)=0.5*(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
           SE(I,J)=SE(I,J)*(GAMMAL(I,J)+GAMMAR(I,J))
           IF (I==J) THEN
              DO R=-M,M
                 SELT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
                 SELT(I,J)=SEGT(I,J)*((1-FERMI(W-PHENERGY*M,MU_L,T_L))*GAMMAL(I,J)+(1-FERMI(W-PHENERGY*M,MU_R,T_R))*GAMMAR(I,J))
                 WX=W+ETA+PHENERGY*M
                 !FILL THE COEFF. MATRIX (BUT NOT THE DIAGONAL)
                 DO L=K,N
                    DO Q=K,N,C
                       IF (I==J) THEN 
                       ELSE
                          A(L,Q)=-(T(L,Q)+SE(L,Q))
                          A(L,:)=A(L,:)/(WX-ENERGY(L)-SE(L,L))
                       END IF
                    END DO
                    K=K+1
                 END DO
                 DO L=1,N
                    A(:,L)=A(L,:)
                    DO Q=1,N
                       IF (L==Q) A(L,Q)=1 !FILL THE DIAGONAL
                    END DO
                 END DO
                 !FILL THE INDEPENDENT TERM
                 DO L=1,N,C+1
                    !THIS IS QUITE TRICKY BECAUSE OF THE ORDERING OF THE TERMS IN THE SYSTEM
                    !TERM(1) IS RELATED WITH ENERGY(1), SE(1,1)
                    !TERM(5) IS RELATED WITH ENERGY(2), SE(2,2)
                    !TERM(9) IS RELATED WITH ENERGY(3),SE(3,3)
                    !AND SO ON
                    IF(L==1) P=1
                    IF(L==N) P=C
                    IF((L.NE.1).AND.(L.NE.N)) THEN
                       P=INT(SQRT(1.*(L*L-(C+2))))
                    END IF
                    TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
                 END DO
                 !CALL THE LINEAR SYSTEM SOLVER TO GET GR
                 CALL SOLVER(A,TERM,X) !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET GR IN X
                 !HERE, WE WANT JUST THE DIAGONAL TERMS (THOSE OF GR(1,1)...)
                 !*****THIS SAVES THE DESIRED TERMS IN THE DIAGONAL*****!
                 Q=1
                 DO L=1,C
                    GR(L,L)=X(Q)
                    Q=Q+1+C
                 END DO
                 !*****************************************************!
              END DO
           ELSE
              !REPEAT FOR NON DIAGONAL TERMS SAME CODE BUT WX=W+ETA
              WX=W+ETA
              SE(I,J)=0.5*(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
              SE(I,J)=SE(I,J)*(GAMMAL(I,J)+GAMMAR(I,J))
              SELT(I,J)=(0,1)*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(I))*PHONON_PART(LAMBDA,PHENERGY,TEMPERATURE(J))
              SELT(I,J)=SELT(I,J)*((1-FERMI(W,MU_L,T_L))*GAMMAL(I,J)+(1-FERMI(W,MU_R,T_R))*GAMMAR(I,J))
              !FILL THE COEFF. MATRIX
              DO L=K,N
                 DO Q=K,N,C
                    IF (I==J) THEN
                    ELSE
                       A(L,Q)=-(T(L,Q)+SE(L,Q))
                       A(L,:)=A(L,:)/(WX-ENERGY(L)-SE(L,L))
                    END IF
                 END DO
                 K=K+1
              END DO
              DO L=1,N
                 A(:,L)=A(L,:)
                 DO Q=1,N
                    IF (L==Q) A(L,Q)=1 !FILL THE DIAGONAL
                 END DO
              END DO
              !FILL THE INDEPENDENT TERM
              DO L=1,N,C+1
                 !THIS IS QUITE TRICKY BECAUSE OF THE ORDERING OF THE TERMS IN THE SYSTEM
                 !TERM(1) IS RELATED WITH ENERGY(1), SE(1,1)
                 !TERM(5) IS RELATED WITH ENERGY(2), SE(2,2)
                 !TERM(9) IS RELATED WITH ENERGY(3),SE(3,3)
                 !AND SO ON
                 IF(L==1) P=1
                 IF(L==N) P=C
                 IF((L.NE.1).AND.(L.NE.N)) THEN
                    P=INT(SQRT(1.*(L*L-(C+2))))
                 END IF
                 TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
              END DO
              !RECALL THE LINEAR SOLVER AND GET AGAIN GR. NOW WE SAVE THE TERMS IN AN AUXILIAR MATRIX
              !BEFORE THE SUSTITUTION
              CALL SOLVER(A,TERM,X) !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET GR IN X
              Q=1
              DO L=1,C
                 AUX(L,:)=X(Q:Q+(C-1))
                 Q=Q+C
              END DO
           END IF
	END DO
     END DO
     !NOW OBTAIN THE RETARDED GREEN FUNCTION
     DO I=1,C
	DO J=1,C
           IF(I==J) THEN
              GR(I,J)=GR(I,J)
           ELSE
              GR(I,J)=AUX(I,J)
           END IF
	END DO
     END DO
     !DRESSED ADVANCED FUNCTION
     SE=-SE  !NOW THE SELF-ENERGY IS CONJUGATED
     DO I=1,C
	DO J=1,C
    !CALCULATE THE LESSER THAN SELF ENERGY
           IF (I==J) THEN
              DO R=-M,M
                 WX=W-ETA+PHENERGY*M
                 !FILL THE COEFF. MATRIX (BUT NOT THE DIAGONAL)
                 DO L=K,N
                    DO Q=K,N,C
                       IF (I==J) THEN 
                       ELSE
                          A(L,Q)=-(T(L,Q)+SE(L,Q))
                          A(L,:)=A(L,:)/(WX-ENERGY(L)-SE(L,L))
                       END IF
                    END DO
                    K=K+1
                 END DO
                 DO L=1,N
                    A(:,L)=A(L,:)
                    DO Q=1,N
                       IF (L==Q) A(L,Q)=1 !FILL THE DIAGONAL
                    END DO
                 END DO
                 !FILL THE INDEPENDENT TERM
                 DO L=1,N,C+1
                    !THIS IS QUITE TRICKY BECAUSE OF THE ORDERING OF THE TERMS IN THE SYSTEM
                    !TERM(1) IS RELATED WITH ENERGY(1), SE(1,1)
                    !TERM(5) IS RELATED WITH ENERGY(2), SE(2,2)
                    !TERM(9) IS RELATED WITH ENERGY(3),SE(3,3)
                    !AND SO ON
                    IF(L==1) P=1
                    IF(L==N) P=C
                    IF((L.NE.1).AND.(L.NE.N)) THEN
                       P=INT(SQRT(1.*(L*L-(C+2))))
                    END IF
                    TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
                 END DO
                 !CALL THE LINEAR SYSTEM SOLVER TO GET GA
                 CALL SOLVER(A,TERM,X) !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET GA IN X
                 !HERE WE WANT JUST THE DIAGONAL TERMS GA(1,1)...

                 !*****THIS SAVES THE DESIRED TERMS IN THE DIAGONAL*****!
                 Q=1
                 DO L=1,C
                    GA(L,L)=X(Q)
                    Q=Q+1+C
                 END DO
                 !*****************************************************!
              END DO
           ELSE
              !REPEAT FOR NON DIAGONAL TERMS SAME CODE BUT WX=W-ETA
              WX=W-ETA

              !FILL THE COEFF. MATRIX
              DO L=K,N
                 DO Q=K,N,C
                    IF (I==J) THEN
                    ELSE
                       A(L,Q)=-(T(L,Q)+SE(L,Q))
                       A(L,:)=A(L,:)/(WX-ENERGY(L)-SE(L,L))
                    END IF
                 END DO
                 K=K+1
              END DO
              DO L=1,N
                 A(:,L)=A(L,:)
                 DO Q=1,N
                    IF (L==Q) A(L,Q)=1 !FILL THE DIAGONAL
                 END DO
              END DO
              !FILL THE INDEPENDENT TERM
              DO L=1,N,C+1
                 !THIS IS QUITE TRICKY BECAUSE OF THE ORDERING OF THE TERMS IN THE SYSTEM
                 !TERM(1) IS RELATED WITH ENERGY(1), SE(1,1)
                 !TERM(5) IS RELATED WITH ENERGY(2), SE(2,2)
                 !TERM(9) IS RELATED WITH ENERGY(3),SE(3,3)
                 !AND SO ON
                 IF(L==1) P=1
                 IF(L==N) P=C
                 IF((L.NE.1).AND.(L.NE.N)) THEN
                    P=INT(SQRT(1.*(L*L-(C+2))))
                 END IF
                 TERM(L)=1./(WX-ENERGY(P)-SE(P,P)) 
              END DO
              !RECALL THE LINEAR SOLVER AND GET AGAIN GA. NOW WE SAVE THE TERMS IN AN AUXILIAR MATRIX
              !BEFORE THE SUSTITUTION
              CALL SOLVER(A,TERM,X) !THIS SHOULD SOLVE THE LINEAR SYSTEM AND WE GET GA IN X

              Q=1
              DO L=1,C
                 AUX(L,:)=X(Q:Q+(C-1))
                 Q=Q+C
              END DO
           END IF
	END DO
     END DO
     !NOW OBTAIN THE ADVANCED GREEN FUNCTION
     DO I=1,C
	DO J=1,C
           IF(I==J) THEN
              GA(I,J)=GA(I,J)
           ELSE
              GA(I,J)=AUX(I,J)
           END IF
	END DO
     END DO
     !NOW WE OBTAIN GLT

     GLT=MATMUL(GR,MATMUL(SELT,GA))

     W=W+WINTERVAL

  END DO

  PRINT*,"--- finished!"

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
    INTEGER::N,N1
    REAL::X
    REAL::OUTPUT,BK
    INTEGER::K,P,M
    N1=ABS(N)
    M=12
    OUTPUT=0
    P=1
    IF(N<0) P=(-1)**N
    DO K=0,M
       BK=1./((2**(2*K+N1))*FACTORIAL(K)*FACTORIAL(K+N1))
       OUTPUT=OUTPUT+BK*((X)**(2*K+N1))
    END DO
    OUTPUT=OUTPUT*P
    BESSEL=OUTPUT
  END FUNCTION BESSEL

  SUBROUTINE SOLVER(MAT, RHS, VEC)

    COMPLEX::MAT(N,N), RHS(N), VEC(N)

    INTEGER::INFO
    INTEGER::IPIV(N)

    PRINT*,"SOLVER()"

    CALL ZGESV( N, 1, MAT, N, IPIV, RHS, N, INFO )

    IF(INFO.NE.0) THEN
       PRINT*,"SOLVER(): INFO=", INFO," --- WRNG!"
    ENDIF

    VEC= RHS

  END SUBROUTINE SOLVER

END PROGRAM linearchain
