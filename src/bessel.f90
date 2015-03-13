FUNCTION bessj(n,x)
  INTEGER n,IACC
  REAL bessj,x,BIGNO,BIGNI
  PARAMETER (IACC=40,BIGNO=1.e10,BIGNI=1.e-10)
  !USES bessj0,bessj1
  INTEGER j,jsum,m
  REAL ax,bj,bjm,bjp,sum,tox,bessj0,bessj1
  if(n.lt.2) PRINT*, 'bad argument n in bessj'
  ax=abs(x)
  if(ax.eq.0.)then
     bessj=0.
  else if(ax.gt.float(n))then
     tox=2./ax
     bjm=bessj0(ax)
     bj=bessj1(ax)
     do 11 j=1,n-1
        bjp=j*tox*bj-bjm
        bjm=bj
        bj=bjp
11   continue
     bessj=bj
  else
     tox=2./ax
     m=2*((n+int(sqrt(float(IACC*n))))/2)
     bessj=0.
     jsum=0
     sum=0.
     bjp=0.
     bj=1.
     do 12 j=m,1,-1
        bjm=j*tox*bj-bjp
        bjp=bj
        bj=bjm
        if(abs(bj).gt.BIGNO)then
           bj=bj*BIGNI
           bjp=bjp*BIGNI
           bessj=bessj*BIGNI
           sum=sum*BIGNI
        endif
        if(jsum.ne.0)sum=sum+bj
        jsum=1-jsum
        if(j.eq.n)bessj=bjp
12      continue
        sum=2.*sum-bj
        bessj=bessj/sum
     endif
     if(x.lt.0..and.mod(n,2).eq.1)bessj=-bessj
     return
  END
!    C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.

FUNCTION bessj1(x)
  REAL bessj1,x
  REAL ax,xx,z
  DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6, &
       s1,s2,s3,s4,s5,s6,y
  SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,&
       s5,s6
  DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0, &
       242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2, &
       s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0, &
       99447.43394d0,376.9991397d0,1.d0/
  DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4, &
       .2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0, &
       -.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
  if(abs(x).lt.8.)then
     y=x**2
     bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+ &
     y*(s4+y*(s5+y*s6)))))
  else
     ax=abs(x)
     z=8./ax
     y=z**2
     xx=ax-2.356194491
     bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y &
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.,x)
  endif
  return
END FUNCTION bessj1
!C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.

FUNCTION bessj0(x)
  REAL bessj0,x
  REAL ax,xx,z
  DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6, &
       s1,s2,s3,s4,s5,s6,y
  SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4, &
       s5,s6
  DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4, &
       -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1, &
       .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
  DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0, &
       651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2, &
       s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0, &
       59272.64853d0,267.8532712d0,1.d0/
  if(abs(x).lt.8.)then
     y=x**2
     bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y &
     *(s4+y*(s5+y*s6)))))
  else
     ax=abs(x)
     z=8./ax
     y=z**2
     xx=ax-.785398164
     bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y &
          *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
  endif
  return
END FUNCTION bessj0
!C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.

FUNCTION bessi(n,x)
      INTEGER n,IACC
      REAL bessi,x,BIGNO,BIGNI
      PARAMETER (IACC=40,BIGNO=1.0e10,BIGNI=1.0e-10)
!CU    USES bessi0
      INTEGER j,m
      REAL bi,bim,bip,tox,bessi0
      if(n.lt.2) PRINT*, 'bad argument n in bessj'
      if (x.eq.0.) then
        bessi=0.
      else
        tox=2.0/abs(x)
        bip=0.0
        bi=1.0
        bessi=0.
        m=2*((n+int(sqrt(float(IACC*n)))))
        do 11 j=m,1,-1
          bim=bip+float(j)*tox*bi
          bip=bi
          bi=bim
          if (abs(bi).gt.BIGNO) then
            bessi=bessi*BIGNI
            bi=bi*BIGNI
            bip=bip*BIGNI
          endif
          if (j.eq.n) bessi=bip
11      continue
        bessi=bessi*bessi0(x)/bi
        if (x.lt.0..and.mod(n,2).eq.1) bessi=-bessi
      endif
      return
END
!C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.

FUNCTION bessi0(x)
      REAL bessi0,x
      REAL ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,1.2067492d0, &
           0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1, &
           0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1, &
           0.2635537d-1,-0.1647633d-1,0.392377d-2/
      
       if (abs(x).lt.3.75) then
        y=(x/3.75)**2
        bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      else
        ax=abs(x)
        y=3.75/ax
        bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
      endif
      return
END
!C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.


      FUNCTION bessi1(x)
      REAL bessi1,x
      REAL ax
      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
      DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,&
     0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1, &
     -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,-0.2895312d-1, &
     0.1787654d-1,-0.420059d-2/
      if (abs(x).lt.3.75) then
        y=(x/3.75)**2
        bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        ax=abs(x)
        y=3.75/ax
        bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
        if(x.lt.0.)bessi1=-bessi1
      endif
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software 6?6>)AY.

