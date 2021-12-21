      Subroutine BetaCI(E0,OneH,ERI,T1,T2,NA,X,Y,R0,R1,R2)
          Implicit None

          Integer, parameter  :: pr = Selected_Real_Kind(15,307)
          Integer, Intent(In) :: NA
          Real (Kind=pr), Intent(In) :: X(NA), Y(NA)
          Real (Kind=pr), Intent(In) :: T1(NA,NA), T2(NA,NA,NA,NA)
          Real (Kind=pr), Intent(In) :: E0(NA), OneH(NA,NA),  ERI(NA,NA,NA,NA)
          Real (Kind=pr), Intent(Out) :: R0
          Real (Kind=pr), Intent(Out) :: R1(NA,NA)
          Real (Kind=pr), Intent(Out) :: R2(NA,NA,NA,NA)
          Real (Kind=pr)              :: Amat(NA,NA,NA,NA)
          Integer :: a, b, c, d, i

          Real (Kind=pr) :: tau0
          Real (Kind=pr) :: tau1
          Real (Kind=pr), dimension(:), allocatable :: tau2
          Real (Kind=pr) :: tau3
          Real (Kind=pr) :: tau4
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau5
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau7
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau8
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau9
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau10
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau11

          ! Pre Processing
          ! Defining the Hamiltonian Matrix Elements

          Real (Kind=pr) ::  h0, h20(na,na), h02(na,na), h11(na,na)
          Real (Kind=pr) ::  h40(na,na,na,na), h04(na,na,na,na)
          Real (Kind=pr) ::  h31(na,na,na,na), h13(na,na,na,na)
          Real (Kind=pr) ::  h221(na,na,na,na)
          Real (Kind=pr) ::  h222(na,na,na,na)
          Real (Kind=pr) ::  scr1(na), scr2(na,na), delK(na,na)

          h0 = 0.0_pr
          delK = 0.0_pr

          do a=1, na
              h0 = h0 + y(a)*y(a)*oneh(a,a)
              delK(a,a) = 1.0_pr
              do b=1, na
                  scr2(a,b) = 0.0_pr
                  do c=1,na
                      scr1(c) = eri(a,c,b,c)
                  end do
                  scr2(a,b) = Sum(y*y*scr1) + oneh(a,b)
                  h0 = h0 + ( y(a)**2 * y(b)**2 * eri(a,b,a,b) )/2
              end do
          end do

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  h20(a,b) = x(a)*x(b)*scr2(a,b)
                  h02(a,b) = -y(a)*y(b)*scr2(a,b)
                  h11(a,b) = x(a)*y(b)*scr2(a,b)
                  do c=1, na
                      do d=1, na
                          h40(a,b,c,d) = eri(a,b,c,d)*x(a)*x(b)*x(c)*x(d)/4.0
                          h04(a,b,c,d) = eri(c,d,a,b)*y(a)*y(b)*y(c)*y(d)/4.0
                          h31(a,b,c,d) = -eri(a,b,c,d)*x(a)*x(b)*y(c)*x(d)/2.0
                          h13(a,b,c,d) = -eri(a,d,b,c)*x(a)*y(b)*y(c)*y(d)/2.0
                          h221(a,b,c,d) = eri(a,b,c,d)*x(a)*x(b)*y(c)*y(d)/4.0
                          h222(a,b,c,d) = eri(a,d,b,c)*x(a)*x(c)*y(b)*y(d)
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp single
          tau0 = Sum(h11*t1)
          !$omp end single nowait

          !$omp single
          tau1 = Sum(h221*t2)
          !$omp end single

          !$omp single
          tau4 = tau0 + tau1
          
          allocate(tau2(1:na))

          do a=1, na
              scr1(a) = t1(a,a)
          end do
          
          tau2 = y**2 - x*y*scr1
          R0 = - (tau4/2) + Sum(e0*tau2)/2 - h0/2

          deallocate(tau2)

          tau3 = Sum(e0*x*y*scr1)
          tau4 = tau4 + tau3
          R1 = tau4*t1/2
          R2 = tau4*t2/2

          allocate(tau5(1:na, 1:na, 1:na, 1:na))
          allocate(tau7(1:na, 1:na, 1:na, 1:na))
          allocate(tau8(1:na, 1:na, 1:na, 1:na))
          allocate(tau9(1:na, 1:na, 1:na, 1:na))
          allocate(tau10(1:na, 1:na, 1:na, 1:na))
          allocate(tau11(1:na, 1:na, 1:na, 1:na))

          tau5 = 0.0_pr
          tau7 = 0.0_pr
          tau8 = 0.0_pr
          tau9 = 0.0_pr
          tau10 = 0.0_pr
          tau11 = 0.0_pr
          Amat = 0.0_pr
          !$omp end single

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  Amat(:,:,a,b) = t2(:,a,:,b)
                  do c=1, na
                      do d=1, na
                          tau5(a,b,c,d) = tau5(a,b,c,d) + h11(a,b)*t1(c,d)
                          tau11(a,b,c,d) = delK(a,d)*t1(c,b)-delK(a,b)*t1(c,d)
                      end do
                  end do
              end do
          end do
          !$omp end do

          !$omp single
          Call A_dot_B(h222,Amat,tau5,na**2,na**2,na**2)
          !$omp end single nowait

          !$omp single
          Call A_dot_Btran(h13,t1,tau7,na**3,na,na)
          !$omp end single nowait

          !$omp single
          Call A_dot_B(h31,t1,tau8,na**3,na,na)
          !$omp end single nowait

          !$omp single
          Call Atran_dot_Btran(-h02,t2,tau9,na,na,na**3)
          !$omp end single nowait

          !$omp single
          Call A_dot_B(h20,t2,tau10,na,na,na**3)
          !$omp end single

          !$omp do schedule(static)
          do b=1, na
              do a=1, na
                  R2(a,b,:,:) = R2(a,b,:,:) - (&
                      tau5(a,:,b,:) - Transpose(tau5(a,:,b,:)) -&
                      tau5(b,:,a,:) + Transpose(tau5(b,:,a,:))&
                      )/2 
                  R2(a,b,:,:) = R2(a,b,:,:) + (&
                      tau7(b,:,:,a) - tau7(a,:,:,b)&
                      ) 
                  R2(a,b,:,:) = R2(a,b,:,:) + (&
                      tau8(a,b,:,:) - Transpose(tau8(a,b,:,:))&
                      )
                  R2(a,b,:,:) = R2(a,b,:,:) - (&
                      tau9(:,a,b,:) - Transpose(tau9(:,a,b,:))&
                      )/2
                  R2(a,b,:,:) = R2(a,b,:,:) - (&
                      tau10(a,b,:,:) - tau10(b,a,:,:)&
                      )/2
                  R2(a,b,:,:) = R2(a,b,:,:) - e0(a)*x(a)*y(a)*(&
                      tau11(a,:,b,:))/2 - e0(b)*x(b)*y(b)*(&
                      Transpose(tau11(b,:,a,:))/2)
              end do
          end do
          !$omp end do

          !$omp single
          Call A_dot_B(-t2,h04,R2,na**2,na**2,na**2)
          !$omp end single nowait

          !$omp single
          Call Atran_dot_B(-h40,t2,R2,na**2,na**2,na**2)
          !$omp end single

          !$omp single
          deallocate(tau5,tau7,tau8,tau9,tau10,tau11)
          R1 = R1 - h11/2 - MatMul(t1,h02)/2 - MatMul(Transpose(h20),t1)/2
          !$omp end single

          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      R1(a, b) = R1(a, b) - ( &
                          e0(c) * x(c) * y(c) * t2(c, a, c, b) / 2&
                      )
                  end do
                  R1(a, b) = R1(a, b) + ( &
                      delK(a, b) * e0(a) * x(a) * y(a) / 2&
                  )
                  R1(a,b) = R1(a,b) + Sum(h31(:,:,:,a)*t2(:,:,:,b))/2 &
                      - Sum(h13(:,:,:,b)*t2(:,a,:,:))/2 &
                      - Sum(h11*t2(:,a,:,b))/2 &
                      - Sum(t1*h222(a,b,:,:))/2
              end do
          end do
          !$omp end do
          !$omp end parallel

          R2 = R2 - 2*h221


      End Subroutine BetaCI

      Subroutine NumberCI(T1, T2, NA, X, Y, R0, R1, R2)
          Implicit None
          Integer, parameter  :: pr = Selected_Real_Kind(15,307)
          Integer, Intent(In) :: NA
          Real (Kind=pr), Intent(In) :: X(NA), Y(NA)
          Real (Kind=pr), Intent(In) :: T1(NA,NA)
          Real (Kind=pr), Intent(In) :: T2(NA,NA,NA,NA)
          Real (Kind=pr), Intent(Out) :: R1(NA,NA)
          Real (Kind=pr), Intent(Out) :: R2(NA,NA,NA,NA)
          Real (Kind=pr), Intent(Out) :: R0
          Integer :: p,q,r,s

          Real (Kind=pr)    ::  t1diag(na)

          R0 = 0.0_pr
          R1 = 0.0_pr
          R2 = 0.0_pr

          t1diag = 0.0_pr

          do p=1, na
              t1diag(p) = t1(p,p)
          end do

          R0 = Sum(x*y*t1diag)

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  R1(p,q) = t1(p,q) * ( x(p)**2 - y(q)**2 )/2.0_pr
                  do r=1, na
                      R1(p, q) = R1(p, q) + ( &
                          -t1(r, r)*t1(p, q)*x(r)*y(r) + t2(r, p, r, q)*x(r)*y(r)&
                      )
                      do s=1, na

                          R2(p,q,r,s) = R2(p,q,r,s) + t2(p,q,r,s)*(&
                              x(p)**2 + x(q)**2 - y(r)**2 - y(s)**2)/2.0_pr -&
                              Sum(x*y*t1diag)*t2(p,q,r,s)
          
                      end do
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel

      End Subroutine NumberCI

      Subroutine A_dot_B(A,B,C,MA,NA,NB)
          ! Subroutine performs the following Matrix Multiplication
          !
          !      C = C + A . B
          ! 
          !  where A is of dimension (MA,NA)
          !    and B is of dimension (MB,NB)
          ! Obviously, NB = NA
           Implicit None
           Integer, parameter  :: pr = Selected_Real_Kind(15,307)
           Integer, Intent(In)           ::  MA,NA,NB
           Real (Kind=pr), Intent(In)    ::  A(MA,NA), B(NA,NB)
           Real (Kind=pr), Intent(InOut) ::  C(MA,NB)

           !C = C + MatMul(A,B)
           Call DGEMM('N','N',MA,NB,NA,1.0_pr,A,MA,B,NA,1.0_pr,C,MA)

           Return

      End Subroutine A_dot_B

      Subroutine A_dot_Btran(A,B,C,MA,NA,MB)
          ! Subroutine performs the following Matrix Multiplication
          !
          !      C = C + A . B^transpose
          ! 
          !  where A is of dimension (MA,NA)
          !    and B is of dimension (MB,NB)
          ! Obviously, NB = NA
           Implicit None
           Integer, parameter  :: pr = Selected_Real_Kind(15,307)
           Integer, Intent(In)           ::  MA,NA,MB
           Real (Kind=pr), Intent(In)    ::  A(MA,NA), B(MB,NA)
           Real (Kind=pr), Intent(InOut) ::  C(MA,MB)

           !C = C + MatMul(A,Transpose(B))
           Call DGEMM('N','T',MA,MB,NA,1.0_pr,A,MA,B,MB,1.0_pr,C,MA)

           Return

      End Subroutine A_dot_Btran 

      Subroutine Atran_dot_B(A,B,C,NA,MA,NB)
          ! Subroutine performs the following Matrix Multiplication
          !
          !      C = C + A^transpose . B
          ! 
          !  where A is of dimension (NA,MA)
          !    and B is of dimension (MB,NB)
          ! Obviously, MB = NA
           Implicit None
           Integer, parameter  :: pr = Selected_Real_Kind(15,307)
           Integer, Intent(In)           ::  MA,NA,NB
           Real (Kind=pr), Intent(In)    ::  A(NA,MA), B(NA,NB)
           Real (Kind=pr), Intent(InOut) ::  C(MA,NB)

           !C = C + MatMul(Transpose(A),B)
           Call DGEMM('T','N',MA,NB,NA,1.0_pr,A,NA,B,NA,1.0_pr,C,MA)

           Return

      End Subroutine Atran_dot_B

      Subroutine Atran_dot_Btran(A,B,C,NA,MA,MB)
          ! Subroutine performs the following Matrix Multiplication
          !
          !      C = C + A^transpose . B
          ! 
          !  where A is of dimension (NA,MA)
          !    and B is of dimension (MB,NB)
          ! Obviously, NB = NA
           Implicit None
           Integer, parameter  :: pr = Selected_Real_Kind(15,307)
           Integer, Intent(In)           ::  MA,NA,MB
           Real (Kind=pr), Intent(In)    ::  A(NA,MA), B(MB,NA)
           Real (Kind=pr), Intent(InOut) ::  C(MA,MB)

           !C = C + MatMul(Transpose(A),Transpose(B))
           Call DGEMM('T','T',MA,MB,NA,1.0_pr,A,NA,B,MB,1.0_pr,C,MA)

           Return

      End Subroutine Atran_dot_Btran
