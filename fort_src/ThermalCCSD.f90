      Subroutine BetaCC(E0, OneH, ERI, S1, S2, NA, X, Y, R0, R1, R2)
          Implicit None
          Integer, parameter  :: pr = Selected_Real_Kind(15,307)
          Integer, Intent(In) :: NA
          Real (Kind=pr), Intent(In) :: E0(Na), OneH(Na,Na), ERI(NA,NA,NA,NA)
          Real (Kind=pr), Intent(In) :: X(NA), Y(NA)
          Real (Kind=pr), Intent(In) :: S1(NA,NA)
          Real (Kind=pr), Intent(In) :: S2(NA,NA,NA,NA)
          Real (Kind=pr), Intent(Out) :: R0
          Real (Kind=pr), Intent(Out) :: R1(NA,NA)
          Real (Kind=pr), Intent(Out) :: R2(NA,NA,NA,NA)
          Integer :: a, b, c, d

          Real (Kind=pr), dimension(:, :), allocatable :: tau0
          Real (Kind=pr), dimension(:), allocatable :: tau1
          Real (Kind=pr), dimension(:, :), allocatable :: tau4
          Real (Kind=pr), dimension(:, :), allocatable :: tau5
          Real (Kind=pr), dimension(:, :), allocatable :: tau6
          Real (Kind=pr), dimension(:, :), allocatable :: tau8
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau11
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau12
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau13
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau14
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau15
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau16
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau17
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau18
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau19
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau20
          Real (Kind=pr), dimension(:, :), allocatable :: tau21
          Real (Kind=pr), dimension(:, :), allocatable :: tau22
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau23

          ! Pre Processing
          ! Defining the Hamiltonian Matrix Elements

          Real (Kind=pr) ::  h0, h20(na,na), h02(na,na), h11(na,na)
          Real (Kind=pr) ::  h40(na,na,na,na), h04(na,na,na,na)
          Real (Kind=pr) ::  h31(na,na,na,na), h13(na,na,na,na)
          Real (Kind=pr) ::  h221(na,na,na,na)
          Real (Kind=pr) ::  h222(na,na,na,na)
          Real (Kind=pr) ::  scr1(na), scr2(na,na)
          Real (Kind=pr) ::  scr4(na,na,na,na)

          h0 = 0.0_pr

          do a=1, na
              h0 = h0 + y(a)*y(a)*oneh(a,a)
              do b=1, na
                  scr2(a,b) = 0.0
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
          !$omp end parallel

          allocate(tau0(1:na, 1:na))
          allocate(tau1(1:na))
          allocate(tau4(1:na, 1:na))
          allocate(tau6(1:na, 1:na))
          allocate(tau21(1:na, 1:na))
          
          tau0 = h11
          tau1 = y**2

          do a=1, na
              tau1(a) = tau1(a) - x(a)*y(a)*s1(a,a)
              do b=1, na
                  tau0(a,b) = tau0(a,b) + 2.0_pr*Sum(s1*h221(a,:,b,:))
                  tau6(a,b) = 2.0_pr*Sum(s1*h13(:,:,a,b))
                  tau4(a,b) = 4.0_pr*Sum(s1*h221(:,a,:,b))
              end do
          end do

          r0 = ( Sum(e0*tau1) - Sum(s1*tau0) )/2.0_pr

          tau21 = -tau6

          do a=1, na
              do b=1, na
                  r1(a,b) = Sum(tau4*s2(a,:,:,b))/2.0_pr
              end do
          end do

          deallocate(tau0, tau1)
          allocate(tau5(1:na, 1:na))
          allocate(tau22(1:na, 1:na))

          tau4 = tau4 + h11
          tau5 = MatMul(Transpose(tau4),s1)
          tau6 = tau6 - tau5 + h02
          tau21 = tau21 + tau5
          tau22 = MatMul(tau4,Transpose(s1))

          do a=1, na
              do b=1, na
                  tau6(a,b) = tau6(a,b) - ( &
                      e0(a) * x(a) * y(a) * s1(a, b)&
                      ) + 2.0_pr*Sum(h221(:,:,a,:)*s2(:,:,:,b))
              end do
          end do

          r1 = r1 - MatMul(s1,tau6)/2.0_pr

          deallocate(tau4, tau5, tau6)
          allocate(tau8(1:na, 1:na))

          scr4 = 0.0_pr

          do a=1, na
              do b=1, na
                  tau8(a,b) = -2.0_pr*Sum(s1*h31(:,a,:,b))
              end do
          end do

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          scr4(a,b,c,d) = Sum(&
                              tau8(:,b)*s2(:,a,c,d) &
                              )/4.0_pr
                      end do
                  end do
              end do
          end do

          tau8 = tau8 + h20

          do a=1, na
              do b=1, na
                  tau8(a,b) = tau8(a,b) + 2.0_pr * Sum(&
                      h221(a,:,:,:)*s2(:,b,:,:))
              end do
          end do

          r1 = r1 - MatMul(Transpose(tau8),s1)/2.0_pr

          deallocate(tau8)
          allocate(tau11(1:na, 1:na, 1:na, 1:na))
          allocate(tau12(1:na, 1:na, 1:na, 1:na))
          allocate(tau13(1:na, 1:na, 1:na, 1:na))
          allocate(tau14(1:na, 1:na, 1:na, 1:na))
          allocate(tau15(1:na, 1:na, 1:na, 1:na))
          allocate(tau16(1:na, 1:na, 1:na, 1:na))
          allocate(tau17(1:na, 1:na, 1:na, 1:na))
          allocate(tau18(1:na, 1:na, 1:na, 1:na))
          allocate(tau19(1:na, 1:na, 1:na, 1:na))
          allocate(tau23(1:na, 1:na, 1:na, 1:na))

          tau11 = 0.0_pr
          tau12 = 0.0_pr
          tau13 = 0.0_pr
          tau14 = 0.0_pr
          tau15 = 0.0_pr
          tau16 = 0.0_pr
          tau17 = 0.0_pr
          tau18 = 0.0_pr
          tau19 = 0.0_pr
          tau23 = 0.0_pr

          tau14 = h13
          tau16 = h04

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau13(a,b,c,d) = Sum(&
                              h221(a,b,:,:)*s2(c,d,:,:))
                          tau11(a,b,c,d) = Sum(&
                              s1(d,:)*h221(a,b,c,:))
                          tau14(a,b,c,d) = tau14(a,b,c,d) + &
                              2.0_pr * Sum(s1(:,d)*h221(:,a,b,c))
                          tau16(a,b,c,d) = tau16(a,b,c,d) + &
                              Sum(s1(:,d)*h13(:,a,b,c))
                          tau17(a,b,c,d) = Sum(&
                              s1(:,d)*h31(:,a,b,c)) + Sum(&
                              h221(a,:,b,:)*s2(:,c,:,d))
                          tau19(a,b,c,d) = - h222(d,c,a,b) + &
                              Sum(s1(d,:)*h13(a,:,b,c))
                          tau23(a,b,c,d) = h13(d,b,c,a) + Sum(&
                              s1(d,:)*h04(:,a,b,c))
                      end do
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel

          do a=1, na
              do b=1, na
                  tau18(a,b,:,:) = Transpose(tau13(a,b,:,:))
                  tau21(a,b) = tau21(a,b) + 2.0_pr*Sum(&
                      h221(:,:,:,a)*s2(:,:,:,b))
                  tau22(a,b) = tau22(a,b) + 2.0_pr*Sum(&
                      h221(:,a,:,:)*s2(:,b,:,:))
              end do
          end do

          tau11 = tau11 + h31
          tau13 = tau13 + 2.0_pr * h40
          tau21 = tau21 - h02
          tau22 = tau22 - h20

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau12(a,b,c,d) = Sum(&
                              s1(d,:)*tau11(a,b,:,c))
                          tau15(a,b,c,d) = 2.0_pr*( &
                              Sum(s1(c,:)*tau14(a,:,b,d))&
                              ) - h222(c,d,a,b)
                      end do
                  end do
                  tau13(a,b,:,:) = tau13(a,b,:,:) + &
                      2.0_pr * Transpose(tau12(a,b,:,:))
                  tau18(a,b,:,:) = tau18(a,b,:,:) + &
                      2.0_pr * tau12(a,b,:,:) + &
                      2.0_pr * Transpose(h40(a,b,:,:))
              end do
          end do
          !$omp end do
          !$omp end parallel

          allocate(tau20(1:na, 1:na, 1:na, 1:na))

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          scr4(a,b,c,d) = scr4(a,b,c,d) - &
                              Sum(s2(:,:,c,d)*tau13(:,:,a,b))/8.0_pr + &
                              Sum(s2(:,b,:,d)*tau15(:,:,a,c))/2.0_pr - &
                              Sum(s2(a,b,:,:)*tau16(:,:,c,d))/4.0_pr + &
                              Sum(s2(:,a,:,d)*tau17(:,:,b,c)) + &
                              Sum(tau21(:,c)*s2(a,b,:,d))/4.0_pr + &
                              Sum(tau22(:,a)*s2(:,b,c,d))/4.0_pr
                          tau20(a,b,c,d) = Sum(&
                              s1(:,b)*tau18(:,c,d,a)) + 2.0_pr*Sum(&
                              s1(d,:)*tau19(c,:,b,a))
                      end do
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          scr4(a,b,c,d) = scr4(a,b,c,d) - Sum(&
                              s1(:,d)*tau20(a,c,:,b))/4.0_pr - Sum(&
                              s1(b,:)*tau23(:,c,d,a))/2.0_pr - Sum(&
                              s1(:,c)*h31(a,b,d,:))/2.0_pr + (&
                              Sum(e0*x*y*s1(a,:)*s2(:,b,c,d)) + &
                              Sum(e0*x*y*s1(:,c)*s2(a,b,:,d)) &
                              )/4.0_pr

                      end do
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel

          deallocate(&
              tau11, tau12, tau13, tau14, tau15, tau16, tau17, tau18, &
              tau19, tau20, tau21, tau22, tau23 &
              )

          r0 = r0 - (h0 + Sum(h221*s2))/2.0_pr
          r1 = r1 - h11/2.0_pr
          scr4 = scr4 - h221/2.0_pr
          r2 = 0.0_pr

          do a=1, na
              r1(a,a) = r1(a,a) + e0(a)*x(a)*y(a)/2
              do b=1, na
                  r1(a,b) = r1(a,b) + Sum(h11*s2(a,:,:,b))/2.0_pr - &
                      Sum(s1*h222(a,b,:,:))/2.0_pr + &
                      Sum(h31(:,:,:,a)*s2(:,:,:,b))/2.0_pr - &
                      Sum(h13(:,:,:,b)*s2(:,a,:,:))/2.0_pr
                  do c=1, na
                      r1(a, b) = r1(a, b) - ( &
                          e0(c) * x(c) * y(c) * s2(c, a, c, b) / 2.0_pr &
                          )
                      do d=1, na
                          r2(a,b,c,d) = scr4(a,b,c,d) - scr4(b,a,c,d) &
                              - scr4(a,b,d,c) + scr4(b,a,d,c)
                      end do
                  end do
              end do
          end do

      End Subroutine BetaCC

      Subroutine NumberCC(S1, S2, NA, X, Y, R0, R1, R2)
          Implicit None
          Integer, parameter  :: pr = Selected_Real_Kind(15,307)
          Integer, Intent(In) :: NA
          Real (Kind=pr), Intent(In) :: X(NA), Y(NA)
          Real (Kind=pr), Intent(In) :: S1(NA,NA)
          Real (Kind=pr), Intent(In) :: S2(NA,NA,NA,NA)
          Real (Kind=pr), Intent(Out) :: R1(NA,NA)
          Real (Kind=pr), Intent(Out) :: R2(NA,NA,NA,NA)
          Real (Kind=pr), Intent(Out) :: R0
          Integer :: p,q,r,s

          Real (Kind=pr)    ::  s1diag(na)

          r0 = 0.0_pr
          r1 = 0.0_pr
          r2 = 0.0_pr

          s1diag = 0.0_pr

          do p=1, na
              s1diag(p) = s1(p,p)
          end do

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do p=1, na
              do q=1, na
                  r1(p,q) = (x(p)**2 - y(q)**2)*s1(p,q)/2.0_pr &
                      - Sum(s1(:,q)*s1(p,:) * x * y)
                  do r=1, na
                      r1(p, q) = r1(p, q) + s2(r, p, r, q)*x(r)*y(r)
                      do s=1, na
                          r2(p,q,r,s) = s2(p,q,r,s)*( &
                              x(p)**2 + x(q)**2 - y(r)**2 - y(s)**2 &
                              )/2.0_pr + &
                              Sum(x*y*s1(:,s)*s2(p,q,:,r)) - &
                              Sum(x*y*s1(:,r)*s2(p,q,:,s)) - &
                              Sum(x*y*s1(p,:)*s2(:,q,r,s)) + &
                              Sum(x*y*s1(q,:)*s2(:,p,r,s))
                      end do
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel

          r0 = Sum(x*y*s1diag)

     End Subroutine NumberCC
