      Subroutine BetaCC(E0, ERI, T1, T2, NA, X, Y, R0, R1, R2)
          Implicit None
          Integer, parameter  :: pr = Selected_Real_Kind(15,307)
          Integer, Intent(In) :: NA
          Real (Kind=pr), Intent(In) :: E0(Na), ERI(NA,NA,NA,NA)
          Real (Kind=pr), Intent(In) :: X(NA), Y(NA)
          Real (Kind=pr), Intent(In) :: T1(NA,NA)
          Real (Kind=pr), Intent(In) :: T2(NA,NA,NA,NA)
          Real (Kind=pr), Intent(Out) :: R0
          Real (Kind=pr), Intent(Out) :: R1(NA,NA)
          Real (Kind=pr), Intent(Out) :: R2(NA,NA,NA,NA)
          Integer :: a, b, c, d, i
      
          Real (Kind=pr), dimension(:, :), allocatable :: tau0
          Real (Kind=pr), dimension(:, :), allocatable :: tau1
          Real (Kind=pr), dimension(:), allocatable :: tau2
          Real (Kind=pr), dimension(:), allocatable :: tau3
          Real (Kind=pr), dimension(:, :), allocatable :: tau4
          Real (Kind=pr), dimension(:, :), allocatable :: tau5
          Real (Kind=pr), dimension(:, :), allocatable :: tau6
          Real (Kind=pr), dimension(:, :, :), allocatable :: tau7
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau8
          Real (Kind=pr), dimension(:, :), allocatable :: tau9
          Real (Kind=pr), dimension(:, :), allocatable :: tau10
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau11
          Real (Kind=pr), dimension(:, :), allocatable :: tau13
          Real (Kind=pr), dimension(:, :), allocatable :: tau15
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau16
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau17
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau18
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau20
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau23
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau24
          Real (Kind=pr), dimension(:, :), allocatable :: tau25
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau30
          Real (Kind=pr), dimension(:, :), allocatable :: tau32
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau37
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau38
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau39
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau42
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau44
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau45
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau47
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau49
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau50
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau51
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau52
          Real (Kind=pr), dimension(:, :, :, :, :), allocatable :: tau53
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau54
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau55
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau56
          Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau57

          ! Pre Processing
          ! Defining the Hamiltonian Matrix Elements

          Real (Kind=pr) ::  h0, h20(na,na), h02(na,na), h11(na,na)
          Real (Kind=pr) ::  h40(na,na,na,na), h04(na,na,na,na)
          Real (Kind=pr) ::  h31(na,na,na,na), h13(na,na,na,na)
          Real (Kind=pr) ::  h221(na,na,na,na)
          Real (Kind=pr) ::  h222(na,na,na,na)
          Real (Kind=pr) ::  scr1(na), f(na,na)
          Real (Kind=pr) ::  scr4(na,na,na,na)
          Real (Kind=pr) ::  xx(na,na), yy(na,na), xy(na,na)

          h0 = Sum(y*y*e0)

          do a=1, na
              do b=1, na
                  f(a,b) = 0.0
                  do c=1,na
                      scr1(c) = eri(a,c,b,c)
                  end do
                  f(a,b) = Sum(y*y*scr1)
                  h0 = h0 + ( y(a)**2 * y(b)**2 * eri(a,b,a,b) )/2
              end do
              f(a,a) = f(a,a) + e0(a)
          end do

          do a=1, na
              do b=1, na
                  h20(a,b) = x(a)*x(b)*f(a,b)
                  h02(a,b) = -y(a)*y(b)*f(a,b)
                  h11(a,b) = x(a)*y(b)*f(a,b)
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

          do a=1, na
              do b=1, na
                  xx(a,b) = x(a)*x(b)
                  yy(a,b) = y(a)*y(b)
                  xy(a,b) = x(a)*y(b)
              end do
          end do

          
          allocate(tau0(1:na, 1:na))
          allocate(tau2(1:na))
          allocate(tau1(1:na, 1:na))
          allocate(tau3(1:na))

          tau1 = h11
          scr4 = h221*t2

          do a=1, na
              do b=1, na
                  tau0(a,b) = Sum(scr4(a,b,:,:)*yy)
              end do
          end do

          do a=1, na
              do b=1, na
                  tau1(a,b) = tau1(a,b) + 2*Sum(xy*t1*h221(a,:,b,:))
              end do
          end do

          tau2 = MatMul(Transpose(tau0),x)
          tau2 = tau2 + MatMul(t1*tau1,y)

          tau3 = y**2
          do a=1, na
              tau3(a) = tau3(a) - xy(a,a)**2 * t1(a,a)
          end do

          r0 = -Sum(tau2*x)/2 + Sum(e0*tau3)/2

          deallocate(tau0,tau1,tau2,tau3)

          allocate(tau4(1:na, 1:na))
          allocate(tau30(1:na, 1:na, 1:na, 1:na))
          allocate(tau5(1:na, 1:na))
          allocate(tau6(1:na, 1:na))
          allocate(tau47(1:na, 1:na, 1:na, 1:na))

          do a=1, na
              do b=1, na
                  tau4(a,b) = 0.0_pr
                  tau5(a,b) = Sum(xy**2 * t1 * eri(a,:,b,:))
                  do c=1, na
                      tau4(a,b) = tau4(a,b) + x(c)**2 * Sum(&
                          yy**2 * t2(c,a,:,:) * Transpose(eri(c,b,:,:))&
                          )
                  end do
              end do
          end do

          r1 = -tau4/4
          tau6 = tau5 + f

          do a=1, na
              do b=1, na
                  r1(a,b) = r1(a,b) - Sum(xy**2 * tau6 * t2(:,a,:,b))/2
                  do c=1, na
                      do d=1,na
                          tau30(a,b,c,d) = Sum(x**2 * tau4(b,:) * t2(:,a,c,d))
                          tau47(a,b,c,d) = -2*Sum(&
                              x**2 * tau5(:,d) * t2(:,a,b,c)&
                              )
                      end do
                  end do
              end do
          end do

          deallocate(tau4,tau5)
          allocate(tau25(1:na, 1:na))
          allocate(tau7(1:na, 1:na, 1:na))
          allocate(tau8(1:na, 1:na, 1:na, 1:na))
          allocate(tau9(1:na, 1:na))
          allocate(tau10(1:na, 1:na))
          allocate(tau11(1:na, 1:na, 1:na, 1:na))
          allocate(tau13(1:na, 1:na))

          tau25 = MatMul(Transpose(t1),xy**2 * tau6)
          tau9 = -f

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  tau9(a, b) = tau9(a, b) + ( &
                      y(a)**2 * e0(a) * t1(b, a)&
                  ) + Sum(y**2 * f(a,:)*t1(b,:))
                  do c=1, na
                      tau7(a,c,b) = -Sum(y**2 * t1(a,:) * eri(c,a,b,:))
                      do d=1, na
                          tau8(a,b,c,d) = 2*t1(d,b)*tau7(a,c,b) + &
                              2*t1(a,b)*eri(c,a,d,b) + Sum(&
                              y**2 * t2(a,d,b,:) * eri(a,c,:,b)&
                              )
                          tau11(a,b,c,d) = -Sum(&
                              x**2 * t2(:,a,b,c)*eri(a,:,b,d)) + &
                              2*t1(a,b)*eri(c,a,d,b)
                      end do
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel

          do a=1, na
              do b=1, na
                  tau10(a,b) = -x(a)**2 * Sum(xy**2 * tau8(:,:,a,b)) + &
                      2*x(a)**2 * tau9(a,b)
                  tau13(b,a) = y(a)**2 * Sum(xy**2 * tau11(:,:,b,a))
              end do
          end do

          r1 = r1 + MatMul( Transpose(tau10),t1 )/4.0_pr

          deallocate(tau7,tau8,tau9,tau10,tau11)
          allocate(tau42(1:na, 1:na, 1:na, 1:na))
          allocate(tau15(1:na, 1:na))

          tau15 = 0.0_pr

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau42(a,b,c,d) = -Sum(tau13(d,:)*t2(a,b,:,c))
                      end do
                  end do
              end do
          end do

          do a=1, na
              tau15(a,a) = 1.0_pr
              do b=1, na
                  tau13(a, b) = tau13(a, b) + ( &
                      2 * y(b)**2 * f(b, a)&
                  )
                  tau15(a, b) = tau15(a, b) - ( &
                      y(a)**2 * t1(a, b)&
                  )
              end do
          end do

          r1 = r1 + MatMul(t1,Transpose(tau13))/4.0_pr

          deallocate(tau13)
          allocate(tau16(1:na, 1:na, 1:na, 1:na))

          do a=1, na
              do b=1, na
                  r1(a, b) = r1(a, b) + ( &
                      e0(a) * tau15(a, b) / 2&
                      )
                  do c=1, na
                      do d=1, na
                          tau16(a,b,c,d) = Sum(&
                              y**2 * t1(a,:) * eri(b,:,c,d)&
                              )
                      end do
                  end do
              end do
          end do

          deallocate(tau15)
          allocate(tau17(1:na, 1:na, 1:na, 1:na))
          allocate(tau18(1:na, 1:na, 1:na, 1:na))

          r2 = 0.0_pr

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          r2(a,b,c,d) = Sum(&
                              y**2 * t1(a,:) * tau16(b,:,d,c)&
                              )/2 
                          tau17(a,b,c,d) = -Sum(&
                              xy**2 * t2(:,a,:,b) * eri(d,:,c,:)&
                              )
                          tau18(a,b,c,d) = Sum(&
                              xy**2 * t2(:,a,:,b) * eri(:,c,:,d)&
                              )
                      end do
                  end do
                  r2(a,b,:,:) = r2(a,b,:,:) + (&
                      Transpose(tau16(a,b,:,:)) - Transpose(tau16(b,a,:,:))&
                      )/2
              end do
          end do
          !$omp end do
          !$omp end parallel

          deallocate(tau16)
          allocate(tau44(1:na, 1:na, 1:na, 1:na))
          allocate(tau49(1:na, 1:na, 1:na, 1:na))

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau42(a,b,c,d) = tau42(a,b,c,d) - 2*Sum(&
                              xy**2 * t2(:,b,:,d) * tau18(a,c,:,:))
                          tau44(a,b,c,d) = Sum(&
                              x**2 * t1(:,a) * tau18(b,c,:,d))
                          tau49(a,b,c,d) = Sum(&
                              y**2 * t1(a,:) * tau18(b,c,d,:))
                      end do
                  end do
              end do
          end do

          do a=1, na
              do b=1, na
                  r2(a,b,:,:) = r2(a,b,:,:) + (&
                      tau44(:,a,:,b) - tau44(:,b,:,a) -&
                      Transpose(tau44(:,a,:,b)) + Transpose(tau44(:,b,:,a))&
                      )/2.0_pr - (&
                      tau49(a,b,:,:) - Transpose(tau49(a,b,:,:)) -&
                      tau49(b,a,:,:) + Transpose(tau49(b,a,:,:)) &
                      )/2.0_pr
                  tau18(a,b,:,:) = tau18(a,b,:,:) + Transpose(eri(:,a,:,b))
              end do
          end do

          deallocate(tau44)
          deallocate(tau49)


          allocate(tau20(1:na, 1:na, 1:na, 1:na))

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau20(a,b,c,d) = Sum(&
                              y**2 * t1(a,:) * tau18(b,c,d,:))
                      end do
                  end do
              end do
          end do

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau17(a,b,c,d) = tau17(a,b,c,d) + Sum(&
                              x**2 * t1(:,b) * tau20(a,c,d,:))
                      end do
                  end do
              end do
          end do

          deallocate(tau18,tau20)
          allocate(tau23(1:na, 1:na, 1:na, 1:na))
          allocate(tau24(1:na, 1:na, 1:na, 1:na))

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  r2(a,b,:,:) = r2(a,b,:,:) + (&
                      tau17(a,:,b,:) - Transpose(tau17(a,:,b,:)) -&
                      tau17(b,:,a,:) + Transpose(tau17(b,:,a,:)) &
                      )/2.0_pr
                  tau25(a, b) = tau25(a, b) + ( &
                      y(b)**2 * f(b,a) + &
                      x(b)**2 * y(b)**2 * e0(b) * t1(b, a)&
                  )
                  do c=1, na
                      do d=1, na
                          tau23(a,b,c,d) = Sum(&
                              yy**2 * t2(a,b,:,:) * eri(c,d,:,:))
                      end do
                      do d=1, na
                          tau24(d,a,b,c) = Sum(&
                              x**2 * t1(:,d) * tau23(a,b,c,:))
                          r2(a, b, c, d) = r2(a, b, c, d) + ( &
                              x(c)**2 * e0(c) + x(d)**2 * e0(d) -&
                              y(a)**2 * e0(a) - y(b)**2 * e0(b) &
                              )* t2(a, b, c, d) / 2.0_pr
                      end do
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel

          deallocate(tau17,tau23)

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau24(a,b,c,d) = tau24(a,b,c,d) + 2*Sum(&
                              tau25(a,:)*t2(b,c,:,d))
                      end do
                  end do
              end do
          end do

          deallocate(tau25)
          allocate(tau32(1:na, 1:na))

          tau32 = MatMul(t1,Transpose(tau6*xy**2))

          do a=1, na
              do b=1, na
                  r2(a,b,:,:) = r2(a,b,:,:) + (&
                      tau24(:,a,b,:) -Transpose(tau24(:,a,b,:))&
                      )/4
                  tau32(a,b)  = tau32(a,b) -x(b)**2 * (&
                      f(b, a) - ( &
                      y(b)**2 * e0(b) * t1(a, b)&
                      ))
              end do
          end do

          deallocate(tau24)
          deallocate(tau6)

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau30(a,b,c,d) = tau30(a,b,c,d) + &
                              2*Sum(tau32(a,:)*t2(:,b,c,d))
                      end do
                  end do
              end do
          end do

          deallocate(tau32)
          allocate(tau37(1:na, 1:na, 1:na, 1:na))
          allocate(tau38(1:na, 1:na, 1:na, 1:na))
          allocate(tau39(1:na, 1:na, 1:na, 1:na))

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau37(a,b,c,d) = Sum(&
                              x**2 * t1(:,a) * eri(b,c,d,:))
                          tau38(a,b,c,d) = Sum(&
                              y**2 * t1(a,:) * eri(b,c,d,:))
                      end do
                      do d=1, na
                          tau39(a,d,c,b) = Sum(&
                              y**2 * t1(d,:) * tau38(a,b,c,:))
                      end do
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau42(a,b,c,d) = tau42(a,b,c,d) - 2*Sum(&
                              x**2 * t1(:,c) * tau39(a,b,:,d))
                      end do
                  end do
              end do
          end do

          allocate(tau55(1:na, 1:na, 1:na, 1:na))
          allocate(tau56(1:na, 1:na, 1:na, 1:na))
          allocate(tau50(1:na, 1:na, 1:na, 1:na))

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  r2(a,b,:,:) = r2(a,b,:,:) + (&
                      tau30(a,b,:,:) - tau30(b,a,:,:) &
                      )/4.0_pr - (&
                      tau37(:,b,a,:) - Transpose(tau37(:,b,a,:))&
                      )/2.0_pr + (&
                      tau42(a,b,:,:) - Transpose(tau42(a,b,:,:))&
                      )/4.0_pr
                  tau55(a,b,:,:) = -Transpose(tau39(:,:,a,b))
                  tau56(a,b,:,:) = -Transpose(tau39(a,b,:,:))
                  do c=1, na
                      do d=1, na
                          tau50(a,b,c,d) = Sum(&
                              x**2 * t1(:,b) * tau38(a,c,:,d))
                      end do
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel

          deallocate(tau30,tau37,tau38,tau39,tau42)
          allocate(tau51(1:na, 1:na, 1:na, 1:na))

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau51(a,b,c,d) = -Sum(&
                              x**2 * t1(:,c) * tau50(a,b,:,d))
                      end do
                  end do
              end do
          end do

          deallocate(tau50)
          allocate(tau45(1:na, 1:na, 1:na, 1:na))

          do a=1, na
              do b=1, na
                  r2(a,b,:,:) = r2(a,b,:,:) + (&
                      Transpose(tau51(a,:,:,b) - tau51(b,:,:,a))&
                      )/2.0_pr 
                  do c=1, na
                      do d=1, na
                          tau45(a,b,c,d)  = Sum(&
                              xx**2 * t2(:,:,a,b) * eri(:,:,c,d))
                      end do
                  end do
              end do
          end do

          deallocate(tau51)

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau47(a,b,c,d) = tau47(a,b,c,d) - Sum(&
                              y**2 * t1(a,:) * tau45(b,c,:,d))
                      end do
                  end do
              end do
          end do

          do a=1, na
              do b=1, na
                  r2(a,b,:,:) = r2(a,b,:,:) - (&
                      tau47(a,:,:,b) - tau47(b,:,:,a)&
                      )/4.0_pr
              end do
          end do

          deallocate(tau45,tau47)
          allocate(tau52(1:na, 1:na, 1:na, 1:na))
          allocate(tau53(1:na, 1:na, 1:na, 1:na, 1:na))

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau52(a,b,c,d) = Sum(&
                              x**2 * t1(:,a) * eri(b,:,c,d))
                      end do
                  end do
              end do
          end do

          do a=1, na
              do b=1, na
                  do c=1, na
                      tau53(a,b,c,:,:) = -2*t1(a,b)*Transpose(tau52(c,a,:,:))
                  end do
              end do
          end do

          deallocate(tau52)

          !$omp parallel default(shared)
          !$omp do schedule(static)
          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          do i=1, na
                              tau53(a,b,c,d,i) = tau53(a,b,c,d,i) + Sum(&
                                  x**2 * t2(:,a,b,c) * eri(a,:,i,d))
                          end do
                      end do
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel

          allocate(tau54(1:na, 1:na, 1:na, 1:na))

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau54(a,b,c,d) = -yy(c,d)**2 * Sum(&
                              x**2 * tau53(:,a,b,d,c)) + 2*(&
                              yy(c,d)**2 * eri(d,c,b,a))
                          tau55(a,b,c,d) = tau55(a,b,c,d) + eri(b,a,d,c)
                      end do
                  end do
              end do
          end do

          deallocate(tau53)

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          r2(a,b,c,d) = r2(a,b,c,d) + Sum(&
                              t2(a,b,:,:)*Transpose(tau54(c,d,:,:)))/8.0_pr +&
                              Sum( xx**2 * t2(:,:,c,d) * tau55(:,:,b,a))/4.0_pr
                      end do
                  end do
              end do
          end do

          deallocate(tau54,tau55)
          tau56 = tau56 + eri

          allocate(tau57(1:na, 1:na, 1:na, 1:na))

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          tau57(a,b,c,d) = -Sum(&
                              x**2 * t1(:,a) * tau56(d,c,b,:))
                      end do
                  end do
              end do
          end do

          deallocate(tau56)

          do a=1, na
              do b=1, na
                  do c=1, na
                      do d=1, na
                          r2(a,b,c,d) = r2(a,b,c,d) + Sum(&
                              x**2 * t1(:,c) * tau57(d,:,b,a))/2.0_pr
                      end do
                  end do
              end do
          end do

          deallocate(tau57)

          r0 = r0 - h0/2.0_pr
          r1 = r1 - f/2.0_pr
          
          do a=1, na
              do b=1, na
          
                  r1(a, b) = r1(a, b) + ( &
                      x(b)**2 * e0(b) * t1(a, b) / 2&
                      ) - Sum( xy**2 * t1 * eri(b,:,a,:) )/2.0_pr
                  do c=1, na
                      r1(a,b) = r1(a,b) + x(c)*x(c)*Sum(&
                          xy**2 * t2(c,:,:,b) * eri(:,c,:,a))/4.0_pr
                      scr1(c) = t2(c,a,c,b)
                  end do
                  r1(a,b) = r1(a,b) - Sum(x**2 * y**2 * e0 * scr1)/2.0_pr
              end do
          end do

          r2 = r2 - eri/2.0_pr

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
                  r1(p,q) = -Sum(s1(:,q)*s1(p,:) * x**2 * y**2)
                  do r=1, na
                      r1(p, q) = r1(p, q) + s2(r, p, r, q)*x(r)**2*y(r)**2
                      do s=1, na
                          r2(p,q,r,s) = -Sum(s1(p,:)*s2(:,q,r,s)*x**2 * y**2) &
                              - Sum(s1(q,:)*s2(p,:,r,s)* x**2 * y**2) &
                              - Sum(s1(:,r)*s2(p,q,:,s)* x**2 * y**2) &
                              - Sum(s1(:,s)*s2(p,q,r,:)* x**2 * y**2)
                      end do
                  end do
              end do
          end do
          !$omp end do
          !$omp end parallel

          r0 = Sum(x**2 * y**2 * s1diag)

     End Subroutine NumberCC
