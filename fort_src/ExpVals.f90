       Subroutine EvalNumber(S1, S2, Z1, Z2, NA, X, Y, NExp)
           Implicit None
      
           Integer, parameter  :: pr = Selected_Real_Kind(15,307)
           Integer, Intent(In) :: NA
           Real (Kind=pr), Intent(In) :: X(NA), Y(NA)
           Real (Kind=pr), Intent(In) :: S1(NA,NA), Z1(NA,NA)
           Real (Kind=pr), Intent(In) :: S2(NA,NA,NA,NA)
           Real (Kind=pr), Intent(In) :: Z2(NA,NA,NA,NA)
           Real (Kind=pr), Intent(Out)   ::  NExp
           Real (Kind=pr)    ::  R1(NA,NA)
           Real (Kind=pr)    ::  R2(NA,NA,NA,NA)
           Real (Kind=pr)    ::  R0
           Integer :: p,q,r,s
      
           Real (Kind=pr)    ::  s1diag(na)
      
           r0 = 0.0_pr
           r1 = 0.0_pr
           r2 = 0.0_pr
      
           s1diag = 0.0_pr
      
           do p=1, na
               s1diag(p) = s1(p,p)
           end do
      
           r0 = Sum(x**2 * y**2 * s1diag) + Sum(y**2)
           
           !$omp parallel default(shared)
           !$omp do schedule(static)
           do p=1, na
               r1(p,p) = x(p)*y(p)
               do q=1, na
                   r1(p,q) = s1(p,q)*x(p)*y(q)*(x(p)**2 - y(q)**2) &
                       - x(p)*y(q)*Sum(s1(:,q)*s1(p,:) * x**2 * y**2)
                   do r=1, na
                       r1(p, q) = r1(p, q) + x(p)*y(q)*s2(r, p, r, q)*x(r)**2*y(r)**2
                       do s=1, na
                           r2(p,q,r,s) = s2(p,q,r,s)*(x(p)**2 + x(q)**2 &
                               - y(r)**2 - y(s)**2) &
                               - Sum(s1(p,:)*s2(:,q,r,s)*x**2 * y**2) &
                               - Sum(s1(q,:)*s2(p,:,r,s)* x**2 * y**2) &
                               - Sum(s1(:,r)*s2(p,q,:,s)* x**2 * y**2) &
                               - Sum(s1(:,s)*s2(p,q,r,:)* x**2 * y**2)
                           r2(p,q,r,s) = r2(p,q,r,s)*x(p)*x(q)*y(r)*y(s)
                       end do
                   end do
               end do
           end do
           !$omp end do
           !$omp end parallel
      
            NExp = r0 + Sum(Z1*r1) + Sum(Z2*r2)
      
       End Subroutine EvalNumber

       Subroutine EvalEnergy(E0, ERI, T1, T2, Z1, Z2, NA, X, Y, ECC)
            Implicit None
            Integer, parameter  :: pr = Selected_Real_Kind(15,307)
            Integer, Intent(In) :: NA
            Real (Kind=pr), Intent(In) :: E0(Na), ERI(NA,NA,NA,NA)
            Real (Kind=pr), Intent(In) :: X(NA), Y(NA)
            Real (Kind=pr), Intent(In) :: T1(NA,NA), Z1(NA,NA)
            Real (Kind=pr), Intent(In) :: T2(NA,NA,NA,NA)
            Real (Kind=pr), Intent(In) :: Z2(NA,NA,NA,NA)
            Real (Kind=pr) :: R0
            Real (Kind=pr) :: R1(NA,NA)
            Real (Kind=pr) :: R2(NA,NA,NA,NA)
            Real (Kind=pr), Intent(Out) :: ECC
            Integer :: a, b, c, d
       

            Real (Kind=pr), dimension(:, :), allocatable :: tau0
            Real (Kind=pr), dimension(:, :), allocatable :: tau1
            Real (Kind=pr), dimension(:, :), allocatable :: tau3
            Real (Kind=pr), dimension(:, :), allocatable :: tau4
            Real (Kind=pr), dimension(:, :), allocatable :: tau5
            Real (Kind=pr), dimension(:, :), allocatable :: tau7
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau8
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau9
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau10
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau11
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau14
            Real (Kind=pr), dimension(:, :), allocatable :: tau17
            Real (Kind=pr), dimension(:, :), allocatable :: tau19
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau21
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau22
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau26
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau27
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau28
            Real (Kind=pr), dimension(:, :), allocatable :: tau31
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau33
            Real (Kind=pr), dimension(:, :), allocatable :: tau34
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau37
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau38
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau39
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau40
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau41
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau42
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau43
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau46
            Real (Kind=pr), dimension(:, :, :, :), allocatable :: tau47

            ! Pre Processing
            ! Defining the Hamiltonian Matrix Elements

            Real (Kind=pr) ::  h0, h20(na,na), h02(na,na), h11(na,na)
            Real (Kind=pr) ::  h40(na,na,na,na), h04(na,na,na,na)
            Real (Kind=pr) ::  h31(na,na,na,na), h13(na,na,na,na)
            Real (Kind=pr) ::  h221(na,na,na,na)
            Real (Kind=pr) ::  h222(na,na,na,na)
            Real (Kind=pr) ::  scr1(na), scr2(na,na)

            h0 = Sum(y*y*e0)

            do a=1, na
                do b=1, na
                    scr2(a,b) = 0.0
                    do c=1,na
                        scr1(c) = eri(a,c,b,c)
                    end do
                    scr2(a,b) = Sum(y*y*scr1)
                    h0 = h0 + ( y(a)**2 * y(b)**2 * eri(a,b,a,b) )/2
                end do
                scr2(a,a) = scr2(a,a) + e0(a)
            end do

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
            
            allocate(tau0(1:na, 1:na))
            allocate(tau1(1:na, 1:na))
            allocate(tau3(1:na, 1:na))
            allocate(tau4(1:na, 1:na))
            allocate(tau5(1:na, 1:na))
            allocate(tau7(1:na, 1:na))
            allocate(tau19(1:na, 1:na))
            allocate(tau31(1:na, 1:na))

            tau0 = h11
            do a=1, na
                do b=1, na
                    tau0(a,b) = tau0(a,b) + 2*Sum(&
                        t1*h221(a,:,b,:))
                    tau1(a,b) = 4*Sum(t1*h221(:,a,:,b))
                    tau3(a,b) = Sum(t1*h31(:,a,:,b))
                    tau7(a,b) = 2*Sum(t1*h13(:,:,a,b))
                end do
            end do
            tau1 = tau1 + h11
            tau4 = MatMul(tau1,Transpose(t1))
            tau5 = tau4 + 2*tau3 - h20
            tau19 = tau4 + 2*tau3
            tau31 = MatMul(Transpose(tau1),t1)

            r0 = Sum(t1*tau0) + h0 + Sum(h221*t2)
            
            deallocate(tau0,tau3,tau4)
            allocate(tau8(1:na, 1:na, 1:na, 1:na))
            allocate(tau9(1:na, 1:na, 1:na, 1:na))
            allocate(tau10(1:na, 1:na, 1:na, 1:na))
            allocate(tau11(1:na, 1:na, 1:na, 1:na))
            allocate(tau37(1:na, 1:na, 1:na, 1:na))
            allocate(tau14(1:na, 1:na, 1:na, 1:na))
            allocate(tau17(1:na, 1:na))
            allocate(tau21(1:na, 1:na, 1:na, 1:na))
            allocate(tau22(1:na, 1:na, 1:na, 1:na))
            allocate(tau26(1:na, 1:na, 1:na, 1:na))
            allocate(tau27(1:na, 1:na, 1:na, 1:na))
            allocate(tau28(1:na, 1:na, 1:na, 1:na))
            allocate(tau33(1:na, 1:na, 1:na, 1:na))
            allocate(tau34(1:na, 1:na))
            allocate(tau38(1:na, 1:na, 1:na, 1:na))
            allocate(tau39(1:na, 1:na, 1:na, 1:na))
            allocate(tau40(1:na, 1:na, 1:na, 1:na))
            allocate(tau41(1:na, 1:na, 1:na, 1:na))
            allocate(tau42(1:na, 1:na, 1:na, 1:na))
            allocate(tau43(1:na, 1:na, 1:na, 1:na))
            allocate(tau46(1:na, 1:na, 1:na, 1:na))
            allocate(tau47(1:na, 1:na, 1:na, 1:na))

            r1 = 0.0_pr

            !$omp parallel default(shared)
            !$omp do schedule(static)
            do a=1, na
                do b=1, na
                    r1(a,b) = Sum(tau1*t2(:,a,:,b)) + Sum(t1*h222(a,b,:,:)) + Sum(&
                        h13(:,:,:,b)*t2(:,a,:,:)) - Sum(&
                        h31(:,:,:,a)*t2(:,:,:,b))
                    tau5(a,b) = tau5(a,b) - 2*Sum(&
                        h221(a,:,:,:) * t2(:,b,:,:))
                    do c=1, na
                        do d=1, na
                            tau8(a,b,c,d) = Sum(t1(:,d)*h222(a,b,:,c))
                            tau10(a,b,c,d) = Sum(&
                                h221(a,:,b,:)*t2(:,c,:,d))
                        end do
                    end do
                end do
            end do
            !$omp end do
            !$omp end parallel

            r1 = r1 - MatMul(Transpose(tau5),t1)
            tau31 = tau31 - tau7 - h02
            tau7 = tau7 + h02
            tau19 = tau19 - h20

            !$omp parallel default(shared)
            !$omp do schedule(static)
            do a=1, na
                do b=1, na
                    tau34(b,a) = Sum(h221(:,:,:,b)*t2(:,:,:,a))
                    tau7(a,b) = tau7(a,b) + 2*Sum(&
                        h221(:,:,a,:)*t2(:,:,:,b))
                    tau17(a,b) = Sum(h221(:,a,:,:)*t2(:,b,:,:))
                    do c=1, na
                        do d=1, na
                            tau9(a,b,c,d) = Sum(&
                                t1(c,:)*tau8(a,b,:,d))
                            tau11(a,b,c,d) = Sum(&
                                t1(:,b)*tau10(:,a,c,d))
                            tau14(a,b,c,d) = 2*Sum(t1(:,d)*h31(:,a,b,c)) +&
                                h222(c,d,a,b)
                            tau21(a,b,c,d) = Sum(h31(:,:,a,b)*t2(:,:,c,d))
                            tau22(a,b,c,d) = Sum(t1(:,d)*h31(a,:,b,c))
                            tau27(a,b,c,d) = Sum(t1(:,d)*h13(:,a,b,c))
                            tau40(a,b,c,d) = Sum(h13(:,a,:,b)*t2(:,c,:,d))
                            tau43(a,b,c,d) = Sum(t1(d,:)*h221(a,b,c,:))
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
                            tau9(a,b,c,d) = tau9(a,b,c,d) + &
                                4*Sum(t1(a,:)*tau11(:,b,c,d)) -&
                                Sum(t2(:,c,:,d)*tau14(:,:,a,b))
                            tau37(a,b,c,d) = 2*Sum(&
                                t2(:,a,:,c)*tau10(:,:,b,d)) + Sum(&
                                tau34(:,d)*t2(a,b,:,c))
                            tau26(a,b,c,d) = -2*Sum(&
                                tau17(:,b)*t2(:,a,c,d)) + Sum(&
                                tau19(:,a)*t2(:,b,c,d))
                            tau21(a,b,c,d) = tau21(a,b,c,d) +&
                                2*Sum(t1(:,c)*tau22(:,a,b,d))
                            tau28(a,b,c,d) = Sum(t1(c,:)*tau27(a,:,b,d))
                            tau38(a,b,c,d) = Sum(h13(a,:,:,b)*t2(c,d,:,:))
                            tau41(a,b,c,d) = Sum(t1(b,:)*tau40(:,a,c,d))
                            tau42(a,b,c,d) = Sum(h221(a,b,:,:)*t2(c,d,:,:)) + &
                                2*Sum(t1(c,:)*tau43(a,b,:,d)) + 2*h40(a,b,c,d)
                        end do
                    end do
                end do
            end do
            !$omp end do
            !$omp end parallel

            r1 = r1 + MatMul(t1,tau7) + h11

            r2 = 4*h221

            !$omp parallel default(shared)
            !$omp do schedule(static)
            do a=1, na
                do b=1, na
                    r2(a,b,:,:) = r2(a,b,:,:)-tau9(a,:,b,:) + Transpose(tau9(a,:,b,:)) &
                        + tau9(b,:,a,:) - Transpose(tau9(b,:,a,:))
                    do c=1, na
                        do d=1, na
                            tau26(a,b,c,d) = tau26(a,b,c,d) + Sum( &
                                t1(b,:) * tau21(:,a,c,d))
                            tau33(a,b,c,d) = 2*Sum(t1(b,:)*tau28(:,a,c,d)) -&
                                Sum(tau31(:,a)*t2(b,c,:,d))
                            tau39(a,b,c,d) = Sum(t1(:,b)*tau38(:,a,c,d))
                            tau46(a,b,c,d) = Sum(t1(:,d)*tau42(:,a,b,c)) - &
                                h31(b,c,d,a)
                            tau47(a,b,c,d) = h13(d,b,c,a) + 2*Sum(&
                                t1(d,:)*h04(:,a,b,c))
                        end do
                    end do
                end do
            end do
            !$omp end do
            !$omp end parallel

            deallocate(tau1,tau5,tau7,tau8,tau9,tau10,tau11,&
                tau14,tau17,tau19,tau21,tau22,tau27,&
                tau28,tau31,tau34,tau38,tau40,tau43)

            !$omp parallel default(shared)
            !$omp do schedule(static)
            do a=1, na
                do b=1, na
                    r2(a,b,:,:) = r2(a,b,:,:) - tau26(a,b,:,:) + &
                        tau26(b,a,:,:) +  tau33(:,a,b,:) - &
                        Transpose(tau33(:,a,b,:)) + 2*tau37(a,b,:,:) - &
                        2*Transpose(tau37(a,b,:,:)) + tau39(:,:,a,b) - &
                        Transpose(tau39(:,:,a,b)) - 2*tau41(:,a,b,:) + &
                        2*tau41(:,b,a,:) + 2*Transpose(tau41(:,a,b,:)) - &
                        2*Transpose(tau41(:,b,a,:))
                    do c=1, na
                        do d=1, na
                            r2(a,b,c,d) = r2(a,b,c,d) + Sum(&
                                t2(:,:,c,d)*tau42(:,:,a,b)) + 2*Sum(&
                                t1(:,d)*tau46(:,a,b,c)) + 2*Sum(&
                                t1(b,:)*tau47(:,c,d,a)) + 2*Sum(&
                                h04(:,:,c,d)*t2(a,b,:,:)) - 2*Sum(&
                                t1(a,:)*h13(b,c,d,:)) + 2*Sum(&
                                t1(:,c)*h31(a,b,d,:))
                        end do
                    end do
                end do
            end do
            !$omp end do
            !$omp end parallel

            deallocate(tau26,tau33,tau37,tau39,tau41,tau42,tau46,tau47)

            ECC = r0 + Sum(Z1*r1) + Sum(Z2*r2)

       End Subroutine EvalEnergy
