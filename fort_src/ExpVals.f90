       Subroutine EvalNum(S1, S2, Z1, Z2, NA, X, Y, NExp)
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
      
       End Subroutine EvalNum

