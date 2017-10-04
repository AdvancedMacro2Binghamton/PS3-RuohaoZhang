% PROGRAM NAME: ps4huggett.m
clear, clc

% PARAMETERS
beta = .9932; %discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)
y_s = [1, b]; % endowment in employment states
PI = [.97 .03; .5 .5]; % transition matrix


% ASSET VECTOR
a_lo = -2; %lower bound of grid points
a_hi = 5;%upper bound of grid points
num_a = 500;

a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR q
q_min = 0.98;
q_max = 1;

% ITERATE OVER ASSET PRICES
aggsav = 1 ;
while abs(aggsav) >= 0.01 ;
    q_guess = (q_min + q_max) / 2;
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, a', q_guess * a);
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons<0)=-Inf;
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(2, num_a);
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol >.0001;
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        E=repmat(PI(1,:)*v_guess, [num_a 1]);
        E(:,:,2)=repmat(PI(2,:)*v_guess, [num_a 1]);
        v=ret+beta*E;
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vmax, pol_indx]=max(v, [], 2);
        v_tol = [max(abs(vmax(:,:,1)' - v_guess(1,:))) ; max(abs(vmax(:,:,2)' - v_guess(2,:)))];
        v_guess=[vmax(:,:,1)';vmax(:,:,2)'];
    end;
    
    % KEEP DECSISION RULE
    pol_indx=permute(pol_indx, [3 1 2]);
    pol_fn = a(pol_indx);
    
    % SET UP INITITAL DISTRIBUTION
    Mu=ones(2, num_a)/(2*num_a);
    % ITERATE OVER DISTRIBUTIONS
    dif=1;
    while dif >.00000001;
          [emp_ind, a_ind, mass] = find(Mu > 0); % find non-zero indices
          MuNew = zeros(size(Mu));
          for ii = 1:length(emp_ind)
              apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); % which a prime does the policy fn prescribe?
              MuNew(:, apr_ind) = MuNew(:, apr_ind) + [PI(emp_ind(ii),1)*Mu(emp_ind(ii),a_ind(ii));PI(emp_ind(ii),2)*Mu(emp_ind(ii),a_ind(ii))];% which mass of households goes to which exogenous state?
          end
          dif=max(max(abs(MuNew-Mu)));
          Mu=MuNew;
    end
    aggsav=Mu(1,:)*a'+Mu(2,:)*a';   
    if aggsav>=0;
       q_min=q_guess;
    else
       q_max=q_guess;
    end
    q_min
    q_max
end
subplot(2,1,1)
plot(a,v_guess)
legend('Employeed','Unemployeed','location','northwest')
title(['Value Function'])
subplot(2,1,2)
plot(a,pol_fn)
legend('Employeed','Unemployeed','location','northwest')
title(['Policy Function'])
Mu=Mu';
pop=[Mu(:,1);Mu(:,2)];
wealth=[a'+y_s(1);a'+y_s(2)];
earning=[repmat(y_s(1),num_a,1);repmat(y_s(2),num_a,1)];
wealth(wealth<0)=0;
subplot(1,2,1)
c_w=gini(pop, wealth,true);
title(['Wealth, gini=',num2str(c_w)])
subplot(1,2,2)
c_E=gini(pop, earning,true);
title(['Earning, gini=',num2str(c_E)])
x=PI(2,1)/(1-PI(1,1)+PI(2,1));
Y=x*y_s(1)+(1-x)*y_s(2);
W=[Y^(1-sigma)/(1-sigma)]/(1-beta);
w=repmat(W,2,num_a);
A=(w./v_guess)';
u=sum(Mu(A<1));
W=[Y^(1-sigma)/(1-sigma)]/(1-beta);
w=repmat(W,2,num_a);
lambda=(w./v_guess).^(1/(1-sigma))-1;
plot(a,lambda')
legend('Employeed','Unemployeed','location','northwest')
title(['consumption equivalent'])
sum(sum(lambda.*Mu'))