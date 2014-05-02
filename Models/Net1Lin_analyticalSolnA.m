function A = Net1Lin_analyticalSolnA(p, input, deltaInput)
% A = Net1Lin_analyticalSolnA(p, input, deltaInput)
% Output: A             is a fxn handle
% Input:  p             is the full parameter vector
%         input         is the initial input (u1)
%         deltaInput    is the change in input (u2-u1)
u1 = input;
u2 = input+deltaInput;
c1 = p(1);
c2 = p(2);
AT = p(8);

t1 = c1.*u1.*AT./(c1.*u1+c2);
t2 = c1.*u2.*AT./(c1.*u2+c2);

A = @(t) (t1-t2).*exp((-c1.*u2-c2).*t)+t2; 