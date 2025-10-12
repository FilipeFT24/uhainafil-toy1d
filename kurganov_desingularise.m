function [uo] = kurganov_desingularise(H, hui) % eq. 2.17
H4_1              = H.^4;
H4_2              = H4_1;
tol1              = eps;
H4_2(H4_2 < tol1) = tol1;
uo                = sqrt(2).*H.*hui./sqrt(H4_1+H4_2);
%                 = sqrt(2).*H.*hui./sqrt(H4_1+max(H4_2, tol1));
end