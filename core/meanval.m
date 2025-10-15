function [X_bar] = meanval(g, X_DOF)
X_bar = sum(X_DOF*g.Wbf', 2)./sum(g.W, 2);
end