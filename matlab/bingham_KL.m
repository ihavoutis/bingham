function d = bingham_KL(B1,B2)
% KL divergence between two binghams
d = bingham_cross_entropy(B1,B2) - bingham_entropy(B1);
