function H = bingham_cross_entropy(B1, B2)
% Computes the cross entropy of two binghams

B1 = bingham_stats(B1);
B2 = bingham_stats(B2);

d = B1.d;
V2 = B2.V';
% compute the full, transposed V1
V1_ft = [B1.mu, B1.V];

% rotate B2 into B1's coordinate frame
  V1_ft_inv = inv(V1_ft);
  A = zeros(d-1, d);
  for i = 1 : d-1
	  for j = 1 : d
		  A(i,j) = dot(V1_ft_inv(j,:), V2(i,:));
		  A(i,j) = A(i,j)*A(i,j);
	  end
  end
  
  %  compute H(B1,B2)
  H = log(B2.F);
  for i = 1 : d-1
	  H_i = A(i,1);
	  for j = 2 : d;
		  H_i = H_i + (A(i,j) - A(i,1)) * (B1.dF(j-1)/B1.F);
	  end
	  H_i = H_i * B2.Z(i);
	  H = H - H_i;
  end

end

function B = bingham_stats(B)
% Compute stats for bingham (F, dF, entropy, mode, scatter matrix)
[B.F B.dF] = bingham_F(B.Z);
B.h = bingham_entropy(B); % Compute entropy
B.mu = bingham_mode(B); % Compute mode
B.S = bingham_scatter(B); % Compute the scatter matrix
end
