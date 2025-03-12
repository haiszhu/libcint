function Vijkl = computeVijkl(nd, Norb, idcoefs, Vmunu, Vijkl)
%

collocation_matrix = zeros(nd,Norb,Norb);
tmpidx = 0;
for i=1:Norb
  for j=1:i
    tmpidx = tmpidx + 1;
    collocation_matrix(:,j,i) = idcoefs(:,tmpidx);
  end
end
for i=1:Norb
  for j=i+1:Norb
    collocation_matrix(:,j,i) = collocation_matrix(:,i,j);
  end
end
collocation_matrix = reshape(collocation_matrix,[nd Norb^2]);
% Vijkl = zeros(Norb,Norb,Norb,Norb);
% for i = 1:Norb
%   for j = 1:Norb
%     % collocation_matrix_ij = collocation_matrix(i,:).*collocation_matrix(j,:);
%     collocation_matrix_ij = collocation_matrix(:,(i-1)*Norb+j);
%     for k = 1:Norb
%       for l = 1:Norb
%         % collocation_matrix_kl = collocation_matrix(k,:).*collocation_matrix(l,:);
%         collocation_matrix_kl = collocation_matrix(:,(k-1)*Norb+l);
%         % Vijkl(i,j,k,l) = sum(Vmunu.*( collocation_matrix_ij(:) ...
%         %                              .*collocation_matrix_kl(:)'),'all');
%         % tensor_prod = collocation_matrix_ij(:).*collocation_matrix_kl(:)';
%         % Vijkl(i,j,k,l) = sum(Vmunu.*tensor_prod,'all');
%         Vijkl(i,j,k,l) = collocation_matrix_kl(:)'*(Vmunu*collocation_matrix_ij(:));
%         % if abs(Vijkl(i,j,k,l) - Vijkl2(i,j,k,l)) > 1e-6
%         %   keyboard
%         % end
%       end
%     end
%   end
% end
Vijkl = zeros(Norb,Norb,Norb,Norb);
for i = 1:Norb
  collocation_matrix_i = collocation_matrix(:,(i-1)*Norb+(1:Norb));
  for k = 1:Norb
    collocation_matrix_k = transpose(collocation_matrix(:,(k-1)*Norb+(1:Norb)));
    Vijkl(i,:,k,:) = (collocation_matrix_k*(Vmunu*collocation_matrix_i))';
  end
end
% keyboard
end