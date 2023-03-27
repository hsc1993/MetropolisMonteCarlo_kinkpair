% rotate_tensor_4th.m
%
% rotate 4th order tensor given rotation matrix Q
%��Cp_ijkl = Q_im Q_jn Q_ks Q_lt C_mnst (sum over m,n,s,t)
%

function Cp = rotate_tensor_4th(C,Q)

Cp = zeros(size(C));

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for m=1:3
                    for n=1:3
                        for s=1:3
                            for t=1:3
                                Cp(i,j,k,l)=Cp(i,j,k,l)+Q(i,m)*Q(j,n)*Q(k,s)*Q(l,t)*C(m,n,s,t);
                            end
                        end
                    end
                end
            end
        end
    end
end
