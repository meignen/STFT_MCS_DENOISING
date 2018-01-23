function s = synth(T,Cs,d)


[na N] = size(T);
nr = size(Cs,1);
s = zeros(N,1);

for j=1:nr
    for k=1:N
        s(k) =s(k)+sum(T(max(1,Cs(j,k)-d):min(na,Cs(j,k)+d),k));
    end
end



end




