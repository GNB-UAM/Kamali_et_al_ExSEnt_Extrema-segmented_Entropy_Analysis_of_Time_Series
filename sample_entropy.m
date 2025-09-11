% Computes sample entropy for a 1D data series.

function samp_ent = sample_entropy(data, m, r)
    
    N = length(data);
    if N <= m+1
        samp_ent = NaN;
        return;
    end
    
    % Construct templates of length m.
    X = zeros(N-m+1, m);
    for i = 1:(N-m+1)
        X(i,:) = data(i:i+m-1);
    end
    
    % Count similar template pairs for m
    B = 0;
    for i = 1:size(X,1)
        for j = i+1:size(X,1)
            if max(abs(X(i,:) - X(j,:))) <= r
                B = B + 1;
            end
        end
    end
    
    % Construct templates of length m+1.
    X1 = zeros(N-m, m+1);
    for i = 1:(N-m)
        X1(i,:) = data(i:i+m);
    end
    
    % Count similar template pairs for m+1.
    A = 0;
    for i = 1:size(X1,1)
        for j = i+1:size(X1,1)
            if max(abs(X1(i,:) - X1(j,:))) <= r
                A = A + 1;
            end
        end
    end
    
    if B == 0 || A == 0
        samp_ent = Inf;
    else
        samp_ent = -log(A / B);
    end

end