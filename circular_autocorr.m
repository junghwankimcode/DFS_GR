function R = circular_autocorr(x)
% x: time sequence
    L = length(x);
    R = zeros(1, L);
    for d = 0:L-1
        R(d+1) = sum(x .* conj(circshift(x, -d)));
    end
   
    R = R ;
    %R = ifft(abs(fft(x)).^2 );  % Circular autocorrelation

end
