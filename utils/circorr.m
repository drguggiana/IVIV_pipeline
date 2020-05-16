function correlation = circorr(a,b)
% calculate circular correlation
% taken from https://www.mathworks.com/matlabcentral/answers/
%230482-how-to-calculate-the-circular-correlation-with-2-sequences-arrays-in-matlab

% remove nans in the inputs
nan_vector1 = isnan(a);
nan_vector2 = isnan(b);

nan_vector = ~(nan_vector1|nan_vector2);

a = a(nan_vector);
b = b(nan_vector);

correlation = ifft(fft(a,5).*conj(fft(b,5)));
