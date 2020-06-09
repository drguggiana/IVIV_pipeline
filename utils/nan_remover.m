function [a_out,b_out] = nan_remover(a,b)
% remove nans from the input vectors

nan_vector1 = isnan(a);
nan_vector2 = isnan(b);

nan_vector = ~(nan_vector1|nan_vector2);

a_out = a(nan_vector);
b_out = b(nan_vector);
