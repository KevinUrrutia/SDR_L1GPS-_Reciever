function [out, register] = shift(register, feedback, output)
%GPS shift register
%input feedback -> which positions to use as feedback (1 indexed)
%      output -> which positions are used for the output
%output -> output of shift register

%calculate the output 
out = zeros(1, size(output, 2));
for i = 1:length(output)
    out(i) = register(output(i));
end

if (size(out, 2) > 1)
    out = mod(sum(out), 2);
else
    out = out(1);
end

element = zeros(1, size(feedback, 2));
for i = 1:length(feedback)
    element(i) = register(feedback(i));
end
fb = mod(sum(element), 2);

%shift to the right
register = circshift(register,1);
register(1) = fb;

end