function m = convertNaN2zero(m)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
m(isnan(m)) = 0;
end

