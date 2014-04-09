
function [bool varargout] = isAdaptive(C0,Cpeak,Css,IO,deltaIO)
% BOOL = isAdaptive(C0,Cpeak,Css,IO,deltaIO)
%
% Returns a boolean for adaptivity. Input parameter "out" is a 4-entry
% vector containing [C_0, C_max, peak_time, C_ss] for the simulation
% following simulation to steady state. 
%
% PARAMS:
%   C0 = initial C concentration
%   Cpeak = concentration of C at peak
%   Css = concentration of C at steady state
%   IO = input added to get to steady state
%   deltaIO = second input - original input (difference)

Z = deltaIO./IO;
S = log10((abs(Cpeak-C0)./C0)./Z);
P = log10(Z.*C0)-log10(abs(Css-C0));

if S>-0.5 && P>1
    bool = true;
else
    bool =false;
end
varargout{1} = S;
varargout{2} = P;