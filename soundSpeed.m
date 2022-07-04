function [ speed ] = soundSpeed( temperature, method )
% soundSpeed calculates the velocity of sound
%
% Name:
%   soundSpeed
% 
% Synopsis:
%   speed = soundSpeed( temperature, method )
% 
% Description:
%   The function calculates the velocity of sound for the
%   given temperature.
%
%   Die Schallgeschwindigkeit wird fuer die uebergebene 
%   Temperatur berechnet und zurueckgegeben. (Formel siehe 
%   Gleichung (3.4) in der Diplomarbeit von Jan Wuerfel S. 61)
% 
%   arguments:
%    - temperature (double) the temerature for which the 
%                           velocity will be calculated
%    - method (String) the name of  the method which is used
%                      valid values are 'JAN','MADER' or 'MARCZAK'
%                      default if no method is given is 'MARCZAK'
%
%   return value:
%    - speed (double) the calculated velocity of sound
% 
if (nargin < 2)
    method = 'MARCZAK';
end

switch lower(method)
   
case 'mader'
   k = [3.14643e-009 -1.478e-006 0.000334199 -0.0580852 5.03711 1402.39];
   speed = polyval(k,temperature);
   
case 'marczak'
   k = [2.78786e-9 -1.398845e-6 3.287156e-4 -5.799136e-2 5.038813 1.402385e3];
   speed = polyval(k,temperature);
   
case 'jan'
   speed = 1557 - 0.0245 .* (74 - temperature).^2;

case 'directthrough'
   speed = temperature;
   
   
otherwise
    error(['Unknown method to compute sound speed: ' method]);
end
