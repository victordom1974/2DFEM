% T2= dissambleT2(T)
%
% For T a set of P2 triangles elements 
% return T2 a set of P1 triangle elements
%
% Specifically, if
% 
% [p1 p2 p3  m23 m31 m12]  
% is a P2 "triangle", with mij the mid-point of pi pj, 
% returns four new triangles
%
% [p1  m12 m31]
% [m12 p2  m23]
% [m13 m23  p3]
% [m12 m23 m13]

function T = disassembleT2(T2) 
T=[T2(:,[1 6 5]);...
   T2(:,[6 2 4]);...
   T2(:,[5 4 3]);...
   T2(:,[6 4 5])];
   
