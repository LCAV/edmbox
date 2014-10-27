function [x,nroot]=dummy(a, b, c)
%-----------------------------------------------------------------------
% Usage:       [x,nroot]=quadfcn(a, b, c)
%              [x nroot]=quadfcn(a, b, c)
%
% Solve a quadratic equation where a, b, and c are real.
%   a*x*x + b*x + c = 0
%
% Public Variables
%   a, b, c     ... coefficients (input)
%   x           ... two complex solutions (output)
%   nroot       ... number of roots (output)
%
% Programming Note:
%   [x,nroot] does NOT mean it is a matrix, for the dimensions of x & root do not have to match.
%   When this function is called without assigning to two variables, only the first variable "x" shows up.
%     e.g., x=quadfcn(1,2,3)
%           x=quadfcn(1,2,3)+1
%   When this function is called by assigning to two variables, both "x" and "nroot" are assigned (but they can have different dimensions).
%     e.g., [x,nroot]=quadfcn(1,2,3)
%   To extract out an element of the function, perform masking with a dot product.
%     e.g., quadfcn(1,2,3)*[ 1 0 ]' ... extract 1st element
%           quadfcn(1,2,3)*[ 0 1 ]' ... extract 2nd element
%-----------------------------------------------------------------------
% Instructor: Nam Sun Wang
%-----------------------------------------------------------------------

      if (a == 0)
        if (b == 0)
%         We have a non-equation; therefore, we have no valid solution
          nroot = 0;
          x = [];
        else
%         We have a linear equation with 1 root.
          nroot = 1;
          x = -c/b;
        end
      else
%     We have a true quadratic equation.  Apply the quadratic formula to find two roots.
        nroot = 2;
        DD = b*b-4*a*c;
        x(1) = (-b+sqrt(DD))/2/a;
        x(2) = (-b-sqrt(DD))/2/a;
      end
