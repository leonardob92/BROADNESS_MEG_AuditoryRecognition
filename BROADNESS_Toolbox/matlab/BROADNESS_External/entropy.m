function e=entropy(x)
%ENTROPY   Entropy of a distribution.
%    E=ENTROPY(X) computes the entropy of the
%    distribution X.
%
%    See also MI.

% Copyright (c) 2008-
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 1999-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2016/03/03 14:57:40 $
% $Revision: 2.6 $
%
% $Log: entropy.m,v $
% Revision 2.6  2016/03/03 14:57:40  marwan
% updated input/output check (because Mathworks is removing downwards compatibility)
% bug in crqad_big fixed (calculation of xcf).
%
% Revision 2.5  2014/09/23 07:11:25  marwan
% Matlab R2014b incompatibility issues fixed
%
% Revision 2.4  2013/01/15 09:08:01  marwan
% bug fix for empty vector
%
% Revision 2.3  2009/03/24 08:32:09  marwan
% copyright address changed
%
% Revision 2.2  2005/04/06 12:58:08  marwan
% changed from log2 to natural log
%
% Revision 2.1  2004/11/10 07:07:51  marwan
% initial import
%
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

narginchk(1,1)
nargoutchk(0,1)

if isempty(x)
    e = 0;
else
    for j=1:size(x,2);
      x(:,j)=x(:,j)./sum(x(:,j));
      x2=x(find(x(:,j)),j);
      e(:,j)=sum(-x2.*log(x2));
    end
end

