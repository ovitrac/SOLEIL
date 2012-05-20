function [pos,larg,amp] = monotone(X,type,zero)
% MONOTONE recherche les segments monotones dans un vecteur
%		ex. iX = monotone(X)
%		options : [pos,larg,amp] = monotonie(X,type,zero)
%			entrées :
%				X = vecteur colonne
%				type = '+' (défaut), '-' ou '0'
%				zero = valeur max de dX assimimée 0
%			sorties :
%				pos = position du segment
%				larg = largeur du segment en indices
%				amplitude de variation du segment

% Thermique 1.0 - 30/04/01 - Olivier Vitrac (source : WoodOx 1.22)

if nargin<2, type='+'; end
if nargin<3, zero = 1e5*eps; end

X = X(:);
dX=diff(X);
S = [(dX>zero)-(dX<-zero)];
switch upper(type)
case '+', indS = find(S>0);
case '-', indS = find(S<0);
case '0', indS = find(S==0);
end
dindS = [2;diff(indS)];
if any(indS)
	il = indS(find([2;diff(indS)]>1));
	ir = indS(find([diff(indS);2]>1))+1;
else
   il = []; ir = [];
end

if ~nargout & any(il)
   hold on
   plot(X,'b-'), plot(il,X(il),'ro',ir,X(ir),'ms')
   stem(il,X(ir)-X(il),'gd')
else
   if  nargout>0, pos = il; end
   if  nargout>1, larg = ir-il +1; end
   if  nargout>2, amp = X(ir)-X(il); end
end
