function ypred=bestTraining_56(x)
%BESTTRAINING_56 This model file was automatically generated by the GPTIPS 2 function gpmodel2mfile on 20-Dec-2021 16:27:46

ypred=(3818342324243663.*(x(:,11).^9).^(1./2))./17592186044416 ...
    - (8367714107802115.*x(:,16))./9007199254740992 ...
    - (4684844483679697.*sin(sin(x(:,2).*x(:,12))))./1125899906842624 ...
    - (4831284340333437.*sin((x(:,10).^2.*x(:,26))./cos(x(:,15))))./4503599627370496 ...
    - (3387877962052739.*exp(exp(x(:,11).^3)))./70368744177664 ...
    - (8649112683782109.*x(:,15))./18014398509481984 ...
    + 2421873864284579./(35184372088832.*cos(x(:,11))) ...
    + (188607530293811.*cot(x(:,12)).^(1./4))./562949953421312 ...
    - (12795601803451.*x(:,1).*x(:,11))./8796093022208 ...
    - (3714273495887619.*x(:,14).*cos(x(:,14)).*sin(x(:,28)))./18014398509481984 ...
    + 8842052222550475./140737488355328;