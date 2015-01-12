function [Tabdel,FreqEqu] = DeltaShift(Mastercurve,CWLF,Texp,Tref,TMes,Shi)

%This function puts in Tabdel the data to trace tan(delta) with respect to
%Temperature, at the frenquency matching damping peak at Texp

global MEST;
global NbMe;


%%%%Reference Frequency research%%%%
MaxDel = Mastercurve(1,5);
FreqMax = Mastercurve(1,1);

for i = 1:NbMe
	if Mastercurve(i,5) > MaxDel
		MaxDel = Mastercurve(i,5);
		FreqMax = Mastercurve(i,1);
	end
end

WCurve = @(T) (CWLF(1)*(Tref-T))/(CWLF(2)+T-Tref);
FreqMax;
Shifting = pchip(TMes,Shi);
%Shiftexp = ppval(Shifting,Texp)
Shiftexp = exp(1)^(WCurve(Texp));
FreqEqu = FreqMax./Shiftexp

%%%%Curve creation%%%%
Tabdel = zeros(NbMe,3);
TDisc = 10:1:1000;
YDisc = zeros(length(TDisc),1);
for i=1:length(TDisc)
	YDisc(i) = WCurve(TDisc(i));
end

for i = 1:NbMe
	Obj = log(Mastercurve(i,1)./FreqEqu);
	mini = abs(YDisc(1)-Obj);
	Indic = 1;
	for j = 1:length(TDisc)
		if abs(YDisc(j)-Obj) < mini
			mini = abs(YDisc(j)-Obj);
			Indic = j;
		end
	end
	Tabdel(i,1) = TDisc(Indic);
	Tabdel(i,2) = Mastercurve(i,5);
	Tabdel(i,3) = exp(-1.10*Tabdel(i,2));
	%disp(['Point ',num2str(i),' on ',num2str(NbMe),' done.']);
end

Tabdel = sortrows(Tabdel);

figure;
plot(Tabdel(:,1),Tabdel(:,2));
xlabel('Temperature');
ylabel('Damping factor');
title(['Damping factor with respect to temperature at FreqEqu = ',num2str(FreqEqu),'Hz']);
xlim([25 150]);
hold off;

Data9010 = dlmread(['Temp-Rest-Std.txt'],'\t');

Verif = figure;
errorbar(Data9010(:,1),Data9010(:,2),Data9010(:,3));
hold all;
scatter(Tabdel(:,1),Tabdel(:,3),'r+');
%title(['Verification of COR']);
xlabel('Temperature (Celsius)','FontSize',15);
ylabel('COR / exp(-1.10*tan \delta)','FontSize',15);
print(Verif,'-depsc','-r300','/run/media/apaloo/Dropbox/Stage M1/Images/Resultats/Anton-9010-1-All/VerifCORDelta.eps');
hold off;


end
