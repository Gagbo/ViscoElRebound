function [Tabdel] = Verif1Hz(Mastercurve,CWLF,FreqEqu,Tref,TMes,Shi)

%This function puts in Tabdel the data to trace tan(delta) with respect to
%Temperature, at the frenquency matching damping peak at Texp

global MEST;
global NbMe;

Balai = dlmread('Tsweep1hz-AP.txt','\t');


%%%%Reference Frequency research%%%%
MaxDel = Mastercurve(1,5);
FreqMax = Mastercurve(1,1);

WCurve = @(T) (CWLF(1)*(Tref-T))/(CWLF(2)+T-Tref);

%%%%Curve creation%%%%
Tabdel = zeros(NbMe,3);
TDisc = 10:.05:1000;
YDisc = zeros(length(TDisc));
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
	Tabdel(i,3) = exp(-1.15*Tabdel(i,2));
	%disp(['Point ',num2str(i),' on ',num2str(NbMe),' done.']);
end

Tabdel = sortrows(Tabdel);

figure;
plot(Tabdel(:,1),Tabdel(:,2));
hold all;
xlabel('Temperature');
ylabel('Damping factor');
title(['Damping factor with respect to temperature at FreqEqu = ',num2str(FreqEqu),'Hz']);
xlim([25 100]);
scatter(Balai(:,1),Balai(:,5),'rx');
hold off;


%%%% This part is already in Delta shift
%Data9010 = dlmread(['Temp-Rest-Std.txt'],'\t');

%figure;
%errorbar(Data9010(:,1),Data9010(:,2),Data9010(:,3));
%hold all;
%scatter(Tabdel(:,1),Tabdel(:,3),'ro');
%title(['Verification of COR']);
%hold off;


end
