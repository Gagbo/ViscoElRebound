function [File,Mastercurve,CWLF,TMes,Shi] = GetWLFv3(filename,Tref)


%%%%Global variables%%%%
global MEST; %Measures done per temperature in the DMA
MEST = 16;
global NbMe;

%%%%File opening and other variables initialisation%%%%
File = dlmread([filename,'.txt'], '\t', 1,1);

[NbMe,NbVa] = size(File);
Freq = File(:,1); % Frequency of the DMA measure
Temp = File(:,2); % Temperature of the DMA measure
Gpri = File(:,3); % Storage Modulus
Gsec = File(:,4); % Loss Modulus
Delt = File(:,5); % Damping Factor
%Torq = File(:,6); % Torque applied

Mastercurve = zeros(NbMe,7);
% Columns
% Frequency associated with Tref for the measure
% Temperature
% Storage Modulus
% Damping Factor
% Frenquency-shift with respect to Tref measures
% base 10 logarithm of the Shift

%%%%Tref indice search%%%%
Iref = 0;
i = 1;
while (i<NbMe) && (Iref==0)
	if Temp(i) == Tref
		Iref = i
		Irefmod = (MEST - 1 + i)/MEST
	end
	i = i+1;
end

%%%%Mastercurve Initialisation%%%%

for i = 0:MEST-1
	for j = 1:5
		Mastercurve(Iref+i,j) = File(Iref+i,j);
	end
	Mastercurve(Iref+i,6) = 1;
	Mastercurve(Iref+i,7) = 0;
end
disp(['Mastercurve initialised']);

%%%%Right shitfting%%%%

Imod = 0;
while Imod < Irefmod-1
	disp(['Right shift n ',num2str(Imod+1),' on ',num2str(Irefmod-1),' started']);
	StRef = (Irefmod-1-Imod)*MEST + 1;
	[Shft,Datashft]=ShiftCurve(Mastercurve(StRef:StRef+MEST-1,:),File(StRef-MEST:StRef-1,:),1);
	Mastercurve(StRef-MEST:StRef-1,:) = Datashft;
	Mastercurve(StRef-MEST:StRef-1,6) = Mastercurve(StRef:StRef+MEST-1,6).*Mastercurve(StRef-MEST:StRef-1,6);
	Mastercurve(StRef-MEST:StRef-1,7) = log10(Mastercurve(StRef-MEST:StRef-1,6));
	Imod = Imod + 1;
end

disp(['Right shift ended']);

%%%%Left Shitfing%%%%

Imod = 0;
while Imod < NbMe/MEST - Irefmod
	disp(['Left shift n ',num2str(Imod+1),' on ',num2str(NbMe/MEST - Irefmod),' started']);
	StRef = (Irefmod - 1 + Imod)*MEST + 1;
	[Shft,Datashft]=ShiftCurve(Mastercurve(StRef:StRef+MEST-1,:),File(StRef+MEST:StRef+2*MEST-1,:),0);
	Mastercurve(StRef+MEST:StRef+2*MEST-1,:) = Datashft;
	Mastercurve(StRef+MEST:StRef+2*MEST-1,6) = Mastercurve(StRef+MEST:StRef+2*MEST-1,6).*Mastercurve(StRef:StRef+MEST-1,6);
	Mastercurve(StRef+MEST:StRef+2*MEST-1,7) = log10(Mastercurve(StRef+MEST:StRef+2*MEST-1,6));
	Imod = Imod + 1;
end

disp(['Left shift ended']);

%%%%Post-treatment%%%%

i=1;
j=1;
TMes = zeros([NbMe/MEST-1,1]);
WLFY = zeros([NbMe/MEST-1,1]);
WLFX = zeros([NbMe/MEST-1,1]);
LShi = zeros([NbMe/MEST-1,1]);
while i < NbMe + 1
	if File(i,2) == Tref
		i = i + MEST;
		continue;
	else
		TMes(j) = File(i,2);
		i = i + MEST;
		j = j + 1;
	end
end

j=1;
for i = 1:NbMe/MEST-1
	WLFX(i) = 1./(TMes(i)-Tref);
end
for i = 1:NbMe/MEST
	if i == Irefmod
		continue;
	else
		WLFY(j) = 1./Mastercurve((i-1)*MEST+1,7);
		LShi(j) = Mastercurve((i-1)*MEST+1,7);
		Shi(j) = Mastercurve((i-1)*MEST+1,6);
		j = j+1;
	end
end

h1 = figure;
%subplot(2,2,1);
%Curve = scatter(Mastercurve(:,1),Mastercurve(:,3),'r+');
set(gca, 'XScale', 'log','YScale', 'log');
hold all;
%title(['Storage modulus (Pa) with respect to frequency at Tref = ',num2str(Tref)]);

I = 1;

while I < NbMe/MEST
     switch I
        case 1
            rgb = [0 0 1];
        case 2
            rgb = [102./255 0 1];
        case 3
            rgb = [1 0 1];
        case 4
            rgb = [1 0 102./255];
        case 5
            rgb = [1 102./255 0];
        case 6
            rgb = [1 1 0];
        case 7
            rgb = [122./255 204./255 41./255];
        case 8
            rgb = [0 1 0];
        case 9
            rgb = [0 1 1];
        case 10
            rgb = [0 153./255 1];
    end
	CurveI = scatter(Mastercurve(1+(I-1)*MEST:I*MEST,1),Mastercurve(1+(I-1)*MEST:I*MEST,3),'x','MarkerFaceColor',rgb,'MarkerEdgeColor',rgb);
	CurveII = plot(Mastercurve(1+(I-1)*MEST:I*MEST,1),Mastercurve(1+(I-1)*MEST:I*MEST,3),'Color',rgb);
	I = I+1;
end
xlabel('Frequency','FontSize',15);
ylabel('G''','FontSize',15);

hold off;
%print(h1,'-depsc','-r300','/path/to/eps');

h2=figure;
%subplot(2,2,2);
I = 1;

while I < NbMe/MEST
    switch I
        case 1
            rgb = [0 0 1];
        case 2
            rgb = [102./255 0 1];
        case 3
            rgb = [1 0 1];
        case 4
            rgb = [1 0 102./255];
        case 5
            rgb = [1 102./255 0];
        case 6
            rgb = [1 1 0];
        case 7
            rgb = [122./255 204./255 41./255];
        case 8
            rgb = [0 1 0];
        case 9
            rgb = [0 1 1];
        case 10
            rgb = [0 153./255 1];
    end
	Curve3 = scatter(Mastercurve(1+(I-1)*MEST:I*MEST,1),Mastercurve(1+(I-1)*MEST:I*MEST,5),'x','MarkerFaceColor',rgb,'MarkerEdgeColor',rgb);
	hold all;
	I = I+1;
end

%title(['tan(delta) with respect to frequency at Tref = ',num2str(Tref)]);
set(gca,'XScale','log');
xlabel('Frequency','FontSize',15);
ylabel('tan(\delta)','FontSize',15);
hold off;
%print(h2,'-depsc','-r300','/path/to/eps');

h3=figure;
%subplot(2,2,3);
TMes(NbMe/MEST) = Tref;
LShi(NbMe/MEST) = 0;
Shi(NbMe/MEST) = 1;
WLFX
WLFY
WLFYlight = zeros(length(WLFY));
WLFXlight = zeros(length(WLFX));
j = 1;
for i=1:length(WLFY)
	if abs(WLFX(i)) > 0
		WLFXlight(j) = WLFX(i);
		WLFYlight(j) = WLFY(i);
		j=j+1;
	else
		continue
	end
end
WLFYlight(j:length(WLFY)) = [];
WLFXlight(j:length(WLFY)) = [];
[FitL,S] = polyfit(WLFXlight,WLFYlight,1)
%FitLplus = fit(WLFXlight,WLFYlight,'poly1'); %%%%To use with the Curve
%fitting toolbox
%FitL = coeffvalues(FitLplus)
Output = polyval(FitL,WLFXlight);
Correlation2 = corrcoef(WLFYlight,Output)
CWLF(1) = -1./FitL(2)
CWLF(2) = FitL(1)./FitL(2)
CWLFTool = zeros(2,1);
CWLFTool1 = 1/0.2925;
CWLFTool2 = 7.2192/0.2925;
funcWLFTool = @(T) (CWLFTool1*(Tref-T)./(CWLFTool2+(T-Tref)));
funcWLF = @(T) (CWLF(1)*(Tref-T)./(CWLF(2)+(T-Tref)));
funcWLFTheo = @(T) (-17.44*(Tref-T)./(51.6+(T-Tref)));
X = Temp(1):.1:Temp(NbMe);
WLF = funcWLF(X);
WLFTool = funcWLFTool(X);
%CWLF(1)=CWLFTool1;
%CWLF(2) = CWLFTool2;
%Curve4 = plot(X,WLF);
hold all;
Curve4bis = plot(X,WLFTool,'color','r');
scatter(TMes,LShi,'r+');
ylim([-50 50]);
xlim([25 95]);
title('Verification WLF Law');
xlabel('Temperature (Celsius)','FontSize',15);
ylabel('log(shift)','FontSize',15);
hold off;
%print(h3,'-depsc','-r300','/path/to/eps');

h4=figure;
%subplot(2,2,4);
scatter(WLFX,WLFY,'bx');
hold all;
FitY = FitL(1).*WLFX+FitL(2);
%Theo = 51.6/(-17.44).*WLFX-1/17.44;
Curve5 = plot(WLFX,-7.2192.*WLFX-0.2925,'color','b');
%Curve6 = plot(WLFX,Theo,'color','g');
%title('Linear fit for WLF verification');
xlabel('1/(T-Tref)','FontSize',15);
ylabel('-1/log(shift)','FontSize',15);
%text('Units','normalized','Position',[0.5 0.85],'String',char('y = A * x + B', ['A = ',num2str(FitL(1))],['B = ',num2str(FitL(2))]));
%print(h,'-dpng','-r1000',['/path/to/',filename,'/',filename,'-',num2str(Tref),'.png']);
hold off;
%print(h4,'-depsc','-r300','/path/to/eps');

figure;
plot(0:0.25:150,funcWLF(0:0.25:150),'color','b');
hold all;
%plot(0:0.25:150,funcWLFTheo(0:0.25:150),'color','g');
plot(0:0.25:150,funcWLFTool(0:0.25:150),'color','r');
title('WLF Function');
hold off;

Mastercurve = sortrows(Mastercurve);
[Tab1Hz] = Verif1Hz(Mastercurve,CWLF,1,Tref,TMes,Shi);
[Tabdel,FreqEqu] = DeltaShift(Mastercurve,CWLF,84.36,Tref,TMes,Shi);

end

function [Fact,Datashft]=ShiftCurve(DataRef,DataMob,Right)

global MEST;

Datashft = zeros(MEST,7);
Targ = floor(MEST/4);

%DataRef(:,1)
%DataRef(:,3)

Ref = pchip(DataRef(:,1),DataRef(:,3));
Mob = pchip(DataMob(:,1),DataMob(:,3));

if Right == 1
	Obj1 = DichotSolve(DataRef,DataMob(1,3),[DataRef(1,1),DataRef(MEST,1)]);
	Obj11 = DichotSolve(DataRef,DataMob(2,3),[DataRef(1,1),DataRef(MEST,1)]);
    Obj12 = DichotSolve(DataRef,DataMob(3,3),[DataRef(1,1),DataRef(MEST,1)]);
    disp(['Dichotomy 1 done']);
	Obj2 = DichotSolve(DataRef,DataMob(Targ,3),[DataRef(1,1),DataRef(MEST,1)]);
    Obj21 = DichotSolve(DataRef,DataMob(Targ+1,3),[DataRef(1,1),DataRef(MEST,1)]);
    Obj22 = DichotSolve(DataRef,DataMob(Targ+2,3),[DataRef(1,1),DataRef(MEST,1)]);
	disp(['Dichotomy 2 done']);
	Obj3 = DichotSolve(DataRef,DataMob(2*Targ,3),[DataRef(1,1),DataRef(MEST,1)]);
	disp(['Dichotomy 3 done']);

	Obj4 = DichotSolve(DataMob,DataRef(3*Targ,3),[DataRef(1,1),DataRef(MEST,1)]);
    Obj41 = DichotSolve(DataMob,DataRef(3*Targ+1,3),[DataRef(1,1),DataRef(MEST,1)]);
    Obj42 = DichotSolve(DataMob,DataRef(3*Targ+2,3),[DataRef(1,1),DataRef(MEST,1)]);
	disp(['Dichotomy 4 done']);

	Obj5 = DichotSolve(DataMob,DataRef(MEST,3),[DataRef(1,1),DataRef(MEST,1)]);
	disp(['Dichotomy 5 done']);
	Shft1 = Obj1 - DataMob(1,1);
	Fact1 = (DataMob(1,1)+Shft1)/DataMob(1,1)
    Shft11 = Obj11 - DataMob(2,1);
	Fact11 = (DataMob(2,1)+Shft11)/DataMob(2,1)
    Shft12 = Obj12 - DataMob(3,1);
	Fact12 = (DataMob(3,1)+Shft12)/DataMob(3,1)
	Shft2 = Obj2-DataMob(Targ,1);
	Fact2 = (DataMob(Targ,1) + Shft2)/DataMob(Targ,1)
    Shft21 = Obj21-DataMob(Targ+1,1);
	Fact21 = (DataMob(Targ+1,1) + Shft21)/DataMob(Targ+1,1)
    Shft22 = Obj22-DataMob(Targ+2,1);
	Fact22 = (DataMob(Targ+2,1) + Shft22)/DataMob(Targ+2,1)
	Shft3 = Obj3 - DataMob(2*Targ,1);
	Fact3 = (DataMob(2*Targ,1) + Shft3)/DataMob(2*Targ,1)
	Shft4 = DataRef(3*Targ,1) - Obj4;
	Fact4 = DataRef(3*Targ,1)/(DataRef(3*Targ,1)-Shft4)
    Shft41 = DataRef(3*Targ+1,1) - Obj41;
	Fact41 = DataRef(3*Targ+1,1)/(DataRef(3*Targ+1,1)-Shft41)
    Shft42 = DataRef(3*Targ+2,1) - Obj42;
	Fact42 = DataRef(3*Targ+2,1)/(DataRef(3*Targ+2,1)-Shft42)
	Shft5 = DataRef(MEST,1) - Obj5;
	Fact5 = DataRef(MEST,1)/(DataRef(MEST,1)-Shft5)
	%Fact = (Fact1*Fact2*Fact3^(3)*Fact4*Fact5)^(1./7);
	Facttab = [Fact1,Fact11,Fact12,Fact2,Fact21,Fact22,Fact3,Fact4,Fact41,Fact42,Fact5];
	Facttab = sort(Facttab);
	Fact = (Facttab(6)*Facttab(7)*Facttab(8)*Facttab(9)*Facttab(10)*Facttab(11))^(1./6);

else
	Obj1 = DichotSolve(DataRef,DataMob(3*Targ,3),[DataRef(1,1),DataRef(MEST,1)]);
    disp(['Dichotomy 1 done']);
	Obj2 = DichotSolve(DataRef,DataMob(MEST,3),[DataRef(1,1),DataRef(MEST,1)]);
    disp(['Dichotomy 2 done']);
	Obj3 = DichotSolve(DataRef,DataMob(2*Targ,3),[DataRef(1,1),DataRef(MEST,1)]);
    disp(['Dichotomy 3 done']);
	Obj4 = DichotSolve(DataMob,DataRef(Targ,3),[DataRef(1,1),DataRef(MEST,1)]);
    disp(['Dichotomy 4 done']);
	Obj5 = DichotSolve(DataMob,DataRef(1,3),[DataRef(1,1),DataRef(MEST,1)]);
	disp(['Dichotomy done']);
	Shft1 = Obj1 - DataMob(3*Targ,1);
	Fact1 = (DataMob(3*Targ,1) + Shft1)/DataMob(3*Targ,1);
	Shft2 = Obj2 - DataMob(MEST,1);
	Fact2 = (DataMob(MEST,1) + Shft2)/DataMob(MEST,1);
	Shft3 = Obj3 - DataMob(2*Targ,1);
	Fact3 = (DataMob(2*Targ,1) + Shft3)/DataMob(2*Targ,1);
	Shft4 = DataRef(Targ,1) - Obj4;
	Fact4 = DataRef(Targ,1)/(DataRef(Targ,1) - Shft4);
	Shft5 = DataRef(1,1) - Obj5;
	Fact5 = DataRef(1,1)/(DataRef(1,1) - Shft5);
	%Fact = (Fact1*Fact2*Fact3^(3)*Fact4*Fact5)^(1./7);
	Facttab = [Fact1,Fact2,Fact3,Fact4,Fact5];
	Facttab = sort(Facttab);
	Fact = (Facttab(2)*Facttab(1))^(1./2);

	
end

%%%%Filling of the Mastercurve w/ modified frequencies%%%%
Datashft(Targ,1) = DataMob(Targ,1)*Fact; %New frequency
Datashft(Targ,2) = DataMob(Targ,2); %Temperature, Gpri, Gsec and damping factor untouched
Datashft(Targ,3) = DataMob(Targ,3);
Datashft(Targ,4) = DataMob(Targ,4);
Datashft(Targ,5) = DataMob(Targ,5);
Datashft(Targ,6) = Fact; %Shift-factor
Datashft(Targ,7) = log10(Datashft(Targ,6));

disp(['log10(Shift-factor) = ',num2str(Datashft(Targ,7))]);


for i = [1:Targ-1,Targ+1:MEST]
	Datashft(i,1) = DataMob(i,1)*Datashft(Targ,6);
	Datashft(i,2) = DataMob(i,2);
	Datashft(i,3) = DataMob(i,3);
	Datashft(i,4) = DataMob(i,4);
	Datashft(i,5) = DataMob(i,5);
	Datashft(i,6) = Datashft(Targ,6);
	Datashft(i,7) = log10(Datashft(i,6));
end

end

function [Answ]=DichotSolve(Data,Const,Domain)

%This function solves f(X) = Const w/ regard to X, where f is represented by Data(:,1)
%and Data(:,2). Domain is the first domain which will be cut by Dichotomy.

Func = pchip(Data(:,1),Data(:,3));
Doma = Domain;

Answ = mean(Doma);
Rech = (Doma(2)-Doma(1))/2;
count = 0;

while abs(ppval(Func,Answ)-Const) > 0.001*Const
	if abs(ppval(Func,Answ+Rech)-Const)>abs(ppval(Func,Answ-Rech)-Const)
		Answ = Answ-Rech/2;
	else
		Answ = Answ+Rech/2;
	end
	Rech = Rech/2;
	if abs(ppval(Func,Answ+Rech)-ppval(Func,Answ-Rech)) < 0.0001 * Const
		Rech = Rech * 2^(2+rand);
    end
    count = count+1;
    if count > 300000
        disp('Giving up the dichotomy');
        break
    end
	%disp([num2str(ppval(Func,Answ)-Const), ' ' , num2str(Const)]);
	%disp([num2str(Answ),' ',num2str(Rech)]);
end
end

