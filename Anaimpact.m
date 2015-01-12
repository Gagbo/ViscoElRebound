function [NewFile,Dura,Peak,Int,Cond] = Anaimpact(filename)

%%%%Variable initialisation%%%%
Seuil = 0.60;
Peak = 0;
PeakInd = 0;
Dura = 0;
Base = 0;
BegI = 1;
BegL = 0;
EndI = 1;
StartFound = 0;
EndFound = 0;
Int = 0;

%%%%Measure conditions reading%%%%
Cond = sscanf(filename,'%d %*[-] %d %*[-] %d %*[-] %d');
%Cond(1) is the date/series
%Cond(2) is the temperature
%Cond(3) is the number of the ball used
%Cond(4) is a boolean. 1 if there was the steel rod.

%%%%File reading%%%%
FileID = fopen(filename);
File = textscan(FileID,'%f %f %*[^\n]','Delimiter',',','HeaderLines',2);
%File = textscan(FileID,'%f %f %*[^\n]','Delimiter',' ','MultipleDelimsAsOne',1,'HeaderLines',3,'BufSize',10000000);
fclose(FileID);
File = [File{1},File{2}];
Siz = length(File);
%First column : time in seconds
%Second column : voltage in volts

%%%%Data analysis%%%%
Peak = File(1,2);
Base = File(1,2);
PeakInd = 1;
StartFound = 0;
EndFound = 0;

for i = 1:Siz
	if File(i,2) > Peak
		PeakInd = i;
		Peak = File(i,2);
	end

	if ~StartFound && (File(i,2)-Base)/Base > Seuil
		BegI = i;
		BegL = File(i,2);
		StartFound = 1;
		disp('Start of impact found');
	end

	if StartFound && ~EndFound && File(i,2) < BegL
		if File(i,1) - File(BegI,1) < 0.00015
			continue
		end
	
		EndI = i;
		EndFound = 1;
		disp('End of impact found');
	end

end

if ~StartFound
	disp('Start of impact not found');
elseif ~EndFound
	disp('End of impact not found');
end

Dura = File(EndI,1) - File(BegI,1);
Int = trapz(File(BegI:EndI,1),File(BegI:EndI,2));

%%%%Creation of Impact Data starting at t=0%%%
NewSiz = EndI-BegI+1;
NewFile = File(BegI:EndI,:);
for i = 1:NewSiz
	NewFile(i,1) = NewFile(i,1) - File(BegI,1);
	NewFile(i,2) = NewFile(i,2) - File(EndI,2);
end

disp([filename,' : the peak is at ',num2str(Peak),' Volts']);
disp([filename,' : the duration of impact is ',num2str(Dura),' seconds']);
disp([filename,' : the area under the impact is ',num2str(Int),' Volts.seconds']);

%%%%Verification of the acquisition of the impact%%%%
h = figure;
Curve = plot(File(:,1),File(:,2),'Linewidth',1,'Color','blue');
hold all;
xlim([File(1,1) File(1+NewSiz,1)]);
ylim([0 1.2*Peak]);
%plot([File(BegI,1),File(BegI,1)],[0,Peak],'Linewidth',1,'Color','red');
%plot([File(EndI,1),File(EndI,1)],[0,Peak],'Linewidth',1,'Color','red');
%title([filename,' : temperature is ', num2str(Cond(2)),' C, ball is n ',num2str(Cond(3))]);
xlabel('Time (s)','FontSize',15);
ylabel('Force (N)','FontSize',15);
%print(h,'-depsc','-r300','/path/to/eps');
hold off;

end
