clc, clear all
clf;
name = 'fig1_data';
figname = sprintf('x%s',name);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf,'DefaultAxesLineWidth',2);
set(gcf,'PaperPositionMode','auto');


for k = 1:2
    if k == 1
        cd DataSet1
    else
        cd DataSet2
    end
sumdata = zeros(13,3);
fid = fopen('accnums.txt');
cells = textscan(fid,'%s');
fclose(fid);
numaccs = length(cells{1});

for i = 1:numaccs
  accnum = cell2mat(cells{1}(i));	 
  filename = ['output/output_' accnum '.txt'];
  fid = fopen(filename,'r');
  tmp = fscanf(fid,'%d');
  fclose(fid);
  data(:,:,i) = reshape(tmp,13,3);
  sumdata = sumdata + data(:,:,i);
end

% make frequency data rather than counts
for j =1:3
 pdata(j,:) = sumdata(:, j)./sum(sumdata(:, j));
end

% keep only those proteins >1%

totals = sum(sumdata');
all = sum(totals);
ikeep = find(totals./all>.01);  

% order the protein types by their % change in frequency
[deltas,order] = sort((pdata(1,ikeep)-pdata(3,ikeep))./pdata(1,ikeep),'descend');

labels = [ 'Term'; 'Port'; 'Head'; 'Injn'; 'Tail'; 'Prot'; 'Tran'; 'Intg'; 'Lyss'; 'Plat'; 'Caps'; 'Lysn'; 'Flip'];
titles = [ 'intact'; 'questn'; 'incomp'];

pc = 100*(pdata(3,ikeep(order))-pdata(1,ikeep(order)))./pdata(1,ikeep(order));

% %statistical significance

 totals = sum(sumdata');
 iqi = sum(sumdata);   
 fiqi = iqi/sum(iqi);
 
 nreps = 10000;
for i = 1:nreps
  for g = 1:length(sumdata)
    r = rand(totals(g),1);
    fdata(1,g) = length(find(r<fiqi(1)));  
    fdata(3,g) = length(find(r>1-fiqi(3))); 
  end
  fpdata(1,:) = fdata(1,:)./sum(fdata(1,:));
  fpdata(3,:) = fdata(3,:)./sum(fdata(3,:));
  pdiffs(:,i) = 100*(fpdata(3,:) - fpdata(1,:))./fpdata(3,:);
end
  
for g=1:length(sumdata)
  tmp = sort(pdiffs(g,:));
  p25(g) = tmp(floor(0.025*nreps));
  p975(g) = tmp(floor(0.975*nreps));
end

%  with p25 and p975 which are computed for all 13 proteins and not shuffled
isig = find(pc<p25(ikeep(order)));
jsig = find(pc>p975(ikeep(order)));


fprintf(1,'Significantly lower:\n');
labels(ikeep(order(isig)),:)
fprintf(1,'Significantly higher:\n');
labels(ikeep(order(jsig)),:)


subplot(2,2,2+k)
grey = 0.5*[1 1 1];
b = bar(pc, 'FaceColor',grey,'EdgeColor',grey,'LineWidth',1);
hold on;

plot(isig, pc(isig)-10, '*r', 'MarkerSize', 10)
hold on
plot(jsig, pc(jsig)+10, '*g', 'MarkerSize', 10)
ylabel({'Percent change'},'fontsize',18,'verticalalignment','bottom','interpreter','latex');
set(gca,'xticklabel',labels(ikeep(order),:));
%set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',16);
if k== 1
    annotation('textbox', [0.05, 0.55, 0.0, 0.0], 'String', "C", 'fontweight','bold', 'FontSize', 18)
else
    annotation('textbox', [0.5, 0.5, 0.0, 0.0], 'String', "D", 'fontweight','bold', 'FontSize', 18)
end
xlim([0.5, 11.5]);   

finaldata = [pdata(1, ikeep(order)); pdata(2, ikeep(order)); pdata(3, ikeep(order))];
subplot(2,2,k)
H = bar(finaldata');
H(1).FaceColor = [0.2 0.2 0.2];
H(2).FaceColor = [0.5 0.5 0.5];
H(3).FaceColor = [0.8 0.8 0.8];
if k ==1
    legend({'Intact', 'Questionable', 'Incomplete'}, 'Location', 'Best', 'fontsize',16,'interpreter','latex');
    legend('boxoff')
end
set(gca,'fontsize',16);
ylabel({'Frequency'}, 'fontsize',18,'verticalalignment','bottom','interpreter','latex');
set(gca,'xticklabel',labels(ikeep(order),:));
%set(gca,'TickLabelInterpreter','latex')
if k==1
    annotation('textbox', [0.05, 0.999999, 0.0, 0.0], 'String', "A", 'fontweight','bold', 'FontSize', 18)
else
    annotation('textbox', [0.5, 0.999999, 0.0, 0.0], 'String', "B", 'fontweight','bold', 'FontSize', 18)
end
cd ..
clear all
end

stampname(13.7268748547337, -0.386146607954113 ,90);

creat(name);

