function plot_state_probabilities(prb,time,data,max_ind,S)

ind=max_ind;

% figure
s1=subplot(5,1,4);
yyaxis left
p1=plot(time,data,'.');
ylabel('Data','fontsize',14)
yyaxis right
p2=plot(time,ind);
ylim([0.9 3.1])
ylabel('State Label','fontsize',14)
datetick
axis tight

s2=subplot(5,1,5);
h=area(time,prb','LineStyle','none');
% h(1).FaceColor = [0 0 1];
% h(2).FaceColor = [1 0.4 0.4];
% h(3).FaceColor = [1 0 0];
datetick
axis tight
xlabel('Date/Time','fontsize',14)
ylabel('State Prb','fontsize',14)

linkaxes([s1 s2], 'x')