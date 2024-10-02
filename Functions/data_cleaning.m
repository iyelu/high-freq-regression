function output = data_cleaning(dir_rawData, dir_output)
% -------------------------------------------------------------------------
% Purpose: Process the raw data, saved in 'dir_rawData', for Models I-IV
%          in the empirical illustrations in the paper
%          "Understanding regressions with observations collected at
%           high-frequency over long span" by Chang, Lu and Park (2024)
% -------------------------------------------------------------------------
% Outputs:
%       'dir_output/M1_TB.csv': 3-month T-bill rate and 10-year T-bond rate
%       'dir_output/M2_TBEU.csv': 3-month T-bill rate and 3-month Eurodollar rate
%       'dir_output/M3_logFX.csv': Log US/UK spot and 3-month forward exchange rates
%       'dir_output/M4_logSP500.csv': Log S&P 500 index future and S&P 500 index
% -------------------------------------------------------------------------
% This version: September 2024 by Ye Lu (ye.lu1@sydney.edu.au)
% -------------------------------------------------------------------------

%% Model I: 3-month T-bill rate and 10-year T-bond rate (1962-Jan-02 to 2019-Dec-31)
% --------
% Raw data downloaded from FRED dataset provided by Federal Reserve Bank
% of St. Louis, saved in
%   'DTB3.csv':  3-month treasury bill market rate, percent, daily
%   'DGS10.csv': 10-year treasury constant maturity rate, percent, daily
% --------
% 3-month T-bill rate
% --------
tb3  = readtable(strcat(dir_rawData,'DTB3.csv'));
data = tb3.DTB3;
% Replace missing values with the data on the previous day
ind = find(isnan(data)); % indice of missing values
while ~isempty(find(isnan(data),1))
    data(ind) = data(ind-1);
    ind = find(isnan(data)); % find again missing values
end
tb3 = timetable(tb3.DATE,data);
% --------
% 10-year T-bond rate
% --------
tb10 = readtable(strcat(dir_rawData,'DGS10.csv'));
data = tb10.DGS10;
% Replace missing values with the data on the previous day
ind = find(isnan(data)); % indice of missing values
while ~isempty(find(isnan(data),1))
    data(ind) = data(ind-1);
    ind = find(isnan(data)); % find again missing values
end
tb10 = timetable(tb10.DATE,data);
% --------
% Save 3-month T-bill and 10-year T-bond rates
% --------
T1 = timetable(tb3.Time,tb3.Var1,tb10.Var1);
T1.Properties.VariableNames = {'TB3m','TB10y'};
writetimetable(T1,strcat(dir_output,'M1_TB.csv'))
output.M1 = T1;

%% Model II: 3-month T-bill rate and 3-month Eurodollar rate (1971-Jan-04 to 2016-Oct-07)
% --------
% Raw data downloaded from FRED dataset provided by Federal Reserve Bank
% of St. Louis, saved in
%   'DTB3_FRED.csv': 3-month treasury bill market rate, percent, daily
%   'DED3_FRED.csv': 3-month Eurodollar deposit rate, percent, daily
% --------
% 3-month Eurodollar ate
% --------
eudol = readtable(strcat(dir_rawData,'DED3.csv'));
data  = eudol.DED3;
% Replace missing values with the data on the previous day
ind = find(isnan(data)); % indice of missing values
while ~isempty(find(isnan(data),1))
    data(ind) = data(ind-1);
    ind = find(isnan(data)); % find again missing values
end
eudol = timetable(eudol.DATE,data);
% --------
% 3-month T-bill rate
% --------
S = timerange('1971-Jan-04','2016-Oct-07', 'closed');
tb3_eudol = tb3(S,:);
% --------
% Save 3-month T-bill and 3-month Eurodollar rates
% --------
T2 = timetable(tb3_eudol.Time, tb3_eudol.Var1, eudol.Var1);
T2.Properties.VariableNames = {'TB3m','Eudollar'};
writetimetable(T2,strcat(dir_output,'M2_TBEU.csv'))
output.M2 = T2;

%% Model III: Log US/UK spot and 3-month forward exchange rates (1979-Jan-02 to 2017-Dec-29)
% --------
% Raw data downloaded from Bank of England Database, saved in
%   'USUKS.xls': US/UK spot exchange rates, daily
%   'USUKF.xls': US/UK forward exchange rates, daily
% --------
% Spot exchange rate
% --------
[~,~,data] = xlsread(strcat(dir_rawData,'USUKS.xls'));
spot = cell2mat(data(3:end,2)); % US/UK spot exchange rates
% Construct datetime variable for the dates
date_spot = cell2mat(data(3:end,1));
d = str2num(date_spot(:,1:2));   % day
m = string(date_spot(:,4:6));    % month
y = str2num(date_spot(:,8:9)); y = [y(y<70)+2000;y(y>70)+1900]; % year
s = cell(9858,1);
for i = 1:9858
    s{i}=sprintf('%d-%s-%d',y(i),m(i),d(i));
end
Time = datetime(s);
% --------
% Forward exchange rate
% --------
[~,~,data] = xlsread(strcat(dir_rawData, 'USUKF.xls'));
frwd3m  = cell2mat(data(3:end,3)); % US/UK 3-month forward exchange rates
% --------
% Save log US/UK spot and 3-month forward exchange rates
% --------
T3 = sortrows(timetable(Time, log(spot), log(frwd3m)));
T3.Properties.VariableNames = {'Spot','Frwd3m'};
writetimetable(T3,strcat(dir_output,'M3_logFX.csv'))
output.M3 = T3;

%% Model IV: Log S&P 500 index future and S&P 500 index (1997-Sep-10 to 2018-Mar-29)
% --------
% Raw S&P 500 index price data downloaded at Yahoo Finance, saved in
%   'SP500.csv': S&P 500 index price, daily
% Raw S&P 500 index futures data downloaded at Investing.com, saved in
%   'SP500F.csv': S&P 500 index futures price, daily
% --------
% S&P 500 daily adjusted close price
% --------
dataC = readtable(strcat(dir_rawData, 'SP500.csv'));
% remove all the variables but the adjusted close price
dataC = removevars(dataC, {'Open', 'High', 'Low', 'Close', 'Volume'});
dataC.Properties.VariableNames = {'Date','data'};
% --------
% S&P 500 futures daily close price
% --------
dataF = readtable(strcat(dir_rawData, 'SP500F.csv'));
% remove all the variables but the close price (first non-date variable)
dataF = removevars(dataF, {'Var3', 'Var4', 'Var5', 'Var6', 'Var7'});
dataF.Properties.VariableNames = {'Date','data'};
% convert text to numeric values
dataF.data = str2double(dataF.data);
% Match the span of SP500 and SP500 futures
dataF = dataF(ismember(dataF.Date, dataC.Date),:);
% check the dates in two data sets match (answer 0: good to go)
sum(dataF.Date ~= dataC.Date);
% --------
% Save log S&P 500 index price and S&P index futures price
% --------
T4 = timetable(dataC.Date, log(dataC.data), log(dataF.data));
T4.Properties.VariableNames = {'SP500','SP500futures'};
writetimetable(T4,strcat(dir_output,'M4_logSP500.csv'))
output.M4 = T4;

end
