%% Loading Data
data = readtable("GDP.csv");
all_gdp = data.GDP;
all_date = data.DATE;

%% Historical Mean Calculation
hm = zeros(303,1);
hm = all_gdp;
for i=1:46 
    hm(i+257) = mean(all_gdp(257:257+i));
end
MSE_hm = mean((hm(258:257+newlength)-all_gdp(258:257+newlength)).^2);

%% Testing Stationarity
% (The following code is commented out, but can be uncommented if needed)
% Y = log(gdp);
% h = adftest(Y,Model="TS",Lags=0:2);

%% Values from 2001 to Present
gdp = data.GDP(218:257);
date = data.DATE(218:257);

%% Plotting Initial Data
figure;
plot(date,gdp, 'r', 'LineWidth',2);
xlabel('Years');
ylabel('GDP in Billions');

%% ACF and PACF Plots
figure;
subplot(2,1,1);
autocorr(log(gdp));
title("ACF");
subplot(2,1,2);
parcorr(log(gdp));
title("PACF");

%% AIC and BIC Computation
[aic,bic]= aicbic(LogLikelihood,2,40);
Model_table = "arima(0,1,0)";
aic_table = aic;
bic_table = bic;

%% Computation of AIC and BIC tests
m=2;
for i = 1:2
    for j = 1:2
        for k=1:2
            arima_model(m-1) = arima(i,j,k);
            [~,~,LoglikehoodE] = estimate(arima_model(m-1),gdp,'display','off');
            [aic_inter,bic_inter] = aicbic(LoglikehoodE,2,250);
            Model_table(m) = ['arima(',num2str(i),',',num2str(j),',',num2str(k),')'];
            aic_table(m) = aic_inter;
            bic_table(m) = bic_inter;
            m=m+1;
        end
    end
end

Comparision = table(Model_table',aic_table',bic_table');
disp(Comparision);

%% AIC Minimum Value
AICMin = aic_table(1);
AICindex = 1;
for i=1:length(aic_table)
    if(aic_table(i)<AICMin)
        AICindex=i;
        AICMin=aic_table(i);
    end
end
disp(aic_table);
disp("Lowest AIC value at:");
disp(AICindex);
disp(Model_table(AICindex));

%% BIC Minimum Value
BICMin = bic_table(1);
BICindex = 1;
for i=1:length(bic_table)
    if(bic_table(i)<BICMin)
        BICindex=i;
        BICMin=bic_table(i);
    end
end
disp(bic_table);
disp("Lowest BIC value at:");
disp(BICindex);
disp(Model_table(BICindex));

%% ARIMA Model Estimations
ARIMA_Theo = arima(0,1,0); %Theoretical ARIMA model
[ARIMA_Theo1,~,LogLikelihood] = estimate(ARIMA_Theo, gdp);

ARIMA_Prac = arima(2,1,2); %Practical ARIMA model obtained after the AIC BIC tests 
[ARIMA_Prac1,~,LogLikelihood2] = estimate(ARIMA_Prac, gdp);

%% RMS Values
residual = infer(ARIMA_Theo1, gdp); 
root_mean_square_error_Theo1 = sqrt(sum(residual.^2))/length(residual);

residual2 = infer(ARIMA_Prac1, gdp );
root_mean_square_error_Prac1 = sqrt(sum(residual2.^2))/length(residual2);

prediction = gdp +residual; 
prediction2 = gdp + residual2;

%% Forecasting the Future
% In this section, we will be forecasting the future data for three months and
% compare it with the original values to get a better understanding.
newlength = 40;

%% ARIMA(0,1,0) Forecasting
[Final_Theo_Forecast,ymse1] = forecast(ARIMA_Theo1, newlength,'Y0', gdp);
lower_theo = Final_Theo_Forecast - 1.96*sqrt(ymse1);
upper_theo = Final_Theo_Forecast + 1.96*sqrt(ymse1);

figure;
h1 = plot(date,gdp,'Color','g', 'LineWidth',1.5);
hold on;
h2 = plot(all_date(258:257+newlength),Final_Theo_Forecast,'LineWidth',2, 'Color','r');
h3 = plot(all_date(258:257+newlength),all_gdp(258:257+newlength));
h4 = plot(all_date(258:257+newlength),lower_theo,'r:','LineWidth',1);
h5 = plot(all_date(258:257+newlength),upper_theo,'r:','LineWidth',1);
legend([h1,h2,h3,h4],"Original Data ","Theoretical Arima - ARIMA(0,1,0)","True ","95% Confidence Interval","Location","Best");
grid on;
hold off;

%% ARIMA(2,1,2) Forecasting
[Final_Prac_Forecast,ymse1] = forecast(ARIMA_Prac1, newlength,'Y0', gdp);
lower_prac = Final_Prac_Forecast - 1.96*sqrt(ymse1);
upper_prac = Final_Prac_Forecast + 1.96*sqrt(ymse1);

figure;
h1 = plot(date,gdp,'Color','g', 'LineWidth',1.5);
hold on;
h2 = plot(all_date(258:257+newlength),Final_Prac_Forecast,'LineWidth',2, 'Color','r');
h3 = plot(all_date(258:257+newlength),all_gdp(258:257+newlength));
h4 = plot(all_date(258:257+newlength),lower_prac,'r:','LineWidth',2);
h5 = plot(all_date(258:257+newlength),upper_prac,'r:','LineWidth',2);
h6 = plot(all_date(258:257+newlength),hm(258:257+newlength),'LineWidth',2, 'Color','m', 'LineStyle','--');
legend([h1,h2,h3,h4, h6],"Original Data","Practical Arima - ARIMA(2,1,2)","Original 2018","95% Confidence Interval", "Historical Mean","Location","Best");
grid on;
hold off;

%% Historical Mean Plot
figure;
hold on;
k1 = plot(all_date,all_gdp,'LineWidth',2, 'Color','r');
k2 = plot(all_date(258:257+newlength),hm(258:257+newlength),'LineWidth',2, 'Color','m', 'LineStyle','--');
legend([k1,k2],"Original Data ","Historical Mean", "Location","Best");
grid on;
hold off;

%% Additional Forecasting for 10 years
sys_2 = arima(2,0,0);
Md1_AR = estimate(sys_2,all_gdp);
yf = forecast(Md1_AR,40,'Y0',all_gdp); % Forecast 200 points ahead
plot(1:length(all_gdp),all_gdp,'b',length(all_gdp):length(all_gdp)+length(yf),[all_gdp(end);yf],'r');
hold on;

%% ARIMA Forecasts
sys_3 = arima(2,1,2);
Md1_ARIMA = estimate(sys_3,all_gdp);
yf = forecast(Md1_ARIMA,40,'Y0',all_gdp);
plot(length(all_gdp):length(all_gdp)+length(yf),[all_gdp(end);yf],'k');

sys_3 = arima(0,1,0);
Md1_ARIMA = estimate(sys_3,all_gdp);
yf = forecast(Md1_ARIMA,40,'Y0',all_gdp);
plot(length(all_gdp):length(all_gdp)+length(yf),[all_gdp(end);yf],'c');
hold off;
