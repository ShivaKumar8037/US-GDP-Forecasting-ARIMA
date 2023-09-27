# Forecasting GDP of USA using ARIMA 📈

## Introduction
Gross domestic product (GDP) is the total monetary or market value of all the finished goods and services produced within a country’s borders in a specific time period. As a broad measure of overall domestic production, it functions as a comprehensive scorecard of a given country’s economic health.

## Personal Note
This was my first time diving into the realm of Financial Econometrics. I had the opportunity to explore MATLAB and work with ARIMA models, both of which were new territories for me. The learning curve was steep, but the journey has been incredibly rewarding.

## Objective (Research Question)
The project aims to forecast the GDP of the United States of America with a dataset taken from the Federal Reserve Economic Data.

## Theory
ARIMA, or autoregressive integrated moving average, is a statistical model used for time series forecasting. It is a popular and widely used method for forecasting a variety of economic and financial data, such as GDP, inflation, and stock prices. Detailed theory about the ARIMA model is available in the project documentation.

## Data
FRED provided the data in `.csv` form. The `GDP.csv` contains 2 columns: `DATE` and `GDP`. The frequency of the GDP is quarterly, i.e., GDP of every quarter from the beginning of January 1947 to 2nd Quarter of 2022.

## Research Method
This project considers data from the 1st Quarter of 2001 to the last Quarter of 2010 and uses this decade of data to forecast the data from 2011 to 2022. This forecasted series is then compared with the true values. The detailed research method, including plots and code snippets, is available in the project documentation.

## Results
The project utilized both Theoretical and Practical ARIMA models for forecasting. A comparative analysis based on Root Mean Square Error (RMSE) showed the effectiveness of the models.

## Conclusion
While the ARIMA model provided good forecasts for the years 2011 to 2022, it's not recommended to solely rely on this model for forecasting GDP indefinitely. The computation of macroeconomic factors, such as GDP, is challenging, and researchers might use other advanced models for better accuracy.


