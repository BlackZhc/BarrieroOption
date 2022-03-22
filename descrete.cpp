#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>

using namespace std;

float risk_free_rate, initial_stock_price, expiration_time, volatility;
float strike_price, barrier_price;
int no_of_trials, no_of_divisions;
vector<double> probability;
vector<double> price;
vector<int> act;

unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator;

double get_uniform()
{
    // http://www.cplusplus.com/reference/random/exponential_distribution/
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double number = distribution(generator);
    return (number);
}

double max(double a, double b) {
    return(b < a) ? a : b;
}

double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a = fabs(z);
    double t = 1.0 / (1.0 + a * p);
    double b = c2 * exp((-z) * (z / 2.0));
    double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
    n = 1.0 - b * n;
    if (z < 0.0) n = 1.0 - n;
    return n;
};

double option_price_put_black_scholes(const double& S,      // spot price
    const double& K,      // Strike (exercise) price,
    const double& r,      // interest rate
    const double& sigma,  // volatility
    const double& time)
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
    double d2 = d1 - (sigma * time_sqrt);
    return K * exp(-r * time) * N(-d2) - S * N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
    const double& K,       // strike (exercise) price,
    const double& r,       // interest rate
    const double& sigma,   // volatility
    const double& time)      // time to maturity
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
    double d2 = d1 - (sigma * time_sqrt);
    return S * N(d1) - K * exp(-r * time) * N(d2);
};


float closed_form_down_and_out_european_call_option()
{
    // I took this formula from Wilmott, Howison and Dewynne, "The Mathematics of Financial Derivatives"
    float K = (2 * risk_free_rate) / (volatility * volatility);
    float A = option_price_call_black_scholes(initial_stock_price, strike_price,
        risk_free_rate, volatility, expiration_time);
    float B = (barrier_price * barrier_price) / initial_stock_price;
    float C = pow(initial_stock_price / barrier_price, -(K - 1));
    float D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
    return (A - D * C);
}

float closed_form_down_and_in_european_put_option()
{
    // just making it easier by renaming the global variables locally
    float S = initial_stock_price;
    float r = risk_free_rate;
    float T = expiration_time;
    float sigma = volatility;
    float H = barrier_price;
    float X = strike_price;

    // Took these formulae from some online reference
    float lambda = (r + ((sigma * sigma) / 2)) / (sigma * sigma);
    float temp = 2 * lambda - 2.0;
    float x1 = (log(S / H) / (sigma * sqrt(T))) + (lambda * sigma * sqrt(T));
    float y = (log(H * H / (S * X)) / (sigma * sqrt(T))) + (lambda * sigma * sqrt(T));
    float y1 = (log(H / S) / (sigma * sqrt(T))) + (lambda * sigma * sqrt(T));
    return (-S * N(-x1) + X * exp(-r * T) * N(-x1 + sigma * sqrt(T)) +
        S * pow(H / S, 2 * lambda) * (N(y) - N(y1)) -
        X * exp(-r * T) * pow(H / S, temp) * (N(y - sigma * sqrt(T)) - N(y1 - sigma * sqrt(T))));
}

float closed_form_down_and_out_european_put_option()
{
    float vanilla_put = option_price_put_black_scholes(initial_stock_price, strike_price,
        risk_free_rate, volatility, expiration_time);
    float put_down_in = closed_form_down_and_in_european_put_option();
    return (vanilla_put - put_down_in);
}




int main(int argc, const char * argv[]) {
   
/*
   expiration_time = 1;
   risk_free_rate = 0.05;
   volatility = 0.25;
   initial_stock_price = 50;
   no_of_divisions  = 1000;
   no_of_trials = 10000;
   strike_price = 40;
   barrier_price = 20;
    */

    
    sscanf(argv[1], "%f", &expiration_time);
    sscanf(argv[2], "%f", &risk_free_rate);
    sscanf(argv[3], "%f", &volatility);
    sscanf(argv[4], "%f", &initial_stock_price);
    sscanf(argv[5], "%f", &strike_price);
    sscanf(argv[6], "%d", &no_of_trials);
    sscanf(argv[7], "%d", &no_of_divisions);
    sscanf(argv[8], "%f", &barrier_price);

   
   double delta_T = expiration_time/((double) no_of_divisions);
   double delta_R = (risk_free_rate - 0.5*pow(volatility,2))*delta_T;
   double delta_SD = volatility*sqrt(delta_T);
    
    
    for (int i = 0; i < no_of_trials; i++){
        double current_stock_price1 = initial_stock_price;
        double current_stock_price2 = initial_stock_price;
        double current_stock_price3 = initial_stock_price;
        double current_stock_price4 = initial_stock_price;
        int act1 = 0; int act2 = 0; int act3 = 0; int act4 = 0;
        
        for (int j = 0; j < no_of_divisions; j++)
        {

            double x = get_uniform();
            double y = get_uniform();
            double a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
            double b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
            
            
            current_stock_price1 = current_stock_price1*exp(delta_R + delta_SD*a);
            if (current_stock_price1 <= barrier_price)
                act1 = 1;
            
            current_stock_price2 = current_stock_price2*exp(delta_R - delta_SD*a);
            if (current_stock_price2 <= barrier_price)
                act2 = 1;
            
            current_stock_price3 = current_stock_price3*exp(delta_R + delta_SD*b);
            if (current_stock_price3 <= barrier_price)
                act3 = 1;
            
            current_stock_price4 = current_stock_price4*exp(delta_R - delta_SD*b);
            if (current_stock_price4 <= barrier_price)
                act4 = 1;
            
        }
        probability.push_back(exp( (-2.0/(volatility*volatility*expiration_time)) *
                                  log(initial_stock_price/barrier_price)*
                                  log(current_stock_price1/barrier_price) ));
        probability.push_back(exp( (-2.0/(volatility*volatility*expiration_time)) *
                                  log(initial_stock_price/barrier_price)*
                                  log(current_stock_price2/barrier_price) ));
        probability.push_back(exp( (-2.0/(volatility*volatility*expiration_time)) *
                                  log(initial_stock_price/barrier_price)*
                                  log(current_stock_price3/barrier_price) ));
        probability.push_back(exp( (-2.0/(volatility*volatility*expiration_time)) *
                                  log(initial_stock_price/barrier_price)*
                                  log(current_stock_price4/barrier_price) ));

        price.push_back(current_stock_price1);
        price.push_back(current_stock_price2);
        price.push_back(current_stock_price3);
        price.push_back(current_stock_price4);
        act.push_back(act1);
        act.push_back(act2);
        act.push_back(act3);
        act.push_back(act4);
        
    }

    
        double call = 0;
        double count = 0;
        for (int i = 0; i < no_of_trials; i++){
            if (act[i] == 0){
                call = call + max(0.0, price[i] - strike_price);
                count = count + 1;}
        }
        call = (call / count)* exp(-risk_free_rate * expiration_time);
        
        double call_adjust = 0;
        for (int i = 0; i < no_of_trials; i++){
            call_adjust = call_adjust + max(0.0, price[i] - strike_price)*(1 - probability[i]);
        }
        call_adjust = (call_adjust / no_of_trials)* exp(-risk_free_rate * expiration_time);
    
    
    
        double put = 0;
        for (int i = 0; i < no_of_trials; i++){
            if (act[i] == 0){
                put = put + max(0.0, strike_price- price[i]);
                }
        }
        put = (put / count)* exp(-risk_free_rate * expiration_time);
        
        double put_adjust = 0;
        for (int i = 0; i < no_of_trials; i++){
            put_adjust = put_adjust + max(0.0, strike_price - price[i])*(1 - probability[i]);
        }
        put_adjust = (put_adjust / no_of_trials)* exp(-risk_free_rate * expiration_time);
        
    

    
    
    cout << "--------------------------------" << endl;
    cout << "Down-and-out Continuous Barrier  Options Pricing via Monte Carlo Simulation" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price= " << strike_price << endl;
    cout << "Barrier Price= " << barrier_price << endl;
    cout << "Number of Trials= " << no_of_trials << endl;
    cout << "Number of Divisions= " << no_of_divisions << endl;
    cout << "--------------------------------" << endl;
    cout << "--------------------------------" << endl;
    cout << "The average Call Price by explixit simulation= " << call << endl;
    cout << "The Call Price using the (1-p)-adjustment term= " << call_adjust << endl;
    cout << "Theoretical Call Price= " << closed_form_down_and_out_european_call_option() << endl;
    cout << "--------------------------------" << endl;
    cout << "--------------------------------" << endl;
    //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    cout << "The average Put Price by explixit simulation= " << put << endl;
    cout << "The Put Price using the (1-p)-adjustment term= " << put_adjust << endl;
    cout << "Theoretical Put Price= " << closed_form_down_and_out_european_put_option() << endl;
    return 0;
    
    }










