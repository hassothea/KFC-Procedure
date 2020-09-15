
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double euclidDist(NumericVector u, NumericVector v)
{
  double d = 0;
  for (int i = 0; i < u.size(); i++)
  {
    d += (u[i] - v[i]) * (u[i] - v[i]);
  }
  return d;
}

// [[Rcpp::export]]

double gklDist(NumericVector u, NumericVector v)
{
  double d = 0;
  for (int i = 0; i < u.size(); i++)
  {
    d += u[i] * log(abs(u[i] / v[i])) - (u[i] - v[i]);
  }
  return d;
}

// [[Rcpp::export]]
double logisticDist(NumericVector u, NumericVector v)
{
  double d = 0;
  double tot_u = std::accumulate(u.begin(), u.end(), 0.0, std::plus<double>());
  double tot_v = std::accumulate(v.begin(), v.end(), 0.0, std::plus<double>());
  double x_i = 0, y_i = 0;
  for (int i = 0; i < u.size(); i++)
  {
    x_i = u[i] / tot_u;
    y_i = v[i] / tot_v;
    d += x_i * log(abs(x_i / y_i)) + (1 - x_i) * log(abs((1 - x_i) / (1 - y_i)));
  }
  return d;
}

// [[Rcpp::export]]
double itakuraDist(NumericVector u, NumericVector v)
{
  double d = 0;
  for (int i = 0; i < u.size(); i++)
  {
    d += u[i] / v[i] - log(abs(u[i] / v[i])) - 1;
  }
  return d;
}

// [[Rcpp::export]]
double expDist(NumericVector u, NumericVector v)
{
  double d = 0;
  for (int i = 0; i < u.size(); i++)
  {
    d += exp(u[i]) - exp(v[i]) - (u[i] - v[i]) * exp(v[i]);
  }
  return d;
}

// [[Rcpp::export]]
double polyDist(NumericVector u, NumericVector v, int deg)
{
  double d = 0;
  if (deg % 2 == 0)
  {
    for (int i = 0; i < u.size(); i++)
    {
      d += std::pow(u[i], deg) - std::pow(v[i], deg) - deg * (u[i] - v[i]) * std::pow(v[i], deg - 1);
    }
  }
  else
  {
    for (int i = 0; i < u.size(); i++)
    {
      if (v[i] > 0)
      {
        d += std::pow(u[i], deg) - std::pow(v[i], deg) - deg * (u[i] - v[i]) * std::pow(v[i], deg - 1);
      }
      else
      {
        d += std::pow(u[i], deg) - std::pow(v[i], deg) + deg * (u[i] - v[i]) * std::pow(v[i], deg - 1);
      }
    }
  }
  return d;
}

