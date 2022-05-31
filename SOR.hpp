#include <bits/stdc++.h>
#define CNT_MAX 200
#define TOL 1e-6
#define omega 1
using namespace std;

double conv(vector<double> a, vector<double> b) {
  double up = 0, down = 0;
  int N = a.size();
  for (int i = 0; i < N; i++) {
    up += (a.at(i) - b.at(i)) * (a.at(i) - b.at(i));
    down += b.at(i) * b.at(i);
  }
  return up / down;
}

vector<double> SOR(vector<vector<double> > A, vector<double> b) {
  int N = b.size();
  vector<double> x(N, nan("")), x_past(N,0);
  for (int k = 0; k < CNT_MAX; k++) {
    for (int i = 0; i < N; i++) {
      double sigma = 0;
      for(int j = 0; j < i; j++) {
	sigma += A.at(i).at(j) * x.at(j);
      }
      for(int j = i+1; j < N; j++) {
	sigma += A.at(i).at(j) * x_past.at(j);
      }
      sigma = (b.at(i) - sigma) / A.at(i).at(i);
      x.at(i) = x_past.at(i) + omega * (sigma - x_past.at(i));
    }
    // cout << x.at(0) << ' ' << x.at(1) << ' ' << x.at(2) << endl;
    swap(x, x_past);
    if (conv(x, x_past) < TOL) break;
  }
  return x_past;
}

