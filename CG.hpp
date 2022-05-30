#include <bits/stdc++.h>
#define tol 1e-2
#define CNT_MAX 200
using namespace std;

double norm2(vector<double> &x, vector<double> &y) {
  double ans = 0;
  for (int j=0; j<x.size(); j++) {
    ans += x.at(j) * y.at(j);
  }
  return ans;
}

// http://www.jicfus.jp/wiki/index.php?Bi-CGSTAB%20%E6%B3%95
void  Bi_CGSTAB(vector<double> A_val, vector<int> A_r, vector<int> A_c, vector<double> b, vector<double> &x) {
  int N = b.size();
  vector<double> r(N), p(N), _r(N), r_past(N);
  for (int j=0; j<N; j++) {
    for (int k=A_r.at(j); k<A_r.at(j+1); k++) {
      b.at(j) -= A_val.at(k) * x.at(A_c.at(k));
    }
    r.at(j) = b.at(j);
    p.at(j) = b.at(j);
    _r.at(j) = b.at(j);
  }
  int cnt = 0;
  
  while (cnt < CNT_MAX && tol < norm2(r,r)) {
    
    cnt++;

    // Line 6
    vector<double> ap(N,0);
    for (int j=0; j<N; j++) {
      for (int k=A_r.at(j); k<A_r.at(j+1); k++) {
	ap.at(j) += A_val.at(k) * p.at(A_c.at(k));
      }
    }
    double alpha = norm2(_r, r) / norm2(_r, ap);

    // Line 7
    vector<double> s(N);
    for (int j=0; j<N; j++) {
      s.at(j) = r.at(j) - alpha * ap.at(j);
    }

    // Line 8
    vector<double> as(N,0);
    for (int j=0; j<N; j++) {
      for (int k=A_r.at(j); k<A_r.at(j+1); k++) {
	as.at(j) += A_val.at(k) * s.at(A_c.at(k));
      }
    }
    double omega = norm2(as, s) / norm2(as, as);

    // Line 9
    for (int j=0; j<N; j++) {
      x.at(j) += omega * s.at(j) + alpha * p.at(j);
    }

    // Line 10
    for (int j=0; j<N; j++) {
      r_past.at(j) = r.at(j);
      r.at(j) = s.at(j) - omega * as.at(j);
    }

    // Line 11
    double beta = alpha / omega * norm2(_r, r) / norm2(_r, r_past);

    // Line 12
    for (int j=0; j<N; j++) {
      p.at(j) = r.at(j) + beta * (p.at(j) - omega * ap.at(j));
    }
    cout << norm2(r,r) << endl;
  }
  cout << cnt << ' ' << norm2(r, r) << endl;
}
