/*

(0,0) (1,0) (2,0) (3,0) (4,0)... (101,0)
(0,1) (1,1) (2,1) (3,1) (4,1)
(0,2) (1,2) (2,2) (3,2) (4,2)
(0,3) (1,3) (2,3) (3,3) (4,3)
(0,4) (1,4) (2,4) (3,4) (4,4)
.                            
.
.
(0,101)

左の壁の温度が15度、右の壁の温度が0度になるように
(0,i)と(4,i) (i=0,1,2,3,...,101) のマスの温度を定める。
左右上下の壁に沿った速度が0になるように
(0,i)と(4,i) (i=0,1,2,3,...,101) のマスの速度を定める。

最初液体は一様に0度とする。

圧力は (0,0) において常に 0 とする。
レイノルズ数:0.71
プラントル数:7100
dt:1e-4
dx=1/100
dy=1/100

 */
#include <bits/stdc++.h>
#include "SOR.hpp"
using namespace std;
#define N 10
// N に合わせる
#define dx 0.1
#define dy 0.1
#define Ra 7100
#define Pr 0.71
#define dt 1e-4
// N に合わせる
#define place(i,j) (i-1)*N+j-2

vector<vector<double> > DIV(vector<vector<double> > &u, vector<vector<double> > &v) {
  vector<vector<double> > ans(N+2, vector<double>(N+2, 0));
  for(int i=N+1; i>=1; i--) {
    for(int j=N+1; j>=1; j--) {
      ans.at(i).at(j) = (u.at(i).at(j)-u.at(i-1).at(j))/dx + (v.at(i).at(j)-v.at(i).at(j-1))/dy;
    }
  }
  return ans;
}

vector<vector<double> > CNVU(vector<vector<double> > &u, vector<vector<double> > &v) {
  vector<vector<double> > ans(N+2, vector<double>(N+2, 0));
  for(int i=N; i>=1; i--) {
    for(int j=N; j>=1; j--) {
      double vave =
	(v.at(i).at(j) +
	 v.at(i).at(j-1) +
	 v.at(i+1).at(j) +
	 v.at(i+1).at(j-1)) / 4;
      ans.at(i).at(j) =
	u.at(i).at(j) * (u.at(i+1).at(j) - u.at(i-1).at(j)) / (2 * dx) -
	abs(u.at(i).at(j)) * (u.at(i-1).at(j) - 2 * u.at(i).at(j) + u.at(i+1).at(j)) / (2 * dx) +
	vave * (u.at(i).at(j+1) - u.at(i).at(j-1)) * (2 * dy) -
	abs(vave) * (u.at(i).at(j-1) - 2 * u.at(i).at(j) + u.at(i).at(j+1)) / (2 * dy);
    }
  }
  return ans;
}
	
vector<vector<double> > DIFU(vector<vector<double> > &u) {
  vector<vector<double> > ans(N+2, vector<double>(N+2, 0));
  for(int i=N; i>=1; i--) {
    for(int j=N; j>=1; j--) {
      ans.at(i).at(j) =
	Pr * ((u.at(i-1).at(j) - 2 * u.at(i).at(j) + u.at(i+1).at(j)) / (dx * dx)
	      + ((u.at(i).at(j-1) - 2 * u.at(i).at(j) + u.at(i).at(j+1)) / (dy * dy)));
    }
  }
  return ans;
}

vector<vector<double> > CNVV(vector<vector<double> > &u, vector<vector<double> > &v) {
  vector<vector<double> > ans(N+2, vector<double>(N+2, 0));
  for(int i=N; i>=1; i--) {
    for(int j=N; j>=1; j--) {
      double uave =
	(u.at(i).at(j) +
	 u.at(i).at(j-1) +
	 u.at(i+1).at(j) +
	 u.at(i+1).at(j-1)) / 4;
      ans.at(i).at(j) =
	u.at(i).at(j) * (v.at(i+1).at(j) - v.at(i-1).at(j)) / (2 * dx) -
	abs(uave) * (v.at(i-1).at(j) - 2 * v.at(i).at(j) + v.at(i+1).at(j)) / (2 * dx) +
	v.at(i).at(j) * (v.at(i).at(j+1) - v.at(i).at(j-1)) * (2 * dy) -
	abs(v.at(i).at(j)) * (v.at(i).at(j-1) - 2 * v.at(i).at(j) + v.at(i).at(j+1)) / (2 * dy);
    }
  }
  return ans;
}
	
vector<vector<double> > DIFV(vector<vector<double> > &v) {
  vector<vector<double> > ans(N+2, vector<double>(N+2, 0));
  for(int i=N; i>=1; i--) {
    for(int j=N; j>=1; j--) {
      ans.at(i).at(j) =
	Pr * ((v.at(i-1).at(j) - 2 * v.at(i).at(j) + v.at(i+1).at(j)) / (dx * dx)
	      + ((v.at(i).at(j-1) - 2 * v.at(i).at(j) + v.at(i).at(j+1)) / (dy * dy)));
    }
  }
  return ans;
}

vector<vector<double> > BUOV(vector<vector<double> > &temp) {
  vector<vector<double> > ans(N+2, vector<double>(N+2, 0));
  for(int i=N; i>=0; i--) {
    for(int j=N; j>=0; j--) {
      ans.at(i).at(j) = Ra * Pr * (temp.at(i).at(j) + temp.at(i).at(j+1)) / 2;
    }
  }
  return ans;
}


vector<vector<double> > CNVT(vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &temp) {
  vector<vector<double> > ans(N+2, vector<double>(N+2, 0));
  for(int i=N; i>=1; i--) {
    for(int j=N; j>=1; j--) {
      double uave = (u.at(i-1).at(j) + u.at(i).at(j)) / 2;
      double vave = (v.at(i).at(j-1) + v.at(i).at(j)) / 2;
      ans.at(i).at(j) =
	uave * (temp.at(i+1).at(j) - temp.at(i-1).at(j)) / (2 * dx) -
	abs(uave) * (temp.at(i-1).at(j) - 2 * temp.at(i).at(j) + temp.at(i+1).at(j)) / (2 * dx) +
	vave * (temp.at(i).at(j+1) - temp.at(i).at(j-1)) / (2 * dy) -
	abs(vave) * (temp.at(i).at(j-1) - 2 * temp.at(i).at(j) + temp.at(i).at(j+1)) / (2 * dy);
    }
  }
  return ans;
}

vector<vector<double> > DIFT(vector<vector<double> > &temp) {
  vector<vector<double> > ans(N+2, vector<double>(N+2, 0));
  for(int i=N; i>=1; i--) {
    for(int j=N; j>=1; j--) {
      ans.at(i).at(j) =
	(temp.at(i-1).at(j) - 2 * temp.at(i).at(j) + temp.at(i+1).at(j)) / (dx * dx) +
	(temp.at(i).at(j-1) - 2 * temp.at(i).at(j) + temp.at(i).at(j+1)) / (dy * dy);
    }
  }
  return ans;
}

double make_b(int i, int j,
	    vector<vector<double> > &div,
	    vector<vector<double> > &cnvu,
	    vector<vector<double> > &difu,
	    vector<vector<double> > &cnvv,
	    vector<vector<double> > &difv,
	    vector<vector<double> > &buov) {
  return (div.at(i).at(j) / dt +
	 (cnvu.at(i-1).at(j) - cnvu.at(i).at(j) + difu.at(i).at(j) - difu.at(i-1).at(j)) / dx +
	 (cnvv.at(i).at(j-1) - cnvv.at(i).at(j) + difv.at(i).at(j) - difv.at(i).at(j-1) + buov.at(i).at(j) - buov.at(i).at(j-1)) / dy) * (dx * dx);
}

void make_A(vector<vector<double> > &A, int i, int j, int cnt) {
  int check = -4;
  if (i-1 > 0 && !(i==2 && j==1)) A.at(cnt).at(place(i-1,j)) = 1;
  else check++;
  
  if (j-1 > 0 && !(i==1 && j==2)) A.at(cnt).at(place(i,j-1)) = 1;
  else check++;
  
  if (i+1 < N+1) A.at(cnt).at(place(i+1,j)) = 1;
  else check++;

  if (j+1 < N+1) A.at(cnt).at(place(i,j+1)) = 1;
  else check++;

  if (i * j == 2) check--;
  A.at(cnt).at(place(i,j)) = check;
  return;
}

// P_{0,0}, P_{0,1},...
//P_{i,j}は上からP.at((i*(N+2))+j)
vector<vector<double> > nextP(vector<vector<double> > &div,
			      vector<vector<double> > &cnvu,
			      vector<vector<double> > &difu,
			      vector<vector<double> > &cnvv,
			      vector<vector<double> > &difv,
			      vector<vector<double> > &buov) {
  vector<vector<double> > A((N-1)*(N+1), vector<double>((N-1)*(N+1),0));
  vector<double> b((N-1)*(N+1));

  // P_{1,1} = 0
  // P_{0,i} = P_{1,i}
  // P_{i,0} = P_{i,1}
  // P_{i,N} = P_{i,N+1}
  // P_{N,i} = P_{N+1,i}
  // P_{i-1,j} + P_{i,j-1} + P_{i+1,j} + P_{i,j+1} - 4 * P_{i,j} = ...
  // P_{1,2} ... P_{N,N} についての方程式
  int cnt = 0;
  for (int i = 1; i <= N; i++) {
    for (int j = 1; j <= N; j++) {
      if (i == 1 && j == 1) continue;
      make_A(A, i, j, cnt);
      b.at(cnt) = make_b(i, j, div, cnvu, difu, cnvv, difv, buov);
      cnt++;
    }
  }
  vector<double> x = SOR(A, b);
  
  vector<vector<double> > ans(N+2, vector<double>(N+2,0));
  // P_{0,2} ... P_{0,N} = P_{1,2} ... P_{1,N}
  for (int i = 2; i <= N; i++) {
    ans.at(0).at(i) = x.at(place(1,i));
  }
  
  for (int i = 1; i <= N; i++) {
    if (i > 1) ans.at(i).at(0) = x.at(place(i,1));
    for (int j = 1; j <= N; j++) {
      if (i * j != 1) ans.at(i).at(j) = x.at(place(i,j));
    }
    ans.at(i).at(N+1) = x.at(place(i,N));
  }

  for (int i = 1; i <= N; i++) {
    ans.at(N+1).at(i) = x.at(place(N,i));
  }
  // 
  return ans;
}

vector<vector<double> > nextU(vector<vector<double> > &p,
			      vector<vector<double> > &u,
			      vector<vector<double> > &cnvu,
			      vector<vector<double> > &difu) {
  vector<vector<double> > ans(N+2, vector<double>(N+2,0));
  for (int i=1; i<N+1; i++) {
    for (int j=0; j<N+2; j++) {
      ans.at(i).at(j) = u.at(i).at(j) + (- (p.at(i+1).at(j) - p.at(i).at(j)) / dx + difu.at(i).at(j) - cnvu.at(i).at(j)) * dt;
    }
  }
  for (int i=0; i<N+2; i++) {
    ans.at(i).at(0) = - ans.at(i).at(1);
    ans.at(i).at(N+1) = - ans.at(i).at(N);
  }
  return ans;
}

vector<vector<double> > nextV(vector<vector<double> > &p,
			      vector<vector<double> > &v,
			      vector<vector<double> > &cnvv,
			      vector<vector<double> > &difv,
			      vector<vector<double> > &buov) {
  vector<vector<double> > ans(N+2, vector<double>(N+2,0));
  for (int i=0; i<N+2; i++) {
    for (int j=1; j<N+1; j++) {
      ans.at(i).at(j) = v.at(i).at(j) + (- (p.at(i).at(j+1) - p.at(i).at(j)) / dy + difv.at(i).at(j) - cnvv.at(i).at(j) + buov.at(i).at(j)) * dt;
    }
  }
  for (int i=0; i<N+2; i++) {
    ans.at(0).at(i) = - ans.at(1).at(i);
    ans.at(N+1).at(i) = - ans.at(N).at(i);
  }
  return ans;
}

vector<vector<double> > nexttemp(vector<vector<double> > &temp,
				 vector<vector<double> > &cnvt,
				 vector<vector<double> > &dift) {
  vector<vector<double> > ans(N+2, vector<double>(N+2,0));
  for (int i=0; i<N+2; i++) {
    for (int j=0; j<N+2; j++) {
      ans.at(i).at(j) = (dift.at(i).at(j) - cnvt.at(i).at(j)) * dt + temp.at(i).at(j);
    }
  }
  for (int i=0; i<N+2; i++) {
    ans.at(0).at(i) = 1 - ans.at(1).at(i);
    ans.at(N+1).at(i) = -1 - ans.at(N).at(i);
  }
  return ans;
}

vector<vector<vector<double> > > ITER(vector<vector<vector<double> > > &uvtemp) {
  vector<vector<double> > u = uvtemp.at(0), v = uvtemp.at(1), temp = uvtemp.at(2);
  vector<vector<double> > cnvu = CNVU(u, v);
  vector<vector<double> > div = DIV(u, v);
  vector<vector<double> > difu = DIFU(u);
  vector<vector<double> > cnvv = CNVV(u, v);
  vector<vector<double> > difv = DIFV(v);
  vector<vector<double> > buov = BUOV(temp);
  vector<vector<double> > cnvt = CNVT(u, v, temp);
  vector<vector<double> > dift = DIFT(temp);
  vector<vector<double> > p = nextP(div, cnvu, difu, cnvv, difv, buov);

  vector<vector<vector<double> > > ans(3, vector<vector<double> >(N+2, vector<double>(N+2)));
  ans.at(0) =  nextU(p, u, cnvu, difu);
  ans.at(1) = nextV(p, v, cnvv, difv, buov);
  ans.at(2) = nexttemp(temp, cnvt, dift);
  return ans;
}

void show_param(vector<vector<vector<double> > > &uvtemp) {
  for (int i = 0; i < N+2; i++) {
    for(int j = N+1; j>=0; j--) {
      cout << fixed << setprecision(2) << uvtemp.at(0).at(j).at(i) << ' ';
    }
    cout << endl;
  }
  cout << endl;
  for (int i = 0; i < N+2; i++) {
    for(int j = N+1; j>=0; j--) {
      cout << fixed << setprecision(2) << uvtemp.at(1).at(j).at(i) << ' ';
    }
    cout << endl;
  }
  cout << endl;
  for (int i = 0; i < N+2; i++) {
    for(int j = N+1; j >=0; j--) {
      cout << fixed << setprecision(2) << uvtemp.at(2).at(j).at(i) << ' ';
    }
    cout << endl;
  }
  return;
}

int main(int argc, char *argv[]) {
  int iter;
  if (argc != 2) {
    cout << "Num of iter!" << endl;
    return 0;
  } else {
    iter = stoi(argv[1]);
  }
  
  vector<vector<vector<double> > > uvtemp(3, vector<vector<double> >(N+2, vector<double>(N+2, 0)));
  for (int i=0; i<N+2; i++) {
    uvtemp.at(2).at(0).at(i) = 0.5;
    uvtemp.at(2).at(1).at(i) = 0.5;
    uvtemp.at(2).at(N).at(i) = -0.5;
    uvtemp.at(2).at(N+1).at(i) = -0.5;
  }

  for (int i = 0; i < iter; i++) {
    // cerr << i << endl;
    cout << "--------------------------------------- NUM :" << i << "-----------------------" << endl;
    show_param(uvtemp);
    cout << endl;
    vector<vector<vector<double> > > uvtemp1 = ITER(uvtemp);
    swap(uvtemp1, uvtemp);
  }
  return 0;
}
