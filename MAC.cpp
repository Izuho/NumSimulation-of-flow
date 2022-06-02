/*
左の壁の温度が -0.5、右の壁の温度が 0.5 になるように

最初液体温度は一様に 0 とする。
壁で流速は 0

圧力は (1,1) において常に 0 とする。
レイノルズ数:0.71
プラントル数:7100
dt:1e-4
dx=1/N
dy=1/N

 */
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <omp.h>

using namespace std;
#define n 4
// = 2^n
#define N 16
// N に合わせる
#define dx 0.0625 // 1 / N
#define dy 0.0625
#define Ra 7100
#define Pr 0.71
#define dt 1e-4
#define _dt 1e4 // 1 / dt
#define place(i,j) (i-1)*N+j-2
#define CNT_MAX 200
#define TOL 1e-6
#define omega 1.9

double conv(vector<double> a, vector<double> b) {
    double up = 0, down = 0;
    int M = a.size();
    for (int i = 0; i < M; i++) {
        up += (a.at(i) - b.at(i)) * (a.at(i) - b.at(i));
        down += b.at(i) * b.at(i);
    }
    return up / down;
}

int A2a(int i, int j) {
    int x = ((i + 1) >> n) + 1, y = ((i + 1) & (N - 1)) + 1;
    int a = ((j + 1) >> n) + 1, b = ((j + 1) & (N - 1)) + 1;
    int ans;
    switch (abs(x - a) + abs(y - b)) {
    case 1:
        ans = 1;
        break;
    case 0:
        ans = -4;
        if (x <= 1) ans += 1;

        if (y <= 1) ans += 1;

        if (x >= N) ans += 1;

        if (y >= N) ans += 1;
        break;
    default:
        ans = 0;
    }
    return ans;
}

vector<double> SOR(vector<double> b) {
    int M = b.size();
    vector<double> x(M, nan("")), x_past(M, 0);
    for (int k = 0; k < CNT_MAX; k++) {
        for (int i = 0; i < M; i++) {
            double sigma = 0;

            #pragma omp parallel for reduction(+:sigma)
            for (int j = max(0, i - N); j < min(M, i + N + 1); j++) {
                if (j < i)      sigma += A2a(i, j) * x.at(j);
                else if (j > i) sigma += A2a(i, j) * x_past.at(j);
                else continue;
            }

            sigma = (b.at(i) - sigma) / A2a(i, i);
            x.at(i) = x_past.at(i) + omega * (sigma - x_past.at(i));
        }
        swap(x, x_past);
        if (conv(x, x_past) < TOL) break;
    }
    return x_past;
}

vector<vector<double> > DIV(vector<vector<double> >& u, vector<vector<double> >& v) {
    vector<vector<double> > ans(N + 2, vector<double>(N + 2, 0));
    for (int i = N + 1; i >= 1; i--) {
        for (int j = N + 1; j >= 1; j--) {
            ans.at(i).at(j) =
                N *
                (u.at(i).at(j) - u.at(i - 1).at(j) + v.at(i).at(j) - v.at(i).at(j - 1));
        }
    }
    return ans;
}

vector<vector<double> > CNVU(vector<vector<double> >& u, vector<vector<double> >& v) {
    vector<vector<double> > ans(N + 2, vector<double>(N + 2, 0));
    for (int i = N; i >= 1; i--) {
        for (int j = N; j >= 1; j--) {
            double vave =
                (v.at(i).at(j) +
                    v.at(i).at(j - 1) +
                    v.at(i + 1).at(j) +
                    v.at(i + 1).at(j - 1)) * 0.25;
            ans.at(i).at(j) =
                0.5 * N *
                (u.at(i).at(j) * (u.at(i + 1).at(j) - u.at(i - 1).at(j)) -
                    abs(u.at(i).at(j)) * (u.at(i - 1).at(j) - 2 * u.at(i).at(j) + u.at(i + 1).at(j)) +
                    vave * (u.at(i).at(j + 1) - u.at(i).at(j - 1)) -
                    abs(vave) * (u.at(i).at(j - 1) - 2 * u.at(i).at(j) + u.at(i).at(j + 1)));
        }
    }
    return ans;
}

vector<vector<double> > DIFU(vector<vector<double> >& u) {
    vector<vector<double> > ans(N + 2, vector<double>(N + 2, 0));
    for (int i = N; i >= 1; i--) {
        for (int j = N; j >= 1; j--) {
            ans.at(i).at(j) =
                Pr * N * N *
                (u.at(i - 1).at(j) + u.at(i + 1).at(j) + u.at(i).at(j - 1) + u.at(i).at(j + 1) - 4 * u.at(i).at(j));
        }
    }
    return ans;
}

vector<vector<double> > CNVV(vector<vector<double> >& u, vector<vector<double> >& v) {
    vector<vector<double> > ans(N + 2, vector<double>(N + 2, 0));
    for (int i = N; i >= 1; i--) {
        for (int j = N; j >= 1; j--) {
            double uave =
                (u.at(i).at(j) +
                    u.at(i).at(j - 1) +
                    u.at(i + 1).at(j) +
                    u.at(i + 1).at(j - 1)) * 0.25;
            ans.at(i).at(j) =
                0.5 * N *
                (uave * (v.at(i + 1).at(j) - v.at(i - 1).at(j)) -
                    abs(uave) * (v.at(i - 1).at(j) - 2 * v.at(i).at(j) + v.at(i + 1).at(j)) +
                    v.at(i).at(j) * (v.at(i).at(j + 1) - v.at(i).at(j - 1)) -
                    abs(v.at(i).at(j)) * (v.at(i).at(j - 1) - 2 * v.at(i).at(j) + v.at(i).at(j + 1)));
        }
    }
    return ans;
}

vector<vector<double> > DIFV(vector<vector<double> >& v) {
    vector<vector<double> > ans(N + 2, vector<double>(N + 2, 0));
    for (int i = N; i >= 1; i--) {
        for (int j = N; j >= 1; j--) {
            ans.at(i).at(j) =
                Pr * N * N *
                (v.at(i - 1).at(j) - 2 * v.at(i).at(j) + v.at(i + 1).at(j) +
                    v.at(i).at(j - 1) - 2 * v.at(i).at(j) + v.at(i).at(j + 1));
        }
    }
    return ans;
}

vector<vector<double> > BUOV(vector<vector<double> >& temp) {
    vector<vector<double> > ans(N + 2, vector<double>(N + 2, 0));
    for (int i = N; i >= 0; i--) {
        for (int j = N; j >= 0; j--) {
            ans.at(i).at(j) = Ra * Pr * 0.5 * (temp.at(i).at(j) + temp.at(i).at(j + 1));
        }
    }
    return ans;
}


vector<vector<double> > CNVT(vector<vector<double> >& u, vector<vector<double> >& v, vector<vector<double> >& temp) {
    vector<vector<double> > ans(N + 2, vector<double>(N + 2, 0));
    for (int i = N; i >= 1; i--) {
        for (int j = N; j >= 1; j--) {
            double uave = (u.at(i - 1).at(j) + u.at(i).at(j)) * 0.5;
            double vave = (v.at(i).at(j - 1) + v.at(i).at(j)) * 0.5;
            ans.at(i).at(j) =
                0.5 * N *
                (uave * (temp.at(i + 1).at(j) - temp.at(i - 1).at(j)) -
                    abs(uave) * (temp.at(i - 1).at(j) - 2 * temp.at(i).at(j) + temp.at(i + 1).at(j)) +
                    vave * (temp.at(i).at(j + 1) - temp.at(i).at(j - 1)) -
                    abs(vave) * (temp.at(i).at(j - 1) - 2 * temp.at(i).at(j) + temp.at(i).at(j + 1)));
        }
    }
    return ans;
}

vector<vector<double> > DIFT(vector<vector<double> >& temp) {
    vector<vector<double> > ans(N + 2, vector<double>(N + 2, 0));
    for (int i = N; i >= 1; i--) {
        for (int j = N; j >= 1; j--) {
            ans.at(i).at(j) =
                N * N *
                (temp.at(i - 1).at(j) - 2 * temp.at(i).at(j) + temp.at(i + 1).at(j) +
                    temp.at(i).at(j - 1) - 2 * temp.at(i).at(j) + temp.at(i).at(j + 1));
        }
    }
    return ans;
}

double make_b(int i, int j,
    vector<vector<double> >& div,
    vector<vector<double> >& cnvu,
    vector<vector<double> >& difu,
    vector<vector<double> >& cnvv,
    vector<vector<double> >& difv,
    vector<vector<double> >& buov) {
    return (div.at(i).at(j) * _dt +
        (cnvu.at(i - 1).at(j) - cnvu.at(i).at(j) + difu.at(i).at(j) - difu.at(i - 1).at(j)) * N +
        (cnvv.at(i).at(j - 1) - cnvv.at(i).at(j) + difv.at(i).at(j) - difv.at(i).at(j - 1) + buov.at(i).at(j) - buov.at(i).at(j - 1)) * N) * (dx * dx);
}

// P_{0,0}, P_{0,1},...
//P_{i,j}は上からP.at((i*(N+2))+j)
vector<vector<double> > nextP(vector<vector<double> >& div,
    vector<vector<double> >& cnvu,
    vector<vector<double> >& difu,
    vector<vector<double> >& cnvv,
    vector<vector<double> >& difv,
    vector<vector<double> >& buov) {
    vector<double> b((N - 1) * (N + 1));

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
            b.at(cnt) = make_b(i, j, div, cnvu, difu, cnvv, difv, buov);
            cnt++;
        }
    }

    vector<double> x = SOR(b);

    vector<vector<double> > ans(N + 2, vector<double>(N + 2, 0));
    // P_{0,2} ... P_{0,N} = P_{1,2} ... P_{1,N}
    for (int i = 2; i <= N; i++) {
        ans.at(0).at(i) = x.at(place(1, i));
    }

    for (int i = 1; i <= N; i++) {
        if (i > 1) ans.at(i).at(0) = x.at(place(i, 1));
        for (int j = 1; j <= N; j++) {
            if (i * j != 1) ans.at(i).at(j) = x.at(place(i, j));
        }
        ans.at(i).at(N + 1) = x.at(place(i, N));
    }

    for (int i = 1; i <= N; i++) {
        ans.at(N + 1).at(i) = x.at(place(N, i));
    }

    return ans;
}

vector<vector<double> > nextU(vector<vector<double> >& p,
    vector<vector<double> >& u,
    vector<vector<double> >& cnvu,
    vector<vector<double> >& difu) {
    vector<vector<double> > ans(N + 2, vector<double>(N + 2, 0));
    for (int i = 1; i < N + 1; i++) {
        for (int j = 0; j < N + 2; j++) {
            ans.at(i).at(j) = u.at(i).at(j) + (-(p.at(i + 1).at(j) - p.at(i).at(j)) * N + difu.at(i).at(j) - cnvu.at(i).at(j)) * dt;
        }
    }
    for (int i = 0; i < N + 2; i++) {
        ans.at(i).at(0) = -ans.at(i).at(1);
        ans.at(i).at(N + 1) = -ans.at(i).at(N);
    }
    return ans;
}

vector<vector<double> > nextV(vector<vector<double> >& p,
    vector<vector<double> >& v,
    vector<vector<double> >& cnvv,
    vector<vector<double> >& difv,
    vector<vector<double> >& buov) {
    vector<vector<double> > ans(N + 2, vector<double>(N + 2, 0));
    for (int i = 0; i < N + 2; i++) {
        for (int j = 1; j < N + 1; j++) {
            ans.at(i).at(j) = v.at(i).at(j) + (-(p.at(i).at(j + 1) - p.at(i).at(j)) / dy + difv.at(i).at(j) - cnvv.at(i).at(j) + buov.at(i).at(j)) * dt;
        }
    }
    for (int i = 0; i < N + 2; i++) {
        ans.at(0).at(i) = -ans.at(1).at(i);
        ans.at(N + 1).at(i) = -ans.at(N).at(i);
    }
    return ans;
}

vector<vector<double> > nexttemp(vector<vector<double> >& temp,
    vector<vector<double> >& cnvt,
    vector<vector<double> >& dift) {
    vector<vector<double> > ans(N + 2, vector<double>(N + 2, 0));
    for (int i = 0; i < N + 2; i++) {
        for (int j = 0; j < N + 2; j++) {
            ans.at(i).at(j) = (dift.at(i).at(j) - cnvt.at(i).at(j)) * dt + temp.at(i).at(j);
        }
    }
    for (int i = 0; i < N + 2; i++) {
        ans.at(0).at(i) = 1 - ans.at(1).at(i);
        ans.at(N + 1).at(i) = -1 - ans.at(N).at(i);
    }
    return ans;
}

vector<vector<vector<double> > > ITER(vector<vector<vector<double> > >& uvtemp) {
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

    vector<vector<vector<double> > > ans(3, vector<vector<double> >(N + 2, vector<double>(N + 2)));
    ans.at(0) = nextU(p, u, cnvu, difu);
    ans.at(1) = nextV(p, v, cnvv, difv, buov);
    ans.at(2) = nexttemp(temp, cnvt, dift);
    return ans;
}

void show_param(vector<vector<vector<double> > >& uvtemp) {
    for (int i = 0; i < N + 2; i++) {
        for (int j = N + 1; j >= 0; j--) {
            cout << fixed << setprecision(2) << uvtemp.at(0).at(j).at(i) << ' ';
        }
        cout << endl;
    }
    cout << endl;
    for (int i = 0; i < N + 2; i++) {
        for (int j = N + 1; j >= 0; j--) {
            cout << fixed << setprecision(2) << uvtemp.at(1).at(j).at(i) << ' ';
        }
        cout << endl;
    }
    cout << endl;
    for (int i = 0; i < N + 2; i++) {
        for (int j = N + 1; j >= 0; j--) {
            cout << fixed << setprecision(2) << uvtemp.at(2).at(j).at(i) << ' ';
        }
        cout << endl;
    }
    return;
}

int main(int argc, char* argv[]) {
    int iter;
    if (argc != 2) {
        cout << "Num of iter!" << endl;
        return 0;
    }
    else {
        iter = stoi(argv[1]);
    }

    vector<vector<vector<double> > > uvtemp(3, vector<vector<double> >(N + 2, vector<double>(N + 2, 0)));
    for (int i = 0; i < N + 2; i++) {
        uvtemp.at(2).at(0).at(i) = 0.5;
        uvtemp.at(2).at(1).at(i) = 0.5;
        uvtemp.at(2).at(N).at(i) = -0.5;
        uvtemp.at(2).at(N + 1).at(i) = -0.5;
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