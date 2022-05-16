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
#define N 100
#define dx 1e-2
#define dy 1e-2
#define Ra 0.71
#define Pe 7100
#define dt 1e-4

vector<vector<double> > DIV(vector<vector<double> > u, vector<vector<double> > v) {
  vector<vector<double> > ans(N+2, vector<double>(N+2, 0));
  for(int i=N+1; i>=1; i--) {
    for(int j=N+1; j>=1; j--) {
      ans.at(i).at(j) = (u.at(i).at(j)-u.at(i-1).at(j))/dx + (v.at(i).at(j)-v.at(i).at(j-1))/dy;
    }
  }
  return ans;
}

vector<vector<double> > CNVU(vector<vector<double> > u, vector<vector<double> > v) {
  vector<vector<double> > ans(N+2, vector<double>(N+2, 0));
  for(int i=N+1; i>=1; i--) {
    for(int j=N+1; j>=1; j--) {
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
	
