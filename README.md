# NumSimulation-of-flow
## 「流れの数値計算と可視化」第3版　平野博之著
- MAC.cpp の実行ファイルの作り方<br>
`g++ -O2 -o MAC MAC.cpp -Wall -isystem $HOME/eigen-3.4.0`
- 実行の仕方 （いつまで計算するか整数を指定）<br>
`./MAC 1000 > outer.dat`
- 可視化の仕方<br>
`ShowFlow.ipynb`に処理をさせれば`test.gif`ができる

## 正しさの保証
等温で回転があるという明らかな場合でのシミュレーションに成功(test-for-check.gif)

## 性能評価の方法
- MacBook Air (Retina, 13-inch, 2018)
- プロセッサ　1.6 GHz デュアルコアIntel Core i5
- メモリ　8 GB 2133 MHz LPDDR3
- `instruments -t "Time Profiler" ./MAC 1000`からの `open instrumentscli0.trace`で解析

## 今後
- 性能評価をするためにpapi performanceを使って、flopsを計算する。
- 理想性能よりも小さいはずなのでスレッド並列化などをして計算時間を減らす。
