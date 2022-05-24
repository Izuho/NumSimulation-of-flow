# NumSimulation-of-flow
## 「流れの数値計算と可視化」第3版　平野博之著
- MAC.cpp の実行ファイルの作り方<br>
`g++ -O2 -o MAC MAC.cpp -Wall -isystem $HOME/eigen-3.4.0`
- 実行の仕方 （いつまで計算するか整数を指定）<br>
`./MAC 1000 > outer.dat`
- 可視化の仕方<br>
`ShowFlow.ipynb`に処理をさせれば`test.gif`ができる

## 今後
- 性能評価をするためにpapi performanceを使って、flopsを計算する。
- 理想性能よりも小さいはずなのでスレッド並列化などをして計算時間を減らす。
- 単純な場合で誤り確認をしておく。
