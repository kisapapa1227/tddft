## 概要
ESIPT(励起状態分子内プロトン移動）を模して、量子化学計算を行います。
ディレクトリー esipt2 の計算は、現最新バージョンの numpy=2.4.0 ではエラーが出ます。
エラーが出る場合は、numpy の別バージョンで実行してみてください。動作確認済み(pyscf=2.11,numpy=1.26,2.0.0,2.1.0,2.2.0,2.3.0)

### 1. 解析手順
計算コードは bin にあります。相対位置は読み替えてください。

### 1.1 解析概要
<img width="1067" height="594" alt="image" src="https://github.com/user-attachments/assets/043ea0ef-8719-4efc-81f5-444e4d710c8b" />

### 1.2 解析手順

① 基底状態の原子の座標(モルブロック)を求めます。 
```
>>python3 ../bin/esipt.py AAP ‘c1(C(=O)C)ccccc1N’ –phase 1
```
第一引数：物質名（出力ファイル名）、第二引数：対象物質のSMILES構造式

解析結果は、物質名.txt に保存されます。

② 可視化して、対象原子の番号を確認します。

分子構造を投影図にするためのパラメータファイルを準備します。
```
>>cat project.prm
AAP # 物質名
10,19,3 #平面を定義する原子の座標。
flip:True # 裏返す。
rotate:240.0 # 回転させる。
```

パラメータファイルを利用して投影図を作成します。
```
>>python3 ../bin/projectMol.py project.prm
```

対象となる原子のモルブロック内での番号を確認します。
```
>>display AAP_proj.svg
```
必要であれば、対象の箇所が見やすいように、パラメータファイルを書き換えて、投影図を作り直してください。

<img width="349" height="448" alt="image" src="https://github.com/user-attachments/assets/3adbe02b-2ff9-4f7c-b48d-68bde5e629f4" />

このサンプルでは、アミノ基の水素とカルボキシル基内の二重結合を持つ酸素との影響を対象とします。図より、対象の水素と酸素は#19と#3になります。

③対象とする原子間距離を指定して、基底エネルギー、吸収スペクトルを計算します。
```
>>python3 ../bin/esipt.py AAP ‘c1(C(=O)C)ccccc1N’ -phase 2 -distance 3,19,2.05 -pt 1 –ms 100
```
オプションの説明
```
-distance #1,#2,#3:原子 #2の座標を移動させ、原子 #1との距離を#3(Å)にする。
-pt #1: ファイルの枝番号を指定する。この場合は AAP_01.tdd が作られる。
-no-op: 分子の位置を移動させるのみで、最適化を行わない場合に指定します。
-ms #1: 基底状態の最適化の最大反復回数を指定する。
```

３回の反復で、吸収スペクトルのグラフ化に十分な精度が得られますが、基底エネルギーのプロファイルを得る場合には、有効数字６桁程度まで求める必要があり、50～の反復が必要になります。
成果物として、AAP-#N.txt(原子の配置位置)、AAP-#N.tdd(吸収強度）が作成され、AAP.erg に基底状態のエネルギーが追記されます。

ディレクトリー esipt2 には、サンプルとして、H(#3)とN(#10)の距離を拘束条件として計算した結果があります。

④原子間距離を変えて、計算を繰り返します(sample3.sh)。
```
>>python3 ../bin/esipt.py AAP ‘c1(C(=O)C)ccccc1N’ –phase 2 –distance 3,19,2.00 –pt 2
...
>>python3 ../bin/esipt.py AAP ‘c1(C(=O)C)ccccc1N’ –phase 2 –distance 3,19,1.00 –pt 22
```
枝番号は1からの連続で指定してください。

⑤解析結果のグラフ化、可視化(動画作成)を行う(sample4.sh)。
```
>>sh ../bin/esiptTxt2Gif.sh project.prm
```
解析手順②で利用した設定ファイル project.prm を利用して、AAP-#n.tdd ファイルより吸収スペクトル、APP-#n.txt ファイルより投影図を作成します。
最終成果物として、AAP_esipt_cb.gif (動画)が作成されます。最適化の過程で分子の座標が動いてしまっています。


![AAP_esipt_cb](https://github.com/user-attachments/assets/1d6d4072-b1af-424a-adf3-a33c6d16a2e6)


ディレクトリー esipt2 では、H(#3)とN(#10)の距離を拘束条件とします。

![AAP_esipt2_cb](https://github.com/user-attachments/assets/52aa8290-cb66-4c8e-b8d5-f1190dbbf53c)

次に、基底状態のエネルギーをグラフ化します。
パラメタファイルを準備します。
```
>> cat plotEnergy.prm
files:../esipt/AAP,../esipt2/AAP
label:'O H','N-H'
output:AAP_erg
ops:[0.9,2.1],60
```

ポテンシャルエネルギーをグラフ化します。
```
>> python3 ../bin/esiptErgPlot.py plotEnergy.prm
```


<img width="640" height="480" alt="AAP_erg" src="https://github.com/user-attachments/assets/ee166833-66bb-42d6-bb60-18bb2707a63c" />

### 二面角で指定する

二面角を拘束条件にすることで、一重結合の回転を模した解析が行えます。また、オプションで –phase 3 を指定すると、最適化を行って基底エネルギーを求めたところで計算を打ち切り、解析時間の短縮ができます（吸収スペクトルは求まりません）。
```
>>python3 ../bin/esipt.py substance_0403_4 'O=C(C1=CC=CC=C1)C2=C(NC(N(CC3=CC=CC=C3)CC4=CC=CC=C4)=O)C=CC=C2‘　–phase 3 -dihedral 1,2,9,10,-133 -pt 1 –ms 100
```
オプション -direhdral #1,#2,#3,#4,#5: 原子 #1,#2,#3で定義される面と#2,#3,#4で定義される面の角度を#5にします。

<img width="439" height="430" alt="image" src="https://github.com/user-attachments/assets/80c9385e-2c0a-4aa3-9a23-4dc917b47baf" />

サンプル(dihedral)では 炭素#2-炭素#9間の一重結合が回転した時のポテンシャルエネルギーの変化を求めています。

ポテンシャルエネルギーをグラフ化します。
```
>>python3 ../bin/esiptErgPlot.py plotEnergy.prm
```
<img width="624" height="468" alt="image" src="https://github.com/user-attachments/assets/8e8b2269-620a-4b8a-a842-e225d2c70a00" />

### 補足

<img width="1128" height="634" alt="image" src="https://github.com/user-attachments/assets/ed372cba-2338-44d4-9ee1-410303399878" />

