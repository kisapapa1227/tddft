## 概要
ESIPT(励起状態分子内プロトン移動）を模擬して量子化学計算を行います。
ディレクトリー esipt2 の計算は、現時点の numpy　2.4.0 ではエラーが発生します。
古いnumpy のバージョンを試してください。動作確認済みの組み合わせ例
- pyscf=2.11
- numpy=1.26.4, 2.0.0, 2.1.0, 2.2.0, 2.3.0

### 1. 解析手順
計算コードは bin にあります。パスは環境に合わせて読み替えてください。

### 1.1 解析概要
<img width="1067" height="594" alt="image" src="https://github.com/user-attachments/assets/043ea0ef-8719-4efc-81f5-444e4d710c8b" />

### 1.2 解析手順

① 基底状態の原子座標(molブロック)を求めます。 
```
python3 ../bin/esipt.py AAP ‘c1(C(=O)C)ccccc1N’ –phase 1
```
- 第1引数：物質名（出力ファイル名）
- 第2引数：対象物質のSMILES構造式(クォートで囲む)
- 出力: AAP.txt (原子配置等を含むテキスト)


② 投影図(可視化)で対象原子の番号を確認します。

- 投影図作成のためのパラメータファイルを準備します。

例
```
# ファイル project.prm
AAP # 物質名
3,19,10 # 投影面を定義する原子の番号(処理時に先頭2個の原子間距離を表示する)。
flip:False # 裏返し（True/False)
rotate:240.0 # 回転角度(度)
```

- 投影図を作成
```
python3 ../bin/projectMol.py project.prm
...
03-19 : 2.00Å O H
```
この処理は投影図を作るとともに、結合している原子対や指定原子間の距離を標準出力に表示します。

- 作成された SVG を表示して、モルブロック内での対象原子の番号を確認します:
```
display AAP_proj.svg
```
必要に応じて、対象の箇所が見やすいように、パラメータファイルを書き換えて、投影図を作り直してください。

<img width="349" height="448" alt="image" src="https://github.com/user-attachments/assets/3adbe02b-2ff9-4f7c-b48d-68bde5e629f4" />

このサンプルでは、「アミノ基の水素」と「カルボキシル基内の炭素と二重結合する酸素」との影響を見るため、図より対象の水素と酸素は#19と#3に対応します。

③原子間距離を指定して基底最適化・吸収スペクトルを計算します。
- 例：
```
python3 ../bin/esipt.py AAP ‘c1(C(=O)C)ccccc1N’ -phase 2 -distance 3,19,2.05 -pt 1 –ms 100
```
オプションの説明
```
-distance #1,#2,#3:原子 #2の座標を移動させ、原子 #1との距離を#3(Å)にする
-pt #1: 出力ファイルの枝番号(例:AAP_01.tdd が作られる)
-no-op: 分子位置を移動するのみで最適化を行わない
-ms #1: 基底状態最適化の最大反復回数を指定(デフォルト 3)
```

運用上の目安：
- 吸収スペクトルの作成に必要な精度であれば数回（例：３回）の反復でも十分なことが多いです。
- 基底エネルギーのプロファイルを精密に求める場合（有効数字6桁程度）は、50回以上の反復が必要になります。

出力物:
- AAP-#N.txt (原子の配置位置)
- AAP-#N.tdd (吸収強度データ）
- AAP.erg (基底エネルギーを追記)

注:ディレクトリー esipt2 には、サンプルとして、H(#3)とN(#10)の距離を拘束した計算結果が含まれています。

④原子間距離を変えて、計算を繰り返します(sample3.sh)。

距離を段階的に変えて繰り返し計算します。例として、投影図で得られた、H(#3)とO(#19)の距離 2.00Å から単結合相当の1.00Åまで近づけていく場:

```
python3 ../bin/esipt.py AAP ‘c1(C(=O)C)ccccc1N’ –phase 2 –distance 3,19,2.00 –pt 2
...
python3 ../bin/esipt.py AAP ‘c1(C(=O)C)ccccc1N’ –phase 2 –distance 3,19,1.00 –pt 22
```
- -ptの枝番号は1からの連続で指定してください。

⑤解析結果のグラフ化、可視化(動画作成)を行う(sample4.sh)。
- 一連の .tdd/.txt ファイルから GIF を作成します:
```
sh ../bin/esiptTxt2Gif.sh project.prm
```
- project.prm (解析手順②で作成したもの)を用いて、.tdd から吸収スペクトル、.txt から投影図を作成し、最終的に AAP_esipt_cb.gif (例)が作成されます。
- 最適化の過程で分子座標が変化するため、アニメーションでは構造が動くことがあることに注意してください。


![AAP_esipt_cb](https://github.com/user-attachments/assets/1d6d4072-b1af-424a-adf3-a33c6d16a2e6)

(補足)
- ディレクトリー esipt2 内のサンプルでは、H(#3)とN(#10)の距離を拘束して計算しています。

![AAP_esipt2_cb](https://github.com/user-attachments/assets/52aa8290-cb66-4c8e-b8d5-f1190dbbf53c)


### 基底状態エネルギー(ポテンシャル曲線)の作図
プロット用のパラメタファイルの例:
```
# ファイル plotEnergy.prm
files:../esipt/AAP,../esipt2/AAP # データのあるディレクトリーをカンマ区切りで指定
label:'O H','N-H' # 凡例のラベル（カンマ区切り）
output:AAP_erg #　出力ファイルの先頭
ops:[0.9,2.1],60 # 横軸 0.9-2.1 Å、縦軸 0-60 kcal/mol
```
- プロット実行:
```
python3 ../bin/esiptErgPlot.py plotEnergy.prm
```


<img width="640" height="480" alt="AAP_erg" src="https://github.com/user-attachments/assets/ee166833-66bb-42d6-bb60-18bb2707a63c" />

### 二面角(ジヘドラル)指定による解析
- 一重結合の回転を模した解析（二面角を拘束）も可能です。-phase 3 を指定すると、基底状態エネルギーの最適化のみで計算を打ち切り、解析時間を短縮できます（このモードでは吸収スペクトルは生成されません）。

例
```
python3 ../bin/esipt.py substance_0403_4 'O=C(C1=CC=CC=C1)C2=C(NC(N(CC3=CC=CC=C3)CC4=CC=CC=C4)=O)C=CC=C2‘　–phase 3 -dihedral 1,2,9,10,-133 -pt 1 –ms 100
```
- オプション -direhdral #1,#2,#3,#4,#angle: 原子 #1,#2,#3で定義される面と、#2,#3,#4で定義される面のなす角を #angle(度) に設定

<img width="439" height="430" alt="image" src="https://github.com/user-attachments/assets/80c9385e-2c0a-4aa3-9a23-4dc917b47baf" />

サンプル(dihedral)では C#2-C#9間の一重結合が回転した時のポテンシャルエネルギーの変化を評価しています。結果はプロットも esiptErgPlot.py を用いて作成します。

プロット用のパラメタファイルの例
```
# ファイル plotEnergy.prm
files:substance_0403_4
label:substance_0403
outpu:tsubstance_0403_erg
ops:[-140,210],20
axis:dihedral angle(deg)
```

```
python3 ../bin/esiptErgPlot.py plotEnergy.prm
```
<img width="624" height="468" alt="image" src="https://github.com/user-attachments/assets/8e8b2269-620a-4b8a-a842-e225d2c70a00" />

### 補足
- 図や生成物（SVG、PNG、GIF 等）はプロジェクト内の対応ファイルを参照してください。
- コマンド中の引用符やダッシュ（– と - の混在）に注意してください。CLI オプションには必ず半角ハイフン（-）を使用してください。
- 計算中に数値的な問題が発生する場合は、NumPy / pyscf のバージョンを変更して試してください（冒頭の互換性を参照）。
<img width="1128" height="634" alt="image" src="https://github.com/user-attachments/assets/ed372cba-2338-44d4-9ee1-410303399878" />

