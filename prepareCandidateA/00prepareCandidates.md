### 概要

骨格分子に官能基を修飾することで、候補物質を作成します。

### 1. 官能基の自動作成、及び、候補物質の作成

計算コードは bin にあるので、相対位置は読み替えてください。

溶解度データセット solubility.train.sdfに含まれる物質の側鎖を取り出し、指定した物質を修飾して、候補物質を作成します。
主に、吸収スペクトルの変化に用いる官能基のスクリーニングを行います。

分子骨格、及び、官能基で修飾する部位を指定する設定ファイルを準備します。
```
>>cat molecularBackbone.prm
furan-　#出力ファイル内の候補物質の頭文字
c1occc1 # 骨格分子のSMILES構造式
c1oc([*:1])cc1 # 官能基で修飾する部位を指定する式
WeighFromTo:0,200 # （水素を除く官能基の）分子量の範囲の指定
AtomNum:3,5 # (水素を除く官能基中の）原子の数の指定
```
下記コマンドで、官能基、及び、候補物質を作成します。
```
>> python3 ../bin/mkCandidatesFromMolecularBackbone.py molecularBackbone.prm
```

官能基で修飾した最大60の候補分子を molecularBackbone.smi に書き込みます。
確認用に、分子構造の画像ファイル(molecularBackbone_?.png/molecularBackbone_A?.png)も作成されます。

```
>>cat molecularBackbone.smi
furan-000 'c1ccoc1'
furan-001 'CC(=O)Sc1ccco1'
furan-002 'CS(=O)c1ccco1'
furan-003 'NS(=O)(=O)c1ccco1'
furan-004 'ClC(Cl)c1ccco1’
...
```
作成した官能基の一覧
<img width="1031" height="773" alt="image" src="https://github.com/user-attachments/assets/ef4560db-cc0a-4c3b-ba08-9943d268cba9" />

候補物質の一覧

<img width="687" height="859" alt="image" src="https://github.com/user-attachments/assets/290d0c21-82c5-4f7c-9942-948f915e35d0" />

### 2. 吸収スペクトルの計算

候補物質を一括処理して、吸収スペクトルを計算するための実行スクリプトを作ります。
```
>>python3 ../bin/mkScript.py molecularBackbone.smi > runAll.sh
```
スクリプトを実行します。
```
>>sh runAll.sh
```

解析結果の可視化を行います。

結果が得られた条件を表示するためのパラメータファイルを作成します。
```
>> sh ../bin/mkPrmForTdds2svg.sh > plot.prm
```
下記の例では、波長250nm 以上に有効な吸収ピークを持つ物質を選別します。ここでは、有効な吸収ピークは最大の吸収ピークとの相対値で10%以上としています。
```
>>python3 ../bin/tdds2svg.py plot.prm -up 250 -r 0.1
```
uv_preA_plot_all.pptx が作られます。

<img width="1037" height="624" alt="image" src="https://github.com/user-attachments/assets/80580323-d7ce-4444-b6c4-02db7d46469a" />

オプションを付けなければ、全候補物資のスペクトルが表示されます。

各々のパラメータファイルはエディターで編集可能なので、必要な場合は、ディレクトリー tddft　内のマニュアルを参照に修正してください。
