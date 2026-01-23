### 概要

骨格分子に自動生成した側鎖を追加することで、候補物質を作成します。

### 1. 官能基の自動作成、及び、候補物質の作成

計算コードは bin にあるので、相対位置は読み替えてください。

溶解度データセット solubility.train.sdfに含まれる物質を利用して、側鎖を作ります。作成した側鎖で対象物質を修飾して、候補物質を作成します。

骨格分子、及び、修飾する部位を指定する設定ファイルを準備します。
```
>>cat molecularBackbone.prm
furan-　# 出力ファイル内の候補物質名の頭文字
c1occc1 # 骨格分子のSMILES構造式
c1oc([*:1])cc1 # 官能基で修飾する部位を指定する式
WeighFromTo:0,200 # （官能基の水素を除く）分子量の範囲の指定
AtomNum:3,5 # (官能基中の水素を除く）原子の数の指定
```
下記コマンドで候補物質を作成します。
```
>> python3 ../bin/mkCandidatesFromMolecularBackbone.py molecularBackbone.prm
```

最大60の候補分子を molecularBackbone.smi に書き込みます。
確認用の画像ファイル(molecularBackbone_?.png/molecularBackbone_A?.png)も作成されます。

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

<img width="800" height="1000" alt="molecularBackbone_1" src="https://github.com/user-attachments/assets/509708ff-a406-425d-b383-cc57671c29a7" />

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

結果を表示するためのパラメータファイルを作成します。
```
>> sh ../bin/mkPrmForTdds2svg.sh > plot.prm
```
有効な吸収ピークを持つ物質を表示します。下記の例では、波長250nm 以上に、正規化した時のピーク高が0.1以上の吸収ピークがあるものを選別します。
```
>>python3 ../bin/makeTddftReport.py plot.prm -up 250 -r 0.1
```
uv_tddft-git_prepareCandidateA_plot_all.pptx が作られます。

<img width="1037" height="624" alt="image" src="https://github.com/user-attachments/assets/80580323-d7ce-4444-b6c4-02db7d46469a" />

オプションを付けなければ、全候補物資のスペクトルが表示されます。

各々のパラメータファイルはエディターで編集可能なので、必要な場合は、ディレクトリー tddft　内のマニュアルを参照に修正してください。
