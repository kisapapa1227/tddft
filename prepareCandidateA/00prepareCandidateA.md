### 概要

溶解度データセット(solubility.train.sdf)に含まれる分子から官能基を生成し、それらで骨格分子を修飾して候補物質を作成します。

計算コードは bin にあります。パスは環境に合わせて読み替えてください。
### 1. 官能基、及び、候補物質の自動作成

骨格分子、及び、修飾する部位を指定する設定ファイルを準備します。

例:
```
#ファイル molecularBackbone.prm:
furan-                # 出力ファイル内の候補物質名の接頭辞
c1occc1               # 骨格分子の SMILES 構造式
c1oc([*:1])cc1        # 官能基で修飾する部位を指定する SMARTS/SMILES 式
WeighFromTo:0,200     # （官能基の水素を除く）分子量の範囲
AtomNum:3,5           # （官能基の水素を除く）原子数の範囲
```
下記コマンドで候補物質を作成します。
```
python3 ../bin/mkCandidatesFromMolecularBackbone.py molecularBackbone.prm
```

- 最大60の候補分子を molecularBackbone.smi に書き出されます。
- 確認用の画像も生成されます(molecularBackbone_?.png/molecularBackbone_A?.png)。

例:
```
# ファイル molecularBackbone.smi
furan-000 'c1ccoc1'
furan-001 'CC(=O)Sc1ccco1'
furan-002 'CS(=O)c1ccco1'
furan-003 'NS(=O)(=O)c1ccco1'
furan-004 'ClC(Cl)c1ccco1’
...
```
- 作成した官能基や候補物質の一覧は生成されたPNGファイルで確認してください。

官能基一覧
<img width="1031" height="773" alt="image" src="https://github.com/user-attachments/assets/ef4560db-cc0a-4c3b-ba08-9943d268cba9" />


候補物質の一覧

<img width="800" height="1000" alt="molecularBackbone_1" src="https://github.com/user-attachments/assets/509708ff-a406-425d-b383-cc57671c29a7" />

### 2. 吸収スペクトルの計算

候補物質に対して一括で吸収スペクトル計算を行います。

2.1 候補物質一覧から実行スクリプトを作ります。

```
python3 ../bin/mkScript.py molecularBackbone.smi > runAll.sh
```

2.2 生成されたスクリプトを実行します。
```
sh runAll.sh
```

解析結果の可視化を行います。

2.3 解析結果の可視化用パラメータファイルを作成します。
```
sh ../bin/mkPrmForTdds2svg.sh > plot.prm
```
2.4 有効な吸収ピークを持つ物質のみを抽出してレポートを作成します。以下は、波長 250 nm 以上で正規化時のピーク高さが 0.1 以上のものを選別する例です。
```
python3 ../bin/makeTddftReport.py plot.prm -up 250 -r 0.1
```
- 上記コマンドの実行で uv_tddft-git_prepareCandidateA_plot_all.pptx が作成されます。
 
（出力例イメージ）

<img width="1037" height="624" alt="image" src="https://github.com/user-attachments/assets/80580323-d7ce-4444-b6c4-02db7d46469a" />

- オプションを付けなければ、全候補物資のスペクトルが表示されます。

### 補足
- 各パラメータファイル（例: molecularBackbone.prm, plot.prm）はテキストエディタで編集が可能です。必要に応じてパラメータを調整してください。
- 詳細は tddft ディレクトリ内のマニュアルを参照してください。
- 実行コマンドやパスは環境によって異なるため、スクリプト実行前に bin/ への参照や Python 環境を確認してください。
