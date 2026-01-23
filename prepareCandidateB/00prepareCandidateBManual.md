### 概要
骨格分子に複数の側鎖を追加して、候補物質を作成します。

- 計算コードは bin にあります。パスは環境に合わせて読み替えてください。
- 分子骨格、修飾する官能基、修飾する場所を指定して、候補物質を作成します。

### 1. 候補物質の生成（組み合わせ）

分子骨格、修飾する部位、官能基を設定するファイルを準備します。修飾箇所と官能基の記述行数が同じ場合は、数学で言う、組み合わせで分子骨格に修飾します。
```
>>cat chalconeDerivative.prm
chalcone-derivative-          #出力ファイル内の候補物質の頭文字
c1ccccc1-C=CC(=O)-c2ccccc(O)2 # 骨格分子のSMARTS/SMILES構造式
c1c([*:1])c([*:2])c([*:3])c([*:4])c1-C=CC(=O)-c2ccccc(O)2 # 官能基で修飾する部位
H,N,S(=O)(=O)                 # :1に結合する官能基
N,S(=O)(=O)                   # :2に結合する官能基
N,S(=O)(=O)                   # :3に結合する官能基
N,S(=O)(=O)N,ON(=O)(=O)       # :4に結合する官能基
```

下記コマンドで候補物質を作成します。
```
python3 ../bin/mkCandidatesWithSpecifiedSideChains.py chalconeDerivative.prm
```
- 候補物質はchalconeDerivative.smi に書き出されます。
- 確認用の画像ファイルも作成されます。


例：
```
#ファイル chalconeDerivative.smi
chalcone-derivative-01010101 'O=C(C=Cc1ccccc1)c1ccccc1O'
chalcone-derivative-02020202 'Nc1ccc(C=CC(=O)c2ccccc2O)c(N)c1N’
...
```

修飾する官能基の一覧

<img width="1200" height="300" alt="chalconeDerivative_A1" src="https://github.com/user-attachments/assets/3f2504c1-3331-4819-8ea2-d5dc26760bc7" />

作成された候補物質の一覧

<img width="800" height="1000" alt="chalconeDerivative_1" src="https://github.com/user-attachments/assets/ec421ed1-39f0-4dc2-96ee-6407a7d106a2" />

### 2. 候補物質の生成（順列）

分子骨格、修飾する部位、官能基を設定するファイルを準備します。官能基の記述が一行の場合は、数学で言う、順列で分子骨格を修飾します。
```
# ファイル furanDerivative.prm
furan-derivative-                   # 出力ファイル内の候補物質の頭文字
c1occc1                             # 骨格分子のSMILES構造式
c1([*:4])oc([*:1])c([*:2])c1([*:3]) # 官能基で修飾する部位
H,N(=O)=O,OCCCC,C(=O)O              # 修飾する官能基
```

下記コマンドで候補物質を作成します。
```
python3 ../bin/mkCandidatesWithSpecifiedSideChains.py franDerivative.prm
# ファイル franDerivative.smi
furan-derivative-01 'c1ccoc1'
furan-derivative-02 'CCCCOc1c([N+](=O)[O-])coc1C(=O)O'
furan-derivative-03 'CCCCOc1occ([N+](=O)[O-])c1C(=O)O'
furan-derivative-04 'CCCCOc1coc(C(=O)O)c1[N+](=O)[O-]’
.......
```

修飾する官能基

<img width="1200" height="150" alt="furanDerivative_A1" src="https://github.com/user-attachments/assets/eb15f392-00ed-4bcb-b2b5-628531877f24" />

作成された候補物質の一覧

<img width="800" height="800" alt="furanDerivative_0" src="https://github.com/user-attachments/assets/35db6836-e89f-497c-8601-8dbb33ebf2d1" />

### 3. 候補物質の一括処理

手順は prepareCandidateA と同一です。

3.1 候補物質一覧から実行スクリプトを作ります。

ただし、この候補物質(furanDerivative.smi)の多くは numpy>1.26.4 で、コード内の optimizerでエラーが出ます。numpy==1.26.4で試してください。

```
python3 ../bin/mkScript.py furanDerivative.smi > runAll.sh
```
3.2 生成されたスクリプトを実行します。
```
sh runAll.sh
```


解析結果の可視化を行います。

3.3 解析結果の可視化用パラメータファイルを作成します。
```
sh ../bin/mkPrmForTdds2svg.sh > plot.prm
```

有効な吸収ピークを持つ物質のみを抽出してレポートを作成します。以下は、波長 250 nm 以上で正規化時のピーク高さが 0.1 以上のものを選別する例です。
```
python3 ../bin/makeTddftReport.py plot.prm -up 250 r 0.1 -num 6 -N
```
- 上記コマンドでuv_prepareCandidateB_plot_all.pptx が作られます。

(出力イメージ)
<img width="970" height="683" alt="image" src="https://github.com/user-attachments/assets/d5bc38ca-6eda-4ceb-a57a-72fe763cb178" />

-オプションを付けなければ、全候補物資のスペクトルが表示されます。

### 補足
- 各パラメータファイルはテキストエディタで編集が可能です。必要に応じてパラメータを調整してください。
- 詳細は tddft ディレクトリ内のマニュアルを参照してください。
- 実行コマンドやパスは環境によって異なるため、スクリプト実行前に bin/ への参照や Python 環境を確認してください。
