## 概要

SMILES記法で記述した物質の吸収スペクトルを化学量子計算で求めます。ディレクトリーにサンプルコードが準備してあります、

### 1. 計算手順

計算コードは bin にあります。相対位置は読み替えてください。

計算は下記コマンドで行います。
```
>>python3 ../bin/tddftSolver.py Indigo_pcm ‘C1=CC=C2C(=C1)C(=C(N2)C3=NC4=CC=CC=C4C3=O)O’
```
第一引数：物質名（出力ファイル名）、第二引数：対象物質のSMILES構造式

解析結果は、物質名.tdd に保存されます。

デフォルトで連続溶媒モデル（Continuum Solvent Model）のPCM（Polarizable Continuum Model）を組み込み、chloroform 相当の効果（4.7113）にしています。 

例えば、methanol 相当にしたい場合は 
```
>>python3 ../bin/tddftSolver.py Indigo 'C1=CC=C2C(=C1)C(=C(N2)C3=NC4=CC=CC=C4C3=O)O' -eps 32.613
```
と、オプションで指定してください。 

あるいは、溶媒モデルを利用したくない場合は、 
```
>>python3 ../bin/tddftSolver.py Indigo 'C1=CC=C2C(=C1)C(=C(N2)C3=NC4=CC=CC=C4C3=O)O' -wo 
```
-wo をオプションを着けてください。 

###  オプションの説明
```
>>python3 ../bin/tddftSolver.py Indigo_pcm ‘C1=CC=C2C(=C1)C(=C(N2)C3=NC4=CC=CC=C4C3=O)O’  options
-mc max_cycle :収束計算の最大反復回数を指定する。
-bs 6311g :基底関数を指定する。
-xc CAMB3LYP :汎関数を指定する。 
-eps 4.7113 :　溶媒モデルで利用する誘電率を指定する。
-wo　:溶媒モデルを利用せずに計算する。
```

### 2. 解析結果の可視化

入力データとして、上の計算手順で作成した tddファイルを利用します。

### 2.1 パワーポイントファイルで作成する。

下記コマンドで設定ファイルを作成します。
```
>>python3 ../bin/preparePlotParameterFile.py > plot.prm
```
カレントディレクトリーに存在する tdd ファイルが対象になります。
```
>> cat plot.prm
files:Indigo_pcm,Indigo
output:all
ops:[100,700,601],0.60,4000
```
1行目が処理をする tdd ファイルのリスト、二行目が出力ファイルの末尾名(uv_tddft-git_tddft_***.pptx)、三行目がグラフ化のパラメータになります。
この例では横軸が100-700nmで601等分(=1nm間隔)、縦軸の吸収強度が0-0.60 a.u.、共振強度から吸収スペクトルを作成する時の分散に関連するパラメータを4000 に設定します。分散に関連するパラメータが大きいと、ピーク分布の半値幅が小さくなります。
必要に応じて、このファイルを書き換えてください。

次のコマンドでパワーポイントのレポートが作成されます。
```
>> python3 ../bin/makeTddftReport.py plot.prm
```

<img width="1130" height="626" alt="image" src="https://github.com/user-attachments/assets/44e127cd-d533-4b3b-ab9d-b4c3f71a4631" />

### オプション
```
-n:吸収スペクトルを正規化（最大値で無次元化）します。
-N:吸収スペクトルを正規化（最大値で無次元化）し、最大値を凡例に書きます。
-s:uv域（220-300nm）の最大値で正規化します。
-up 280: 閾値 threshold (280nm) 以上に吸収ピークの存在するの分布のみ描画します。
-r 0.1: 上の吸収ピークで、threshold 以下の
-num 5: 一グラフ（一ページ）に描画する級数スペクトルの数を決めます。 
-no_table: 表を描きません。
```

### 2.2 エクセルファイルで作成する。
エクセルファイルでレポートを作成する場合は、下記コマンドを利用します。
```
>> python3 ../bin/makeTddft2Excel.py –d .
```
### オプション

```
-d input_directory: データ(tdd file)のあるディレクトリーを指定する。
-n tag_name:データを書き込むタグの名前を指定します。指定されない場合は、データのディレクトリー名になります。
-k key: key で指定した文字列をファイル名に含む tdd ファイルのみ処理します。
-range x0,x1,num: 横軸の最小値、最大値、分割数を指定します。デフォルトは100,600,501 です。
-stdev stdev:共振強度から吸収スペクトルを作成する時の分散に関連するパラメータを指定します。デフォルト値は4000です。 
-o outputFileName:出力するファイル名を指定します。デフォルトは、uv_all.xlsxです。
```
-o オプションで、存在する出力ファイルを指定すると、追記します。

<img width="1960" height="787" alt="image" src="https://github.com/user-attachments/assets/fc7c0fa5-18f7-4f93-902f-b730649be558" />

### 3. 補足説明
### 3.1 tddftSolver.py の概略
<img width="1139" height="639" alt="image" src="https://github.com/user-attachments/assets/1df647c2-a794-43de-9de6-8c09547dbe1e" />

<img width="1141" height="641" alt="image" src="https://github.com/user-attachments/assets/7cdf5cbf-1d6c-4ab2-b75f-bf519de9280c" />

### 3.2 基底関数の比較
<img width="1149" height="638" alt="image" src="https://github.com/user-attachments/assets/30843572-56af-4246-8610-cb9f4ca58943" />

### 3.3 汎関数の比較
<img width="1138" height="639" alt="image" src="https://github.com/user-attachments/assets/0449ae52-d05d-4e32-8e5c-3c3defb8d4f4" />
<img width="1144" height="421" alt="image" src="https://github.com/user-attachments/assets/ad5c77ce-a87a-4e64-9701-ef3a704fbd55" />

### 3.4 nstatの比較
<img width="1137" height="636" alt="image" src="https://github.com/user-attachments/assets/1b057e2b-d9f3-4a9d-90a7-5cb823f40ea0" />
<img width="1127" height="635" alt="image" src="https://github.com/user-attachments/assets/ae0982bf-124c-4ac0-8b66-2d3b708282c0" />

### 3.4 stddevの比較
<img width="1160" height="641" alt="image" src="https://github.com/user-attachments/assets/d6f4d439-ca5b-4e02-a51a-472376b10000" />

