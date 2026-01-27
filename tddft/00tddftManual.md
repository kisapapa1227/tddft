## 概要

SMILES記法で記述した物質の吸収スペクトルを量子化学計算で求めます。ディレクトリーにサンプルコードが準備してあります。

### 1. 計算手順

計算コードは bin にあります。パスは環境に合わせて読み替えてください。

計算は下記コマンドで行います。
```
python3 ../bin/tddftSolver.py Indigo_pcm ‘C1=CC=C2C(=C1)C(=C(N2)C3=NC4=CC=CC=C4C3=O)O’
```
- 第1引数：物質名（出力ファイル名）
- 第2引数：対象物質のSMILES構造式(クォートで囲む)

解析結果は、物質名.tdd に保存されます。

デフォルトで連続溶媒モデル（Continuum Solvent Model）のPCM（Polarizable Continuum Model）を組み込み、chloroform 相当の誘電率（ε=4.7113）にしています。 

methanol 相当の誘電率（ε=4.7113)に変更する例： 
```
python3 ../bin/tddftSolver.py Indigo 'C1=CC=C2C(=C1)C(=C(N2)C3=NC4=CC=CC=C4C3=O)O' -eps 32.613
```

溶媒モデルを利用しない場合は -wo オプションを付けます。 
```
python3 ../bin/tddftSolver.py Indigo 'C1=CC=C2C(=C1)C(=C(N2)C3=NC4=CC=CC=C4C3=O)O' -wo 
```

###  オプションの説明
```
python3 ../bin/tddftSolver.py Indigo_pcm ‘C1=CC=C2C(=C1)C(=C(N2)C3=NC4=CC=CC=C4C3=O)O’  options
-mc max_cycle :収束計算の最大反復回数を指定する。
-bs 6311g :基底関数を指定する。
-xc CAMB3LYP :汎関数を指定する。 
-eps 4.7113 :　溶媒モデルで利用する誘電率を指定する。
-wo　:溶媒モデルを利用せずに計算する。
```

### 2. 解析結果の可視化

入力データとして、上の計算手順で作成した .tddファイルを利用します。

### 2.1 パワーポイントレポートの作成
- プロット設定ファイルを生成
```
python3 ../bin/preparePlotParameterFile.py > plot.prm
```
カレントディレクトリーに存在する .tdd ファイルが対象になります。

- plot.prm の例
```
files:Indigo_pcm,Indigo
output:all
ops:[100,700,601],0.60,4000
```
- 1行目:処理をする .tdd ファイルのリスト
- 2行目:出力ファイルの末尾名(出力ファイル名は”ディレクトリー名_末尾名.pptx")
- 3行目:グラフ化のパラメータ(例の意味)
   - 横軸範囲:100-700nm、分割数 601(=1nm間隔)
   - 縦軸(吸収強度)の上限　0.60 a.u.
   - 共振強度から吸収スペクトルを作成する時の分散に関連するパラメータ : 4000 (デフォルト)。値によりピーク幅が変化します（補足説明 3.4 stddev の比較を参照）。

PowerPointを作成:
```
 python3 ../bin/makeTddftReport.py plot.prm
```

<img width="1130" height="626" alt="image" src="https://github.com/user-attachments/assets/44e127cd-d533-4b3b-ab9d-b4c3f71a4631" />

### オプション
```
-n:吸収スペクトルを正規化（最大値で無次元化）する
-N:正規化した上で、最大値を凡例に表示する
-s:uv域（220-300nm）の最大値で正規化する
-up 280: 閾値(例:280nm) 以上に吸収ピークの存在する分布のみ描画
-r 0.1: 上の条件で、正規化時のピーク高が閾値(例：10%) 以上の分布のみ描画
-num 5: グラフ（ページ）あたりに描画する吸収スペクトルの数を指定 
-no_table: 表を出力しない
```

### 2.2 Excel レポートの作成
実行例：
```
 python3 ../bin/makeTddft2Excel.py –d .
```
### オプション

```
-d input_directory: データ(.tdd ファイル)のあるディレクトリーを指定
-n tag_name:データを書き込むシート名を指定（デフォルトはデータディレクトリー名)
-k key: key をファイル名に含む .tdd ファイルのみ処理
-range x0,x1,num: 横軸の最小値、最大値、分割数を指定(デフォルトは100,600,501)
-stddev stddev:共振強度から吸収スペクトルを作成する時の分散に関連するパラメータを指定(デフォルトは4000) 
-o outputFileName:出力するファイル名を指定(デフォルトはuv_all.xlsx)。既存ファイルを指定した場合は追記します。
```

<img width="1960" height="787" alt="image" src="https://github.com/user-attachments/assets/fc7c0fa5-18f7-4f93-902f-b730649be558" />

### 3. 参考資料

### 3.1 tddftSolver.py の概略
<img width="1139" height="639" alt="image" src="https://github.com/user-attachments/assets/1df647c2-a794-43de-9de6-8c09547dbe1e" />

<img width="1141" height="641" alt="image" src="https://github.com/user-attachments/assets/7cdf5cbf-1d6c-4ab2-b75f-bf519de9280c" />

### 3.2 基底関数の比較
<img width="1142" height="626" alt="image" src="https://github.com/user-attachments/assets/32697224-638b-47f0-9762-ad1c5c50bcc2" />

### 3.3 汎関数の比較
<img width="1153" height="657" alt="image" src="https://github.com/user-attachments/assets/1eebd96b-28ef-4f2f-86ca-5f8a971369e8" />

### 3.4 nstatの比較
<img width="1137" height="636" alt="image" src="https://github.com/user-attachments/assets/1b057e2b-d9f3-4a9d-90a7-5cb823f40ea0" />
<img width="1127" height="635" alt="image" src="https://github.com/user-attachments/assets/ae0982bf-124c-4ac0-8b66-2d3b708282c0" />

### 3.4 stddevの比較
<img width="1160" height="641" alt="image" src="https://github.com/user-attachments/assets/d6f4d439-ca5b-4e02-a51a-472376b10000" />

### 補足・注意事
- ファイルパスやコマンド中の引用符（'...'）などは環境に合わせて正しく指定してください。
- plot.prm や他の設定ファイルは必要に応じて手動で編集できます。
- 計算条件や対象物質によっては数値計算でエラーが生じることがあります。その場合はライブラリのバージョンを変えるなどして試してください。

