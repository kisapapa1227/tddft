## 概要
分子軌道計算を行い、電子密度分布(HOMO/LUMO等)を可視化します。

### 1. 計算手順
計算コードは bin にあるので、パスは環境に合わせて読み替えてください。

実行例:
```
python3 ../bin/getHomoLumo.py aniline ‘c1ccccc1N’
```
- 第1引数：物質名（出力ファイル名）
- 第2引数：対象物質のSMILES構造式（クォートで囲む）

出力:
- 解析結果本体:物質名_ext_mf.pkl
- 計算に使用したパラメータ:物質名.erg

デフォルトで連続溶媒モデル（Continuum Solvent Model）のPCM（Polarizable Continuum Model）を組み込み、chloroform 相当の誘電率（ε=4.7113）にしています。

methanol 相当の誘電率（ε=4.7113)に変更する例：
```
python3 ../bin/getHomoLumo.py aniline ‘c1ccccc1N’ -eps 32.613 
```

溶媒モデルを利用しない場合は -wo オプションを付けます。
```
python3 ../bin/getHomoLumo.py aniline ‘c1ccccc1N’ -wo 
```

### オプション
```
-mc max_cycle :収束計算の最大反復回数を指定する。
-bs 6311g :基底関数を指定する。
-xc CAMB3LYP :汎関数を指定する。 
-eps 4.7113 :　溶媒モデルで利用する誘電率を指定する。
-wo　:溶媒モデルを利用せずに計算する。
```

### 2. 結果の可視化

可視化は JupyterLab（ローカル）または Google Colab（クラウド）のどちらかで行えます。事前に必要パッケージをインストールしてください（例：pip install jupyterlab pyscf numpy rdkit 等）。

### 2.1 jupyterlab を利用する場合
- jupyterlab を起動
```
jupyter lab
```
- ブラウザーでhttp://127.0.0.1:888/labにアクセスする。
- メニューから File → Open From Path を選び、bin/showHomoLumo.ipynb を読み込みます。
- ノートブック中のセルを順に実行します（パッケージのインストールセルと解析用コードセルの2箇所）。
- ファイル選択ダイアログが出たら、物質名_ext_mf.pkl を選択して読み込みます。

### 2.2 colab を利用する場合
- ブラウザーでhttps://colab.research.google.com/にアクセスする。
- メニューから「ファイル」→「ノートブックをアップロード」で bin/showHomoLumo.ipynb をアップロードします。
- 物質名_ext_mf.pkl を Colab にアップロードします（左サイドバーのファイルタブからアップロードできます）。
- ノートブック中の冒頭付近にある datafile = "物質名_ext_mf.pkl" の行を書き換え、アップロードしたファイル名に合わせます。
- ノートブックのセルを順に実行します（パッケージのインストールセルと解析用コードセルの2箇所）。

### 画面の操作(可視化)
- 表示したい分子軌道を選択します。
- マウス操作：
   - 右ボタンドラッグ：回転
   - Shift + 右ボタンドラッグ：拡大／縮小
- スクリーンショットを使って画像を保存してください。
<img width="1038" height="855" alt="image" src="https://github.com/user-attachments/assets/ac291820-04b5-4524-bad7-afb37374bf94" />
