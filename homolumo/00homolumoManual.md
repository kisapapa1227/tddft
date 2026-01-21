## 概要
分子軌道計算を行い、電子密度分布を表示します。

### 1. 計算手順
計算コードは bin にあるので、相対位置は読み替えてください。

計算は下記コマンドで行います。
```
>>python3 ../bin/getHomoLumo.py aniline ‘c1ccccc1N’
```
第一引数：物質名（出力ファイル名）、第二引数：対象物質のSMILES構造式

解析結果は物質名_ext_mf.pkl に保存されます。また、計算で使用したパラメータは物質名.erg に保存されます。

デフォルトで連続溶媒モデル（Continuum Solvent Model）のPCM（Polarizable Continuum Model）を組み込み、chloroform 相当の効果（4.7113）にしています。 

例えば、methanol 相当にしたい場合は 
```
>>python3 ../bin/getHomoLumo.py aniline ‘c1ccccc1N’　-eps 32.613 
```
と、オプションで指定してください。 

あるいは、溶媒モデルを利用したくない場合は、 
```
>>python3 ../bin/getHomoLumo.py aniline ‘c1ccccc1N’　-wo 
```
-wo をオプションを着けてください。 

### オプション
```
-mc max_cycle :収束計算の最大反復回数を指定する。
-bs 6311g :基底関数を指定する。
-xc CAMB3LYP :汎関数を指定する。 
-eps 4.7113 :　溶媒モデルで利用する誘電率を指定する。
-wo　:溶媒モデルを利用せずに計算する。
```

### 2. 結果の可視化

jupyterlab をインストール(pip install jupyterlab)するか、 colabを利用してください。

### 2.1 jupyterlab を利用する場合
1. jupyterlab を起動して(>>jupyter lab)、任意のブラウザーでhttp://127.0.0.1:888/labに接続する。
2. File->Open from Pathで、bin/showHomoLumo.ipynb を読み込む。
3. コードを実行する(パッケージのインストールとpythonのコードの二箇所)。
4. ファイル選択用のウィンドーが現れるので、物質名_ext_mf.pkl　を読み込む。 

### 2.2 colab を利用する場合
1. https://colab.research.google.com/に接続する。
2. ファイル->ノートブックをアップロードで、 bin/showHomoLumo.ipynb を読み込む。
3. 物質名_ext_mf.pkl　をアップロードする。
4. python コードの冒頭付近のdatafile=“物質名_ext_mf.pkl” を書き換える。
5. コードを実行する(パッケージのインストールとpythonのコードの二箇所)。

### 画面の操作
表示したい分子軌道を選択した後、マウスの右ボタンで回転、シフト＋右ボタンで拡大縮小をして、スクリーンショットで画像を保存する。
<img width="1038" height="855" alt="image" src="https://github.com/user-attachments/assets/ac291820-04b5-4524-bad7-afb37374bf94" />
