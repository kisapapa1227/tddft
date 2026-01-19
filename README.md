このプロジェクトでは、オープンソースの量子化学計算ツール「PySCF(Python Simulations of Chemistry Framework)」を利用して、物質の吸収スペクトルを求めるための python の解析コード、及び、前処理として候補物質の作成、後処理として解析結果をグラフ化、レポート作成の手順を提供します。

手順は下記の項目に分類されます。各ディレクトリーにサンプルコード、解析例、マニュアルがあります。、

1. 量子化学計算、及び、可視化、レポートの作成
    1) 吸収スペクトルを計算する(tddft)。
    2) 分子軌道 (HOMO/LUMO)を計算する(homoLumo)。
    3) 励起状態分子内プロトン移動 ESIPT (Excited state intramolecular proton transfer)を模して、吸収スペクトルを計算する(esipt)。 

1. 候補物質の作成（前処理）
    1) 修飾する官能基を探す。データベースより官能基のリストを作成し、それらを SMILES で指定した分子骨格に修飾した候補物質を作成する(pre1)。
    2) 官能基の組み合わせを探す。複数の有効な官能基を、骨格となる分子の複数箇所に修飾した時の候補物質の吸収スペクトルを比較する(pre2)。

```
conda create -n tddft-git python=3.12
pip install numpy matplot pyscf geometric python-pptx rdkit Chart Sequence scikit-learn xgboost
apt install imageMagick
```

python3.10 以降を利用する場合は、collections を利用しているパッケージに関して、下記を参考にコードの書き換えを行ってください。
```
in envs/tddft-git/lib/python3.12/site-packages/pptx/chart/chart.py series.py
from collection import Sequence -> from collection.abc import Sequence
'''


