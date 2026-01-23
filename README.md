このプロジェクトでは、オープンソースの量子化学計算ツール「PySCF（Python Simulations of Chemistry Framework）」を用いて、物質の吸収スペクトルを求めるための Python 解析コードを提供します。あわせて、前処理（候補物質の作成）および後処理（解析結果の可視化、報告書作成）の手順も含まれます。

各手順は以下のカテゴリに分類され、各ディレクトリにサンプルコード、解析例、マニュアルを配置しています

主な処理

- 吸収スペクトルを計算する（tddft）。
- 分子軌道（HOMO／LUMO）を計算する（homolumo）。
- 励起状態での分子内プロトン移動（ESIPT: Excited State Intramolecular Proton Transfer）を模擬して吸収スペクトルを計算する（esipt）。

候補物質の作成（前処理）

- 官能基候補の収集：データベースから官能基のリストを作成し、SMILES で指定した分子骨格に修飾して候補物質を生成します（prepareCandidateA）。
- 官能基組合せの検討：複数の有効な官能基を骨格の複数箇所に修飾した場合の候補物質について吸収スペクトルを比較します（prepareCandidateB）。


動作確認環境（例）

- conda create -n tddft-git python=3.12
- pip install pyscf==3.11.0 numpy==2.4.0 matplotlib geometric python-pptx rdkit Chart Sequence scikit-learn xgboost
- apt install ImageMagick

Python 3.10 以降を使用する場合の注意
- 一部のパッケージが古い import 形式（collections の Sequence を直接 import する書き方）を使っているため、Python 3.10 以降ではコードの修正が必要になる場合があります。例として、python-pptx 内のファイルを以下のように修正してください。
  
修正例:

```
Python
#### 変更前
from collections import Sequence

#### 変更後
from collections.abc import Sequence
```
（該当ファイルの例） envs/tddft-git/lib/python3.12/site-packages/pptx/chart/chart.py または envs/tddft-git/lib/python3.12/site-packages/pptx/chart/series.py

補足（トラブルシューティング）

- 対象物質や計算条件によっては、NumPy の行列演算でエラーが出ることがあります。その場合、古いバージョンの NumPy（例：numpy 1.26.4 など）へ戻すことで解決することがあるため、バージョンを切り替えて試してください。

