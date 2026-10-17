# OklchSV — S/P 比（v0.1）

[English](README.md)

薄明視の表示モデル OklchSV で使う、暗所視/明所視の輝度比（S/P）の計算部分です。
論文 arXiv:2609.24819（CC BY 4.0、<https://arxiv.org/abs/2609.24819>）に対応します。

含むもの:

- **閉形式 S/P**（線形 sRGB から。論文で公開した 6 個の重み）
- **分光 S/P**（放射輝度から。V(λ)・V′(λ) の積分、白で正規化）
- 閉形式の前提になる **D65 への色順応**（線形 sRGB ⇄ XYZ、Bradford）

**変換本体（OklchSV の薄明視変換）は含みません。**

閉形式は「D65 に順応した**線形** sRGB」で評価してください（論文の条件 C2）。
普通の画像（ガンマ付きの sRGB）は `srgb_to_linear` で線形に戻してから渡します。
等色関数・視感効率関数の表は同梱していません。入手先は [data/README.md](data/README.md) にあります。

使い方・再現スクリプト・期待値は [README.md](README.md)（英語）を見てください。

解説記事（Zenn）: <https://zenn.dev/nekotrack>

ライセンス: MIT
