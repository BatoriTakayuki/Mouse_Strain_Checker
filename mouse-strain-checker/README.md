# Mouse Strain Checker

マウスの公開RNA-seqのSRX番号から、**RNAで観測できた既知SNPに適合するマウス系統候補**を、染色体別・区間別に推定する独立ツールです。KO_Checkerへの変更・追加はありません。

**研究用プロトタイプ v0.1.0。合成データとモックによる31件のテストが合格。外部ツールを使う追加1件は依存実行環境がないためスキップ。実SRX、HISAT2/SAMtoolsの実行、Dockerビルド、ブラウザUIの実動作はこの作成環境では未検証です。判定閾値は精度評価で調整する必要があります。**

## 方針は妥当か

既知SNPを利用する方法は妥当です。ただし「染色体全体がどの系統由来かをRNA-seqだけで確定する」という目的には制約があります。

| 問題 | 実装した対処 | 残る制約 |
|---|---|---|
| 未発現領域ではSNPを観測できない | 無観測区間を明示し、補間しない | RNAが読めた領域の候補しか評価できない |
| 染色体内で由来が切り替わる | 10 Mb区間と染色体集計を併記。区間の候補が食い違えば染色体の確定的ラベルを出さない | HMM・相推定・組換え位置の精密推定は未実装。小さい導入断片は見逃しうる |
| 1遺伝子の高発現や近接SNPに偏る | パネルSNPを既定1 kb間隔で間引き、1座位の有効深度を20に制限。染色体ラベルには原則2つ以上の判定区間を要求 | 遺伝子内・連鎖による相関は残る。独立な遺伝子座の数とは同じでない |
| C57BL/6J参照へのマッピングバイアス | 結果に制約を明示。HISAT2 SNPグラフインデックスの作成オプションを用意 | 自動setupは線形参照。バイアスを自動除去した実装ではない |
| アレル特異的発現・インプリンティング | ベータ二項の過分散モデル、2系統モデルを併用 | 片アレル発現が強いとF1でも単一系統に見える可能性がある |
| C57BL/6Jと6NJ、129亜系統などが近い | 候補と次点を区別できるSNP数・スコア差を要求。不足時は判定不能 | 参照にない亜系統を同定できない。「N系」全体の同定でもない |
| エラー、マルチマップ、ペアの重複計数 | MAPQ、BQ、NH、SAMフラグで除外し、同じRG/QNAMEのペアは1回に集計 | PCR重複は既にduplicateフラグが付いたもののみ除く。UMI処理は未実装 |
| RNA編集・indel近傍・誤アラインメント | 既知のPASS SNPのみ利用し、第3塩基が多い座位を除外 | 既知A/G SNPに重なるRNA編集は除ききれない。ブラックリスト・indel近傍除外は未実装 |
| 参照にない系統でも何かが1位になる | 適合不足・情報不足を別の状態として出力 | 候補パネル外の系統を必ず検出できる保証はない |

**使いどころ:** 公共RNA-seqの背景系統の確認、B6/129混合の疑いの探索、KO近傍の残存ドナー領域候補の探索。DNA上の系統を確定する場合は、DNA SNPタイピングやWGS等と照合してください。混在、F1、汚染はこの出力だけでは区別できません。

## Dockerで起動

このフォルダで実行します。Docker DesktopまたはDocker Engine + Composeが必要です。

```bash
docker compose build
docker compose run --rm strain-checker setup --out /data/resources --threads 4
docker compose up -d
```

ブラウザで **http://localhost:8502** を開き、SRX番号を入力して「解析開始」を押します。参照データは `/data/resources`、結果は `/data/results` が初期値です。複数SRXは別々に解析します。同じSRXに属する複数SRRは、ENA上で同じ生物試料に属することを確認して統合します。

初回のsetupはGRCm39のFASTA、Mouse Genomes Project v8のSNP VCFを取得し、パネルとHISAT2インデックスを作成します。**大容量のダウンロード・ディスク・計算時間が必要です。** 実行時間とメモリ量は本環境では計測していません。FASTQ、BAM、SRA変換用一時ファイルも保存するため、参照ファイルだけの容量で見積もらないでください。

停止は `docker compose down`。解析中断は画面のボタン、またはコンテナ停止で行います。参照準備はCLIで実施します。UIは1つのタブで使用してください。ブラウザ再接続後のジョブ再接続・永続的なジョブ管理は未実装です。

### 既定の候補系統

MGP VCF内の `129S1_SvImJ, C57BL_6NJ, BALB_cJ, DBA_2J, FVB_NJ, CAST_EiJ, PWK_PhJ, WSB_EiJ` と、参照配列から追加する `C57BL_6J` の9候補です。その他を対象にする場合、初回に指定します。

```bash
docker compose run --rm strain-checker setup --out /data/resources_custom --threads 4 --strains 129S1_SvImJ,129P2_OlaHsd,129S5SvEvBrd,C57BL_6NJ,BALB_cJ,DBA_2J
```

GUIの参照データ欄を `/data/resources_custom` に変更してください。候補系統はVCFのsample名と厳密に一致させます。任意の「129」や「BALB/c」などの総称への変換は行いません。

`C57BL_6J` は「参照FASTA由来なので選択SNPでREF」という仮定で追加します。VCFで測定された実試料のGTではありません。欠損GTをREF扱いする処理はありません。

## CLIで実行

```bash
docker compose run --rm strain-checker run --resources /data/resources --srx SRX123456 --out /data/results --threads 4
```

`SRX123456` は入力例です。解析対象の実際の番号に置き換えてください。ERX/DRXにも対応します。対応対象はENAメタデータでマウス・Illumina・RNA-Seq・TRANSCRIPTOMICに該当するデータです。10x等のバーコード・UMI付きライブラリ用の処理はありません。通常のbulk RNA-seqを想定しています。

FASTQ取得はENAのHTTPSとMD5検証を優先します。FASTQ情報が不足している場合は `prefetch` + `fasterq-dump` を使います。ENAの通信エラーが続く場合には `--download-method sra` を指定してください。ペアエンドの孤立リードは解析から除外し、ログに記録します。

マッピングは `hisat2 --dta -k 5` とSAMtoolsによる座標ソートです。複数配置を探索した上で、系統推定に使うのは **NH=1かつMAPQ 30以上の配置**です。`-k 1` で多重配置を隠す運用はしません。

### 既存HISAT2 BAMを利用

```bash
docker compose run --rm strain-checker run --resources /data/resources --bam /data/existing/sample.bam --bam-assembly GRCm39 --out /data/results_existing
```

座標ソート済みBAMが必要です。NHタグのない読み取りは除外されます。BAMの染色体名・長さを照合し、assemblyの明示を要求しますが、同じ染色体長を持つ異なる参照の利用を完全に検出するものではありません。既存BAMが対応するFASTAでマッピングされたことを利用者側でも確認してください。

### 区間を細かくする

```bash
docker compose run --rm strain-checker run --resources /data/resources --srx SRX123456 --out /data/results_1Mb --window-size 1000000
```

区間を狭くすると残存ドナー断片を調べやすくなる一方、SNP不足の判定不能区間が増えます。閾値を下げてラベルを増やすことは精度保証にはなりません。

## 判定の中身

1. PASSのbiallelic SNPを選択。各候補のGTが高品質のホモ接合または半数体の座位のみ残す。FIがあれば全候補でFI=1、GQがあれば20以上を要求する。
2. 選んだ候補全てでGTが揃った共通座位を使う。欠損値・ヘテロ接合・低品質は座位ごと除外し、候補間の比較条件を揃える。
3. FASTAとVCFの座標系・REFを確認。候補間で差のないSNPを除外し、既定1 kb以上離して間引く。
4. 既定でBQ≥25、MAPQ≥30、NH=1。unmapped/secondary/supplementary/QC-fail/marked-duplicateを除外する。
5. 同じリードグループとQNAMEのペアが同じSNPを読む場合は1断片として数える。高品質塩基が食い違うペアはその座位で両方除く。
6. REF+ALTの断片数5以上、その他の塩基の比率10%以下の座位を使う。pileup深度上限に達した座位は除く。
7. 各系統のホモ接合モデルと常染色体上の2系統モデルを比較する。2系統の異なる座位での期待ALT比率は0.5。ベータ二項モデルで過分散を許し、1座位の有効深度を20に制限する。
8. 染色体は30 SNP、区間は10 SNP、次点と異なる診断SNPは10、座位平均スコア差は0.1以上を既定とする。適合損失が大きい場合も判定を保留する。染色体には原則2区間以上の局所的な支持が必要。

これらは**暫定的な工学的閾値**です。スコア差・適合損失はp値、事後確率、校正済みの確信度ではありません。スコアは高いほど適合します。候補を増やすと共通高品質座位が減り、近縁候補が増えると判定保留も増えます。

X、Y、MTは単一候補モデルのみで別扱いです。X不活化、性別、PAR、ミトコンドリアの共通ハプロタイプ、NUMT等を明示的にモデル化していません。性染色体・MTの候補ラベルから両親の由来を推定しないでください。

## 出力

各SRXのディレクトリ内に保存します。

| ファイル | 内容 |
|---|---|
| `report.html` | 区間バーと染色体別の表。ブラウザで開ける |
| `chromosomes.csv` | 染色体別の候補、次点、判定、SNP数、診断SNP数、スコア差、観測・判定区間数 |
| `windows.csv` | 固定幅区間別の同情報。無観測区間も含む。start/endは1-based inclusive |
| `candidate_scores.csv` | 各染色体の全候補モデルの順位と平均スコア |
| `sites.tsv.gz` | 観測されたパネル座位ごとのREF/ALT/その他の断片数、利用可否。全パネル座位は参照側のsites.bedにある |
| `evidence.json` | パネル情報、閾値、QC集計。NH欠損などの数はSNPに重なったentry数で、固有リード数ではない |
| `ena_runs.json` | 公開データのメタデータ・URL・チェックサム（SRX入力時） |
| `complete.json` | 解析完了マーカー、入力BAM、ツール版、参照情報 |
| `cache/` | FASTQ/SRA/BAM、マッピングログ、再利用用情報 |

| status | 読み方 |
|---|---|
| `candidate_supported` | 観測SNPでは単一候補が暫定基準を満たす |
| `mixed_or_heterozygous` | 2系統モデルに適合。F1/混在/汚染を区別する意味ではない |
| `regional_discordance` | 区間間で支持候補が異なる。染色体全体のcallは空欄 |
| `localized_evidence` | 支持が一部区間に限られる。染色体全体のcallは空欄 |
| `ambiguous` | 近縁候補等を区別する証拠が不足 |
| `poor_fit` | 候補モデルとの適合が悪い |
| `insufficient_sites` / `no_observed_snps` | 有効SNP不足 / 無観測 |

`candidate` は保留時も参考として出る単なる最上位候補です。**判定結果として扱うのはstatusとcallの組み合わせ**です。`complete.json` のないディレクトリは解析が完了していません。

同じ入力で再実行すると、MD5検証済みFASTQと、入力・参照・ツール版が一致する完成BAMを再利用します。系統スコアは再計算します。同じ出力先の並行実行には対応していません。初回setupの途中でパネルやインデックスの作成が失敗した場合、途中ディレクトリを別名に移すか、新しい保存先で実行してください。

## 自分の参照データ / mm10を利用

GRCm38またはGRCm39の**同じassemblyのFASTAと多系統VCF**を用意します。FASTAは非圧縮です。mm10用VCFをmm39のBAMに適用することはできません。染色体名のchr接頭辞変換と座標変換は異なります。

```bash
python -m mouse_strain samples --vcf strains.vcf.gz
python -m mouse_strain panel --fasta reference.fa --vcf strains.vcf.gz --assembly GRCm38 --strains 129S1_SvImJ,C57BL_6NJ --add-reference C57BL_6J --out resources38/panel
python -m mouse_strain build-index --fasta reference.fa --assembly GRCm38 --out resources38/hisat2/genome --threads 4
python -m mouse_strain run --fasta reference.fa --panel resources38/panel --index resources38/hisat2/genome --srx SRX123456 --out results38
```

FIを持たないVCFは標準で拒否します。独立に品質を確認したVCFの場合のみ `panel --allow-no-fi` を付けてください。GQがVCF自体にない場合、GQ閾値は適用できません。ヘッダーにassemblyの記述がなく、REFも偶然一致する誤ったassemblyは検出しきれないため、配布元も確認してください。

既存インデックスを推測だけで登録する機能は設けていません。`build-index` がFASTAのSHA-256を含むsidecarを作り、解析時にパネルと照合します。SNPグラフを使う場合はHISAT2形式のSNPファイルを同じ参照座標で準備し、`build-index --snp genome.snp ...` を使ってください。グラフ作成のメモリ負荷とSNPファイルの適切さは別途確認が必要です。

## ローカル環境 / テスト

Python 3.10以上、HISAT2、SAMtools 1.13以上が必要です。SRA取得にはSRA Toolkitも必要です。解析本体はPython標準ライブラリのみです。GUIのみ追加パッケージを使います。

```bash
python -m pip install -r requirements-web.txt
python -m mouse_strain web
python -m unittest discover -s tests -v
```

Docker内での検証コマンドです。外部ツールが揃っていれば、追加のスモークテストで人工FASTQ→HISAT2→BAM→染色体候補まで実行します。実SRXや実試料の精度検証の代用ではありません。

```bash
docker compose run --rm --entrypoint python strain-checker -m unittest discover -s /app/tests -v
```

ローカル起動のGUIは http://localhost:8501 です。Dockerとローカルでは初期データ保存先が異なるので、GUIのパス欄を確認してください。

## 精度検証で次に行うこと

1. 同じ組織または類似条件の既知純系SRXで、対象候補を含む・含まないパネルを比較する。B6J/6NJ、複数129亜系統の取り違え率を測る。
2. 既知F1で2系統モデルの感度を評価する。異なる純系RNAの計算上の混合も別に作り、F1との区別ができない点を確認する。
3. 同じ試料をサブサンプリングし、有効SNP数・区間幅・判定保留率の関係を評価する。保存先を分け、学習に使わない試料で閾値を評価する。
4. 線形参照とHISAT2 SNPグラフで結果の変化を比較し、組織、発現、ASEによる偏りを確認する。
5. 由来が既知の交雑・コンジェニック試料をDNA遺伝子型と比較し、残存ドナー断片をどの大きさまで検出できるか検証する。

## 出典

- [HISAT2公式マニュアル](https://daehwankimlab.github.io/hisat2/manual/)
- [Mouse Genomes Project](https://www.sanger.ac.uk/data/mouse-genomes-project/)
- [SNPsplit genome preparation: MGP v8、FI、GRCm39](https://felixkrueger.github.io/SNPsplit/genome_prep/genome_preparation/)
- [ENA programmatic access](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access.html)
- [SAMtools mpileup](https://www.htslib.org/doc/samtools-mpileup.html)
- [NCBI prefetch / fasterq-dump](https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump)
- [Keane et al., Mouse genomic variation and its effect on phenotypes and gene regulation](https://pmc.ncbi.nlm.nih.gov/articles/PMC3276836/)
- [Munger et al., RNA-Seq Alignment to Individualized Genomes Improves Transcript Abundance Estimates in Multiparent Populations](https://academic.oup.com/genetics/article/198/1/59/6073398)
