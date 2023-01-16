# 2022年度GPU講習会2におけるハンズオン資料
これは[2022年度GPU講習会2](https://www.cfca.nao.ac.jp/content/gpu_workshop_2022_2)の[流体力学計算のFortran+OpenACCによる実装例](https://qiita.com/takiwaki_tomoya/items/af059b7fb8877f0e3d79)の資料です。

## ログイン
vpnを繋ぎ、`ssh`で`g00.cfca.nao.ac.jp`にログインします。
    
    ssh -XY <username>@g00.cfca.nao.ac.jp

## ソースファイル
`/cfca-work/gpuws01/gpuhydro`にファイルが用意されているのでコピーします。
    
    cd /cfca-work/username
    cp -r /cfca-work/gpuws01/gpuhydro .
    cd gpuhydro/KHf90openacc

OpenACC化されたお手本のファイルとして`main_ori.f90`が用意されています。
また、GPU化されてない練習用フォートランのファイルは`main_pra.f90`です。
こちらのファイルを開いて、[流体力学計算のFortran+OpenACCによる実装例]を参考にOpenACC化してみましょう。

## コンパイル
まずはフォートランコンパイラを使えるようにするにするため、以下のコマンドを実行します。

    module load nvhpc
    
お手本のファイルをコンパイルしたいときは以下のようにしてください。
    
    make kh_ori.x
    
 練習のファイルをコンパイルしたいときは以下です。
    
    make kh_pra.x
    
## 実行
OpenACC化できたら、テストしてみましょう。まずはお手本コードの場合。
    
    sbatch sj_ori.sh
    
次に練習コードの場合。
    
    sbatch sj_pra.sh
    
