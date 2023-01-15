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

OpenACC化されたファイルとして`main.f90`が用意されています。
