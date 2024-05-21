# QUMASUNマニュアル：コンパイルおよび実行

## はじめに

QUMASUNは密度汎関数理論(DFT)に基づく電子状態計算コードです。波動関数を実空間グリッドデータとして表現し、固有状態を求めます。

また、電子状態の時間発展(TDDFT)も扱うことができます。

MPIによる並列化が施されており、空間領域分割およびk点並列に対応しています。

## 必要な環境

以下のコンパイラ及びライブラリが必要です。

- C++コンパイラ（C++17以降推奨）
- MPI
- BLAS
- LAPACK
- ScaLAPACK（optional）
- FFTW

実際に動作確認しているコンパイラは次のものです。

- Microsoft Visual Studio
- Intel OneAPI (icpx, Intel MKL, Intel MPI)
- NECSDK (nc++, ただしscalapackリンク時にnfortが必要)

以下のシステムで動作確認しています。

- Intel Xeon Gold, (x86-64, skylake, cascadelake, icelake), (LINUX)
- Intel Core i7, (x86-64), (Windows)
- AMD Ryzen, (x86-64), (Windows)
- NEC SX-Aurora TSUBASA, (Vector Engine, VE10), (LINUX)

上記以外のシステムでも、適切にmakefileを設定すれば動作すると思います。

## ファイル構成

配布ファイルを解凍してください。解凍されたディレクトリの中身は次のような構成になっています。
```
qumasunYYYYMMDD/
|---source/      #ソースファイル
    |---makefile      #Intel OneAPI用のmakefileサンプル
    |---makefile.nec  #NECSDK環境用のmakefileサンプル
|---input/       #計算条件を記した入力ファイル
|---PP/          #計算に必要な擬ポテンシャルファイル
|---job.sh       #Slurmのジョブ投入スクリプトサンプル
|---jobW16k222_nec.sh   #NEC環境でのジョブ投入スクリプトサンプル
|---jobW128k222_nec.sh   #NEC環境でのジョブ投入スクリプトサンプル
```
ディレクトリ構成は適宜変えることが出来ますが、実行時に`input`および`PP`以下のファイルが読み込めるようにしておいてください。実効が失敗し、メッセージの出ないときには、これらのファイルのパスが間違っていることが多いです。

## コンパイル方法

コンパイルするには、`source/`ディレクトリ以下に移動し`make`を行います。
```
cd ~/fuga/hoge/qumasunYYYYMMDD/source/
make clean
make install
```
Intel OneAPI環境の場合は上記でコンパイルできます。NECSDK環境の場合には、`makefile.nec`を`makefile`とリネームしてから上記を実行してください。また、これ以外の環境の場合は、適宜makefileを編集してください。

makefileの中で記述するコンパイルコマンドとしてはMPIのC++コマンドを使ってください。以下にIntel OneAPI環境の場合の例を示します。
```
#コンパイルコマンドの例(Intel OneAPIの例)
CXX = mpiicpx -DUSE_MPI -DUSE_SCALAPACK \\
   -std=c++17 -O3 -ipo -qmkl -xCORE-AVX512 -qopt-zmm-usage=high
```
ここで、`-D`オプションによってディレクティブ定義を加えています。以下に`-D`で指定可能なディレクティブの説明を列挙します。

- `-DUSE_MPI`: MPI並列を有効化します。殆どの場合にこの記載は必須です。
- `-DUSE_SCALAPACK`: ScaLAPACKを利用する場合に指定します。ScaLAPACKを利用できる場合は`-DUSE_SCALAPACK`を記載することで固有値問題の解法にScaLAPACKを利用します。そうでない場合はLAPACKを使って固有値問題を解きます。
- `-D_NEC`: NECSDK環境でコンパイルする場合に記載してください。他の環境では記載しないでください。
- `-DASLFFTW_COMPLEX_LEGACY`: NECSDK環境でコンパイルする場合に記載してください。他の環境では記載しないでください。
- `-DLAPACK_NO_HEADER`: LAPACKのヘッダーファイル`lapack.h`が提供されない環境ではこのオプションを記載してください。それでも上手くコンパイルできない場合は、環境に合わせて`source/wrap_lapack.h`および`source/wrap_scalapack.h`を適宜書き換えてください。

makeに成功すると、実行ファイルとして`qumasun.exe`が生成されます。続けて`make install`を実行することで、`qumasunYYYYMMDD/`直下にコピーされます。

## 実行方法

コンパイルが成功し、上記のディレクトリ構成に対して、`qumasunYYYYMMDD/qumasun.exe`が存在するものとします。コマンド上で`qumasunYYYYMMDD/`直下に移動してください。

### 空間領域分割による実行

実行は`mpirun`もしくは`mpiexec`で実行します。ただし、Slurmでの`srun`などジョブ実行環境の指定に従ってください。`mpirun`を利用する場合に空間分割領域並列を行う場合は、次のように指示します。
```
mpi -n $NP ./qumasun.exe (入力ファイル) -ddm $DX $DY $DZ
```
ここで、`$NP, $DX, $DY, $DZ`は
```
$NP : 全プロセス数
$DX : x方向の空間分割数
$DY : y方向の空間分割数
$DZ : z方向の空間分割数
```
を意味し、実際の実行時には数値を記入します。

具体例として、計算条件を記した`input/W16k222.txt`を読み込んで計算を実行する場合は次のようになります。
```
mpirun -n 16 ./qumasun.exe input/W16k222.txt -ddm 2 2 4 
```
この例では、16プロセスの並列計算を行いますが、`-ddm 2 2 4`に従って空間領域をx, y, z方向にそれぞれ$2\times 2 \times 4$分割した空間領域分割を行って計算します。

また、全プロセス数`$NP`と空間分割数`$DX, $DY, $DZ`の間には
```
$NP = $DX * $DY * $DZ
```
の関係が成立している必要があります。ただし、次のk点並列を行う場合はこれとは異なります。

### 空間領域分割とk点並列による実行

DFTではk点サンプリングと呼ばれる手法を利用することが度々あります。これは、周期性のある構造に対して、そのうちの小さな領域だけを切り取った系を計算をする代わりに、逆空間(k空間)において波数だけをずらしたミラーとなる系を複数種類計算することで、元の広い空間と同様の解を得る方法です。逆空間で波数だけのズレたミラーの系たちは電子密度を共有しますが、それ以外の部分では独立に計算することが可能です。これを並列化によって同時実行することを、ここではk点並列と呼びます。

例として入力ファイル`input/W16k222.txt`の中には、次のような項目が記載されています。
```
System.KpointSample   2 2 2
```
これは、空間のx, y, z方向に対応したそれぞれの逆空間で$2\times 2\times 2$倍のミラーの系を作って計算することを指示しています。$2\times 2\times 2=8$ですので`M=8`個のミラーの系が生成されます。それに合わせて、k点並列も最大8並列まで実行可能です。

また、k点並列は空間領域分割との併用が可能です。この場合の実行方法も上記と同じです（最後の`-kpoint`オプションは記述しなくても現状では自動で判断されます）。
```
mpi -n $NP ./qumasun.exe (入力ファイル) -ddm $DX $DY $DZ -kpoint
```
ただし、k点並列数を`$NK`とすると、
```
$NP = $DX * $DY * $DZ * $NK
```
となるように、自動的に`$NK`が決められます。

例えば上記の`input/W16k222.txt`について、今度は次のように256プロセスを利用して実行することにします。
```
mpirun -n 128 ./qumasun.exe input/W16k222.txt -ddm 2 2 4 -kpoint
```
すると、`$NK=128/(2*2*4)=8`と設定されます。ここで`input/W16k222.txt`では`M=8`個のミラーが生成されますので、ちょうど`M=$NK`となっています。よって、全128プロセスは、16プロセスごとの8組に分割され、各組が各ミラー系を1つずつ担当します。各組では$2\times 2\times 4$空間領域分割として計算が実行されます。

ここで、k点並列数`$NK`はk点サンプル数(ミラーの系の数)`M`と一致している必要はなく、
```
M >= $NK >= 1
```
であれば大丈夫です。例えば、上記の計算を128プロセスではなく64プロセスで実行した場合は
```
mpirun -n 64 ./qumasun.exe input/W16k222.txt -ddm 2 2 4 -kpoint
```
k点並列数は`$NK=4`となります。この場合は16プロセスごとの4組に分割され、各組はミラーの系を2つずつ担当することになります。さらに約数である必要もなく、`$NK=3`の場合は、ミラーの系を3つ担当する組が２つ、ミラーの系を2つ担当する組が1つとなります。（ただし、ScaLAPACKおよびその内部のBLACSの実装次第ではこのようなアンバランスなk点並列はエラーとなる可能性があります）。このような理由から、計算効率が高くなるのは`$NK`がk点サンプル数`M`の約数となる場合です。

### 空間領域分割と状態分割による実行(TDDFTのみ)

空間領域分割だけでは、空間的な担当場所は違えど、全てのプロセスが全ての波動関数についての計算を行います。TDDFT計算（Ehrenfest MD含む）の場合には、担当する波動関数がプロセスごとに異なるように分割することができます。例えば10本の波動関数に対して、状態分割数2の場合、0番プロセスは波動関数0番から4番を、1番プロセスは波動関数5番から9番を担当します。これを状態分割と呼ぶことにします。

状態分割は、波動関数の内積（H行列やS行列などの作成）を計算する必要があるDFTのSCF計算では非常に遅くなるため実装していませんが、それらを必要としないTDDFT計算においては、効率の良い並列化が期待できます。

状態分割を行うには、実行時引数に`-sd 状態分割数`と指定します。以下に例を示します。
```
mpirun -n 128 ./qumasun.exe input/graphene_tddft.txt -ddm 2 2 2 -sd 16
```
これは、空間領域分割として$2\times 2\times 2$に分割し、かつ、状態分割として16個のプロセスグループに分割しています。併せて128プロセス必要です。このように、空間領域分割と併用できます。

また、次のようにすると、空間領域分割、状態分割、k点分割(この場合はspinのみ)を併用できます。
```
mpirun -n 128 ./qumasun.exe input/graphene_tddft.txt -ddm 2 2 2 -sd 8
```
ここでは状態分割数が8となっており、空間分割数と合わせても64プロセスとなります。よって、$(128/64=2)$分割がk点サンプリングに対して分割されます。TDDFTではk点サンプリングは未対応ですが、spin分極を有効にした計算においては、これはup状態とdown状態を別のプロセスで担当することに相当します。逆に言えば、spin分極を無効にした計算ではこのようなk点分割はエラーとなります。

## 実行結果

無事に計算が実行されると、標準出力の最後に次のような表示がみられます。
```
init:mu = -0.156688, -0.927815, 0.614440
final:mu = 0.517816, 0.517816, 0.517816
total rho = 192.000000
diff_abs_rho = 0.014845
Etot = -1109.431566  (test: 2668.030761, 395.323776, 395.323776, -891.498529)
Ekin = 288.875113
E_PP_local = -312.396941
  Eext_nucl_density = -312.396941 (non-use)
  Eext_nucl_point = -389.626374 (non-use)
E_PP_nonlocal = -12.324437
Ehart = 67.169687
Exc = -244.727654
  Ex = -222.458759
  Ec = -22.268895
  VxRho = -296.611679 (non-use)
  VcRho = -24.918978 (non-use)
Enn(Vlocal*rho_nucl-self) = -896.027334
  Enn(Vlocal*rho_nucl) = 390.794971
  Enn(self) = 1286.822305
  Enn(Vlocal) = 608.728008 (non-use)
  Enn(simpleCoulomb) = 2881.434993 (non-use)
correct_Coulomb_Vlocal = 0.000000
delta E_tot = -0.000280962373836
\int |rho-rho_prev| dr / Ne = 0.000077315157900
End SCF Step 100 ==============================

Note: definitions of energies
  Etot = Ekin + E_PP_local + Ehart + Exc + E_PP_nonlocal + Enn(Vlocal-self) + correct_Coulomb_Vlocal
  E_PP_local = \int{\rho V'_local}dV
  Eext_nucl_density = \int{\rho_nucl V_hart}dV  (non-use)
  Eext_nucl_point = \int{Q_n \delta(x-R_n) V_hart}dV  (non-use)
  E_PP_nonlocal = <psi|V_nonlocal|psi>
  Ehart = \int{\rho V_hart}dV
  Exc = Ex + Ec
  Enn(Vlocal-self) = Enn(Vlocal) - Enn(self)
  Enn(Vlocal) = \int{Q_n \delta(x-R_n) V'_local}dV
  Enn(self) = Q_n V_local(r=0),
    which is self interaction of nuclei estimated as an independent state
  Enn(direct) = QQ/r,
    where r is distance of two nuclei. (non-use)
  correct_Coulomb_Vlocal = QQ/r - Q*V_local(r),
    which is correction of two nuclei when r < core cutoff.
  V_local is on radial grid of an independent atom.
  \rho_nucl is corresponding charge of V_local,
    and is converted from radial grid to Cartesian grid.
  V'_local on Cartesian grid is the solution of the Poisson equation from rho_nucl,
    and is including the effect of periodic boundary.


Force =========================
Ftot:  0.00018   0.00002  -0.00088
(略)

Peak performance in GEMM=====================
  DGEMM1(x2HS)   : 8890.358817 GFLOPS  (71.588571 GFLOPS/root_proc)
  DGEMM2(x2x)    : 17835.539324 GFLOPS  (140.216640 GFLOPS/root_proc)
  ZGEMM3(HS2HS)  : 12807.690324 GFLOPS  (99.263422 GFLOPS/root_proc)

Calculation time====================
Total time      : 132.294099 [s]
mInitialState   : 5.183982 [s] / 1 [call]
mPrepareCore    : 0.841639 [s] / 1 [call]
InitialDensity  : 1.020854 [s] / 1 [call]

Calculation time: 125.247624 [s]
mSetOccupancy   : 3.928345 [s] / 101 [call]
mSetDensity     : 0.220120 [s] / 100 [call]
mSetPotential   : 6.578511 [s] / 101 [call]
mGetTotalEnergy : 7.448295 [s] / 100 [call]
EigenSolver     :--
--K-operation   : 16.340480 [s] / 39900 [call]
--V-operation   : 0.333452 [s] / 39900 [call]
--PPNonlocal    : 26.066257 [s] / 39900 [call]
--GEMM1(x2HS)   : 30.539525 [s] / 1668 [call]
--GEMM2(x2x)    : 8.338256 [s] / 1784 [call]
--GEMM3(HS2HS)  : 0.174727 [s] / 2168 [call]
--LAPACK        : 21.043631 [s] / 199 [call]
--MPI           : 3.549496 [s] / 796 [call]
--transpose     : 0.000000 [s] / 0 [call]
--others        : 0.682547 [s] / 798 [call]
Count_num 1.347456 [s] / 39900 [call]
Find_pp   0.000000 [s] / 0 [call]
Allocate  1.117222 [s] / 239400 [call]
Cut_psi   5.594827 [s] / 199500 [call]
Inner     4.397001 [s] / 199500 [call]
MPI       5.390747 [s] / 39900 [call]
Add local 4.171034 [s] / 199500 [call]
Paste     3.603377 [s] / 199500 [call]
```

このうち、`Etot`がDFTで求めた全エネルギーになります。
```
Etot = -8875.540690   (test: 109443.075174, 3171.256575, 3171.256575, -7123.321869)
```
この行の最初の数値`-8875.540690`が全エネルギーです。括弧内の数値は開発向けの数値ですので無視してください（将来は記載が消える予定です）。

この例では
```
\int |rho-rho_prev| dr / Ne = 0.000018471065783
```
の行で表示される値が十分に小さくない為、まだ収束に至る前に、設定したステップ数（この例では100）に達して計算が終了したことになります。

## ライセンス

Copyright (c) 2024, Atsushi M. Ito Released under the MIT license 

https://opensource.org/licenses/mit-license.php
