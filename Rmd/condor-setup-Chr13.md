Setting up for use of Condor CHTC with Chr 13 hotspot
================
Frederick Boehm
06/08/2019

## Goals

Goal here is to document the setup & organization of files that I need
for using Condor CHTC to do my two-dimensional scans as part of the
analyses of the four Keller hotspots. Note that they reported 5
hotspots.

## Setting

Recall that I have downloaded a Keller 2018 .RData file from DataDryad
that contains many of their results. The .RData file,
“data/Attie\_DO378\_eQTL\_viewer\_v1.Rdata” annotates 139 traits as
belonging to the chromosome 2 hotspot.

Following Karl’s suggestion, I want to look at not only the local trait
that Keller identifies as causal for the hotspot, but also many other
local-ish traits that have reasonable univariate LOD
scores.

## Preparing RDS files for Chr 13 analysis

``` r
load("~/Box Sync/attie/keller2018-chr2-hotspot-chtc/data-to-ignore/Attie_DO378_eQTL_viewer_v1.Rdata")
ls()
```

    ## [1] "dataset.clinical.phenotypes" "dataset.islet.hotspots"     
    ## [3] "dataset.islet.modules"       "dataset.islet.rnaseq"       
    ## [5] "ensembl.version"             "genoprobs"                  
    ## [7] "K"                           "map"                        
    ## [9] "markers"

``` r
library(tidyverse)
```

    ## ── Attaching packages ──────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.1       ✔ purrr   0.3.2  
    ## ✔ tibble  2.1.1       ✔ dplyr   0.8.0.1
    ## ✔ tidyr   0.8.3       ✔ stringr 1.4.0  
    ## ✔ readr   1.3.1       ✔ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
dataset.islet.rnaseq$expr %>%
  dim
```

    ## [1]   378 21771

``` r
colnames(dataset.islet.rnaseq$expr) %>%
  head
```

    ## [1] "ENSMUSG00000000001" "ENSMUSG00000000028" "ENSMUSG00000000037"
    ## [4] "ENSMUSG00000000049" "ENSMUSG00000000056" "ENSMUSG00000000058"

We need to get the names of the nonlocal traits that map to Chromosome
13 hotspot.

``` r
library(qtl2)
```

    ## 
    ## Attaching package: 'qtl2'

    ## The following object is masked from 'package:readr':
    ## 
    ##     read_csv

We want the ensembl ids for the traits that map to Chr 13 hotspot.
First, define the boundaries of the hotspot. The supplemental table in
Keller et al. 2018 lists only the midpoint of the hotspot. I believe
that the half-width is 2 Mb. Yes, this is indeed the half-width per the
Rmd file that Dan Gatti prepared.

``` r
hot_chrom <- 13
hot_center <- 112.5# for chr13
hot_lower <- hot_center - 2
hot_upper <- hot_center + 2
```

Now, filter to get those traits that map to the hotspot. Note that we’re
likely to get local traits, too, in this step.

``` r
dataset.islet.rnaseq$lod.peaks %>%
  filter(chrom == hot_chrom, pos <= hot_upper, pos >= hot_lower) %>%
  left_join(y = dataset.islet.rnaseq$annots, by = c("annot.id" = "gene_id")) %>%
  filter(!duplicated(annot.id)) %>%
  filter(lod >= 7.18)
```

    ##               annot.id    marker.id chrom      pos        lod
    ## 1   ENSMUSG00000038552 13_112640753    13 112.6408   7.652434
    ## 2   ENSMUSG00000039328 13_111608670    13 111.6087   8.752663
    ## 3   ENSMUSG00000039697 13_111832170    13 111.8322   9.195440
    ## 4   ENSMUSG00000040043 13_111977823    13 111.9778   7.188884
    ## 5   ENSMUSG00000040151 13_111532266    13 111.5323   7.812328
    ## 6   ENSMUSG00000040451 13_111712277    13 111.7123   8.807171
    ## 7   ENSMUSG00000040481 13_111664743    13 111.6647   7.715262
    ## 8   ENSMUSG00000006390 13_112506522    13 112.5065   8.983694
    ## 9   ENSMUSG00000040722 13_111425219    13 111.4252   9.624795
    ## 10  ENSMUSG00000042348 13_114143712    13 114.1437  31.669348
    ## 11  ENSMUSG00000042379 13_113610203    13 113.6102  57.443471
    ## 12  ENSMUSG00000042417 13_112978421    13 112.9784  10.415082
    ## 13  ENSMUSG00000042426 13_112830925    13 112.8309  32.136340
    ## 14  ENSMUSG00000042607 13_110749910    13 110.7499   7.680622
    ## 15  ENSMUSG00000045761 13_112251959    13 112.2520  10.910632
    ## 16  ENSMUSG00000046688 13_112459835    13 112.4598  13.554697
    ## 17  ENSMUSG00000047789 13_112689915    13 112.6899  35.950422
    ## 18  ENSMUSG00000048329 13_112468416    13 112.4684   7.380910
    ## 19  ENSMUSG00000049191 13_112506522    13 112.5065   8.947704
    ## 20  ENSMUSG00000049985 13_112459835    13 112.4598  93.090172
    ## 21  ENSMUSG00000052310 13_111832170    13 111.8322  10.402367
    ## 22  ENSMUSG00000054321 13_112509520    13 112.5095   7.553424
    ## 23  ENSMUSG00000010307 13_111951232    13 111.9512   7.478484
    ## 24  ENSMUSG00000058248 13_112472706    13 112.4727   7.590095
    ## 25  ENSMUSG00000059714 13_111712277    13 111.7123   9.959086
    ## 26  ENSMUSG00000059866 13_111357170    13 111.3572   8.238779
    ## 27  ENSMUSG00000061878 13_112530500    13 112.5305  10.054673
    ## 28  ENSMUSG00000061981 13_112205697    13 112.2057  11.109564
    ## 29  ENSMUSG00000066441 13_111357170    13 111.3572   8.583333
    ## 30  ENSMUSG00000066800 13_112216066    13 112.2161   8.865916
    ## 31  ENSMUSG00000070565 13_111951232    13 111.9512   9.968589
    ## 32  ENSMUSG00000015652 13_111425219    13 111.4252   9.186866
    ## 33  ENSMUSG00000071286 13_113510482    13 113.5105   7.227139
    ## 34  ENSMUSG00000071379 13_112527503    13 112.5275   7.846047
    ## 35  ENSMUSG00000071856 13_111608670    13 111.6087   7.587528
    ## 36  ENSMUSG00000072683 13_112506522    13 112.5065   7.354954
    ## 37  ENSMUSG00000074151 13_112211516    13 112.2115   8.892326
    ## 38  ENSMUSG00000016018 13_112833286    13 112.8333  13.792436
    ## 39  ENSMUSG00000016458 13_112216066    13 112.2161  10.886499
    ## 40  ENSMUSG00000016528 13_111608670    13 111.6087   7.308811
    ## 41  ENSMUSG00000081232 13_113102831    13 113.1028   8.172086
    ## 42  ENSMUSG00000083816 13_111425219    13 111.4252   8.375781
    ## 43  ENSMUSG00000085204 13_111832170    13 111.8322  41.336950
    ## 44  ENSMUSG00000085240 13_112506522    13 112.5065   8.197156
    ## 45  ENSMUSG00000086111 13_111900621    13 111.9006  39.076117
    ## 46  ENSMUSG00000017929 13_112205697    13 112.2057  10.026268
    ## 47  ENSMUSG00000018001 13_112828563    13 112.8286   7.747569
    ## 48  ENSMUSG00000018500 13_111712277    13 111.7123   9.120063
    ## 49  ENSMUSG00000091370 13_111772224    13 111.7722   9.818949
    ## 50  ENSMUSG00000018862 13_112506522    13 112.5065   8.850742
    ## 51  ENSMUSG00000018965 13_112509520    13 112.5095   9.377684
    ## 52  ENSMUSG00000097110 13_112870831    13 112.8708   7.331156
    ## 53  ENSMUSG00000097273 13_111883727    13 111.8837   7.658997
    ## 54  ENSMUSG00000097926 13_112211516    13 112.2115   7.390634
    ## 55  ENSMUSG00000098188 13_112509520    13 112.5095   8.459124
    ## 56  ENSMUSG00000020083 13_113258281    13 113.2583   7.664903
    ## 57  ENSMUSG00000020607 13_112459835    13 112.4598  11.175418
    ## 58  ENSMUSG00000020834 13_112665334    13 112.6653   7.993359
    ## 59  ENSMUSG00000021133 13_112468416    13 112.4684   9.464521
    ## 60  ENSMUSG00000021360 13_112828563    13 112.8286   9.250351
    ## 61  ENSMUSG00000000826 13_112530500    13 112.5305   7.917426
    ## 62  ENSMUSG00000001120 13_112205697    13 112.2057   7.949503
    ## 63  ENSMUSG00000001227 13_111712277    13 111.7123   8.873786
    ## 64  ENSMUSG00000021754 13_111772224    13 111.7722  17.416967
    ## 65  ENSMUSG00000021756 13_112468416    13 112.4684 122.115897
    ## 66  ENSMUSG00000021758 13_112640753    13 112.6408  24.192736
    ## 67  ENSMUSG00000021759 13_112828563    13 112.8286  13.141513
    ## 68  ENSMUSG00000021760 13_113240258    13 113.2403  14.186191
    ## 69  ENSMUSG00000021763 13_113249270    13 113.2493  20.332189
    ## 70  ENSMUSG00000021764 13_114143712    13 114.1437  20.031153
    ## 71  ENSMUSG00000022016 13_112472706    13 112.4727   8.354470
    ## 72  ENSMUSG00000022500 13_111532266    13 111.5323  10.288415
    ## 73  ENSMUSG00000001507 13_112640753    13 112.6408   8.621933
    ## 74  ENSMUSG00000024036 13_112455544    13 112.4555  13.755156
    ## 75  ENSMUSG00000024697 13_111368009    13 111.3680  10.815753
    ## 76  ENSMUSG00000024985 13_112464125    13 112.4641   9.433954
    ## 77  ENSMUSG00000025234 13_112530500    13 112.5305   9.123671
    ## 78  ENSMUSG00000025949 13_112530500    13 112.5305   9.416698
    ## 79  ENSMUSG00000026070 13_112472706    13 112.4727  13.345952
    ## 80  ENSMUSG00000026072 13_112530500    13 112.5305   7.330709
    ## 81  ENSMUSG00000001946 13_112530500    13 112.5305   7.814057
    ## 82  ENSMUSG00000002428 13_112205697    13 112.2057  10.156422
    ## 83  ENSMUSG00000026417 13_112764754    13 112.7648  10.474036
    ## 84  ENSMUSG00000026475 13_112640753    13 112.6408   7.428318
    ## 85  ENSMUSG00000026672 13_113249270    13 113.2493   8.805851
    ## 86  ENSMUSG00000026880 13_112509520    13 112.5095   8.162196
    ## 87  ENSMUSG00000026887 13_112506522    13 112.5065   9.988480
    ## 88  ENSMUSG00000026934 13_110609526    13 110.6095   7.667522
    ## 89  ENSMUSG00000027111 13_112472706    13 112.4727  13.565496
    ## 90  ENSMUSG00000027293 13_112205697    13 112.2057   9.447003
    ## 91  ENSMUSG00000027339 13_111608670    13 111.6087   9.735459
    ## 92  ENSMUSG00000027340 13_112530500    13 112.5305   8.090962
    ## 93  ENSMUSG00000027540 13_112828563    13 112.8286  10.387917
    ## 94  ENSMUSG00000028001 13_112472706    13 112.4727   7.950234
    ## 95  ENSMUSG00000002578 13_111608670    13 111.6087   8.236996
    ## 96  ENSMUSG00000028149 13_112205697    13 112.2057  19.085265
    ## 97  ENSMUSG00000029178 13_111772224    13 111.7722   7.720119
    ## 98  ENSMUSG00000029352 13_111532266    13 111.5323  11.477034
    ## 99  ENSMUSG00000029426 13_112506522    13 112.5065  10.884415
    ## 100 ENSMUSG00000029708 13_111608670    13 111.6087  10.023750
    ## 101 ENSMUSG00000030199 13_111532266    13 111.5323   7.421417
    ## 102 ENSMUSG00000030259 13_112689915    13 112.6899   7.311877
    ## 103 ENSMUSG00000030283 13_112472706    13 112.4727  11.073552
    ## 104 ENSMUSG00000030522 13_112281956    13 112.2820   8.478469
    ## 105 ENSMUSG00000030523 13_112281956    13 112.2820   7.490721
    ## 106 ENSMUSG00000031555 13_112472706    13 112.4727  11.728945
    ## 107 ENSMUSG00000004040 13_112828563    13 112.8286   7.502323
    ## 108 ENSMUSG00000032030 13_112472706    13 112.4727  11.071497
    ## 109 ENSMUSG00000032338 13_110537457    13 110.5375   7.242060
    ## 110 ENSMUSG00000032369 13_111951232    13 111.9512   8.496910
    ## 111 ENSMUSG00000032515 13_112530500    13 112.5305   7.606507
    ## 112 ENSMUSG00000032572 13_112216066    13 112.2161   9.203811
    ## 113 ENSMUSG00000032578 13_111883727    13 111.8837   8.708278
    ## 114 ENSMUSG00000032583 13_111368009    13 111.3680   7.297706
    ## 115 ENSMUSG00000032727 13_112764754    13 112.7648   8.502076
    ## 116 ENSMUSG00000033392 13_112509520    13 112.5095   9.738075
    ## 117 ENSMUSG00000033740 13_112506522    13 112.5065  11.480646
    ## 118 ENSMUSG00000034006 13_112640753    13 112.6408   7.968233
    ## 119 ENSMUSG00000034247 13_112216066    13 112.2161   9.841057
    ## 120 ENSMUSG00000034616 13_112459835    13 112.4598   7.426829
    ## 121 ENSMUSG00000035049 13_112665334    13 112.6653   9.304589
    ## 122 ENSMUSG00000035266 13_112530500    13 112.5305   7.339308
    ## 123 ENSMUSG00000035493 13_111832170    13 111.8322   7.817699
    ## 124 ENSMUSG00000037426 13_112281956    13 112.2820   8.086635
    ## 125 ENSMUSG00000037447 13_112066663    13 112.0667  10.059070
    ##            symbol chr      start        end strand     middle
    ## 1           Fndc4   5  31.292255  31.295877     -1  31.294066
    ## 2          Rnf122   8  31.111820  31.131482      1  31.121651
    ## 3           Ncoa7  10  30.645584  30.803107     -1  30.724346
    ## 4           Rbms2  10 128.129470 128.180297     -1 128.154884
    ## 5          Hs2st1   3 144.431103 144.570181     -1 144.500642
    ## 6           Sgms1  19  32.122727  32.389714     -1  32.256220
    ## 7            Bptf  11 107.033081 107.132127     -1 107.082604
    ## 8          Elovl1   4 118.428093 118.432953      1 118.430523
    ## 9          Scamp5   9  57.441328  57.468024     -1  57.454676
    ## 10          Arl15  13 113.794505 114.157461      1 113.975983
    ## 11           Esm1  13 113.209659 113.218098      1 113.213878
    ## 12           Ccno  13 112.987802 112.990778      1 112.989290
    ## 13          Dhx29  13 112.927730 112.969431      1 112.948580
    ## 14           Asb4   6   5.383386   5.433022      1   5.408204
    ## 15        Fam179a  17  71.673261  71.729669      1  71.701465
    ## 16           Tifa   3 127.789872 127.832135      1 127.811004
    ## 17        Slc38a9  13 112.660766 112.738743      1 112.699754
    ## 18         Mfsd6l  11  68.556186  68.558244      1  68.557215
    ## 19          Rgag4   X 102.066544 102.071304     -1 102.068924
    ## 20        Ankrd55  13 112.288451 112.384002      1 112.336226
    ## 21        Slc39a1   3  90.248172  90.253612      1  90.250892
    ## 22          Taf4b  18  14.783245  14.900359      1  14.841802
    ## 23        Tmem86a   7  47.050640  47.054776      1  47.052708
    ## 24          Kcnh1   1 192.190784 192.510159      1 192.350472
    ## 25          Flot1  17  35.823230  35.832791      1  35.828010
    ## 26          Tnip2   5  34.496087  34.513997     -1  34.505042
    ## 27          Sphk1  11 116.530925 116.536674      1 116.533800
    ## 28          Flot2  11  78.037931  78.060434      1  78.049182
    ## 29          Rdh11  12  79.174337  79.192293     -1  79.183315
    ## 30         Rnasel   1 153.749426 153.764221      1 153.756824
    ## 31         Rasal2   1 157.135183 157.412595     -1 157.273889
    ## 32         Steap1   5   5.736317   5.749326     -1   5.742822
    ## 33         Sowahc  10  59.221953  59.226431      1  59.224192
    ## 34         Hpcal1  12  17.690814  17.791926      1  17.741370
    ## 35            Mcc  18  44.425061  44.812182     -1  44.618622
    ## 36         Gm7457   6 142.814231 142.833494      1 142.823862
    ## 37          Nlrc5   8  94.472763  94.527272      1  94.500018
    ## 38        Skiv2l2  13 112.867780 112.927380     -1 112.897580
    ## 39            Wt1   2 105.126529 105.173616      1 105.150072
    ## 40       Mapkapk2   1 131.053704 131.097543     -1 131.075624
    ## 41        Gm14373   X   5.699727   5.701884     -1   5.700806
    ## 42        Gm13033   4 132.884545 132.886626      1 132.885586
    ## 43        Gm15327  13 111.809141 111.811039      1 111.810090
    ## 44        Gm12119  11  33.695894  33.697908      1  33.696901
    ## 45        Gm15326  13 111.867936 111.874730      1 111.871333
    ## 46        B4galt5   2 167.298444 167.349183     -1 167.323814
    ## 47          Cyth3   5 143.622447 143.710250      1 143.666348
    ## 48        Adora2b  11  62.248984  62.266453      1  62.257718
    ## 49  5730435O14Rik  10  41.568016  41.576476      1  41.572246
    ## 50          Otop3  11 115.334731 115.346927      1 115.340829
    ## 51          Ywhah   5  33.018816  33.027966      1  33.023391
    ## 52        Gm26726   X  73.402211  73.411436     -1  73.406824
    ## 53  C630004M23Rik   2  10.047838  10.049519      1  10.048678
    ## 54        Gm26575  18   5.572824   5.575055     -1   5.573940
    ## 55         Sowahc  10  59.221953  59.226434      1  59.224194
    ## 56  2010107G23Rik  10  62.107655  62.143920     -1  62.125788
    ## 57         Fam84a  12  14.147599  14.152038     -1  14.149818
    ## 58         Dhrs13  11  78.032280  78.037866      1  78.035073
    ## 59  4933426M11Rik  12  80.790532  80.880832      1  80.835682
    ## 60          Gcnt2  13  40.859768  40.960891      1  40.910330
    ## 61         Dnajc5   2 181.520485 181.555133      1 181.537809
    ## 62          Pcbp3  10  76.761857  76.961887     -1  76.861872
    ## 63         Sema6b  17  56.123085  56.140343     -1  56.131714
    ## 64         Map3k1  13 111.746428 111.808993     -1 111.777710
    ## 65          Il6st  13 112.464070 112.510086      1 112.487078
    ## 66           Ddx4  13 112.598333 112.652310     -1 112.625322
    ## 67         Ppap2a  13 112.800894 112.867881      1 112.834388
    ## 68           Gpx8  13 113.042763 113.046388     -1 113.044576
    ## 69       BC067074  13 113.317084 113.379711      1 113.348398
    ## 70         Ndufs4  13 114.287795 114.388094     -1 114.337944
    ## 71         Akap11  14  78.492246  78.536860     -1  78.514553
    ## 72          Litaf  16  10.959275  11.066157     -1  11.012716
    ## 73          Itga3  11  95.044474  95.076801     -1  95.060638
    ## 74        Slc37a1  17  31.295483  31.350696      1  31.323090
    ## 75          Gna14  19  16.435667  16.610818      1  16.523242
    ## 76         Tcf7l2  19  55.741810  55.933654      1  55.837732
    ## 77          Arih1   9  59.388258  59.486618     -1  59.437438
    ## 78        Pikfyve   1  65.186750  65.274012      1  65.230381
    ## 79         Il18r1   1  40.465552  40.500854      1  40.483203
    ## 80          Il1r1   1  40.225080  40.316198      1  40.270639
    ## 81           Esam   9  37.528078  37.538319      1  37.533198
    ## 82           Hltf   3  20.057811  20.118490      1  20.088150
    ## 83           Pigr   1 130.826684 130.852249      1 130.839466
    ## 84          Rgs16   1 153.740349 153.745468      1 153.742908
    ## 85           Optn   2   5.020642   5.064051     -1   5.042346
    ## 86           Stom   2  35.313986  35.336976     -1  35.325481
    ## 87           Mrrf   2  36.136389  36.190647      1  36.163518
    ## 88           Lhx3   2  26.200212  26.208289     -1  26.204250
    ## 89          Itga6   2  71.745616  71.858416      1  71.802016
    ## 90           Ehd4   2 120.089175 120.154606     -1 120.121890
    ## 91         Rassf2   2 131.989415 132.030258     -1 132.009836
    ## 92        Slc23a2   2 132.052496 132.145108     -1 132.098802
    ## 93          Ptpn1   2 167.932057 167.979385      1 167.955721
    ## 94            Fga   3  83.026153  83.033615      1  83.029884
    ## 95          Ikzf4  10 128.630843 128.645991     -1 128.638417
    ## 96       Rap1gds1   3 138.925906 139.075199     -1 139.000552
    ## 97           Klf3   5  64.803523  64.830129      1  64.816826
    ## 98         Crybb3   5 113.075839 113.081584     -1 113.078712
    ## 99         Scarb2   5  92.443873  92.505608     -1  92.474740
    ## 100          Gcc1   6  28.416091  28.428390     -1  28.422240
    ## 101          Etv6   6 134.035700 134.270158      1 134.152929
    ## 102        Rassf8   6 145.746748 145.821079      1 145.783914
    ## 103       St8sia1   6 142.821545 142.964452     -1 142.892998
    ## 104        Mtmr10   7  64.287653  64.340407      1  64.314030
    ## 105         Trpm1   7  64.153835  64.269775      1  64.211805
    ## 106         Adam9   8  24.949611  25.016922     -1  24.983266
    ## 107         Stat3  11 100.885098 100.939540     -1 100.912319
    ## 108          Cul5   9  53.614582  53.670014     -1  53.642298
    ## 109          Hcn4   9  58.823512  58.860955      1  58.842234
    ## 110        Plscr1   9  92.250057  92.271975      1  92.261016
    ## 111        Csrnp1   9 119.971166 119.977250     -1 119.974208
    ## 112        Col6a4   9 105.989454 106.096783     -1 106.043118
    ## 113          Cish   9 107.296026 107.302784      1 107.299405
    ## 114         Mon1a   9 107.888129 107.903125      1 107.895627
    ## 115         Mier3  13 111.686178 111.718596      1 111.702387
    ## 116        Clasp2   9 113.812586 113.919697      1 113.866142
    ## 117          St18   1   6.487231   6.860940      1   6.674086
    ## 118         Pqlc1  18  80.253292  80.292725      1  80.273008
    ## 119       Plekhm1  11 103.364275 103.412687     -1 103.388481
    ## 120          Ssh3  19   4.261668   4.269172     -1   4.265420
    ## 121         Rrp12  19  41.862852  41.896153     -1  41.879502
    ## 122          Helq   5 100.762145 100.798598     -1 100.780372
    ## 123         Tgfbi  13  56.609603  56.639339      1  56.624471
    ## 124        Depdc5   5  32.863701  32.994231      1  32.928966
    ## 125        Arid5a   1  36.307733  36.324029      1  36.315881
    ##     nearest.marker.id              biotype         module hotspot
    ## 1          5_31207687       protein_coding     lightgreen   chr13
    ## 2          8_31104305       protein_coding     lightgreen   chr13
    ## 3         10_30682488       protein_coding      royalblue   chr13
    ## 4        10_128145937       protein_coding      royalblue   chr13
    ## 5         3_144525411       protein_coding           grey   chr13
    ## 6         19_32252908       protein_coding      royalblue   chr13
    ## 7        11_107030096       protein_coding      royalblue   chr13
    ## 8         4_118424020       protein_coding      royalblue   chr13
    ## 9          9_57459238       protein_coding      royalblue   chr13
    ## 10       13_114014695       protein_coding          brown    <NA>
    ## 11       13_113222235       protein_coding          brown    <NA>
    ## 12       13_112978421       protein_coding           grey    <NA>
    ## 13       13_112927057       protein_coding          white    <NA>
    ## 14          6_5410616       protein_coding      turquoise   chr13
    ## 15        17_71718856       protein_coding           grey   chr13
    ## 16        3_127934331       protein_coding      royalblue   chr13
    ## 17       13_112689915       protein_coding            tan    <NA>
    ## 18        11_68640199       protein_coding      royalblue   chr13
    ## 19        X_102059421       protein_coding     lightgreen   chr13
    ## 20       13_112281956       protein_coding           grey    <NA>
    ## 21         3_90257553       protein_coding      royalblue   chr13
    ## 22        18_14917630       protein_coding     lightgreen   chr13
    ## 23         7_47035721       protein_coding          brown   chr13
    ## 24        1_192380604       protein_coding         yellow   chr13
    ## 25        17_35860095       protein_coding      royalblue   chr13
    ## 26         5_34342197       protein_coding      royalblue   chr13
    ## 27       11_116536263       protein_coding      royalblue   chr13
    ## 28        11_78062799       protein_coding      royalblue   chr13
    ## 29        12_78979554       protein_coding          green   chr13
    ## 30        1_153756004       protein_coding        magenta   chr13
    ## 31        1_157166180       protein_coding            red   chr13
    ## 32          5_5728638       protein_coding          brown   chr13
    ## 33        10_59200975           pseudogene     lightgreen    <NA>
    ## 34        12_17661910       protein_coding         yellow   chr13
    ## 35        18_44532060       protein_coding           grey   chr13
    ## 36        6_142824658 processed_transcript           grey   chr13
    ## 37         8_94563209       protein_coding        magenta   chr13
    ## 38       13_112903651       protein_coding      turquoise    <NA>
    ## 39        2_105152903       protein_coding      royalblue   chr13
    ## 40        1_130914982       protein_coding     lightgreen   chr13
    ## 41          X_5710894           pseudogene           grey   chr13
    ## 42        4_132886054           pseudogene     lightgreen   chr13
    ## 43       13_111832170 processed_transcript           grey    <NA>
    ## 44        11_33693427            antisense           grey   chr13
    ## 45       13_111883727              lincRNA           grey    <NA>
    ## 46        2_167322761       protein_coding      royalblue   chr13
    ## 47        5_143699438       protein_coding      royalblue   chr13
    ## 48        11_62213989       protein_coding      royalblue   chr13
    ## 49        10_41568532 processed_transcript           grey   chr13
    ## 50       11_115309983       protein_coding darkolivegreen   chr13
    ## 51         5_33008992       protein_coding     lightgreen   chr13
    ## 52         X_73407505            antisense      royalblue   chr13
    ## 53         2_10066961       protein_coding      royalblue    <NA>
    ## 54         18_5499717       sense_intronic           grey   chr13
    ## 55        10_59200975       protein_coding      royalblue   chr13
    ## 56        10_62123963       protein_coding           blue   chr13
    ## 57        12_14201843       protein_coding      royalblue   chr13
    ## 58        11_78062799       protein_coding     lightgreen   chr13
    ## 59        12_80754775       protein_coding      royalblue   chr13
    ## 60        13_40855519       protein_coding darkolivegreen    <NA>
    ## 61        2_181340625       protein_coding            red   chr13
    ## 62        10_76671270       protein_coding          black   chr13
    ## 63        17_56080483       protein_coding      royalblue   chr13
    ## 64       13_111772224       protein_coding      royalblue    <NA>
    ## 65       13_112472706       protein_coding      royalblue    <NA>
    ## 66       13_112640753       protein_coding           grey    <NA>
    ## 67       13_112833286       protein_coding          brown    <NA>
    ## 68       13_113079634       protein_coding          brown    <NA>
    ## 69       13_113344694       protein_coding    greenyellow    <NA>
    ## 70       13_114297395       protein_coding         salmon    <NA>
    ## 71        14_78445459       protein_coding      darkgreen   chr13
    ## 72        16_10992283       protein_coding      royalblue   chr13
    ## 73        11_95049781       protein_coding      royalblue   chr13
    ## 74        17_31326113       protein_coding      royalblue   chr13
    ## 75        19_16522825       protein_coding          brown   chr13
    ## 76        19_55848648       protein_coding      royalblue   chr13
    ## 77         9_59422066       protein_coding      darkgreen   chr13
    ## 78         1_65285084       protein_coding           pink   chr13
    ## 79         1_40487515       protein_coding      royalblue   chr13
    ## 80         1_40264675       protein_coding    lightyellow   chr13
    ## 81         9_37470492       protein_coding    lightyellow   chr13
    ## 82         3_20112430       protein_coding           grey   chr13
    ## 83        1_130884612       protein_coding        magenta   chr13
    ## 84        1_153748053       protein_coding         yellow   chr13
    ## 85          2_5037185       protein_coding darkolivegreen   chr13
    ## 86         2_35314728       protein_coding      royalblue   chr13
    ## 87         2_36215697       protein_coding     lightgreen   chr13
    ## 88         2_26244470       protein_coding      royalblue   chr13
    ## 89         2_71876558       protein_coding      royalblue   chr13
    ## 90        2_120096292       protein_coding     lightgreen   chr13
    ## 91        2_131822707       protein_coding      royalblue   chr13
    ## 92        2_132240855       protein_coding     lightgreen   chr13
    ## 93        2_167928557       protein_coding      royalblue   chr13
    ## 94         3_83032966       protein_coding      royalblue   chr13
    ## 95       10_128651252       protein_coding      royalblue   chr13
    ## 96        3_138988923       protein_coding      royalblue   chr13
    ## 97         5_64791826       protein_coding         grey60   chr13
    ## 98        5_113102018       protein_coding      royalblue   chr13
    ## 99         5_92505264       protein_coding     lightgreen   chr13
    ## 100        6_28433183       protein_coding     lightgreen   chr13
    ## 101       6_134185817       protein_coding      royalblue   chr13
    ## 102       6_145827741       protein_coding            red   chr13
    ## 103       6_142893152       protein_coding           grey   chr13
    ## 104        7_64311475       protein_coding           grey   chr13
    ## 105        7_64166021       protein_coding      royalblue   chr13
    ## 106        8_25064849       protein_coding      royalblue   chr13
    ## 107      11_100915707       protein_coding     lightgreen   chr13
    ## 108        9_53658586       protein_coding      turquoise   chr13
    ## 109        9_58980535       protein_coding   midnightblue   chr13
    ## 110        9_92257673       protein_coding      royalblue   chr13
    ## 111       9_119975235       protein_coding     lightgreen   chr13
    ## 112       9_106023952       protein_coding           grey   chr13
    ## 113       9_107283384       protein_coding      royalblue   chr13
    ## 114       9_107891367       protein_coding     orangered4   chr13
    ## 115      13_111712277       protein_coding      turquoise    <NA>
    ## 116       9_114004686       protein_coding     lightgreen   chr13
    ## 117         1_6673421       protein_coding     lightgreen   chr13
    ## 118       18_80242500       protein_coding     orangered4   chr13
    ## 119      11_103418350       protein_coding      royalblue   chr13
    ## 120        19_4267368       protein_coding         purple   chr13
    ## 121       19_41902281       protein_coding            red   chr13
    ## 122       5_100814844       protein_coding            red   chr13
    ## 123       13_56536955       protein_coding     darkorange    <NA>
    ## 124        5_32935215       protein_coding           grey   chr13
    ## 125        1_36365564       protein_coding      royalblue   chr13

The above tibble has local and nonlocal traits.

``` r
local_annot <- dataset.islet.rnaseq$lod.peaks %>%
  filter(chrom == hot_chrom, pos <= hot_upper, pos >= hot_lower) %>%
  left_join(y = dataset.islet.rnaseq$annots, by = c("annot.id" = "gene_id")) %>%
  filter(!duplicated(annot.id)) %>%
  filter(lod >= 7.18) %>%
  mutate(local_tx = chr ==hot_chrom) %>%
  filter(local_tx)
```

``` r
nonlocal_annot <- dataset.islet.rnaseq$lod.peaks %>%
  filter(chrom == hot_chrom, pos <= hot_upper, pos >= hot_lower) %>%
  left_join(y = dataset.islet.rnaseq$annots, by = c("annot.id" = "gene_id")) %>%
  filter(!duplicated(annot.id)) %>%
  filter(lod >= 7.18) %>% #genomewide threshold
  mutate(local_tx = chr == hot_chrom) %>%
  filter(!local_tx) %>% 
  filter(!is.na(hotspot))
```

``` r
dim(nonlocal_annot)
```

    ## [1] 104  16

We’ll use only those nonlocal traits that Keller et al. annotated as
mapping to a hotspot.

We also want to look at the local traits for Chr 13 hotspot. We have 19
local traits. We want to identify a subset of those for two-dimensional
scans. Let’s first limit to those near the hotspot, ie, within 5Mb of
either end.

``` r
local_annot2 <- local_annot %>%
  filter(middle <= hot_upper + 5, middle >= hot_lower - 5) %>%
  arrange(desc(lod))
```

Recall that *Pdx1* is the gene that Keller et al. identified as the key
mediator for Chr 5 hotspot. Its LOD is only 15, which seems to indicate
that I should be looking at a broad set of local traits, not just those,
eg., with LOD \> 40\!

Also, I’m limiting consideration of local genes to those with middle
within 5Mb of the hotspot.

Now, I make the two expression trait matrices that I need for use with
Condor.

``` r
hot_local <- dataset.islet.rnaseq$expr[, colnames(dataset.islet.rnaseq$expr) %in% local_annot2$annot.id]
hot_nonlocal <- dataset.islet.rnaseq$expr[, colnames(dataset.islet.rnaseq$expr) %in% nonlocal_annot$annot.id]
saveRDS(hot_local, "../data/Chr13hot_local.rds") # change this
saveRDS(hot_nonlocal, "../data/Chr13hot_nonlocal.rds") # change this
saveRDS(K$`13`, "../data/Chr13_kinship.rds") # change this - both
saveRDS(genoprobs$`13`, "../data/Chr13_aprobs.rds") # change this - both
```

Lastly, we need to identify the scan region for Chr 13.

We consider the physical
map.

``` r
(start_index <- which(map$`13`== map$`13`[map$`13` > hot_lower][1]) - 50) # change this - all 3!
```

    ## 13_110537457 
    ##         2814

``` r
(stop_index <- which(map$`13`== map$`13`[map$`13` > hot_upper][1]) + 50)
```

    ## 13_114554864 
    ##         3039

``` r
# change this - all 3!
```

``` r
(nsnp <- stop_index - start_index + 1)
```

    ## 13_114554864 
    ##          226

``` r
dim(local_annot2)
```

    ## [1] 17 16

## References

Keller et al. 2018. Genetics. Genetic drivers of pancreatic islet
function.

## Session info

``` r
devtools::session_info()
```

    ## ─ Session info ──────────────────────────────────────────────────────────
    ##  setting  value                       
    ##  version  R version 3.5.3 (2019-03-11)
    ##  os       macOS Mojave 10.14.5        
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  ctype    en_US.UTF-8                 
    ##  tz       America/Chicago             
    ##  date     2019-06-08                  
    ## 
    ## ─ Packages ──────────────────────────────────────────────────────────────
    ##  package     * version date       lib source                        
    ##  assertthat    0.2.1   2019-03-21 [1] CRAN (R 3.5.2)                
    ##  backports     1.1.4   2019-04-10 [1] CRAN (R 3.5.2)                
    ##  bit           1.1-14  2018-05-29 [1] CRAN (R 3.5.0)                
    ##  bit64         0.9-7   2017-05-08 [1] CRAN (R 3.5.0)                
    ##  blob          1.1.1   2018-03-25 [1] CRAN (R 3.5.0)                
    ##  broom         0.5.1   2018-12-05 [1] CRAN (R 3.5.0)                
    ##  callr         3.2.0   2019-03-15 [1] CRAN (R 3.5.3)                
    ##  cellranger    1.1.0   2016-07-27 [1] CRAN (R 3.5.0)                
    ##  cli           1.1.0   2019-03-19 [1] CRAN (R 3.5.3)                
    ##  colorspace    1.4-1   2019-03-18 [1] CRAN (R 3.5.3)                
    ##  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)                
    ##  data.table    1.12.0  2019-01-13 [1] CRAN (R 3.5.2)                
    ##  DBI           1.0.0   2018-05-02 [1] CRAN (R 3.5.0)                
    ##  desc          1.2.0   2018-05-01 [1] CRAN (R 3.5.0)                
    ##  devtools      2.0.2   2019-04-08 [1] CRAN (R 3.5.2)                
    ##  digest        0.6.19  2019-05-20 [1] CRAN (R 3.5.2)                
    ##  dplyr       * 0.8.0.1 2019-02-15 [1] CRAN (R 3.5.2)                
    ##  evaluate      0.14    2019-05-28 [1] CRAN (R 3.5.3)                
    ##  forcats     * 0.4.0   2019-02-17 [1] CRAN (R 3.5.2)                
    ##  fs            1.3.1   2019-05-06 [1] CRAN (R 3.5.2)                
    ##  generics      0.0.2   2018-11-29 [1] CRAN (R 3.5.0)                
    ##  ggplot2     * 3.1.1   2019-04-07 [1] CRAN (R 3.5.3)                
    ##  glue          1.3.1   2019-03-12 [1] CRAN (R 3.5.2)                
    ##  gtable        0.3.0   2019-03-25 [1] CRAN (R 3.5.2)                
    ##  haven         2.1.0   2019-02-19 [1] CRAN (R 3.5.2)                
    ##  hms           0.4.2   2018-03-10 [1] CRAN (R 3.5.0)                
    ##  htmltools     0.3.6   2017-04-28 [1] CRAN (R 3.5.0)                
    ##  httr          1.4.0   2018-12-11 [1] CRAN (R 3.5.0)                
    ##  jsonlite      1.6     2018-12-07 [1] CRAN (R 3.5.0)                
    ##  knitr         1.23    2019-05-18 [1] CRAN (R 3.5.2)                
    ##  lattice       0.20-38 2018-11-04 [1] CRAN (R 3.5.3)                
    ##  lazyeval      0.2.2   2019-03-15 [1] CRAN (R 3.5.3)                
    ##  lubridate     1.7.4   2018-04-11 [1] CRAN (R 3.5.0)                
    ##  magrittr      1.5     2014-11-22 [1] CRAN (R 3.5.0)                
    ##  memoise       1.1.0   2017-04-21 [1] CRAN (R 3.5.0)                
    ##  modelr        0.1.4   2019-02-18 [1] CRAN (R 3.5.2)                
    ##  munsell       0.5.0   2018-06-12 [1] CRAN (R 3.5.0)                
    ##  nlme          3.1-137 2018-04-07 [1] CRAN (R 3.5.3)                
    ##  pillar        1.3.1   2018-12-15 [1] CRAN (R 3.5.0)                
    ##  pkgbuild      1.0.3   2019-03-20 [1] CRAN (R 3.5.2)                
    ##  pkgconfig     2.0.2   2018-08-16 [1] CRAN (R 3.5.0)                
    ##  pkgload       1.0.2   2018-10-29 [1] CRAN (R 3.5.0)                
    ##  plyr          1.8.4   2016-06-08 [1] CRAN (R 3.5.0)                
    ##  prettyunits   1.0.2   2015-07-13 [1] CRAN (R 3.5.0)                
    ##  processx      3.3.1   2019-05-08 [1] CRAN (R 3.5.2)                
    ##  ps            1.3.0   2018-12-21 [1] CRAN (R 3.5.0)                
    ##  purrr       * 0.3.2   2019-03-15 [1] CRAN (R 3.5.3)                
    ##  qtl2        * 0.18    2019-02-08 [1] local                         
    ##  R6            2.4.0   2019-02-14 [1] CRAN (R 3.5.2)                
    ##  Rcpp          1.0.1.3 2019-05-29 [1] Github (RcppCore/Rcpp@6062d56)
    ##  readr       * 1.3.1   2018-12-21 [1] CRAN (R 3.5.0)                
    ##  readxl        1.3.0   2019-02-15 [1] CRAN (R 3.5.2)                
    ##  remotes       2.0.4   2019-04-10 [1] CRAN (R 3.5.2)                
    ##  rlang         0.3.4   2019-04-07 [1] CRAN (R 3.5.3)                
    ##  rmarkdown     1.13    2019-05-22 [1] CRAN (R 3.5.2)                
    ##  rprojroot     1.3-2   2018-01-03 [1] CRAN (R 3.5.0)                
    ##  RSQLite       2.1.1   2018-05-06 [1] CRAN (R 3.5.0)                
    ##  rstudioapi    0.10    2019-03-19 [1] CRAN (R 3.5.2)                
    ##  rvest         0.3.2   2016-06-17 [1] CRAN (R 3.5.0)                
    ##  scales        1.0.0   2018-08-09 [1] CRAN (R 3.5.0)                
    ##  sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.5.0)                
    ##  stringi       1.4.3   2019-03-12 [1] CRAN (R 3.5.2)                
    ##  stringr     * 1.4.0   2019-02-10 [1] CRAN (R 3.5.2)                
    ##  testthat      2.1.1   2019-04-23 [1] CRAN (R 3.5.3)                
    ##  tibble      * 2.1.1   2019-03-16 [1] CRAN (R 3.5.3)                
    ##  tidyr       * 0.8.3   2019-03-01 [1] CRAN (R 3.5.2)                
    ##  tidyselect    0.2.5   2018-10-11 [1] CRAN (R 3.5.0)                
    ##  tidyverse   * 1.2.1   2017-11-14 [1] CRAN (R 3.5.0)                
    ##  usethis       1.5.0   2019-04-07 [1] CRAN (R 3.5.2)                
    ##  withr         2.1.2   2018-03-15 [1] CRAN (R 3.5.0)                
    ##  xfun          0.7     2019-05-14 [1] CRAN (R 3.5.2)                
    ##  xml2          1.2.0   2018-01-24 [1] CRAN (R 3.5.0)                
    ##  yaml          2.2.0   2018-07-25 [1] CRAN (R 3.5.0)                
    ## 
    ## [1] /Library/Frameworks/R.framework/Versions/3.5/Resources/library
