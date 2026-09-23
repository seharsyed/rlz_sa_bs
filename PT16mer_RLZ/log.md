
some random 100
```bash
 ./pt16_v2 --reference ../../RLZ_DATA/cleaned-references/GCA_000005845.2_ASM584v2.cleaned --suffix-array ../../RLZ_DATA/cleaned-references/GCA_000005845.2_ASM584v2.cleaned.sa --filenames ../../RLZ_DATA/cleaned_input_list_100_first.txt --results results/pt16_v2_ecoli100.csv
 ```
lex 100 first with the same reference as above

```bash
 ./pt16 --reference ../../RLZ_DATA/cleaned-references/GCA_000005845.2_ASM584v2.cleaned --suffix-array ../../RLZ_DATA/cleaned-references/GCA_000005845.2_ASM584v2.cleaned.sa --filenames ../../RLZ_DATA/input-list-100-lex.txt --results results/pt16_v2_ecoli100.csv
 ```

using diff reference
 GCA_001012175.1_CFSAN026787_02.0.txt.sa
 ```bash
 ./pt16_v2 --reference ../../RLZ_DATA/cleaned-references/GCA_000005845.2_ASM584v2.cleaned --suffix-array ../../RLZ_DATA/cleaned-references/GCA_000005845.2_ASM584v2.cleaned.sa --filenames ../../RLZ_DATA/cleaned_input_list_100_first.txt --results results/pt16_v2_ecoli100.csv
 ```
 ```bash
 ./pt16_v2 --reference ../../RLZ_DATA/cleaned-references/GCA_000005845.2_ASM584v2.cleaned --suffix-array ../../RLZ_DATA/cleaned-references/GCA_000005845.2_ASM584v2.cleaned.sa --filenames ../../RLZ_DATA/input-list-100-lex.txt --results results/pt16_v2_ecoli100.csv
 ```

yet another reff
 ../../RLZ_DATA/cleaned-references/GCA_000179135.1_ASM17913v1.txt.sa
 ```bash
 ./pt16_v2_non_interleaved --reference ../../RLZ_DATA/cleaned_inputs/GCA_000179135.1_ASM17913v1.txt --suffix-array ../../RLZ_DATA/cleaned-references/GCA_000179135.1_ASM17913v1.txt.sa --filenames ../../RLZ_DATA/cleaned_input_list_100_first.txt --results results/pt16_v2_ecoli100.csv
 ```
 ```bash
 ./pt16_v2 --reference ../../RLZ_DATA/cleaned-references/GCA_000005845.2_ASM584v2.cleaned --suffix-array ../../RLZ_DATA/cleaned-references/GCA_000005845.2_ASM584v2.cleaned.sa --filenames ../../RLZ_DATA/ --results results/pt16_v2_ecoli100.csv
 ```

100 first
```text
../../RLZ_DATA/cleaned_inputs/GCA_000006665.1_ASM666v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000007445.1_ASM744v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000008865.2_ASM886v2.txt
../../RLZ_DATA/cleaned_inputs/GCA_000009565.2_ASM956v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000010245.1_ASM1024v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000010385.1_ASM1038v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000010485.1_ASM1048v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000010745.1_ASM1074v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000010765.1_ASM1076v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000013265.1_ASM1326v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000013305.1_ASM1330v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000014845.1_ASM1484v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000017745.1_ASM1774v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000017765.1_ASM1776v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000017985.1_ASM1798v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000019385.1_ASM1938v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000019425.1_ASM1942v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000019645.1_ASM1964v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000021125.1_ASM2112v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000022225.1_ASM2222v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000022345.1_ASM2234v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000022665.2_ASM2266v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000023365.1_ASM2336v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000023665.1_ASM2366v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000025165.1_ASM2516v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000025745.1_ASM2574v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000026245.1_ASM2624v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000026265.1_ASM2626v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000026285.2_ASM2628v2.txt
../../RLZ_DATA/cleaned_inputs/GCA_000026305.1_ASM2630v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000026325.2_ASM2632v2.txt
../../RLZ_DATA/cleaned_inputs/GCA_000026345.1_ASM2634v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000026545.1_ASM2654v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000027125.1_ASM2712v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000091005.1_ASM9100v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000146735.1_ASM14673v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000147755.2_ASM14775v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000147855.3_ASM14785v3.txt
../../RLZ_DATA/cleaned_inputs/GCA_000148365.1_ASM14836v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000148605.1_ASM14860v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000155005.1_ASM15500v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000155125.1_ASM15512v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000157115.2_Escherichia_sp_3_2_53FAA_V2.txt
../../RLZ_DATA/cleaned_inputs/GCA_000158395.1_ASM15839v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000159295.1_ASM15929v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000163155.1_ASM16315v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000163175.1_ASM16317v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000163195.1_ASM16319v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000163215.1_ASM16321v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000163235.1_ASM16323v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164195.1_ASM16419v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164215.1_ASM16421v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164235.1_ASM16423v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164255.1_ASM16425v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164275.1_ASM16427v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164295.1_ASM16429v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164315.1_ASM16431v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164335.1_ASM16433v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164355.1_ASM16435v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164375.1_ASM16437v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164415.1_ASM16441v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164435.1_ASM16443v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164455.1_ASM16445v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164475.1_ASM16447v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164495.1_ASM16449v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164515.1_ASM16451v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164535.1_ASM16453v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164555.1_ASM16455v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164575.1_ASM16457v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164595.1_ASM16459v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000164615.1_ASM16461v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000165655.2_ASM16565v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000166535.3_ASM16653v3.txt
../../RLZ_DATA/cleaned_inputs/GCA_000166555.2_ASM16655v2.txt
../../RLZ_DATA/cleaned_inputs/GCA_000166575.2_ASM16657v2.txt
../../RLZ_DATA/cleaned_inputs/GCA_000166595.2_ASM16659v2.txt
../../RLZ_DATA/cleaned_inputs/GCA_000166615.2_ASM16661v2.txt
../../RLZ_DATA/cleaned_inputs/GCA_000167815.1_ASM16781v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000167835.1_ASM16783v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000167855.1_ASM16785v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000167875.1_ASM16787v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000167895.3_ASM16789v3.txt
../../RLZ_DATA/cleaned_inputs/GCA_000176695.2_ASM17669v2.txt
../../RLZ_DATA/cleaned_inputs/GCA_000176815.1_ASM17681v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000178315.1_ASM17831v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000178695.1_ASM17869v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000178715.1_ASM17871v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000178735.1_ASM17873v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000178755.1_ASM17875v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000178775.1_ASM17877v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000178795.1_ASM17879v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000179075.1_ASM17907v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000179095.1_ASM17909v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000179115.1_ASM17911v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000179135.1_ASM17913v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000179155.1_ASM17915v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000179175.1_ASM17917v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000179795.1_ASM17979v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000181735.1_ASM18173v1.txt
../../RLZ_DATA/cleaned_inputs/GCA_000181755.1_ASM18175v1.txt
```