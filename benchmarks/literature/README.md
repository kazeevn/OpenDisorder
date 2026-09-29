# Literature cases: structures disputed as "new"

`python run.py --db <index>` regenerates `results.csv`. Verdicts: **known** = `pristine` / `disordered_realization`; **suspicious** = `disordered_family` / `substitutional_ordering` / `sibling_ordering`; **novel** = no COD match. The papers compare against ICSD; many of the phases they name are not in COD, so "novel" here means *no COD evidence*, not a rebuttal of the paper.

## Cheetham&Seshadri2024

| formula             | id         | known phase according to the paper                                                | verdict    | kind                    | COD match         |
|:--------------------|:-----------|:----------------------------------------------------------------------------------|:-----------|:------------------------|:------------------|
| K2SmLuB24           | 00edc694b5 | LaB6 (hexaborides) [Pm-3m] (Kiessling 1950)                                       | novel      | none                    |                   |
| K3Nd(AsO4)2         | 98bccb878a | K3Nd(PO4)2 (phosphate analogue) [P2_1/m] (Hong & Chinn 1976 (ref 8))              | novel      | none                    |                   |
| DyHo2Rh9            | 00e7694a54 | (Dy,Ho)Rh3 ~ CeNi3 [P6_3/mmc (CeNi3)] (Cromer&Olsen 1959; YRh2Si ordered variant) | novel      | none                    |                   |
| TbHo2Tm9Co4         | 00e086512f | Ln3Co e.g. Ho3Co [Pnma] (Buschow 1969 (ref 12))                                   | suspicious | substitutional_ordering | 1524134 (Co Ho3)  |
| Tb16Ho(ErIr4)3      | 0057fbfd06 | Ln5Ir3 e.g. Tb5Ir3 [P4/ncc] (Le Roy 1980 (ref 13))                                | novel      | none                    |                   |
| Dy6Y2Ho11Lu(Cd3Ru)2 | 05e70f1374 | none named (HEA-like) [] ()                                                       | novel      | none                    |                   |
| Dy4Tc3B8Ru          | 7c5484aa4a | B2DyRu [Pnma] (ICSD 100774)                                                       | novel      | none                    |                   |
| Y3SnPd3Pb2          | cce1e6134a | PbPtGd [P-62m] (ICSD 230854)                                                      | novel      | none                    |                   |
| Hf4Zr11Nb(Ir2N)4    | ebc757e2d1 | Zr4Ir2N (Ti2Ni-filled, Zr/Hf/Nb disordered) [Fd-3m] (ICSD 640826)                 | novel      | none                    |                   |
| Tb4Er4Tc(B8Os)3     | 6bffec4d39 | B6Er2Os [Pban] (ICSD 603495)                                                      | novel      | none                    |                   |
| Ta16AlV3Te16        | 7035244016 | FeTa4Te4 [Pbam] (ICSD 79796)                                                      | novel      | none                    |                   |
| La2Ga4Pt3Rh         | cace2fa555 | Ce2Pt1.5Ga6.5 [I4/mmm] (ICSD 621150)                                              | novel      | none                    |                   |
| SrTm(Cu2P)4         | f595e58e43 | CaCu4P2 [P4_2/mnm] (ICSD 62553)                                                   | novel      | none                    |                   |
| Tb9DySi5P           | a7c8bcfdf9 | Gd5Si3 [P6_3/mcm] (ICSD 99641)                                                    | suspicious | substitutional_ordering | 1536977 (Dy5 Si3) |
| AcTl2P2S7           | 093058048a | La2P4S14Rb4 (RbLa2(P2S6)(PS4)2) [P2_1/c] (ICSD 25121)                             | novel      | none                    |                   |
| Nb10GaOs4           | d95b844297 | Ge3.9Nb20.4Rh5.9 (sigma phase) [P4_2/mnm] (ICSD 230694)                           | novel      | none                    |                   |
| K3Nd(AsO4)2         | 0418958d3f | K3Nd(PO4)2 (phosphate analogue) [P2_1/m] (Hong & Chinn 1976 (ref 8))              | novel      | none                    |                   |

## Leeman2024

| formula            | id         | known phase according to the paper                                                      | verdict    | kind                    | COD match                                  |
|:-------------------|:-----------|:----------------------------------------------------------------------------------------|:-----------|:------------------------|:-------------------------------------------|
| Ba2ZrSnO6          | mp-1228067 | BaZr0.5Sn0.5O3 (+BaSnO3,SnO2,ZrO2) [Pm-3m] (ICSD 43137, 188149)                         | suspicious | substitutional_ordering | 1532743 (Ba O3 Zr)                         |
| Ba6Na2Ta2V2O17     | mp-1214664 | BaNaVO4 + NaVO3 [P-3m1; C2/c] (ICSD 130002, 2103)                                       | novel      | none                    |                                            |
| Ba6Na2V2Sb2O17     | mp-1214658 | BaNaVO4 + NaVO3 [P-3m1; C2/c] (ICSD 130002, 2103)                                       | novel      | none                    |                                            |
| Ba9Ca3La4(Fe4O15)2 | mp-1228537 | Ba4.5Ca1.5La2Fe4O15 [P6_3mc] (ICSD 72336)                                               | known      | disordered_realization  | 2002594 (Ba4.5 Ca1.5 Fe4 La2 O15)          |
| CaCo(PO3)4         | mp-1045787 | Co2P2O7 + Ca3(PO4)2 / M2(PO3)4 type [P2_1/c] (ICSD 2830, 923)                           | suspicious | substitutional_ordering | 1542153 (Co2 O12 P4)                       |
| CaFe2P2O9          | mp-1040941 | CaFe2(PO4)2O crocobelonite [] (Britvin 2023 (ref 44))                                   | novel      | none                    |                                            |
| CaGd2Zr(GaO3)4     | mp-686296  | Ca0.95Zr0.95Gd2.05Ga4.05O12 [Ia-3d] (ICSD 202850)                                       | novel      | none                    |                                            |
| CaMn(PO3)4         | mp-1045779 | Mn2P4O12 [C2/c] (ICSD 412558)                                                           | suspicious | substitutional_ordering | 1534609 (Mn2 O12 P4)                       |
| CaNi(PO3)4         | mp-1045813 | CaNi3(P2O7)2 + Ca2P2O7 / Ni2(PO3)4 [P2_1/c; P4_1] (ICSD 74046, 14313)                   | novel      | none                    |                                            |
| FeSb3Pb4O13        | mp-1224890 | Pb2(Fe0.5Sb1.5)O6.5 [Fd-3m] (ICSD 60805; Cascales 1985)                                 | novel      | none                    |                                            |
| Hf2Sb2Pb4O13       | mp-1224490 | Pb2(HfSb)O6.5 + PbHfO3 [Fd-3m] (ICSD 62723, 174110; Cascales 1986)                      | novel      | none                    |                                            |
| InSb3(PO4)6        | mp-1224667 | Sb0.5In0.5P2O7; Sb(Sb0.5In0.5)(PO4)3 [P2_1/n; Pna2_1] (ICSD 166834, 166835)             | novel      | none                    |                                            |
| InSb3Pb4O13        | mp-1223746 | Pb2(In0.5Sb1.5)O6.5 + In2O3 [Fd-3m] (ICSD 14388/41119; Cascales 1985)                   | novel      | none                    |                                            |
| K2TiCr(PO4)3       | mp-1224541 | K2Ti2(PO4)3 + Cr2O3 [P2_13] (ICSD 202888, 626479)                                       | known      | disordered_realization  | 2103459 (Cr K2 O12 P3 Ti)                  |
| K4MgFe3(PO4)5      | mp-532755  | K4MgFe3(PO4)5 [P-42_1c] (ICSD 161484; Hidouri 2008)                                     | novel      | none                    |                                            |
| K4TiSn3(PO5)4      | mp-1224290 | K(Ti0.25Sn0.75)OPO4 [Pna2_1] (ICSD 250088)                                              | known      | disordered_realization  | 1537558 (K O5 P Sn0.75 Ti0.25)             |
| KBaGdWO6           | mp-1523079 | Ba2GdWO6 + Gd2O3 [Fm-3m] (ICSD 138973, 40473)                                           | novel      | none                    |                                            |
| KBaPrWO6           | mp-1523149 | Ba11W4O23 + Pr6O11 [Fd-3m] (ICSD 418207)                                                | novel      | none                    |                                            |
| KMn3O6             | mp-1016190 | K0.3MnO2 / K0.48Mn1.94O5.18 [P6_3/mmc] (ICSD 240249; Kim 1999)                          | novel      | none                    |                                            |
| KNa2Ga3(SiO4)3     | mp-1211711 | NaGaSiO4 [P2_1/n] (ICSD 46861)                                                          | novel      | none                    |                                            |
| KNaP6(PbO3)8       | mp-1223429 | (Na,K)Pb4(PO4)3 + Pb8O5(PO4)2 [P6_3/m] (ICSD 182501, 98702)                             | suspicious | substitutional_ordering | 1528728 (K2 O24 P6 Pb8)                    |
| KNaTi2(PO5)2       | mp-1211611 | Na0.5(Na,K)0.5TiOPO4 (Na-doped KTP) [Pna2_1] (ICSD 59284; Norberg 2003)                 | suspicious | disordered_family       | 2103915 (K1.886 Na0.114 O10 P2 Ti2)        |
| KPr9(Si3O13)2      | mp-1223421 | Pr9.33Si6O26-type (EVA) [] ()                                                           | novel      | none                    |                                            |
| Mg3MnNi3O8         | mp-1222170 | Mg3MnNi3O8 / (Ni1-xMgx)6MnO8 murdochite [Fm-3m] (Taguchi 1995 (ICSD); ICSD 84517, 9866) | suspicious | substitutional_ordering | 9013975 (Mn Ni6 O8)                        |
| Mg3NiO4            | mp-1099253 | MgNiO2 / MgO-NiO solid solution [Fm-3m] (ICSD 290603)                                   | suspicious | substitutional_ordering | 1010093 (Ni O)                             |
| MgCuP2O7           | mp-1041741 | Mg2P2O7 / (Mg,Cu)2P2O7 disordered [C2/m] (ICSD 20295)                                   | novel      | none                    |                                            |
| MgNi(PO3)4         | mp-1045786 | Mg3(PO4)2, Ni3(PO4)2 / M2(PO3)4 [P2_1/n] (ICSD 31005, 153159, 158525)                   | novel      | none                    |                                            |
| MgTi2NiO6          | mp-1221952 | (Ni0.5Mg0.5)TiO3 [R-3] (ICSD 171583)                                                    | known      | disordered_realization  | 9007388 (Mg0.5 Ni0.5 O3 Ti)                |
| MgTi4(PO4)6        | mp-1222070 | Mg0.5Ti2(PO4)3 [R-3c] (ICSD 74287; Barth 1993)                                          | novel      | none                    |                                            |
| MgV4Cu3O14         | mp-1222158 | (Cu1.5Mg0.5)V2O7 [C2/c] (ICSD 69731)                                                    | known      | disordered_realization  | 2002572 (Cu1.5 Mg0.5 O7 V2)                |
| Mn2VPO7            | mp-1210613 | Mn2P2O7 / Mn2V2O7 [C2/m] (ICSD 47136)                                                   | known      | disordered_realization  | 1539969 (Mn2 O7 P1.09 V0.91)               |
| Mn4Zn3(NiO6)2      | mp-1222033 | (Zn0.759Mn0.241)(Mn1.35Ni0.65)O4 + NiO [Fd-3m] (ICSD 92223; Guillemet-Fritsch 2000)     | known      | disordered_realization  | 1526211 (Mn1.338 Ni0.65 O4 Zn1.012)        |
| Mn7(P2O7)4         | mp-778008  | Mn2(PO3)4 [C2/c] (ICSD 145534)                                                          | novel      | none                    |                                            |
| MnAgO2             | mp-996995  | MnAgO2 [C2/m] (ICSD 139006; Griesemer 2021)                                             | novel      | none                    |                                            |
| Na3Ca18Fe(PO4)14   | mp-725491  | Na3Ca18Fe(PO4)14 (whitlockite) [R3c] (ICSD 85103)                                       | novel      | none                    |                                            |
| Na7Mg7Fe5(PO4)12   | mp-1173791 | Na2Mg2Fe(PO4)3 (alluaudite) [C2/c] (ICSD 138263)                                        | known      | disordered_realization  | 2201754 (Fe1.21 Mg1.79 Na1.79 O12 P3)      |
| NaCaMgFe(SiO3)4    | mp-1221075 | (Ca,Na)(Mg,Fe)Si2O6 pyroxene [C2/c] (ICSD 75294, 417169)                                | known      | disordered_realization  | 9014598 (Ca0.58 Fe0.5 Mg0.5 Na0.42 O6 Si2) |
| NaMnFe(PO4)2       | mp-1173592 | Na-Mn-Fe alluaudite [C2/c] (ICSD 138263 (Daidouh 2002))                                 | known      | disordered_realization  | 9003657 (Fe1.5 Mn1.538 Na1.56 O12 P3)      |
| Sn2Sb2Pb4O13       | mp-1219056 | Pb2SnSbO6.5 + SnO2 [Fd-3m] (ICSD 62722, 9163)                                           | novel      | none                    |                                            |
| Y3In2(GaO4)3       | mp-1207946 | Y3In2Ga3O12 [Ia-3d] (ICSD 54664; Li & Zhong 2022)                                       | novel      | none                    |                                            |
| Zn2Cr3FeO8         | mp-1215741 | ZnFe(2-x)CrxO4 / Zn(FeCrO4) [Fd-3m] (ICSD 167362, 81207; Oles 1970)                     | suspicious | substitutional_ordering | 1011278 (Cr2 O4 Zn)                        |
| Zn3Ni4(SbO6)2      | mp-1216023 | Zn-Ni-Sb spinel e.g. Zn6NiSb2O12 (Gama 2003) [Fd-3m] (ICSD 109465/109468/173996)        | suspicious | substitutional_ordering | 1007077 (O32.004 Sb5.334 Zn18.669)         |
| Zr2Sb2Pb4O13       | mp-1215826 | Pb2(ZrSb)O6.5 [Fd-3m] (ICSD 62721; Cascales 1986)                                       | novel      | none                    |                                            |

## Juelsholt2026

| formula              | id           | known phase according to the paper   | verdict    | kind                   | COD match                    |
|:---------------------|:-------------|:-------------------------------------|:-----------|:-----------------------|:-----------------------------|
| TaCr2O6              | TaCr2O6      | Ta0.5Cr0.5O2 rutile (COD 9011856)    | suspicious | disordered_family      | 9011856 (Cr0.5 O2 Ta0.5)     |
| Ta0.66666Cr1.33334O4 | Ta1/3Cr2/3O2 | Ta0.5Cr0.5O2 rutile (COD 9011856)    | suspicious | disordered_family      | 2002668 (Cr0.761 O6 Ta2.033) |
| TaCrO4               | mp-753467    | Ta0.5Cr0.5O2 rutile                  | known      | disordered_realization | 9011856 (Cr0.5 O2 Ta0.5)     |
| TaCrO4               | mp-756340    | Ta0.5Cr0.5O2 rutile                  | known      | disordered_realization | 9011856 (Cr0.5 O2 Ta0.5)     |
| TaCrO4               | agm003743920 | Ta0.5Cr0.5O2 rutile                  | known      | disordered_realization | 9011856 (Cr0.5 O2 Ta0.5)     |
| Ta2CrO6              | mp-21512     | Ta2CrO6 (monoclinic above x_Ta 0.65) | known      | pristine               | 1001091 (Cr O6 Ta2)          |
| Mn5Fe3O8             | Mn5Fe3O8     | (Mn,Fe)O rock salt                   | suspicious | disordered_family      | 4344257 (Fe0.437 Mn0.556 O)  |

## Zeni2025 SI

| formula   | id      | known phase according to the paper   | verdict   | kind   | COD match   |
|:----------|:--------|:-------------------------------------|:----------|:-------|:------------|
| Cr2MoO6   | Cr2MoO6 | Cr2MoO6 (Collomb et al.)             | novel     | none   |             |
