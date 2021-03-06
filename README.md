# Variant-Normalization

[Variant Normalization](https://genome.sph.umich.edu/wiki/Variant_Normalization)



The normalization of a variant representation in VCF consists of two parts: **parsimony** and **left alignment** pertaining to the nature of a variant's length and position respectively.

- Parsimony means representing a variant in as few nucleotides as possible without reducing the length of any allele to 0.

- Left aligning a variant means shifting the start position of that variant to the left till it is no longer possible to do so.

- A variant is normalized if and only if it is parsimonious and left-aligned

## National Center for Clinical Laboratories(NCCL) 3'rule

> 所有变异都应按照靠近基因转录方向末端（3’ 端）的位置来表述（3’rule）
> 
> 正链基因（以正链为模板转录，出现在重复序列上的突变根据转录本上的 3’rule 发生在 3’ 端，在基因组层面，也按照转录本上 3’ 端书写，即突变发生在重复序列的右侧）
> 
> 负链基因（以负链为模板转录，出现在重复序列上的突变根据转录本上的 3’rule，发生在 3’ 端，在基因组层面，也按照转录本上 3’ 端书写，即突变发生在重复序列的左侧）

### Example of Deletion variant normalization

Raw variant :

| Chr  | Pos      | Ref   | Alt  | Gene | Strand | Type     | Transcript  |
| ---- | -------- | ----- | ---- | ---- | ------ | -------- | ----------- |
| 3    | 10183820 | CCCTA | C    | VHL  | +      | Deletion | NM_000551.4 |

Equivalent deletions :

![](VHL.c.292_295del.PNG)

Variant Normalization result :

| Chr  | Pos      | Ref  | Alt  | Gene | Strand | Type     | Transcript  |
| ---- | -------- | ---- | ---- | ---- | ------ | -------- | ----------- |
| 3    | 10183823 | TACC | -    | VHL  | +      | Deletion | NM_000551.4 |



### Example of Insertion variant normalization

Raw variant :

| Chr  | Pos      | Ref  | Alt   | Gene  | Strand | Type      | Transcript  |
| ---- | -------- | ---- | ----- | ----- | ------ | --------- | ----------- |
| X    | 48649699 | C    | CTACT | GATA1 | +      | Insertion | NM_002049.4 |

Equivalent insertions :

![](GATA1.c.184_187dup.PNG)

Variant Normalization result:

| Chr  | Pos      | Ref  | Alt  | Gene  | Strand | Type      | Transcript  |
| ---- | -------- | ---- | ---- | ----- | ------ | --------- | ----------- |
| X    | 48649703 | -    | TACT | GATA1 | +      | Insertion | NM_002049.4 |

