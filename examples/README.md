# Examples

## Full

| Filetype  | Keyword | Description                                                               | Type     |
| ---       | ---     | ---                                                                       | ---      |
| BAM       | `bam`   | Aligned reads                                                             | Required |
| SV VCF    | `vcf`   | SV VCF produced by an external caller [_used to filter reads for GRIDSS_] | Optional |

```text
id        subject_name   sample_name          sample_type  filetype  filepath
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        bam       /path/to/tumor_bam/sample_one_tumor.bam
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        sv        /path/to/tumor_sv_vcf/sample_one_tumor.vcf.gz
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       bam       /path/to/normal_bam/sample_one_normal.bam
STWO-1    SUBJECT_TWO    SAMPLE_TWO_TUMOR-1   tumor        bam       /path/to/tumor_bam/sample_two_tumor_one.bam
STWO-1    SUBJECT_TWO    SAMPLE_TWO_NORMAL    normal       bam       /path/to/normal_bam/sample_two_normal.bam
STWO-2    SUBJECT_TWO    SAMPLE_TWO_TUMOR-2   tumor        bam       /path/to/tumor_bam/sample_two_tumor_two.bam
STWO-2    SUBJECT_TWO    SAMPLE_TWO_NORMAL    normal       bam       /path/to/normal_bam/sample_two_normal.bam
STRHEE-1  SUBJECT_THREE  SAMPLE_THREE_TUMOR   tumor        bam       /path/to/tumor_bam/sample_three_tumor.bam
STRHEE-1  SUBJECT_THREE  SAMPLE_THREE_TUMOR   tumor        sv        /path/to/tumor_sv_vcf/sample_three_tumor.vcf.gz
STRHEE-1  SUBJECT_THREE  SAMPLE_THREE_NORMAL  normal       bam       /path/to/normal_bam/sample_three_normal.bam
```

## GRIDSS

See [Full section](#full)

## PURPLE

| Filetype                      | Keyword          | Description                 | Type     |
| ---                           | ---              | ---                         | ---      |
| AMBER directory               | `amber_dir`      | AMBER output directory      | Required |
| COBALT directory              | `cobalt_dir`     | COBALT output directory     | Required |
| GRIPSS SV VCF (hard filtered) | `gripss_hard_sv` | Hard filtered GRIPSS SV VCF | Required |
| GRIPSS SV VCF (soft filtered) | `gripss_soft_sv` | Soft filtered GRIPSS SV VCF | Required |
| SNV/MNV and INDEL VCF         | `smlv`           | Small SNV/MNV VCF           | Optional |

```text
id        subject_name   sample_name          sample_type  filetype         filepath
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        amber_dir        /path/to/amber_dir/
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        cobalt_dir       /path/to/cobalt_dir/
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        gripss_hard_sv   /path/to/tumor_gripss_hard_sv/sample_one_tumor.vcf.gz
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        gripss_soft_sv   /path/to/tumor_gripss_soft_sv/sample_one_tumor.vcf.gz
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        smlv             /path/to/tumor_smlv_vcf/sample_one_tumor.vcf.gz
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       smlv             /path/to/normal_smlv_vcf/sample_one_normal.vcf.gz
```

## LINX

| Filetype                      | Keyword          | Description                                   | Type     |
| ---                           | ---              | ---                                           | ---      |
| PURPLE directory              | `purple_dir`     | PURPLE output directory [_LINX somatic_]      | Required |
| GRIPSS SV VCF (hard filtered) | `gripss_hard_sv` | Hard filtered GRIPSS SV VCF [_LINX germline_] | Required |

```text
id        subject_name   sample_name          sample_type  filetype         filepath
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        purple_dir       /path/to/purple_dir/
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       gripss_hard_sv   /path/to/normal_gripss_hard_sv/sample_one_normal.vcf.gz
```

## GRIDSS-PURPLE-LINX

| Filetype              | Keyword | Description                                                               | Type     |
| ---                   | ---     | ---                                                                       | ---      |
| BAM                   | `bam`   | Aligned reads                                                             | Required |
| SV VCF                | `vcf`   | SV VCF produced by an external caller [_used to filter reads for GRIDSS_] | Optional |
| SNV/MNV and INDEL VCF | `smlv`  | Small SNV/MNV VCF                                                         | Optional |

```text
id        subject_name   sample_name          sample_type  filetype  filepath
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        bam       /path/to/tumor_bam/sample_one_tumor.bam
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        sv        /path/to/tumor_sv_vcf/sample_one_tumor.vcf.gz
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       bam       /path/to/normal_bam/sample_one_normal.bam
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        smlv      /path/to/tumor_smlv_vcf/sample_one_tumor.vcf.gz
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       smlv      /path/to/normal_smlv_vcf/sample_one_normal.vcf.gz
```


## LILAC

| Filetype              | Keyword      | Description             | Type     |
| ---                   | ---          | ---                     | ---      |
| BAM                   | `bam`        | Aligned reads           | Required |
| PURPLE directory      | `purple_dir` | PURPLE output directory | Required |

```text
id        subject_name   sample_name          sample_type  filetype     filepath
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        bam          /path/to/tumor_bam/sample_one_tumor.bam
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        purple_dir   /path/to/purple_dir/
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       bam          /path/to/normal_bam/sample_one_normal.bam
```

## TEAL

| Filetype              | Keyword      | Description             | Type     |
| ---                   | ---          | ---                     | ---      |
| BAM                   | `bam`        | Aligned reads           | Required |
| COBALT directory      | `cobalt_dir` | COBALT output directory | Required |
| PURPLE directory      | `purple_dir` | PURPLE output directory | Required |

```text
id        subject_name   sample_name          sample_type  filetype     filepath
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        bam          /path/to/tumor_bam/sample_one_tumor.bam
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        purple_dir   /path/to/purple_dir/
SONE-1    SUBJECT_ONE    SAMPLE_ONE_TUMOR     tumor        cobalt_dir   /path/to/cobalt_dir/
SONE-1    SUBJECT_ONE    SAMPLE_ONE_NORMAL    normal       bam          /path/to/normal_bam/sample_one_normal.bam
```
