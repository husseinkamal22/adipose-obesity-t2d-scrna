# Legacy / exploratory scripts

These scripts are **previous workflow versions**. They are kept for provenance
only and are **not** part of the final reproducible pipeline.

| File | What it was |
| ---- | ----------- |
| `02_downstream_analysis.R` | Early single-sample downstream pass (annotation + composition) for BRI-1456 only. |
| `02_downstream_analysis_v1_pilot.R` | Larger exploratory pilot that mixed annotation, composition, inflammation scoring and object saving in one file. |

**Do not run these scripts.** They predate the numbered pipeline in `../`:

```
01_primary_BRI145*.R  ->  02_annotation_and_composition.R
                       ->  03_inflammation_scoring.R
                       ->  04_cellchat_analysis.R
                       ->  05_cellchat_results.R
                       ->  06_final_figures.R
```

They are not executed by the pipeline, they are not validated, and they are not
required to reproduce any published figure or result table.
