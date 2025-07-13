# OGtoberfest
Benchmarking Orthogroups and Orthologs

## Preprocessing

Example commands:
- **OrthoBench**
```bash
    ogtoberfest orthogroups_preprocess -i OrthoBench/predicted_orthogroups --database OrthoBench/proteomes -t 8
```

- **Simulations**
```bash
    ogtoberfest orthogroups_preprocess -i Sim1k/predicted_orthogroups --database Sim1k/proteome_files --use-id -t 8
```

## Benchmarking
### Orthogroups Benchmarking

Example commands:
- **OrthoBench**
```bash
    ogtoberfest orthogroups_benchmark \
        -i OrthoBench/preprocessed_predicted_orthogroups \
        --database OrthoBench/proteomes \
        -t 8
```

- **Simulations**
```bash
    ogtoberfest orthogroups_benchmark \
    -i Sim1k/preprocessed_predicted_orthogroups \
    --database Sim1k/proteome_files \
    --refog Sim1k/reference_orthogroups.txt \
    --refog-stats Sim1k/reference_orthogroups_stats.txt \
    --use-id \
    -t 8
```
