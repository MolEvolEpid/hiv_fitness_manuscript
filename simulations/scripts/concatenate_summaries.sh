
cat output/*/simulations/*_tr_summary_stats.csv > summaries/tr_summary_stats.csv
find . -type f -name '*_counts.csv' -print0 | xargs -0 -I {} bash -c 'awk -v filename="{}" "{print filename \",\" \$0}" {}' > summaries/counts.csv
find . -type f -name '*_fitness.csv' -print0 | xargs -0 -I {} bash -c 'awk -v filename="{}" "{print filename \",\" \$0}" {}' > summaries/fitness.csv
