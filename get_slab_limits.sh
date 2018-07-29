more ../openquake-inslab-builder/slab_demo_data/sum_contours_clean_remove_duplicates_minlon_102E_maxlon_114E.in | grep 300 > java_slab_source_model300.xml.limits.txt
more ../openquake-inslab-builder/slab_demo_data/sum_contours_clean_remove_duplicates_minlon_102E_maxlon_114E.in | grep -P '\-20\t' > java_slab_source_model20.xml.limits.txt
sort -t, -nk1 java_slab_source_model300.xml.limits.txt > java_slab_source_model300sort.xml.limits.txt
sort -t, -k1,1rn java_slab_source_model20.xml.limits.txt > java_slab_source_model20sort.xml.limits.txt
cat java_slab_source_model20sort.xml.limits.txt java_slab_source_model300sort.xml.limits.txt > data/java_slab_source_model.xml.limits.txt
tr -s '\t' ',' < data/java_slab_source_model.xml.limits.txt > data/java_slab_source_model.xml.limits.txt.tmp
mv data/java_slab_source_model.xml.limits.txt.tmp data/java_slab_source_model.xml.limits.txt
