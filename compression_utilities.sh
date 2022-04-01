# to extract subtranches
for i in *.tar; do tar -xvf "$i"; done

# to extract the molecules
find . -name '*.tar.gz' -execdir tar -xzvf '{}' \;

# to compute number of the molecules in a tranch
find . -name '*.pdbqt' | wc -l;