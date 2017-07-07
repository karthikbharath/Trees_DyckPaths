#!/usr/bin/env bash
# Author Anju Kambadur
# Please change:
# --dot-in-dir <where-you-have-your-files>
# Currently it's set to my directory:
# --dot-in-dir /Users/pkambadu/Downloads/${method}

if [ ! -d results ]; then
  mkdir results
fi

for method in Mean Median Centroid Complete; do
  for smpl_rate in 0.1 0.2 0.3 0.4 0.5; do
    cat <<EOF > options.txt;
--lca-seed 1234
--verbosity 0
--use-random-weights false
--measure-lca true
--percent-lca-samples ${smpl_rate}
--gen-method Graphviz
--dot-in-dir /Users/pkambadu/Downloads/${method}
--sample-leaves true
EOF
  echo "Running ${method} with ${smpl_rate} samples"
  ./harness --response-file options.txt \
            > results/result.${method}.${smpl_rate}.txt \
            2> results/result.${method}.${smpl_rate}.err
  done
done

if [ -f options.txt ]; then
  rm options.txt
fi
