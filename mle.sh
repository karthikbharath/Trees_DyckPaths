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
  cat <<EOF > options.txt;
--mle-seed 1234
--verbosity 1
--use-random-weights false
--measure-mle true
--gen-method Graphviz
--dot-in-dir ~/Desktop/imagedata_results/${method}
EOF
  echo "Running ${method}"
  ./harness --response-file options.txt \
            > results/result.${method}.txt 
done

if [ -f options.txt ]; then
  rm options.txt
fi
