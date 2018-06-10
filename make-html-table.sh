#!/bin/bash -e

cat <<EOF
<div style="font-family: courier;">
  <table class="rounded table">
    <tr>
      <th>blastn</th>
      <th>bwa mem</th>
      <th>bwa aln</th>
      <th>bwa aln -l</th>
      <th>bowtie2</th>
    </tr>
    <tr>
EOF

for alg in blastn bwa-mem bwa-aln bwa-aln-l bowtie2
do
    echo "Running $alg" >&2
    echo '    <td>'
    echo -n '      <pre>'
    script=match-$alg.sh
    $script "$@"
    echo '</pre>'
    echo '    </td>'
done

echo '    </tr>'
echo '  </table>'
echo '</div>'
