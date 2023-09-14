#!/bin/bash

filename=${1:0: -4}
dssp -i $filename".pdb" -o $filename".dssp"
gawk 'NR>28' $filename".dssp" > $filename"_dssp.txt"
#sed -i '1,28d' $filename".dssp"
#sed -i '/!/d' $filename".dssp"

cat $filename"_dssp.txt" |while read line;do
  v=($line)
  sschar=${v[4]}
  nline_char=${v[0]}
  nres_char=${v[1]}
  chain_char=${v[2]}
  res_char=${v[3]}
  if [ "$sschar" != H ] && [ "$sschar" != B ] && [ "$sschar" != G ] && [ "$sschar" != E ] && [ "$sschar" != I ] && [ "$sschar" != T ] && [ "$sschar" != S ]; then
    sschar=-
  fi
  echo -e $nline_char $nres_char $chain_char $res_char $sschar "\n" >> $filename"_dssp.tmp"
done

sed -i '/^$/d' $filename"_dssp.tmp"
mv $filename"_dssp.tmp" $filename"_dssp.txt"

#rm aKNL2_ranked_0.dssp



#loop pra cada frame
#   abrir um frame do dcd com o vmd
#   criar arquivo pdb temporário
#   rodar o dssp
#   loop pra cada resíduo
#      ler coluna 14 desse resíduo
#      trocar letra por número
#      escrever um número numa linha em novo arquivo
#   end loop
#   excluir arquivos dssp e pdb temporários
#end loop




#inicia o vmd com o arquivo de estrutura .psf. Escreva isso:
#vmd -dispdev text -psf ../1fsf_sol.psf -e selectPRO.tcl
#set k 32
# loop para carregar, ler, e reescrever cada dcd
#for {set l 28} {$l < 32} {incr l 1} {
#loop pra cada 400 frames
#for {set i 0} {$i < 4} {incr i 1} {
#set fs [expr {$i*4000}]
#set fe [expr {$i*4000 + 3999}]
#set n [expr {$i + 46}]
#animate read dcd ../1fsf_dyn.31.dcd beg $fs end $fe waitfor all top
#set sel [atomselect top "protein"]
#animate write dcd 1fsf_dyn.$n.pro.dcd waitfor all sel $sel
#$sel delete
#animate delete all
#}
#set k n
#}

#exit
