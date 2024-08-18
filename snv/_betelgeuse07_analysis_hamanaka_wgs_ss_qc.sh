#(venv) [hamanaka@canopus wgs]$ ls /betelgeuse07/analysis/ncgm/*/*/DA0*/QCmetrics/DA0*.WgsMetrics.chrY.txt > tmp.y.txt
#(venv) [hamanaka@canopus wgs]$ less -SN tmp.y.2.txt 
#(venv) [hamanaka@canopus wgs]$ head -n1 tmp.y.txt
#/betelgeuse07/analysis/ncgm/data01/G011/DA0000000348/QCmetrics/DA0000000348.WgsMetrics.chrY.txt
#(venv) [hamanaka@canopus wgs]$ less -SN /betelgeuse07/analysis/ncgm/data01/G011/DA0000000348/QCmetrics/DA0000000348.WgsMetrics.chrY.txt
#[1]+  終了                  less tmp.y.txt | xargs -I{} sh -c 'echo -n {}; head -n3 {}|tail -n1' > tmp.y.2.txt
#(venv) [hamanaka@canopus wgs]$ sed -n2,2p < /betelgeuse07/analysis/ncgm/data01/G011/DA0000000348/QCmetrics/DA0000000348.WgsMetrics.chrY.txt|less -SN
#(venv) [hamanaka@canopus wgs]$ sed -p2,2n < /betelgeuse07/analysis/ncgm/data01/G011/DA0000000348/QCmetrics/DA0000000348.WgsMetrics.chrY.txt|less -SN
#(venv) [hamanaka@canopus wgs]$ sed -n 2p < /betelgeuse07/analysis/ncgm/data01/G011/DA0000000348/QCmetrics/DA0000000348.WgsMetrics.chrY.txt|less -SN


#ls /betelgeuse07/analysis/ncgm/yokohama-005/DA0*/QCmetrics/DA0*.WgsMetrics.autosome.txt > tmp.a.txt
#ls /betelgeuse07/analysis/ncgm/yokohama-005/DA0*/QCmetrics/DA0*.WgsMetrics.chrX.txt > tmp.x.txt
#ls /betelgeuse07/analysis/ncgm/yokohama-005/DA0*/QCmetrics/DA0*.WgsMetrics.chrY.txt > tmp.y.txt
find /mira03/NCGM -name 'DA0*.WgsMetrics.autosome.txt' > tmp.a.txt
find /mira03/NCGM -name 'DA0*.WgsMetrics.chrX.txt'     > tmp.x.txt
find /mira03/NCGM -name 'DA0*.WgsMetrics.chrY.txt'     > tmp.y.txt
cat <(sed -n 2p < /betelgeuse07/analysis/ncgm/data01/G011/DA0000000348/QCmetrics/DA0000000348.WgsMetrics.autosome.txt) <(less tmp.a.txt | xargs -I{} sh -c 'echo -n {}; head -n3 {}|tail -n1') > tmp.a2.txt
cat <(sed -n 2p < /betelgeuse07/analysis/ncgm/data01/G011/DA0000000348/QCmetrics/DA0000000348.WgsMetrics.chrX.txt)     <(less tmp.x.txt | xargs -I{} sh -c 'echo -n {}; head -n3 {}|tail -n1') > tmp.x2.txt
cat <(sed -n 2p < /betelgeuse07/analysis/ncgm/data01/G011/DA0000000348/QCmetrics/DA0000000348.WgsMetrics.chrY.txt)     <(less tmp.y.txt | xargs -I{} sh -c 'echo -n {}; head -n3 {}|tail -n1') > tmp.y2.txt
paste <(cut -f1,2 tmp.a2.txt) <(cut -f2 tmp.x2.txt) <(cut -f2 tmp.y2.txt)|awk 'NR==1{print $0}NR>1{print $1,$3/$2,$4/$2}'|sort -k2,2n > tmp.res.txt


