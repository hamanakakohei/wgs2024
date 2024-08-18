# denovo file wo kireinisuru
import pandas as pd
DENOVO = 'denovoqc.anno.all.txt'
OUT = 'denovoqc.anno.all.filtered.txt'

dt = pd.read_table(DENOVO, low_memory = False)
dt['gt_gt']              = dt['gt'   ].map(lambda x: x.split(':')[0])
dt['gt_ad']              = dt['gt'   ].map(lambda x: int(x.split(':')[1].split(',')[1]))
dt['gt_dp']              = dt['gt'   ].map(lambda x: int(x.split(':')[2]))
dt['gt_allele_ratio']    = dt['gt_ad'].astype('int')    / dt['gt_dp'   ].astype('int')
dt['gt.fa_gt']           = dt['gt.fa'].map(lambda x: x.split(':')[0])
dt['gt.fa_ad']           = dt['gt.fa'].map(lambda x: int(x.split(':')[1].split(',')[1]))
dt['gt.fa_dp']           = dt['gt.fa'].map(lambda x: int(x.split(':')[2]))
dt['gt.fa_allele_ratio'] = dt['gt.fa_ad'].astype('int') / dt['gt.fa_dp'].astype('int')
dt['gt.mo_gt']           = dt['gt.mo'].map(lambda x: x.split(':')[0])
dt['gt.mo_ad']           = dt['gt.mo'].map(lambda x: int(x.split(':')[1].split(',')[1]))
dt['gt.mo_dp']           = dt['gt.mo'].map(lambda x: int(x.split(':')[2]))
dt['gt.mo_allele_ratio'] = dt['gt.mo_ad'].astype('int') / dt['gt.mo_dp'].astype('int')

dt.\
  query('gt_dp > 9').\
  query('`gt.fa_dp` > 9').\
  query('`gt.mo_dp` > 9').\
  query('gt_ad > 4').\
  query('gt_allele_ratio > 0.2').\
  query('`gt.fa_ad` == 0 or `gt.mo_ad` == 0').\
  to_csv(OUT, sep = '\t', index = False)
  #query('FILTER == "PASS"').\
  #query('`gt.fa_ad` > 9').\
  #query('`gt.mo_ad` > 9').\
  #query('gt.fa_gt == "0/0"').\
  #query('gt.mo_gt == "0/0"').\
  #query('gt_gt    == "0/1" or gt_gt == "1/0" or gt_gt == "1/1"').\


# each family
DENOVOALL1=denovoqc.anno.all.txt
DENOVOALL2=denovoqc.anno.all.filtered.txt

PED=/betelgeuse07/analysis/hamanaka/wgs/sample1211.plus40trio.ped
less $PED | awk '$NF==2{print $2}' | grep -e 00603 -e 00640 |  while read PATIENT; do
  awk -v PATIENT="${PATIENT}" 'NR==1 || $3==PATIENT' $DENOVOALL1 > denovoqc.anno.all.$PATIENT.txt
  awk -v PATIENT="${PATIENT}" 'NR==1 || $3==PATIENT' $DENOVOALL2 > denovoqc.anno.all.filtered.$PATIENT.txt
done
