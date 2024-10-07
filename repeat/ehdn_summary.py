import argparse
import pandas as pd


def select_one_family(df,proband,others,controls):
  members = [proband]+others
  df = df[df['counts'].str.contains(proband)]
  counts4member_d = dict(zip(members,[[] for i in range(len(members))]))
  mafAll_l = []
  mafCtl_l = []
  for idx in range(len(df)):
    count4member_d = dict(zip(members,[float(0)]*len(members)))
    count4all = df.iloc[idx]['counts']
    for count4each in count4all.split(','):
      SAMPLE = count4each.split(':')[0]
      COUNT  = float(count4each.split(':')[1])
      if SAMPLE in members:
        count4member_d[SAMPLE] = COUNT
    for member in members:
      counts4member_d[member] = counts4member_d[member]+[count4member_d[member]]
    count4proband = count4member_d[proband]
    mafAll = len([i for i in count4all.split(',') if float(i.split(':')[1]) > count4proband])
    mafCtl = len([i for i in count4all.split(',') if float(i.split(':')[1]) > count4proband and i.split(':')[0] in controls])
    mafAll_l = mafAll_l + [mafAll]
    mafCtl_l = mafCtl_l + [mafCtl]
  df['countInAll'] = mafAll_l
  df['countInCtl'] = mafCtl_l
  for member in members:
    df[member] = counts4member_d[member]
  return df.drop('top_case_zscore',axis=1)


def rename_motif_proband(motif_proband,proband,others):
  members = [proband]+others
  motif_proband = motif_proband.rename({'high_case_counts':'high_case_counts_motif','counts':'counts_motif','countInAll':'countInAll_motif','countInCtl':'countInCtl_motif'},axis=1)
  for member in members:
    motif_proband = motif_proband.rename({member:member+'_motif'},axis=1)
  return motif_proband


def run():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument("-proband", type=str, help="") 
  parser.add_argument("-others", type = str, help="")
  parser.add_argument("-controls", type = str, help="")
  parser.add_argument("-locus_summary", type = str, help="")
  parser.add_argument("-motif_summary", type = str, help="")
  parser.add_argument("-out", type = str, help="")
  args = parser.parse_args()
  proband = args.proband
  out = args.out
  with open(args.others,'r') as f:
    others__l = [i.strip() for i in f.readlines()]    
  with open(args.controls,'r') as f:
    controls__l = [i.strip() for i in f.readlines()]    
  locus = pd.read_table(args.locus_summary)
  motif = pd.read_table(args.motif_summary)
  locus_proband = select_one_family(locus,proband,others__l,controls__l)
  motif_proband = select_one_family(motif,proband,others__l,controls__l)
  motif_proband = rename_motif_proband(motif_proband,proband,others__l)
  pd.merge(locus_proband,motif_proband,how='outer').to_csv(out,sep='\t',index=False)


if __name__ == "__main__":
  run()

