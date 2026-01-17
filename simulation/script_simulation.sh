sim_dim=2
deme_size=100
mut_rate=0.6
adv_rate="0.1"
s_coef=0.4
repl=1
path="./csv"
rd=10
birth_rate=0.4
death_rate=0.3
push_prop=0.1
mig_rate=0
punch_diameter=3
puch_density=0.1
spacing=3
title="example"
snapname="snap_example"
python 3DTumorSimulPush.py \
  $sim_dim \
  $deme_size \
  $mut_rate \
  $adv_rate \
  $s_coef \
  $repl \
  $path \
  $rd \
  $birth_rate \
  $death_rate \
  $push_prop \
  $mig_rate \
  $punch_diameter \
  $puch_density \
  $spacing \
  $title \
  $snapname
