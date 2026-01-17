#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l mem=120g
#PBS -l walltime=720:00:00
#PBS -q large
#PBS -V
#PBS -o myout.txt
#PBS -e myerr.txt

cd $PBS_O_WORKDIR
source /project/yangy/packages/mambaforge/bin/activate r_env


seeds=$(seq 1 100)
base_s_ceof=0
base_death_rate=0.3
base_adv_rate=0
base_push_prob=1
base_mut_rate=0.6
base_sample_diameter=16
base_vaf_cutoff=0

# 创建 output 目录
mkdir -p output
mkdir -p csv 


# 并行运行
export R_SCRIPT="Simdata_chenbz_yy.R"


#Rscript $R_SCRIPT ${base_s_ceof} ${base_death_rate} ${base_adv_rate} ${base_push_prob} ${base_mut_rate} ${base_sample_diameter} ${base_vaf_cutoff} ${seeds[@]}


# #基准组合
#parallel --jobs 20 Rscript $R_SCRIPT {1} {2} {3} {4} {5} {6} {7} {8} ::: ${base_s_ceof} ::: ${base_death_rate} ::: ${base_adv_rate} ::: ${base_push_prob} ::: ${base_mut_rate} ::: ${base_sample_diameter} ::: ${base_vaf_cutoff} ::: ${seeds[@]}


# #自然选择
# s_ceof=$(seq 0 0.2 0.6)
# s_ceof=("0.4" "0.6")
# adv_rate=0.000001
# parallel --jobs 20 Rscript $R_SCRIPT {1} {2} {3} {4} {5} {6} {7} {8} ::: ${s_ceof[@]} ::: ${base_death_rate} ::: ${adv_rate} ::: ${base_push_prob} ::: ${base_mut_rate} ::: ${base_sample_diameter} ::: ${base_vaf_cutoff} ::: ${seeds[@]}

# #克隆比例
# s_ceof=0.4
# adv_rate=$(awk 'BEGIN { for (i = 0.00001; i >= 0.0000001; i /= 10) print i}')
# #adv_rate=0.000001
# parallel --jobs 20 Rscript $R_SCRIPT {1} {2} {3} {4} {5} {6} {7} {8} ::: ${s_ceof} ::: ${base_death_rate} ::: ${adv_rate[@]} ::: ${base_push_prob} ::: ${base_mut_rate} ::: ${base_sample_diameter} ::: ${base_vaf_cutoff} ::: ${seeds[@]}

# #死亡率
# death_rate=$(seq 0.2 0.1 0.4)
# parallel --jobs 20 Rscript $R_SCRIPT {1} {2} {3} {4} {5} {6} {7} {8} ::: ${base_s_ceof} ::: ${death_rate[@]} ::: ${base_adv_rate} ::: ${base_push_prob} ::: ${base_mut_rate} ::: ${base_sample_diameter} ::: ${base_vaf_cutoff} ::: ${seeds[@]}

# #push比例
# push_prob=$(seq 0 0.2 0.8)
push_prob="0.5"
parallel --jobs 20 Rscript $R_SCRIPT {1} {2} {3} {4} {5} {6} {7} {8} ::: ${base_s_ceof} ::: ${base_death_rate} ::: ${base_adv_rate} ::: ${push_prob[@]} ::: ${base_mut_rate} ::: ${base_sample_diameter} ::: ${base_vaf_cutoff} ::: ${seeds[@]}

# # # # #mutation rate
# mut_rate=("0.3" "1.2")
# parallel --jobs 20 Rscript $R_SCRIPT {1} {2} {3} {4} {5} {6} {7} ::: ${base_s_ceof} ::: ${base_death_rate} ::: ${base_adv_rate} ::: ${base_push_prob} ::: ${mut_rate[@]} ::: ${base_sample_diameter} ::: ${seeds[@]}

# #sample diameter
# sample_diameter=$(seq 8 20 80)
# parallel --jobs 20 Rscript $R_SCRIPT {1} {2} {3} {4} {5} {6} {7} {8} ::: ${base_s_ceof} ::: ${base_death_rate} ::: ${base_adv_rate} ::: ${base_push_prob} ::: ${base_mut_rate} ::: ${sample_diameter[@]} ::: ${base_vaf_cutoff} ::: ${seeds[@]}

# #vaf
# base_vaf_cutoff=$(seq 0.01 0.02 0.1)
# parallel --jobs 20 Rscript $R_SCRIPT {1} {2} {3} {4} {5} {6} {7} {8} ::: ${base_s_ceof} ::: ${base_death_rate} ::: ${base_adv_rate} ::: ${base_push_prob} ::: ${base_mut_rate} ::: ${base_sample_diameter} ::: ${base_vaf_cutoff[@]} ::: ${seeds[@]}

# mut_rate
# mut_rate=("0.6" "1.2")
# parallel --jobs 20 Rscript $R_SCRIPT {1} {2} {3} {4} {5} {6} {7} {8} ::: ${base_s_ceof} ::: ${base_death_rate} ::: ${base_adv_rate} ::: ${base_push_prob} ::: ${mut_rate[@]} ::: ${base_sample_diameter} ::: ${base_vaf_cutoff[@]} ::: ${seeds[@]}


# #自然选择2
# s_ceof=$(seq 0.2 0.2 0.8)
# adv_rate=0.000001
# parallel --jobs 20 Rscript $R_SCRIPT {1} {2} {3} {4} {5} {6} {7} {8} ::: ${s_ceof[@]} ::: ${base_death_rate} ::: ${adv_rate} ::: 0 ::: ${base_mut_rate} ::: ${base_sample_diameter} ::: ${base_vaf_cutoff} ::: ${seeds[@]}

