import os
from collections import defaultdict  
import csv

"""
call loci are the regions generated from looking at the underlying read depth so as to call snps under reasonablye
even size compute chunks. here size means the coverage, and the compute time

assess loci are simpler, they are just chunks over which we want to compute stats, they could, for instance be overlapping 
"""

call_loci = [l.rstrip().replace(":","_") for l in open("./get_regions/input_loci_MINUS_SATELLITE.txt").readlines()]

assess_loci = []

chunk_size, chunk_slide = 1000000, 1000000
for l in open("/net/eichler/vol7/home/psudmant/genomes/fastas/human_1kg_v37/human_1kg_v37.fasta.fai"):
    contig, length = l.split("\t", 3)[:2]
    if not "MT" in contig and not "GL" in contig:
        length = int(length)
        s = 1
        e = min(s+chunk_size, length)
        while s<e:
            assess_loci.append("%s_%d-%d"%(contig, s, e))
            s=e
            e = min(s+chunk_size, length)

rule all:
    input: "merged_SNPs/C_TEAM.vcf.gz", "merged_SNPs/C_TEAM.vcf.gz.tbi", "./analyses/snp_counts.summary"
    params: sge_opts="-l mfree=2G -N run_all"

rule get_stats:
    input: expand("per_locus_analyses/{assess_locus}.txt", assess_locus=assess_loci)
    params: sge_opts="-l mfree=2G -N index_VCF"
    output: "./analyses/snp_counts.summary","./analyses/snp_counts.along_chr.dataframe"
    run:
        FOUT_sum = open(output[0],'w')
        FOUT_chrsum = open(output[1],'w')
        
        FOUT_sum.write("indiv\tcontinental_population\tpopulation\thets\thoms\n")
        FOUT_chrsum.write("indiv\tcontinental_population\tpopulation\tcontig\tstart\tend\thets\thoms\n")

        sum_inf_by_indiv = defaultdict(lambda: defaultdict(int)) 
            
        for f in input:
            if "X" in f or "Y" in f: continue
            contig, s_e = f.split("/")[-1].split(".")[0].split("_") 
            s,e = s_e.split("-")
            
            inf_by_indiv = defaultdict(lambda: defaultdict(int))
            for ldict in csv.DictReader(open(f),delimiter="\t"):
                t=ldict['TYPE'] 
                for k,v in ldict.items():
                    if k!="TYPE":
                        inf_by_indiv[k][t]=int(v)
            for indiv, inf in inf_by_indiv.items():
                cont_pop, pop = indiv.split("_",3)[:2]
                FOUT_chrsum.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n"%(indiv,cont_pop, pop, contig, s,e, inf["HET"], inf["HOM"]))
                sum_inf_by_indiv[indiv]['HOM'] += inf['HOM']
                sum_inf_by_indiv[indiv]['HET'] += inf['HET']
        
        for indiv, inf in sum_inf_by_indiv.items():
            cont_pop, pop = indiv.split("_",3)[:2]
            FOUT_sum.write("%s\t%s\t%s\t%d\t%d\n"%(indiv,cont_pop, pop, inf["HET"], inf["HOM"]))

        FOUT_sum.close()
        FOUT_chrsum.close()
rule compute_stat:
    input: "merged_SNPs/C_TEAM.vcf.gz", "merged_SNPs/C_TEAM.vcf.gz.tbi" 
    output: "per_locus_analyses/{assess_locus}.txt"
    params: loc = "{assess_locus}", sge_opts="-l mfree=4G -N assess_{assess_locus}"
    run:
        loc = params.loc.replace("_",":")
        shell("""vcffilter -f "DP > 10" -f "QUAL > 20" -f "TYPE = snp" -r %s {input[0]} | vcfhethomcount - > {output[0]}"""%loc)

rule index_VCF:
    input: "merged_SNPs/C_TEAM.vcf.gz"
    output: "merged_SNPs/C_TEAM.vcf.gz.tbi"
    params: sge_opts="-l mfree=2G -N index_VCF"
    shell: "tabix -p vcf {input[0]}"

rule merge_VCFs:
    input: expand("C_team_SNPs/chunked/{locus}.vcf",locus=call_loci)
    output: "merged_SNPs/C_TEAM.vcf.gz"
    params: sge_opts="-l mfree=10G -N merge_VCF"
    run: 
        a_chunk = "C_team_SNPs/chunked/%s.vcf"%(call_loci[0])
        chunks = "\n".join(["C_team_SNPs/chunked/%s.vcf"%(vcf) for vcf in call_loci])
        F = open("./tmp_files",'w')
        F.write(chunks)
        F.close()
        shell("""((head -100 %s | grep "^#"; cat ./tmp_files | while read f; do echo $f >&2; cat $f | egrep -v "^#"; done) | vcfstreamsort -w 10000 | vcfuniq ) | bgzip -c > {output[0]}"""%(a_chunk))
        os.unlink("tmp_files")

rule call_SNPs:
    output: "C_team_SNPs/chunked/{locus}.vcf"
    params: loc = "{locus}", sge_opts="-l mfree=40G -N loc_{locus}"
    run:
        loc = params.loc.replace("_",":")
        print(loc)
        print(params.loc)
        #shell("freebayes -C 5 --report-monomorphic -L ./bamlist.txt -f ~/genomes/fastas/hg19_1kg_phase2_reference/human_g1k_v37.fasta -r %s >{output[0]}"%loc)
        shell("freebayes -C 5 -L ./bamlist.txt -f ~/genomes/fastas/hg19_1kg_phase2_reference/human_g1k_v37.fasta -r %s >{output[0]}"%loc)
