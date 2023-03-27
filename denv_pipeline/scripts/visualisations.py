import os
import matplotlib.pyplot as plt
import csv
import matplotlib.patches as mpatches
import numpy as np



def prepare_for_plots(final_serotype_calls):

    serotype_dict = {}
    with open(final_serotype_calls) as f:
        data = csv.DictReader(f, delimiter="\t")
        for l in data:
            serotype_dict[l['sample_id']] = l['serotype']

    colour_dict = {"DENV1": "#2E9680",
                "DENV2": "#E1A900",
                "DENV3":"#DC810B",
                "DENV4":"#75B7D5"}

    patch_list = []
    for serotype,hexc in colour_dict.items():
        patch_list.append(mpatches.Patch(color=hexc, label=serotype))


    return serotype_dict, colour_dict, patch_list

def variant_plot(results_dir, variants_summary_file, serotype_dict, colour_dict, patch_list):
    
    variant_num = {}
    with open(variants_summary_file) as f:
        data = csv.DictReader(delimiter="\t")
        for l in data:
            variant_num[l['sample_id']] = float(l['variant_count'])

    variant_num = {k:v for k,v in sorted(variant_num.items(), key=lambda x: x[1], reverse=True)}

    fig, ax = plt.subplots(1,1, figsize=(20,10))

    x = []
    y = []
    colours = []

    for sample, variants in variant_num.items():
        x.append(sample)
        y.append(variants)
        colours.append(colour_dict[serotype_dict[sample]])

    plt.xticks(rotation=90, size=15)
    plt.yticks(size=15)
    plt.scatter(x,y,color=colours,s=70)

    plt.xlabel("Sample ID", size=20)
    plt.ylabel("Variant number", size=20)

    plt.legend(handles=patch_list,fontsize=15,frameon=False)

    plt.savefig(os.path.join(results_dir, "variant_plot.pdf"), bbox_inches="tight")


def ct_plot(results_dir, ct_file, ct_column, id_column, final_serotype_calls, serotype_dict, colour_dict, patch_list):

    ct_dict = {}

    with open(ct_file) as f:
        data = csv.DictReader(f)
        for l in data:
            if l[id_column] in serotype_dict:
                if type(l[ct_column]) == float and not np.isnan(l[ct_column]):
                    ct_dict[l[id_column]] = float(l[ct_column])
                else:
                    ct_dict[l[id_column]] = 45
                
                
    coverage_dict = {}   
    with open(final_serotype_calls) as f:
        data = csv.DictReader(f, delimiter="\t")
        for l in data:
            if l['coverage_trimmed'] != "NA":
                cov = l['coverage_trimmed']
            else:
                cov = l['coverage_untrimmed']
            coverage_dict[l['sample_id']] = float(cov)


    fig, ax = plt.subplots(1,1, figsize=(20,10))

    x = []
    y = []
    colours = []

    for sample, ct in ct_dict.items():
        x.append(ct)
        y.append(coverage_dict[sample])
        colours.append(colour_dict[serotype_dict[sample]])
        
    plt.xticks(rotation=90, size=15)
    plt.yticks(size=15)
    plt.scatter(x,y,color=colours,s=70)

    plt.xlabel("Ct value", size=20)
    plt.ylabel("Coverage", size=20)

    plt.legend(handles=patch_list,fontsize=15,frameon=False)

    plt.savefig(os.path.join(results_dir, "ct_plot.pdf"), bbox_inches="tight")
