import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import csv
import matplotlib.patches as mpatches
from matplotlib.colors import rgb2hex
import numpy as np



def prepare_for_plots(final_serotype_calls):

    virus_dict = {}
    all_viruses = set()
    with open(final_serotype_calls) as f:
        data = csv.DictReader(f, delimiter="\t")
        for l in data:
            virus_dict[l['sample_id']] = l['serotype_called']
            all_viruses.add(l['serotype_called'])

    colour_dict = {}
    custom_cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["#567CBE","#D58A80", "#ADB3D9"], len(all_viruses))

    lst = custom_cmap(range(len(all_viruses)))
    for i,j in enumerate(all_viruses):
        colour_dict[j] = rgb2hex(lst[i])

    patch_list = []
    for serotype,hexc in colour_dict.items():
        patch_list.append(mpatches.Patch(color=hexc, label=serotype))
 
    return virus_dict, colour_dict, patch_list

def variant_plot(results_dir, variants_summary_file, virus_dict, colour_dict, patch_list):
    
    variant_num = {}
    with open(variants_summary_file) as f:
        data = csv.DictReader(f, delimiter="\t")
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
        colours.append(colour_dict[virus_dict[sample]])

    plt.xticks(rotation=90, size=15)
    plt.yticks(size=15)
    plt.scatter(x,y,color=colours,s=70)

    plt.xlabel("Sample ID", size=20)
    plt.ylabel("Variant number", size=20)

    plt.legend(handles=patch_list,fontsize=15,frameon=False)

    plt.savefig(os.path.join(results_dir, "variant_plot.pdf"), bbox_inches="tight")


def ct_plot(results_dir, ct_file, ct_column, id_column, final_serotype_calls, virus_dict, colour_dict, patch_list):

    ct_dict = {}

    with open(ct_file) as f:
        data = csv.DictReader(f)
        for l in data:
            if l[id_column] in virus_dict:
                try:
                    value = float(l[ct_column])
                    if np.isnan(value):
                        value = 45
                except ValueError:
                    value = 45
                    
                ct_dict[l[id_column]] = value
                
                
    coverage_dict = {}   
    with open(final_serotype_calls) as f:
        data = csv.DictReader(f, delimiter="\t")
        for l in data:
            if l['coverage_trimmed'] != "0" and l['coverage_trimmed'] != "NA":
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
        colours.append(colour_dict[virus_dict[sample]])
        
    plt.xticks(rotation=90, size=15)
    plt.yticks(size=15)
    plt.scatter(x,y,color=colours,s=70)

    plt.xlabel("Ct value", size=20)
    plt.ylabel("Coverage", size=20)

    plt.legend(handles=patch_list,fontsize=15,frameon=False)

    plt.savefig(os.path.join(results_dir, "ct_plot.pdf"), bbox_inches="tight")

