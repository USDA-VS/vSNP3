#!/usr/bin/env python

__version__ = "3.15"

import os

class Step2_Summary:

    def html_step2_summary(self, summary_dictionary, working_directory='.'):

        st = summary_dictionary.get("st", "empty")
        starting_vcf_number = summary_dictionary.get("starting_vcf_number", "empty")
        fixed_vcf_number_removed = summary_dictionary.get("fixed_vcf_number_removed", "empty")
        reference = summary_dictionary.get("reference", "empty")
        remove_from_analysis_count = summary_dictionary.get("remove_from_analysis_count", "empty")
        all_vcf_boolen = summary_dictionary.get("all_vcf_boolen", "empty")
        subset_boolen = summary_dictionary.get("subset_boolen", "empty")
        keep_count = summary_dictionary.get("keep_count", "empty")
        samples_groups_dict = summary_dictionary.get("samples_groups_dict", {"No_groups": "No_groups"})
        samples_groups_dict = dict(sorted(samples_groups_dict.items()))
        malformed = summary_dictionary.get("malformed", "empty")
        starttime = summary_dictionary.get("starttime", "empty")
        endtime = summary_dictionary.get("endtime", "empty")
        runtime = summary_dictionary.get("runtime", "empty")

        htmlfile = open(f'{working_directory}/vSNP_step2_summary-{st}.html', 'at')
        
        #MAKE HTML FILE:
        print("<html>\n<head><style> table { font-family: arial, sans-serif; border-collapse: collapse; width: 40%; } td, th { border: 1px solid #dddddd; padding: 4px; text-align: left; font-size: 11px; } </style></head>\n<body style=\"font-size:12px;\">", file=htmlfile)
        print(f"<h2>Script ran using <u>{reference}</u> variables</h2>", file=htmlfile)
        print(f"<h4>{starting_vcf_number} VCF files in this run<br>", file=htmlfile)
        if subset_boolen:
            print(f"\n<h4>Subset ran, {keep_count:,} VCF files used in analysis. </h4>", file=htmlfile)
        print(f"{remove_from_analysis_count} VCF files in this run removed from analysis<br>", file=htmlfile)
        print(f"{fixed_vcf_number_removed} VCF files in this run corrupt and therefore removed</h4>", file=htmlfile)

        #OPTIONS
        if all_vcf_boolen:
            print("\n<h4>All_VCFs is available</h4>", file=htmlfile)

        #TIME
        print(f"\n<h4>Start time: {starttime} <br>", file=htmlfile)
        print(f"End time: {endtime} <br>", file=htmlfile)
        print(f"Total run time: {runtime}: </h4>", file=htmlfile)

        # ERROR LIST
        if len(malformed) < 1:
            print("<h2>No corrupt files found</h2>", file=htmlfile)
        else:
            print("\n<h2>Corrupt files removed</h2>", file=htmlfile)
            for i in malformed:
                print(f"{i} <br>", file=htmlfile)
            print("<br>", file=htmlfile)

        #GROUPING TABLE
        samples_count_listed_groups = len(samples_groups_dict)
        print(f'<h2>Groupings with {samples_count_listed_groups:,} listed:</h2>', file=htmlfile)
        print("<table>", file=htmlfile)
        print("<tr align=\"left\"><th>Sample Name</th><tr>", file=htmlfile)

        samples_count_listed_groups = len(samples_groups_dict)
        for key, value in samples_groups_dict.items():
            print("<tr>", file=htmlfile)
            print(f"<td>{key}</td>", end='\t', file=htmlfile)
            for group in value:
                print(f"<td>{group}</td>", end='\t', file=htmlfile)
            print("</tr>", file=htmlfile)
        print("</table>", file=htmlfile)

        try:
            print("\n<h2>Program versions:</h2>", file=htmlfile)
            versions = os.popen('conda list biopython | grep -v "^#"; \
            conda list numpy | egrep -v "^#|numpydoc"; \
            conda list pandas | grep -v "^#"; \
            conda list pysam | grep -v "^#"; \
            conda list raxml | grep -v "^#"').read()
            versions = versions.split('\n')
            for i in versions:
                print("%s<br>" % i, file=htmlfile)
        except:
            pass

        print("</body>\n</html>", file=htmlfile)
        htmlfile.close()