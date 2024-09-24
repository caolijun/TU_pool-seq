import sys
with open(sys.argv[1]) as rc:
    num_of_pop = int((len(next(rc).split())-9)/2)
    for line in rc:
        tabs = line.strip().split()
        if tabs[3] == "2":
            maj_alleles = tabs[7]
            min_alleles = tabs[8]
            maj_allele, min_allele = tabs[4].split("/")
            pos_str = []
            for i in range(0, len(maj_alleles)):
                # if (maj_alleles[i] == maj_allele or maj_alleles[i] == min_allele) and min_alleles[i] == "N":
                if maj_alleles[i] == min_allele:
                    pos_str.append(tabs[9+i+num_of_pop].split("/")[0] + "," + tabs[9+i].split("/")[0])
                elif maj_alleles[i] == maj_allele:
                    pos_str.append(tabs[9+i].split("/")[0] + "," + tabs[9+i+num_of_pop].split("/")[0])
                #
                # elif maj_alleles[i] == maj_allele and min_alleles[i] == min_allele:
                #     pos_str.append(tabs[9+i].split("/")[0] + "," + tabs[9+i+num_of_pop].split("/")[0])
                # elif maj_alleles[i] == min_allele and min_alleles[i] == "N":
                #     pos_str.append(tabs[9+i].split("/")[0] + "," + tabs[9+i+num_of_pop].split("/")[0])
                elif maj_alleles[i] == "N" and min_alleles[i] == "N":
                    pos_str.append(tabs[9+i+num_of_pop].split("/")[0] + "," + tabs[9+i].split("/")[0])
            if len(pos_str) == num_of_pop:
                print(tabs[0] + "\t" + tabs[1] + "\t" + "\t".join(pos_str))
