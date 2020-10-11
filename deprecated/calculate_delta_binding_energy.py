import argparse


class Variant:
    def __init__(self, info, chain_ids=None):
        items = info.split(',')
        self.energy_change = items[5]
        self.rmsd = items[13]
        
        if chain_ids:
            self.coevolution = items[0]
            self.conservative = items[1]
            self.count = items[2]
            self.date_first = items[3]
            self.date_last = items[4]
            self.id_1 = items[6]
            self.id_2 = items[7]
            self.location_1 = items[10]
            self.location_2 = items[11]
            for chain_id in chain_ids:
                exec('self.monomer_' + chain_id + '=list()')
            if len(items[14]) > 0:
                for mutation in items[14].split(';'):
                    for chain_id in chain_ids:
                        if mutation[0] == chain_id:
                            exec('self.monomer_' + chain_id + '.append(mutation)')
        else:
            self.is_in_PDB = items[8]
            self.layer = items[9]
            self.substitutions = list()
            if len(items[14]) > 0:
                for mutation in items[-2].split(';'):
                    self.substitutions.append(mutation[1:])
            self.tag = items[15]
    
    def match(self, monomer_dict, chain_id):
        mutations_str = eval("';'.join(sorted(self.monomer_" + chain_id + "))")
        monomer = monomer_dict.get(mutations_str)
        if monomer:
            exec('self.monomer_' + chain_id + '=monomer')


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--complex', type=str, required=True)
    parser.add_argument('-m', '--monomers', type=str, nargs='*', required=True)
    parser.add_argument('-id', '--chain_ids', type=str, nargs='*', required=True)
    parser.add_argument('-o', '--output_name', type=str)
    return parser.parse_args()

def read_monomer_info_from_csv(monomer_csv_list, chain_id_list):
    chain_dict = dict()
    for i, monomer_csv in enumerate(monomer_csv_list):
        monomer_dict = dict()
        with open(monomer_csv, 'r') as p_mono_csv:
            not_first_line = False
            for line in p_mono_csv:
                if not_first_line: # Skip the first line of the .csv file.
                    substitutions = line.split(',')[-2].split(';')
                    key = ';'.join(sorted(substitutions))
                    monomer = Variant(line[:-1])
                    monomer_dict[key]=monomer
                not_first_line = True
        chain_dict[chain_id_list[i]] = monomer_dict
    return chain_dict

def write_summary_csv(chain_dict, cplx_csv, smry_csv):
    with open(smry_csv, 'w+') as p_smry_csv:
        p_smry_csv.write('coevolution,conservative,count,date_first,date_last,id_1,id_2,')
        p_smry_csv.write('location_1,location_2,multimer_name,multimer_energy,multimer_rmsd,')
        for chain_name in chain_dict.keys():
            p_smry_csv.write('monomer_' + chain_name + '_name,monomer_' + chain_name + \
                '_mutations,monomer_' + chain_name + '_is_in_PDB,monomer_' + chain_name + \
                '_layer,monomer_' + chain_name + '_energy,monomer_' + chain_name + '_rmsd,')
        p_smry_csv.write('delta_binding_energy\n')

    with open(cplx_csv, 'r') as p_cplx_csv:
        not_first_line = False
        for line in p_cplx_csv:
            if not_first_line: # Skip the first line of the .csv file.
                with open(smry_csv, 'a+') as p_smry_csv:
                    cplx = Variant(line[:-1], chain_dict.keys())
                    for chain_name, monomer_info_dict in chain_dict.items():
                        cplx.match(monomer_info_dict, chain_name)
                    # Write the multimer protein data to the new summary.csv file.
                    p_smry_csv.write(cplx.coevolution + ',' + cplx.conservative + ',' + \
                        cplx.count + ',' + cplx.date_first + ',' + cplx.date_last + ',' + \
                        cplx.id_1 + ',' + cplx.id_2 + ',' + cplx.location_1 + ',' + \
                        cplx.location_2 + ',')
                    # Concatenate monomer names to get the multimer name.
                    monomer_list = list()
                    for chain_name in chain_dict.keys():
                        monomer_list.append(eval('cplx.monomer_' + chain_name + '.tag'))
                    p_smry_csv.write('-'.join(sorted(monomer_list)) + ',' + \
                        cplx.energy_change + ',' + cplx.rmsd + ',')
                    # Write monomers protein data to the new summary.csv file.
                    if cplx.energy_change == 'NA':
                        delta_binding_energy = 0
                    else:
                        delta_binding_energy = float(cplx.energy_change)
                    for chain_name in chain_dict.keys():
                        p_smry_csv.write(eval('cplx.monomer_' + chain_name + '.tag') + ',')
                        
                        mutation_list = list()
                        for mutation in eval('cplx.monomer_' + chain_name + '.substitutions'):
                            mutation_list.append(mutation)
                        p_smry_csv.write(';'.join(mutation_list) + ',')

                        p_smry_csv.write(eval('cplx.monomer_' + chain_name + '.is_in_PDB') + ',')
                        p_smry_csv.write(eval('cplx.monomer_' + chain_name + '.layer') + ',')

                        monomer_energy = eval('cplx.monomer_' + chain_name + '.energy_change')
                        p_smry_csv.write(monomer_energy + ',')
                        if monomer_energy == 'NA':
                            monomer_energy = 0
                        else:
                            monomer_energy = float(monomer_energy)
                        delta_binding_energy -= monomer_energy

                        p_smry_csv.write(eval('cplx.monomer_' + chain_name + '.rmsd') + ',')
                    p_smry_csv.write(str(delta_binding_energy) + '\n')
            not_first_line = True


if __name__ == "__main__":
    args = parse_arguments()
    chain_info_dict = read_monomer_info_from_csv(args.monomers, args.chain_ids)
    if args.output_name:
        write_summary_csv(chain_info_dict, args.complex, args.output_name)
    else:
        write_summary_csv(chain_info_dict, args.complex, args.complex[:-4] + '_summary.csv')
