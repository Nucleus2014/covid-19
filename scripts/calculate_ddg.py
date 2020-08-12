import argparse


class Variant:
    def __init__(self, info, chain_ids=None):
        items = info.split(',')
        if items[5] == 'NA':
            self.energy_change = 0
        else:
            self.energy_change = float(items[5])
        self.id_1 = items[6]
        self.id_2 = items[7]
        if len(items[-2]) > 0:
            self.substitutions = items[-2].split(';')
        else:
            self.substitutions = list()
        if chain_ids:
            for chain_id in chain_ids:
                exec('self.monomer_' + chain_id + '=list()')
            for mutation in self.substitutions:
                for chain_id in chain_ids:
                    if mutation[0] == chain_id:
                        exec('self.monomer_' + chain_id + '.append(mutation)')
    
    def match(self, monomer_dict, chain_id):
        mutations_str = eval("';'.join(sorted(self.monomer_" + chain_id + "))")
        monomer = monomer_dict.get(mutations_str)
        if monomer:
            exec('self.monomer_' + chain_id + '=monomer')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--complex', type=str, required=True)
    parser.add_argument('-m', '--monomers', type=str, nargs='*', required=True)
    parser.add_argument('-id', '--chain_ids', type=str, nargs='*', required=True)
    args = parser.parse_args()

    for i, chain_id in enumerate(args.chain_ids):
        exec('monomer_' + chain_id + '_dict=dict()')
        with open(args.monomers[i], 'r') as monomer_csv:
            pass_first_line = False
            for line in monomer_csv:
                if pass_first_line:
                    substitutions = line.split(',')[-2].split(';')
                    key = ';'.join(sorted(substitutions))
                    monomer = Variant(line[:-1])
                    exec('monomer_' + chain_id + '_dict[key]=monomer')
                pass_first_line = True

    with open(args.complex[:-4] + '_ddG.csv', 'w+') as ddg_csv:
        ddg_csv.write('id_1,id_2,mutations,complex_energy,')
        for chain_id in args.chain_ids:
            ddg_csv.write('monomer_' + chain_id + '_energy,')
        ddg_csv.write('ddG\n')

    with open(args.complex, 'r') as cplx_csv:
        pass_first_line = False
        for line in cplx_csv:
            if pass_first_line:
                with open(args.complex[:-4] + '_ddG.csv', 'a+') as ddg_csv:
                    cplx = Variant(line, args.chain_ids)
                    for chain_id in args.chain_ids:
                        exec("cplx.match(monomer_" + chain_id + "_dict,'" + chain_id + "')")

                    ddg_csv.write(cplx.id_1 + ',' + cplx.id_2 + ',' + ';'.join(cplx.substitutions) + \
                        ',' + str(cplx.energy_change) + ',')
                    
                    ddg = cplx.energy_change
                    for chain_id in args.chain_ids:
                        exec('monomer_' + chain_id + '_energy=cplx.monomer_' + chain_id + '.energy_change')
                        exec("ddg_csv.write(str(monomer_" + chain_id + "_energy) + ',')")
                        exec('ddg-=monomer_' + chain_id + '_energy')
                    ddg_csv.write(str(ddg) + '\n')
            pass_first_line = True
