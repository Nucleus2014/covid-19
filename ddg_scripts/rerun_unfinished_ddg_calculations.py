import os
from shutil import copyfile

for mut in filter(lambda x: x.endswith('.mut'), os.listdir()):
    # read mutation sets from mut file
    with open(mut, 'r') as pmut:
        mut_lines = pmut.readlines()
    # reform the mut file
    variants = list()
    variant = list()
    for line in mut_lines[1:]:
        if not len(list(filter(lambda x: x != '', line.split(' ')))) == 3:
            variants.append(variant)
            variant = list()
        else:
            variant.append(line)
    # since there is not a new number line after the last line, 
    # manually append the last variant to variants.
    variants.append(variant)
    # read ddg file
    with open(mut[:-3] + 'ddg', 'r+') as pddg:
        ddg_lines = pddg.readlines()
    # calculate finished variants
    completed_vars = int(len(ddg_lines) / 3)
    if completed_vars != len(variants):
        print(mut)
        if len(ddg_lines) < 6:
            # rerun the whole calculation
            print("Discard " + mut[:-3] + 'ddg and recalculate it.')
            os.remove(mut[:-3] + 'ddg')
            copyfile(mut, mut[:-4] + '-2.mut')
        else:
            # truncate ddg file if unfinished
            if len(ddg_lines) % 3 == 2:
                with open(mut[:-3] + 'ddg', 'a+') as pddg:
                    pddg.write(ddg_lines[-1])
                completed_vars += 1
                if completed_vars == len(variants):
                    print(mut + " is finished!")
                    continue
            elif len(ddg_lines) % 3 == 1:
                os.remove(mut[:-3] + 'ddg')
                with open(mut[:-3] + 'ddg', 'w+') as pddg:
                    pddg.writelines(ddg_lines[:-1])
            # write a new mut file for recalculation
            total_point_mutations = 0
            for variant in variants[completed_vars:]:
                total_point_mutations += len(variant)
            with open(mut[:-4] + '-2.mut', 'w+') as pmut2:
                pmut2.write('total ' + str(total_point_mutations) + '\n')
                for variant in variants[completed_vars:]:
                    pmut2.write(str(len(variant)) + '\n')
                    pmut2.writelines(variant)
