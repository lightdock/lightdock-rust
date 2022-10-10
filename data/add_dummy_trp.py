import os

"""Simple script to add a new atom type to the DFIRE scoring function"""

params = open("DCparams").readlines()
output_file_name = "DCParams.new"

dfire_atom_types = [
    'CYSN', 'CYSCA', 'CYSC', 'CYSO', 'CYSCB', 'CYSSG', 'METN', 'METCA', 'METC', 'METO', 'METCB',
    'METCG', 'METSD', 'METCE', 'PHEN', 'PHECA', 'PHEC', 'PHEO', 'PHECB', 'PHECG', 'PHECD1',
    'PHECD2', 'PHECE1', 'PHECE2', 'PHECZ', 'ILEN', 'ILECA', 'ILEC', 'ILEO', 'ILECB', 'ILECG1',
    'ILECG2', 'ILECD1', 'LEUN', 'LEUCA', 'LEUC', 'LEUO', 'LEUCB', 'LEUCG', 'LEUCD1', 'LEUCD2',
    'VALN', 'VALCA', 'VALC', 'VALO', 'VALCB', 'VALCG1', 'VALCG2', 'TRPN', 'TRPCA', 'TRPC', 'TRPO',
    'TRPCB', 'TRPCG', 'TRPCD1', 'TRPCD2', 'TRPNE1', 'TRPCE2', 'TRPCE3', 'TRPCZ2', 'TRPCZ3',
    'TRPCH2', 'TYRN', 'TYRCA', 'TYRC', 'TYRO', 'TYRCB', 'TYRCG', 'TYRCD1', 'TYRCD2', 'TYRCE1',
    'TYRCE2', 'TYRCZ', 'TYROH', 'ALAN', 'ALACA', 'ALAC', 'ALAO', 'ALACB', 'GLYN', 'GLYCA', 'GLYC',
    'GLYO', 'THRN', 'THRCA', 'THRC', 'THRO', 'THRCB', 'THROG1', 'THRCG2', 'SERN', 'SERCA', 'SERC',
    'SERO', 'SERCB', 'SEROG', 'GLNN', 'GLNCA', 'GLNC', 'GLNO', 'GLNCB', 'GLNCG', 'GLNCD', 'GLNOE1',
    'GLNNE2', 'ASNN', 'ASNCA', 'ASNC', 'ASNO', 'ASNCB', 'ASNCG', 'ASNOD1', 'ASNND2', 'GLUN',
    'GLUCA', 'GLUC', 'GLUO', 'GLUCB', 'GLUCG', 'GLUCD', 'GLUOE1', 'GLUOE2', 'ASPN', 'ASPCA', 'ASPC',
    'ASPO', 'ASPCB', 'ASPCG', 'ASPOD1', 'ASPOD2', 'HISN', 'HISCA', 'HISC', 'HISO', 'HISCB', 'HISCG',
    'HISND1', 'HISCD2', 'HISCE1', 'HISNE2', 'ARGN', 'ARGCA', 'ARGC', 'ARGO', 'ARGCB', 'ARGCG',
    'ARGCD', 'ARGNE', 'ARGCZ', 'ARGNH1', 'ARGNH2', 'LYSN', 'LYSCA', 'LYSC', 'LYSO', 'LYSCB',
    'LYSCG', 'LYSCD', 'LYSCE', 'LYSNZ', 'PRON', 'PROCA', 'PROC', 'PROO', 'PROCB', 'PROCG', 'PROCD',
    'MMBBJ'
]

new_type = "MMYDU"

dfire_atom_types.append(new_type)

num_types = len(set([x[:3] for x in dfire_atom_types]))

z_trp = [10.0, 10.0, 10.0, 10.0, -0.58913422, 0.487351896, 0.154414342, -0.340395132, -0.039789371, -0.189618084, 0.154676419, 0.344108203, 0.345315762, 0.0736003874, -0.0313145271, -0.00786236502, 0.0147609806, -0.0265503577, 0.0097722418, 0.0]

with open(output_file_name, 'w') as out:
    i = 0
    for x in dfire_atom_types:
        for y in dfire_atom_types:
            for z in range(20):
                if x == new_type or y == new_type:
                    p = z_trp[z]
                    out.write(str(z) + '\t' + x + '-' + y + '\t' + str(p) + os.linesep)
                else:
                    out.write(str(z) + '\t' + x + '-' + y + '\t' + str(params[i]))
                    i += 1


